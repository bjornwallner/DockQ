#!/usr/bin/env python

import sys
import os
import pickle
import traceback
import itertools
from functools import lru_cache 
from argparse import ArgumentParser

import Bio.PDB
import numpy as np
from Bio import Align
from Bio.SeqUtils import seq1
from Bio.SVDSuperimposer import SVDSuperimposer

# fallback in case the cython version doesn't work, though it will be slower
try:
    from .operations import residue_distances, get_fnat_stats
except:
    print(
        "WARNING: It looks like cython is not working, falling back on native python. This will make DockQ slower"
    )
    from operations_nocy import residue_distances, get_fnat_stats


def parse_args():
    parser = ArgumentParser(
        description="DockQ - Quality measure for protein-protein docking models"
    )
    parser.add_argument("model", metavar="<model>", type=str, help="path to model file")
    parser.add_argument(
        "native", metavar="<native>", type=str, help="path to native file"
    )
    parser.add_argument(
        "-capri_peptide",
        default=False,
        action="store_true",
        help="use version for capri_peptide (DockQ cannot not be trusted for this setting)",
    )
    parser.add_argument(
        "-short", default=False, action="store_true", help="short output"
    )
    parser.add_argument(
        "-verbose", default=False, action="store_true", help="talk a lot!"
    )
    parser.add_argument(
        "-quiet", default=False, action="store_true", help="keep quiet!"
    )
    parser.add_argument(
        "-useCA", default=False, action="store_true", help="use CA instead of backbone"
    )
    parser.add_argument(
        "-mmcif_model",
        default=False,
        action="store_true",
        help="The model is in mmCIF format",
    )
    parser.add_argument(
        "-mmcif_native",
        default=False,
        action="store_true",
        help="The native is in mmCIF format",
    )
    parser.add_argument(
        "-no_needle",
        default=False,
        action="store_true",
        help="Do not align native and model using sequence alignments, but use the numbering of residues instead",
    )
    parser.add_argument(
        "-perm1",
        default=False,
        action="store_true",
        help="use all chain1 permutations to find maximum DockQ (number of comparisons is n! = 24, if combined with -perm2 there will be n!*m! combinations",
    )
    parser.add_argument(
        "-perm2",
        default=False,
        action="store_true",
        help="use all chain2 permutations to find maximum DockQ (number of comparisons is n! = 24, if combined with -perm1 there will be n!*m! combinations",
    )
    parser.add_argument(
        "-model_chain1",
        metavar="model_chain1",
        type=str,
        nargs="+",
        help="pdb chain order to group together partner 1",
    )
    parser.add_argument(
        "-model_chain2",
        metavar="model_chain2",
        type=str,
        nargs="+",
        help="pdb chain order to group together partner 2 (complement to partner 1 if undef)",
    )
    parser.add_argument(
        "-native_chain1",
        metavar="native_chain1",
        type=str,
        nargs="+",
        help="pdb chain order to group together from native partner 1",
    )
    parser.add_argument(
        "-native_chain2",
        metavar="native_chain2",
        type=str,
        nargs="+",
        help="pdb chain order to group together from native partner 2 (complement to partner 1 if undef)",
    )

    return parser.parse_args()


def calc_DockQ(
    sample_model,
    ref_model,
    ref_model_original,
    group1,
    group2,
    nat_group1,
    nat_group2,
    use_CA_only=False,
    capri_peptide=False,
):
    atom_for_sup = ["CA", "C", "N", "O"] if not use_CA_only else ["CA"]
    fnat_threshold = 4.0 if capri_peptide else 5.0
    interface_threshold = 8.0 if capri_peptide else 10.0

    # total number of native contacts is calculated on untouched native structure
    ref_res_distances = get_residue_distances(
        ref_model_original, nat_group1, nat_group2
    )
    nat_total = np.nonzero(np.asarray(ref_res_distances) < fnat_threshold**2)[
        0
    ].shape[0]
    
    if nat_total == 0:
        # if the native has no interface between the two chain groups
        # nothing to do here
        print(nat_group1, nat_group2)
        return None

    sample_res_distances = get_residue_distances(sample_model, group1, group2)
    ref_res_distances = get_residue_distances(ref_model, nat_group1, nat_group2)

    assert (
        sample_res_distances.shape == ref_res_distances.shape
    ), f"Native and model have incompatible sizes ({sample_res_distances.shape} != {ref_res_distances.shape})"

    nat_correct, nonnat_count, _, model_total = get_fnat_stats(
        sample_res_distances, ref_res_distances, threshold=fnat_threshold
    )
    # avoids divide by 0 errors
    fnat = nat_total and nat_correct / nat_total or 0
    fnonnat = model_total and nonnat_count / model_total or 0

    if capri_peptide:
        ref_res_distances = get_residue_distances(
            ref_model, nat_group1, nat_group2, all_atom=False
        )
    # Get interfacial atoms from reference, and corresponding atoms from sample
    interacting_pairs = get_interacting_pairs(
        # working with squared thresholds to avoid using sqrt in distance calculations
        ref_res_distances,
        threshold=interface_threshold**2,
    )

    # get a copy of each structure, then only keep backbone atoms
    sample_model_backbone = sample_model
    ref_model_backbone = ref_model
    set_common_backbone_atoms(
        sample_model_backbone, ref_model_backbone, atom_types=atom_for_sup
    )

    sample_interface_atoms, ref_interface_atoms = get_interface_atoms(
        interacting_pairs,
        sample_model_backbone,
        ref_model_backbone,
        group1,
        group2,
        nat_group1,
        nat_group2,
    )

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(sample_interface_atoms, ref_interface_atoms)

    irms = super_imposer.rms

    # assign which group of chains constitutes the receptor, then the other is the ligand
    ref_group1_size = np.sum([len(ref_model_original[chain]) for chain in nat_group1])
    ref_group2_size = np.sum([len(ref_model_original[chain]) for chain in nat_group2])
    receptor_chains = (
        (nat_group1, group1)
        if ref_group1_size > ref_group2_size
        else (nat_group2, group2)
    )
    ligand_chains = (
        (nat_group1, group1)
        if ref_group1_size <= ref_group2_size
        else (nat_group2, group2)
    )
    class1, class2 = (
        ("receptor", "ligand")
        if ref_group1_size > ref_group2_size
        else ("ligand", "receptor")
    )
    receptor_atoms_native = [
        atom
        for chain in receptor_chains[0]
        for atom in ref_model_backbone[chain].get_atoms()
    ]
    receptor_atoms_sample = [
        atom
        for chain in receptor_chains[1]
        for atom in sample_model_backbone[chain].get_atoms()
    ]

    # Set to align on receptor
    super_imposer.set_atoms(receptor_atoms_native, receptor_atoms_sample)
    super_imposer.apply(sample_model_backbone.get_atoms())

    coord1 = np.array(
        [
            atom.coord
            for chain in ligand_chains[0]
            for atom in ref_model_backbone[chain].get_atoms()
        ]
    )
    coord2 = np.array(
        [
            atom.coord
            for chain in ligand_chains[1]
            for atom in sample_model_backbone[chain].get_atoms()
        ]
    )

    sup = SVDSuperimposer()
    Lrms = sup._rms(
        coord1, coord2
    )  # using the private _rms function which does not superimpose

    DockQ = (
        float(fnat)
        + 1 / (1 + (irms / 1.5) * (irms / 1.5))
        + 1 / (1 + (Lrms / 8.5) * (Lrms / 8.5))
    ) / 3
    info = {}
    info["DockQ"] = DockQ
    info["irms"] = irms
    info["Lrms"] = Lrms
    info["fnat"] = fnat
    info["nat_correct"] = nat_correct
    info["nat_total"] = nat_total

    info["fnonnat"] = fnonnat
    info["nonnat_count"] = nonnat_count
    info["model_total"] = model_total

    info["chain1"] = " ".join(group1)
    info["chain2"] = " ".join(group2)
    info["len1"] = ref_group1_size
    info["len2"] = ref_group2_size
    info["class1"] = class1
    info["class2"] = class2

    return info


def align_model_to_native(
    model_structure, native_structure, model_chain, native_chain, use_numbering=False
):
    """
    Function to align two PDB structures. This can be done by sequence (default) or by
    numbering. If the numbering is used, then each residue number from the pdb structure
    is converted to a unique character. Then the two vectors of character are aligned
    as if they were two sequences
    """
    
    if use_numbering:
        model_numbering = []
        native_numbering = []

        for residue in model_structure[model_chain].get_residues():
            resn = int(residue.id[1])
            model_numbering.append(resn)

        for residue in native_structure[native_chain].get_residues():
            resn = int(residue.id[1])
            native_numbering.append(resn)
        # if the samllest resn is negative, it will be used to shift all numbers so they start from 0
        # the minimum offset is 45 to avoid including the "-" character that is reserved for gaps
        min_resn = max(45, -min(model_numbering + native_numbering))

        model_sequence = "".join([chr(resn + min_resn) for resn in model_numbering])
        native_sequence = "".join([chr(resn + min_resn) for resn in native_numbering])

    else:
        model_sequence = "".join(
            seq1(residue.get_resname())
            for residue in model_structure[model_chain].get_residues()
        )

        native_sequence = "".join(
            seq1(residue.get_resname())
            for residue in native_structure[native_chain].get_residues()
        )

    aligner = Align.PairwiseAligner()
    aligner.match = 5
    aligner.mismatch = 0
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aln = aligner.align(model_sequence, native_sequence)[0]
    return aln


def format_alignment(aln):
    alignment = {}
    formatted_aln = aln.format().split("\n")
    alignment["seqA"] = formatted_aln[0]
    alignment["matches"] = formatted_aln[1]
    alignment["seqB"] = formatted_aln[2]
    return alignment


def remove_extra_chains(model, chains_to_keep):
    chains = [chain.id for chain in model.get_chains()]

    chains_to_remove = set(chains).difference(set(chains_to_keep))
    for chain in chains_to_remove:
        model.detach_child(chain)


def remove_hetatms(model):
    chains = [chain.id for chain in model.get_chains()]
    residues_to_delete = []

    for chain in chains:
        residues = model[chain].get_residues()

        for res in residues:
            if res.id[0] != " ":
                residues_to_delete.append(res.get_full_id())
    for _, _, chain, res in residues_to_delete:
        model[chain].detach_child(res)


def fix_chain_residues(model, chain, alignment, invert=False):
    residues = model[chain].get_residues()
    residues_to_delete = []

    seqA = alignment["seqA"] if not invert else alignment["seqB"]
    seqB = alignment["seqB"] if not invert else alignment["seqA"]
    for (aligned_residue_A, aligned_residue_B) in zip(seqA, seqB):
        if aligned_residue_A != "-":
            residue = next(residues)
        if aligned_residue_B == "-":  # gap residue: remove from structure
            residues_to_delete.append(residue.get_full_id())

    for _, _, _, res in residues_to_delete:
        model[chain].detach_child(res)


def list_atoms_per_residue(model, group):
    n_atoms_per_residue = []
    for chain in group:
        for residue in model[chain].get_residues():
            # important to remove duplicate atoms (e.g. alternates) at this stage, remove also hydrogens
            atom_ids = set(
                [a.id for a in residue.get_unpacked_list() if a.element != "H"]
            )
            n_atoms_per_residue.append(len(atom_ids))
    return np.array(n_atoms_per_residue).astype(int)


#@lru_cache
def get_residue_distances(model, group1, group2, all_atom=True):
    if all_atom:
        # get information about how many atoms correspond to each amino acid in each group of chains
        n_atoms_per_res_group1 = list_atoms_per_residue(model, group1)
        n_atoms_per_res_group2 = list_atoms_per_residue(model, group2)
        model_A_atoms = np.asarray(
            [
                atom.get_coord()
                for chain in group1
                for res in model[chain].get_residues()
                for atom in res.get_atoms()
                if atom.element != "H"
            ]
        )
        model_B_atoms = np.asarray(
            [
                atom.get_coord()
                for chain in group2
                for res in model[chain].get_residues()
                for atom in res.get_atoms()
                if atom.element != "H"
            ]
        )

    else:  # distances were already between CBs only
        model_A_atoms = np.asarray(
            [
                res["CB"].get_coord() if "CB" in res else res["CA"].get_coord()
                for chain in group1
                for res in model[chain].get_residues()
            ]
        )
        model_B_atoms = np.asarray(
            [
                res["CB"].get_coord() if "CB" in res else res["CA"].get_coord()
                for chain in group2
                for res in model[chain].get_residues()
            ]
        )

        n_atoms_per_res_group1 = np.ones(model_A_atoms.shape[0]).astype(int)
        n_atoms_per_res_group2 = np.ones(model_B_atoms.shape[0]).astype(int)

    model_res_distances = residue_distances(
        model_A_atoms, model_B_atoms, n_atoms_per_res_group1, n_atoms_per_res_group2
    )
    return model_res_distances


def get_interacting_pairs(distances, threshold):
    return np.nonzero(np.asarray(distances) < threshold)


# @profile
def get_interface_atoms(
    interacting_pairs,
    model_backbone,
    ref_backbone,
    model_group1,
    model_group2,
    ref_group1,
    ref_group2,
):
    ref_interface = []
    mod_interface = []

    ref_residues_group1 = [
        res for chain in ref_group1 for res in ref_backbone[chain].get_residues()
    ]
    ref_residues_group2 = [
        res for chain in ref_group2 for res in ref_backbone[chain].get_residues()
    ]

    mod_residues_group1 = [
        res for chain in model_group1 for res in model_backbone[chain].get_residues()
    ]
    mod_residues_group2 = [
        res for chain in model_group2 for res in model_backbone[chain].get_residues()
    ]

    # remove duplicate residues
    interface_residues_group1 = set(interacting_pairs[0])
    interface_residues_group2 = set(interacting_pairs[1])

    for i in interface_residues_group1:
        ref_interface += [atom for atom in ref_residues_group1[i].get_atoms()]
        mod_interface += [atom for atom in mod_residues_group1[i].get_atoms()]

    for j in interface_residues_group2:
        ref_interface += [atom for atom in ref_residues_group2[j].get_atoms()]
        mod_interface += [atom for atom in mod_residues_group2[j].get_atoms()]

    return mod_interface, ref_interface


def set_common_backbone_atoms(model, reference, atom_types=["CA", "C", "N", "O"]):
    # model and reference should have the same number of amino acids and be aligned
    for mod_res, ref_res in zip(model.get_residues(), reference.get_residues()):
        mod_atoms = [atom for atom in mod_res.get_atoms()]
        ref_atoms = [atom for atom in ref_res.get_atoms()]

        atom_ids_in_mod_res = [atm.id for atm in mod_atoms]
        atom_ids_in_ref_res = [atm.id for atm in ref_atoms]

        atom_ids_in_ref_and_mod_res = (
            set(atom_ids_in_mod_res)
            .intersection(atom_types)
            .intersection(atom_ids_in_ref_res)
        )
        # whatever atom is not in the shared list, remove it from the both structures
        for atom_id in set(atom_ids_in_mod_res).difference(atom_ids_in_ref_and_mod_res):
            mod_res.detach_child(atom_id)

        for atom_id in set(atom_ids_in_ref_res).difference(atom_ids_in_ref_and_mod_res):
            ref_res.detach_child(atom_id)


#@lru_cache
def run_on_groups(
    model_structure,
    native_structure,
    group1,
    group2,
    nat_group1,
    nat_group2,
    no_needle=False,
    use_CA_only=False,
    capri_peptide=False,
):
    remove_extra_chains(model_structure, chains_to_keep=group1 + group2)
    native_structure_original = pickle.loads(pickle.dumps(native_structure, -1))
    remove_extra_chains(native_structure, chains_to_keep=nat_group1 + nat_group2)

    # realign each model chain against the corresponding native chain
    for model_chain, native_chain in zip(group1 + group2, nat_group1 + nat_group2):
        aln = align_model_to_native(
            model_structure,
            native_structure,
            model_chain,
            native_chain,
            use_numbering=no_needle,
        )
        alignment = format_alignment(aln)
        fix_chain_residues(model_structure, model_chain, alignment)
        fix_chain_residues(native_structure, native_chain, alignment, invert=True)
    info = calc_DockQ(
        model_structure,
        native_structure,
        native_structure_original,
        group1,
        group2,
        nat_group1,
        nat_group2,
        use_CA_only=use_CA_only,
        capri_peptide=capri_peptide,
    )
    return info


def run_on_all_native_interfaces(
    model_structure,
    native_structure,
    chain_map={"A": "A", "B": "B"},
    no_needle=False,
    use_CA_only=False,
    capri_peptide=False,
):
    """Given a native-model chain map, finds all non-null native interfaces and runs DockQ for each native-model pair of interfaces"""
    results_dic = {}
    native_chains = [c.id for c in native_structure]
    for chain_pair in itertools.combinations(native_chains, 2):
        interface_size = np.sum(
            np.asarray(
                get_residue_distances(
                    native_structure, (chain_pair[0]), (chain_pair[1])
                )
            )
            < 25.0
        )

        if (
            interface_size > 0
            and chain_pair[0] in chain_map
            and chain_pair[1] in chain_map
        ):
            model_structure_this = pickle.loads(pickle.dumps(model_structure, -1))
            native_structure_this = pickle.loads(pickle.dumps(native_structure, -1))
            info = run_on_groups(
                model_structure_this,
                native_structure_this,
                (chain_map[chain_pair[0]]),
                (chain_map[chain_pair[1]]),
                (chain_pair[0]),
                (chain_pair[1]),
            )
            results_dic[chain_pair] = info

    return results_dic


def run_DockQ(
    path_to_model,
    path_to_native,
    group1=["A"],
    group2=["B"],
    nat_group1=["A"],
    nat_group2=["B"],
    model_is_mmcif=False,
    native_is_mmcif=False,
    no_needle=False,
    use_CA_only=False,
    capri_peptide=False,
):
    model = load_PDB(path_to_model, is_mmcif=model_is_mmcif)
    native = load_PDB(path_to_native, is_mmcif=native_is_mmcif)

    return run_on_groups(
        model,
        native,
        group1,
        group2,
        nat_group1,
        nat_group2,
        no_needle,
        use_CA_only,
        capri_peptide,
    )


def load_PDB(path, n_model=0, is_mmcif=False):

    if not is_mmcif:
        pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    else:
        pdb_parser = Bio.PDB.MMCIFParser(QUIET=True)

    try:
        structure = pdb_parser.get_structure("-", path)
        model = structure[n_model]
    except Exception as e:
        print("ERROR: is the file in the correct format? (.pdb, .mmcif)")
        if not is_mmcif:
            print("       (use -mmcif_model or -mmcif_native with mmCIF inputs)")
        print(traceback.format_exc())
        sys.exit(1)
    return model


def group_model_chains(model_structure, native_structure, model_chains, native_chains):
    alignment_targets = itertools.product(model_chains, native_chains)
    native_chain_clusters = {chain:[] for chain in native_chains}

    for model_chain, native_chain in alignment_targets:
        aln = align_model_to_native(model_structure, native_structure, model_chain, native_chain)
        alignment = format_alignment(aln)
        if "." not in alignment["matches"] and ("-" not in alignment["seqA"] or "-" not in alignment["seqB"]):
            # 100% sequence identity, 100% coverage of native sequence in model sequence
            native_chain_clusters[native_chain].append(model_chain)
    return native_chain_clusters

# @profile
def main():
    args = parse_args()

    bio_ver = 1.79
    if float(Bio.__version__) < bio_ver:
        print(
            f"WARNING: Biopython version {Bio.__version__} is older than the recommended version {bio_ver}"
        )

    native_structure = load_PDB(args.native, is_mmcif=args.mmcif_native)
    model_structure = load_PDB(args.model, is_mmcif=args.mmcif_model)
    remove_hetatms(native_structure)
    remove_hetatms(model_structure)

    best_info = ""
    info = {}
    model_chains = [c.id for c in model_structure]
    native_chains = [c.id for c in native_structure]

    if (len(model_chains) > 2 or len(native_chains) > 2) and (
        not args.model_chain1 and not args.native_chain1
    ):
        print(
            "Multi-chain model need sets of chains to group\nuse -native_chain1 and/or -model_chain1 if you want a different mapping than 1-1"
        )
        print(f"Model chains: {' '.join(model_chains)}")
        print(f"Native chains: {' '.join(native_chains)}")
        sys.exit()
    if len(model_chains) < 2 or len(native_chains) < 2:
        print("Need at least two chains in the two inputs\n")
        sys.exit()

    # Some of these might be None
    group1 = [model_chains[0]] if len(model_chains) == 2 else args.model_chain1
    group2 = args.model_chain2
    nat_group1 = args.native_chain1
    nat_group2 = args.native_chain2
    
    # at this stage either group1 or nat_group1 are not None
    if not nat_group1:  # then the user has set group1. Try to follow the same mapping
        if (
            model_chains == native_chains
        ):  # easier case: the chains have the same naming between native/model, so just copy them
            nat_group1 = group1
            # use complement to nat_group1 if group2 hasn't been decided yet
            nat_group2 = group2
        else:  # otherwise, group the chains by however many where in either model group
            nat_group1 = native_chains[: len(group1)]
            nat_group2 = (
                native_chains[len(group1) : len(group1) + len(group2)]
                if group2
                else None
            )

    if not group1:  # viceversa, the user has set nat_group1
        if model_chains == native_chains:
            group1 = nat_group1
            group2 = nat_group2
        else:
            group1 = model_chains[: len(nat_group1)]
            group2 = (
                model_chains[len(nat_group1) : len(nat_group1) + len(nat_group2)]
                if nat_group2
                else None
            )

    if not group2:  # no group2 set yet, use the complement to group1
        group2 = [chain for chain in model_chains if chain not in group1]
    if not nat_group2:
        nat_group2 = [chain for chain in native_chains if chain not in nat_group1]

    if not args.perm1 and not args.perm2:
        info = run_on_groups(
            model_structure,
            native_structure,
            group1,
            group2,
            nat_group1,
            nat_group2,
            args.no_needle,
            args.useCA,
            args.capri_peptide,
        )
    else:  # permute chains and run on a for loop
        best_dockq = -1
        best_mapping = None

        native_chain_clusters = group_model_chains(model_structure, native_structure, model_chains, native_chains)
        all_mappings = itertools.product(*native_chain_clusters.values())
        
        # remove mappings where the same model chain is present more than once
        all_mappings = [element for element in all_mappings if len(set(element)) == len(element)]
        combo_dockq = -1
        for mapping in all_mappings:
            chain_map = {native_chain:mapping[i] for i, native_chain in enumerate(native_chains)}

            result_this_mapping = run_on_all_native_interfaces(
                model_structure,
                native_structure,
                chain_map=chain_map,
                no_needle=args.no_needle,
                use_CA_only=args.useCA,
                capri_peptide=args.capri_peptide,
            )
            #print(chain_map)
            #print(result_this_mapping.values())
            total_dockq = sum([result["DockQ"] for result in result_this_mapping.values()])
            print(result_this_mapping, total_dockq)
            if total_dockq > best_dockq:
                best_dockq = total_dockq
                best_result = result_this_mapping

        info["model"] = args.model
        info["native"] = args.native
        info["best_dockq"] = best_dockq
        info["best_result"] = best_result

    print_results(info, args.short, args.capri_peptide)


def print_results(info, short=False, capri_peptide=False):
    if short:
        capri_peptide_str = "-capri_peptide" if capri_peptide else ""
        print(
            f"DockQ{capri_peptide_str} {info['DockQ']:.3f} Fnat {info['fnat']:.3f} iRMS {info['irms']:.3f} LRMS {info['Lrms']:.3f} Fnonnat {info['fnonnat']:.3f} {info['model']} {info['native']} {info['best']}"
        )

    else:
        if capri_peptide:
            print("****************************************************************")
            print("*                DockQ-CAPRI peptide                           *")
            print("*   Do not trust any thing you read....                        *")
            print("*   OBS THE DEFINITION OF Fnat and iRMS are different for      *")
            print("*   peptides in CAPRI                                          *")
            print("*                                                              *")
            print("*   For the record:                                            *")
            print("*   Definition of contact <4A all heavy atoms (Fnat)           *")
            print("*   Definition of interface <8A CB (iRMS)                      *")
            print("*   For comments, please email: bjorn.wallner@.liu.se          *")
            print("****************************************************************")
        else:
            print("****************************************************************")
            print("*                       DockQ                                  *")
            print("*   Scoring function for protein-protein docking models        *")
            print("*   Statistics on CAPRI data:                                  *")
            print("*    0.00 <= DockQ <  0.23 - Incorrect                         *")
            print("*    0.23 <= DockQ <  0.49 - Acceptable quality                *")
            print("*    0.49 <= DockQ <  0.80 - Medium quality                    *")
            print("*            DockQ >= 0.80 - High quality                      *")
            print("*   Ref: S. Basu and B. Wallner, DockQ: A quality measure for  *")
            print("*   protein-protein docking models                             *")
            print("*                            doi:10.1371/journal.pone.0161879  *")
            print("*   For the record:                                            *")
            print("*   Definition of contact <5A (Fnat)                           *")
            print("*   Definition of interface <10A all heavy atoms (iRMS)        *")
            print("*   For comments, please email: bjorn.wallner@.liu.se          *")
            print("*                                                              *")
            print("****************************************************************")
        print(f"Model  : {info['model']}")
        print(f"Native : {info['native']}")
        if "best_dockq" in info:
            print(info["best_result"])
            print(info["best_dockq"])
        else:
            print(
                f"Number of equivalent residues in chain {info['chain1']} {info['len1']} ({info['class1']})"
            )
            print(
                f"Number of equivalent residues in chain {info['chain2']} {info['len2']} ({info['class2']})"
            )
            print(
                f"Fnat {info['fnat']:.3f} {info['nat_correct']} correct of {info['nat_total']} native contacts"
            )
            print(
                f"Fnonnat {info['fnonnat']:.3f} {info['nonnat_count']} non-native of {info['model_total']} model contacts"
            )
            print(f"iRMS {info['irms']:.3f}")
            print(f"LRMS {info['Lrms']:.3f}")

            peptide_disclaimer = (
                " DockQ not reoptimized for CAPRI peptide evaluation"
                if capri_peptide
                else ""
            )
            print(f"DockQ {info['DockQ']:.3f}{peptide_disclaimer}")


if __name__ == "__main__":
    main()
