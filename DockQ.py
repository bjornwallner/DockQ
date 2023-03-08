#!/usr/bin/env python

import sys
import os
import pickle
import itertools
import numpy as np
from argparse import ArgumentParser
import Bio.PDB
from Bio import pairwise2
from Bio.SeqUtils import seq1
from Bio.SVDSuperimposer import SVDSuperimposer


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
        "-mmcif_model", default=False, action="store_true", help="The model is in mmCIF format"
    )
    parser.add_argument(
        "-mmcif_native", default=False, action="store_true", help="The native is in mmCIF format"
    )
    parser.add_argument(
        "-no_needle",
        default=False,
        action="store_true",
        help="do not use global alignment to fix residue numbering between native and model during chain permutation (use only in case needle is not installed, and the residues between the chains are identical",
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


def capri_class(fnat, iRMS, LRMS, capri_peptide=False):
    if capri_peptide:

        if fnat < 0.2 or (LRMS > 5.0 and iRMS > 2.0):
            return "Incorrect"
        elif (
            (fnat >= 0.2 and fnat < 0.5)
            and (LRMS <= 5.0 or iRMS <= 2.0)
            or (fnat >= 0.5 and LRMS > 2.0 and iRMS > 1.0)
        ):
            return "Acceptable"
        elif (
            (fnat >= 0.5 and fnat < 0.8)
            and (LRMS <= 2.0 or iRMS <= 1.0)
            or (fnat >= 0.8 and LRMS > 1.0 and iRMS > 0.5)
        ):
            return "Medium"
        elif fnat >= 0.8 and (LRMS <= 1.0 or iRMS <= 0.5):
            return "High"
        else:
            return "Undef"
    else:

        if fnat < 0.1 or (LRMS > 10.0 and iRMS > 4.0):
            return "Incorrect"
        elif (
            (fnat >= 0.1 and fnat < 0.3)
            and (LRMS <= 10.0 or iRMS <= 4.0)
            or (fnat >= 0.3 and LRMS > 5.0 and iRMS > 2.0)
        ):
            return "Acceptable"
        elif (
            (fnat >= 0.3 and fnat < 0.5)
            and (LRMS <= 5.0 or iRMS <= 2.0)
            or (fnat >= 0.5 and LRMS > 1.0 and iRMS > 1.0)
        ):
            return "Medium"
        elif fnat >= 0.5 and (LRMS <= 1.0 or iRMS <= 1.0):
            return "High"
        else:
            return "Undef"


def capri_class_DockQ(DockQ, capri_peptide=False):
    if capri_peptide:
        return "Undef for capri_peptides"

    (c1, c2, c3) = (0.23, 0.49, 0.80)
    if DockQ < c1:
        return "Incorrect"
    elif DockQ >= c1 and DockQ < c2:
        return "Acceptable"
    elif DockQ >= c2 and DockQ < c3:
        return "Medium"
    elif DockQ >= c3:
        return "High"
    else:
        return "Undef"


#@profile
def calc_DockQ(
    sample_model,
    ref_model,
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

    sample_res_distances = get_residue_distances(sample_model, group1, group2)
    ref_res_distances = get_residue_distances(ref_model, nat_group1, nat_group2)

    nat_correct, nonnat_count, nat_total, model_total = get_fnat_stats(
        sample_res_distances, ref_res_distances, thr=fnat_threshold
    )
    fnat = nat_correct / nat_total
    fnonnat = nonnat_count / model_total

    # get a copy of each structure, then only keep backbone atoms. This is faster than copy.deepcopy()
    sample_model_backbone = pickle.loads(pickle.dumps(sample_model, -1))
    ref_model_backbone = pickle.loads(pickle.dumps(ref_model, -1))
    set_common_backbone_atoms(
        sample_model_backbone, ref_model_backbone, atom_types=atom_for_sup
    )

    if capri_peptide:
        ref_res_distances = get_residue_distances(
            ref_model, nat_group1, nat_group2, all_atom=False
        )
    # Get interfacial atoms from reference, and corresponding atoms from sample
    interacting_pairs = get_interacting_pairs(
        ref_res_distances, thr=interface_threshold
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
    ref_group1_size = np.sum([len(ref_model[chain]) for chain in nat_group1])
    ref_group2_size = np.sum([len(ref_model[chain]) for chain in nat_group2])
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
    receptor_chain_rms = super_imposer.rms

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
    model_sequence = "".join(
        seq1(residue.get_resname())
        for residue in model_structure[model_chain].get_residues()
    )

    native_sequence = "".join(
        seq1(residue.get_resname())
        for residue in native_structure[native_chain].get_residues()
    )
    alignment = {}
    if use_numbering:
        model_numbering = []
        native_numbering = []

        for residue in model_structure[model_chain].get_residues():
            resn = int(residue.id[1])
            model_numbering.append(resn)

        for residue in native_structure[native_chain].get_residues():
            resn = int(residue.id[1])
            native_numbering.append(resn)

        start = min(native_numbering + model_numbering)
        end = max(native_numbering + model_numbering)
        alignment["seqA"] = []
        alignment["seqB"] = []
        model_sequence = iter(model_sequence)
        native_sequence = iter(native_sequence)

        for i in range(start, end + 1):
            if i in model_numbering:
                next_model_res = next(model_sequence)
            else:
                next_model_res = "-"

            if i in native_numbering:
                next_native_res = next(native_sequence)
            else:
                next_native_res = "-"

            if next_model_res != "-" or next_native_res != "-":
                alignment["seqA"].append(next_model_res)
                alignment["seqB"].append(next_native_res)

    else:
        aln = pairwise2.align.localms(
            model_sequence, native_sequence, match=5, mismatch=0, open=-10, extend=-1
        )[0]
        alignment["seqA"] = aln.seqA
        alignment["seqB"] = aln.seqB

    return alignment


def remove_extra_chains(model, chains_to_keep):
    for chain in model.get_chains():
        if chain.id not in chains_to_keep:
            model.detach_child(chain.id)
    return model


def fix_chain_residues(model, chain, alignment, invert=False):
    residues = model[chain].get_residues()
    residues_to_delete = []
    start = False
    seqA = alignment["seqA"] if not invert else alignment["seqB"]
    seqB = alignment["seqB"] if not invert else alignment["seqA"]
    for (aligned_residue_A, aligned_residue_B) in zip(seqA, seqB):
        if aligned_residue_A != "-":
            residue = next(residues)
        if aligned_residue_B == "-":  # gap residue: remove from structure
            residues_to_delete.append(residue.get_full_id())

    for _, _, _, res in residues_to_delete:
        model[chain].detach_child(res)


def get_distances_across_chains(model, group1, group2, all_atom=True):
    if all_atom:
        model_A_atoms = np.asarray(
            [
                atom.get_coord()
                for chain in group1
                for res in model[chain].get_residues()
                for atom in res.get_atoms()
            ]
        )
        model_B_atoms = np.asarray(
            [
                atom.get_coord()
                for chain in group2
                for res in model[chain].get_residues()
                for atom in res.get_atoms()
            ]
        )
    else:
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

    distances = np.sqrt(
        ((model_A_atoms[:, None] - model_B_atoms[None, :]) ** 2).sum(-1)
    )
    return distances


# @profile
def atom_distances_to_residue_distances(atom_distances, group_dic1, group_dic2):
    res_distances = np.zeros((len(group_dic1), len(group_dic2)))

    cum_i_atoms = 0
    for i, i_atoms in enumerate(group_dic1):
        cum_j_atoms = 0
        for j, j_atoms in enumerate(group_dic2):
            res_distances[i, j] = atom_distances[
                cum_i_atoms : cum_i_atoms + i_atoms, cum_j_atoms : cum_j_atoms + j_atoms
            ].min()
            cum_j_atoms += j_atoms
        cum_i_atoms += i_atoms
    return res_distances


def list_atoms_per_residue(model, group):
    n_atoms_per_residue = []
    for chain in group:
        for residue in model[chain].get_residues():
            # important to remove duplicate atoms (e.g. alternates) at this stage
            atom_ids = set([a.id for a in residue.get_unpacked_list()])
            n_atoms_per_residue.append(len(atom_ids))
    return n_atoms_per_residue


# @profile
def get_residue_distances(structure, group1, group2, all_atom=True):
    # get information about how many atoms correspond to each amino acid in each group of chains
    n_atoms_per_res_group1 = list_atoms_per_residue(structure, group1)
    n_atoms_per_res_group2 = list_atoms_per_residue(structure, group2)

    model_atom_distances = get_distances_across_chains(
        structure, group1, group2, all_atom=all_atom
    )

    if all_atom:
        model_res_distances = atom_distances_to_residue_distances(
            model_atom_distances, n_atoms_per_res_group1, n_atoms_per_res_group2
        )
    else:  # distances were already between CBs only
        model_res_distances = model_atom_distances
    return model_res_distances


def get_fnat_stats(model_res_distances, native_res_distances, thr=5.0):
    native_contacts = native_res_distances < thr
    model_contacts = model_res_distances < thr
    n_native_contacts = np.sum(native_contacts)
    n_model_contacts = np.sum(model_contacts)
    n_shared_contacts = np.sum(model_contacts * native_contacts)
    n_non_native_contacts = np.sum(model_contacts * (1 - native_contacts))
    return (
        n_shared_contacts,
        n_non_native_contacts,
        n_native_contacts,
        n_model_contacts,
    )


def get_interacting_pairs(distances, thr=5.0):
    return np.nonzero(distances < thr)


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


def run_on_groups(model_structure, native_structure, group1, group2, nat_group1, nat_group2, no_needle=False, use_CA_only=False, capri_peptide=False):
    model_structure = remove_extra_chains(
        model_structure, chains_to_keep=group1 + group2
    )
    native_structure = remove_extra_chains(
        native_structure, chains_to_keep=nat_group1 + nat_group2
    )

    # realign each model chain against the corresponding native chain
    for model_chain, native_chain in zip(
        group1 + group2, nat_group1 + nat_group2
    ):
        alignment = align_model_to_native(
            model_structure,
            native_structure,
            model_chain,
            native_chain,
            use_numbering=no_needle,
        )
        fix_chain_residues(model_structure, model_chain, alignment)
        fix_chain_residues(
            native_structure, native_chain, alignment, invert=True
        )
    info = calc_DockQ(
        model_structure,
        native_structure,
        group1,
        group2,
        nat_group1,
        nat_group2,
        use_CA_only=use_CA_only,
        capri_peptide=capri_peptide,
    )
    return info


def load_PDB(path, n_model=0, is_mmcif=False):

    if not is_mmcif:
        pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    else:
        pdb_parser = Bio.PDB.MMCIFParser(QUIET=True)

    try:
        structure = pdb_parser.get_structure("-", path)
        model = structure[n_model]
    except Exception as e:
        print(e)
    return model

#@profile
def main():
    args = parse_args()

    bio_ver = 1.79
    if float(Bio.__version__) < bio_ver:
        print(
            f"WARNING: Biopython version {Bio.__version__} is older than the recommended version {bio_ver}"
        )

    native_structure = load_PDB(args.native, is_mmcif=args.mmcif_native)
    model_structure = load_PDB(args.model, is_mmcif=args.mmcif_model)

    best_info = ""

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
    if not nat_group1: # then the user has set group1. Try to follow the same mapping
        if model_chains == native_chains: # easier case: the chains have the same naming between native/model, so just copy them
            nat_group1 = group1
            # use complement to nat_group1 if group2 hasn't been decided yet
            nat_group2 = group2
        else: # otherwise, group the chains by however many where in either model group
            nat_group1 = native_chains[: len(group1)]
            nat_group2 = native_chains[len(group1) :]
            
    if not group1: # viceversa, the user has set nat_group1
        if model_chains == native_chains:
            group1 = nat_group1
            group2 = nat_group2
        else:
            group1 = model_chains[: len(nat_group1)]
            group2 = model_chains[len(nat_group1) :]
            
    if not group2: # no group2 set yet, use the complement to group1
        group2 = [chain for chain in model_chains if chain not in group1]
    if not nat_group2:
        nat_group2 = [chain for chain in native_chains if chain not in nat_group1]
        
    if not args.perm1 and not args.perm2:
        info = run_on_groups(model_structure, native_structure, group1, group2, nat_group1, nat_group2, args.no_needle, args.useCA, args.capri_peptide)
    else: # permute chains and run on a for loop
        pe = 0
        best_DockQ = -1
        best_g1 = []
        best_g2 = []

        iter_perm1 = itertools.combinations(group1, len(group1))
        iter_perm2 = itertools.combinations(group2, len(group2))
        if args.perm1:
            iter_perm1 = itertools.permutations(group1)
        if args.perm2:
            iter_perm2 = itertools.permutations(group2)

        combos1 = []
        combos2 = []
        for g1 in iter_perm1:  # _temp:
            combos1.append(g1)
        for g2 in iter_perm2:
            combos2.append(g2)

        for g1 in combos1:
            for g2 in combos2:
                pe = pe + 1
        pe_tot = pe
        pe = 1
        if args.verbose:
            print(
                f"Starting chain order permutation search (number of permutations: {pe_tot})"
            )

        for g1 in combos1:
            for g2 in combos2:

                model_structure_this = pickle.loads(
                    pickle.dumps(model_structure, -1)
                )
                model_structure_this = remove_extra_chains(
                    model_structure_this, chains_to_keep=g1 + g2
                )
                test_info = run_on_groups(model_structure_this, native_structure, g1, g2, nat_group1, nat_group2, args.no_needle, args.useCA, args.capri_peptide)

                if not args.quiet:
                    print(
                        f"{pe}/{pe_tot} {''.join(g1)} -> {''.join(g2)} {test_info['DockQ']}"
                    )

                if test_info["DockQ"] > best_DockQ:
                    best_DockQ = test_info["DockQ"]
                    info = test_info
                    best_g1 = g1
                    best_g2 = g2
                    best_info = f"Best score ({best_DockQ}) found for model -> native, chain1: {''.join(best_g1)} -> {''.join(nat_group1)} chain2: {''.join(best_g2)} -> {''.join(nat_group2)}"

                    if args.verbose:
                        print(best_info)
                    if not args.quiet:
                        print(f"Current best: {best_DockQ}")
                pe = pe + 1
        if not args.quiet:
            print(best_info)

    info["model"] = args.model
    info["native"] = args.native
    info["best"] = best_info
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
        if "best" in info:
            print(info["best"])
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
        peptide_suffix = ""
        if capri_peptide:
            peptide_suffix = "_peptide"

        peptide_disclaimer = (
            " DockQ not reoptimized for CAPRI peptide evaluation"
            if capri_peptide
            else ""
        )
        print(f"DockQ {info['DockQ']:.3f}{peptide_disclaimer}")


if __name__ == "__main__":
    main()
