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
    parser.add_argument(
        "model", metavar="<model>", type=str, help="path to model file"
    )
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
        "-skip_check",
        default=False,
        action="store_true",
        help="skip initial check fo speed up on two chain examples",
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
def calc_DockQ(sample_model, ref_model, group1, group2, nat_group1, nat_group2, use_CA_only=False, capri_peptide=False):
    atom_for_sup = ["CA", "C", "N", "O"] if not use_CA_only else ["CA"]
    threshold = 4.0 if capri_peptide else 5.0
    interface_threshold = 8.0 if capri_peptide else 10.0
    all_atom = not capri_peptide
    
    sample_res_distances = get_residue_distances(sample_model, group1, group2)
    ref_res_distances = get_residue_distances(ref_model, nat_group1, nat_group2)

    nat_correct, nonnat_count, nat_total, model_total = get_fnat(sample_res_distances, ref_res_distances, thr=threshold)
    fnat = nat_correct / nat_total
    fnonnat = nonnat_count / model_total

    # get a copy of each structure containing shared backbone atoms
    sample_model_backbone = pickle.loads(pickle.dumps(sample_model, -1))
    ref_model_backbone = pickle.loads(pickle.dumps(ref_model, -1))
    set_common_backbone_atoms(sample_model_backbone, ref_model_backbone, atom_types=atom_for_sup)

    # Get interfacial atoms from reference, and corresponding atoms from sample
    interacting_pairs = get_interacting_pairs(ref_res_distances, thr=interface_threshold)
    sample_interface_atoms, ref_interface_atoms = get_interface_atoms(interacting_pairs,
                                                                      sample_model_backbone, 
                                                                      ref_model_backbone, 
                                                                      group1, 
                                                                      group2, 
                                                                      nat_group1, 
                                                                      nat_group2)
 
    super_imposer = Bio.PDB.Superimposer()
    #print([atom.coord for atom in ref_interface_atoms])
    super_imposer.set_atoms(sample_interface_atoms, ref_interface_atoms)
    #super_imposer.apply(sample_model_backbone.get_atoms())

    irms = super_imposer.rms
    
    ref_group1_size = np.sum([len(ref_model[chain]) for chain in nat_group1])
    ref_group2_size = np.sum([len(ref_model[chain]) for chain in nat_group2])

    receptor_chains = (nat_group1, group1) if ref_group1_size > ref_group2_size else (nat_group2, group2)
    ligand_chains = (nat_group1, group1) if ref_group1_size <= ref_group2_size else (nat_group2, group2)

    class1, class2 = ("receptor", "ligand") if ref_group1_size > ref_group2_size else ("ligand", "receptor")
        
    receptor_atoms_native = [atom for chain in receptor_chains[0] for atom in ref_model_backbone[chain].get_atoms()]
    receptor_atoms_sample = [atom for chain in receptor_chains[1] for atom in sample_model_backbone[chain].get_atoms()]
    # Set to align on receptor
    super_imposer.set_atoms(receptor_atoms_native, receptor_atoms_sample)
    super_imposer.apply(sample_model_backbone.get_atoms())
    receptor_chain_rms = super_imposer.rms

    coord1 = np.array([atom.coord for chain in ligand_chains[0] for atom in ref_model_backbone[chain].get_atoms()])
    coord2 = np.array([atom.coord for chain in ligand_chains[1] for atom in sample_model_backbone[chain].get_atoms()])

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


def align_sequence_features(model_sequence, native_sequence, native_numbering):
    """Realigns sequence-specific features between PDB structures and antigen sequences
    For example, if from_sequence (the sequence for which features have
    been calculated) is missing an amino acid but contains more amino acids
    at the C terminal than the to_sequence, the features will be assigned
    as in this schematic:
    
    native_numbering            23456789
    native_sequence_aligned    -WEKLAPTG
    model_sequence_aligned     DWEKLAPT-
    model_numbering            12345678-1
    Wherever it is not possible to assign features (missing AA in native_sequence)
    a pseudovalue of -1 will be assigned instead
    """
    alignment = pairwise2.align.localms(
        model_sequence, native_sequence, match=5, mismatch=0, open=-10, extend=-1
    )[0]
    model_sequence_aligned = np.array(list(alignment.seqA))
    native_sequence_aligned = np.array(list(alignment.seqB))

    model_indexes = np.where(model_sequence_aligned != "-")[0]
    native_indexes = np.where(native_sequence_aligned != "-")[0]
    alignment_length = np.sum(
        (model_sequence_aligned != "-") & (native_sequence_aligned != "-")
    )

    aligned_numbering = np.zeros(len(model_indexes), dtype=int) - 1

    aligned_numbering[native_indexes[:alignment_length]] = native_numbering[
        model_indexes[:alignment_length]
    ]

    return aligned_numbering

def align_model_to_native(model_structure, native_structure, model_chain, native_chain):
    model_sequence = "".join(
        seq1(residue.get_resname()) for residue in model_structure[model_chain].get_residues()
    )
    
    native_sequence = "".join(
        seq1(residue.get_resname()) for residue in native_structure[native_chain].get_residues()
    )

    alignment = pairwise2.align.localms(
        model_sequence, native_sequence, match=5, mismatch=0, open=-10, extend=-1
    )[0]

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
    seqA = alignment.seqA if not invert else alignment.seqB
    seqB = alignment.seqB if not invert else alignment.seqA
    for (aligned_residue_A, aligned_residue_B) in zip(seqA, seqB):
        if aligned_residue_A != "-":
            residue = next(residues)
        if aligned_residue_B == "-": # gap residue: remove from structure
            residues_to_delete.append(residue.get_full_id())

    for _, _, _, res in residues_to_delete:
        model[chain].detach_child(res)


def get_distances_across_chains(model, group1, group2, all_atom=True):
    
    if all_atom:
        model_A_atoms = np.asarray([atom.get_coord() for chain in group1 for res in model[chain].get_residues() for atom in res.get_atoms()])
        model_B_atoms = np.asarray([atom.get_coord() for chain in group2 for res in model[chain].get_residues() for atom in res.get_atoms()])
    else:
        model_A_atoms = np.asarray([res["CB"].get_coord() if "CB" in res else res["CA"].get_coord() for chain in group1 for res in model[chain].get_residues()])
        model_B_atoms = np.asarray([res["CB"].get_coord() if "CB" in res else res["CA"].get_coord() for chain in group2 for res in model[chain].get_residues()])

    distances = np.sqrt(((model_A_atoms[:, None] - model_B_atoms[None, :]) ** 2).sum(-1))
    return distances

#@profile
def group_atom_into_res_distances(atom_distances, group_dic1, group_dic2):
    res_distances = np.zeros((len(group_dic1), len(group_dic2)))
    
    cum_i_atoms = 0
    for i, i_atoms in enumerate(group_dic1):
        cum_j_atoms = 0
        for j, j_atoms in enumerate(group_dic2):
            res_distances[i, j] = atom_distances[cum_i_atoms:cum_i_atoms+i_atoms, cum_j_atoms:cum_j_atoms+j_atoms].min()
            cum_j_atoms += j_atoms    
        cum_i_atoms += i_atoms
    return res_distances


def get_group_dictionary(model, group):
    n_atoms_per_residue = []
    for chain in group:
        for residue in model[chain].get_residues():
            n_atoms_per_residue.append(len(residue.get_unpacked_list()))
    return n_atoms_per_residue


#@profile
def get_residue_distances(structure, group1, group2):
    # get information about how many atoms correspond to each amino acid in each group of chains
    model_group_dic1 = get_group_dictionary(structure, group1)
    model_group_dic2 = get_group_dictionary(structure, group2)
    
    model_atom_distances = get_distances_across_chains(structure, group1, group2)
    
    model_res_distances = group_atom_into_res_distances(model_atom_distances, model_group_dic1, model_group_dic2)
    return model_res_distances


def get_fnat(model_res_distances, native_res_distances, thr=5.0):

    native_contacts = native_res_distances < thr
    model_contacts = model_res_distances < thr
    n_native_contacts = np.sum(native_contacts)
    n_model_contacts = np.sum(model_contacts)
    n_shared_contacts = np.sum(model_contacts * native_contacts)
    n_non_native_contacts = np.sum(model_contacts * (1 - native_contacts))
    return (n_shared_contacts, n_non_native_contacts, n_native_contacts, n_model_contacts)


def get_interacting_pairs(distances, thr=0.5):
    return np.nonzero(distances < thr)


def get_interface_atoms(interacting_pairs, model_backbone, ref_backbone, model_group1, model_group2, ref_group1, ref_group2):
    ref_interface = []
    mod_interface = []

    ref_residues_group1 = [res for chain in ref_group1 for res in ref_backbone[chain].get_residues()]
    ref_residues_group2 = [res for chain in ref_group2 for res in ref_backbone[chain].get_residues()]
    
    mod_residues_group1 = [res for chain in model_group1 for res in model_backbone[chain].get_residues()]
    mod_residues_group2 = [res for chain in model_group2 for res in model_backbone[chain].get_residues()]
    
    # get the native interfacial residues, along with the corresponding model residues, and select common backbone atoms
    for i, j in zip(interacting_pairs[0], interacting_pairs[1]):
        ref_res1_atoms = ref_residues_group1[i].get_atoms()
        mod_res1_atoms = mod_residues_group1[i].get_atoms()
        
        ref_res2_atoms = ref_residues_group2[j].get_atoms()
        mod_res2_atoms = mod_residues_group2[j].get_atoms()
        
        for atom1, atom2 in zip(ref_res1_atoms, mod_res1_atoms):
            if atom1 not in ref_interface:
                ref_interface.append(atom1)
            if atom2 not in mod_interface:
                mod_interface.append(atom2)
            
        for atom1, atom2 in zip(ref_res2_atoms, mod_res2_atoms):
            if atom1 not in ref_interface:
                ref_interface.append(atom1)
            if atom2 not in mod_interface:
                mod_interface.append(atom2)
                
    return list(mod_interface), list(ref_interface)


def set_common_backbone_atoms(model, reference, atom_types=["CA", "C", "N", "O"]):
    # model and reference should have the same number of amino acids and be aligned
    for mod_res, ref_res in zip(model.get_residues(), reference.get_residues()):
        mod_atoms = [atom for atom in mod_res.get_atoms()]
        ref_atoms = [atom for atom in ref_res.get_atoms()]
        
        atom_ids_in_mod_res = [atm.id for atm in mod_atoms]
        atom_ids_in_ref_res = [atm.id for atm in ref_atoms]
        
        atom_ids_in_ref_and_mod_res = set(atom_ids_in_mod_res).intersection(atom_types).intersection(atom_ids_in_ref_res)
        
        # whatever atom is not in the shared list, remove it from the both structures
        for atom_id in set(atom_ids_in_mod_res).difference(atom_ids_in_ref_and_mod_res):
            mod_res.detach_child(atom_id)
        
        for atom_id in set(atom_ids_in_ref_res).difference(atom_ids_in_ref_and_mod_res):
            ref_res.detach_child(atom_id)

#@profile
def main():
    args = parse_args()

    bio_ver = 1.61
    if float(Bio.__version__) < bio_ver:
        print(
            "Biopython version (%s) is too old need at least >=%.2f"
            % (Bio.__version__, bio_ver)
        )
        sys.exit()

    exec_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    fix_numbering = exec_path + "/scripts/fix_numbering.pl"
    model = args.model
    native = args.native
    use_CA_only = args.useCA
    capri_peptide = args.capri_peptide

    # Start the parser
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)

    # Get the structures
    native_structure = pdb_parser.get_structure("reference", native)[0]
    model_structure = pdb_parser.get_structure("model", model)[0]

    model_chains = []
    native_chains = []
    best_info = ""
    if not args.skip_check:
        native_chains = [c.id for c in native_structure]
        model_chains = [c.id for c in model_structure]

    files_to_clean = []

    if (len(model_chains) > 2 or len(native_chains) > 2) and (
        args.model_chain1 == None and args.native_chain1 == None
    ):
        print(
            "Multi-chain model need sets of chains to group\nuse -native_chain1 and/or -model_chain1 if you want a different mapping than 1-1"
        )
        print("Model chains  : " + str(model_chains))
        print("Native chains : " + str(native_chains))
        sys.exit()
    if not args.skip_check and (len(model_chains) < 2 or len(native_chains) < 2):
        print("Need at least two chains in the two inputs\n")
        sys.exit()
    
    group1 = model_chains[0]
    group2 = model_chains[1]
    
    nat_group1 = native_chains[0]
    nat_group2 = native_chains[1]
    if len(model_chains) > 2 or len(native_chains) > 2:
        if args.model_chain1 != None:
            group1 = args.model_chain1
            nat_group1 = group1
            if args.model_chain2 != None:
                group2 = args.model_chain2
            else:
                # will use the complement from group1
                group2 = []
                for c in model_chains:
                    if c not in group1:
                        group2.append(c)
            nat_group1 = group1
            nat_group2 = group2

        if args.native_chain1 != None:
            nat_group1 = args.native_chain1
            if args.native_chain2 != None:
                nat_group2 = args.native_chain2
            else:
                # will use the complement from group1
                nat_group2 = []
                for c in native_chains:
                    if c not in nat_group1:
                        nat_group2.append(c)

        if args.model_chain1 == None:
            group1 = nat_group1
            group2 = nat_group2

        pe = 0
        if args.perm1 or args.perm2:
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
                    "Starting chain order permutation search (number of permutations: "
                    + str(pe_tot)
                    + ")"
                )
            for g1 in combos1:
                for g2 in combos2:
                    model_renum = make_two_chain_pdb_perm(model, g1, g2)
                    model_fixed = model_renum
                    if not args.no_needle:

                        fix_numbering_cmd = (
                            fix_numbering
                            + " "
                            + model_renum
                            + " "
                            + native
                            + " > /dev/null"
                        )
                        model_fixed = model_renum + ".fixed"
                        os.system(fix_numbering_cmd)
                        os.remove(model_renum)
                        if not os.path.exists(model_fixed):
                            print(
                                "If you are sure the residues are identical you can use the options -no_needle"
                            )
                            sys.exit()
                    test_dict = calc_DockQ(model_structure, native_structure,  group1, group2, nat_group1, nat_group2, use_CA_only)
                    os.remove(model_fixed)
                    if not args.quiet:
                        print(
                            str(pe)
                            + "/"
                            + str(pe_tot)
                            + " "
                            + "".join(g1)
                            + " -> "
                            + "".join(g2)
                            + " "
                            + str(test_dict["DockQ"])
                        )
                    if test_dict["DockQ"] > best_DockQ:
                        best_DockQ = test_dict["DockQ"]
                        info = test_dict
                        best_g1 = g1
                        best_g2 = g2
                        best_info = (
                            "Best score ( "
                            + str(best_DockQ)
                            + " ) found for model -> native, chain1:"
                            + "".join(best_g1)
                            + " -> "
                            + "".join(nat_group1)
                            + " chain2:"
                            + "".join(best_g2)
                            + " -> "
                            + "".join(nat_group2)
                        )

                        if args.verbose:
                            print(best_info)
                        if not args.quiet:
                            print("Current best " + str(best_DockQ))
                    pe = pe + 1
            if not args.quiet:
                print(best_info)
        else:

            model_structure = remove_extra_chains(model_structure, chains_to_keep=group1 + group2)
            native_structure = remove_extra_chains(native_structure, chains_to_keep=nat_group1 + nat_group2)
        
            # realign each model chain against the corresponding native chain
            for model_chain, native_chain in zip(group1 + group2, nat_group1 + nat_group2):
                alignment = align_model_to_native(model_structure, native_structure, model_chain, native_chain)
                fix_chain_residues(model_structure, model_chain, alignment)
                fix_chain_residues(native_structure, native_chain, alignment, invert=True)
            info = calc_DockQ(model_structure, native_structure,  group1, group2, nat_group1, nat_group2, use_CA_only=use_CA_only, capri_peptide=capri_peptide)

    else:
        model_structure = remove_extra_chains(model_structure, chains_to_keep=group1 + group2)
        native_structure = remove_extra_chains(native_structure, chains_to_keep=nat_group1 + nat_group2)
    
        # realign each model chain against the corresponding native chain
        for model_chain, native_chain in zip(group1 + group2, nat_group1 + nat_group2):
            alignment = align_model_to_native(model_structure, native_structure, model_chain, native_chain)
            fix_chain_residues(model_structure, model_chain, alignment)
            fix_chain_residues(native_structure, native_chain, alignment, invert=True)
        info = calc_DockQ(model_structure, native_structure,  group1, group2, nat_group1, nat_group2, use_CA_only=use_CA_only, capri_peptide=capri_peptide)

    irms = info["irms"]
    Lrms = info["Lrms"]
    fnat = info["fnat"]
    DockQ = info["DockQ"]
    fnonnat = info["fnonnat"]

    if args.short:
        if capri_peptide:
            print(
                (
                    "DockQ-capri_peptide %.3f Fnat %.3f iRMS %.3f LRMS %.3f Fnonnat %.3f %s %s %s"
                    % (DockQ, fnat, irms, Lrms, fnonnat, model, native_in, best_info)
                )
            )
        else:
            print(
                (
                    "DockQ %.3f Fnat %.3f iRMS %.3f LRMS %.3f Fnonnat %.3f %s %s %s"
                    % (DockQ, fnat, irms, Lrms, fnonnat, model, native, best_info)
                )
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
            print("*   Reference: Sankar Basu and Bjorn Wallner, DockQ: A quality *")
            print("*   measure for protein-protein docking models, submitted      *")
            print("*                                                              *")
            print("*   For the record:                                            *")
            print("*   Definition of contact <5A (Fnat)                           *")
            print("*   Definition of interface <10A all heavy atoms (iRMS)        *")
            print("*   For comments, please email: bjorn.wallner@.liu.se          *")
            print("*                                                              *")
            print("****************************************************************")
        print(("Model  : %s" % model))
        print(("Native : %s" % native))
        if len(best_info):
            print(best_info)
        print(
            "Number of equivalent residues in chain "
            + info["chain1"]
            + " "
            + str(info["len1"])
            + " ("
            + info["class1"]
            + ")"
        )
        print(
            "Number of equivalent residues in chain "
            + info["chain2"]
            + " "
            + str(info["len2"])
            + " ("
            + info["class2"]
            + ")"
        )
        print(
            (
                "Fnat %.3f %d correct of %d native contacts"
                % (info["fnat"], info["nat_correct"], info["nat_total"])
            )
        )
        print(
            (
                "Fnonnat %.3f %d non-native of %d model contacts"
                % (info["fnonnat"], info["nonnat_count"], info["model_total"])
            )
        )
        print(("iRMS %.3f" % irms))
        print(("LRMS %.3f" % Lrms))
        peptide_suffix = ""
        if capri_peptide:
            peptide_suffix = "_peptide"

        peptide_disclaimer = ""
        if capri_peptide:
            peptide_disclaimer = "DockQ not reoptimized for CAPRI peptide evaluation"
        print(("DockQ {:.3f} {}".format(DockQ, peptide_disclaimer)))

    for f in files_to_clean:
        os.remove(f)


if __name__ == "__main__":
    main()
