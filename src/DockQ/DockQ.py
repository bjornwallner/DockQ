#!/usr/bin/env python

import sys
import gzip
import hashlib
import traceback
import itertools
from functools import lru_cache, wraps
from argparse import ArgumentParser

import Bio.PDB
import numpy as np
from Bio import Align
from Bio.SeqUtils import seq1
from Bio.SVDSuperimposer import SVDSuperimposer
from .parsers import PDBParser

# fallback in case the cython version doesn't work, though it will be slower
try:
    from .operations import residue_distances, get_fnat_stats
except ImportError:
    print(
        "WARNING: It looks like cython is not working, \
        falling back on native python. This will make DockQ slower"
    )
    from operations_nocy import residue_distances, get_fnat_stats


def parse_args():
    parser = ArgumentParser(
        description="DockQ - Quality measure for \
        protein-protein docking models"
    )
    parser.add_argument("model", metavar="<model>", type=str, help="path to model file")
    parser.add_argument(
        "native", metavar="<native>", type=str, help="path to native file"
    )
    parser.add_argument(
        "--capri_peptide",
        default=False,
        action="store_true",
        help="use version for capri_peptide \
        (DockQ cannot not be trusted for this setting)",
    )
    parser.add_argument(
        "--short", default=False, action="store_true", help="short output"
    )
    parser.add_argument(
        "--verbose", "-v", default=False, action="store_true", help="talk a lot!"
    )
    parser.add_argument(
        "--use_CA", "-ca", default=False, action="store_true", help="use CA instead of backbone"
    )
    parser.add_argument(
        "--no_needle",
        default=False,
        action="store_true",
        help="Do not align native and model using sequence alignments, but use the numbering of residues instead",
    )
    parser.add_argument(
        "--optDockQF1",
        default=False,
        action="store_true",
        help="optimize on DockQ_F1 instead of DockQ",
    )
    parser.add_argument(
        "--allowed_mismatches",
        default=0,
        type=int,
        help="number of allowed mismatches when mapping model sequence to native sequence.",
    )
    parser.add_argument(
        "--mapping",
        default=None,
        metavar="MODELCHAINS:NATIVECHAINS",
        help="""Specify a chain mapping between model and native structure.
            If the native contains two chains "H" and "L"
            while the model contains two chains "A" and "B",
            and chain A is a model of native chain
            H and chain B is a model of native chain L,
            the flag can be set as: '--mapping AB:HL'.
            This can also help limit the search to specific native interfaces.
            For example, if the native is a tetramer (ABCD) but the user is only interested
            in the interface between chains B and C, the flag can be set as: '--mapping :BC'
            or the equivalent '--mapping *:BC'."""
             )

    return parser.parse_args()


def open_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename) as f:
            yield f
    else:
        with open(filename) as f:
            yield f


@lru_cache
def get_aligned_residues(chainA, chainB, alignment):
    aligned_resA = []
    aligned_resB = []
    resA = chainA.get_residues()
    resB = chainB.get_residues()

    if alignment[0] == alignment[2]:
        return tuple(resA), tuple(resB)

    for A, match, B in zip(*alignment):
        if A != "-":
            rA = next(resA)
        if B != "-":
            rB = next(resB)

        if match == "|":
            aligned_resA.append(rA)
            aligned_resB.append(rB)

    return tuple(aligned_resA), tuple(aligned_resB)


@lru_cache
def list_atoms_per_residue2(chain, what):
    n_atoms_per_residue = []
    atoms_ids = []
    for residue in chain:
        # important to remove duplicate atoms (e.g. alternates) at this stage
        atom_ids = tuple([a.id for a in residue.get_unpacked_list()])
        atoms_ids.append(atom_ids)
        n_atoms_per_residue.append(len(set(atom_ids)))
    return tuple(n_atoms_per_residue), tuple(atoms_ids)


@lru_cache
def get_atoms(chain, what):
    atoms = np.array(
            [
                atom.get_coord()
                for res in chain
                for atom in res.get_atoms()
            ]
        )
    return atoms


@lru_cache
def get_residue_distances(chain1, chain2, what, all_atom=True):
    if all_atom:
        # how many atoms per aligned amino acid
        n_atoms_per_res_chain1 = list_atoms_per_residue(chain1, what)
        n_atoms_per_res_chain2 = list_atoms_per_residue(chain2, what)
        model_A_atoms = np.asarray(
            [
                atom.get_coord()
                for res in chain1
                for atom in res.get_atoms()
            ]
        )
        model_B_atoms = np.asarray(
            [
                atom.get_coord()
                for res in chain2
                for atom in res.get_atoms()
            ]
        )

    else:  # distances were already between CBs only
        model_A_atoms = np.asarray(
            [
                res["CB"].get_coord() if "CB" in res else res["CA"].get_coord()
                for res in chain1
            ]
        )
        model_B_atoms = np.asarray(
            [
                res["CB"].get_coord() if "CB" in res else res["CA"].get_coord()
                for res in chain2
            ]
        )

        n_atoms_per_res_chain1 = np.ones(model_A_atoms.shape[0]).astype(int)
        n_atoms_per_res_chain2 = np.ones(model_B_atoms.shape[0]).astype(int)

    model_res_distances = residue_distances(
        model_A_atoms, model_B_atoms, n_atoms_per_res_chain1, n_atoms_per_res_chain2
    )
    return model_res_distances


@lru_cache
def get_aligned_atoms(n_atoms_per_res_chainA, n_atoms_per_res_chainB, alignment):

    if alignment[0] == alignment[2]: # no need to align
        return np.ones(sum(list(n_atoms_per_res_chainA))).astype(bool), np.ones(sum(list(n_atoms_per_res_chainB))).astype(bool), [i for i in range(len(alignment[0]))], [i for i in range(len(alignment[1]))]
    atom_index_A = np.zeros(sum(list(n_atoms_per_res_chainA))).astype(bool)
    atom_index_B = np.zeros(sum(list(n_atoms_per_res_chainB))).astype(bool)
    res_index_A = []
    res_index_B = []

    n_atoms_per_res_chainA = iter(n_atoms_per_res_chainA)
    n_atoms_per_res_chainB = iter(n_atoms_per_res_chainB)
    countA = iA = 0
    countB = iB = 0

    for A, match, B in zip(*alignment):
        if A != "-":
            nA = next(n_atoms_per_res_chainA)
            countA += nA
            iA += 1
        if B != "-":
            nB = next(n_atoms_per_res_chainB)
            countB += nB
            iB += 1

        if match == "|":
            atom_index_A[countA-nA:countA] = 1
            atom_index_B[countB-nB:countB] = 1
            res_index_A.append(iA-1)
            res_index_B.append(iB-1)

    return atom_index_A, atom_index_B, res_index_A, res_index_B


class NpWrapper:
    def __init__(self, x: np.array) -> None:
        self.values = x
        # here you can use your own hashing function
        self.h = hashlib.sha224(x.view()).hexdigest()

    def __hash__(self) -> int:
        return hash(self.h)

    def __eq__(self, __value: object) -> bool:
        return __value.h == self.h


def dist_memoizer(get_atom_distances):
    @lru_cache()
    def cached_wrapper(np1, np2, n_atoms, what):
        return get_atom_distances(np1.values, np2.values, n_atoms, what)

    @wraps(get_atom_distances)
    def wrapper(x: np.array, y: np.array, n_atoms, what):
        np1 = NpWrapper(x)
        np2 = NpWrapper(y)
        return cached_wrapper(np1, np2, n_atoms, what)

    return wrapper


@dist_memoizer
def get_atom_distances(chain1_atoms, chain2_atoms, n_atoms_per_res, what):
    n_atoms_per_res_chain1, n_atoms_per_res_chain2 = n_atoms_per_res

    model_res_distances = residue_distances(
        chain1_atoms, chain2_atoms, np.array(n_atoms_per_res_chain1).astype(int), np.array(n_atoms_per_res_chain2).astype(int)
    )
    return model_res_distances


def sup_memoizer(superimpose):
    @lru_cache()
    def cached_wrapper(np1, np2):
        return superimpose(np1.values, np2.values)

    @wraps(superimpose)
    def wrapper(x: np.array, y: np.array):
        np1 = NpWrapper(x)
        np2 = NpWrapper(y)
        return cached_wrapper(np1, np2)

    return wrapper


@sup_memoizer
def superimpose(from_atoms, to_atoms):
    super_imposer = SVDSuperimposer()
    # Set to align on receptor
    super_imposer.set(to_atoms, from_atoms)
    super_imposer.run()
    rot, tran = super_imposer.get_rotran()
    return rot, tran, super_imposer


# @profile
def calc_DockQ2(
    model_chains,
    native_chains,
    alignments,
    use_CA_only=False,
    capri_peptide=False,
):
    atom_for_sup = ("CA", "C", "N", "O") if not use_CA_only else ("CA")
    fnat_threshold = 4.0 if capri_peptide else 5.0
    interface_threshold = 8.0 if capri_peptide else 10.0

    native_atoms_1 = get_atoms(native_chains[0], "native")
    native_atoms_2 = get_atoms(native_chains[1], "native")
    native_n_atoms_1, native_atom_ids_1 = list_atoms_per_residue2(native_chains[0], "native")
    native_n_atoms_2, native_atom_ids_2 = list_atoms_per_residue2(native_chains[1], "native")

    native_distances = get_atom_distances(native_atoms_1, native_atoms_2, (native_n_atoms_1, native_n_atoms_2), "native")
    nat_total = np.nonzero(np.asarray(native_distances) < fnat_threshold**2)[
        0
    ].shape[0]

    if nat_total == 0:
        return None

    model_atoms_1 = get_atoms(model_chains[0], "model")
    model_atoms_2 = get_atoms(model_chains[1], "model")
    model_n_atoms_1, model_atom_ids_1 = list_atoms_per_residue2(model_chains[0], "model")
    model_n_atoms_2, model_atom_ids_2 = list_atoms_per_residue2(model_chains[1], "model")

    model_atom_index_1, native_atom_index_1, model_res_index_1, native_res_index_1 = get_aligned_atoms(
        model_n_atoms_1, native_n_atoms_1, tuple(alignments[0])
    )
    model_atom_index_2, native_atom_index_2, model_res_index_2, native_res_index_2 = get_aligned_atoms(
        model_n_atoms_2, native_n_atoms_2, tuple(alignments[1])
    )

    if len(model_res_index_1) != len(model_chains[0]):
        # filter model, native atoms by alignment atom indexes
        model_atoms_1 = model_atoms_1[model_atom_index_1]
        model_n_atoms_1 = tuple([n_atoms for i, n_atoms in enumerate(model_n_atoms_1) if i in model_res_index_1])
        model_atom_ids_1 = tuple([atom_ids for i, atom_ids in enumerate(model_atom_ids_1) if i in model_res_index_1])

    if len(model_res_index_2) != len(model_chains[1]):
        model_atoms_2 = model_atoms_2[model_atom_index_2]
        model_n_atoms_2 = tuple([n_atoms for i, n_atoms in enumerate(model_n_atoms_2) if i in model_res_index_2])
        model_atom_ids_2 = tuple([atom_ids for i, atom_ids in enumerate(model_atom_ids_2) if i in model_res_index_2])

    model_distances = get_atom_distances(model_atoms_1, model_atoms_2, (model_n_atoms_1, model_n_atoms_2), "model")
    # only recalculate distances if model was missing residues compared to native
    if len(native_res_index_1) != len(native_chains[0]) or len(native_res_index_2) != len(native_chains[1]):
        native_atoms_1 = native_atoms_1[native_atom_index_1]
        native_atoms_2 = native_atoms_1[native_atom_index_2]
        native_n_atoms_1 = tuple([n_atoms for i, n_atoms in enumerate(native_n_atoms_1) if i in native_res_index_1])
        native_n_atoms_2 = tuple([n_atoms for i, n_atoms in enumerate(native_n_atoms_2) if i in native_res_index_2])
        native_atom_ids_1 = tuple([atom_ids for i, atom_ids in enumerate(native_atom_ids_1) if i in native_res_index_1])
        native_atom_ids_2 = tuple([atom_ids for i, atom_ids in enumerate(native_atom_ids_2) if i in native_res_index_2])
        native_distances = get_atom_distances(native_atoms_1, native_atoms_2, (native_n_atoms_1, native_n_atoms_2), "native")

    assert (
        model_distances.shape == native_distances.shape
    ), f"Native and model have incompatible sizes ({model_distances.shape} != {native_distances.shape})"

    nat_correct, nonnat_count, _, model_total = get_fnat_stats(
        model_distances, native_distances, threshold=fnat_threshold
    )

    # avoids divide by 0 errors
    fnat = nat_total and nat_correct / nat_total or 0
    fnonnat = model_total and nonnat_count / model_total or 0

    # filter out backbone atoms now
    model_index_1, native_index_1, n_atoms_1 = filter_atoms(model_atom_ids_1, native_atom_ids_1, filter=atom_for_sup)
    model_index_2, native_index_2, n_atoms_2 = filter_atoms(model_atom_ids_2, native_atom_ids_2, filter=atom_for_sup)
    model_atoms_1 = model_atoms_1[model_index_1]
    model_atoms_2 = model_atoms_2[model_index_2]
    native_atoms_1 = native_atoms_1[native_index_1]
    native_atoms_2 = native_atoms_2[native_index_2]

    # assign which chains constitute the receptor, ligand
    native1_size = len(native_chains[0])
    native2_size = len(native_chains[1])
    receptor_atoms = (
        (model_atoms_1, native_atoms_1)
        if native1_size > native2_size
        else (model_atoms_2, native_atoms_2)
    )
    ligand_atoms = (
        (model_atoms_1, native_atoms_1)
        if native1_size <= native2_size
        else (model_atoms_2, native_atoms_2)
    )
    receptor_n_atoms = (
        n_atoms_1
        if native1_size > native2_size
        else n_atoms_2
    )
    ligand_n_atoms = (
        n_atoms_1
        if native1_size <= native2_size
        else n_atoms_2
    )

    class1, class2 = (
        ("receptor", "ligand")
        if native1_size > native2_size
        else ("ligand", "receptor")
    )

    # Set to align on receptor
    rot, tran, super_imposer = superimpose(receptor_atoms[0], receptor_atoms[1])
    rotated_ligand_atoms = np.dot(ligand_atoms[0], rot) + tran

    Lrms = super_imposer._rms(
        ligand_atoms[1], rotated_ligand_atoms
    )  # using the private _rms function which does not superimpose

    # Get interfacial atoms from reference, and corresponding atoms from sample
    interacting_pairs = get_interacting_pairs(
        # working with squared thresholds to avoid using sqrt
        native_distances,
        threshold=interface_threshold**2,
    )

    if class1 == "ligand":
        interacting_pairs = interacting_pairs[::-1]

    receptor_interface_index, _, _ = filter_atoms(receptor_n_atoms, None, filter=interacting_pairs[0])
    ligand_interface_index, _, _ = filter_atoms(ligand_n_atoms, None, filter=interacting_pairs[1])

    model_interface_atoms = np.vstack((receptor_atoms[0][receptor_interface_index], ligand_atoms[0][ligand_interface_index]))
    native_interface_atoms = np.vstack((receptor_atoms[1][receptor_interface_index], ligand_atoms[1][ligand_interface_index]))

    _, _, super_imposer = superimpose(model_interface_atoms, native_interface_atoms)
    irms = super_imposer.get_rms()

    info = {}

    info["DockQ_F1"] = dockq_formula(f1(nat_correct, nonnat_count, nat_total), irms, Lrms)
    info["DockQ"] = dockq_formula(fnat, irms, Lrms)
    info["irms"] = irms
    info["Lrms"] = Lrms
    info["fnat"] = fnat
    info["nat_correct"] = nat_correct
    info["nat_total"] = nat_total

    info["fnonnat"] = fnonnat
    info["nonnat_count"] = nonnat_count
    info["model_total"] = model_total
    
    info["len1"] = native1_size
    info["len2"] = native2_size
    info["class1"] = class1
    info["class2"] = class2
    return info


@lru_cache
def filter_atoms(model_info, native_info, filter=("CA", "C", "N", "O", "P")):
    model_bools = []
    native_bools = []
    if "CA" in filter:  # backbone filter
        model_atom_ids = model_info
        native_atom_ids = native_info

        n_atoms_per_res = []
        for model_ids, native_ids in zip(model_atom_ids, native_atom_ids):
            model_select = [1 if atom in filter and atom in native_ids else 0 for atom in list(dict.fromkeys(model_ids))]
            native_select = [1 if atom in filter and atom in model_ids else 0 for atom in list(dict.fromkeys(native_ids))]
            model_bools.extend(model_select)
            native_bools.extend(native_select)
            n_atoms_per_res.append(sum(native_select))
            if sum(model_select) != sum(native_select):
                print(model_ids, model_select)
                print(native_ids, native_select)
                print("ERROR!")
    else:  # filter by index
        n_atoms_per_res = model_info

        for i, n_atoms in enumerate(n_atoms_per_res):
            if i in filter:
                model_bools.extend([1 for i in range(n_atoms)])
            else:
                model_bools.extend([0 for i in range(n_atoms)])
        native_bools = model_bools
    model_bools = np.array(model_bools).astype(bool)
    native_bools = np.array(native_bools).astype(bool)
    return model_bools, native_bools, tuple(n_atoms_per_res)


# @profile
def calc_DockQ(
    sample_chains,
    ref_chains,
    alignments,
    use_CA_only=False,
    capri_peptide=False,
):
    atom_for_sup = ("CA", "C", "N", "O", "P") if not use_CA_only else ("CA", "P")
    fnat_threshold = 4.0 if capri_peptide else 5.0
    interface_threshold = 8.0 if capri_peptide else 10.0
    clash_threshold=2.0
    # total number of native contacts is calculated on untouched native structure
    ref_res_distances = get_residue_distances(ref_chains[0], ref_chains[1], "ref")
    nat_total = np.nonzero(np.asarray(ref_res_distances) < fnat_threshold**2)[
        0
    ].shape[0]

    if nat_total == 0:
        # if the native has no interface between the two chain groups
        # nothing to do here
        return None

    aligned_sample_1, aligned_ref_1 = get_aligned_residues(
        sample_chains[0], ref_chains[0], alignments[0]
    )
    aligned_sample_2, aligned_ref_2 = get_aligned_residues(
        sample_chains[1], ref_chains[1], alignments[1]
    )

    sample_res_distances = get_residue_distances(aligned_sample_1, aligned_sample_2, "sample")
    ref_res_distances = get_residue_distances(aligned_ref_1, aligned_ref_2, "ref")

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
            ref_chains[0], ref_chains[1], "ref", all_atom=False
        )
    # Get interfacial atoms from reference, and corresponding atoms from sample
    interacting_pairs = get_interacting_pairs(
        # working with squared thresholds to avoid using sqrt
        ref_res_distances,
        threshold=interface_threshold**2,
    )

    # get a copy of each structure, then only keep backbone atoms
    sample_interface_atoms, ref_interface_atoms = get_interface_atoms(
        interacting_pairs,
        (aligned_sample_1, aligned_sample_2),
        (aligned_ref_1, aligned_ref_2),
        atom_types=atom_for_sup,
    )
    super_imposer = SVDSuperimposer()
    super_imposer.set(sample_interface_atoms, ref_interface_atoms)
    super_imposer.run()
    irms = super_imposer.get_rms()

    # assign which chains constitute the receptor, ligand
    ref_group1_size = len(ref_chains[0])
    ref_group2_size = len(ref_chains[1])
    receptor_chains = (
        (aligned_ref_1, aligned_sample_1)
        if ref_group1_size > ref_group2_size
        else (aligned_ref_2, aligned_sample_2)
    )
    ligand_chains = (
        (aligned_ref_1, aligned_sample_1)
        if ref_group1_size <= ref_group2_size
        else (aligned_ref_2, aligned_sample_2)
    )
    class1, class2 = (
        ("receptor", "ligand")
        if ref_group1_size > ref_group2_size
        else ("ligand", "receptor")
    )

    receptor_atoms_native, receptor_atoms_sample = np.asarray(
        get_atoms_per_residue(receptor_chains, what="receptor", atom_types=atom_for_sup)
    )
    ligand_atoms_native, ligand_atoms_sample = np.asarray(
        get_atoms_per_residue(ligand_chains, what="ligand", atom_types=atom_for_sup)
    )
    # Set to align on receptor
    super_imposer.set(receptor_atoms_native, receptor_atoms_sample)
    super_imposer.run()

    rot, tran = super_imposer.get_rotran()
    rotated_sample_atoms = np.dot(ligand_atoms_sample, rot) + tran

    Lrms = super_imposer._rms(
        ligand_atoms_native, rotated_sample_atoms
    )  # using the private _rms function which does not superimpose

    info = {}

    info["DockQ_F1"] = dockq_formula(f1(nat_correct, nonnat_count, nat_total), irms, Lrms)
    info["DockQ"] = dockq_formula(fnat, irms, Lrms)
    info["irms"] = irms
    info["Lrms"] = Lrms
    info["fnat"] = fnat
    info["nat_correct"] = nat_correct
    info["nat_total"] = nat_total

    info["fnonnat"] = fnonnat
    info["nonnat_count"] = nonnat_count
    info["model_total"] = model_total
    info["clashes"]=np.nonzero(np.asarray(sample_res_distances) < clash_threshold**2)[0].shape[0]
    info["len1"] = ref_group1_size
    info["len2"] = ref_group2_size
    info["class1"] = class1
    info["class2"] = class2

    return info


def f1(tp, fp, p):
    return 2 * tp / (tp + fp + p)


def dockq_formula(fnat, irms, Lrms):
    return (
        float(fnat)
        + 1 / (1 + (irms / 1.5) * (irms / 1.5))
        + 1 / (1 + (Lrms / 8.5) * (Lrms / 8.5))
    ) / 3


@lru_cache
def align_chains(model_chain, native_chain, use_numbering=False):
    """
    Function to align two PDB structures. This can be done by sequence (default) or by
    numbering. If the numbering is used, then each residue number from the pdb structure
    is converted to a unique character. Then the two vectors of character are aligned
    as if they were two sequences
    """

    if use_numbering:
        model_numbering = []
        native_numbering = []

        for residue in model_chain.get_residues():
            resn = int(residue.id[1])
            model_numbering.append(resn)

        for residue in native_chain.get_residues():
            resn = int(residue.id[1])
            native_numbering.append(resn)
        # if the samllest resn is negative, it will be used to shift all numbers so they start from 0
        # the minimum offset is 45 to avoid including the "-" character that is reserved for gaps
        min_resn = max(45, -min(model_numbering + native_numbering))

        model_sequence = "".join([chr(resn + min_resn) for resn in model_numbering])
        native_sequence = "".join([chr(resn + min_resn) for resn in native_numbering])

    else:
        model_sequence = [residue.get_resname() for residue in model_chain.get_residues()]
        native_sequence = [residue.get_resname() for residue in native_chain.get_residues()]
        model_sequence = "".join(
            seq1(r) if len(r) == 3 else r[:-1] if (len(r) == 2) else r for r in model_sequence
        )

        native_sequence = "".join(
            seq1(r) if len(r) == 3 else r[:-1] if len(r) == 2 else r for r in native_sequence
        )

    aligner = Align.PairwiseAligner()
    aligner.match = 5
    aligner.mismatch = 0
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -0.5
    aln = aligner.align(model_sequence, native_sequence)[0]
    return aln


def format_alignment(aln):
    alignment = {}
    try:
        alignment["seqA"] = aln[0, :]
        alignment["matches"] = "".join(["|" if aa1 == aa2 else
                                        " " if (aa1 == "-" or aa2 == "-") else
                                        "." for aa1, aa2 in
                                        zip(aln[0, :], aln[1, :])])
        alignment["seqB"] = aln[1, :]
    except NotImplementedError:
        formatted_aln = aln.format().split("\n")
        alignment["seqA"] = formatted_aln[0]
        alignment["matches"] = formatted_aln[1]
        alignment["seqB"] = formatted_aln[2]

    return alignment


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

    for chain in chains:
        if not model[chain]:
            model.detach_child(chain)


def remove_h(model):
    chains = [chain.id for chain in model.get_chains()]
    atoms_to_delete = []

    for chain in chains:
        residues = model[chain].get_residues()
        for residue in residues:
            for atom in residue.get_atoms():
                if atom.element == "H":
                    atoms_to_delete.append(atom.get_full_id())

    for _, _, chain, res, atom in atoms_to_delete:
        model[chain][res].detach_child(atom[0])


@lru_cache
def list_atoms_per_residue(chain, what):
    n_atoms_per_residue = []

    for residue in chain:
        # important to remove duplicate atoms (e.g. alternates) at this stage
        atom_ids = set([a.id for a in residue.get_unpacked_list()])
        n_atoms_per_residue.append(len(atom_ids))
    return np.array(n_atoms_per_residue).astype(int)


@lru_cache
def get_atoms_per_residue(
    chains,
    what,
    atom_types=("CA", "C", "N", "O", "P"),
):
    chain1, chain2 = chains
    atoms1 = (
        [
            atom.coord
            for res1, res2 in zip(chain1, chain2)
            for atom in res1.get_atoms()
            if atom.id in atom_types
            and atom.id in [a.id for a in res2.get_atoms()]
        ]
    )

    atoms2 = (
        [
            atom.coord
            for res1, res2 in zip(chain1, chain2)
            for atom in res2.get_atoms()
            if atom.id in atom_types
            and atom.id in [a.id for a in res1.get_atoms()]
        ]
    )
    return atoms1, atoms2


def get_interacting_pairs(distances, threshold):
    interacting_pairs = np.nonzero(np.asarray(distances) < threshold)
    return tuple(interacting_pairs[0]), tuple(interacting_pairs[1])


@lru_cache
# @profile
def get_interface_atoms(
    interacting_pairs,
    model_chains,
    ref_chains,
    atom_types=[],
):
    ref_interface = []
    mod_interface = []

    ref_residues_group1 = [res for res in ref_chains[0]]
    ref_residues_group2 = [res for res in ref_chains[1]]

    mod_residues_group1 = [res for res in model_chains[0]]
    mod_residues_group2 = [res for res in model_chains[1]]
    # remove duplicate residues
    interface_residues_group1 = set(interacting_pairs[0])
    interface_residues_group2 = set(interacting_pairs[1])

    for i in interface_residues_group1:
        ref_atoms = [atom for atom in ref_residues_group1[i].get_atoms()]
        mod_atoms = [atom for atom in mod_residues_group1[i].get_atoms()]
        ref_atoms_ids = [atom.id for atom in ref_atoms]
        mod_atoms_ids = [atom.id for atom in mod_atoms]
        ref_interface += [
            atom.coord
            for atom in ref_atoms
            if atom.id in atom_types and atom.id in mod_atoms_ids
        ]
        mod_interface += [
            atom.coord
            for atom in mod_atoms
            if atom.id in atom_types and atom.id in ref_atoms_ids
        ]

    for j in interface_residues_group2:
        ref_atoms = [atom for atom in ref_residues_group2[j].get_atoms()]
        mod_atoms = [atom for atom in mod_residues_group2[j].get_atoms()]
        ref_atoms_ids = [atom.id for atom in ref_atoms]
        mod_atoms_ids = [atom.id for atom in mod_atoms]
        ref_interface += [
            atom.coord
            for atom in ref_atoms
            if atom.id in atom_types and atom.id in mod_atoms_ids
        ]
        mod_interface += [
            atom.coord
            for atom in mod_atoms
            if atom.id in atom_types and atom.id in ref_atoms_ids
        ]

    return np.asarray(mod_interface), np.asarray(ref_interface)


@lru_cache
# @profile
def run_on_chains(
    model_chains,
    native_chains,
    no_needle=False,
    use_CA_only=False,
    capri_peptide=False,
):
    # realign each model chain against the corresponding native chain
    alignments = []
    for model_chain, native_chain in zip(model_chains, native_chains):
        aln = align_chains(
            model_chain,
            native_chain,
            use_numbering=no_needle,
        )
        alignment = format_alignment(aln)
        alignments.append(tuple(alignment.values()))

    info = calc_DockQ(
        model_chains,
        native_chains,
        alignments=tuple(alignments),
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
    """Given a native-model chain map, finds all non-null native interfaces
    and runs DockQ for each native-model pair of interfaces"""
    results_dic = {}
    native_chain_ids = list(chain_map.keys())

    for chain_pair in itertools.combinations(native_chain_ids, 2):
        native_chains = tuple([native_structure[chain] for chain in chain_pair])
        model_chains = tuple(
            [
                model_structure[chain]
                for chain in [chain_map[chain_pair[0]], chain_map[chain_pair[1]]]
            ]
        )
        if len(set(model_chains)) < 2:
            continue
        if chain_pair[0] in chain_map and chain_pair[1] in chain_map:
            info = run_on_chains(
                model_chains,
                native_chains,
                no_needle=no_needle,
                use_CA_only=use_CA_only,
                capri_peptide=capri_peptide,
            )
            if info:
                info["chain1"], info["chain2"] = chain_map[chain_pair[0]], chain_map[chain_pair[1]]
                results_dic[chain_pair] = info

    return results_dic


#@profile
def load_PDB(path, chains=[], n_model=0):
    try:
        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure("-", (gzip.open if path.endswith(".gz") else open)(path, "rt"), chains=chains)
        model = structure[n_model]
    except Exception:
        pdb_parser = Bio.PDB.MMCIFParser(QUIET=True)
        structure = pdb_parser.get_structure("-", (gzip.open if path.endswith(".gz") else open)(path, "rt"))
        model = structure[n_model]

    #remove_hetatms(model)
    remove_h(model)
    return model


def group_model_chains(model_structure, native_structure, model_chains, native_chains,allowed_mismatches=0):
    alignment_targets = itertools.product(model_chains, native_chains)
    native_chain_clusters = {chain: [] for chain in native_chains}
    mismatches_={} # for diagnostics
    for model_chain, native_chain in alignment_targets:
        aln = align_chains(model_structure[model_chain], native_structure[native_chain])
        alignment = format_alignment(aln)
        mismatches=alignment["matches"].count('.')
        if mismatches > 0 and mismatches < 10:
            mismatches_[(model_chain,native_chain)]=mismatches

        if mismatches <= allowed_mismatches: 
            # 100% sequence identity, 100% coverage of native sequence in model sequence
            native_chain_clusters[native_chain].append(model_chain)
    chains_without_match=[native_chain for native_chain in native_chain_clusters if len(native_chain_clusters[native_chain])==0]
    if len(chains_without_match):
        print(f'For these chains {chains_without_match} no match was found between model and native, try increasing the --allowed_mismatches from {allowed_mismatches}')
        print(f'Current number of alignments with 1-10 mismatches: {mismatches_}')   
    return native_chain_clusters


def format_mapping(mapping_str):
    mapping = None
    model_chains = None
    native_chains = None
    if not mapping_str:
        return mapping, model_chains, native_chains

    model_mapping, native_mapping = mapping_str.split(":")
    if not native_mapping:
        print("When using --mapping, native chains must be set (e.g. ABC:ABC or :ABC)")
        sys.exit()
    else:
        # :ABC or *:ABC only use those natives chains, permute model chains
        if not model_mapping or model_mapping == "*":
            native_chains = [chain for chain in native_mapping]
        elif len(model_mapping) == len(native_mapping):
            # ABC*:ABC* fix the first part of the mapping, try all other combinations
            mapping = {nm: mm for nm, mm in zip(native_mapping, model_mapping) if nm != "*" and mm != "*"}
            if model_mapping[-1] != "*" and native_mapping[-1] != "*":
                # ABC:ABC use the specific mapping
                model_chains = [chain for chain in model_mapping]
                native_chains = [chain for chain in native_mapping]
    return mapping, model_chains, native_chains

def format_mapping_string(chain_mapping, model_first=True):
    chain1=""
    chain2=""
    if model_first:
        mapping=sorted([(b,a) for a,b in chain_mapping.items()])
    else:
        mapping=sorted(chain_mapping.items())   
    for model_chain,native_chain in mapping: #sorted([(b,a) for a,b in chain_mapping.items()]):
        chain1+=model_chain
        chain2+=native_chain
  
    return f'{chain1}:{chain2}'
#@profile
def main():
    args = parse_args()
    initial_mapping, model_chains, native_chains = format_mapping(args.mapping)

    model_structure = load_PDB(args.model, chains=model_chains)
    native_structure = load_PDB(args.native, chains=native_chains)

    info = {}
    model_chains = [c.id for c in model_structure] if not model_chains else model_chains
    native_chains = [c.id for c in native_structure] if not native_chains else native_chains

    if len(model_chains) < 2 or len(native_chains) < 2:
        print("Need at least two chains in the two inputs\n")
        sys.exit()

    # permute chains and run on a for loop
    best_dockq = -1
    best_result = None
    best_mapping = None
    ref_switch=False
    native_chain_clusters = group_model_chains(
        model_structure, native_structure, model_chains, native_chains,args.allowed_mismatches
    )
    all_mappings = itertools.product(*[cluster for cluster in native_chain_clusters.values() if cluster])

    # remove mappings where the same model chain is present more than once
    # only if the mapping is supposed to be 1-1
    if len(model_chains) == len(native_chains):
        all_mappings = [
            element for element in all_mappings if len(set(element)) == len(element)
        ]

    for mapping in all_mappings:
        chain_map = {
            native_chain: mapping[i] for i, native_chain in enumerate(native_chains)
        }

        if initial_mapping and not initial_mapping.items() <= chain_map.items():
            continue

        result_this_mapping = run_on_all_native_interfaces(
            model_structure,
            native_structure,
            chain_map=chain_map,
            no_needle=args.no_needle,
            use_CA_only=args.use_CA,
            capri_peptide=args.capri_peptide,
        )
        if args.optDockQF1:
            total_dockq = sum(
                [result["DockQ_F1"] for result in result_this_mapping.values()]
                )
        else:
             total_dockq = sum(
                [result["DockQ"] for result in result_this_mapping.values()]
                )
        if total_dockq > best_dockq:
            best_dockq = total_dockq
            best_result = result_this_mapping
            best_mapping = chain_map

    info["model"] = args.model
    info["native"] = args.native
    info["best_dockq"] = best_dockq
    info["best_result"] = best_result
    info["GlobalDockQ"] = best_dockq/len(best_result)
    info['best_mapping']=best_mapping
    info['best_mapping_str']=f'{format_mapping_string(best_mapping,model_first=not ref_switch)} {best_mapping}' 
    print_results(info, args.short, args.verbose, args.capri_peptide)


def print_results(info, short=False, verbose=False, capri_peptide=False):
    if short:
        capri_peptide_str = "-capri_peptide" if capri_peptide else ""
        print(f"Total DockQ over {len(info['best_result'])} native interfaces: {info['GlobalDockQ']:.3f} with {info['best_mapping_str']} model:native mapping")
        print(info['best_result'])
        for chains, results in info["best_result"].items():
            print(
                f"DockQ{capri_peptide_str} {results['DockQ']:.3f} DockQ_F1 {results['DockQ_F1']:.3f} Fnat {results['fnat']:.3f} iRMS {results['irms']:.3f} LRMS {results['Lrms']:.3f} Fnonnat {results['fnonnat']:.3f} {info['native']} {chains[0]} {chains[1]} -> {info['model']} {results['chain1']} {results['chain2']}"
            )

    else:
        print_header(verbose, capri_peptide)
        print(f"Model  : {info['model']}")
        print(f"Native : {info['native']}")
        print(f"Total DockQ over {len(info['best_result'])} native interfaces: {info['best_dockq']:.3f}")
        items = ["DockQ_F1", "DockQ", "irms", "Lrms", "fnat"]

        for chains, results in info["best_result"].items():
            print(f"Native chains: {chains[0]}, {chains[1]}")
            print(f"\tModel chains: {results['chain1']}, {results['chain2']}")
            print("\n".join([f"\t{item}: {results[item]:.3f}" for item in items]))


def print_header(verbose=False, capri_peptide=False):
    reference = "*   Ref: S. Basu and B. Wallner, DockQ: A quality measure for  *\n" \
                "*   protein-protein docking models                             *\n" \
                "*                            doi:10.1371/journal.pone.0161879  *\n" \
                "*   For comments, please email: bjorn.wallner@.liu.se          *"
    if not capri_peptide:
        header = "****************************************************************\n" \
                 "*                       DockQ                                  *\n" \
                 "*   Scoring function for protein-protein docking models        *\n" \
                 "*   Statistics on CAPRI data:                                  *\n" \
                 "*    0.00 <= DockQ <  0.23 - Incorrect                         *\n" \
                 "*    0.23 <= DockQ <  0.49 - Acceptable quality                *\n" \
                 "*    0.49 <= DockQ <  0.80 - Medium quality                    *\n" \
                 "*            DockQ >= 0.80 - High quality                      *"
    else:
        header = "****************************************************************\n" \
                 "*                DockQ-CAPRI peptide                           *\n" \
                 "*   Do not trust any thing you read....                        *\n" \
                 "*   OBS THE DEFINITION OF Fnat and iRMS are different for      *\n" \
                 "*   peptides in CAPRI                                          *\n" \
                 "*                                                              *"

    if verbose:
        notice = "*   For the record:                                            *\n" \
                 f"*   Definition of contact <{'5A' if not capri_peptide else '4A'} (Fnat)                           *\n" \
                 f"*   Definition of interface <{'10A all heavy atoms (iRMS)       ' if not capri_peptide else '8A CB (iRMS)                     '} *\n" \
                 "****************************************************************"
    else:
        notice = "****************************************************************"

    print(header)
    print(reference)
    print(notice)


if __name__ == "__main__":
    main()
