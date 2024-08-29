#!/usr/bin/env python

import sys
import json
import gzip
import math
import logging
import itertools
import traceback
from collections import Counter
from argparse import ArgumentParser
from functools import lru_cache, partial

import numpy as np
from Bio import Align
from Bio.SVDSuperimposer import SVDSuperimposer
from parallelbar import progress_map

# fallback in case the cython version doesn't work, though it will be slower
try:
    from .operations import residue_distances, get_fnat_stats
    from .parsers import PDBParser, MMCIFParser
    from .constants import *
except ImportError:
    logging.warning(
        """Invoking DockQ as script rather than binary. This will slow down computations."""
    )
    from operations_nocy import residue_distances, get_fnat_stats
    from parsers import PDBParser, MMCIFParser
    from constants import *


def parse_args():
    parser = ArgumentParser(
        description="DockQ - Quality measure for \
        protein-protein docking models"
    )
    parser.add_argument("model", metavar="<model>", type=str, help="Path to model file")
    parser.add_argument(
        "native", metavar="<native>", type=str, help="Path to native file"
    )
    parser.add_argument(
        "--capri_peptide",
        default=False,
        action="store_true",
        help="use version for capri_peptide \
        (DockQ cannot not be trusted for this setting)",
    )
    parser.add_argument(
        "--small_molecule",
        help="If the docking pose of a small molecule should be evaluated",
        action="store_true",
    )
    parser.add_argument(
        "--short", default=False, action="store_true", help="Short output"
    )
    parser.add_argument(
        "--json",
        default=None,
        metavar="out.json",
        help="Write outputs to a chosen json file",
    )
    parser.add_argument(
        "--verbose", "-v", default=False, action="store_true", help="Verbose output"
    )
    parser.add_argument(
        "--no_align",
        action="store_true",
        help="Do not align native and model using sequence alignments, but use the numbering of residues instead",
    )
    parser.add_argument(
        "--n_cpu",
        default=8,
        type=int,
        metavar="CPU",
        help="Number of cores to use",
    )
    parser.add_argument(
        "--max_chunk",
        default=512,
        type=int,
        metavar="CHUNK",
        help="Maximum size of chunks given to the cores, actual chunksize is min(max_chunk,combos/cpus)",
    )
    parser.add_argument(
        "--allowed_mismatches",
        default=0,
        type=int,
        help="Number of allowed mismatches when mapping model sequence to native sequence.",
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
            or the equivalent '--mapping *:BC'.""",
    )

    return parser.parse_args()


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
def get_residue_distances(chain1, chain2, what, all_atom=True):
    if all_atom:
        # how many atoms per aligned amino acid
        n_atoms_per_res_chain1 = list_atoms_per_residue(chain1, what)
        n_atoms_per_res_chain2 = list_atoms_per_residue(chain2, what)
        model_A_atoms = np.asarray(
            [atom.coord for res in chain1 for atom in res.get_atoms()]
        )
        model_B_atoms = np.asarray(
            [atom.coord for res in chain2 for atom in res.get_atoms()]
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


def calc_sym_corrected_lrmsd(
    sample_chains,
    ref_chains,
    alignments,
):
    import networkx as nx

    is_het_sample_0 = bool(sample_chains[0].is_het)
    is_het_sample_1 = bool(sample_chains[1].is_het)

    if is_het_sample_0 and not is_het_sample_1:
        sample_ligand = sample_chains[0]
        sample_receptor = sample_chains[1]
        ref_ligand = ref_chains[0]
        ref_receptor = ref_chains[1]
        receptor_alignment = alignments[1]
    elif not is_het_sample_0 and is_het_sample_1:
        sample_ligand = sample_chains[1]
        sample_receptor = sample_chains[0]
        ref_ligand = ref_chains[1]
        ref_receptor = ref_chains[0]
        receptor_alignment = alignments[0]
    else:
        return  # both ligands, no lrmsd

    ref_res_distances = get_residue_distances(ref_receptor, ref_ligand, "ref")
    receptor_interface, _ = get_interacting_pairs(
        # working with squared thresholds to avoid using sqrt
        ref_res_distances,
        threshold=INTERFACE_THRESHOLD ** 2,
    )
    if not receptor_interface:
        return
    aligned_sample_receptor, aligned_ref_receptor = get_aligned_residues(
        sample_receptor, ref_receptor, receptor_alignment
    )

    sample_interface_atoms, ref_interface_atoms = subset_atoms(
        aligned_sample_receptor,
        aligned_ref_receptor,
        atom_types=BACKBONE_ATOMS,
        residue_subset=receptor_interface,
        what="receptor",
    )

    sample_ligand_atoms_ids = [atom.id for atom in sample_ligand.get_atoms()]
    sample_ligand_atoms_ele = [atom.element for atom in sample_ligand.get_atoms()]

    ref_ligand_atoms_ids = [atom.id for atom in ref_ligand.get_atoms()]
    ref_ligand_atoms_ele = [atom.element for atom in ref_ligand.get_atoms()]

    sample_ligand_atoms = np.array(
        [
            atom.coord
            for atom in sample_ligand.get_atoms()
            if atom.id in ref_ligand_atoms_ids
        ]
    )
    ref_ligand_atoms = np.array(
        [
            atom.coord
            for atom in ref_ligand.get_atoms()
            if atom.id in sample_ligand_atoms_ids
        ]
    )

    # Set to align on receptor interface
    super_imposer = SVDSuperimposer()
    super_imposer.set(
        np.asarray(ref_interface_atoms), np.asarray(sample_interface_atoms)
    )
    super_imposer.run()
    rot, tran = super_imposer.get_rotran()

    sample_rotated_ligand_atoms = np.dot(sample_ligand_atoms, rot) + tran

    sample_graph = create_graph(sample_ligand_atoms, sample_ligand_atoms_ele)
    ref_graph = create_graph(ref_ligand_atoms, ref_ligand_atoms_ele)

    min_lrms = float("inf")
    best_mapping = None

    for isomorphism in nx.vf2pp_all_isomorphisms(sample_graph, ref_graph):
        model_i = list(isomorphism.keys())
        native_i = list(isomorphism.values())

        lrms = super_imposer._rms(
            sample_rotated_ligand_atoms[model_i], ref_ligand_atoms[native_i]
        )
        if lrms < min_lrms:
            best_mapping = isomorphism
            min_lrms = lrms
    dockq = dockq_formula(0, 0, min_lrms)
    info = {
        "DockQ": dockq,
        "LRMSD": min_lrms,
        "mapping": best_mapping,
        "is_het": sample_ligand.is_het,
    }
    return info


# @profile
def calc_DockQ(
    sample_chains,
    ref_chains,
    alignments,
    capri_peptide=False,
    low_memory=False,
):

    fnat_threshold = FNAT_THRESHOLD if not capri_peptide else FNAT_THRESHOLD_PEPTIDE
    interface_threshold = (
        INTERFACE_THRESHOLD if not capri_peptide else INTERFACE_THRESHOLD_PEPTIDE
    )

    # total number of native contacts is calculated on untouched native structure
    ref_res_distances = get_residue_distances(ref_chains[0], ref_chains[1], "ref")
    nat_total = np.nonzero(np.asarray(ref_res_distances) < fnat_threshold ** 2)[
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

    sample_res_distances = get_residue_distances(
        aligned_sample_1, aligned_sample_2, "sample"
    )

    if ref_res_distances.shape != sample_res_distances.shape:
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
            aligned_ref_1, aligned_ref_2, "ref", all_atom=False
        )
    # Get interfacial atoms from reference, and corresponding atoms from sample
    interacting_pairs = get_interacting_pairs(
        # working with squared thresholds to avoid using sqrt
        ref_res_distances,
        threshold=interface_threshold ** 2,
    )

    sample_interface_atoms1, ref_interface_atoms1 = subset_atoms(
        aligned_sample_1,
        aligned_ref_1,
        atom_types=BACKBONE_ATOMS,
        residue_subset=interacting_pairs[0],
    )
    sample_interface_atoms2, ref_interface_atoms2 = subset_atoms(
        aligned_sample_2,
        aligned_ref_2,
        atom_types=BACKBONE_ATOMS,
        residue_subset=interacting_pairs[1],
    )

    sample_interface_atoms = np.asarray(
        sample_interface_atoms1 + sample_interface_atoms2
    )
    ref_interface_atoms = np.asarray(ref_interface_atoms1 + ref_interface_atoms2)

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

    receptor_atoms_native, receptor_atoms_sample = subset_atoms(
        receptor_chains[0],
        receptor_chains[1],
        atom_types=BACKBONE_ATOMS,
        what="receptor",
    )
    ligand_atoms_native, ligand_atoms_sample = subset_atoms(
        ligand_chains[0], ligand_chains[1], atom_types=BACKBONE_ATOMS, what="ligand"
    )
    # Set to align on receptor
    super_imposer.set(
        np.asarray(receptor_atoms_native), np.asarray(receptor_atoms_sample)
    )
    super_imposer.run()

    rot, tran = super_imposer.get_rotran()
    rotated_sample_atoms = np.dot(np.asarray(ligand_atoms_sample), rot) + tran

    lrms = super_imposer._rms(
        np.asarray(ligand_atoms_native), rotated_sample_atoms
    )  # using the private _rms function which does not superimpose

    info = {}
    F1 = f1(nat_correct, nonnat_count, nat_total)
    info["DockQ"] = dockq_formula(fnat, irms, lrms)
    if low_memory:
        return info

    info["F1"] = F1
    info["iRMSD"] = irms
    info["LRMSD"] = lrms
    info["fnat"] = fnat
    info["nat_correct"] = nat_correct
    info["nat_total"] = nat_total

    info["fnonnat"] = fnonnat
    info["nonnat_count"] = nonnat_count
    info["model_total"] = model_total
    info["clashes"] = np.nonzero(
        np.asarray(sample_res_distances) < CLASH_THRESHOLD ** 2
    )[0].shape[0]
    info["len1"] = ref_group1_size
    info["len2"] = ref_group2_size
    info["class1"] = class1
    info["class2"] = class2
    info["is_het"] = False

    return info


def f1(tp, fp, p):
    return 2 * tp / (tp + fp + p)


def dockq_formula(fnat, irms, lrms):
    return (
        float(fnat)
        + 1 / (1 + (irms / 1.5) * (irms / 1.5))
        + 1 / (1 + (lrms / 8.5) * (lrms / 8.5))
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
        # if the smallest resn is negative, it will be used to shift all numbers so they start from 0
        # the minimum offset is 45 to avoid including the "-" character that is reserved for gaps
        min_resn = max(45, -min(model_numbering + native_numbering))

        model_sequence = "".join([chr(resn + min_resn) for resn in model_numbering])
        native_sequence = "".join([chr(resn + min_resn) for resn in native_numbering])

    else:
        model_sequence = model_chain.sequence
        native_sequence = native_chain.sequence

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
        alignment["matches"] = "".join(
            [
                "|" if aa1 == aa2 else " " if (aa1 == "-" or aa2 == "-") else "."
                for aa1, aa2 in zip(aln[0, :], aln[1, :])
            ]
        )
        alignment["seqB"] = aln[1, :]
    except NotImplementedError:
        formatted_aln = aln.format().split("\n")
        alignment["seqA"] = formatted_aln[0]
        alignment["matches"] = formatted_aln[1]
        alignment["seqB"] = formatted_aln[2]

    return alignment


@lru_cache
def list_atoms_per_residue(chain, what):
    n_atoms_per_residue = []

    for residue in chain:
        # important to remove duplicate atoms (e.g. alternates) at this stage
        atom_ids = set([a.id for a in residue.get_unpacked_list()])
        n_atoms_per_residue.append(len(atom_ids))
    return np.array(n_atoms_per_residue).astype(int)


def get_interacting_pairs(distances, threshold):
    interacting_pairs = np.nonzero(np.asarray(distances) < threshold)
    return tuple(interacting_pairs[0]), tuple(interacting_pairs[1])


@lru_cache
def subset_atoms(
    mod_chain,
    ref_chain,
    atom_types,
    residue_subset=None,
    what="",
):
    mod_atoms = []
    ref_atoms = []

    mod_residues = [res for res in mod_chain]
    ref_residues = [res for res in ref_chain]

    # remove duplicate residues
    residue_subset = set(residue_subset) if residue_subset else range(len(mod_residues))

    for i in residue_subset:
        mod_res_atoms = list(mod_residues[i].get_atoms())
        ref_res_atoms = list(ref_residues[i].get_atoms())
        mod_res_atoms_ids = [atom.id for atom in mod_res_atoms]
        ref_res_atoms_ids = [atom.id for atom in ref_res_atoms]

        for atom_type in atom_types:
            try:
                mod_i = mod_res_atoms_ids.index(atom_type)
                ref_i = ref_res_atoms_ids.index(atom_type)
                mod_atoms += [mod_res_atoms[mod_i].coord]
                ref_atoms += [ref_res_atoms[ref_i].coord]
            except:
                continue

    return mod_atoms, ref_atoms


@lru_cache
def run_on_chains(
    model_chains,
    native_chains,
    no_align=False,
    capri_peptide=False,
    small_molecule=True,
    low_memory=False,
):
    # realign each model chain against the corresponding native chain
    alignments = []
    for model_chain, native_chain in zip(model_chains, native_chains):
        aln = align_chains(
            model_chain,
            native_chain,
            use_numbering=no_align,
        )
        alignment = format_alignment(aln)
        alignments.append(tuple(alignment.values()))

    if not small_molecule:
        info = calc_DockQ(
            model_chains,
            native_chains,
            alignments=tuple(alignments),
            capri_peptide=capri_peptide,
            low_memory=low_memory,
        )
    else:
        info = calc_sym_corrected_lrmsd(
            model_chains,
            native_chains,
            alignments=tuple(alignments),
        )
    return info


def create_graph(atom_list, atom_ids):
    import networkx as nx

    G = nx.Graph()

    for i, atom_i in enumerate(atom_list):
        cr_i = COVALENT_RADIUS[atom_ids[i]]
        for j, atom_j in enumerate(atom_list):
            cr_j = COVALENT_RADIUS[atom_ids[j]]
            distance = np.linalg.norm(atom_i - atom_j)
            threshold = (cr_i + cr_j + BOND_TOLERANCE) if i != j else 1
            if distance < threshold:  # Adjust threshold as needed
                G.add_edge(i, j)

    return G


def run_on_all_native_interfaces(
    model_structure,
    native_structure,
    chain_map={"A": "A", "B": "B"},
    no_align=False,
    capri_peptide=False,
    low_memory=False,
):
    """Given a native-model chain map, finds all non-null native interfaces
    and runs DockQ for each native-model pair of interfaces"""
    result_mapping = dict()
    native_chain_ids = list(chain_map.keys())

    for chain_pair in itertools.combinations(native_chain_ids, 2):
        native_chains = tuple([native_structure[chain] for chain in chain_pair])
        model_chains = tuple(
            [
                model_structure[chain]
                for chain in [chain_map[chain_pair[0]], chain_map[chain_pair[1]]]
            ]
        )

        small_molecule = native_chains[0].is_het or native_chains[1].is_het

        if len(set(model_chains)) < 2:
            continue
        if chain_pair[0] in chain_map and chain_pair[1] in chain_map:
            info = run_on_chains(
                model_chains,
                native_chains,
                no_align=no_align,
                capri_peptide=capri_peptide,
                small_molecule=small_molecule,
                low_memory=low_memory,
            )
            if info:
                info["chain1"], info["chain2"] = (
                    chain_map[chain_pair[0]],
                    chain_map[chain_pair[1]],
                )
                info["chain_map"] = chain_map  # diagnostics
                result_mapping["".join(chain_pair)] = info
    total_dockq = sum([result["DockQ"] for result in result_mapping.values()])
    return result_mapping, total_dockq


def load_PDB(path, chains=[], small_molecule=False, n_model=0):
    try:
        pdb_parser = PDBParser(QUIET=True)
        model = pdb_parser.get_structure(
            "-",
            (gzip.open if path.endswith(".gz") else open)(path, "rt"),
            chains=chains,
            parse_hetatms=small_molecule,
            model_number=n_model,
        )
    except Exception:
        pdb_parser = MMCIFParser(QUIET=True)
        model = pdb_parser.get_structure(
            "-",
            (gzip.open if path.endswith(".gz") else open)(path, "rt"),
            chains=chains,
            parse_hetatms=small_molecule,
            auth_chains=not small_molecule,
            model_number=n_model,
        )
    model.id = path
    return model


def group_chains(
    query_structure, ref_structure, query_chains, ref_chains, allowed_mismatches=0
):
    reverse_map = False
    mismatch_dict = {}  # for diagnostics
    # this might happen e.g. when modelling only part of a large homomer
    if len(query_chains) < len(ref_chains):
        query_structure, ref_structure = ref_structure, query_structure
        query_chains, ref_chains = ref_chains, query_chains
        reverse_map = True

    alignment_targets = itertools.product(query_chains, ref_chains)
    chain_clusters = {chain: [] for chain in ref_chains}

    for query_chain, ref_chain in alignment_targets:
        try:
            qc = query_structure[query_chain]
        except KeyError:
            logging.error(
                f"""The specified model chain {query_chain} is not found in the PDB structure.
This is possibly due to using the wrong chain identifier in --mapping,
or forgetting to specify --small_molecule if this is a HETATM chain.
If working with mmCIF files, make sure you use the right chain identifier.
            """
            )
            print(traceback.format_exc())
            sys.exit(1)
        try:
            rc = ref_structure[ref_chain]
        except KeyError:
            logging.error(
                f"""The specified native chain {ref_chain} is not found in the PDB structure.
This is possibly due to using the wrong chain identifier in --mapping,
or forgetting to specify --small_molecule if this is a HETATM chain.
If working with mmCIF files, make sure you use the right chain identifier.
            """
            )
        het_qc = qc.is_het
        het_rc = rc.is_het

        if het_qc is None and het_rc is None:
            aln = align_chains(
                qc,
                rc,
                use_numbering=False,
            )
            alignment = format_alignment(aln)
            n_mismatches = alignment["matches"].count(".")

            if 0 < n_mismatches < 10:
                mismatch_dict[(query_chain, ref_chain)] = n_mismatches

            if n_mismatches <= allowed_mismatches:
                # 100% sequence identity, 100% coverage of native sequence in model sequence
                chain_clusters[ref_chain].append(query_chain)
        elif het_qc and het_rc and het_qc == het_rc:
            chain_clusters[ref_chain].append(query_chain)
    chains_without_match = [
        chain for chain in chain_clusters if not chain_clusters[chain]
    ]

    if mismatch_dict:
        logging.warning(
            f"""Some chains have a limited number of sequence mismatches and are treated as non-homologous. 
Try increasing the --allowed_mismatches for the following: {", ".join(f"Model chain {c[1]}, native chain {c[0]}: {m} mismatches" for c, m in mismatch_dict.items())}
if they should be treated as homologous."""
        )

    if chains_without_match:
        logging.error(
            f"For chains {chains_without_match} no identical corresponding chain was found between in the native."
        )
        sys.exit(1)

    return chain_clusters, reverse_map


def format_mapping(mapping_str, small_molecule=None):
    mapping = dict()
    model_chains = None
    native_chains = None
    if not mapping_str:
        return mapping, model_chains, native_chains

    model_mapping, native_mapping = mapping_str.split(":")
    if not native_mapping:
        logging.error(
            "When using --mapping, native chains must be set (e.g. ABC:ABC or :ABC)"
        )
        sys.exit()
    else:
        # :ABC or *:ABC only use those natives chains, permute model chains
        if not model_mapping or model_mapping == "*":
            native_chains = [chain for chain in native_mapping]
        elif len(model_mapping) == len(native_mapping):
            # ABC*:ABC* fix the first part of the mapping, try all other combinations
            mapping = {
                nm: mm
                for nm, mm in zip(native_mapping, model_mapping)
                if nm != "*" and mm != "*"
            }
            if model_mapping[-1] != "*" and native_mapping[-1] != "*":
                # ABC:ABC use the specific mapping
                model_chains = [chain for chain in model_mapping]
                native_chains = [chain for chain in native_mapping]
    return mapping, model_chains, native_chains


def format_mapping_string(chain_mapping):
    chain1 = ""
    chain2 = ""

    # Sorting might change LRMSD since the definition of receptor/ligand for equal length depends on order
    mapping = [(b, a) for a, b in chain_mapping.items()]
    for (
        model_chain,
        native_chain,
    ) in mapping:
        chain1 += model_chain
        chain2 += native_chain

    return f"{chain1}:{chain2}"


def product_without_dupl(*args, repeat=1):
    pools = [tuple(pool) for pool in args] * repeat
    result = [[]]
    for pool in pools:
        result = [
            x + [y] for x in result for y in pool if y not in x
        ]  # here we added condition

    for prod in result:
        yield tuple(prod)


def count_chain_combinations(chain_clusters):
    try:
        clusters = [tuple(li) for li in chain_clusters.values()]
        number_of_combinations = np.prod(
            [
                int(math.factorial(len(a)) / math.factorial(len(a) - b))
                for a, b in Counter(clusters).items()
            ]
        )
    except ValueError:
        logging.error(
            """Couldn't find a match between each model-native chain specified in the mapping.
Make sure that all chains in your model have a homologous chain in the native, or specify the right subset of chains with --mapping"""
        )
        sys.exit()
    return number_of_combinations


def get_all_chain_maps(
    chain_clusters,
    initial_mapping,
    reverse_map,
    model_chains_to_combo,
    native_chains_to_combo,
):
    all_mappings = product_without_dupl(
        *[cluster for cluster in chain_clusters.values() if cluster]
    )
    for mapping in all_mappings:
        chain_map = {key: value for key, value in initial_mapping.items()}
        if reverse_map:
            chain_map.update(
                {
                    mapping[i]: model_chain
                    for i, model_chain in enumerate(model_chains_to_combo)
                }
            )
        else:
            chain_map.update(
                {
                    native_chain: mapping[i]
                    for i, native_chain in enumerate(native_chains_to_combo)
                }
            )
        yield chain_map


def get_chain_map_from_dockq(result):
    chain_map = {}
    for ch1, ch2 in result:
        chain_map[ch1] = result[ch1, ch2]["chain1"]
        chain_map[ch2] = result[ch1, ch2]["chain2"]
    return chain_map


# @profile
def main():
    args = parse_args()

    initial_mapping, model_chains, native_chains = format_mapping(
        args.mapping, args.small_molecule
    )
    model_structure = load_PDB(
        args.model, chains=model_chains, small_molecule=args.small_molecule
    )
    native_structure = load_PDB(
        args.native, chains=native_chains, small_molecule=args.small_molecule
    )
    # check user-given chains are in the structures
    model_chains = [c.id for c in model_structure] if not model_chains else model_chains
    native_chains = (
        [c.id for c in native_structure] if not native_chains else native_chains
    )

    if len(model_chains) < 2 or len(native_chains) < 2:
        print("Need at least two chains in the two inputs\n")
        sys.exit()

    # permute chains and run on a for loop
    best_dockq = -1
    best_result = None
    best_mapping = None

    model_chains_to_combo = [
        mc for mc in model_chains if mc not in initial_mapping.values()
    ]
    native_chains_to_combo = [
        nc for nc in native_chains if nc not in initial_mapping.keys()
    ]

    chain_clusters, reverse_map = group_chains(
        model_structure,
        native_structure,
        model_chains_to_combo,
        native_chains_to_combo,
        args.allowed_mismatches,
    )
    chain_maps = get_all_chain_maps(
        chain_clusters,
        initial_mapping,
        reverse_map,
        model_chains_to_combo,
        native_chains_to_combo,
    )

    num_chain_combinations = count_chain_combinations(chain_clusters)
    # copy iterator to use later
    chain_maps, chain_maps_ = itertools.tee(chain_maps)

    low_memory = num_chain_combinations > 100
    run_chain_map = partial(
        run_on_all_native_interfaces,
        model_structure,
        native_structure,
        no_align=args.no_align,
        capri_peptide=args.capri_peptide,
        low_memory=low_memory,
    )

    if num_chain_combinations > 1:
        cpus = min(num_chain_combinations, args.n_cpu)
        chunk_size = min(args.max_chunk, max(1, num_chain_combinations // cpus))

        # for large num_chain_combinations it should be possible to divide the chain_maps in chunks
        result_this_mappings = progress_map(
            run_chain_map,
            chain_maps,
            total=num_chain_combinations,
            n_cpu=cpus,
            chunk_size=chunk_size,
        )

        for chain_map, (result_this_mapping, total_dockq) in zip(
            chain_maps_, result_this_mappings
        ):

            if total_dockq > best_dockq:
                best_dockq = total_dockq
                best_result = result_this_mapping
                best_mapping = chain_map

        if low_memory:  # retrieve the full output by rerunning the best chain mapping
            best_result, total_dockq = run_on_all_native_interfaces(
                model_structure,
                native_structure,
                chain_map=best_mapping,
                no_align=args.no_align,
                capri_peptide=args.capri_peptide,
                low_memory=False,
            )

    else:  # skip multi-threading for single jobs (skip the bar basically)
        best_mapping = next(chain_maps)
        best_result, best_dockq = run_chain_map(best_mapping)

    if not best_result:
        logging.error(
            "Could not find interfaces in the native model. Please double check the inputs or select different chains with the --mapping flag."
        )
        sys.exit(1)

    info = dict()
    info["model"] = args.model
    info["native"] = args.native
    info["best_dockq"] = best_dockq
    info["best_result"] = best_result
    info["GlobalDockQ"] = best_dockq / len(best_result)
    info["best_mapping"] = best_mapping
    info["best_mapping_str"] = f"{format_mapping_string(best_mapping)}"

    if args.json:
        with open(args.json, "w") as fp:
            json.dump(info, fp)

    print_results(
        info, args.short, args.verbose, args.capri_peptide, args.small_molecule
    )


def print_results(
    info, short=False, verbose=False, capri_peptide=False, small_molecule=False
):

    score = (
        "DockQ-small_molecules"
        if small_molecule
        else "DockQ-capri_peptide"
        if capri_peptide
        else "DockQ"
    )
    if short:
        print(
            f"Total {score} over {len(info['best_result'])} native interfaces: {info['GlobalDockQ']:.3f} with {info['best_mapping_str']} model:native mapping"
        )
        for chains, results in info["best_result"].items():
            reported_measures = (
                [
                    "DockQ",
                    "iRMSD",
                    "LRMSD",
                    "fnat",
                    "fnonnat",
                    "F1",
                    "clashes",
                ]
                if not results["is_het"]
                else ["LRMSD"]
            )
            hetname = f" ({results['is_het']})" if results["is_het"] else ""
            score_str = " ".join(
                [
                    f"{item} {results[item]:.3f}"
                    if item != "clashes"
                    else f"{item} {results[item]}"
                    for item in reported_measures
                ]
            )
            print(
                f"{score_str} mapping {results['chain1']}{results['chain2']}:{chains[0]}{chains[1]}{hetname} {info['model']} {results['chain1']} {results['chain2']} -> {info['native']} {chains[0]} {chains[1]}"
            )
    else:
        print_header(verbose, capri_peptide)
        print(f"Model  : {info['model']}")
        print(f"Native : {info['native']}")
        print(
            f"Total {score} over {len(info['best_result'])} native interfaces: {info['GlobalDockQ']:.3f} with {info['best_mapping_str']} model:native mapping"
        )
        for chains, results in info["best_result"].items():
            reported_measures = (
                [
                    "DockQ",
                    "iRMSD",
                    "LRMSD",
                    "fnat",
                    "fnonnat",
                    "F1",
                    "clashes",
                ]
                if not results["is_het"]
                else ["LRMSD"]
            )
            hetname = f" ({results['is_het']})" if results["is_het"] else ""
            print(f"Native chains: {chains[0]}, {chains[1]}{hetname}")
            print(f"\tModel chains: {results['chain1']}, {results['chain2']}")
            print(
                "\n".join(
                    [
                        f"\t{item}: {results[item]:.3f}"
                        if item != "clashes"
                        else f"\t{item}: {results[item]}"
                        for item in reported_measures
                    ]
                )
            )


def print_header(verbose=False, capri_peptide=False):
    reference = (
        "*   Ref: Mirabello and Wallner, 'DockQ v2: Improved automatic  *\n"
        "*   quality measure for protein multimers, nucleic acids       *\n"
        "*   and small molecules'                                       *\n"
        "*                                                              *\n"
        "*   For comments, please email: bjorn.wallner@.liu.se          *"
    )

    header = (
        "****************************************************************\n"
        "*                       DockQ                                  *\n"
        "*   Docking scoring for biomolecular models                    *\n"
        "*   DockQ score legend:                                        *\n"
        "*    0.00 <= DockQ <  0.23 - Incorrect                         *\n"
        "*    0.23 <= DockQ <  0.49 - Acceptable quality                *\n"
        "*    0.49 <= DockQ <  0.80 - Medium quality                    *\n"
        "*            DockQ >= 0.80 - High quality                      *"
    )

    if verbose:
        notice = (
            "*   For the record:                                            *\n"
            f"*   Definition of contact <{'5A' if not capri_peptide else '4A'} (Fnat)                           *\n"
            f"*   Definition of interface <{'10A all heavy atoms (iRMSD)      ' if not capri_peptide else '8A CB (iRMSD)                    '} *\n"
            "****************************************************************"
        )
    else:
        notice = "****************************************************************"

    print(header)
    print(reference)
    print(notice)


if __name__ == "__main__":
    main()
