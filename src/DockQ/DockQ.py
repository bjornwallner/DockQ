#!/usr/bin/env python

import sys
import gzip
import hashlib
import warnings
import traceback
import itertools
import math
from functools import lru_cache, wraps, partial
from argparse import ArgumentParser
from tqdm import tqdm
from parallelbar import progress_map

import Bio.PDB
import numpy as np
from Bio import Align
from Bio.SeqUtils import seq1
from Bio.SVDSuperimposer import SVDSuperimposer

# fallback in case the cython version doesn't work, though it will be slower
try:
    from .operations import residue_distances, get_fnat_stats
    from .parsers import PDBParser, MMCIFParser
except ImportError:
    warnings.warn(
        """WARNING: It looks like cython is not working,
         falling back on native python. This will make DockQ slower"""
    )
    from operations_nocy import residue_distances, get_fnat_stats
    from parsers import PDBParser, MMCIFParser


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
        "--no_align",
        default=False,
        action="store_true",
        help="Do not align native and model using sequence alignments, but use the numbering of residues instead",
    )
    parser.add_argument(
        "--n_cpu",
        default=8,
        type=int,
        metavar="n_cpu",
        help="Number of cores to use",
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
            [atom.get_coord() for res in chain1 for atom in res.get_atoms()]
        )
        model_B_atoms = np.asarray(
            [atom.get_coord() for res in chain2 for atom in res.get_atoms()]
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


# @profile
def calc_DockQ(
    sample_chains,
    ref_chains,
    alignments,
    capri_peptide=False,
    low_memory=False,
):
    atom_for_sup = ("CA", "C", "N", "O", "P")
    fnat_threshold = 4.0 if capri_peptide else 5.0
    interface_threshold = 8.0 if capri_peptide else 10.0
    clash_threshold = 2.0
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
        threshold=interface_threshold ** 2,
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

    info["DockQ_F1"] = dockq_formula(
        f1(nat_correct, nonnat_count, nat_total), irms, Lrms
    )
    info["DockQ"] = dockq_formula(fnat, irms, Lrms)
    if low_memory:
        return info
    info["irms"] = irms
    info["Lrms"] = Lrms
    info["fnat"] = fnat
    info["nat_correct"] = nat_correct
    info["nat_total"] = nat_total

    info["fnonnat"] = fnonnat
    info["nonnat_count"] = nonnat_count
    info["model_total"] = model_total
    info["clashes"] = np.nonzero(
        np.asarray(sample_res_distances) < clash_threshold ** 2
    )[0].shape[0]
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
def align_chains(model_chain, native_chain, use_numbering=False, verbose=False):
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
        custom_map = {"MSE": "M", "CME": "C"}
        model_sequence = [
            residue.get_resname() for residue in model_chain.get_residues()
        ]
        native_sequence = [
            residue.get_resname() for residue in native_chain.get_residues()
        ]
        model_sequence = "".join(
            seq1(r, custom_map=custom_map)
            if len(r) == 3
            else r[:-1]
            if (len(r) == 2)
            else r
            for r in model_sequence
        )

        native_sequence = "".join(
            seq1(r, custom_map=custom_map)
            if len(r) == 3
            else r[:-1]
            if len(r) == 2
            else r
            for r in native_sequence
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


@lru_cache
def get_atoms_per_residue(
    chains,
    what,
    atom_types=("CA", "C", "N", "O", "P"),
):
    chain1, chain2 = chains
    atoms1 = [
        atom.coord
        for res1, res2 in zip(chain1, chain2)
        for atom in res1.get_atoms()
        if atom.id in atom_types and atom.id in [a.id for a in res2.get_atoms()]
    ]

    atoms2 = [
        atom.coord
        for res1, res2 in zip(chain1, chain2)
        for atom in res2.get_atoms()
        if atom.id in atom_types and atom.id in [a.id for a in res1.get_atoms()]
    ]
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
    no_align=False,
    capri_peptide=False,
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

    info = calc_DockQ(
        model_chains,
        native_chains,
        alignments=tuple(alignments),
        capri_peptide=capri_peptide,
        low_memory=False,
    )
    return info


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
                no_align=no_align,
                capri_peptide=capri_peptide,
                low_memory=False,
            )
            if info:
                info["chain1"], info["chain2"] = (
                    chain_map[chain_pair[0]],
                    chain_map[chain_pair[1]],
                )
                results_dic[chain_pair] = info
            

    return results_dic


# @profile
def load_PDB(path, chains=[], n_model=0):
    try:
        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure(
            "-",
            (gzip.open if path.endswith(".gz") else open)(path, "rt"),
            chains=chains,
        )
        model = structure[n_model]
    except Exception:
        pdb_parser = MMCIFParser(QUIET=True)
        structure = pdb_parser.get_structure(
            "-",
            (gzip.open if path.endswith(".gz") else open)(path, "rt"),
            chains=None,
        )
        model = structure[n_model]

    # remove_h(model)
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
        aln = align_chains(query_structure[query_chain], ref_structure[ref_chain])
        alignment = format_alignment(aln)
        n_mismatches = alignment["matches"].count(".")

        if n_mismatches > 0 and n_mismatches < 10:
            mismatch_dict[(query_chain, ref_chain)] = n_mismatches

        if n_mismatches <= allowed_mismatches:
            # 100% sequence identity, 100% coverage of native sequence in model sequence
            chain_clusters[ref_chain].append(query_chain)

    chains_without_match = [
        chain for chain in chain_clusters if not chain_clusters[chain]
    ]

    if chains_without_match:
        print(
            f"For these chains {chains_without_match} no match was found between model and native, try increasing the --allowed_mismatches from {allowed_mismatches}"
        )
        print(f"Current number of alignments with 1-10 mismatches: {mismatch_dict}")

    return chain_clusters, reverse_map


def format_mapping(mapping_str):
    mapping = dict()
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

    # mapping = sorted([(b, a) for a, b in chain_mapping.items()])
    # Sorting might change LRMSD since the definition of receptor/ligand for equal length depends on order
    mapping = [(b, a) for a, b in chain_mapping.items()]
    for (
        model_chain,
        native_chain,
    ) in mapping:  # sorted([(b,a) for a,b in chain_mapping.items()]):
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
    # result = set(list(map(lambda x: tuple(sorted(x)), result))) # to remove symmetric duplicates
    for prod in result:
        yield tuple(prod)

def count_chain_combinations(chain_clusters):
    counts = {}
    for chain in chain_clusters:
        chains = tuple(chain_clusters[chain])
        if chains not in counts:
            counts[chains]=0
        counts[chains]+=1
    number_of_combinations=np.prod([math.factorial(a) for a in counts.values()])
    return number_of_combinations
    #combos=itertools.product(*[itertools.permutations(chains) for chains in set([tuple(ch) for ch in chain_clusters.values()])])
    
    #return(number_of_combinations,counts)
    #set(chain_clusters.values())



def get_all_chain_maps(chain_clusters,initial_mapping,reverse_map,model_chains_to_combo,native_chains_to_combo):
    all_mappings = product_without_dupl(
        *[cluster for cluster in chain_clusters.values() if cluster]
    )
    for mapping in all_mappings:   
        chain_map = {key:value for key, value in initial_mapping.items()}
        if reverse_map:
            chain_map.update(
                {
                    mapping[i]: model_chain
                    for i, model_chain in enumerate(model_chains_to_combo)
                }
            )
        else:
            chain_map.update({
                native_chain: mapping[i] for i, native_chain in enumerate(native_chains_to_combo)
            })
        yield(chain_map)



# @profile
def main():
    args = parse_args()
    initial_mapping, model_chains, native_chains = format_mapping(args.mapping)
    model_structure = load_PDB(args.model, chains=model_chains)
    native_structure = load_PDB(args.native, chains=native_chains)

    # check user-given chains are in the structures
    model_chains = [c.id for c in model_structure] if not model_chains else model_chains
    native_chains = (
        [c.id for c in native_structure] if not native_chains else native_chains
    )
    info = {}
    if len(model_chains) < 2 or len(native_chains) < 2:
        print("Need at least two chains in the two inputs\n")
        sys.exit()

    # permute chains and run on a for loop
    best_dockq = -1
    best_result = None
    best_mapping = None

    model_chains_to_combo = [mc for mc in model_chains if mc not in initial_mapping.values()]
    native_chains_to_combo = [nc for nc in native_chains if nc not in initial_mapping.keys()]

    chain_clusters, reverse_map = group_chains(
        model_structure,
        native_structure,
        model_chains_to_combo,
        native_chains_to_combo,
        args.allowed_mismatches,
    )

    num_chain_combinations=count_chain_combinations(chain_clusters)
    chain_maps=get_all_chain_maps(chain_clusters,initial_mapping,reverse_map,model_chains_to_combo,native_chains_to_combo)
    
    low_memory=num_chain_combinations>100
    run_chain_map=partial(run_on_all_native_interfaces, 
                          model_structure, 
                          native_structure, 
                          no_align=args.no_align,  
                          capri_peptide=args.capri_peptide, 
                          low_memory=low_memory) ##args: chain_map
    
    if num_chain_combinations>1: 
        #chunk_size=max(1,num_chain_combinations // args.n_cpu)
        #I suspect large chunk_size will result in large input arguments to the workers.
        chunk_size=512
        #for large num_chain_combinations it should be possible to divide the chain_maps in chunks
        result_this_mappings=progress_map(run_chain_map,chain_maps, total=num_chain_combinations,n_cpu=args.n_cpu, chunk_size=chunk_size)
        #get a fresh iterator
        chain_maps=get_all_chain_maps(chain_clusters,initial_mapping,reverse_map,model_chains_to_combo,native_chains_to_combo)
        for chain_map,result_this_mapping in zip(chain_maps,result_this_mappings):
            total_dockq = sum(
                [result["DockQ_F1" if args.optDockQF1 else "DockQ"] for result in result_this_mapping.values()]
            )
         
            if total_dockq > best_dockq:
                best_dockq = total_dockq
                best_result = result_this_mapping
                best_mapping = chain_map

    else: #skip multi-threading for single jobs (skip the bar basically)
       # result_this_mappings=[run_chain_map(chain_map) for chain_map in chain_maps]
        for chain_map in chain_maps:
            result_this_mapping=run_chain_map(chain_map)
            total_dockq = sum(
                [result["DockQ_F1" if args.optDockQF1 else "DockQ"] for result in result_this_mapping.values()]
            )
            if total_dockq > best_dockq:
                best_dockq = total_dockq
                best_result = result_this_mapping
                best_mapping = chain_map
        
    
    
    if low_memory: #retrieve the full output by reruning the best chain mapping
        best_result=run_on_all_native_interfaces(
            model_structure,
            native_structure,
            best_mapping,
            args.no_align,
            args.capri_peptide,
            low_memory=False,
        )

    info["model"] = args.model
    info["native"] = args.native
    info["best_dockq"] = best_dockq
    info["best_result"] = best_result
    info["GlobalDockQ"] = best_dockq / len(best_result)
    info["best_mapping"] = best_mapping
    info["best_mapping_str"] = f"{format_mapping_string(best_mapping)}"
    print_results(info, args.short, args.verbose, args.capri_peptide)


def print_results(info, short=False, verbose=False, capri_peptide=False):
    if short:
        capri_peptide_str = "-capri_peptide" if capri_peptide else ""
        print(
            f"Total DockQ over {len(info['best_result'])} native interfaces: {info['GlobalDockQ']:.3f} with {info['best_mapping_str']} model:native mapping"
        )
        # print(info["best_result"])
        for chains, results in info["best_result"].items():
            print(
                f"DockQ{capri_peptide_str} {results['DockQ']:.3f} DockQ_F1 {results['DockQ_F1']:.3f} Fnat {results['fnat']:.3f} iRMS {results['irms']:.3f} LRMS {results['Lrms']:.3f} Fnonnat {results['fnonnat']:.3f} clashes {results['clashes']} mapping {results['chain1']}{results['chain2']}:{chains[0]}{chains[1]} {info['model']} {results['chain1']} {results['chain2']} -> {info['native']} {chains[0]} {chains[1]}"
            )
    else:
        print_header(verbose, capri_peptide)
        print(f"Model  : {info['model']}")
        print(f"Native : {info['native']}")
        print(
            f"Total DockQ over {len(info['best_result'])} native interfaces: {info['best_dockq']:.3f}"
        )
        items = ["DockQ_F1", "DockQ", "irms", "Lrms", "fnat"]

        for chains, results in info["best_result"].items():
            print(f"Native chains: {chains[0]}, {chains[1]}")
            print(f"\tModel chains: {results['chain1']}, {results['chain2']}")
            print("\n".join([f"\t{item}: {results[item]:.3f}" for item in items]))


def print_header(verbose=False, capri_peptide=False):
    reference = (
        "*   Ref: S. Basu and B. Wallner, DockQ: A quality measure for  *\n"
        "*   protein-protein docking models                             *\n"
        "*                            doi:10.1371/journal.pone.0161879  *\n"
        "*   For comments, please email: bjorn.wallner@.liu.se          *"
    )
    if not capri_peptide:
        header = (
            "****************************************************************\n"
            "*                       DockQ                                  *\n"
            "*   Scoring function for protein-protein docking models        *\n"
            "*   Statistics on CAPRI data:                                  *\n"
            "*    0.00 <= DockQ <  0.23 - Incorrect                         *\n"
            "*    0.23 <= DockQ <  0.49 - Acceptable quality                *\n"
            "*    0.49 <= DockQ <  0.80 - Medium quality                    *\n"
            "*            DockQ >= 0.80 - High quality                      *"
        )
    else:
        header = (
            "****************************************************************\n"
            "*                DockQ-CAPRI peptide                           *\n"
            "*   Do not trust any thing you read....                        *\n"
            "*   OBS THE DEFINITION OF Fnat and iRMS are different for      *\n"
            "*   peptides in CAPRI                                          *\n"
            "*                                                              *"
        )

    if verbose:
        notice = (
            "*   For the record:                                            *\n"
            f"*   Definition of contact <{'5A' if not capri_peptide else '4A'} (Fnat)                           *\n"
            f"*   Definition of interface <{'10A all heavy atoms (iRMS)       ' if not capri_peptide else '8A CB (iRMS)                     '} *\n"
            "****************************************************************"
        )
    else:
        notice = "****************************************************************"

    print(header)
    print(reference)
    print(notice)


if __name__ == "__main__":
    main()
