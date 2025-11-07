![CI status](https://github.com/bjornwallner/DockQ/actions/workflows/main.yml/badge.svg)

# DockQ
**A Quality Measure for Protein, Nucleic Acids and Small Molecule Docking Models**

## Installation

With `pip`:

```
pip install DockQ
```

Or, if you want the latest commit: clone the repository and run `pip` within the repo directory:

```
git clone https://github.com/bjornwallner/DockQ/
cd DockQ
pip install .
```

## Quick start

After installing DockQ with `pip`, the `DockQ` binary will be in your path. Just run DockQ with:

`DockQ <model> <native>`

**Example**

When running DockQ on model/native complexes with one or more interfaces, you will get a result for each interface. Results are computed to maximise the average DockQ across all interfaces:

```

$ DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb

****************************************************************
*                       DockQ                                  *
*   Docking scoring for biomolecular models                    *
*   DockQ score legend:                                        *
*    0.00 <= DockQ <  0.23 - Incorrect                         *
*    0.23 <= DockQ <  0.49 - Acceptable quality                *
*    0.49 <= DockQ <  0.80 - Medium quality                    *
*            DockQ >= 0.80 - High quality                      *
*   Ref: Mirabello and Wallner, 'DockQ v2: Improved automatic  *
*   quality measure for protein multimers, nucleic acids       *
*   and small molecules'                                       *
*                                                              *
*   For comments, please email: bjorn.wallner@.liu.se          *
****************************************************************
Model  : examples/1A2K_r_l_b.model.pdb
Native : examples/1A2K_r_l_b.pdb
Total DockQ over 3 native interfaces: 0.653 with BAC:ABC model:native mapping
Native chains: A, B
	Model chains: B, A
	DockQ: 0.994
	iRMSD: 0.000
	LRMSD: 0.000
	fnat: 0.983
	fnonnat: 0.008
	F1: 0.987
	clashes: 0
Native chains: A, C
	Model chains: B, C
	DockQ: 0.511
	iRMSD: 1.237
	LRMSD: 6.864
	fnat: 0.333
	fnonnat: 0.000
	F1: 0.500
	clashes: 0
Native chains: B, C
	Model chains: A, C
	DockQ: 0.453
	iRMSD: 2.104
	LRMSD: 8.131
	fnat: 0.500
	fnonnat: 0.107
    F1: 0.641
	clashes: 0
```

A more compact output option is available with the flag `--short`:

```
$ DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --short

Total DockQ over 3 native interfaces: 0.653 with BAC:ABC model:native mapping
DockQ 0.994 iRMSD 0.000 LRMSD 0.000 fnat 0.983 fnonnat 0.008 F1 0.987 clashes 0 mapping BA:AB examples/1A2K_r_l_b.model.pdb B A -> examples/1A2K_r_l_b.pdb A B
DockQ 0.511 iRMSD 1.237 LRMSD 6.864 fnat 0.333 fnonnat 0.000 F1 0.500 clashes 0 mapping BC:AC examples/1A2K_r_l_b.model.pdb B C -> examples/1A2K_r_l_b.pdb A C
DockQ 0.453 iRMSD 2.104 LRMSD 8.131 fnat 0.500 fnonnat 0.107 F1 0.641 clashes 0 mapping AC:BC examples/1A2K_r_l_b.model.pdb A C -> examples/1A2K_r_l_b.pdb B C
```

A dictionary containing the complete results data can be dumped in JSON format with the `--json filename.json` flag:

```
$ DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --json filename.json
```

## Model/Native chain mapping

By default, DockQ will try to find the optimal mapping between interfaces found in the native and in the model.

The simplest case is when a homodimer has been modelled. Then, the native interface "AB" (between native chains A and B) could be compared to
either the model AB interface, but also BA, as changing the order of the chains will generally change the results. If the user runs:

`DockQ homodimer_model.pdb homodimer_native.pdb`

the software will report the mapping (AB -> AB or AB -> BA) with highest DockQ score.

If the user wishes to enforce a certain mapping, the flag `--mapping` can be used. This is useful, for example, when model/native contain a large number of homologous chains,
as it will speed up computations.

**Complete mapping**

The user defines the complete mapping between native and model chains with: `--mapping MODELCHAINS:NATIVECHAINS`. For example, in the previous case, two possible mappings can be:

* `--mapping AB:AB` (native chain A corresponds to model chain A, native B to model B),
* `--mapping AB:BA` (native chain A corresponds to model chain B, native B to model A).

The pair before the colon `:` defines the chain order in the model, the pair after defines the order in the native.

**Partial mapping**

If the user wishes to fix part of the mapping and let DockQ optimize the rest, wildcards can be used. For example, if a tetramer has chains `ABCD` in the model and `WXYZ` in the native,
the user might use:

```
--mapping A*:W*
```

where the wildcard `*` indicates that DockQ should optimize the mapping between BCD and XYZ while keeping A -> W fixed. Multiple chains can be fixed:

```
--mapping AD*:WY*
```

**Limiting the search to subsets of native interfaces**

If the user is interested one or more specific interfaces in the native, while the rest should be ignored, the following can be used:

```
--mapping :WX
```

which is equivalent to:


```
--mapping *:WX
```

Then DockQ will find the interface in the model that best matches the WX interface in the native. Multiple native interfaces can be included:

```
--mapping :WXY
--mapping *:WXY
```

## Scoring small molecule docking poses

Small molecules in PDB or mmCIF files can be scored and the mapping optimized in the same way as for proteins. Just add the flag `--small_molecules`:

```
# Compare docking of hemoglobin chains (chain A and B in native) as well as HEM and PO4 groups (chains E, F, G)
$ DockQ examples/1HHO_hem.cif examples/2HHB_hem.cif --small_molecule --mapping :ABEFG --short

Total DockQ-small_molecules over 4 native interfaces: 0.659 with ABDCF:ABEFG model:native mapping
DockQ 0.950 iRMSD 0.455 LRMSD 1.451 fnat 0.964 fnonnat 0.070 clashes 0.000 F1 0.946 DockQ_F1 0.945 mapping AB:AB examples/1HHO_hem.cif A B -> examples/2HHB_hem.cif A B
LRMSD 0.585 mapping AD:AE (HEM) examples/1HHO_hem.cif A D -> examples/2HHB_hem.cif A E
LRMSD 28.096 mapping BC:BF (PO4) examples/1HHO_hem.cif B C -> examples/2HHB_hem.cif B F
LRMSD 1.311 mapping BF:BG (HEM) examples/1HHO_hem.cif B F -> examples/2HHB_hem.cif B G
```

Only LRMSD is reported for small molecules.

**NB: Small molecules must be in the PDB/mmCIF files. They also must have separate chain identifiers** (the `label_asym_id` field is used in mmCIF formatted files).

## Scoring DNA/RNA poses

Interfaces involving nucleic acids are seamlessly scored along with protein interfaces. The DockQ score is calculated for protein-NA or NA-NA interfaces in the same way as for protein-protein interfaces (two DockQ scores are reported for double helix chains).

## Other uses

Run DockQ with `-h/--help` to see a list of the available flags:

```
usage: DockQ [-h] [--capri_peptide] [--small_molecule] [--short] [--verbose] [--no_align] [--n_cpu CPU] [--max_chunk CHUNK] [--optDockQF1] [--allowed_mismatches ALLOWED_MISMATCHES] [--mapping MODELCHAINS:NATIVECHAINS] <model> <native>

DockQ - Quality measure for protein-protein docking models

positional arguments:
  <model>               Path to model file
  <native>              Path to native file

optional arguments:
  -h, --help            show this help message and exit
  --capri_peptide       use version for capri_peptide (DockQ cannot not be trusted for this setting)
  --small_molecule      If the docking pose of a small molecule should be evaluated
  --short               Short output
  --verbose, -v         Verbose output
  --no_align            Do not align native and model using sequence alignments, but use the numbering of residues instead
  --n_cpu CPU           Number of cores to use
  --max_chunk CHUNK     Maximum size of chunks given to the cores, actual chunksize is min(max_chunk,combos/cpus)
  --optDockQF1          Optimize on DockQ_F1 instead of DockQ
  --allowed_mismatches ALLOWED_MISMATCHES
                        Number of allowed mismatches when mapping model sequence to native sequence.
  --mapping MODELCHAINS:NATIVECHAINS
                        Specify a chain mapping between model and native structure. If the native contains two chains "H" and "L" while the model contains two chains "A" and "B", and chain A is a model of native chain H and chain B is a model
                        of native chain L, the flag can be set as: '--mapping AB:HL'. This can also help limit the search to specific native interfaces. For example, if the native is a tetramer (ABCD) but the user is only interested in the
                        interface between chains B and C, the flag can be set as: '--mapping :BC' or the equivalent '--mapping *:BC'.
```

**Importing DockQ as python module**

Once DockQ is installed with pip, it can also be used as a module in your python code (note that in this case the chain mapping has to be provided in the form native:model, not model:native):

```{python}
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

model = load_PDB("examples/1A2K_r_l_b.model.pdb")
native = load_PDB("examples/1A2K_r_l_b.pdb")

# native:model chain map dictionary for two interfaces
chain_map = {"A":"A", "B":"B"}
# returns a dictionary containing the results and the total DockQ score
run_on_all_native_interfaces(model, native, chain_map=chain_map)

({('A', 'B'): {'DockQ_F1': 0.9437927182141027,
   'DockQ': 0.9425398964102757,
   'iRMSD': 0.3753064373774967,
   'LRMSD': 0.5535111803522507,
   'fnat': 0.8907563025210085,
   'nat_correct': 106,
   'nat_total': 119,
   'fnonnat': 0.1016949152542373,
   'nonnat_count': 12,
   'model_total': 118,
   'clashes': 0,
   'len1': 124,
   'len2': 124,
   'class1': 'ligand',
   'class2': 'receptor',
   'chain1': 'A',
   'chain2': 'B',
   'chain_map': {'A': 'A', 'B': 'B'}}},
 0.9425398964102757)
```

**Merge multiple chains in a single receptor or ligand**

See issue [#33](https://github.com/bjornwallner/DockQ/issues/33) if you want to merge multiple chains (e.g. heavy and light antibody chains) into a single receptor or ligand

## Legend of outputs

* `DockQ`: DockQ score as defined in the [original DockQ paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0161879)
* `iRMSD`: RMSD of interfacial residues
* `LRMSD`: Ligand RMSD
* `fnat`: Fraction of retrieved native contacts (same as Recall or TPR)
* `fnonnat`: Fraction of predicted contacts that are not native (same as FPR)
* `F1`: Harmonic mean of the precision and recall in predicted interfacial contacts
* `clashes`: number of clashing interfacial residues in model (distance between residues below 2Å)

## Citing DockQ

If you use DockQ, please cite: [10.1093/bioinformatics/btae586](https://academic.oup.com/bioinformatics/article/40/10/btae586/7796530)
