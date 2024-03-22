![CI status](https://github.com/bjornwallner/DockQ/actions/workflows/main.yml/badge.svg)

# DockQ: A Quality Measure for Protein-Protein Docking Models

## Installation

Clone the repository, then install the necessary libraries with `pip`:

```
git clone https://github.com/bjornwallner/DockQ/
cd DockQ
pip install .
```

## Quick start:

After installing DockQ with `pip`, the `DockQ` binary will be in your path. Just run DockQ with:

`DockQ <model> <native>`

**Example**

When running DockQ on model/native complexes with one or more interfaces, you will get a result for each interface. Results are computed to maximise the average DockQ across all interface:

```

$ DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb

****************************************************************
*                       DockQ                                  *
*   Scoring function for protein-protein docking models        *
*   Statistics on CAPRI data:                                  *
*    0.00 <= DockQ <  0.23 - Incorrect                         *
*    0.23 <= DockQ <  0.49 - Acceptable quality                *
*    0.49 <= DockQ <  0.80 - Medium quality                    *
*            DockQ >= 0.80 - High quality                      *
*   Ref: S. Basu and B. Wallner, DockQ: A quality measure for  *
*   protein-protein docking models                             *
*                            doi:10.1371/journal.pone.0161879  *
*   For comments, please email: bjorn.wallner@.liu.se          *
****************************************************************
Model  : examples/1A2K_r_l_b.model.pdb
Native : examples/1A2K_r_l_b.pdb
Total DockQ over 3 native interfaces: 1.959
Native chains: A, B
        Model chains: B, A
        DockQ_F1: 0.996
        DockQ: 0.994
        irms: 0.000
        Lrms: 0.000
        fnat: 0.983
Native chains: A, C
        Model chains: B, C
        DockQ_F1: 0.567
        DockQ: 0.511
        irms: 1.237
        Lrms: 6.864
        fnat: 0.333
Native chains: B, C
        Model chains: A, C
        DockQ_F1: 0.500
        DockQ: 0.453
        irms: 2.104
        Lrms: 8.131
        fnat: 0.500
```

A more compact output option is available with the flag `--short`:

```
$ DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --short

DockQ 0.994 DockQ_F1 0.996 Fnat 0.983 iRMS 0.000 LRMS 0.000 Fnonnat 0.008 clashes 0 mapping BA:AB examples/1A2K_r_l_b.model.pdb B A -> examples/1A2K_r_l_b.pdb A B
DockQ 0.511 DockQ_F1 0.567 Fnat 0.333 iRMS 1.237 LRMS 6.864 Fnonnat 0.000 clashes 0 mapping BC:AC examples/1A2K_r_l_b.model.pdb B C -> examples/1A2K_r_l_b.pdb A C
DockQ 0.453 DockQ_F1 0.500 Fnat 0.500 iRMS 2.104 LRMS 8.131 Fnonnat 0.107 clashes 0 mapping AC:BC examples/1A2K_r_l_b.model.pdb A C -> examples/1A2K_r_l_b.pdb B C

```

**Other uses**

Run DockQ with `-h/--help` to see a list of the available flags:

```
bash$ DockQ -h

usage: DockQ [-h] [--capri_peptide] [--short] [--verbose] [--use_CA] [--no_align] [--optDockQF1] [--allowed_mismatches ALLOWED_MISMATCHES] [--mapping MODELCHAINS:NATIVECHAINS]
             <model> <native>

DockQ - Quality measure for protein-protein docking models

positional arguments:
  <model>               path to model file
  <native>              path to native file

optional arguments:
  -h, --help            show this help message and exit
  --capri_peptide       use version for capri_peptide (DockQ cannot not be trusted for this setting)
  --short               short output
  --verbose, -v         talk a lot!
  --use_CA, -ca         use CA instead of backbone
  --no_align            Do not align native and model using sequence alignments, but use the numbering of residues instead
  --optDockQF1          optimize on DockQ_F1 instead of DockQ
  --allowed_mismatches ALLOWED_MISMATCHES
                        number of allowed mismatches when mapping model sequence to native sequence.
  --mapping MODELCHAINS:NATIVECHAINS
                        Specify a chain mapping between model and native structure. If the native contains two chains "H" and "L" while the model contains two chains "A" and "B",
                        and chain A is a model of native chain H and chain B is a model of native chain L, the flag can be set as: '--mapping AB:HL'. This can also help limit the
                        search to specific native interfaces. For example, if the native is a tetramer (ABCD) but the user is only interested in the interface between chains B and
                        C, the flag can be set as: '--mapping :BC' or the equivalent '--mapping *:BC'.

```

