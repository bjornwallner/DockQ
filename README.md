# DockQ
Requires python packages: `numpy` and `Biopython`

Installation
```
git clone https://github.com/bjornwallner/DockQ/
cd DockQ
make
```
Install (i) `numpy` (a prerequisite to install 'Biopython') and (ii) `Biopython` 

- Numpy: http://www.scipy.org/install.html
- Biopython version >=1.64: http://biopython.org/wiki/Download#Installation_Instructions

Quick start for two interacting partners (two-chain-models) run with

`./DockQ.py <model> <native>`

To fix the residue numbering, in case there are inconsistencies, missing residues or small sequence differences between `model` and `native`

`scripts/fix_numbering.pl model.pdb native.pdb`

will output a file, `model.pdb.fixed`, with a residue numbering corresponding to the that in the `native.pdb` based on the sequences from the two pdb files.


Example

```

bash$ ./DockQ.py examples/model.pdb examples/native.pdb
***********************************************************
*                       DockQ                             *
*   Scoring function for protein-protein docking models   *
*   Statistics on CAPRI data:                             *
*    0    <  DockQ <  0.23 - Incorrect                    *
*    0.23 <= DockQ <  0.49 - Acceptable quality           *
*    0.49 <= DockQ <  0.80 - Medium quality               *
*            DockQ >= 0.80 - High quality                 *
*   Reference: Sankar Basu and Bjorn Wallner, DockQ:...   *
*   For comments, please email: bjornw@ifm.liu.se         *
***********************************************************

Number of equivalent residues in chain A 1492 (receptor)
Number of equivalent residues in chain B 912 (ligand)
Fnat 0.533 32 correct of 60 native contacts
Fnonnat 0.238 10 non-native of 42 model contacts
iRMS 1.232
LRMS 1.516
CAPRI Medium
DockQ_CAPRI Medium
DockQ 0.700

```

Help page
```
bash$ ./DockQ.py -h
usage: DockQ.py [-h] [-short] [-verbose] [-useCA] [-skip_check] [-no_needle]
                [-perm1] [-perm2]
                [-model_chain1 model_chain1 [model_chain1 ...]]
                [-model_chain2 model_chain2 [model_chain2 ...]]
                [-native_chain1 native_chain1 [native_chain1 ...]]
                [-native_chain2 native_chain2 [native_chain2 ...]]
                <model> <native>

DockQ - Quality measure for protein-protein docking models

positional arguments:
  <model>               path to model file
  <native>              path to native file

optional arguments:
  -h, --help            show this help message and exit
  -short                short output
  -verbose              talk a lot!
  -useCA                use CA instead of backbone
  -skip_check           skip initial check fo speed up on two chain examples
  -no_needle            do not use global alignment to fix residue numbering
                        between native and model during chain permutation (use
                        only in case needle is not installed, and the residues
                        between the chains are identical
  -perm1                use all chain1 permutations to find maximum DockQ
                        (number of comparisons is n! = 24, if combined with
                        -perm2 there will be n!*m! combinations
  -perm2                use all chain2 permutations to find maximum DockQ
                        (number of comparisons is n! = 24, if combined with
                        -perm1 there will be n!*m! combinations
  -model_chain1 model_chain1 [model_chain1 ...]
                        pdb chain order to group together partner 1
  -model_chain2 model_chain2 [model_chain2 ...]
                        pdb chain order to group together partner 2
                        (complement to partner 1 if undef)
  -native_chain1 native_chain1 [native_chain1 ...]
                        pdb chain order to group together from native partner
                        1
  -native_chain2 native_chain2 [native_chain2 ...]
                        pdb chain order to group together from native partner
                        2 (complement to partner 1 if undef)


					
```


##### Multi-chain functionality

For targets with more than two interacting chains. For instance a
dimer interacting with a partner. There are options to control which
chains to group together and also in which order to combine
them. There are also options to try all possible chain combinations
(`-perm1` and `-perm2`), this is important if for instance a homo
oligomer is interacting asymmetrically with a third partner, or if
there are symmetries that make multiple solution possibly correct.

For this mode to work if there are missing residues the global
alignment program `needle` from the [EMBOSS
package](http://emboss.sourceforge.net/download/) needs to be in your
path, i.e `which needle` should return the location.

This mode is illustrated by a homodimer that is interacting with a
third partner asymmetrically (A,B) <-> C (1A2K from docking benchmark
5.0).

The following commands will put the chains A,B as one partner and the
remaining chain, C, as the second partner. It will assume the chain
naming is the same in the model protein:

`./DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb -native_chain1 A B -model_chain1 A B -native_chain2 C -model_chain2 C`

Assuming the chains are the same in the model and native it is enough to just specify one set chains to group and the second group will be formed from the complement using the the remaining chains.

`./DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb -native_chain1 A B`
(chain C is remaining)

these are also equivalent:

`./DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb -native_chain1 C`
(chain AB is remaining)

`./DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb -native_chain1 A B -model_chain1 A B`

This will reverse the relative chain order of AB, comparing modelBA with nativeAB interacting with chain C:

`./DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb -native_chain1 A B -model_chain1 B A` (observe how the score increases)


To try all permutations for model_chain1, observe at the reverse
mapping BA -> AB gets the best score:

```
bash$ ./DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb -native_chain1 A B -perm1
1/2 AB -> C 0.00972962403319
Current best 0.00972962403319
2/2 BA -> C 0.476267208558
Current best 0.476267208558
Best score ( 0.476267208558 ) found for model -> native, chain1:BA -> AB chain2:C -> C
***********************************************************
*                       DockQ                             *
*   Scoring function for protein-protein docking models   *
*   Statistics on CAPRI data:                             *
*    0.00 <= DockQ <  0.23 - Incorrect                    *
*    0.23 <= DockQ <  0.49 - Acceptable quality           *
*    0.49 <= DockQ <  0.80 - Medium quality               *
*            DockQ >= 0.80 - High quality                 *
*   Reference: Sankar Basu and Bjorn Wallner, DockQ:...   *
*   For comments, please email: bjornw@ifm.liu.se         *
***********************************************************
Model  : examples/1A2K_r_l_b.model.pdb
Native : examples/1A2K_r_l_b.pdb
Best score ( 0.476267208558 ) found for model -> native, chain1:BA -> AB chain2:C -> C
Number of equivalent residues in chain A 248 (receptor)
Number of equivalent residues in chain B 196 (ligand)
Fnat 0.491 26 correct of 53 native contacts
Fnonnat 0.103 3 non-native of 29 model contacts
iRMS 1.988
LRMS 7.300
CAPRI Medium
DockQ_CAPRI Acceptable
DockQ 0.476


```


To try all permutations for model_chain1 and model_chain2 (ok only 1 chain in this example:-):

`./DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb -native_chain1 A B -perm1 -perm2`

For a dimer interacting with one partner this is only 2 (2!\*1!),
however for larger oligomers the number of possible permutations
explodes. For two tetramers the number will be 576 (4!\*4!)

Multimeric biological assemblies of proteins (of higher order than
that of dimers) are also found in nature to interact: e.g., PDB IDS:
3J07, 4IXZ, 4IY7, 4IYO, 3IYW, 2VQ0, 1E57, 2IZN etc. particularly
common in viral envelopes / capsids. We choose an instance of two
interacting tetramers (1EXB) to further demonstrate the multi-chain
functionality of DockQ

Tetramer example (24 combinations):

`./DockQ.py examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb -native_chain1 A B C D -perm1`

Tetramer example with all possible permutations (576 combinations):

`./DockQ.py examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb -native_chain1 A B C D -perm1 -perm2`


