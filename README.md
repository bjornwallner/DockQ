# DockQ
Requires `Biopython` and `numpy` 

Installation
```
git clone https://github.com/bjornwallner/DockQ/
cd DockQ
make
```
Install `Biopython` and `numpy` 
- Biopython version >=1.64: http://biopython.org/wiki/Download#Installation_Instructions
- Numpy: http://www.scipy.org/install.html


Quick start for two interacting partners (two-chain-models) run with

`./DockQ.py <model> <native>`

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
vpn-83$ ./DockQ.py -h
usage: DockQ.py [-h] [-short] [-verbose] [-useCA] [-skip_check] [-perm1]
                [-perm2] [-chain1 chain1 [chain1 ...]]
                [-chain2 chain2 [chain2 ...]]
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
  -perm1                use all chain1 permutations to find maximum DockQ
                        (number of comparisons is n! = 24, if combined with
                        -perm2 there will be n!*m! combinations
  -perm2                use all chain2 permutations to find maximum DockQ
                        (number of comparisons is n! = 24, if combined with
                        -perm1 there will be n!*m! combinations
  -chain1 chain1 [chain1 ...]
                        pdb chain order to group together partner 1
  -chain2 chain2 [chain2 ...]
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
(`-perm1` and `-perm2`), this is important if for instance a partner is interacting
asymmetrically with a dimer. This should only be performed for
symmetric oligomers, where multiple solution are correct.

For this mode to work if there are missing residues the global
alignment program `needle` from the [EMBOSS
package](http://emboss.sourceforge.net/download/) needs to be in your
path.

To illustrate this mode we are using two dimers that are
interacting (A,B) -> (L H) that we are aligning to itself.

This command will put the chains A,B as one partner and the
remaining L H as the second partner. It will assume the chain
naming is the same in the model protein:

`./DockQ.py examples/dimer_dimer.pdb examples/dimer_dimer.model.pdb -native_chain1 A B`

this will be the same:

`./DockQ.py examples/dimer_dimer.pdb examples/dimer_dimer.model.pdb -native_chain1 A B -chain1 A B`

This will reverse the relative chain order of AB:

`./DockQ.py examples/dimer_dimer.pdb examples/dimer_dimer.model.pdb -native_chain1 A B -chain1 B A`

This will reverse the relative chain order of AB and LH

`./DockQ.py examples/dimer_dimer.pdb examples/dimer_dimer.model.pdb -native_chain1 A B -native_chain2 L H -chain1 B A -chain2 H L`

To try all permutations for chain1:

`./DockQ.py examples/dimer_dimer.pdb examples/dimer_dimer.model.pdb -native_chain1 A B -perm1`

To try all permutations for chain1 and chain2:

`./DockQ.py examples/dimer_dimer.pdb examples/dimer_dimer.model.pdb -native_chain1 A B -perm1 -perm2`

for two dimer this is only 4 (2!\*2!), however for a two tetramers
interacting the number will be 576 (4!\*4!)

24 combinations:

`./DockQ.py examples/tetramer_tetramer.pdb examples/tetramer_tetramer.pdb -native_chain1 A B C D -perm1`

576 combinations:

`./DockQ.py examples/tetramer_tetramer.pdb examples/tetramer_tetramer.pdb -native_chain1 A B C D -perm1 -perm2`


