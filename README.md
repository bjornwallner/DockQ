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


Quick start
For two interacting partners (two-chain-models) run with
`./DockQ.py <model> <native>`

`vpn-77$ ./DockQ.py -h
usage: DockQ.py [-h] [-short] [-verbose] [-useCA] [-skip_check] [-perm]
                [-chain1 chain1 [chain1 ...]] [-chain2 chain2 [chain2 ...]]
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
  -perm                 use all chain permutations to find maximum DockQ
                        (number of comparisons is n!*m! = 24*24 = 576 for two
                        tetramers interacting)
  -chain1 chain1 [chain1 ...]
                        chains to group together partner 1
  -chain2 chain2 [chain2 ...]
                        chains to group together partner 2 (complement to
                        partner 1 if undef
  -native_chain1 native_chain1 [native_chain1 ...]
                        chains to group together from native partner 1
  -native_chain2 native_chain2 [native_chain2 ...]
                        chains to group together from native partner 2
                        (complement to partner 1 if undef)
			
`

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

Multi-chain functionality


