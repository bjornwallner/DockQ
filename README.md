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




Run with
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
*    0.40 <= DockQ <  0.79 - Medium quality               *
*            DockQ >= 0.79 - High quality                 *
*   Reference: Sankar Basu and Bjorn Wallner, DockQ:...   *
*   For comments, please email: bjornw@ifm.liu.se         *
***********************************************************

Number of equivalent residues in chain A 1492 (receptor)
Number of equivalent residues in chain B 912 (ligand)
Fnat 0.533 32 correct of 60 native contacts
Fnonnat 0.238 10 non-native of 42 model contacts
iRMS 1.304
LRMS 1.516
CAPRI Medium
DockQ_CAPRI Medium
DockQ 0.691


