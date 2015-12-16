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
*** DockQ and other docking quality measures *** 
Number of equivalent residues in chain A 1492 (receptor)
Number of equivalent residues in chain B 912 (ligand)
Fnat 0.533333 32 correct of 60 native contacts
Fnonnat 0.238095 10 non-native of 42 model contacts
iRMS 1.30439307514
LRMS 1.51595045277
CAPRI Medium
DockQ_CAPRI Medium
DockQ 0.690639443815

