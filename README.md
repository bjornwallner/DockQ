# DockQ
Requires `Biopython` and `numpy` 

Installation
```
git clone https://github.com/bjornwallner/DockQ/
cd DockQ
make
```
Install `Biopython` and `numpy` 
- Biopython: http://biopython.org/wiki/Download#Installation_Instructions
- Numpy: http://www.scipy.org/install.html




Run with
`./DockQ.py <model> <native>`

Example
```
bash$ ./DockQ.py model.pdb native.pdb
*** DockQ and other docking quality measures *** 
Number of equivalent residues in chain A 373 (receptor)
Number of equivalent residues in chain B 228 (ligand)
Fnat 0.533333 32 correct of 60
iRMS 1.15759298174
LRMS 1.50152873613
CAPRI Medium
DockQ 0.70993637399
