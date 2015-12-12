# DockQ
Requires `Biopython` and `numpy` 

Compile using `make`

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
