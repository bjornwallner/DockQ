#!/bin/bash
set -euo pipefail

# Test that cython version behaves the same as nocython
DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb
python src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb
DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align
python src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align
# Multiple interfaces
DockQ examples/dimer_dimer.model.pdb examples/dimer_dimer.pdb  --short
# Test on structures with slightly different sequences
DockQ examples/model.pdb examples/native.pdb --allowed_mismatches 1
# Test various mapping strategies
DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping AB*:BA*
DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping :ABC
DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping ABCDEFGH:BADCFEHG
# Test that cif parsing behaves same as pdb parsing
DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --mapping DH:AE
DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB.cif.gz
