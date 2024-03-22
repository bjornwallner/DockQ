#!/bin/bash
set -euo pipefail

DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb
python src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb
DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align
python src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align
DockQ examples/model2.pdb examples/native2.pdb
DockQ examples/dimer_dimer.model.pdb examples/dimer_dimer.pdb  --short
DockQ examples/model.pdb examples/native.pdb --allowed_mismatches 1
DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping AB*:BA*
DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping :ABC
DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping ABCDEFGH:BADCFEHG
