#!/bin/bash
set -euo pipefail

# Test that cython version behaves the same as nocython
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb
coverage run -a --source=src/ src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align
coverage run -a --source=src/ src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align
# Multiple interfaces
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/dimer_dimer.model.pdb examples/dimer_dimer.pdb  --short
# Test on structures with slightly different sequences
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/model.pdb examples/native.pdb --allowed_mismatches 1
# Test various mapping strategies
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping AB*:BA*
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping :ABC
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping ABCDEFGH:BADCFEHG
# Test that cif parsing behaves same as pdb parsing
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --mapping DH:AE
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB.cif.gz --mapping DH:AE

coverage report
