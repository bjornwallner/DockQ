#!/bin/bash
set -euo pipefail

which DockQ

# Test that cython version behaves the same as nocython
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb > test
diff test testdata/1A2K.dockq
coverage run -a --source='DockQ.DockQ' src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb > test
diff test testdata/1A2K.dockq
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align > test
diff test testdata/1A2K.dockq
coverage run -a --source='DockQ.DockQ' src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align > test
diff test testdata/1A2K.dockq

# Multiple interfaces
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/dimer_dimer.model.pdb examples/dimer_dimer.pdb  --short > test
diff test testdata/dimer_dimer.dockq

# Test on structures with slightly different sequences
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/model.pdb examples/native.pdb --allowed_mismatches 1 > test
diff test testdata/model.dockq

# Test various mapping strategies
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping AB*:BA* > test
diff test testdata/1EXB_AB.BA.dockq
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping :ABC > test
diff test testdata/1EXB_.ABC.dockq
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping ABCDEFGH:BADCFEHG > test
diff test testdata/1EXB_ABCDEFGH.BADCFEHG.dockq

# Test that cif parsing behaves same as pdb parsing
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --mapping DH:AE > test
diff test testdata/1EXB_DH.AE.dockq
coverage run -a --source='DockQ.DockQ' -m DockQ.DockQ examples/1EXB_r_l_b.model.pdb examples/1EXB.cif.gz --mapping DH:AE > test
diff test testdata/1EXB_DH.AE_cif.dockq

coverage report
