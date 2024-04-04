#!/bin/bash
set -euo pipefail

rm -f test .coverage

if command -v coverage &> /dev/null; then
    binary="coverage run -a -m DockQ.DockQ"
else
    binary="DockQ"
fi

$binary examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align > test
diff test testdata/1A2K.dockq
$binary examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb > test
diff test testdata/1A2K.dockq

# Multiple interfaces
$binary examples/dimer_dimer.model.pdb examples/dimer_dimer.pdb  --short > test
diff test testdata/dimer_dimer.dockq

# Test on structures with slightly different sequences
$binary examples/model.pdb examples/native.pdb --allowed_mismatches 1 > test
diff test testdata/model.dockq

# lowmem test
$binary examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short > test
diff test testdata/1EXB.dockq

# Test various mapping strategies
$binary examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping AB*:BA* > test
diff test testdata/1EXB_AB.BA.dockq
$binary examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping :ABC > test
diff test testdata/1EXB_.ABC.dockq
$binary examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping ABCDEFGH:BADCFEHG > test
diff test testdata/1EXB_ABCDEFGH.BADCFEHG.dockq

# Test that cif parsing behaves same as pdb parsing
$binary examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --mapping DH:AE > test
diff test testdata/1EXB_DH.AE.dockq
$binary examples/1EXB_r_l_b.model.pdb examples/1EXB.cif.gz --mapping DH:AE > test
diff test testdata/1EXB_DH.AE_cif.dockq

# Peptide measures
$binary examples/6qwn-assembly1.cif.gz examples/6qwn-assembly2.cif.gz --capri_peptide > test
diff test testdata/6q2n_peptide.dockq

# Test that cython version behaves the same as nocython
python src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb > test
diff test testdata/1A2K.dockq
python src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align > test
diff test testdata/1A2K.dockq
