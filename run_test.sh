#!/bin/bash
set -euo pipefail

rm -f test .coverage

PYTHON=${1:-"python"}

if command -v coverage &> /dev/null; then
    binary="coverage run --parallel-mode -m DockQ.DockQ "
else
    binary="$PYTHON -m DockQ.DockQ"
fi

$binary examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb > test
diff <(grep -v "*" test) <(grep -v "*" testdata/1A2K.dockq)

$binary examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align > test
diff <(grep -v "*" test) <(grep -v "*" testdata/1A2K.dockq)

# Multiple interfaces
$binary examples/dimer_dimer.model.pdb examples/dimer_dimer.pdb --short > test
diff test testdata/dimer_dimer.dockq

# Test on structures with slightly different sequences
$binary examples/model.pdb examples/native.pdb --allowed_mismatches 1 > test
diff <(grep -v "*" test) <(grep -v "*" testdata/model.dockq)

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
diff <(grep -v "*" test) <(grep -v "*" testdata/1EXB_DH.AE.dockq)
$binary examples/1EXB_r_l_b.model.pdb examples/1EXB.cif.gz --mapping DH:AE > test
diff <(grep -v "*" test) <(grep -v "*" testdata/1EXB_DH.AE_cif.dockq)

# Peptide measures
$binary examples/6qwn-assembly1.cif.gz examples/6qwn-assembly2.cif.gz --capri_peptide > test
diff <(grep -v "*" test) <(grep -v "*" testdata/6q2n_peptide.dockq)

# Small molecule test
$binary examples/1HHO_hem.cif examples/2HHB_hem.cif --small_molecule --mapping :ABEFG > test

# Test that cython version behaves the same as nocython
$PYTHON src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb > test
diff <(grep -v "*" test) <(grep -v "*" testdata/1A2K.dockq)
$PYTHON src/DockQ/DockQ.py examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align > test
diff <(grep -v "*" test) <(grep -v "*" testdata/1A2K.dockq)

if command -v coverage &> /dev/null; then
	coverage combine
fi
