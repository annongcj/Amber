#!/bin/bash

# Testing MMPBSA_py package on mutations with igb = 66 and PB
# this unit testing has been added just for GBNSR6 testing, other options of MMPBSA_py are not implemented yet

OUTPUT_NAME=TEST_OUTPUT.dat

echo "================================================="
echo "Testing MMPBSA_py with igb = $1"

$AMBERHOME/bin/MMPBSA.py -O -i mmpbsa.in \
  -sp ras-raf_solvated.prmtop \
  -cp rasraf.prmtop \
  -rp ras.prmtop \
  -lp raf.prmtop \
  -y prod_compact.mdcrd \
  -mc rasraf_mutant.prmtop \
  -mr ras_mutant.prmtop \
  -o  $OUTPUT_NAME >> /dev/null

if test -f "$OUTPUT_NAME"; then
  echo "MMPBSA_py with igb = $1 passed successfully"
else
  echo "MMPBSA_py with igb = $1 failed!"
fi

echo "================================================="