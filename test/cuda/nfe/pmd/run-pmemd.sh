#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../common.sh

set_SANDER pmemd.cuda_$1
set_SANDER_CMD 1

JUNK="mdout mdinfo restrt pmd.txt inpcrd prmtop sander.out logfile"

#
# remove the junk
#

/bin/rm -rf ${JUNK} junk.*

#
# prepare files
#

cp -p ../prmtop .
cp -p ../inpcrd.4 ./inpcrd

#
# run SANDER
#

${SANDER_CMD} -O -AllowSmallBox > sander.out
CheckError $?

if [ "`basename $SANDER`" = "pmemd.cuda_DPFP" -o "`basename $SANDER`" = "pmemd.cuda_DPFP.MPI" ]; then
  ../../../dacdif -r 1.0e-4 save.pmemd/mdout mdout
  ../../../dacdif -r 1.0e-5 save.pmemd/pmd.txt pmd.txt
else
  ../../../dacdif -r 1.0e-2 save.pmemd/mdout mdout
  ../../../dacdif -r 1.0e-3 save.pmemd/pmd.txt pmd.txt
fi

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK} profile_mpi
