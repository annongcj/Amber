#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../common.sh

set_SANDER pmemd.cuda_$1.MPI
set_SANDER_CMD 4

JUNK="mdout.? restrt.? pmd.?.txt groups prmtop inpcrd.* sander.out *.log *.type logfile.*"

#
# remove the junk
#

/bin/rm -rf ${JUNK} junk.*

#
# prepare files
#

cp -p ../prmtop ../inpcrd.* ../groups .

#
# run SANDER
#

${SANDER_CMD} -ng 4 -rem 3 -groupfile groups > sander.out 2>&1
CheckError $?

if [ "`basename $SANDER`" = "pmemd.cuda_DPFP" -o "`basename $SANDER`" = "pmemd.cuda_DPFP.MPI" ]; then
  ../../../dacdif -r 1.0e-4 save/rem.log rem.log

  for i in 1 2 3 4; do
    ../../../dacdif -r 1.0e-5 save/pmd.$i.txt pmd.$i.txt
    ../../../dacdif -r 1.0e-4 save/mdout.$i mdout.$i
  done
else
  ../../../dacdif -r 5.0e-2 save/rem.log rem.log

  for i in 1 2 3 4; do
    ../../../dacdif -r 1.0e-3 save/pmd.$i.txt pmd.$i.txt
    ../../../dacdif -r 1.0e-2 save/mdout.$i mdout.$i
  done
fi

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK} profile_mpi
