#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../common.sh

#
# Test to make sure there are at most six processors.
#
numprocs=`${DO_PARALLEL} echo "Testing number of processors" | wc -l`
if [ "$numprocs" -gt 2 ] ; then
	echo "This test requires 2 processors. Skipping."
	exit
fi


set_SANDER pmemd.cuda_$1
set_SANDER_CMD 1

JUNK="mdout mdinfo restrt work.txt sander.out logfile"

#
# remove the junk
#

/bin/rm -rf ${JUNK} junk.*

#
# run SANDER
#

${SANDER_CMD} > sander.out 2>&1
CheckError $?

if [ "`basename $SANDER`" = "pmemd.cuda_DPFP" -o "`basename $SANDER`" = "pmemd.cuda_DPFP.MPI" ]; then
  ../../../dacdif -r 1.0e-4 save.pmemd/mdout mdout
  ../../../dacdif -r 1.0e-5 save.pmemd/work.DPFP.txt work.txt
else
  ../../../dacdif -r 1.0e-2 save.pmemd/mdout mdout
  ../../../dacdif -r 1.0e-3 save.pmemd/work.DPFP.txt work.txt
fi

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK} profile_mpi
