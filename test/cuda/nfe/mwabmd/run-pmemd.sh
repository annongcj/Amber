#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../common.sh

#
# don't run if there is no netCDF
#

set_NC_UTILS

if test $? -ne 0; then
   exit 0;
fi

set_SANDER pmemd.cuda_$1.MPI
set_SANDER_CMD 4

JUNK="*umbrella.?.nc mdout.? restrt.? abmd.?.txt inpcrd.?"
JUNK="${JUNK} mdcrd.? groups prmtop sander.out *.log logfile.* *.ncdump"

#
# remove the junk
#

/bin/rm -rf ${JUNK} junk.*

#
# prepare files
#

cp -p ../prmtop ../groups ../inpcrd.? .

#
# run SANDER
#

${SANDER_CMD} -ng 4 -rem 0 -groupfile groups > sander.out 2>&1
CheckError $?

if [ "`basename $SANDER`" = "pmemd.cuda_DPFP" -o "`basename $SANDER`" = "pmemd.cuda_DPFP.MPI" ]; then
  for i in 1 2 3 4; do
    ../../../dacdif -r 1.0e-5 save.pmemd/abmd.$i.txt abmd.$i.txt
    ../../../dacdif -r 1.0e-4 save.pmemd/mdout.$i mdout.$i
  done

  for i in 1 2 3 4; do
    do_ncdump '%14.10f ' umbrella.$i.nc umbrella.$i.ncdump
    ../../../dacdif -r 1.0e-5 save.pmemd/umbrella.$i.ncdump umbrella.$i.ncdump
  done

  for i in 1 2 3 4; do
    do_ncdump '%14.10f ' wt_umbrella.$i.nc wt_umbrella.$i.ncdump
    ../../../dacdif -r 1.0e-5 save.pmemd/wt_umbrella.$i.ncdump wt_umbrella.$i.ncdump
  done
else
  for i in 1 2 3 4; do
    ../../../dacdif -r 1.0e-3 save.pmemd/abmd.$i.txt abmd.$i.txt
    ../../../dacdif -r 1.0e-2 save.pmemd/mdout.$i mdout.$i
  done

  for i in 1 2 3 4; do
    do_ncdump '%14.10f ' umbrella.$i.nc umbrella.$i.ncdump
    ../../../dacdif -r 1.0e-3 save.pmemd/umbrella.$i.ncdump umbrella.$i.ncdump
  done

  for i in 1 2 3 4; do
    do_ncdump '%14.10f ' wt_umbrella.$i.nc wt_umbrella.$i.ncdump
    ../../../dacdif -r 1.0e-3 save.pmemd/wt_umbrella.$i.ncdump wt_umbrella.$i.ncdump
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
