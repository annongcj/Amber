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

JUNK="bbmd.log bbmd.?.txt umbrella.?.nc mdout.? inpcrd.? wt_umbrella.?.nc"
JUNK="${JUNK} mt19937.nc* restrt.? mdcrd.? groups prmtop sander.out logfile.* *.ncdump"

#
# remove the junk
#

/bin/rm -rf ${JUNK} junk.*

#
# prepare files
#

cp -p ../prmtop ../inpcrd.? ../groups ../mt19937.nc .

#
# run SANDER
#

${SANDER_CMD} -ng 4 -groupfile groups > sander.out 2>&1
CheckError $?

if [ "`basename $SANDER`" = "pmemd.cuda_DPFP" -o "`basename $SANDER`" = "pmemd.cuda_DPFP.MPI" ]; then
  ../../../dacdif -r 1.0e-5 save.pmemd/bbmd.log bbmd.log

  do_ncdump '%14.10f ' save.pmemd/mt19937.nc mt19937.save.ncdump
  do_ncdump '%14.10f ' mt19937.nc mt19937.ncdump
  ../../../dacdif -r 1.0e-4 mt19937.save.ncdump mt19937.ncdump

  for i in 1 2 3 4; do
    ../../../dacdif -f -r 1.0e-5 save.pmemd/bbmd.$i.txt bbmd.$i.txt
    ../../../dacdif -r 1.0e-4 save.pmemd/mdout.$i mdout.$i
  done

  for i in 1 2 3; do
    do_ncdump '%14.10f ' umbrella.$i.nc umbrella.$i.ncdump
    ../../../dacdif -f -r 1.0e-3 save.pmemd/umbrella.$i.ncdump umbrella.$i.ncdump
  done

  do_ncdump '%14.10f ' wt_umbrella.1.nc wt_umbrella.1.ncdump
  ../../../dacdif -r 1.0e-3 save.pmemd/wt_umbrella.1.ncdump wt_umbrella.1.ncdump
else
  ../../../dacdif -r 1.0e-3 save.pmemd/bbmd.log bbmd.log

  do_ncdump '%14.10f ' save.pmemd/mt19937.nc mt19937.save.ncdump
  do_ncdump '%14.10f ' mt19937.nc mt19937.ncdump
  ../../../dacdif -r 1.0e-3 mt19937.save.ncdump mt19937.ncdump

  for i in 1 2 3 4; do
    ../../../dacdif -r 1.0e-3 save.pmemd/bbmd.$i.txt bbmd.$i.txt
    ../../../dacdif -r 1.0e-2 save.pmemd/mdout.$i mdout.$i
  done

  for i in 1 2 3; do
    do_ncdump '%14.10f ' umbrella.$i.nc umbrella.$i.ncdump
    ../../../dacdif -r 1.0e-3 save.pmemd/umbrella.$i.ncdump umbrella.$i.ncdump
  done

  do_ncdump '%14.10f ' wt_umbrella.1.nc wt_umbrella.1.ncdump
  ../../../dacdif -r 1.0e-3 save.pmemd/wt_umbrella.1.ncdump wt_umbrella.1.ncdump
fi

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK}  profile_mpi
