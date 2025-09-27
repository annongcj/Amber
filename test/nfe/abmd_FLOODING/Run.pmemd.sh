#!/bin/bash
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

set_SANDER pmemd.MPI
set_SANDER_CMD 1

JUNK="mdout mdinfo restrt monitor.txt *umbrella.nc inpcrd prmtop logfile sander.out"

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

${SANDER_CMD} > sander.out 2>&1
CheckError $?

../../dacdif -t 1 save.pmemd/mdout mdout
../../dacdif -t 1 save.pmemd/monitor.txt monitor.txt

do_ncdump '%14.10f' save.pmemd/umbrella.nc umbrella.save.ncdump
do_ncdump '%14.10f' umbrella.nc umbrella.ncdump
../../dacdif umbrella.save.ncdump umbrella.ncdump

do_ncdump '%14.10f' save.pmemd/wt_umbrella.nc wt_umbrella.save.ncdump
do_ncdump '%14.10f' wt_umbrella.nc wt_umbrella.ncdump
../../dacdif wt_umbrella.save.ncdump wt_umbrella.ncdump

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK} *.ncdump profile_mpi
