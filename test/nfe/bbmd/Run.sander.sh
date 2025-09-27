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

set_SANDER sander.MPI
set_SANDER_CMD 4

JUNK="bbmd.log bbmd.?.txt umbrella.?.nc mdout.? inpcrd.? wt_umbrella.?.nc"
JUNK="${JUNK} mt19937.nc restrt.? mdcrd.? groups prmtop sander.out"

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

../../dacdif -t 1 save.sander/bbmd.log bbmd.log

../../dacdif save.sander/mt19937.nc mt19937.nc

for i in 1 2 3 4; do
    ../../dacdif -t 2 save.sander/bbmd.$i.txt bbmd.$i.txt
    ../../dacdif -t 1 save.sander/mdout.$i mdout.$i
done

for i in 1 2 3; do
   do_ncdump '%14.10f' save.sander/umbrella.$i.nc umbrella.$i.save.ncdump
   do_ncdump '%14.10f' umbrella.$i.nc umbrella.$i.ncdump
   ../../dacdif umbrella.$i.save.ncdump umbrella.$i.ncdump
done

do_ncdump '%14.10f ' save.sander/wt_umbrella.1.nc wt_umbrella.1.save.ncdump
do_ncdump '%14.10f ' wt_umbrella.1.nc wt_umbrella.1.ncdump
../../dacdif -r 1.0e-4 wt_umbrella.1.save.ncdump wt_umbrella.1.ncdump

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK} *.ncdump profile_mpi
