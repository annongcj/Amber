#!/bin/sh
#TEST-PROGRAM pmemd
#TEST-DESCRIP ASM simulation of conformational change in alanine dipeptide
#TEST-PURPOSE Test and demonstrate the above
#TEST-STATE First version (Apr 2024)

. ../../program_error.sh

if [ -z "$TESTsander" ]; then
   TESTsander="${AMBERHOME}/bin/sander.MPI"
fi

if [ -z "$DO_PARALLEL" ]; then
   echo "This is a parallel test, DO_PARALLEL must be set!"
   exit 0
fi

if [ "`basename $TESTsander`" -ne "pmemd.MPI" ]; then
   echo "This test is for pmemd.MPI!"
   exit 0
fi

numprocs=`$DO_PARALLEL ../../numprocs`

if [ "$numprocs" -lt 24 ]; then
   echo "This test requires at least 24 processors!"
   exit 0
fi

numgroups=$((numprocs/2))

mkdir results
./in.sh $numprocs
cp STOP_STRING results/
$DO_PARALLEL $TESTsander -ng $numgroups -groupfile string.groupfile || error

/bin/rm -rf results [1-9]*
exit 0
