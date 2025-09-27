#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP ASM simulation of SN2 reaction in water
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

numprocs=`$DO_PARALLEL ../../numprocs`

if [ "$numprocs" -lt 16 ]; then
   echo "This test requires at least 16 processors!"
   exit 0
fi

mkdir results
./in.sh $numprocs
cp STOP_STRING results/
$DO_PARALLEL $TESTsander -ng $numprocs -groupfile string.groupfile || error

/bin/rm -rf results [1-9]*
exit 0
