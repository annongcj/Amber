#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

if test -z "${AMBERHOME}"; then
    export AMBERHOME="`pwd`/../../"
fi

if test -z "${TESTsander}"; then
    export TESTsander="${AMBERHOME}/bin/pmemd.cuda_$1.MPI"
fi

TESTS_2="abmd_ANALYSIS abmd_FLOODING abmd_UMBRELLA smd smd2"

# 4*n mpi threads
TESTS_4="abmd_ANALYSIS abmd_FLOODING abmd_UMBRELLA smd smd2 bbmd mwabmd premd abremd"

TESTS="${TESTS_2} ${TESTS_4}"

for t in ${TESTS_4} ; do
    echo ">>>>>>> doing '$t'"
    (cd $t ; ./run-pmemd.sh ; cd ..)
done
