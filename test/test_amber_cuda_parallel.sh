#!/bin/sh

#This script is for testing the cuda implmentation via the test
#Makefiles. It supports the use of 2 command line flags to control
#testing:
#
# $1 = PREC_MODEL  - default = DPFP
#
#Defaults
# PREC_MODEL = DPFP
#
#Example use:
#
# ./test_amber_cuda_parallel.x         (This will run with defaults)
# ./test_amber_cuda_parallel.x DPFP (This selects the hybrid DPFP
#                                    precision model to use)
#
# ./test_amber_cuda_parallel.x SPFP
#                            (This will run the hybrid SPFP precision model)
#                            This assumes each MPI thread is running on a
#                            separate node. There is currently no logic to
#                            be able to run the test cases on multiple GPUs
#                            within a single node without using the automatic
#                            selection algorithm. If you want to run a complex
#                            setup then you will need to run the test cases
#                            manually specifying the command line for each mpi
#                            task in your mpirun script. Check your MPI manual
#                            for details on how to accomplish this.

export AMBERHOME=$PMEMDHOME
. `dirname $0`/../test/test_check.sh
check_environment `dirname \`pwd\`` parallel

script_dir_relative=`dirname $0`
script_dir=`cd $script_dir_relative && pwd`

if [ $# -lt 1 ]; then
    echo "Using default PREC_MODEL = DPFP"
    PREC_MODEL="DPFP"
else
    PREC_MODEL=$1
fi

if [ ! -d $script_dir/cuda ]; then
    echo "CUDA tests not found; You probably only have AmberTools"
    exit 0
fi

if [ ! "${PREC_MODEL}" = "DPFP" ] && [ ! "${PREC_MODEL}" = "SPFP" ] ; then
    echo "Precision model options are: DPFP (double-precision, testing mode) "
    echo "and SPFP (single-precision hybrid model, production mode)."
    exit 0
fi

date_string=`date +%Y-%m-%d_%H-%M-%S`
logdir="${AMBERHOME}/logs/test_amber_cuda_parallel"
logprefix="${logdir}/${date_string}"
logfile="${logprefix}.log"
difffile="${logprefix}.diff"

mkdir -p ${logdir}
(make -k -f Makefile test.parallel.cuda PREC_MODEL=$PREC_MODEL 2>&1) | tee ${logfile}
(make -k -f Makefile finished.cuda.parallel 2>&1) | tee -a ${logfile}

passed_count=`grep --text PASS ${logfile} | wc -l`
questionable_count=`grep --text "FAILURE:" ${logfile} | wc -l`
questionable_ignored_count=`grep --text -e "FAILURE: (ignored)" ${logfile} | wc -l`
error_count=`grep --text "Program [Ee]rror" ${logfile} | wc -l`

echo "${passed_count} file comparisons passed" | tee -a ${logfile}
if [ ${questionable_count} -eq 0 ]; then
    echo "${questionable_count} file comparisons failed" | tee -a ${logfile}
else
    echo "${questionable_count} file comparisons failed (${questionable_ignored_count}" \
         "of which can be ignored)" | tee -a ${logfile}
fi
echo "${error_count} tests experienced errors" | tee -a ${logfile}

echo "Test log file saved as ${logfile}" | tee -a ${logfile}

if [ -f TEST_FAILURES.diff ]; then
   mv TEST_FAILURES.diff ${difffile}
   echo "Test diffs file saved as ${difffile}" | tee -a ${logfile}
else
   echo "No test diffs to save!" | tee -a ${logfile}
fi

if [ ${questionable_count} -gt ${questionable_ignored_count} -o ${error_count} -gt 0 ]; then
    # Make sure the exit code reflects the error count if we have comparison failures that
    # are not intended to be ignored
    exit 1
fi
