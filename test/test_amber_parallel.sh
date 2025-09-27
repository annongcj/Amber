#!/bin/sh

# $1 represents the mode of execution

export AMBERHOME=$PMEMDHOME
. `dirname $0`/../test/test_check.sh
check_environment `dirname \`pwd\`` parallel
numprocs=`${DO_PARALLEL} ./numprocs`
if [ -z "$numprocs" ] ; then
    echo "Error: Invocation of numprocs with DO_PARALLEL did not produce output."
    echo "       Verify that DO_PARALLEL, which is \"$DO_PARALLEL\", is correct."
    echo "       Exiting."
    exit
fi
# Attempt to prune spurious output from mpi:
if [ `echo "$numprocs" |wc -w` -gt 1 ] ; then
    numprocs=`${DO_PARALLEL} ./numprocs | tail -1`
fi
if [ "$numprocs" -eq 1 ] ; then
    echo "Error: Single processor parallel jobs are not supported. Exiting."
    exit
fi

date_string=`date +%Y-%m-%d_%H-%M-%S`
logdir="${AMBERHOME}/logs/test_amber_parallel"
logprefix="${logdir}/${date_string}"
logfile="${logprefix}.log"
difffile="${logprefix}.diff"

mkdir -p ${logdir}
if [ $# -lt 1 ]; then
  (make -k -f Makefile test.parallel.pmemd 2>&1) | tee ${logfile}
else
  (make -k -f Makefile $1 2>&1) | tee -a ${logfile}
fi

(make -k -f Makefile finished.parallel 2>&1) | tee -a ${logfile}

passed_count=`grep --text PASS ${logfile} | wc -l`
questionable_count=`grep --text "FAILURE:" ${logfile} | wc -l`
questionable_ignored_count=`grep --text -e "FAILURE: (ignored)" ${logfile} | wc -l`
error_count=`grep --text "Program error" ${logfile} | wc -l`

echo "${passed_count} file comparisons passed" | tee -a ${logfile}
if [ ${questionable_count} -eq 0 ]; then
    echo "${questionable_count} file comparisons failed" | tee -a ${logfile}
else
    echo "${questionable_count} file comparisons failed (${questionable_ignored_count}" \
         "of which can be ignored)" | tee -a ${logfile}
fi
echo "${error_count} tests experienced an error" | tee -a ${logfile}

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
