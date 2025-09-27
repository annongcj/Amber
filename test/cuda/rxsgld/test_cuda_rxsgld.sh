#!/bin/bash

if [ -z "$DO_PARALLEL" ] ; then
  echo "CUDA RXSGLD tests require that DO_PARALLEL is set."
  exit 0
fi

if [ -z "$N_TEST_THREADS" ] ; then
  if [ ! -f "../../numprocs" ] ; then
    echo "Error: '../../numprocs' is missing."
    exit 1
  fi
  N_TEST_THREADS=`$DO_PARALLEL ../../numprocs`
  export N_TEST_THREADS
fi

# Default to DPFP if no precision model
if [ -z "$PREC_MODEL" ] ; then
  echo "No precision model specified. Defaulting to DPFP"
  PREC_MODEL="DPFP"
  export PREC_MODEL
fi

if [ $N_TEST_THREADS -eq 4 ] ; then
  make test.4rep 
else
  echo "CUDA RXSGLD tests can only be run with  4 threads ($N_TEST_THREADS currently)."
fi

exit 0
