# Provides functions common to CUDA REMD tests.
# Source this at the top of CUDA REMD test scripts.

# Common variables
PROGRAM_TYPE=0 # 0 = pmemd.cuda, 1 = pmemd.MPI, 2 = sander.MPI
IG=-71277      # Default random seed for GPU
CLEAN_ONLY=0   # If 1, exit with no error after CleanFiles called

# Environment variables
# PREC_MODEL     # GPU precision model
# N_TEST_THREADS # Current number of threads that will be used

# ------------------------------------------------------------------------------
# CleanFiles <file list>
CleanFiles() {
  while [ ! -z "$1" ] ; do
    if [ -f "$1" ] ; then
      /bin/rm -f $1
    fi
    shift
  done
  if [ $CLEAN_ONLY -eq 1 ] ; then
    exit 0
  fi
}

# ------------------------------------------------------------------------------
# Title <test title>
Title() {
  if [ $CLEAN_ONLY -eq 0 ] ; then
    echo "--------------------------------------------------------------------------------"
    echo "$1"
  fi
}

# ------------------------------------------------------------------------------
# CheckError <error value> <message>
CheckError() {
  if [ $1 -ne 0 ] ; then
    echo "  $2:  Program error"
    exit 1
  fi
}

# ------------------------------------------------------------------------------
# CheckNumProcessors <num groups>
CheckNumProcessors() {
  NGROUPS=$1
  # Get current number of threads
  if [ -z "$N_TEST_THREADS" ] ; then
    N_TEST_THREADS=`$DO_PARALLEL ../../../numprocs`
  fi
  if [ $PROGRAM_TYPE -eq 0 ] ; then
    # GPU. Must be exact.
    if [ $N_TEST_THREADS -ne $NGROUPS ] ; then
      echo "This test requires $NGROUPS threads ($N_TEST_THREADS currently)."
      exit 0
    fi
  else
    # pmemd.MPI requires at least two threads per group
    if [ $PROGRAM_TYPE -eq 1 ] ; then
      NGROUPS=$(( $NGROUPS * 2 ))
    fi
    if [ `expr $N_TEST_THREADS % $NGROUPS` -ne 0 ] ; then
      echo "This test requires a multiple of $NGROUPS threads ($N_TEST_THREADS currently)."
      exit 0
    fi
  fi
}

# ------------------------------------------------------------------------------
PREC_MODEL=$1
#$2 = NETCDF

# If 'clean' was specified instead of a precision model do not worry about setup
if [ "$PREC_MODEL" = 'clean' ] ; then
  CLEAN_ONLY=1
else
  if [ -z "$DO_PARALLEL" ] ; then
    echo "This test must be run with DO_PARALLEL set."
    exit 0
  fi

  # If the precision model has not been set, default to DPFP
  if [ -z "$PREC_MODEL" ] ; then
    echo "No precision model specified. Defaulting to DPFP."
    PREC_MODEL="DPFP"
  fi

  # Set up default binary if necessary
  if [ -z "$TESTsander" ] ; then
    TESTsander=${AMBERHOME}/bin/pmemd.cuda_$PREC_MODEL.MPI
  fi
  if [ ! -f "$TESTsander" ] ; then
    echo "  Binary $TESTsander missing:  Program error"
    exit 1
  fi

  # Set IG for CPU/GPU (debugging only)
  BASENAME=`basename $TESTsander`
  if [ -z `echo $BASENAME | grep cuda` ] ; then
    IG=71277
    if [ -z `echo $BASENAME | grep sander` ] ; then
      PROGRAM_TYPE=1
    else
      PROGRAM_TYPE=2
    fi
  fi
fi
