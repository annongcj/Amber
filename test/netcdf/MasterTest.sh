# All test scripts in subdirectories should source this

# CleanFiles(): For every arg passed to the function, check for the file and rm it
CleanFiles() {
  while [[ ! -z $1 ]] ; do
    #for RMFILE in `find . -name "$1"` ; do
    if [[ -e $1 ]] ; then
      #echo "  Cleaning $1"
      rm $1
    fi
    #done
    shift
  done
  # If only cleaning requested no run needed, exit now
  if [[ $CLEAN -eq 1 ]] ; then
    exit 0
  fi
}

SetCpptraj() {
  CPPTRAJ="$AMBERHOME/bin/cpptraj"
  if [ ! -e "$CPPTRAJ" ] ; then
    CPPTRAJ='../../../bin/cpptraj'
    if [ ! -e "$CPPTRAJ" ] ; then
      echo "Warning: This check requires cpptraj to perform NetCDF to ASCII conversion."
      return 1
    fi
  fi
  return 0
}

CheckError() {
  if [[ $1 -ne 0 ]] ; then
    echo "${0}: Program error"
    exit 1
  fi
}

#==============================================================================

# If the first argument is "clean" then no set-up is required. Script will
# exit when CleanFiles is called from sourcing script.
CLEAN=0
if [[ $1 = "clean" ]] ; then
  CLEAN=1
fi

# Get options, check environment
CPPTRAJ=""
if [[ $CLEAN -eq 0 ]] ; then
  # Check for TESTsander
  if [[ -z $TESTsander ]] ; then
    TESTsander="$AMBERHOME/bin/sander"
  fi
  # Check for DO_PARALLEL
  if [[ -z $DO_PARALLEL ]] ; then
    DO_PARALLEL=" "
  fi
fi

