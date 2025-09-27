# Provides a function to check that the Amber environment is sufficient before
# running tests that will ultimately fail.  This prevents spewing garbage and
# replaces it with a succinct error message.

# These environment variables are checked:  PMEMDHOME, DO_PARALLEL, and # TESTsander

# Usage: check_environment actual_Amber_home_path serial | parallel

check_environment() {
   # check_environment should be called with 2 arguments. The first should be
   # the value of PMEMDHOME. The second should be either "serial" or "parallel".
   # The third argyment is optional and tells us if the python tests are being
   # skipped or not.

   # First check that PMEMDHOME is set
   if [ -z "$PMEMDHOME" ]; then
      echo ""
      echo "Error: PMEMDHOME is not defined !"
      echo ""
      echo "Set PMEMDHOME to $1 and re-run the tests:"
      echo "The best way to do this is via one of the initialization scripts"
      echo "\$PMEMDHOME/amber.sh and \$PMEMDHOME/amber.csh which can be sourced"
      echo "now or from your shell resource file (e.g., ~/.bashrc or ~/.cshrc)."
      echo ""
      exit 1
   fi

   # For parallel simulations, check that DO_PARALLEL is not empty. For serial
   # simulations, empty DO_PARALLEL and warn if it was set
   if [ "$2" = "serial" ]; then
      if [ ! -z "$DO_PARALLEL" ]; then
         echo ""
         echo "Warning: DO_PARALLEL is set to \"$DO_PARALLEL\"."
         echo "This environment variable is being unset for serial testing."
         echo ""
         unset DO_PARALLEL
      fi
   fi

   if [ "$2" = "parallel" ]; then
      if [ -z "$DO_PARALLEL" ]; then
         echo ""
         echo "Error: DO_PARALLEL must be set for parallel tests!"
         echo ""
         exit 1
      fi
      echo ""
      echo "Tests being run with DO_PARALLEL=\"$DO_PARALLEL\"."
      echo ""
   fi

   # Unset TESTsander, since that's used by most tests
   if [ ! -z "$TESTsander" ]; then
      echo ""
      echo "Warning: TESTsander set to \"$TESTsander\"."
      echo "Unsetting TESTsander for tests."
      echo ""
      unset TESTsander
   fi

}
