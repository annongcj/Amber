# Source this script to define the environment variables necessary to use Amber.
# This script must be located in the Amber root folder!

# Amber was configured on @CONFIGURATION_TIMESTAMP@

set myname = 'amber.csh'

# Incomplete path coverage, but may catch novice user attempts:
if ( "$0" == "$myname" || "$0" == ./"$myname" ) then
        echo "Warning:  $myname is not a script; it should be sourced not executed!"
        echo "          Use it like this:  source $myname"
endif

if ( "$0" !~ '*csh' ) then
        echo "Warning:  $myname is a C shell source file!"
        echo "          Your shell does not appear to be a C shell:  $0"
endif

# Get path used for this source file (credit scott brozell).
set invocationpath = `echo $_ | cut -d' ' -f2- | sed "s@$myname.*@@"`
if ( "$invocationpath" == '' ) then
        set invocationpath = '.'
endif

setenv PMEMDHOME `cd "$invocationpath" >&! /dev/null; pwd`
setenv PATH "$PMEMDHOME/bin:@EXTRA_PATH_PART@$PATH"

