#!/bin/sh
#TEST-PROGRAM pmemd.cuda 
#TEST-DESCRIP Test replica exchange calc, NVT, should match CPU exactly. 
#TEST-PURPOSE regression, basic
#TEST-STATE   partially documented

#$1 = PREC_MODEL
#$2 = NETCDF

# Common setup. Sets PREC_MODEL, IG, and TESTsander if necessary
. ../CudaRemdCommon.sh
. ../../../program_error.sh

# Clean any previous test output.
CleanFiles rem.log rem.out.00? \
           rem.r* mdinfo* mdcrd* reminfo.* rem.type groupfile logfile.00? rem.in.00?

Title "Two replica PME NVT REMD test."

# Check that number of processors is appropriate
numprocs=`$DO_PARALLEL ../../../numprocs`
if [ "$numprocs" -ne 2 -a "$numprocs" -ne 4 ] ; then
  echo "This test requires 2 or 4 GPUs."
  exit 0
fi

log=rem.log
output=rem.out

TOP=../../../rem_wat/ala3.sol.top
CRD=../../../rem_wat/mdrestrt

CUT='5.5'
for REP in 000 001 ; do
  if [ "$REP" = '000' ] ; then
    TEMP0='300.0'
    P0='1.0'
  else
    TEMP0='304.0'
    P0='1.0'
  fi
        #ntb = 2, ntp = 1,
  cat > rem.in.$REP <<EOF
Title Line                
 &cntrl
        imin = 0, ntx = 5, nstlim = 100,
        ntc = 2, ntf = 2, tol=0.0000001, ntt = 1, dt = 0.002,
        ntb = 1, ntp = 0,
        barostat = 1, pres0 = $P0, taup = 1.0, irest = 1,
        ntwx = 100, ntwe = 0, ntwr = 100, ntpr = 50,
        cut = $CUT,
        ntr = 0, ibelly = 0, temp0 = $TEMP0, tempi = $TEMP0,
        nscm = 500, iwrap = 1,
        nsnb = 20,
        tautp = 0.1, offset = 0.09,
        numexchg = 2, 
        ntave = 0, ig=71277,
 &end
 &ewald
   nfft1 = 32, nfft2 = 32, nfft3 = 32,
 &end
EOF
  echo "-O -rem 1 -remlog $log -i rem.in.$REP -p $TOP -c $CRD -o ./$output.$REP -inf reminfo.$REP -r ./rem.r.$REP" >> groupfile
done

$DO_PARALLEL $TESTsander -O -ng 2 -groupfile groupfile < /dev/null || error

../../../dacdif $log.save $log
../../../dacdif rem.out.000.save rem.out.000
../../../dacdif rem.out.001.save rem.out.001

CleanFiles rem.r* mdinfo* mdcrd* reminfo.* rem.type groupfile logfile.00? rem.in.00?
exit 0
