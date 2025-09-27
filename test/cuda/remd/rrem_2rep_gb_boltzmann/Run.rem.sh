#!/bin/bash
#TEST-PROGRAM pmemd.cuda
#TEST-DESCRIP Two replica GB Boltzmann Reservoir REMD.
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

#$1 = PREC_MODEL
#$2 = NETCDF

# Common setup. Sets PREC_MODEL, IG, and TESTsander if necessary
. ../CudaRemdCommon.sh
Title "Two replica GB Boltzmann Reservoir REMD test."

# Remove any previous test output
CleanFiles rem.in.00? mdinfo.00? rem.out.00? mdcrd.00? rst7.00? groupfile \
           rem.log rem.type logfile.00?

# Check that number of processors is appropriate
CheckNumProcessors 2

TOP=tc5b.top
CRD=tc5b.crd
i=0
for TEMP0 in "319.8" "340.9" ; do
  EXT=`printf "%03i" $i`
  cat > rem.in.$EXT <<EOF
Trp-cage GB REMD
&cntrl
   imin = 0, nstlim = 100, dt = 0.002,
   ntx = 5, irest = 1, ig = $IG,
   ntwx = 500, ntwe = 0, ntwr = 500, ntpr = 100,
   ioutfm = 1, 
   ntt = 1, tautp = 1, tol=0.000001, temp0 = $TEMP0,
   ntc = 2, ntf = 2, ntb = 0,
   cut = 9999.0, nscm = 500,
   igb = 8, 
   numexchg = 10,
&end
EOF
  echo "-O -i rem.in.$EXT -p $TOP -c $CRD -o rem.out.$EXT -x mdcrd.$EXT -r rst7.$EXT" >> groupfile
  i=$(( $i + 1 ))
done

$DO_PARALLEL $TESTsander -O -ng 2 -groupfile groupfile -rem 1 -remlog rem.log -rremd 1 -reservoir reservoir.nc
CheckError $? "${0}"

DACDIF=../../../dacdif
DIFFOPTS=""
if [ "$PREC_MODEL" != "DPFP" ] ; then
  DIFFOPTS="-r 0.00004"
fi
$DACDIF rem.out.000.save rem.out.000 $DIFFOPTS
$DACDIF rem.out.001.save rem.out.001 $DIFFOPTS
$DACDIF rem.log.save rem.log

# Cleanup
CleanFiles rem.in.00? mdinfo.00? mdcrd.00? rst7.00? groupfile rem.type logfile.00?
exit 0







