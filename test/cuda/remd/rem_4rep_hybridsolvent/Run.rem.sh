#!/bin/bash
#TEST-PROGRAM pmemd.cuda
#TEST-DESCRIP Four replica Hybrid Solvent REMD.
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

#$1 = PREC_MODEL
#$2 = NETCDF

# Common setup. Sets PREC_MODEL, IG, and TESTsander if necessary
. ../CudaRemdCommon.sh
Title "Four replica Hybrid Solvent REMD test."

# Remove any previous test output
CleanFiles rem.in.00? mdinfo.00? rem.out.00? mdcrd.00? rst7.00? groupfile \
           rem.log rem.type logfile.00?

# Check that number of processors is appropriate
CheckNumProcessors 4

TOP=hbsp.top
CRD=hbsp.crd
i=0
for TEMP0 in "287.0" "300.0" "313.5" "327.7" ; do
  EXT=`printf "%03i" $i`
  cat > rem.in.$EXT <<EOF
HBSP GB REMD
&cntrl
   imin = 0, nstlim = 100, dt = 0.004,
   ntx = 5, irest = 1, ig = $IG, ntxo = 2,
   ntwx = 500, ntwe = 0, ntwr = 500, ntpr = 100,
   ioutfm = 1, 
   ntt = 1, tautp = 1, tol = 0.000001, temp0 = $TEMP0,
   ntc = 2, ntf = 2, ntb = 1,
   cut = 8.0, nscm = 500, iwrap = 1,
   hybridgb = 8, numwatkeep = 150,
   numexchg = 4,
&end
EOF
  echo "-O -i rem.in.$EXT -p $TOP -c $CRD -o rem.out.$EXT -x mdcrd.$EXT -r rst7.$EXT" >> groupfile
  i=$(( $i + 1 ))
done

$DO_PARALLEL $TESTsander -O -ng 4 -groupfile groupfile -rem 1 -remlog rem.log
CheckError $? "${0}"

DACDIF=../../../dacdif
DIFFOPTS=""
if [ "$PREC_MODEL" != "DPFP" ] ; then
  DIFFOPTS="-r 0.00004"
fi
$DACDIF rem.out.000.save rem.out.000 $DIFFOPTS
$DACDIF rem.out.001.save rem.out.001 $DIFFOPTS
$DACDIF rem.out.002.save rem.out.002 $DIFFOPTS
$DACDIF rem.out.003.save rem.out.003 $DIFFOPTS
$DACDIF rem.log.save rem.log

# Cleanup
CleanFiles rem.in.00? mdinfo.00? mdcrd.00? rst7.00? groupfile rem.type logfile.00?
exit 0







