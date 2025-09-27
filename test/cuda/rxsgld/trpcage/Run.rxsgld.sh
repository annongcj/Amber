#!/bin/bash
#TEST-PROGRAM pmemd.cuda
#TEST-DESCRIP Four replica GB REMD.
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

#$1 = PREC_MODEL
#$2 = NETCDF

# Common setup. Sets PREC_MODEL, IG, and TESTsander if necessary
. ../../remd/CudaRemdCommon.sh
Title "Four replica GB RXSGLD test."

# Remove any previous test output
CleanFiles rxsgld.in.00? mdinfo.00? rxsgld.out.00? mdcrd.00? rst7.00? groupfile \
           rxsgld.log rxsgld.type logfile.00?

# Check that number of processors is appropriate
CheckNumProcessors 4

TOP=prmtop
CRD=inpcrd
i=0
for SGFT in "0.0" "0.5" "0.8" "1.0" ; do
  EXT=`printf "%03i" $i`
  cat > rxsgld.in.$EXT <<EOF
TRPCAGE  GB RXSGLD
&cntrl
   imin = 0, nstlim = 10, dt = 0.002,
   ntx = 5, irest = 1, ig = $IG,
   ntwx = 100, ntwe = 0, ntwr = 500, ntpr = 10,
   ioutfm = 0, ntxo = 1,
   ntt = 3, tautp = 5.0, tempi = 0.0, temp0 = 300.0 ,
   ntc = 2, tol = 0.000001, ntf = 2, ntb = 0,
   cut = 9999.0, nscm = 10, gamma_ln=1.0,
   igb = 5, offset = 0.09,
   numexchg = 20,
   isgld=1,tsgavg=0.2,sgft=$SGFT, ig=71277,
&end
EOF
  echo "-O -i rxsgld.in.$EXT -p $TOP -c $CRD -o rxsgld.out.$EXT -x rxsgld.trj.$EXT -r rst7.$EXT" >> groupfile
  i=$(( $i + 1 ))
done

$DO_PARALLEL $TESTsander -O -ng 4 -groupfile groupfile -rem 1 -remlog rxsgld.log
CheckError $? "${0}"

DACDIF=../../../dacdif
DIFFOPTS=""
if [ "$PREC_MODEL" != "DPFP" ] ; then
  DIFFOPTS="-r 0.00004"
fi
$DACDIF rxsgld.out.000.GPU_$1 rxsgld.out.000 $DIFFOPTS
$DACDIF rxsgld.out.001.GPU_$1 rxsgld.out.001 $DIFFOPTS
$DACDIF rxsgld.out.002.GPU_$1 rxsgld.out.002 $DIFFOPTS
$DACDIF rxsgld.out.003.GPU_$1 rxsgld.out.003 $DIFFOPTS
$DACDIF rxsgld.log.GPU_$1 rxsgld.log

# Cleanup
CleanFiles rxsgld.in.00? mdinfo.00? mdcrd.00? rst7.00? groupfile rxsgld.type logfile.00? rxsgld.trj.00?

exit 0







