#!/bin/bash
#TEST-PROGRAM pmemd.cuda
#TEST-DESCRIP Four replica GB HREMD.
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

#$1 = PREC_MODEL
#$2 = NETCDF

# Common setup. Sets PREC_MODEL, IG, and TESTsander if necessary
. ../CudaRemdCommon.sh
Title "Four replica GB H-REMD test."

# Remove any previous test output
CleanFiles mdin mdinfo.00? rem.out.00? mdcrd.00? restrt.00? groupfile \
           rem.log rem.type logfile.00?

# Check that number of processors is appropriate
CheckNumProcessors 4

cat > mdin <<EOF
GB HREMD
&cntrl
   imin = 0, nstlim = 100, dt = 0.002,
   ntx = 5, irest = 1, ig = $IG,
   ntwx = 500, ntwe = 0, ntwr = 1000, ntpr = 50,
   ioutfm = 0, ntxo = 1,
   ntt = 1, tautp = 5.0, tempi = 0.0, temp0 = $TEMP0,
   ntc = 2, tol = 0.000001, ntf = 2, ntb = 0,
   cut = 9999.0, nscm = 500,
   igb = 5, offset = 0.09,
   numexchg = 8,
&end
EOF

cat > groupfile <<EOF
# 1.00 ASH + 0.00 ASP
-O -p ala2_ash_1.0.parm7 -c ala2_ash.heat.rst7 -o rem.out -suffix 000
# 0.67 ASH + 0.33 ASP
-O -p ala2_ash_0.67.parm7 -c ala2_ash.heat.rst7 -o rem.out -suffix 001
# 0.33 ASH + 0.67 ASP
-O -p ala2_ash_0.33.parm7 -c ala2_ash.heat.rst7 -o rem.out -suffix 002
# 0.00 ASH + 1.00 ASP
-O -p ala2_ash_0.0.parm7 -c ala2_ash.heat.rst7 -o rem.out -suffix 003
EOF

$DO_PARALLEL $TESTsander -O -ng 4 -groupfile groupfile -rem 3 -remlog rem.log
CheckError $? "${0}"

DACDIF=../../../dacdif
DIFFOPTS=""
if [ "$PREC_MODEL" != "DPFP" ] ; then
  DIFFOPTS="-r 0.0003"
fi
$DACDIF rem.out.000.save rem.out.000 $DIFFOPTS
$DACDIF rem.out.001.save rem.out.001 $DIFFOPTS
$DACDIF rem.out.002.save rem.out.002 $DIFFOPTS
$DACDIF rem.out.003.save rem.out.003 $DIFFOPTS
$DACDIF rem.log.save rem.log $DIFFOPTS

# Cleanup
CleanFiles mdin mdinfo.00? mdcrd.00? restrt.00? groupfile rem.type logfile.00?
exit 0







