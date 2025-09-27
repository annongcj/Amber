#!/bin/bash
#TEST-PROGRAM pmemd.cuda
#TEST-DESCRIP Four replica GB Multi-D REMD.
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

#$1 = PREC_MODEL
#$2 = NETCDF

# Common setup. Sets PREC_MODEL, IG, and TESTsander if necessary
. ../CudaRemdCommon.sh
Title "Four replica GB Multi-dimensional REMD test."

# Remove any previous test output
CleanFiles remd.mdin? mdinfo.00? rem.out.00? mdcrd.00? restrt.00? groupfile \
           rem.log rem.type logfile.00? remd.dim

# Check that number of processors is appropriate
CheckNumProcessors 4

TOP1=ala4_ash.parm7
TOP2=ala4_asp.parm7

# Create MD input files
idx=1
for TEMP0 in "300.0" "350.0" ; do
  cat > remd.mdin$idx <<EOF
GB MREMD, $TEMP0 K
&cntrl
   imin = 0, nstlim = 10, dt = 0.001,
   ntx = 5, irest = 1, ig = $IG,
   ntwx = 500, ntwe = 0, ntwr = 1000, ntpr = 100,
   ntt = 1, tautp = 1.0, tempi = 0.0, temp0 = $TEMP0,
   ntc = 2, tol = 0.000001, ntf = 2, ntb = 0,
   cut = 9999.0, nscm = 500,
   igb = 5, offset = 0.09,
   numexchg = 10,
&end
EOF
  idx=$(( $idx + 1 ))
done
IN1=remd.mdin1
IN2=remd.mdin2

# Create groupfile
cat > groupfile <<EOF
# T=300K, ASH
-O -i $IN1 -p $TOP1 -c inpcrd.001 -o rem.out.001 -r restrt.001 -remlog rem.log
# T=350K, ASH
-O -i $IN2 -p $TOP1 -c inpcrd.002 -o rem.out.002 -r restrt.002 -remlog rem.log
# T=300K, ASP
-O -i $IN1 -p $TOP2 -c inpcrd.003 -o rem.out.003 -r restrt.003 -remlog rem.log
# T=350K, ASP
-O -i $IN2 -p $TOP2 -c inpcrd.004 -o rem.out.004 -r restrt.004 -remlog rem.log
EOF

# Create REMD dimensions file.
cat > remd.dim <<EOF
Temperature REMD
&multirem
   exch_type = 'TEMPERATURE',
   group(1,:) = 1,2,
   group(2,:) = 3,4,
   desc = 'Temperature exchange from 300K to 350K'
/
Hamiltonian REMD
&multirem
   exch_type='HAMILTONIAN',
   group(1,:) = 1,3,
   group(2,:) = 2,4,
   desc = 'Protonated ASP to Deprotonated ASP mutation'
/
EOF

$DO_PARALLEL $TESTsander -O -ng 4 -groupfile groupfile -remd-file remd.dim
CheckError $? "${0}"

DACDIF=../../../dacdif
DIFFOPTS=""
if [ "$PREC_MODEL" != "DPFP" ] ; then
  DIFFOPTS="-r 0.00004"
fi
$DACDIF rem.out.001.save rem.out.001 $DIFFOPTS
$DACDIF rem.out.002.save rem.out.002 $DIFFOPTS
$DACDIF rem.out.003.save rem.out.003 $DIFFOPTS
$DACDIF rem.out.004.save rem.out.004 $DIFFOPTS
$DACDIF rem.log.1.save rem.log.1 $DIFFOPTS
$DACDIF rem.log.2.save rem.log.2 $DIFFOPTS

# Cleanup
CleanFiles remd.mdin? mdinfo.00? mdcrd.00? restrt.00? groupfile rem.type \
           logfile.00? remd.dim
exit 0







