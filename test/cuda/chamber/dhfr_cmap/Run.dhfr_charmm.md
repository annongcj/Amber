#!/bin/csh -f
#TEST-PROGRAM pmemd.cuda
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

#$1 = PREC_MODEL
#$2 = NETCDF

if( ! $?DO_PARALLEL ) then
  setenv DO_PARALLEL " "
  if( $?TESTsander ) then
      set sander = $TESTsander
  else
      set sander = ${AMBERHOME}/bin/pmemd.cuda_$1
  endif
else
  if( $?TESTsander ) then
      set sander = $TESTsander
  else
      set sander = ${AMBERHOME}/bin/pmemd.cuda_$1.MPI
  endif
endif

cat > mdin <<EOF
 short md
 &cntrl
   ntx=1, irest=0,
   imin=0,nstlim=20,
   dt=0.002,ntc=2,ntf=2,
   ntt=1,tempi=300.0,temp0=300.0,
   ntpr=1,igb=1,cut=9999.0,ntwx=0,
   ntwr=0,ntwe=0,ntb=0, ig=71277,
 /
EOF

set output = mdout.dhfr_charmm_md

set output_save = $output.GPU_$1

touch dummy
$DO_PARALLEL $sander -O -i mdin -o $output <dummy || goto error
#Use different dacdif tolerances based on precision model.
if ( "$1" == "DPFP" ) then
  #7 sig figs
  ../../../dacdif -r 1.0e-7 $output_save $output
else
  #5 sig figs
  ../../../dacdif -r 1.0e-5 $output_save $output
endif

/bin/rm -f mdin logfile mdinfo dummy mdcrd restrt
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)







