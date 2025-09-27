#!/bin/csh -f
#TEST-PROGRAM sander, pmemd
#TEST-DESCRIP microcanonical ensemble, periodic boundaries, shake, water
#TEST-PURPOSE regression, basic, pedagogical
#TEST-STATE   ok

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

cat > nvein << EOF
Constant energy dynamics of benzene
 &cntrl
  imin=0, ntx=5, irest=1,
  nstlim=10, dt=0.001,
  ntc=2, ntf=2,
  ntpr=1, ntwx=0, ntwr=0,
  cut=9.0, ntb=1, ntp=0, ntt=0, 
  nmropt=0, ioutfm=0, tol=0.0000001,
 /
EOF
$DO_PARALLEL $sander -O -i nvein -p ../../../virtual_sites/BenzenePi/prmtop_type8 \
                     -c ../../../virtual_sites/BenzenePi/eqrst -o nve.type8 || goto error
		     
grep "Etot" nve.type8 > etot.type8

$AMBERHOME/test/dacdif etot.type8_$1.save etot.type8

/bin/rm -f nvein restrt mdinfo

exit(0)

error:
echo "  ${0}:  Program error"
exit(1)
