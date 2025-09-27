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
Constant volume and energy dynamics of bromobenzene
 &cntrl
  imin=0, ntx=5, irest=1,
  nstlim=10, dt=0.001,
  ntc=2, ntf=2,
  ntpr=1, ntwx=0, ntwr=0,
  cut=10.5, ntb=1, ntp=0, ntt=0, 
  nmropt=0, ioutfm=0, tol=0.0000001,
 &end
 &ewald
  nfft1 = 48, nfft2 = 48, nfft3 = 48,
 &end
EOF
$DO_PARALLEL $sander -O -i nvein -p ../../../virtual_sites/BromoBenzene/prmtop_type4 \
                     -c ../../../virtual_sites/BromoBenzene/eqrst -o nve.type4 || goto error

grep "Etot" nve.type4 > etot.type4

$AMBERHOME/test/dacdif etot.type4_$1.save etot.type4

$DO_PARALLEL $sander -O -i nvein -p ../../../virtual_sites/BromoBenzene/prmtop_type9 \
                     -c ../../../virtual_sites/BromoBenzene/eqrst -o nve.type9 || goto error
		     
grep "Etot" nve.type9 > etot.type9

$AMBERHOME/test/dacdif etot.type9_$1.save etot.type9

/bin/rm -f nvein nve.type4 nve.type9 mdinfo

exit(0)

error:
echo "  ${0}:  Program error"
exit(1)
