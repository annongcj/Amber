#!/bin/csh -f
#TEST-PROGRAM sander, pmemd
#TEST-DESCRIP microcanonical ensemble, periodic boundaries, shake, water
#TEST-PURPOSE regression, basic, pedagogical
#TEST-STATE   ok

if( $?TESTsander ) then
   set sander = $TESTsander
endif

if( ! $?DO_PARALLEL ) then
    setenv DO_PARALLEL " "
endif

cat > nvein << EOF
Constant volume and energy dynamics of tip4p
 &cntrl
  imin=0, ntx=5, irest=1,
  nstlim=10, dt=0.001,
  ntc=2, ntf=2,
  ntpr=1, ntwx=0, ntwr=0,
  cut=9.0, ntb=1, ntp=0, ntt=0, 
  nmropt=0, ioutfm=0, tol=0.0000001,
 /
EOF
$DO_PARALLEL $sander -O -i nvein -p prmtop_type5 -c eqrst -o nve.type5 \
                     -r nverst || goto error
		     
grep "Etot" nve.type5 > etot.type5
  
$AMBERHOME/test/dacdif etot.type5.save etot.type5

$DO_PARALLEL $sander -O -i nvein -p prmtop_type6 -c eqrst -o nve.type6 \
                     -r nverst || goto error

grep "Etot" nve.type6 > etot.type6

$AMBERHOME/test/dacdif etot.type6.save etot.type6

$DO_PARALLEL $sander -O -i nvein -p prmtop_type7 -c eqrst -o nve.type7 \
                     -r nverst || goto error

grep "Etot" nve.type7 > etot.type7

$AMBERHOME/test/dacdif etot.type7.save etot.type7

/bin/rm -f nvein restrt mdinfo nve.type5 nve.type6 nve.type7
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)
