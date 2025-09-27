#!/bin/csh -f

if( ! $?TESTsander ) set TESTsander = "${AMBERHOME}/bin/pmemd.MPI"

if( ! $?DO_PARALLEL ) then
   echo "This test must be run in parallel - skipping."
   exit(0)
endif

$DO_PARALLEL $TESTsander \
    -i ti_ligands.in \
    -c ti.ncrst \
    -p ti.parm7 \
 -O -o ti001.out \
  -inf ti001.info \
    -e ti001.en \
    -r ti001.ncrst \
    -x ti001.nc \
    -l ti001.log

../dacdif -r 1.e-6 ti001.out.save ti001.out

/bin/rm -f ti001.info ti001.en ti001.ncrst ti001.nc ti001.log
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)

