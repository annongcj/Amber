#!/bin/bash -f

if [ -f groupfile ];
then
  rm remd.groupfile groupfile
fi

nrep=`wc temperatures.dat | awk '{print $1}'`
echo $nrep
count=0
for TEMP in `cat temperatures.dat`
do
  let COUNT+=1
  REP=`printf "%03d" $COUNT`
  echo "TEMPERATURE: $TEMP K ==> FILE: ala5.mdin.$REP"
  sed "s/XXXXX/$TEMP/g" ala5.mdin > temp
  sed "s/RANDOM_NUMBER/$RANDOM/g" temp > ala5.mdin.$REP
  cp ala5.inpcrd ala5.inpcrd.$REP
  echo "-O -rem 1 -remlog ala5.rem.log1 -i ala5.mdin.$REP -o ala5.mdout1.$REP -c ala5.inpcrd.$REP -r ala5.rst1.$REP -x ala5.mdcrd1.$REP -inf ala5.mdinfo1.$REP -p ala5.prmtop" >> remd.groupfile

  rm -f temp
done
echo "#" >> groupfile

echo "N REPLICAS  = $nrep"
echo " Done."
