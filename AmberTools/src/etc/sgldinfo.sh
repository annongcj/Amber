#!/usr/bin/env csh
# This script print SGLD properties from sander simulation output
# csh
if ($#argv != 1 )then
  echo " Usage: sgldinfo.sh sander_output "
  exit
endif

set z = $argv[1]

awk -v sanderfile=$z  'BEGIN{t=-1;nt=0; sgfg=0;fsgldg=0;\
while(getline <sanderfile){if($1=="Local"&&$2=="averaging")tlavg=$4;\
if($1=="sgfti:"&&$3=="psgldg:"){sgft=$2;psgldg=$4;}\
if($1=="sgffi:"&&NF==2)sgff=$2;\
if($3=="sgfgi="&&$8=="fsgldg="){sgfg=$4;fsgldg=$9;}\
if($1=="NSTEP"){if(t!=$6){nt++;t=$6;}}} \
close(sanderfile);sta=nt/2;t=-1;nt=0;m=0;} \
{  if($1=="NSTEP"){skip=0;if(t==$6)skip=1; \
else {nt++;t=$6;temp=$9;if(nt<sta)skip=1;}} \
if(skip)next; \
if($1=="SGLD:"||$1=="SGMD:"){gamm=$2;templf=$3;temphf=$4;\
epotlf=$5;epothf=$6;epotllf=$7;wsg=$8; \
avgtemp+=temp;avggamm+=gamm;avgtlf+=templf;avgthf+=temphf;avgelf+=epotlf; avgehf+=epothf; \
avgellf+=epotllf;avgwsg+=wsg;m++;  \
}} \
END{if(m==0){printf(" This is not a SGMD/SGLD simulation: %-s\n",  \
sanderfile); exit;} \
avgtemp/=m;avggamm/=m;avgtlf/=m;avgthf/=m;avgelf/=m; avgehf/=m; \
avgellf/=m;avgwsg/=m;     \
printf("  SGLD properties from SANDER output: %-40s \n",sanderfile); \
printf("  local averaging time: %10.4f \n",tlavg); \
printf("  momentum guiding factor sgft, psgldg: %10.4f %10.4f\n",sgft,psgldg); \
printf("  force guiding factor sgff: %10.4f \n",sgff); \
printf("  SGLD-GLE momentum guiding factor sgfg, fsgldg: %10.4f %10.4f\n\n",\
sgfg,fsgldg); \
printf("   NPRINT         TIME       TEMP    TSGAVG   \n");    \
printf(" %8d   %10.2f   %8.2f  %8.4f \n\n",nt,t,avgtemp,tlavg); \
printf("      GAMM   TEMPLF   TEMPHF       EPOTLF   EPOTHF      EPOTLLF       WTSG\n"); \
printf("%10.4f %8.2f %8.2f %12.2f %8.2f %12.2f %10.4f\n",\
avggamm,avgtlf,avgthf,avgelf,avgehf,avgellf,avgwsg); \
}'   $z



