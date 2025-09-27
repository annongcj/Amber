#!/usr/bin/env csh
# This script calculate SGLD weight from sander simulation output
# csh
if ($#argv != 1 &&$#argv != 2 &&$#argv != 3 )then
  echo " Usage: sgldwt.sh sander_output  "
  exit
endif

set z = $argv[1]


awk -v sanderfile=$z  'BEGIN{t=-1;nt=0; sgfg=0;fsgldg=0;\
while(getline <sanderfile){if($1=="Local"&&$2=="averaging")tlavg=$4;\
if($1=="sgfti:"&&$3=="psgldg:"){sgft=$2;psgldg=$4;}\
if($1=="sgffi:"&&NF==2)sgff=$2;\
if($3=="sgfgi="&&$8=="fsgldg="){sgfg=$4;fsgldg=$9;}\
if($1=="NSTEP"){if(t!=$6){nt++;t=$6;}}} \
close(sanderfile);sta=nt/2;t=-1;nt=0;m=0; \
printf("      NOUT         TIME        WSG             WT   \n");    \
} \
{  if($1=="NSTEP"){skip=0;if(t==$6)skip=1; \
else {nt++;t=$6;temp=$9;if(nt<sta)skip=1;}} \
if($1=="SGLD:"||$1=="SGMD:"){gamm=$2;templf=$3;temphf=$4;\
epotlf=$5;epothf=$6;epotllf=$7;wsg=$8; \
wt=exp(wsg); printf("%10d %12.4f %10.4f %14.6f\n",nt,t,wsg,wt);\
if(skip)next; \
avgtemp+=temp;avggamm+=gamm;avgtlf+=templf;avgthf+=temphf;avgelf+=epotlf; avgehf+=epothf; \
avgellf+=epotllf;avgwsg+=wsg;m++;  \
wtgamm+=gamm*wt;wttlf+=templf*wt;wtthf+=temphf*wt;wtelf+=epotlf*wt; \
wtehf+=epothf*wt; wtellf+=epotllf*wt;wtwsg+=wt;  \
}} \
END{avgtemp/=m;avggamm/=m;avgtlf/=m;avgthf/=m;avgelf/=m; avgehf/=m; \
avgellf/=m;avgwsg/=m;     \
wtgamm/=wtwsg;wttlf/=wtwsg;wtthf/=wtwsg;wtelf/=wtwsg; wtehf/=wtwsg; \
wtellf/=wtwsg;     \
printf("\n  SGLD properties from SANDER output: %-40s \n",sanderfile); \
printf("  local averaging time: %10.4f \n",tlavg); \
printf("  momentum guiding factor sgft, psgldg: %10.4f %10.4f\n",sgft,psgldg); \
printf("  force guiding factor sgff: %10.4f \n",sgff); \
printf("  SGLD-GLE momentum guiding factor sgfg, fsgldg: %10.4f %10.4f\n\n",\
sgfg,fsgldg); \
printf("   NPRINT         TIME       TEMP    TSGAVG   \n");    \
printf(" %8d   %10.2f   %8.2f  %8.4f \n\n",nt,t,avgtemp,tlavg); \
printf("  Average properties before reweighting \n"); \
printf("      GAMM   TEMPLF   TEMPHF       EPOTLF   EPOTHF      EPOTLLF       WTSG\n"); \
printf("%10.4f %8.2f %8.2f %12.2f %8.2f %12.2f %10.4f\n",\
avggamm,avgtlf,avgthf,avgelf,avgehf,avgellf,avgwsg); \
printf("\n  Average properties after reweighting \n"); \
printf("      GAMM   TEMPLF   TEMPHF       EPOTLF   EPOTHF      EPOTLLF       WTSG\n"); \
printf("%10.4f %8.2f %8.2f %12.2f %8.2f %12.2f %10.4f\n",\
wtgamm,wttlf,wtthf,wtelf,wtehf,wtellf,wtwsg); \
}'   $z



