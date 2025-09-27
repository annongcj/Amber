#include "copyright.i"
{
  unsigned int pos       = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  int atomsize=cSim.coarseMaxAtoms;
  int offset=3;
  const PMEFloat sixth = (PMEFloat)0.16666666666667;
  const PMEFloat twelvth = (PMEFloat)0.083333333333333;
#ifdef TI
  PMEFloat& tiWeight_ele1 = cSim.TIItemWeight[LambdaSchedule::TypeEleCC][0];
  PMEFloat& tiWeight_vdw1 = cSim.TIItemWeight[LambdaSchedule::TypeVDW][0];
  PMEFloat& tiWeight_ele2 = cSim.TIItemWeight[LambdaSchedule::TypeEleCC][1];
  PMEFloat& tiWeight_vdw2 = cSim.TIItemWeight[LambdaSchedule::TypeVDW][1];
#endif
  //Cycle through array that is allocated.
  while(pos < cSim.mcressize) 
  {
    double coarse_ene=0.0;
    int watatomid=cSim.pRandomWater[pos];
    double lj_ene=0.0;
    double ele_ene=0.0;
    PMEFloat deltax = cSim.pCoarseTrans[pos];
    PMEFloat deltay = cSim.pCoarseTrans[pos+cSim.mcressize];
    PMEFloat deltaz = cSim.pCoarseTrans[pos+2*cSim.mcressize];
    for(int i=0; i<cSim.waterAtomCount; i++)
    {
      int atm_i=watatomid+i;
      int img_i=cSim.pImageAtomLookup[atm_i];
      PMEFloat crd_x=cSim.pImageX[img_i];
      PMEFloat crd_y=cSim.pImageY[img_i];
      PMEFloat crd_z=cSim.pImageZ[img_i];
#ifdef IS_ORTHOG
      PMEFloat fx = cSim.recipf[0][0]*crd_x;
      PMEFloat fy = cSim.recipf[1][1]*crd_y;
      PMEFloat fz = cSim.recipf[2][2]*crd_z;
#else
      PMEFloat fx = cSim.recipf[0][0]*crd_x + cSim.recipf[1][0]*crd_y + cSim.recipf[2][0]*crd_z;
      PMEFloat fy =                           cSim.recipf[1][1]*crd_y + cSim.recipf[2][1]*crd_z;
      PMEFloat fz =                                                     cSim.recipf[2][2]*crd_z;
#endif
      fx = fx - round(fx);                     
      fy = fy - round(fy);
      fz = fz - round(fz);
      if(fx < 0.0) fx+=1.0;
      if(fy < 0.0) fy+=1.0;
      if(fz < 0.0) fz+=1.0;
      int gridx=fx*cSim.coarseMaxXvxl;
      int gridy=fy*cSim.coarseMaxYvxl;
      int gridz=fz*cSim.coarseMaxZvxl;
      for(int x=-1*offset+gridx; x<=gridx+offset; x++)
      {
        int bkt_x= x%cSim.coarseMaxXvxl;
        if(bkt_x < 0) bkt_x+=cSim.coarseMaxXvxl;
        for(int y=-1*offset+gridy; y<=offset+gridy; y++)
        {
          int bkt_y=y%cSim.coarseMaxYvxl;
          if(bkt_y < 0) bkt_y+=cSim.coarseMaxYvxl;
          for(int z=-1*offset+gridz; z<=offset+gridz; z++)
          {
            int bkt_z=z%cSim.coarseMaxZvxl;
            if(bkt_z < 0) bkt_z+=cSim.coarseMaxZvxl;
            int cur_voxel_flat=bkt_z * atomsize + bkt_y * (cSim.coarseMaxZvxl) * atomsize + bkt_x * (cSim.coarseMaxZvxl) * (cSim.coarseMaxYvxl) * atomsize;
            int bkt_atm_cnt=cSim.pCoarseGrid[0 + cur_voxel_flat];
            for(int atom=1; atom < bkt_atm_cnt+1; atom++)
            {
              int atm_j=cSim.pCoarseGrid[atom+cur_voxel_flat]-1;
              int img_j=cSim.pImageAtomLookup[atm_j];
              if(atm_j-atm_i+i >= cSim.waterAtomCount || atm_j-atm_i+i < 0) //Check we're not computing against self
              {
#ifdef IS_ORTHOG
                PMEFloat fxj = cSim.recipf[0][0]*cSim.pImageX[img_j];
                PMEFloat fyj = cSim.recipf[1][1]*cSim.pImageY[img_j];
                PMEFloat fzj = cSim.recipf[2][2]*cSim.pImageZ[img_j];
#else
                PMEFloat fxj = cSim.recipf[0][0]*cSim.pImageX[img_j] + cSim.recipf[1][0]*cSim.pImageY[img_j] + cSim.recipf[2][0]*cSim.pImageZ[img_j];
                PMEFloat fyj =                                         cSim.recipf[1][1]*cSim.pImageY[img_j] + cSim.recipf[2][1]*cSim.pImageZ[img_j];
                PMEFloat fzj =                                                                                 cSim.recipf[2][2]*cSim.pImageZ[img_j];
#endif
                fxj = fxj - round(fxj);
                fyj = fyj - round(fyj);
                fzj = fzj - round(fzj);
                fxj = (fxj < (PMEFloat)1.0 ? fxj : (PMEFloat)0.0);
                fyj = (fyj < (PMEFloat)1.0 ? fyj : (PMEFloat)0.0);
                fzj = (fzj < (PMEFloat)1.0 ? fzj : (PMEFloat)0.0);
                PMEFloat xi=fx-fxj-round(fx-fxj);
                PMEFloat yi=fy-fyj-round(fy-fyj);
                PMEFloat zi=fz-fzj-round(fz-fzj);
#ifdef IS_ORTHOG
                PMEFloat dx = xi*cSim.ucellf[0][0];
                PMEFloat dy = yi*cSim.ucellf[1][1];
                PMEFloat dz = zi*cSim.ucellf[2][2];
#else
                PMEFloat dx = xi*cSim.ucellf[0][0] + yi*cSim.ucellf[0][1] + zi*cSim.ucellf[0][2];
                PMEFloat dy =                        yi*cSim.ucellf[1][1] + zi*cSim.ucellf[1][2];
                PMEFloat dz =                                               zi*cSim.ucellf[2][2];
#endif
                PMEFloat delr2=dx*dx+dy*dy+dz*dz;             
                if(delr2 <= cSim.cut2)
                {
#ifdef TI
                  double timult_ele=1.0;
                  double timult_vdw=1.0;
                  if(cSim.pTIList[atm_j] > 0)
                  {
                    timult_ele = tiWeight_ele1;
                    timult_vdw = tiWeight_vdw1;
                    //printf("Pre R1: %i %i %i %f %f\n",atm_j,cSim.pTIList[atm_j],cSim.pTIList[atm_j+cSim.stride],timult_ele,timult_vdw);
                  }
                  if(cSim.pTIList[atm_j + cSim.stride] > 0)
                  {
                    timult_ele = tiWeight_ele2;
                    timult_vdw = tiWeight_vdw2;
                    //printf("Pre R2: %i %i %i %f %f\n",atm_j,cSim.pTIList[atm_j],cSim.pTIList[atm_j+cSim.stride],timult_ele,timult_vdw);
                  }
#endif
                  unsigned int index = (cSim.pImageLJID[img_i]) * cSim.LJTypes + cSim.pImageLJID[img_j];
                  lj_ene=0;
#ifndef use_DPFP
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
                  PMEFloat2 term = cSim.pLJTerm[index];
#  else
                  PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#  endif
#else
                  PMEFloat2 term = cSim.pLJTerm[index];
#endif
                  PMEFloat rinv  = rsqrt(delr2);
                  PMEFloat r = delr2 * rinv;
                  PMEFloat r2inv = rinv * rinv;
                  PMEFloat r6inv = r2inv * r2inv * r2inv;
                  PMEFloat f6 = term.y * r6inv;
                  PMEFloat f12 = term.x * r6inv * r6inv;
                  lj_ene = f12*twelvth - f6*sixth;
                  PMEFloat qiqj=cSim.pImageCharge[img_i]*cSim.pImageCharge[img_j];
#  ifdef use_DPFP
                  PMEFloat swtch = erfc(cSim.ew_coeffSP * r) * rinv;
#  else
                  PMEFloat swtch = fasterfc(r) * rinv;
#  endif
                  PMEFloat b0=qiqj*swtch; 
                  ele_ene = b0;
#ifdef TI
                  ele_ene *= timult_ele;
                  lj_ene *= timult_vdw;
#endif
                  coarse_ene-=lj_ene;
                  coarse_ene-=ele_ene;
                }
              }
            }
          }
        }
      }

      crd_x= crd_x - deltax;
      crd_y= crd_y - deltay;
      crd_z= crd_z - deltaz;
#ifdef IS_ORTHOG
      fx = cSim.recipf[0][0]*crd_x;
      fy = cSim.recipf[1][1]*crd_y;
      fz = cSim.recipf[2][2]*crd_z;
#else
      fx = cSim.recipf[0][0]*crd_x + cSim.recipf[1][0]*crd_y + cSim.recipf[2][0]*crd_z;
      fy =                           cSim.recipf[1][1]*crd_y + cSim.recipf[2][1]*crd_z;
      fz =                                                     cSim.recipf[2][2]*crd_z;
#endif
      fx = fx - round(fx);                     
      fy = fy - round(fy);
      fz = fz - round(fz);
      if(fx < 0.0) fx+=1.0;
      if(fy < 0.0) fy+=1.0;
      if(fz < 0.0) fz+=1.0;
      gridx=fx*cSim.coarseMaxXvxl;
      gridy=fy*cSim.coarseMaxYvxl;
      gridz=fz*cSim.coarseMaxZvxl;
      // Cycle through nearby grid atoms 
      for(int x=-1*offset+gridx; x<=gridx+offset; x++)
      {
        int bkt_x= x%cSim.coarseMaxXvxl;
        if(bkt_x < 0) bkt_x+=cSim.coarseMaxXvxl;
        for(int y=-1*offset+gridy; y<=offset+gridy; y++)
        {
          int bkt_y=y%cSim.coarseMaxYvxl;
          if(bkt_y < 0) bkt_y+=cSim.coarseMaxYvxl;
          for(int z=-1*offset+gridz; z<=offset+gridz; z++)
          {
            int bkt_z=z%cSim.coarseMaxZvxl;
            if(bkt_z < 0) bkt_z+=cSim.coarseMaxZvxl;
            int cur_voxel_flat=bkt_z * atomsize + bkt_y * (cSim.coarseMaxZvxl) * atomsize + bkt_x * (cSim.coarseMaxZvxl) * (cSim.coarseMaxYvxl) * atomsize;
            int bkt_atm_cnt=cSim.pCoarseGrid[0 + cur_voxel_flat];
            for(int atom=1; atom < bkt_atm_cnt+1; atom++)
            {
              int atm_j=cSim.pCoarseGrid[atom+cur_voxel_flat]-1;
              int img_j=cSim.pImageAtomLookup[atm_j];
              if(atm_j-atm_i+i >= cSim.waterAtomCount || atm_j-atm_i+i < 0) //Check we're not computing against self
              {
#ifdef IS_ORTHOG
                PMEFloat fxj = cSim.recipf[0][0]*cSim.pImageX[img_j];
                PMEFloat fyj = cSim.recipf[1][1]*cSim.pImageY[img_j];
                PMEFloat fzj = cSim.recipf[2][2]*cSim.pImageZ[img_j];
#else
                PMEFloat fxj = cSim.recipf[0][0]*cSim.pImageX[img_j] + cSim.recipf[1][0]*cSim.pImageY[img_j] + cSim.recipf[2][0]*cSim.pImageZ[img_j];
                PMEFloat fyj =                                         cSim.recipf[1][1]*cSim.pImageY[img_j] + cSim.recipf[2][1]*cSim.pImageZ[img_j];
                PMEFloat fzj =                                                                                 cSim.recipf[2][2]*cSim.pImageZ[img_j];
#endif
                fxj = fxj - round(fxj);
                fyj = fyj - round(fyj);
                fzj = fzj - round(fzj);
                fxj = (fxj < (PMEFloat)1.0 ? fxj : (PMEFloat)0.0);
                fyj = (fyj < (PMEFloat)1.0 ? fyj : (PMEFloat)0.0);
                fzj = (fzj < (PMEFloat)1.0 ? fzj : (PMEFloat)0.0);
                PMEFloat xi=fx-fxj-round(fx-fxj);
                PMEFloat yi=fy-fyj-round(fy-fyj);
                PMEFloat zi=fz-fzj-round(fz-fzj);
#ifdef IS_ORTHOG
                PMEFloat dx = xi*cSim.ucellf[0][0];
                PMEFloat dy = yi*cSim.ucellf[1][1];
                PMEFloat dz = zi*cSim.ucellf[2][2];
#else
                PMEFloat dx = xi*cSim.ucellf[0][0] + yi*cSim.ucellf[0][1] + zi*cSim.ucellf[0][2];
                PMEFloat dy =                        yi*cSim.ucellf[1][1] + zi*cSim.ucellf[1][2];
                PMEFloat dz =                                               zi*cSim.ucellf[2][2];
#endif
                PMEFloat delr2=dx*dx+dy*dy+dz*dz;             
                if(delr2 <= cSim.cut2)
                {
#ifdef TI
                  double timult_ele=1.0;
                  double timult_vdw=1.0;
                  if(cSim.pTIList[atm_j] > 0)
                  {
                    timult_ele = tiWeight_ele1;
                    timult_vdw = tiWeight_vdw1;
                    //printf("Post R1: %i %i %i %i %f %f\n",atm_j,cSim.pTIList[atm_j],cSim.pTIList[atm_j+cSim.stride],cSim.pSCList[atm_j],timult_ele,timult_vdw);
                  }
                  if(cSim.pTIList[atm_j + cSim.stride] > 0)
                  {
                    timult_ele = tiWeight_ele2;
                    timult_vdw = tiWeight_vdw2;
                    //printf("Post R2: %i %i %i %i %f %f\n",atm_j,cSim.pTIList[atm_j],cSim.pTIList[atm_j+cSim.stride],cSim.pSCList[atm_j],timult_ele,timult_vdw);
                  }
#endif
                  unsigned int index = (cSim.pImageLJID[img_i]) * cSim.LJTypes + cSim.pImageLJID[img_j];
#ifndef use_DPFP
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
                  PMEFloat2 term = cSim.pLJTerm[index];
#  else
                  PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#  endif
#else
                  PMEFloat2 term = cSim.pLJTerm[index];
#endif
                  PMEFloat rinv  = rsqrt(delr2);
                  PMEFloat r = delr2 * rinv;
                  PMEFloat r2inv = rinv * rinv;
                  PMEFloat r6inv = r2inv * r2inv * r2inv;
                  PMEFloat f6 = term.y * r6inv;
                  PMEFloat f12 = term.x * r6inv * r6inv;
                  lj_ene = f12*twelvth - f6*sixth;
                  PMEFloat qiqj=cSim.pImageCharge[img_i]*cSim.pImageCharge[img_j];
#  ifdef use_DPFP
                  PMEFloat swtch = erfc(cSim.ew_coeffSP * r) * rinv;
#  else
                  PMEFloat swtch = fasterfc(r) * rinv;
#  endif
                  PMEFloat b0=qiqj*swtch; 
                  ele_ene = b0;
#ifdef TI
                  ele_ene *= timult_ele;
                  lj_ene *= timult_vdw;
#endif
                  coarse_ene+=lj_ene+ele_ene;
                }
              }
            }
          }
        }
      }
    }
    cSim.pCoarseGridEne[pos] = coarse_ene;
    pos+=increment;
  }
}
