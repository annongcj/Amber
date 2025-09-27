#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// Feb 2014, by Scott Le Grand and Ross C. Walker
//---------------------------------------------------------------------------------------------
{
#define NMR_TIMED

#ifdef NMR_ENERGY
  const int THREADS_PER_BLOCK = NMRFORCES_THREADS_PER_BLOCK;
  volatile __shared__ PMEAccumulator sE[THREADS_PER_BLOCK];
#  ifdef NMR_AFE
  volatile __shared__ PMEAccumulator sDVDL[THREADS_PER_BLOCK];
  volatile __shared__ PMEAccumulator sESCRestR1[THREADS_PER_BLOCK];
  volatile __shared__ PMEAccumulator sESCRestR2[THREADS_PER_BLOCK];
  __shared__  PMEDouble sLambda[3];

  // This is 1, 1-lambda, lambda
  sLambda[0] = (PMEDouble)1.0;
  sLambda[1] = cSim.AFElambda[0];
  sLambda[2] = cSim.AFElambda[1];

  //these are needed for dvdl
  __shared__  PMEDouble sTISigns[2];
  sTISigns[0] = -1.0;
  sTISigns[1] = 1.0;

#  endif
#endif
  int pos = blockIdx.x*blockDim.x + threadIdx.x + cSim.NMRTorsionOffset;

  // Calculate COM distance restraints
  if (pos < cSim.NMRCOMDistanceOffset) {
#ifdef NMR_ENERGY
    PMEAccumulator eCOMdistance = 0;
#ifdef NMR_AFE
    PMEAccumulator dvdl_COMdistance        = 0;
    PMEAccumulator scCOMdistance[2]            = {0, 0};
#endif
#endif
    pos -= cSim.NMRTorsionOffset;
    if (pos < cSim.NMRCOMDistances) {
#ifdef NMR_NEIGHBORLIST
      int2 atom = cSim.pImageNMRCOMDistanceID[pos];
#ifdef NMR_AFE
      int TIx = cSim.pImageTIRegion[atom.x];
      int TIy = cSim.pImageTIRegion[atom.y];
#endif
#else
      int2 atom = cSim.pNMRCOMDistanceID[pos];
#ifdef NMR_AFE
      int TIx = cSim.pTIRegion[atom.x];
      int TIy = cSim.pTIRegion[atom.y];
#endif
#endif

#ifdef NMR_AFE
      bool CVCOMdistance = ((((TIx >> 1) > 0) || ((TIy >> 1) > 0))
                           && (!((TIx | TIy) & 0x1)));
      int SCCOMdistance = ((TIx | TIy) & 0x1);
      int TICOMdistanceRegion = ((TIx | TIy) > 0) ? ((TIx | TIy) >> 1) : 0;
      PMEDouble CVlambda  = (CVCOMdistance) ? sLambda[TICOMdistanceRegion] : 1.0;
      PMEDouble lambda = (1 - SCCOMdistance) * CVlambda;
#endif

      PMEDouble2 rmstot;
      rmstot.x = 0.0;
      rmstot.y = 0.0;
      PMEDouble2 R1R2;
      PMEDouble2 R3R4;
      PMEDouble2 K2K3;
#ifdef NMR_TIMED
      int2 Step = cSim.pNMRCOMDistanceStep[pos];
      int Inc   = cSim.pNMRCOMDistanceInc[pos];

      // Skip restraint if not active
      if ((step < Step.x) || ((step > Step.y) && (Step.y > 0))) {
        goto exit4;
      }

      // Read timed data
      PMEDouble2 R1R2Slp = cSim.pNMRCOMDistanceR1R2Slp[pos];
      PMEDouble2 R1R2Int = cSim.pNMRCOMDistanceR1R2Int[pos];
      PMEDouble2 R3R4Slp = cSim.pNMRCOMDistanceR3R4Slp[pos];
      PMEDouble2 R3R4Int = cSim.pNMRCOMDistanceR3R4Int[pos];
      PMEDouble2 K2K3Slp = cSim.pNMRCOMDistanceK2K3Slp[pos];
      PMEDouble2 K2K3Int = cSim.pNMRCOMDistanceK2K3Int[pos];

      // Calculate increment
      double dstep = step - (double)((step - Step.x) % abs(Inc));

      // Calculate restraint values
      R1R2.x = dstep*R1R2Slp.x + R1R2Int.x;
      R1R2.y = dstep*R1R2Slp.y + R1R2Int.y;
      R3R4.x = dstep*R3R4Slp.x + R3R4Int.x;
      R3R4.y = dstep*R3R4Slp.y + R3R4Int.y;
      if (Inc > 0) {
        K2K3.x = dstep*K2K3Slp.x + K2K3Int.x;
        K2K3.y = dstep*K2K3Slp.y + K2K3Int.y;
      }
      else {
        int nstepu = (step - Step.x) / abs(Inc);
        K2K3.x     = K2K3Int.x * pow(K2K3Slp.x, nstepu);
        K2K3.y     = K2K3Int.y * pow(K2K3Slp.y, nstepu);
      }
#else
      R1R2 = cSim.pNMRCOMDistanceR1R2[pos];
      R3R4 = cSim.pNMRCOMDistanceR3R4[pos];
      K2K3 = cSim.pNMRCOMDistanceK2K3[pos];
#endif
      PMEDouble atomIX;
      PMEDouble atomIY;
      PMEDouble atomIZ;
      PMEDouble atomJX;
      PMEDouble atomJY;
      PMEDouble atomJZ;
#ifdef NODPTEXTURE
#  ifdef NMR_NEIGHBORLIST
      if (atom.x > 0) {
        atomIX = cSim.pImageX[atom.x];
        atomIY = cSim.pImageY[atom.x];
        atomIZ = cSim.pImageZ[atom.x];
      }
      else {

        // Find COM
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        int2 COMgrps   = cSim.pImageNMRCOMDistanceCOMGrp[pos * 2];
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
          int2 COMatms = cSim.pImageNMRCOMDistanceCOM[ip];
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass = cSim.pImageMass[j];
            xtot     = xtot + (cSim.pImageX[j] * rmass);
            ytot     = ytot + (cSim.pImageY[j] * rmass);
            ztot     = ztot + (cSim.pImageZ[j] * rmass);
            rmstot.x = rmstot.x + rmass;
          }
        }
       atomIX = xtot / rmstot.x;
       atomIY = ytot / rmstot.x;
       atomIZ = ztot / rmstot.x;
      }
      if (atom.y > 0) {
        atomJX = cSim.pImageX[atom.y];
        atomJY = cSim.pImageY[atom.y];
        atomJZ = cSim.pImageZ[atom.y];
      }
      else {

        // Find the center of mass
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot.y       = 0.0;
        int2 COMgrps   = cSim.pImageNMRCOMDistanceCOMGrp[pos * 2 + 1];
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
          int2 COMatms = cSim.pImageNMRCOMDistanceCOM[ip];
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass = cSim.pImageMass[j];
            xtot     = xtot + cSim.pImageX[j] * rmass;
            ytot     = ytot + cSim.pImageY[j] * rmass;
            ztot     = ztot + cSim.pImageZ[j] * rmass;
            rmstot.y = rmstot.y + rmass;
          }
        }
        atomJX = xtot / rmstot.y;
        atomJY = ytot / rmstot.y;
        atomJZ = ztot / rmstot.y;
      }
#  else  // NMR_NEIGHBORLIST
      if (atom.x > 0) {
        atomIX = cSim.pAtomX[atom.x];
        atomIY = cSim.pAtomY[atom.x];
        atomIZ = cSim.pAtomZ[atom.x];
      }
      else {

        // Find the center of mass
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot.x = 0.0;
        int2 COMgrps = cSim.pNMRCOMDistanceCOMGrp[pos * 2];
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
          int2 COMatms = cSim.pNMRCOMDistanceCOM[ip];
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass = cSim.pAtomMass[j];
            xtot     = xtot + cSim.pAtomX[j] * rmass;
            ytot     = ytot + cSim.pAtomY[j] * rmass;
            ztot     = ztot + cSim.pAtomZ[j] * rmass;
            rmstot.x = rmstot.x + rmass;
          }
        }
        atomIX = xtot / rmstot.x;
        atomIY = ytot / rmstot.x;
        atomIZ = ztot / rmstot.x;
      }
      if (atom.y > 0) {
        atomJX = cSim.pAtomX[atom.y];
        atomJY = cSim.pAtomY[atom.y];
        atomJZ = cSim.pAtomZ[atom.y];
      }
      else {

        // Find the center of mass
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot.y = 0.0;
        int2 COMgrps = cSim.pNMRCOMDistanceCOMGrp[pos * 2 + 1];
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
          int2 COMatms = cSim.pNMRCOMDistanceCOM[ip];
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass = cSim.pAtomMass[j];
            xtot     = xtot + cSim.pAtomX[j] * rmass;
            ytot     = ytot + cSim.pAtomY[j] * rmass;
            ztot     = ztot + cSim.pAtomZ[j] * rmass;
            rmstot.y = rmstot.y + rmass;
          }
        }
        atomJX = xtot / rmstot.y;
        atomJY = ytot / rmstot.y;
        atomJZ = ztot / rmstot.y;
      }
#  endif // NMR_NEIGHBORLIST
#else  // NODPTEXTURE
      if (atom.x > 0) {
#ifdef NMR_NEIGHBORLIST
        int2 iatomIX = tex1Dfetch<int2>(cSim.texImageX, atom.x);
        int2 iatomIY = tex1Dfetch<int2>(cSim.texImageX, atom.x + cSim.stride);
        int2 iatomIZ = tex1Dfetch<int2>(cSim.texImageX, atom.x + cSim.stride2);
#else
        int2 iatomIX = tex1Dfetch<int2>(cSim.texAtomX, atom.x);
        int2 iatomIY = tex1Dfetch<int2>(cSim.texAtomX, atom.x + cSim.stride);
        int2 iatomIZ = tex1Dfetch<int2>(cSim.texAtomX, atom.x + cSim.stride2);
#endif
        atomIX = __hiloint2double(iatomIX.y, iatomIX.x);
        atomIY = __hiloint2double(iatomIY.y, iatomIY.x);
        atomIZ = __hiloint2double(iatomIZ.y, iatomIZ.x);
      }
      else {

        // Find the center of mass
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot.x = 0.0;
#  ifdef NMR_NEIGHBORLIST
        int2 COMgrps = cSim.pImageNMRCOMDistanceCOMGrp[pos * 2];
#  else
        int2 COMgrps = cSim.pNMRCOMDistanceCOMGrp[pos * 2];
#  endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#  ifdef NMR_NEIGHBORLIST
          int2 COMatms = cSim.pImageNMRCOMDistanceCOM[ip];
#  else
          int2 COMatms = cSim.pNMRCOMDistanceCOM[ip];
#  endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
#  ifdef NMR_NEIGHBORLIST
            PMEDouble rmass = cSim.pImageMass[j];
#  else
            PMEDouble rmass = cSim.pAtomMass[j];
#  endif
#  ifdef NMR_NEIGHBORLIST
            int2 iatomIX = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomIY = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomIZ = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
#  else
            int2 iatomIX = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomIY = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomIZ = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
#  endif
            atomIX = __hiloint2double(iatomIX.y, iatomIX.x);
            atomIY = __hiloint2double(iatomIY.y, iatomIY.x);
            atomIZ = __hiloint2double(iatomIZ.y, iatomIZ.x);
            xtot     = xtot + atomIX * rmass;
            ytot     = ytot + atomIY * rmass;
            ztot     = ztot + atomIZ * rmass;
            rmstot.x = rmstot.x + rmass;
          }
        }
        atomIX = xtot / rmstot.x;
        atomIY = ytot / rmstot.x;
        atomIZ = ztot / rmstot.x;
      }
      if (atom.y > 0) {
#  ifdef NMR_NEIGHBORLIST
        int2 iatomJX = tex1Dfetch<int2>(cSim.texImageX, atom.y);
        int2 iatomJY = tex1Dfetch<int2>(cSim.texImageX, atom.y + cSim.stride);
        int2 iatomJZ = tex1Dfetch<int2>(cSim.texImageX, atom.y + cSim.stride2);
#  else
        int2 iatomJX = tex1Dfetch<int2>(cSim.texAtomX, atom.y);
        int2 iatomJY = tex1Dfetch<int2>(cSim.texAtomX, atom.y + cSim.stride);
        int2 iatomJZ = tex1Dfetch<int2>(cSim.texAtomX, atom.y + cSim.stride2);
#  endif
        atomJX = __hiloint2double(iatomJX.y, iatomJX.x);
        atomJY = __hiloint2double(iatomJY.y, iatomJY.x);
        atomJZ = __hiloint2double(iatomJZ.y, iatomJZ.x);
      }
      else {

        // Find the center of mass
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot.y       = 0.0;
#  ifdef NMR_NEIGHBORLIST
        int2 COMgrps = cSim.pImageNMRCOMDistanceCOMGrp[pos * 2 + 1];
#  else
        int2 COMgrps = cSim.pNMRCOMDistanceCOMGrp[pos * 2 + 1];
#  endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#  ifdef NMR_NEIGHBORLIST
          int2 COMatms = cSim.pImageNMRCOMDistanceCOM[ip];
#  else
          int2 COMatms = cSim.pNMRCOMDistanceCOM[ip];
#  endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++ ) {
#  ifdef NMR_NEIGHBORLIST
            PMEDouble rmass = cSim.pImageMass[j];
#  else
            PMEDouble rmass = cSim.pAtomMass[j];
#  endif
#  ifdef NMR_NEIGHBORLIST
            int2 iatomJX = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomJY = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomJZ = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
#  else
            int2 iatomJX = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomJY = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomJZ = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
#  endif
            atomJX = __hiloint2double(iatomJX.y, iatomJX.x);
            atomJY = __hiloint2double(iatomJY.y, iatomJY.x);
            atomJZ = __hiloint2double(iatomJZ.y, iatomJZ.x);
            xtot     = xtot + atomJX * rmass;
            ytot     = ytot + atomJY * rmass;
            ztot     = ztot + atomJZ * rmass;
            rmstot.y = rmstot.y + rmass;
          }
        }
        atomJX = xtot / rmstot.y;
        atomJY = ytot / rmstot.y;
        atomJZ = ztot / rmstot.y;
      }
#endif // NODPTEXTURE
      PMEDouble rij;
      PMEDouble xweight = cSim.pNMRCOMDistanceWeights[pos*5];
      PMEDouble yweight = cSim.pNMRCOMDistanceWeights[pos*5+1];
      PMEDouble zweight = cSim.pNMRCOMDistanceWeights[pos*5+2];

      // Direction ranges from 0 to 3
      // If Direction is 0 then 
      PMEDouble direction = cSim.pNMRCOMDistanceWeights[pos*5+3];

      // Initial neg or pos -1 or 1
      //PMEDouble posneg = cSim.pNMRCOMDistanceWeights[pos*5+4];
      PMEDouble xij = atomIX - atomJX;
      PMEDouble yij = atomIY - atomJY;
      PMEDouble zij = atomIZ - atomJZ;
      cSim.pNMRCOMDistanceXYZ[pos*3]  = xij;
      cSim.pNMRCOMDistanceXYZ[pos*3+1]= yij;
      cSim.pNMRCOMDistanceXYZ[pos*3+2]= zij;
      if (direction == 0) {
        rij = sqrt(xij*xij*xweight*xweight + yij*yij*yweight*yweight +
                   zij*zij*zweight*zweight);
      }
      else {
        //rij = (xij*xweight + yij*yweight + zij*zweight)*posneg;
        rij = (xij*xweight + yij*yweight + zij*zweight);
      }
      PMEDouble df;
#ifdef NMR_ENERGY
      PMEDouble e;
#endif
      if (rij < R1R2.x) {
        PMEDouble dif = R1R2.x - R1R2.y;
        df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
        e             = df*(rij - R1R2.x) + K2K3.x*dif*dif;
#endif
      }
      else if (rij < R1R2.y) {
        PMEDouble dif = rij - R1R2.y;
        df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
        e             = K2K3.x * dif * dif;
#endif
      }
      else if (rij < R3R4.x) {
        df            = (PMEDouble)0.0;
#ifdef NMR_ENERGY
        e             = (PMEDouble)0.0;
#endif
      }
      else if (rij < R3R4.y) {
        PMEDouble dif = rij - R3R4.x;
        df            = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
        e             = K2K3.y * dif * dif;
#endif
      }
      else {
        PMEDouble dif = R3R4.y - R3R4.x;
        df            = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
        e             = df*(rij - R3R4.y) + K2K3.y*dif*dif;
#endif
      }
      if (cSim.bJar) {
        double fold  = cSim.pNMRJarData[2];
        double work  = cSim.pNMRJarData[3];
        double first = cSim.pNMRJarData[4];
        double fcurr = (PMEDouble)-2.0 * K2K3.x * (rij - R1R2.y);
        if (first == (PMEDouble)0.0) {
          fold                = -fcurr;
          cSim.pNMRJarData[4] = (PMEDouble)1.0;
        }
        work                += (PMEDouble)0.5 * (fcurr + fold) * cSim.drjar;
        cSim.pNMRJarData[0]  = R1R2.y;
        cSim.pNMRJarData[1]  = rij;
        cSim.pNMRJarData[2]  = fcurr;
        cSim.pNMRJarData[3]  = work;
      }
#ifdef NMR_ENERGY
#  ifdef NMR_AFE
      scCOMdistance[TICOMdistanceRegion - SCCOMdistance] += ENERGYSCALE * e * SCCOMdistance;
      if (CVCOMdistance){
        dvdl_COMdistance += llrint(ENERGYSCALE * e * sTISigns[TICOMdistanceRegion - 1]);
      }
      e            *= lambda;
#  endif
      eCOMdistance                   += llrint(ENERGYSCALE * e);
#endif
      df            *= (PMEDouble)1.0 / rij;
      PMEDouble dfx  = df * xij * xweight * xweight;
      PMEDouble dfy  = df * yij * yweight * yweight;
      PMEDouble dfz  = df * zij * zweight * zweight;

      PMEAccumulator ifx;
      PMEAccumulator ify;
      PMEAccumulator ifz;
      if (atom.y > 0) {
        ifx = llrint(dfx * FORCESCALE);
        ify = llrint(dfy * FORCESCALE);
        ifz = llrint(dfz * FORCESCALE);
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.y],
                  llitoulli(ifx));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.y],
                  llitoulli(ify));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.y],
                  llitoulli(ifz));
      }
      else {
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps = cSim.pImageNMRCOMDistanceCOMGrp[pos * 2 + 1];
#else
        int2 COMgrps = cSim.pNMRCOMDistanceCOMGrp[pos * 2 + 1];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#ifdef NMR_NEIGHBORLIST
          int2 COMatms           = cSim.pImageNMRCOMDistanceCOM[ip];
#else
          int2 COMatms           = cSim.pNMRCOMDistanceCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
            PMEDouble dcomdx = rmass / rmstot.y;
            PMEDouble fx     = dfx * dcomdx;
            PMEDouble fy     = dfy * dcomdx;
            PMEDouble fz     = dfz * dcomdx;
            ifx = llrint(fx * FORCESCALE);
            ify = llrint(fy * FORCESCALE);
            ifz = llrint(fz * FORCESCALE);
            atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[j],
                      llitoulli(ifx));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[j],
                      llitoulli(ify));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[j],
                      llitoulli(ifz));
          }
        }
      }
      if (atom.x > 0) {
        ifx = llrint(dfx * FORCESCALE);
        ify = llrint(dfy * FORCESCALE);
        ifz = llrint(dfz * FORCESCALE);
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.x],
                  llitoulli(-ifx));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.x],
                  llitoulli(-ify));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.x],
                  llitoulli(-ifz));
      }
      else {
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps = cSim.pImageNMRCOMDistanceCOMGrp[pos * 2];
#else
        int2 COMgrps = cSim.pNMRCOMDistanceCOMGrp[pos * 2];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#ifdef NMR_NEIGHBORLIST
          int2 COMatms = cSim.pImageNMRCOMDistanceCOM[ip];
#else
          int2 COMatms = cSim.pNMRCOMDistanceCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++ ) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass     = cSim.pImageMass[j];
#else
            PMEDouble rmass     = cSim.pAtomMass[j];
#endif
            PMEDouble dcomdx = rmass / rmstot.x;
            PMEDouble fx     = dfx * dcomdx;
            PMEDouble fy     = dfy * dcomdx;
            PMEDouble fz     = dfz * dcomdx;
            ifx = llrint(fx * FORCESCALE);
            ify = llrint(fy * FORCESCALE);
            ifz = llrint(fz * FORCESCALE);
            atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[j],
                      llitoulli(-ifx));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[j],
                      llitoulli(-ify));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[j],
                      llitoulli(-ifz));
          }
        }
      }
    }
    exit4: {}
#ifdef NMR_ENERGY
    eCOMdistance += __SHFL(WARP_MASK, eCOMdistance, threadIdx.x ^ 1);
    eCOMdistance += __SHFL(WARP_MASK, eCOMdistance, threadIdx.x ^ 2);
    eCOMdistance += __SHFL(WARP_MASK, eCOMdistance, threadIdx.x ^ 4);
    eCOMdistance += __SHFL(WARP_MASK, eCOMdistance, threadIdx.x ^ 8);
    eCOMdistance += __SHFL(WARP_MASK, eCOMdistance, threadIdx.x ^ 16);
#ifdef AMBER_PLATFORM_AMD_WARP64
    eCOMdistance += __SHFL(WARP_MASK, eCOMdistance, threadIdx.x ^ 32);
#endif
    // Write out energies
    if ((threadIdx.x & GRID_BITS_MASK) == 0) {
      __syncthreads();
      atomicAdd(cSim.pENMRCOMDistance, llitoulli(eCOMdistance)); 
    }
#ifdef NMR_AFE
    sDVDL[threadIdx.x]  = dvdl_COMdistance;
    sESCRestR1[threadIdx.x]  = scCOMdistance[0];
    sESCRestR2[threadIdx.x]  = scCOMdistance[1];
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 1];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 1];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 1];
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 2];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 2];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 2];
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 4];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 4];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 4];
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 8];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 8];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 8];
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 16];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 16];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 16];
#ifdef AMBER_PLATFORM_AMD_WARP64
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 32];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 32];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 32];
#endif
    if ((threadIdx.x & GRID_BITS_MASK) == 0) {
      atomicAdd(cSim.pDVDL, llitoulli(sDVDL[threadIdx.x]));
      atomicAdd(cSim.pESCNMRCOMDistanceR1, llitoulli(sESCRestR1[threadIdx.x]));
      atomicAdd(cSim.pESCNMRCOMDistanceR2, llitoulli(sESCRestR2[threadIdx.x]));
    }
#endif
#endif
  }


#ifdef R6AV
  // Calculate r6av distance restraints
  else if (pos >= cSim.NMRCOMDistanceOffset && pos < cSim.NMRr6avDistanceOffset) {
#ifdef NMR_ENERGY
    PMEAccumulator er6avdistance = 0;
#ifdef NMR_AFE
    PMEAccumulator dvdl_r6avdistance  = 0;
    PMEAccumulator scr6avdistance[2]    = {0, 0};
#endif
#endif
    pos -= cSim.NMRCOMDistanceOffset;
    if (pos < cSim.NMRr6avDistances) {
#ifdef NMR_NEIGHBORLIST
      int2 atom = cSim.pImageNMRr6avDistanceID[pos];
#ifdef NMR_AFE
      int TIx = cSim.pImageTIRegion[atom.x];
      int TIy = cSim.pImageTIRegion[atom.y];
#endif
#else
      int2 atom = cSim.pNMRr6avDistanceID[pos];
#ifdef NMR_AFE
      int TIx = cSim.pTIRegion[atom.x];
      int TIy = cSim.pTIRegion[atom.y];
#endif
#endif

#ifdef NMR_AFE
      bool CVr6avdistance = ((((TIx >> 1) > 0) || ((TIy >> 1) > 0))
                             && (!((TIx | TIy) & 0x1)));
      int SCr6avdistance = ((TIx | TIy) & 0x1);
      int TIr6avdistanceRegion = ((TIx | TIy) > 0) ? ((TIx | TIy) >> 1) : 0;
      PMEDouble CVlambda  = (CVr6avdistance)
                            ? sLambda[TIr6avdistanceRegion] : 1.0;
      PMEDouble lambda = (1 - SCr6avdistance) * CVlambda;
#endif
      
      int const maxpr = 2048;
      PMEDouble xij[maxpr][3];
      int2 ijat0[maxpr];
      PMEDouble rm8[maxpr];
      PMEDouble fsum;
      PMEDouble rm6bar;
      PMEDouble rij;
      int nsum;
      PMEDouble2 R1R2;
      PMEDouble2 R3R4;
      PMEDouble2 K2K3;
      // ifxyz option, each xyz weight will be 0 or 1
      
      PMEDouble xweight = cSim.pNMRr6avDistanceWeights[pos*5];
      PMEDouble yweight = cSim.pNMRr6avDistanceWeights[pos*5+1];
      PMEDouble zweight = cSim.pNMRr6avDistanceWeights[pos*5+2];

#ifdef NMR_TIMED
      int2 Step = cSim.pNMRr6avDistanceStep[pos];
      int Inc = cSim.pNMRr6avDistanceInc[pos];

      // Skip restraint if not active
      if ((step < Step.x) || ((step > Step.y) && (Step.y > 0))) {
        goto exit5;
      }

      // Read timed data
      PMEDouble2 R1R2Slp = cSim.pNMRr6avDistanceR1R2Slp[pos];
      PMEDouble2 R1R2Int = cSim.pNMRr6avDistanceR1R2Int[pos];
      PMEDouble2 R3R4Slp = cSim.pNMRr6avDistanceR3R4Slp[pos];
      PMEDouble2 R3R4Int = cSim.pNMRr6avDistanceR3R4Int[pos];
      PMEDouble2 K2K3Slp = cSim.pNMRr6avDistanceK2K3Slp[pos];
      PMEDouble2 K2K3Int = cSim.pNMRr6avDistanceK2K3Int[pos];

      // Calculate increment
      double dstep = step - (double)((step - Step.x) % abs(Inc));

      // Calculate restraint values
      R1R2.x = dstep*R1R2Slp.x + R1R2Int.x;
      R1R2.y = dstep*R1R2Slp.y + R1R2Int.y;
      R3R4.x = dstep*R3R4Slp.x + R3R4Int.x;
      R3R4.y = dstep*R3R4Slp.y + R3R4Int.y;
      if (Inc > 0) {
        K2K3.x = dstep*K2K3Slp.x + K2K3Int.x;
        K2K3.y = dstep*K2K3Slp.y + K2K3Int.y;
      }
      else {
        int nstepu = (step - Step.x) / abs(Inc);
        K2K3.x     = K2K3Int.x * pow(K2K3Slp.x, nstepu);
        K2K3.y     = K2K3Int.y * pow(K2K3Slp.y, nstepu);
      }
#else
      R1R2 = cSim.pNMRr6avDistanceR1R2[pos];
      R3R4 = cSim.pNMRr6avDistanceR3R4[pos];
      K2K3 = cSim.pNMRr6avDistanceK2K3[pos];
#endif
      PMEDouble atomIX;
      PMEDouble atomIY;
      PMEDouble atomIZ;
      PMEDouble atomJX;
      PMEDouble atomJY;
      PMEDouble atomJZ;
      int ipr = 0;
      rm6bar = (PMEDouble)0.0;
#ifdef NODPTEXTURE
#  ifdef NMR_NEIGHBORLIST
      if (atom.y > 0) {
        atomJX = cSim.pImageX[atom.y];
        atomJY = cSim.pImageY[atom.y];
        atomJZ = cSim.pImageZ[atom.y];

        // Find r6av
        int2 r6avgrps = cSim.pImageNMRr6avDistancer6avGrp[pos * 2];
        for (int ip = r6avgrps.x; ip < r6avgrps.y; ip++) {
          int2 r6avatms = cSim.pImageNMRr6avDistancer6av[ip];
          for (int j = r6avatms.x; j < r6avatms.y + 1; j++ ) {
            xij[ipr][0] = cSim.pImageX[j] - atomJX;
            xij[ipr][1] = cSim.pImageY[j] - atomJY;
            xij[ipr][2] = cSim.pImageZ[j] - atomJZ;
            ijat0[ipr].x   = j;
            ijat0[ipr].y   = atom.y;
            //ifxyz is normally 1, only zero if the user specified 
            PMEDouble rij2 = xij[ipr][0]*xij[ipr][0]*xweight +
                             xij[ipr][1]*xij[ipr][1]*yweight +
                             xij[ipr][2]*xij[ipr][2]*zweight;
            PMEDouble rm2 = (PMEDouble)1.0/rij2;
            //rm6bar        = rm6bar + pow(rm2, 3);
            rm6bar        = rm6bar + (rm2*rm2*rm2);
            //rm8[ipr]      = pow(rm2, 4);
            rm8[ipr]      = rm2*rm2*rm2*rm2;
            ipr = ipr + 1;
          }
        }
        nsum   = ipr;
        fsum   = __uint2double_rn (nsum);
        rm6bar = rm6bar/fsum;
        rij    = pow(rm6bar,-(PMEDouble)1.0/(PMEDouble)6.0);
      }
      else if (atom.x > 0) {
        atomIX = cSim.pImageX[atom.x];
        atomIY = cSim.pImageY[atom.x];
        atomIZ = cSim.pImageZ[atom.x];

        // Find r6av
        int2 r6avgrps = cSim.pImageNMRr6avDistancer6avGrp[pos * 2 + 1];
        for (int ip = r6avgrps.x; ip < r6avgrps.y; ip++) {
          int2 r6avatms = cSim.pImageNMRr6avDistancer6av[ip];
          for (int j = r6avatms.x; j < r6avatms.y + 1; j++) {
            xij[ipr][0] = atomIX - cSim.pImageX[j];
            xij[ipr][1] = atomIY - cSim.pImageY[j];
            xij[ipr][2] = atomIZ - cSim.pImageZ[j];
            ijat0[ipr].x   = atom.x;
            ijat0[ipr].y   = j;
            PMEDouble rij2 = xij[ipr][0]*xij[ipr][0]*xweight +
                             xij[ipr][1]*xij[ipr][1]*yweight +
                             xij[ipr][2]*xij[ipr][2]*zweight;
            PMEDouble rm2 = (PMEDouble)1.0/rij2;
            rm6bar        = rm6bar + (rm2*rm2*rm2);
            rm8[ipr]      = rm2*rm2*rm2*rm2;
            ipr = ipr + 1;
          }
        }
        nsum   = ipr;
        fsum   = __uint2double_rn (nsum);
        rm6bar = rm6bar/fsum;
        rij    = pow(rm6bar,-(PMEDouble)1.0/(PMEDouble)6.0);
      }
      else {

        // Find r6av of 2 groups
        int2 r6avgrpsi = cSim.pImageNMRr6avDistancer6avGrp[pos * 2];
        for (int ip = r6avgrpsi.x; ip < r6avgrpsi.y; ip++) {
          int2 r6avatmsi = cSim.pImageNMRr6avDistancer6av[ip];
          for (int i = r6avatmsi.x; i < r6avatmsi.y + 1; i++) {
            int2 r6avgrpsj = cSim.pImageNMRr6avDistancer6avGrp[pos * 2 + 1];
            for (int jp = r6avgrpsj.x; jp < r6avgrpsj.y; jp++) {
              int2 r6avatmsj          = cSim.pImageNMRr6avDistancer6av[jp];
              for (int j = r6avatmsj.x; j < r6avatmsj.y + 1; j++) {
                xij[ipr][0] = cSim.pImageX[i] - cSim.pImageX[j];
                xij[ipr][1] = cSim.pImageY[i] - cSim.pImageY[j];
                xij[ipr][2] = cSim.pImageZ[i] - cSim.pImageZ[j];
                ijat0[ipr].x   = i;
                ijat0[ipr].y   = j;
                PMEDouble rij2 = xij[ipr][0]*xij[ipr][0]*xweight +
                             xij[ipr][1]*xij[ipr][1]*yweight +
                             xij[ipr][2]*xij[ipr][2]*zweight;
                PMEDouble rm2 = (PMEDouble)1.0/rij2;
                rm6bar        = rm6bar + (rm2*rm2*rm2);
                rm8[ipr]      = rm2*rm2*rm2*rm2;
                ipr = ipr + 1;
              }
            }
          }
        }
        nsum   = ipr;
        fsum   = __uint2double_rn (nsum);
        rm6bar = rm6bar/fsum;
        rij    = pow(rm6bar, -(PMEDouble)1.0/(PMEDouble)6.0);
      }
#  else  // NMR_NEIGHBORLIST
      if (atom.y > 0) {
        atomJX = cSim.pAtomX[atom.y];
        atomJY = cSim.pAtomY[atom.y];
        atomJZ = cSim.pAtomZ[atom.y];

        // Find r6av
        int2 r6avgrps = cSim.pNMRr6avDistancer6avGrp[pos * 2];
        for (int ip = r6avgrps.x; ip < r6avgrps.y; ip++) {
          int2 r6avatms = cSim.pNMRr6avDistancer6av[ip];
          for (int j = r6avatms.x; j < r6avatms.y + 1; j++) {
            xij[ipr][0] = cSim.pAtomX[j] - atomJX;
            xij[ipr][1] = cSim.pAtomY[j] - atomJY;
            xij[ipr][2] = cSim.pAtomZ[j] - atomJZ;
            ijat0[ipr].x   = j;
            ijat0[ipr].y   = atom.y;
            PMEDouble rij2 = xij[ipr][0]*xij[ipr][0]*xweight +
                             xij[ipr][1]*xij[ipr][1]*yweight +
                             xij[ipr][2]*xij[ipr][2]*zweight;
            PMEDouble rm2 = (PMEDouble)1.0/rij2;
            rm6bar        = rm6bar + (rm2*rm2*rm2);
            rm8[ipr]      = rm2*rm2*rm2*rm2;
            ipr = ipr + 1;
          }
        }
        nsum   = ipr;
        fsum   = __uint2double_rn (nsum);
        rm6bar = rm6bar/fsum;
        rij    = pow(rm6bar, -(PMEDouble)1.0/(PMEDouble)6.0);
      }
      else if (atom.x > 0) {
        atomIX = cSim.pAtomX[atom.x];
        atomIY = cSim.pAtomY[atom.x];
        atomIZ = cSim.pAtomZ[atom.x];

        // Find r6av
        int2 r6avgrps = cSim.pNMRr6avDistancer6avGrp[pos * 2 + 1];
        for (int ip = r6avgrps.x; ip < r6avgrps.y; ip++) {
          int2 r6avatms = cSim.pNMRr6avDistancer6av[ip];
          for (int j = r6avatms.x; j < r6avatms.y + 1; j++) {
            xij[ipr][0] = atomIX - cSim.pAtomX[j];
            xij[ipr][1] = atomIY - cSim.pAtomY[j];
            xij[ipr][2] = atomIZ - cSim.pAtomZ[j];
            ijat0[ipr].x   = atom.x;
            ijat0[ipr].y   = j;
            PMEDouble rij2 = xij[ipr][0]*xij[ipr][0]*xweight +
                             xij[ipr][1]*xij[ipr][1]*yweight +
                             xij[ipr][2]*xij[ipr][2]*zweight;
            PMEDouble rm2 = (PMEDouble)1.0/rij2;
            rm6bar        = rm6bar + (rm2*rm2*rm2);
            rm8[ipr]      = rm2*rm2*rm2*rm2;
            ipr = ipr + 1;
          }
        }
        nsum   = ipr;
        fsum   = __uint2double_rn (nsum);
        rm6bar = rm6bar/fsum;
        rij    = pow(rm6bar, -(PMEDouble)1.0/(PMEDouble)6.0);
      }
      else {

        // r6av of 2 groups
        int2 r6avgrpsi = cSim.pNMRr6avDistancer6avGrp[pos * 2];
        for (int ip = r6avgrpsi.x; ip < r6avgrpsi.y; ip++) {
          int2 r6avatmsi = cSim.pNMRr6avDistancer6av[ip];
          for (int i = r6avatmsi.x; i < r6avatmsi.y + 1; i++) {
            int2 r6avgrpsj = cSim.pNMRr6avDistancer6avGrp[pos * 2 + 1];
            for (int jp = r6avgrpsj.x; jp < r6avgrpsj.y; jp++) {
              int2 r6avatmsj = cSim.pNMRr6avDistancer6av[jp];
              for (int j = r6avatmsj.x; j < r6avatmsj.y + 1; j++) {
                xij[ipr][0] = cSim.pAtomX[i] - cSim.pAtomX[j];
                xij[ipr][1] = cSim.pAtomY[i] - cSim.pAtomY[j];
                xij[ipr][2] = cSim.pAtomZ[i] - cSim.pAtomZ[j];
                ijat0[ipr].x   = i;
                ijat0[ipr].y   = j;
                PMEDouble rij2 = xij[ipr][0]*xij[ipr][0]*xweight +
                             xij[ipr][1]*xij[ipr][1]*yweight +
                             xij[ipr][2]*xij[ipr][2]*zweight;
                PMEDouble rm2 = (PMEDouble)1.0/rij2;
                rm6bar        = rm6bar + (rm2*rm2*rm2);
                rm8[ipr]      = rm2*rm2*rm2*rm2;
                ipr = ipr + 1;
              }
            }
          }
        }
        nsum   = ipr;
        fsum   = __uint2double_rn (nsum);
        rm6bar = rm6bar/fsum;
        rij    = pow(rm6bar, -(PMEDouble)1.0/(PMEDouble)6.0);
      }
#  endif // NMR_NEIGHBORLIST
#else  // NODPTEXTURE
      if (atom.y > 0) {
#  ifdef NMR_NEIGHBORLIST
        int2 iatomJX = tex1Dfetch<int2>(cSim.texImageX, atom.y);
        int2 iatomJY = tex1Dfetch<int2>(cSim.texImageX, atom.y + cSim.stride);
        int2 iatomJZ = tex1Dfetch<int2>(cSim.texImageX, atom.y + cSim.stride2);
#  else
        int2 iatomJX = tex1Dfetch<int2>(cSim.texAtomX, atom.y);
        int2 iatomJY = tex1Dfetch<int2>(cSim.texAtomX, atom.y + cSim.stride);
        int2 iatomJZ = tex1Dfetch<int2>(cSim.texAtomX, atom.y + cSim.stride2);
#  endif
        atomJX = __hiloint2double(iatomJX.y, iatomJX.x);

        atomJY = __hiloint2double(iatomJY.y, iatomJY.x);
        atomJZ = __hiloint2double(iatomJZ.y, iatomJZ.x);

        // Find r6av
#  ifdef NMR_NEIGHBORLIST
        int2 r6avgrps = cSim.pImageNMRr6avDistancer6avGrp[pos * 2];
#  else
        int2 r6avgrps = cSim.pNMRr6avDistancer6avGrp[pos * 2];
#  endif
        for (int ip = r6avgrps.x; ip < r6avgrps.y; ip++) {
#  ifdef NMR_NEIGHBORLIST
          int2 r6avatms = cSim.pImageNMRr6avDistancer6av[ip];
#  else
          int2 r6avatms = cSim.pNMRr6avDistancer6av[ip];
#  endif
          for (int j = r6avatms.x; j < r6avatms.y + 1; j++) {
#  ifdef NMR_NEIGHBORLIST
            int2 iatomIX = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomIY = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomIZ = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
#  else
            int2 iatomIX = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomIY = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomIZ = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
#  endif
            atomIX = __hiloint2double(iatomIX.y, iatomIX.x);
            atomIY = __hiloint2double(iatomIY.y, iatomIY.x);
            atomIZ = __hiloint2double(iatomIZ.y, iatomIZ.x);
            xij[ipr][0] = atomIX - atomJX;
            xij[ipr][1] = atomIY - atomJY;
            xij[ipr][2] = atomIZ - atomJZ;
            ijat0[ipr].x   = j;
            ijat0[ipr].y   = atom.y;
            PMEDouble rij2 = xij[ipr][0]*xij[ipr][0]*xweight +
                             xij[ipr][1]*xij[ipr][1]*yweight +
                             xij[ipr][2]*xij[ipr][2]*zweight;
            PMEDouble rm2 = (PMEDouble)1.0/rij2;
            rm6bar        = rm6bar + (rm2*rm2*rm2);
            rm8[ipr]      = rm2*rm2*rm2*rm2;
            ipr = ipr + 1;
          }
        }
        nsum   = ipr;
        fsum   = __uint2double_rn (nsum);
        rm6bar = rm6bar/fsum;
        rij    = pow(rm6bar, -(PMEDouble)1.0/(PMEDouble)6.0);
      }
      else if (atom.x > 0) {
#  ifdef NMR_NEIGHBORLIST
        int2 iatomIX = tex1Dfetch<int2>(cSim.texImageX, atom.x);
        int2 iatomIY = tex1Dfetch<int2>(cSim.texImageX, atom.x + cSim.stride);
        int2 iatomIZ = tex1Dfetch<int2>(cSim.texImageX, atom.x + cSim.stride2);
#  else
        int2 iatomIX = tex1Dfetch<int2>(cSim.texAtomX, atom.x);
        int2 iatomIY = tex1Dfetch<int2>(cSim.texAtomX, atom.x + cSim.stride);
        int2 iatomIZ = tex1Dfetch<int2>(cSim.texAtomX, atom.x + cSim.stride2);
#  endif
        atomIX = __hiloint2double(iatomIX.y, iatomIX.x);
        atomIY = __hiloint2double(iatomIY.y, iatomIY.x);
        atomIZ = __hiloint2double(iatomIZ.y, iatomIZ.x);

        // Find r6av
#  ifdef NMR_NEIGHBORLIST
        int2 r6avgrps = cSim.pImageNMRr6avDistancer6avGrp[pos * 2 + 1];
#  else
        int2 r6avgrps = cSim.pNMRr6avDistancer6avGrp[pos * 2 + 1];
#  endif
        for (int ip = r6avgrps.x; ip < r6avgrps.y; ip++) {
#  ifdef NMR_NEIGHBORLIST
          int2 r6avatms = cSim.pImageNMRr6avDistancer6av[ip];
#  else
          int2 r6avatms = cSim.pNMRr6avDistancer6av[ip];
#  endif
          for (int j = r6avatms.x; j < r6avatms.y + 1; j++) {
#  ifdef NMR_NEIGHBORLIST
            int2 iatomJX = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomJY = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomJZ = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
#  else
            int2 iatomJX = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomJY = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomJZ = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
#  endif
            atomJX    = __hiloint2double(iatomJX.y, iatomJX.x);
            atomJY    = __hiloint2double(iatomJY.y, iatomJY.x);
            atomJZ    = __hiloint2double(iatomJZ.y, iatomJZ.x);
            xij[ipr][0] = atomIX - atomJX;
            xij[ipr][1] = atomIY - atomJY;
            xij[ipr][2] = atomIZ - atomJZ;
            ijat0[ipr].x   = atom.x;
            ijat0[ipr].y   = j;
            PMEDouble rij2 = xij[ipr][0]*xij[ipr][0]*xweight +
                             xij[ipr][1]*xij[ipr][1]*yweight +
                             xij[ipr][2]*xij[ipr][2]*zweight;
            PMEDouble rm2 = (PMEDouble)1.0/rij2;
            rm6bar        = rm6bar + (rm2*rm2*rm2);
            rm8[ipr]      = rm2*rm2*rm2*rm2;
            ipr = ipr + 1;
          }
        }
        nsum     = ipr;
        fsum     = __uint2double_rn (nsum);
        rm6bar   = rm6bar/fsum;
        rij      = pow(rm6bar,-(PMEDouble)1.0/(PMEDouble)6.0);
      }
      else {

        // r6av of 2 groups
#  ifdef NMR_NEIGHBORLIST
        int2 r6avgrpsi = cSim.pImageNMRr6avDistancer6avGrp[pos * 2];
#  else
        int2 r6avgrpsi = cSim.pNMRr6avDistancer6avGrp[pos * 2];
#  endif
        for (int ip = r6avgrpsi.x; ip < r6avgrpsi.y; ip++) {
#  ifdef NMR_NEIGHBORLIST
          int2 r6avatmsi = cSim.pImageNMRr6avDistancer6av[ip];
#  else
          int2 r6avatmsi = cSim.pNMRr6avDistancer6av[ip];
#  endif
          for (int i = r6avatmsi.x; i < r6avatmsi.y + 1; i++ ) {
#  ifdef NMR_NEIGHBORLIST
            int2 r6avgrpsj = cSim.pImageNMRr6avDistancer6avGrp[pos * 2 + 1];
#  else
            int2 r6avgrpsj = cSim.pNMRr6avDistancer6avGrp[pos * 2 + 1];
#  endif
            for (int jp = r6avgrpsj.x; jp < r6avgrpsj.y; jp++) {
#  ifdef NMR_NEIGHBORLIST
              int2 r6avatmsj = cSim.pImageNMRr6avDistancer6av[jp];
#  else
              int2 r6avatmsj = cSim.pNMRr6avDistancer6av[jp];
#  endif
              for (int j = r6avatmsj.x; j < r6avatmsj.y + 1; j++) {
#  ifdef NMR_NEIGHBORLIST
                int2 iatomIX   = tex1Dfetch<int2>(cSim.texImageX, i);
                int2 iatomIY   = tex1Dfetch<int2>(cSim.texImageX, i + cSim.stride);
                int2 iatomIZ   = tex1Dfetch<int2>(cSim.texImageX, i + cSim.stride2);
#  else
                int2 iatomIX   = tex1Dfetch<int2>(cSim.texAtomX, i);
                int2 iatomIY   = tex1Dfetch<int2>(cSim.texAtomX, i + cSim.stride);
                int2 iatomIZ   = tex1Dfetch<int2>(cSim.texAtomX, i + cSim.stride2);
#  endif
                atomIX         = __hiloint2double(iatomIX.y, iatomIX.x);
                atomIY         = __hiloint2double(iatomIY.y, iatomIY.x);
                atomIZ         = __hiloint2double(iatomIZ.y, iatomIZ.x);
#  ifdef NMR_NEIGHBORLIST
                int2 iatomJX   = tex1Dfetch<int2>(cSim.texImageX, j);
                int2 iatomJY   = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
                int2 iatomJZ   = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
#  else
                int2 iatomJX   = tex1Dfetch<int2>(cSim.texAtomX, j);
                int2 iatomJY   = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
                int2 iatomJZ   = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
#  endif
                atomJX         = __hiloint2double(iatomJX.y, iatomJX.x);
                atomJY         = __hiloint2double(iatomJY.y, iatomJY.x);
                atomJZ         = __hiloint2double(iatomJZ.y, iatomJZ.x);
                xij[ipr][0]    = atomIX - atomJX;
                xij[ipr][1]    = atomIY - atomJY;
                xij[ipr][2]    = atomIZ - atomJZ;
                ijat0[ipr].x   = i;
                ijat0[ipr].y   = j;
                PMEDouble rij2 = xij[ipr][0]*xij[ipr][0]*xweight +
                                 xij[ipr][1]*xij[ipr][1]*yweight +
                                 xij[ipr][2]*xij[ipr][2]*zweight;
                PMEDouble rm2 = (PMEDouble)1.0/rij2;
                rm6bar        = rm6bar + (rm2*rm2*rm2);
                rm8[ipr]      = rm2*rm2*rm2*rm2;
                ipr = ipr + 1;
              }
            }
          }
        }
        nsum   = ipr;
        fsum   = __uint2double_rn (nsum);
        rm6bar = rm6bar/fsum;
        rij    = pow(rm6bar, (PMEDouble)-1.0/(PMEDouble)6.0);
      }
#endif // NODPTEXTURE
      PMEDouble df;
#ifdef NMR_ENERGY
      PMEDouble e;
#endif
      if (rij < R1R2.x) {
        PMEDouble dif = R1R2.x - R1R2.y;
        df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
        e             = df * (rij - R1R2.x) + K2K3.x * dif * dif;
#  ifdef NMR_AFE
        e            *= lambda;
#  endif
#endif
      }
      else if (rij < R1R2.y) {
        PMEDouble dif = rij - R1R2.y;
        df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
        e             = K2K3.x * dif * dif;
#  ifdef NMR_AFE
        e            *= lambda;
#  endif
#endif
      }
      else if (rij < R3R4.x) {
        df            = (PMEDouble)0.0;
#ifdef NMR_ENERGY
        e             = (PMEDouble)0.0;
#endif
      }
      else if (rij < R3R4.y) {
        PMEDouble dif = rij - R3R4.x;
        df            = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
        e             = K2K3.y * dif * dif;
#  ifdef NMR_AFE
        e            *= lambda;
#  endif
#endif
      }
      else {
        PMEDouble dif = R3R4.y - R3R4.x;
        df            = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
        e             = df * (rij - R3R4.y) + K2K3.y * dif * dif;
#  ifdef NMR_AFE
        e            *= lambda;
#  endif
#endif
      }
      if (cSim.bJar) {
        double fold  = cSim.pNMRJarData[2];
        double work  = cSim.pNMRJarData[3];
        double first = cSim.pNMRJarData[4];
        double fcurr = (PMEDouble)-2.0 * K2K3.x * (rij - R1R2.y);
        if (first == (PMEDouble)0.0) {
          fold = -fcurr;
          cSim.pNMRJarData[4] = (PMEDouble)1.0;
        }
        work                += (PMEDouble)0.5 * (fcurr + fold) * cSim.drjar;
        cSim.pNMRJarData[0]  = R1R2.y;
        cSim.pNMRJarData[1]  = rij;
        cSim.pNMRJarData[2]  = fcurr;
        cSim.pNMRJarData[3]  = work;
      }

#ifdef NMR_ENERGY
#  ifdef NMR_AFE
      scr6avdistance[TIr6avdistanceRegion - SCr6avdistance] += ENERGYSCALE * e *
                                                               SCr6avdistance;
      if (CVr6avdistance){
        dvdl_r6avdistance += llrint(ENERGYSCALE * e * sTISigns[TIr6avdistanceRegion - 1]);
      }
      e            *= lambda;
#  endif
      er6avdistance += llrint(ENERGYSCALE * e);
#endif
      PMEDouble fact = df * (rij / (fsum * rm6bar));
      int i3, j3;
      PMEDouble fact2;
      PMEDouble dfx;
      PMEDouble dfy;
      PMEDouble dfz;
      PMEAccumulator ifx;
      PMEAccumulator ify;
      PMEAccumulator ifz;
      for (int ipr = 0; ipr < nsum; ipr++) {
        i3    = ijat0[ipr].x;
        j3    = ijat0[ipr].y;
        fact2 = fact * rm8[ipr];
        // weight are only 0 or 1, 
        dfx   = fact2 * xij[ipr][0] * xweight;
        dfy   = fact2 * xij[ipr][1] * yweight;
        dfz   = fact2 * xij[ipr][2] * zweight;
        ifx = llrint(dfx * FORCESCALE);
        ify = llrint(dfy * FORCESCALE);
        ifz = llrint(dfz * FORCESCALE);
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[i3],
                  llitoulli(-ifx));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[i3],
                  llitoulli(-ify));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[i3],
                  llitoulli(-ifz));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[j3], llitoulli(ifx));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[j3], llitoulli(ify));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[j3], llitoulli(ifz));
      }
    }
    exit5: {}
#ifdef NMR_ENERGY
    er6avdistance += __SHFL(WARP_MASK, er6avdistance, threadIdx.x ^ 1);
    er6avdistance += __SHFL(WARP_MASK, er6avdistance, threadIdx.x ^ 2);
    er6avdistance += __SHFL(WARP_MASK, er6avdistance, threadIdx.x ^ 4);
    er6avdistance += __SHFL(WARP_MASK, er6avdistance, threadIdx.x ^ 8);
    er6avdistance += __SHFL(WARP_MASK, er6avdistance, threadIdx.x ^ 16);
#ifdef AMBER_PLATFORM_AMD_WARP64
    er6avdistance += __SHFL(WARP_MASK, er6avdistance, threadIdx.x ^ 32);
#endif
    // Write out energies
    if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    __syncthreads();
    atomicAdd(cSim.pENMRr6avDistance, llitoulli(er6avdistance)); 
    }
#  ifdef NMR_AFE
    sDVDL[threadIdx.x]  = dvdl_r6avdistance;
    sESCRestR1[threadIdx.x]  = scr6avdistance[0];
    sESCRestR2[threadIdx.x]  = scr6avdistance[1];
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 1];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 1];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 1];
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 2];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 2];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 2];
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 4];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 4];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 4];
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 8];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 8];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 8];
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 16];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 16];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 16];
#ifdef AMBER_PLATFORM_AMD_WARP64
    sDVDL[threadIdx.x] += sDVDL[threadIdx.x ^ 32];
    sESCRestR1[threadIdx.x] += sESCRestR1[threadIdx.x ^ 32];
    sESCRestR2[threadIdx.x] += sESCRestR2[threadIdx.x ^ 32];
#endif
    if ((threadIdx.x & GRID_BITS_MASK) == 0) {
      atomicAdd(cSim.pDVDL, llitoulli(sDVDL[threadIdx.x]));
      atomicAdd(cSim.pESCNMRr6avDistanceR1, llitoulli(sESCRestR1[threadIdx.x]));
      atomicAdd(cSim.pESCNMRr6avDistanceR2, llitoulli(sESCRestR2[threadIdx.x]));
    }
#  endif
#endif
  }
#endif

 //Calculate COM Angle restraints
  else if (pos >= cSim.NMRr6avDistanceOffset && pos < cSim.NMRCOMAngleOffset) {
#ifdef NMR_ENERGY
    PMEAccumulator eCOMangle         = 0;
#endif
    pos                             -= cSim.NMRr6avDistanceOffset;
    if (pos < cSim.NMRCOMAngles){
#ifdef NMR_NEIGHBORLIST
      int2 atom1                   = cSim.pImageNMRCOMAngleID1[pos];
      int  atom2                   = cSim.pImageNMRCOMAngleID2[pos];
#else
      int2 atom1                   = cSim.pNMRCOMAngleID1[pos];
      int  atom2                   = cSim.pNMRCOMAngleID2[pos];
#endif

      PMEDouble2 rmstot1;
      PMEDouble  rmstot2;
      rmstot1.x = 0.0;
      rmstot1.y = 0.0;
      rmstot2   = 0.0;
      PMEDouble2 R1R2;
      PMEDouble2 R3R4;
      PMEDouble2 K2K3;
#ifdef NMR_TIMED
      int2 Step                    = cSim.pNMRCOMAngleStep[pos];
      int  Inc                     = cSim.pNMRCOMAngleInc[pos];
      // skip restraint if not active
      if ((step < Step.x) || ((step > Step.y) && (Step.y > 0))) {
        goto exit6;
      }
      // Read timed data
      PMEDouble2 R1R2Slp = cSim.pNMRCOMAngleR1R2Slp[pos];
      PMEDouble2 R1R2Int = cSim.pNMRCOMAngleR1R2Int[pos];
      PMEDouble2 R3R4Slp = cSim.pNMRCOMAngleR3R4Slp[pos];
      PMEDouble2 R3R4Int = cSim.pNMRCOMAngleR3R4Int[pos];
      PMEDouble2 K2K3Slp = cSim.pNMRCOMAngleK2K3Slp[pos];
      PMEDouble2 K2K3Int = cSim.pNMRCOMAngleK2K3Int[pos];

      // Calculate increment
      double dstep = step - (double)((step - Step.x) % abs(Inc));

      // Calculate restraint values
      R1R2.x = dstep * R1R2Slp.x + R1R2Int.x;
      R1R2.y = dstep * R1R2Slp.y + R1R2Int.y;
      R3R4.x = dstep * R3R4Slp.x + R3R4Int.x;
      R3R4.y = dstep * R3R4Slp.y + R3R4Int.y;
      if (Inc > 0) {
        K2K3.x = dstep * K2K3Slp.x + K2K3Int.x;
        K2K3.y = dstep * K2K3Slp.y + K2K3Int.y;
      }
      else {
        int nstepu = (step - Step.x) / abs(Inc);
        K2K3.x     = K2K3Int.x * pow(K2K3Slp.x, nstepu);
        K2K3.y     = K2K3Int.y * pow(K2K3Slp.y, nstepu);
      }
#else
      R1R2 = cSim.pNMRCOMAngleR1R2[pos];
      R3R4 = cSim.pNMRCOMAngleR3R4[pos];
      K2K3 = cSim.pNMRCOMAngleK2K3[pos];
#endif
      PMEDouble atomIX;
      PMEDouble atomIY;
      PMEDouble atomIZ;
      PMEDouble atomJX;
      PMEDouble atomJY;
      PMEDouble atomJZ;
      PMEDouble atomKX;
      PMEDouble atomKY;
      PMEDouble atomKZ;
#ifdef NODPTEXTURE
  #ifdef NMR_NEIGHBORLIST
      //COM_1: NL
      if (atom1.x > 0) {
        atomIX = cSim.pImageX[atom1.x];
        atomIY = cSim.pImageY[atom1.x];
        atomIZ = cSim.pImageZ[atom1.x];
      }
      else {
      // find COM
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot1.x = 0.0;
        int2 COMgrps   = cSim.pImageNMRCOMAngleCOMGrp[pos * 3];
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++)  {
          int2 COMatms = cSim.pImageNMRCOMAngleCOM[ip];
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass = cSim.pImageMass[j];
            xtot      = xtot + (cSim.pImageX[j] * rmass);
            ytot      = ytot + (cSim.pImageY[j] * rmass);
            ztot      = ztot + (cSim.pImageZ[j] * rmass);
            rmstot1.x = rmstot1.x + rmass;
          }
        }
        atomIX = xtot / rmstot1.x;
        atomIY = ytot / rmstot1.x;
        atomIZ = ztot / rmstot1.x;
      }
      // COM_2: NL
      if (atom1.y > 0) {
        atomJX = cSim.pImageX[atom1.y];
        atomJY = cSim.pImageY[atom1.y];
        atomJZ = cSim.pImageZ[atom1.y];
      }
      else {
        // find the center of mass
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot1.y      = 0.0;
        int2 COMgrps   = cSim.pImageNMRCOMAngleCOMGrp[pos * 3 + 1];
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
          int2 COMatms = cSim.pImageNMRCOMAngleCOM[ip];
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass = cSim.pImageMass[j];
            xtot      = xtot + (cSim.pImageX[j] * rmass);
            ytot      = ytot + (cSim.pImageY[j] * rmass);
            ztot      = ztot + (cSim.pImageZ[j] * rmass);
            rmstot1.y = rmstot1.y + rmass;
          }
        }
        atomJX      = xtot / rmstot1.y;
        atomJY      = ytot / rmstot1.y;
        atomJZ      = ztot / rmstot1.y;
      }
      // COM_3: NL
      if (atom2 > 0) {
        atomKX = cSim.pImageX[atom2];
        atomKY = cSim.pImageY[atom2];
        atomKZ = cSim.pImageZ[atom2];
      }
      // if COM_3 has multiple atoms:
      else {
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot2        = 0.0;
        int2 COMgrps   = cSim.pImageNMRCOMAngleCOMGrp[pos * 3 + 2];
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
          int2 COMatms = cSim.pImageNMRCOMAngleCOM[ip];
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass = cSim.pImageMass[j];
            xtot    = xtot + (cSim.pImageX[j] * rmass);
            ytot    = ytot + (cSim.pImageY[j] * rmass);
            ztot    = ztot + (cSim.pImageZ[j] * rmass);
            rmstot2 = rmstot2 + rmass;
          }
        }
        atomKX = xtot / rmstot2;
        atomKY = ytot / rmstot2;
        atomKZ = ztot / rmstot2;
      }
  #else
      //COM_1: noNeighborList (NL)
      if (atom1.x > 0) {
        atomIX = cSim.pAtomX[atom1.x];
        atomIY = cSim.pAtomY[atom1.x];
        atomIZ = cSim.pAtomZ[atom1.x];
      }
      else {
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot1.x = 0.0;
        int2 COMgrps = cSim.pNMRCOMAngleCOMGrp[pos * 3];
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
          int2 COMatms = cSim.pNMRCOMAngleCOM[ip];
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass      = cSim.pAtomMass[j];
            xtot      = xtot + cSim.pAtomX[j] * rmass;
            ytot      = ytot + cSim.pAtomY[j] * rmass;
            ztot      = ztot + cSim.pAtomZ[j] * rmass;
            rmstot1.x = rmstot1.x + rmass;
          }
        }
        atomIX = xtot / rmstot1.x;
        atomIY = ytot / rmstot1.x;
        atomIZ = ztot / rmstot1.x;
      }
      //COM_2: noNL
      if (atom1.y > 0) {
        atomJX = cSim.pAtomX[atom1.y];
        atomJY = cSim.pAtomY[atom1.y];
        atomJZ = cSim.pAtomZ[atom1.y];
      }
      else {
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot1.y  = 0.0;
        int2 COMgrps = cSim.pNMRCOMAngleCOMGrp[pos * 3 + 1];
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
          int2 COMatms = cSim.pNMRCOMAngleCOM[ip];
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass      = cSim.pAtomMass[j];
            xtot                 = xtot + cSim.pAtomX[j] * rmass;
            ytot                 = ytot + cSim.pAtomY[j] * rmass;
            ztot                 = ztot + cSim.pAtomZ[j] * rmass;
            rmstot1.y            = rmstot1.y + rmass;
          }
        }
        atomJX      = xtot / rmstot1.y;
        atomJY      = ytot / rmstot1.y;
        atomJZ      = ztot / rmstot1.y;
      }
      //COM_3: noNL
      if (atom2 > 0) {
        atomKX      = cSim.pAtomX[atom2];
        atomKY      = cSim.pAtomY[atom2];
        atomKZ      = cSim.pAtomZ[atom2];
      }
      else {
        PMEDouble xtot  = 0.0;
        PMEDouble ytot  = 0.0;
        PMEDouble ztot  = 0.0;
        rmstot2         = 0.0;
        int2 COMgrps = cSim.pNMRCOMAngleCOMGrp[pos * 3 + 2];
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
          int2 COMatms       = cSim.pNMRCOMAngleCOM[ip];
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass = cSim.pAtomMass[j];
            xtot    = xtot + cSim.pAtomX[j] * rmass;
            ytot    = ytot + cSim.pAtomY[j] * rmass;
            ztot    = ztot + cSim.pAtomZ[j] * rmass;
            rmstot2 = rmstot2 + rmass;
          }
        }
        atomKX = xtot / rmstot2;
        atomKY = ytot / rmstot2;
        atomKZ = ztot / rmstot2;
      }
  #endif
#else // NODPTEXTURE

      //COM_1: SPXP
      if (atom1.x > 0) {
#ifdef NMR_NEIGHBORLIST
        int2 iatomIX  = tex1Dfetch<int2>(cSim.texImageX, atom1.x);
        int2 iatomIY  = tex1Dfetch<int2>(cSim.texImageX, atom1.x + cSim.stride);
        int2 iatomIZ  = tex1Dfetch<int2>(cSim.texImageX, atom1.x + cSim.stride2);
#else
        int2 iatomIX  = tex1Dfetch<int2>(cSim.texAtomX, atom1.x);
        int2 iatomIY  = tex1Dfetch<int2>(cSim.texAtomX, atom1.x + cSim.stride);
        int2 iatomIZ  = tex1Dfetch<int2>(cSim.texAtomX, atom1.x + cSim.stride2);
#endif
        atomIX = __hiloint2double(iatomIX.y, iatomIX.x);
        atomIY = __hiloint2double(iatomIY.y, iatomIY.x);
        atomIZ = __hiloint2double(iatomIZ.y, iatomIZ.x);
      }
      else {
        // Find the center of mass
        PMEDouble xtot = 0.0;
        PMEDouble ytot = 0.0;
        PMEDouble ztot = 0.0;
        rmstot1.x  = 0.0;
  #ifdef NMR_NEIGHBORLIST
        int2 COMgrps = cSim.pImageNMRCOMAngleCOMGrp[pos * 3];
  #else
        int2 COMgrps = cSim.pNMRCOMAngleCOMGrp[pos * 3];
  #endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
  #ifdef NMR_NEIGHBORLIST
          int2 COMatms        = cSim.pImageNMRCOMAngleCOM[ip];
  #else
          int2 COMatms        = cSim.pNMRCOMAngleCOM[ip];
  #endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
  #ifdef NMR_NEIGHBORLIST
            PMEDouble rmass = cSim.pImageMass[j];
  #else
            PMEDouble rmass  = cSim.pAtomMass[j];
  #endif
  #ifdef NMR_NEIGHBORLIST
            int2 iatomIX = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomIY = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomIZ = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
  #else
            int2 iatomIX = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomIY = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomIZ = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
  #endif
            atomIX = __hiloint2double(iatomIX.y, iatomIX.x);
            atomIY = __hiloint2double(iatomIY.y, iatomIY.x);
            atomIZ = __hiloint2double(iatomIZ.y, iatomIZ.x);
            xtot = xtot + atomIX * rmass;
            ytot = ytot + atomIY * rmass;
            ztot = ztot + atomIZ * rmass;
            rmstot1.x        = rmstot1.x + rmass;
          }
        }
        atomIX = xtot / rmstot1.x;
        atomIY = ytot / rmstot1.x;
        atomIZ = ztot / rmstot1.x;
      }
      //COM_2: SPXP
      if (atom1.y > 0) {
  #ifdef NMR_NEIGHBORLIST
        int2 iatomJX = tex1Dfetch<int2>(cSim.texImageX, atom1.y);
        int2 iatomJY = tex1Dfetch<int2>(cSim.texImageX, atom1.y + cSim.stride);
        int2 iatomJZ = tex1Dfetch<int2>(cSim.texImageX, atom1.y + cSim.stride2);
  #else
        int2 iatomJX = tex1Dfetch<int2>(cSim.texAtomX, atom1.y);
        int2 iatomJY = tex1Dfetch<int2>(cSim.texAtomX, atom1.y + cSim.stride);
        int2 iatomJZ = tex1Dfetch<int2>(cSim.texAtomX, atom1.y + cSim.stride2);
  #endif
        atomJX = __hiloint2double(iatomJX.y, iatomJX.x);
        atomJY = __hiloint2double(iatomJY.y, iatomJY.x);
        atomJZ = __hiloint2double(iatomJZ.y, iatomJZ.x);
      }
      else {
        PMEDouble xtot  = 0.0;
        PMEDouble ytot  = 0.0;
        PMEDouble ztot  = 0.0;
        rmstot1.y       = 0.0;
  #ifdef NMR_NEIGHBORLIST
        int2 COMgrps = cSim.pImageNMRCOMAngleCOMGrp[pos * 3 + 1];
  #else
        int2 COMgrps = cSim.pNMRCOMAngleCOMGrp[pos * 3 + 1];
  #endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
  #ifdef NMR_NEIGHBORLIST
          int2 COMatms = cSim.pImageNMRCOMAngleCOM[ip];
  #else
          int2 COMatms = cSim.pNMRCOMAngleCOM[ip];
  #endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
  #ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
  #else
            PMEDouble rmass  = cSim.pAtomMass[j];
  #endif
  #ifdef NMR_NEIGHBORLIST
            int2 iatomJX     = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomJY     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomJZ     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
  #else
            int2 iatomJX     = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomJY     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomJZ     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
  #endif
            atomJX           = __hiloint2double(iatomJX.y, iatomJX.x);
            atomJY           = __hiloint2double(iatomJY.y, iatomJY.x);
            atomJZ           = __hiloint2double(iatomJZ.y, iatomJZ.x);
            xtot             = xtot + atomJX * rmass;
            ytot             = ytot + atomJY * rmass;
            ztot             = ztot + atomJZ * rmass;
            rmstot1.y        = rmstot1.y + rmass;
          }
        }
        atomJX       = xtot / rmstot1.y;
        atomJY       = ytot / rmstot1.y;
        atomJZ       = ztot / rmstot1.y;
      }
      //COM_3: SPXP
      if (atom2 > 0) {
  #ifdef NMR_NEIGHBORLIST
        int2 iatomKX             = tex1Dfetch<int2>(cSim.texImageX, atom2);
        int2 iatomKY             = tex1Dfetch<int2>(cSim.texImageX, atom2 + cSim.stride);
        int2 iatomKZ             = tex1Dfetch<int2>(cSim.texImageX, atom2 + cSim.stride2);
  #else
        int2 iatomKX             = tex1Dfetch<int2>(cSim.texAtomX, atom2);
        int2 iatomKY             = tex1Dfetch<int2>(cSim.texAtomX, atom2 + cSim.stride);
        int2 iatomKZ             = tex1Dfetch<int2>(cSim.texAtomX, atom2 + cSim.stride2);
  #endif
        atomKX                   = __hiloint2double(iatomKX.y, iatomKX.x);
        atomKY                   = __hiloint2double(iatomKY.y, iatomKY.x);
        atomKZ                   = __hiloint2double(iatomKZ.y, iatomKZ.x);
      }
      else {
        PMEDouble xtot           = 0.0;
        PMEDouble ytot           = 0.0;
        PMEDouble ztot           = 0.0;
        rmstot2                  = 0.0;
  #ifdef NMR_NEIGHBORLIST
        int2 COMgrps = cSim.pImageNMRCOMAngleCOMGrp[pos * 3 + 2];
  #else
        int2 COMgrps = cSim.pNMRCOMAngleCOMGrp[pos * 3 + 2];
  #endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
  #ifdef NMR_NEIGHBORLIST
          int2 COMatms = cSim.pImageNMRCOMAngleCOM[ip];
  #else
          int2 COMatms = cSim.pNMRCOMAngleCOM[ip];
  #endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++ ) {
  #ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
  #else
            PMEDouble rmass  = cSim.pAtomMass[j];
  #endif
  #ifdef NMR_NEIGHBORLIST
            int2 iatomKX     = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomKY     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomKZ     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
  #else
            int2 iatomKX     = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomKY     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomKZ     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
  #endif
            atomKX   = __hiloint2double(iatomKX.y, iatomKX.x);
            atomKY   = __hiloint2double(iatomKY.y, iatomKY.x);
            atomKZ   = __hiloint2double(iatomKZ.y, iatomKZ.x);
            xtot     = xtot + atomKX * rmass;
            ytot     = ytot + atomKY * rmass;
            ztot     = ztot + atomKZ * rmass;
            rmstot2  = rmstot2 + rmass;
          }
        }
        atomKX = xtot / rmstot2;
        atomKY = ytot / rmstot2;
        atomKZ = ztot / rmstot2;
      }
#endif
      PMEDouble xij                = atomIX - atomJX;
      PMEDouble xkj                = atomKX - atomJX;
      PMEDouble yij                = atomIY - atomJY;
      PMEDouble ykj                = atomKY - atomJY;
      PMEDouble zij                = atomIZ - atomJZ;
      PMEDouble zkj                = atomKZ - atomJZ;
      PMEDouble rij2               = xij * xij  +  yij * yij  +  zij * zij ;
      PMEDouble rkj2               = xkj * xkj  +  ykj * ykj  +  zkj * zkj ;
      PMEDouble rij                = sqrt( rij2 );
      PMEDouble rkj                = sqrt( rkj2 );
      PMEDouble rdenom             = rij * rkj;
      PMEDouble cst                = min(   1.0, max(  -1.0, ( xij * xkj  +  yij * ykj  +  zij * zkj )  /  rdenom  )   );
      PMEDouble theta              = acos(cst);
      PMEDouble df;
#ifdef NMR_ENERGY
      PMEDouble e;
#endif
      if (theta < R1R2.x) {
        PMEDouble dif = R1R2.x - R1R2.y;
        df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
        e             = df * (theta - R1R2.x) + K2K3.x * dif * dif;
#endif
      }
      else if (theta < R1R2.y) {
        PMEDouble dif = theta - R1R2.y;
        df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
        e             = K2K3.x * dif * dif;
#endif
      }
      else if (theta < R3R4.x) {
        df            = (PMEDouble)0.0;
#ifdef NMR_ENERGY
        e             = (PMEDouble)0.0;
#endif
      }
      else if (theta < R3R4.y) {
        PMEDouble dif  = theta - R3R4.x;
        df             = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
        e              = K2K3.y * dif * dif;
#endif
      }
      else {
        PMEDouble dif  = R3R4.y - R3R4.x;
        df             = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
        e              = df * ( theta - R3R4.y) + K2K3.y * dif * dif;
#endif
      }
      if (cSim.bJar) {
        double fold             = cSim.pNMRJarData[2];
        double work             = cSim.pNMRJarData[3];
        double first            = cSim.pNMRJarData[4];
        double fcurr            = (PMEDouble)-2.0 * K2K3.x * ( theta - R1R2.y);
        if (first == (PMEDouble)0.0) {
          fold                    = -fcurr;
          cSim.pNMRJarData[4]     = (PMEDouble)1.0;
        }
        work                   += (PMEDouble)0.5 * ( fcurr + fold ) * cSim.drjar;
        cSim.pNMRJarData[0]     = R1R2.y;
        cSim.pNMRJarData[1]     = theta;
        cSim.pNMRJarData[2]     = fcurr;
        cSim.pNMRJarData[3]     = work;
      }
#ifdef NMR_ENERGY
      eCOMangle                 += llrint(ENERGYSCALE  * e);
      //printf("eCOMangle:",eCOMangle,"\n");
#endif

      PMEDouble snt                = sin(theta);
      if ( abs(snt) < (PMEDouble)1.0e-14 )
        snt                 = (PMEDouble)1.0e-14;
      PMEDouble st                 = -df / snt;
      PMEDouble sth                = st * cst;
      PMEDouble cik                = st / rdenom;
      PMEDouble cii                = sth / rij2;
      PMEDouble ckk                = sth / rkj2;
      PMEDouble fx1                = cii * xij - cik * xkj;
      PMEDouble fy1                = cii * yij - cik * ykj;
      PMEDouble fz1                = cii * zij - cik * zkj;
      PMEDouble fx2                = ckk * xkj - cik * xij;
      PMEDouble fy2                = ckk * ykj - cik * yij;
      PMEDouble fz2                = ckk * zkj - cik * zij;
      //COM stuff:
      PMEAccumulator       ifx1;
      PMEAccumulator       ify1;
      PMEAccumulator       ifz1;
      PMEAccumulator       ifx2;
      PMEAccumulator       ify2;
      PMEAccumulator       ifz2;
      //COM_1:
      if (atom1.x > 0) {
        PMEAccumulator     ifx1    = llrint(  fx1 * FORCESCALE  );
        PMEAccumulator     ify1    = llrint(  fy1 * FORCESCALE  );
        PMEAccumulator     ifz1    = llrint(  fz1 * FORCESCALE  );
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.x],
                  llitoulli(ifx1));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.x],
                  llitoulli(ify1));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.x],
                  llitoulli(ifz1));
      }
      else {
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps = cSim.pImageNMRCOMAngleCOMGrp[pos * 3];
#else
        int2 COMgrps = cSim.pNMRCOMAngleCOMGrp[pos * 3];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#ifdef NMR_NEIGHBORLIST
          int2 COMatms = cSim.pImageNMRCOMAngleCOM[ip];
#else
          int2 COMatms = cSim.pNMRCOMAngleCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
            PMEDouble dcomdx = rmass / rmstot1.x;
            PMEDouble kfx1   = fx1 * dcomdx;
            PMEDouble kfy1   = fy1 * dcomdx;
            PMEDouble kfz1   = fz1 * dcomdx;
            ifx1   = llrint(kfx1 * FORCESCALE);
            ify1   = llrint(kfy1 * FORCESCALE);
            ifz1   = llrint(kfz1 * FORCESCALE);
            atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[j],
                      llitoulli(ifx1));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[j],
                      llitoulli(ify1));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[j],
                      llitoulli(ifz1));
          }
        }
      }
      //COM_2:
      if (atom1.y > 0) {
        PMEAccumulator      ifx1   = llrint(  fx1 * FORCESCALE  );
        PMEAccumulator      ify1   = llrint(  fy1 * FORCESCALE  );
        PMEAccumulator      ifz1   = llrint(  fz1 * FORCESCALE  );
        PMEAccumulator      ifx2   = llrint(  fx2 * FORCESCALE  );
        PMEAccumulator      ify2   = llrint(  fy2 * FORCESCALE  );
        PMEAccumulator      ifz2   = llrint(  fz2 * FORCESCALE  );
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.y],
                  llitoulli( -ifx1 - ifx2 ));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.y],
                  llitoulli( -ify1 - ify2 ));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.y],
                  llitoulli( -ifz1 - ifz2 )   );
      }
      else {
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps  = cSim.pImageNMRCOMAngleCOMGrp[pos * 3 + 1];
#else
        int2 COMgrps  = cSim.pNMRCOMAngleCOMGrp[pos * 3 + 1];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#ifdef NMR_NEIGHBORLIST
          int2 COMatms  = cSim.pImageNMRCOMAngleCOM[ip];
#else
          int2 COMatms  = cSim.pNMRCOMAngleCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
            PMEDouble dcomdx = rmass / rmstot1.y;
            PMEDouble kfx1   = fx1 * dcomdx;
            PMEDouble kfy1   = fy1 * dcomdx;
            PMEDouble kfz1   = fz1 * dcomdx;
            PMEDouble kfx2   = fx2 * dcomdx;
            PMEDouble kfy2   = fy2 * dcomdx;
            PMEDouble kfz2   = fz2 * dcomdx;
            ifx1   = llrint(kfx1 * FORCESCALE);
            ify1   = llrint(kfy1 * FORCESCALE);
            ifz1   = llrint(kfz1 * FORCESCALE);
            ifx2   = llrint(kfx2 * FORCESCALE);
            ify2   = llrint(kfy2 * FORCESCALE);
            ifz2   = llrint(kfz2 * FORCESCALE);
            atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[j],
                      llitoulli( -ifx1 - ifx2 ));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[j],
                      llitoulli( -ify1 - ify2 ));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[j],
                      llitoulli( -ifz1 - ifz2 ));
          }
        }
      }
      //COM_3:
      if (atom2 > 0) {
        PMEAccumulator  ifx2   = llrint(  fx2 * FORCESCALE  );
        PMEAccumulator  ify2   = llrint(  fy2 * FORCESCALE  );
        PMEAccumulator   ifz2   = llrint(  fz2 * FORCESCALE  );
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom2],
                  llitoulli( ifx2 ));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom2],
                  llitoulli( ify2 ));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom2],
                  llitoulli( ifz2 ));
      }
      else {
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps  = cSim.pImageNMRCOMAngleCOMGrp[pos * 3 + 2];
#else
        int2 COMgrps  = cSim.pNMRCOMAngleCOMGrp[pos * 3 + 2];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#ifdef NMR_NEIGHBORLIST
          int2 COMatms        = cSim.pImageNMRCOMAngleCOM[ip];
#else
          int2 COMatms        = cSim.pNMRCOMAngleCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
            PMEDouble dcomdx = rmass / rmstot2;
            PMEDouble kfx2   = fx2 * dcomdx;
            PMEDouble kfy2   = fy2 * dcomdx;
            PMEDouble kfz2   = fz2 * dcomdx;
            ifx2   = llrint(kfx2 * FORCESCALE);
            ify2   = llrint(kfy2 * FORCESCALE);
            ifz2   = llrint(kfz2 * FORCESCALE);
            atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[j],
                      llitoulli(ifx2));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[j],
                      llitoulli(ify2));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[j],
                      llitoulli(ifz2));
          }
        }
      }
    }   // end of " if ( pos < cSim.NMRCOMAngles ) " at line 1410
  exit6: {}
# ifdef NMR_ENERGY
    eCOMangle += __SHFL(WARP_MASK, eCOMangle, threadIdx.x ^ 1);
    eCOMangle += __SHFL(WARP_MASK, eCOMangle, threadIdx.x ^ 2);
    eCOMangle += __SHFL(WARP_MASK, eCOMangle, threadIdx.x ^ 4);
    eCOMangle += __SHFL(WARP_MASK, eCOMangle, threadIdx.x ^ 8);
    eCOMangle += __SHFL(WARP_MASK, eCOMangle, threadIdx.x ^ 16);
#ifdef AMBER_PLATFORM_AMD_WARP64
    eCOMangle += __SHFL(WARP_MASK, eCOMangle, threadIdx.x ^ 32);
#endif
    // Write out energies
    if ((threadIdx.x & GRID_BITS_MASK) == 0) {
      atomicAdd(cSim.pENMRCOMAngle, llitoulli(eCOMangle)); 
    }
#endif
  }             // end of " else if ( pos < cSim.NMRCOMAngleOffset)

//Calculate COM Torsion restraints
  else if (pos >= cSim.NMRCOMAngleOffset && pos < cSim.NMRCOMTorsionOffset) {
#ifdef NMR_ENERGY
    PMEAccumulator eCOMtorsion       = 0;
#endif
    pos                             -= cSim.NMRCOMAngleOffset;
    if ( pos < cSim.NMRCOMTorsions) {
#ifdef NMR_NEIGHBORLIST
      int4 atom   = cSim.pImageNMRCOMTorsionID1[pos];
#else
      int4 atom  = cSim.pNMRCOMTorsionID1[pos];
#endif
      PMEDouble4 rmstot;
      rmstot.x = 0.0;
      rmstot.y = 0.0;
      rmstot.z = 0.0;
      rmstot.w = 0.0;
      PMEDouble2 R1R2;
      PMEDouble2 R3R4;
      PMEDouble2 K2K3;
#ifdef NMR_TIMED
      int2 Step                    = cSim.pNMRCOMTorsionStep[pos];
      int  Inc                     = cSim.pNMRCOMTorsionInc[pos];
      if ((step < Step.x) || ((step > Step.y) && (Step.y > 0))) {
        goto exit7;
      }
      // read timed data
      PMEDouble2 R1R2Slp           = cSim.pNMRCOMTorsionR1R2Slp[pos];
      PMEDouble2 R1R2Int           = cSim.pNMRCOMTorsionR1R2Int[pos];
      PMEDouble2 R3R4Slp           = cSim.pNMRCOMTorsionR3R4Slp[pos];
      PMEDouble2 R3R4Int           = cSim.pNMRCOMTorsionR3R4Int[pos];
      PMEDouble2 K2K3Slp           = cSim.pNMRCOMTorsionK2K3Slp[pos];
      PMEDouble2 K2K3Int           = cSim.pNMRCOMTorsionK2K3Int[pos];

      // Calculate increment
      double dstep                 = step - (double)((step - Step.x) % abs(Inc));

      // Calculate restraint values
      R1R2.x                       = dstep * R1R2Slp.x + R1R2Int.x;
      R1R2.y                       = dstep * R1R2Slp.y + R1R2Int.y;
      R3R4.x                       = dstep * R3R4Slp.x + R3R4Int.x;
      R3R4.y                       = dstep * R3R4Slp.y + R3R4Int.y;
      if (Inc > 0) {
        K2K3.x                   = dstep * K2K3Slp.x + K2K3Int.x;
        K2K3.y                   = dstep * K2K3Slp.y + K2K3Int.y;
      }
      else {
        int nstepu               = (step - Step.x) / abs(Inc);
        K2K3.x                   = K2K3Int.x * pow(K2K3Slp.x, nstepu);
        K2K3.y                   = K2K3Int.y * pow(K2K3Slp.y, nstepu);
      }
#else
      R1R2                         = cSim.pNMRCOMTorsionR1R2[pos];
      R3R4                         = cSim.pNMRCOMTorsionR3R4[pos];
      K2K3                         = cSim.pNMRCOMTorsionK2K3[pos];
#endif
      PMEDouble atomIX;
      PMEDouble atomIY;
      PMEDouble atomIZ;
      PMEDouble atomJX;
      PMEDouble atomJY;
      PMEDouble atomJZ;
      PMEDouble atomKX;
      PMEDouble atomKY;
      PMEDouble atomKZ;
      PMEDouble atomLX;
      PMEDouble atomLY;
      PMEDouble atomLZ;
#ifdef NODPTEXTURE
  #ifdef NMR_NEIGHBORLIST
      //COM_1: NL
      if ( atom.x > 0 ) {
         atomIX      = cSim.pImageX[atom.x];
         atomIY      = cSim.pImageY[atom.x];
         atomIZ      = cSim.pImageZ[atom.x];
      }
      else {
        PMEDouble xtot          = 0.0;
        PMEDouble ytot          = 0.0;
        PMEDouble ztot          = 0.0;
        rmstot.x                = 0.0;
        int2 COMgrps            = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4];
        for ( int ip = COMgrps.x; ip < COMgrps.y; ip++ ) {
          int2 COMatms      = cSim.pImageNMRCOMTorsionCOM[ip];
          for ( int j = COMatms.x; j < COMatms.y + 1; j++ ){
            PMEDouble rmass      = cSim.pImageMass[j];
            xtot                 = xtot + cSim.pImageX[j] * rmass;
            ytot                 = ytot + cSim.pImageY[j] * rmass;
            ztot                 = ztot + cSim.pImageZ[j] * rmass;
            rmstot.x             = rmstot.x + rmass;
          }
        }
        atomIX      = xtot / rmstot.x;
        atomIY      = ytot / rmstot.x;
        atomIZ      = ztot / rmstot.x;
      }
      //COM_2: NL
      if ( atom.y > 0 ) {
        atomJX      = cSim.pImageX[atom.y];
        atomJY      = cSim.pImageY[atom.y];
        atomJZ      = cSim.pImageZ[atom.y];
      }
      else {
        PMEDouble xtot          = 0.0;
        PMEDouble ytot          = 0.0;
        PMEDouble ztot          = 0.0;
        rmstot.y                = 0.0;
        int2 COMgrps            = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 1];
        for ( int ip = COMgrps.x; ip < COMgrps.y; ip++ )  {
          int2 COMatms      = cSim.pImageNMRCOMTorsionCOM[ip];
          for ( int j = COMatms.x; j < COMatms.y + 1; j++ ) {
            PMEDouble rmass      = cSim.pImageMass[j];
            xtot                 = xtot + cSim.pImageX[j] * rmass;
            ytot                 = ytot + cSim.pImageY[j] * rmass;
            ztot                 = ztot + cSim.pImageZ[j] * rmass;
            rmstot.y             = rmstot.y + rmass;
          }
        }
        atomJX      = xtot / rmstot.y;
        atomJY      = ytot / rmstot.y;
        atomJZ      = ztot / rmstot.y;
      }
      //COM_3: NL
      if ( atom.z > 0 ) {
        atomKX      = cSim.pImageX[atom.z];
        atomKY      = cSim.pImageY[atom.z];
        atomKZ      = cSim.pImageZ[atom.z];
      }
      else {
        PMEDouble xtot           = 0.0;
        PMEDouble ytot           = 0.0;
        PMEDouble ztot           = 0.0;
        rmstot.z                 = 0.0;
        int2 COMgrps             = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 2];
        for ( int ip = COMgrps.x; ip < COMgrps.y; ip++ ) {
          int2 COMatms       = cSim.pImageNMRCOMTorsionCOM[ip];
          for ( int j = COMatms.x; j < COMatms.y + 1; j++ ) {
            PMEDouble rmass       = cSim.pImageMass[j];
            xtot                  = xtot + cSim.pImageX[j] * rmass;
            ytot                  = ytot + cSim.pImageY[j] * rmass;
            ztot                  = ztot + cSim.pImageZ[j] * rmass;
            rmstot.z              = rmstot.z + rmass;
          }
        }
        atomKX       = xtot / rmstot.z;
        atomKY       = ytot / rmstot.z;
        atomKZ       = ztot / rmstot.z;
      }
      //COM_4: NL
      if ( atom.w > 0 ) {
        atomLX      = cSim.pImageX[atom.w];
        atomLY      = cSim.pImageY[atom.w];
        atomLZ      = cSim.pImageZ[atom.w];
      }
      else {
        PMEDouble xtot           = 0.0;
        PMEDouble ytot           = 0.0;
        PMEDouble ztot           = 0.0;
        rmstot.w                 = 0.0;
        int2 COMgrps             = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 3];
        for ( int ip = COMgrps.x; ip < COMgrps.y; ip++ ) {
          int2 COMatms       = cSim.pImageNMRCOMTorsionCOM[ip];
          for ( int j = COMatms.x; j < COMatms.y + 1; j++ ) {
            PMEDouble rmass       = cSim.pImageMass[j];
            xtot                  = xtot + cSim.pImageX[j] * rmass;
            ytot                  = ytot + cSim.pImageY[j] * rmass;
            ztot                  = ztot + cSim.pImageZ[j] * rmass;
            rmstot.w              = rmstot.w + rmass;
          }
        }
        atomLX       = xtot / rmstot.w;
        atomLY       = ytot / rmstot.w;
        atomLZ       = ztot / rmstot.w;
      }
#else
      //COM_1: noNL
      if ( atom.x > 0 ) {
        atomIX      = cSim.pAtomX[atom.x];
        atomIY      = cSim.pAtomY[atom.x];
        atomIZ      = cSim.pAtomZ[atom.x];
      }
      else {
        PMEDouble xtot          = 0.0;
        PMEDouble ytot          = 0.0;
        PMEDouble ztot          = 0.0;
        rmstot.x                = 0.0;
        int2 COMgrps            = cSim.pNMRCOMTorsionCOMGrp[pos * 4];

        for ( int ip = COMgrps.x; ip < COMgrps.y; ip++ ) {
          int2 COMatms      = cSim.pNMRCOMTorsionCOM[ip];
          for ( int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass      = cSim.pAtomMass[j];
            xtot                 = xtot + cSim.pAtomX[j] * rmass;
            ytot                 = ytot + cSim.pAtomY[j] * rmass;
            ztot                 = ztot + cSim.pAtomZ[j] * rmass;
            rmstot.x             = rmstot.x + rmass;
          }
        }
        atomIX      = xtot / rmstot.x;
        atomIY      = ytot / rmstot.x;
        atomIZ      = ztot / rmstot.x;
      }
      //COM_2: noNL
      if ( atom.y > 0 ) {
        atomJX      = cSim.pAtomX[atom.y];
        atomJY      = cSim.pAtomY[atom.y];
        atomJZ      = cSim.pAtomZ[atom.y];
      }
      else {
        PMEDouble xtot          = 0.0;
        PMEDouble ytot          = 0.0;
        PMEDouble ztot          = 0.0;
        rmstot.y                = 0.0;
        int2 COMgrps            = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 1];
        for ( int ip = COMgrps.x; ip < COMgrps.y; ip++) {
          int2 COMatms      = cSim.pNMRCOMTorsionCOM[ip];
          for ( int j = COMatms.x; j < COMatms.y + 1; j++) {
            PMEDouble rmass      = cSim.pAtomMass[j];
            xtot                 = xtot + cSim.pAtomX[j] * rmass;
            ytot                 = ytot + cSim.pAtomY[j] * rmass;
            ztot                 = ztot + cSim.pAtomZ[j] * rmass;
            rmstot.y             = rmstot.y + rmass;
          }
        }
        atomJX      = xtot / rmstot.y;
        atomJY      = ytot / rmstot.y;
        atomJZ      = ztot / rmstot.y;
      }
      //COM_3: noNL
      if ( atom.z > 0 ) {
        atomKX      = cSim.pAtomX[atom.z];
        atomKY      = cSim.pAtomY[atom.z];
        atomKZ      = cSim.pAtomZ[atom.z];
      }
      else {
        PMEDouble xtot           = 0.0;
        PMEDouble ytot           = 0.0;
        PMEDouble ztot           = 0.0;
        rmstot.z                 = 0.0;
        int2 COMgrps             = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 2];
        for ( int ip = COMgrps.x; ip < COMgrps.y; ip++ ) {
          int2 COMatms       = cSim.pNMRCOMTorsionCOM[ip];
          for ( int j = COMatms.x; j < COMatms.y + 1; j++ ) {
            PMEDouble rmass       = cSim.pAtomMass[j];
            xtot                  = xtot + cSim.pAtomX[j] * rmass;
            ytot                  = ytot + cSim.pAtomY[j] * rmass;
            ztot                  = ztot + cSim.pAtomZ[j] * rmass;
            rmstot.z              = rmstot.z + rmass;
          }
        }
        atomKX       = xtot / rmstot.z;
        atomKY       = ytot / rmstot.z;
        atomKZ       = ztot / rmstot.z;
      }
      //COM_4: noNL
      if ( atom.w > 0 ) {
        atomLX      = cSim.pAtomX[atom.w];
        atomLY      = cSim.pAtomY[atom.w];
        atomLZ      = cSim.pAtomZ[atom.w];
      }
      else {
        PMEDouble xtot           = 0.0;
        PMEDouble ytot           = 0.0;
        PMEDouble ztot           = 0.0;
        rmstot.w                 = 0.0;
        int2 COMgrps             = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 3];
        for ( int ip = COMgrps.x; ip < COMgrps.y; ip++ ) {
          int2 COMatms       = cSim.pNMRCOMTorsionCOM[ip];
          for ( int j = COMatms.x; j < COMatms.y + 1; j++ ) {
            PMEDouble rmass      = cSim.pAtomMass[j];
            xtot                 = xtot + cSim.pAtomX[j] * rmass;
            ytot                 = ytot + cSim.pAtomY[j] * rmass;
            ztot                 = ztot + cSim.pAtomZ[j] * rmass;
            rmstot.w             = rmstot.w + rmass;
          }
        }
        atomLX       = xtot / rmstot.w;
        atomLY       = ytot / rmstot.w;
        atomLZ       = ztot / rmstot.w;
      }
#endif
#else
      //COM_1: SPXP
      if (atom.x > 0) {
#ifdef NMR_NEIGHBORLIST
        int2 iatomIX = tex1Dfetch<int2>(cSim.texImageX, atom.x);
        int2 iatomIY = tex1Dfetch<int2>(cSim.texImageX, atom.x + cSim.stride);
        int2 iatomIZ = tex1Dfetch<int2>(cSim.texImageX, atom.x + cSim.stride2);
#else
        int2 iatomIX = tex1Dfetch<int2>(cSim.texAtomX, atom.x);
        int2 iatomIY = tex1Dfetch<int2>(cSim.texAtomX, atom.x + cSim.stride);
        int2 iatomIZ = tex1Dfetch<int2>(cSim.texAtomX, atom.x + cSim.stride2);
#endif
        atomIX       = __hiloint2double(iatomIX.y, iatomIX.x);
        atomIY       = __hiloint2double(iatomIY.y, iatomIY.x);
        atomIZ       = __hiloint2double(iatomIZ.y, iatomIZ.x);
      }
      else {
        PMEDouble xtot           = 0.0;
        PMEDouble ytot           = 0.0;
        PMEDouble ztot           = 0.0;
        rmstot.x                 = 0.0;
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps             = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4];
#else
        int2 COMgrps             = cSim.pNMRCOMTorsionCOMGrp[pos * 4];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {

#ifdef NMR_NEIGHBORLIST
          int2 COMatms        = cSim.pImageNMRCOMTorsionCOM[ip];
#else
          int2 COMatms        = cSim.pNMRCOMTorsionCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++ ) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
#ifdef NMR_NEIGHBORLIST
            int2 iatomIX     = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomIY     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomIZ     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
#else
            int2 iatomIX     = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomIY     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomIZ     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
#endif

            atomIX           = __hiloint2double(iatomIX.y, iatomIX.x);
            atomIY           = __hiloint2double(iatomIY.y, iatomIY.x);
            atomIZ           = __hiloint2double(iatomIZ.y, iatomIZ.x);

            xtot             = xtot + atomIX * rmass;
            ytot             = ytot + atomIY * rmass;
            ztot             = ztot + atomIZ * rmass;
            rmstot.x         = rmstot.x + rmass;
          }
        }
        atomIX       = xtot / rmstot.x;
        atomIY       = ytot / rmstot.x;
        atomIZ       = ztot / rmstot.x;
      }
      //COM_2: noNL
      if (atom.y > 0) {
#ifdef NMR_NEIGHBORLIST
        int2 iatomJX = tex1Dfetch<int2>(cSim.texImageX, atom.y);
        int2 iatomJY = tex1Dfetch<int2>(cSim.texImageX, atom.y + cSim.stride);
        int2 iatomJZ = tex1Dfetch<int2>(cSim.texImageX, atom.y + cSim.stride2);
#else
        int2 iatomJX = tex1Dfetch<int2>(cSim.texAtomX, atom.y);
        int2 iatomJY = tex1Dfetch<int2>(cSim.texAtomX, atom.y + cSim.stride);
        int2 iatomJZ = tex1Dfetch<int2>(cSim.texAtomX, atom.y + cSim.stride2);
#endif
        atomJX       = __hiloint2double(iatomJX.y, iatomJX.x);
        atomJY       = __hiloint2double(iatomJY.y, iatomJY.x);
        atomJZ       = __hiloint2double(iatomJZ.y, iatomJZ.x);
      }
      else {
        PMEDouble xtot           = 0.0;
        PMEDouble ytot           = 0.0;
        PMEDouble ztot           = 0.0;
        rmstot.y                 = 0.0;
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps             = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 1];
#else
        int2 COMgrps             = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 1];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {

#ifdef NMR_NEIGHBORLIST
          int2 COMatms        = cSim.pImageNMRCOMTorsionCOM[ip];
#else
          int2 COMatms        = cSim.pNMRCOMTorsionCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++ ) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
#ifdef NMR_NEIGHBORLIST
            int2 iatomJX     = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomJY     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomJZ     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
#else
            int2 iatomJX     = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomJY     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomJZ     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
#endif

            atomJX           = __hiloint2double(iatomJX.y, iatomJX.x);
            atomJY           = __hiloint2double(iatomJY.y, iatomJY.x);
            atomJZ           = __hiloint2double(iatomJZ.y, iatomJZ.x);

            xtot             = xtot + atomJX * rmass;
            ytot             = ytot + atomJY * rmass;
            ztot             = ztot + atomJZ * rmass;
            rmstot.y         = rmstot.y + rmass;
          }
        }
        atomJX       = xtot / rmstot.y;
        atomJY       = ytot / rmstot.y;
        atomJZ       = ztot / rmstot.y;
      }
      //COM_3: SPXP
      if (atom.z > 0) {
#ifdef NMR_NEIGHBORLIST
        int2 iatomKX = tex1Dfetch<int2>(cSim.texImageX, atom.z);
        int2 iatomKY = tex1Dfetch<int2>(cSim.texImageX, atom.z + cSim.stride);
        int2 iatomKZ = tex1Dfetch<int2>(cSim.texImageX, atom.z + cSim.stride2);
#else
        int2 iatomKX = tex1Dfetch<int2>(cSim.texAtomX, atom.z);
        int2 iatomKY = tex1Dfetch<int2>(cSim.texAtomX, atom.z + cSim.stride);
        int2 iatomKZ = tex1Dfetch<int2>(cSim.texAtomX, atom.z + cSim.stride2);
#endif
        atomKX       = __hiloint2double(iatomKX.y, iatomKX.x);
        atomKY       = __hiloint2double(iatomKY.y, iatomKY.x);
        atomKZ       = __hiloint2double(iatomKZ.y, iatomKZ.x);
      }
      else {

        PMEDouble xtot           = 0.0;
        PMEDouble ytot           = 0.0;
        PMEDouble ztot           = 0.0;
        rmstot.z                 = 0.0;
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps             = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 2];
#else
        int2 COMgrps             = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 2];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {

#ifdef NMR_NEIGHBORLIST
          int2 COMatms        = cSim.pImageNMRCOMTorsionCOM[ip];
#else
          int2 COMatms        = cSim.pNMRCOMTorsionCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++ ) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
#ifdef NMR_NEIGHBORLIST
            int2 iatomKX     = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomKY     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomKZ     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
#else
            int2 iatomKX     = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomKY     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomKZ     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
#endif
            atomKX           = __hiloint2double(iatomKX.y, iatomKX.x);
            atomKY           = __hiloint2double(iatomKY.y, iatomKY.x);
            atomKZ           = __hiloint2double(iatomKZ.y, iatomKZ.x);
            xtot             = xtot + atomKX * rmass;
            ytot             = ytot + atomKY * rmass;
            ztot             = ztot + atomKZ * rmass;
            rmstot.z         = rmstot.z + rmass;
          }
        }
        atomKX       = xtot / rmstot.z;
        atomKY       = ytot / rmstot.z;
        atomKZ       = ztot / rmstot.z;
      }
      //COM_4: SPXP
      if (atom.w > 0) {
#ifdef NMR_NEIGHBORLIST
        int2 iatomLX = tex1Dfetch<int2>(cSim.texImageX, atom.w);
        int2 iatomLY = tex1Dfetch<int2>(cSim.texImageX, atom.w + cSim.stride);
        int2 iatomLZ = tex1Dfetch<int2>(cSim.texImageX, atom.w + cSim.stride2);
#else
        int2 iatomLX = tex1Dfetch<int2>(cSim.texAtomX, atom.w);
        int2 iatomLY = tex1Dfetch<int2>(cSim.texAtomX, atom.w + cSim.stride);
        int2 iatomLZ = tex1Dfetch<int2>(cSim.texAtomX, atom.w + cSim.stride2);
#endif
        atomLX       = __hiloint2double(iatomLX.y, iatomLX.x);
        atomLY       = __hiloint2double(iatomLY.y, iatomLY.x);
        atomLZ       = __hiloint2double(iatomLZ.y, iatomLZ.x);
      }
      else {
        PMEDouble xtot           = 0.0;
        PMEDouble ytot           = 0.0;
        PMEDouble ztot           = 0.0;
        rmstot.w                 = 0.0;
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps             = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 3];
#else
        int2 COMgrps             = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 3];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#ifdef NMR_NEIGHBORLIST
          int2 COMatms        = cSim.pImageNMRCOMTorsionCOM[ip];
#else
          int2 COMatms        = cSim.pNMRCOMTorsionCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++ ) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
#ifdef NMR_NEIGHBORLIST
            int2 iatomLX     = tex1Dfetch<int2>(cSim.texImageX, j);
            int2 iatomLY     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride);
            int2 iatomLZ     = tex1Dfetch<int2>(cSim.texImageX, j + cSim.stride2);
#else
            int2 iatomLX     = tex1Dfetch<int2>(cSim.texAtomX, j);
            int2 iatomLY     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride);
            int2 iatomLZ     = tex1Dfetch<int2>(cSim.texAtomX, j + cSim.stride2);
#endif
            atomLX           = __hiloint2double(iatomLX.y, iatomLX.x);
            atomLY           = __hiloint2double(iatomLY.y, iatomLY.x);
            atomLZ           = __hiloint2double(iatomLZ.y, iatomLZ.x);

            xtot             = xtot + atomLX * rmass;
            ytot             = ytot + atomLY * rmass;
            ztot             = ztot + atomLZ * rmass;
            rmstot.w         = rmstot.w + rmass;
          }
        }
        atomLX       = xtot / rmstot.w;
        atomLY       = ytot / rmstot.w;
        atomLZ       = ztot / rmstot.w;
      }
#endif
      // calculate the torsion:
      PMEDouble xij    = atomIX - atomJX;
      PMEDouble xkj    = atomKX - atomJX;
      PMEDouble xkl    = atomKX - atomLX;
      PMEDouble yij    = atomIY - atomJY;
      PMEDouble ykj    = atomKY - atomJY;
      PMEDouble ykl    = atomKY - atomLY;
      PMEDouble zij    = atomIZ - atomJZ;
      PMEDouble zkj    = atomKZ - atomJZ;
      PMEDouble zkl    = atomKZ - atomLZ;
     // calculate ij .cross. jk, and kl .cross. jk
      PMEDouble dx     =  yij * zkj  -  zij * ykj;
      PMEDouble dy     =  zij * xkj  -  xij * zkj;
      PMEDouble dz     =  xij * ykj  -  yij * xkj;
      PMEDouble gx     =  zkj * ykl  -  ykj * zkl;
      PMEDouble gy     =  xkj * zkl  -  zkj * xkl;
      PMEDouble gz     =  ykj * xkl  -  xkj * ykl;
      // calculate the magnitude of above vectors and their dot product
      PMEDouble fxi    =  dx * dx  +  dy * dy  +  dz * dz + tm24;
      PMEDouble fyi    =  gx * gx  +  gy * gy  +  gz * gz + tm24;
      PMEDouble ct     =  dx * gx  +  dy * gy  +  dz * gz;
     // branch if linear torsion
      PMEDouble z1                 = rsqrt( fxi );
      PMEDouble z2                 = rsqrt( fyi );
      PMEDouble z11                = z1 * z1;
      PMEDouble z22                = z2 * z2;
      PMEDouble z12                = z1 * z2;
      ct                          *= z1 * z2;
      // Cannot let ct be exactly 1. or -1., because this results
      // in sin(tau) = 0., and we divide by sin(tau) [sphi] in calculating
      //  the derivatives
      ct                           = max((PMEDouble)-0.9999999999999, min(ct, (PMEDouble)0.9999999999999));
      PMEDouble ap                 = acos(ct);
      PMEDouble s                  = xkj*(dz * gy - dy * gz)  +  ykj*(dx * gz - dz * gx)  +  zkj*(dy * gx - dx * gy);
      if ( s < (PMEDouble)0.0 ){
        ap = -ap;
      }
      ap                           = PI - ap;
      PMEDouble sphi, cphi;
      faster_sincos(ap, &sphi, &cphi);

      // reflect the handedness of the torsion to get the SMALLEST angle...
      PMEDouble apmean = ( R1R2.y + R3R4.x ) * (PMEDouble)0.5;
      if ( ap - apmean > PI ){
        ap  -= (PMEDouble)2.0 * (PMEDouble)PI;
      }
      if ( apmean - ap > PI ) {
        ap                     += (PMEDouble)2.0 * (PMEDouble)PI;
      }
      PMEDouble df;
#ifdef NMR_ENERGY
      PMEDouble e;
#endif
      if ( ap < R1R2.x) {
        PMEDouble dif            = R1R2.x - R1R2.y;
        df                       = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
        e                        = df * (ap - R1R2.x) + K2K3.x * dif * dif;
#endif
      }
      else if ( ap < R1R2.y) {
        PMEDouble dif            =  ap - R1R2.y;
        df                       = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
        e                        = K2K3.x * dif * dif;
#endif
      }
      else if ( ap <=  R3R4.x) {
        df                       = (PMEDouble)0.0;
#ifdef NMR_ENERGY
        e                        = (PMEDouble)0.0;
#endif
      }
      else if ( ap < R3R4.y) {
        PMEDouble dif            =  ap - R3R4.x;
        df                       = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
        e                        = K2K3.y * dif * dif;
#endif
      }
      else {
        PMEDouble dif            = R3R4.y - R3R4.x;
        df                       = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
        e                        = df * (ap - R3R4.y) + K2K3.y * dif * dif;
#endif
      }
      df                              *= -1.0 / sphi;
       // Steered MD for COMTorsion... para Adrian...
      if (cSim.bJar) {
        double fold              = cSim.pNMRJarData[2];
        double work              = cSim.pNMRJarData[3];
        double first             = cSim.pNMRJarData[4];
        double fcurr             = (PMEDouble)-2.0 * K2K3.x * (ap - R1R2.y);
        if (first == (PMEDouble)0.0) {
           fold                 = -fcurr;
           cSim.pNMRJarData[4]  = (PMEDouble)1.0;
        }
        work                    += (PMEDouble)0.5 * (fcurr + fold) * cSim.drjar;
        cSim.pNMRJarData[0]      = R1R2.y;
        cSim.pNMRJarData[1]      = ap;
        cSim.pNMRJarData[2]      = fcurr;
        cSim.pNMRJarData[3]      = work;
      }
#ifdef NMR_ENERGY
      eCOMtorsion      += llrint(ENERGYSCALE * e);
#endif
      PMEDouble dcdx   =-gx * z12 - cphi * dx * z11;
      PMEDouble dcdy   =-gy * z12 - cphi * dy * z11;
      PMEDouble dcdz   =-gz * z12 - cphi * dz * z11;
      PMEDouble dcgx   = dx * z12 + cphi * gx * z22;
      PMEDouble dcgy   = dy * z12 + cphi * gy * z22;
      PMEDouble dcgz   = dz * z12 + cphi * gz * z22;
      PMEDouble dr1    = df * ( dcdz * ykj - dcdy * zkj);
      PMEDouble dr2    = df * ( dcdx * zkj - dcdz * xkj);
      PMEDouble dr3    = df * ( dcdy * xkj - dcdx * ykj);
      PMEDouble dr4    = df * ( dcgz * ykj - dcgy * zkj);
      PMEDouble dr5    = df * ( dcgx * zkj - dcgz * xkj);
      PMEDouble dr6    = df * ( dcgy * xkj - dcgx * ykj);
      PMEDouble drx    = df * (-dcdy * zij + dcdz * yij + dcgy * zkl -  dcgz * ykl);
      PMEDouble dry    = df * ( dcdx * zij - dcdz * xij - dcgx * zkl +  dcgz * xkl);
      PMEDouble drz    = df * (-dcdx * yij + dcdy * xij + dcgx * ykl -  dcgy * xkl);
      // declare COM stuff
      PMEAccumulator idr1;
      PMEAccumulator idr2;
      PMEAccumulator idr3;
      PMEAccumulator idr4;
      PMEAccumulator idr5;
      PMEAccumulator idr6;
      PMEAccumulator idrx;
      PMEAccumulator idry;
      PMEAccumulator idrz;
      //COM_1
      if (atom.x > 0) {
        idr1           = llrint(dr1 * FORCESCALE);
        idr2           = llrint(dr2 * FORCESCALE);
        idr3           = llrint(dr3 * FORCESCALE);
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.x],
                  llitoulli(-idr1));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.x],
                  llitoulli(-idr2));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.x],
                  llitoulli(-idr3));
      }
      else {
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps             = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4];
#else
        int2 COMgrps             = cSim.pNMRCOMTorsionCOMGrp[pos * 4];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#ifdef NMR_NEIGHBORLIST
          int2 COMatms        = cSim.pImageNMRCOMTorsionCOM[ip];
#else
          int2 COMatms        = cSim.pNMRCOMTorsionCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
            PMEDouble dcomdx = rmass / rmstot.x;
            PMEDouble kdr1   = dr1 * dcomdx;
            PMEDouble kdr2   = dr2 * dcomdx;
            PMEDouble kdr3   = dr3 * dcomdx;
            idr1             = llrint(kdr1 * FORCESCALE);
            idr2             = llrint(kdr2 * FORCESCALE);
            idr3             = llrint(kdr3 * FORCESCALE);
            atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[j],
                      llitoulli(-idr1));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[j],
                      llitoulli(-idr2));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[j],
                      llitoulli(-idr3));
          }
        }
      }
      //COM_2
      if (atom.y > 0) {
        idr1           = llrint(dr1 * FORCESCALE);
        idr2           = llrint(dr2 * FORCESCALE);
        idr3           = llrint(dr3 * FORCESCALE);
        idrx           = llrint(drx * FORCESCALE);
        idry           = llrint(dry * FORCESCALE);
        idrz           = llrint(drz * FORCESCALE);
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.y],
                  llitoulli(-idrx + idr1));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.y],
                  llitoulli(-idry + idr2));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.y],
                  llitoulli(-idrz + idr3));
      }
      else {
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps             = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 1];
#else
        int2 COMgrps             = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 1];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#ifdef NMR_NEIGHBORLIST
          int2 COMatms        = cSim.pImageNMRCOMTorsionCOM[ip];
#else
          int2 COMatms        = cSim.pNMRCOMTorsionCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
            PMEDouble dcomdx = rmass / rmstot.y;
            PMEDouble kdr1   = dr1 * dcomdx;
            PMEDouble kdr2   = dr2 * dcomdx;
            PMEDouble kdr3   = dr3 * dcomdx;
            PMEDouble kdrx   = drx * dcomdx;
            PMEDouble kdry   = dry * dcomdx;
            PMEDouble kdrz   = drz * dcomdx;
            idr1             = llrint(kdr1 * FORCESCALE);
            idr2             = llrint(kdr2 * FORCESCALE);
            idr3             = llrint(kdr3 * FORCESCALE);
            idrx             = llrint(kdrx * FORCESCALE);
            idry             = llrint(kdry * FORCESCALE);
            idrz             = llrint(kdrz * FORCESCALE);
            atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[j],
                      llitoulli(-idrx + idr1));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[j],
                      llitoulli(-idry + idr2));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[j],
                      llitoulli(-idrz + idr3));
          }
        }
      }
      //COM_3
      if (atom.z > 0) {
        idr4           = llrint(dr4 * FORCESCALE);
        idr5           = llrint(dr5 * FORCESCALE);
        idr6           = llrint(dr6 * FORCESCALE);
        idrx           = llrint(drx * FORCESCALE);
        idry           = llrint(dry * FORCESCALE);
        idrz           = llrint(drz * FORCESCALE);
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.z], llitoulli(idrx + idr4));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.z], llitoulli(idry + idr5));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.z], llitoulli(idrz + idr6));
      }
      else {
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps             = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 2];
#else
        int2 COMgrps             = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 2];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#ifdef NMR_NEIGHBORLIST
          int2 COMatms        = cSim.pImageNMRCOMTorsionCOM[ip];
#else
          int2 COMatms        = cSim.pNMRCOMTorsionCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
            PMEDouble dcomdx = rmass / rmstot.z;
            PMEDouble kdr4   = dr4 * dcomdx;
            PMEDouble kdr5   = dr5 * dcomdx;
            PMEDouble kdr6   = dr6 * dcomdx;
            PMEDouble kdrx   = drx * dcomdx;
            PMEDouble kdry   = dry * dcomdx;
            PMEDouble kdrz   = drz * dcomdx;
            idr4             = llrint(kdr4 * FORCESCALE);
            idr5             = llrint(kdr5 * FORCESCALE);
            idr6             = llrint(kdr6 * FORCESCALE);
            idrx             = llrint(kdrx * FORCESCALE);
            idry             = llrint(kdry * FORCESCALE);
            idrz             = llrint(kdrz * FORCESCALE);
            atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[j],
                      llitoulli(idrx + idr4));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[j],
                      llitoulli(idry + idr5));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[j],
                      llitoulli(idrz + idr6));
          }
        }
      }
      //COM_4
      if (atom.w > 0) {
        idr4           = llrint(dr4 * FORCESCALE);
        idr5           = llrint(dr5 * FORCESCALE);
        idr6           = llrint(dr6 * FORCESCALE);
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.w],
                  llitoulli(-idr4));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.w],
                  llitoulli(-idr5));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.w],
                  llitoulli(-idr6));
      }
      else {
#ifdef NMR_NEIGHBORLIST
        int2 COMgrps             = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 3];
#else
        int2 COMgrps             = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 3];
#endif
        for (int ip = COMgrps.x; ip < COMgrps.y; ip++) {
#ifdef NMR_NEIGHBORLIST
          int2 COMatms        = cSim.pImageNMRCOMTorsionCOM[ip];
#else
          int2 COMatms        = cSim.pNMRCOMTorsionCOM[ip];
#endif
          for (int j = COMatms.x; j < COMatms.y + 1; j++) {
#ifdef NMR_NEIGHBORLIST
            PMEDouble rmass  = cSim.pImageMass[j];
#else
            PMEDouble rmass  = cSim.pAtomMass[j];
#endif
            PMEDouble dcomdx = rmass / rmstot.z;
            PMEDouble kdr4   = dr4 * dcomdx;
            PMEDouble kdr5   = dr5 * dcomdx;
            PMEDouble kdr6   = dr6 * dcomdx;
            idr4             = llrint(kdr4 * FORCESCALE);
            idr5             = llrint(kdr5 * FORCESCALE);
            idr6             = llrint(kdr6 * FORCESCALE);
            atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[j],
                      llitoulli(-idr4));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[j],
                      llitoulli(-idr5));
            atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[j],
                      llitoulli(-idr6));
          }
        }
      }
    }   // end of " if ( pos < cSim.NMRCOMTorsions ) "
    exit7: {}
#ifdef NMR_ENERGY
    eCOMtorsion += __SHFL(WARP_MASK, eCOMtorsion, threadIdx.x ^ 1);
    eCOMtorsion += __SHFL(WARP_MASK, eCOMtorsion, threadIdx.x ^ 2);
    eCOMtorsion += __SHFL(WARP_MASK, eCOMtorsion, threadIdx.x ^ 4);
    eCOMtorsion += __SHFL(WARP_MASK, eCOMtorsion, threadIdx.x ^ 8);
    eCOMtorsion += __SHFL(WARP_MASK, eCOMtorsion, threadIdx.x ^ 16);
#ifdef AMBER_PLATFORM_AMD_WARP64
    eCOMtorsion += __SHFL(WARP_MASK, eCOMtorsion, threadIdx.x ^ 32);
#endif
    // Write out energies
    if ((threadIdx.x & GRID_BITS_MASK) == 0) {
      atomicAdd(cSim.pENMRCOMTorsion, llitoulli(eCOMtorsion)); 
    }
#endif
  } // end of " else if ( pos < cSim.NMRCOMTorsionOffset )
}// END PROGRAM kCNF.h
