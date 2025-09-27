# Platform-specific info should be found in config.h

include ../../../config.h

CPP_OBJS = cuda_info.o gpu.o \
          base_gpuContext.o base_simulationConst.o gpuMemoryInfo.o \
          gputypes.o gpuContext.o \
          gti_simulationConst.o gti_gpuContext.o \
          gti_gpu.o gti_f95.o gti_cuda.o \
          gti_schdule_functions.o \
          base_xrayHostContext.o base_xrayDevConstants.o mdl_xray.o \
          matrix.o bondRemap.o

CU_OBJS = kForcesUpdate.o kCalculateLocalForces.o kCalculateGBBornRadii.o \
          kCalculatePMENonbondEnergy.o kCalculateGBNonbondEnergy1.o kNLRadixSort.o \
          kCalculateGBNonbondEnergy2.o kShake.o kNeighborList.o kPMEInterpolation.o \
          kCalculateEFieldEnergy.o kDataTransfer.o kCalculateNEBForces.o kXray.o \
          kCalculateCoarseGridEnergy.o

GTI_OBJS = gti_energies_kernels.o gti_general_kernels.o \
          gti_NBList_kernels.o gti_nonBond_kernels.o

.SUFFIXES: .F90 .cpp .cu .o

# CUDA compilations flags
kForcesUpdate.cu.CU_FLAGS                   =
kCalculateLocalForces.cu.CU_FLAGS           =
kCalculateGBBornRadii.cu.CU_FLAGS           =
kCalculateGBNonbondEnergy1.cu.CU_FLAGS      =
kCalculateGBNonbondEnergy2.cu.CU_FLAGS      =
kRandom.cu.CU_FLAGS                         =
kShake.cu.CU_FLAGS                          =
kPMEInterpolation.cu.CU_FLAGS               =
kNeighborList.cu.CU_FLAGS                   =
kCalculatePMENonbondEnergy.cu.CU_FLAGS      =
kCalculateAMDWeights.cu.CU_FLAGS            =
kCalculateGAMDWeights.cu.CU_FLAGS           =
kCalculateEFieldEnergy.cu.CU_FLAGS          =
kCalculateCoarseGridEnergy.cu.CU_FLAGS          =

cuda.o: $(CU_OBJS) $(CPP_OBJS) $(GTI_OBJS)
	$(NVCC) -dlink -o cuda.o $(CPP_OBJS) $(CU_OBJS) $(GTI_OBJS)

cuda.a: cuda.o
	ar rvs cuda.a cuda.o  $(CPP_OBJS) $(CU_OBJS) $(GTI_OBJS)

cuda.so: $(CU_OBJS) $(CPP_OBJS) $(GTI_OBJS)
ifeq (,$(findstring nvcc,$(NVCC)))
	$(NVCC) --hip-link -fgpu-rdc -shared -o cuda.so $(CPP_OBJS) $(CU_OBJS) $(GTI_OBJS)
else
	$(NVCC) -shared -o cuda.o $(CPP_OBJS) $(CU_OBJS) $(GTI_OBJS)
endif

.F90.o:
	$(PMEMD_F90) -fPIC $(AMBERBUILDFLAGS) $(PREC_MODEL) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(PMEMD_CU_INCLUDES) -c $*.F90

.cpp.o:
	$(PMEMD_CC) -fPIC $(AMBERBUILDFLAGS) $(PREC_MODEL) $(PMEMD_COPTFLAGS) $(PMEMD_CU_DEFINES) $(PMEMD_CU_INCLUDES) -c $*.cpp

.cu.o:  gputypes.h
ifeq (,$(findstring nvcc,$(NVCC)))
	$(NVCC) -fgpu-rdc -D __HIPCC_RDC__ $(EXTRAFLAGS) $(AMBERBUILDFLAGS) $(PREC_MODEL) $(PMEMD_CU_DEFINES) $(PMEMD_CU_INCLUDES) $(PMEMD_CU_EXTRA) $($*.CU_FLAGS) -c $*.cu
else
	$(NVCC) $(EXTRAFLAGS) $(AMBERBUILDFLAGS) $(PREC_MODEL) $(PMEMD_CU_DEFINES) $(PMEMD_CU_INCLUDES) $(PMEMD_CU_EXTRA) $($*.CU_FLAGS) -c $*.cu
endif

clean:
	rm -f *.o *.linkinfo cuda.a *.mod


## CUDA dependencies
kNeighborList.o:	kBNL.h kReImageCoord.h
kForcesUpdate.o:	kNTPKernels.h kNTPCalls.h kU.h kRandom.h kLocalToGlobalLoop.h \
                        kLocalToGlobalKernel.h kOrientForcesLoop.h kOrientForcesKernel.h \
                        kRV.h kAFECalls.h kAFEKernels.h kSGLD.h
kCalculatePMENonbondEnergy.o:	kNLCPNE.h kNLCINE.h
kCalculateGBNonbondEnergy1.o:	kCalculateGBNonbondEnergy1.h
kCalculateGBNonbondEnergy2.o:	kCalculateGBNonbondEnergy2.h
kCalculateGBBornRadii.o:	kCalculateGBBornRadii.h
kCalculateLocalForces.o:	kCNF.h kBWU.h \
                                kBWU_bond.h kBWU_angl.h kBWU_dihe.h kBWU_nb14.h kBWU_cnst.h \
                                kBWU_cmap.h kBWU_cimp.h kBWU_urey.h kBWU_nmr2.h kBWU_nmr3.h \
                                kBWU_nmr4.h kBWU_EnergyReduction.h
kPMEInterpolation.o:	kReImageCoord.h kPGS.h kPSSE.h
kCalculateEFieldEnergy.o: kEFE.h
kShake.o:	kShake.h kRattle.h
kDataTransfer.o:                kShuttle.h
kCalculateNEBForces.o:                 kShuttle.h kNEB.h kRMSFit.h kNorm.h
$(CU_OBJS):	gpu.h gputypes.h simulationConst.h
gti_cuda.o:     gti_cuda.cuh gti_general_kernels.cuh gti_NBList_kernels.cuh gti_nonBond_kernels.cuh gti_energies_kernels.cuh
gti_energies_kernels.o: gti_simulationConst.h gti_kernelAttribute.cuh gti_energies_kernels.cuh gti_utils.cuh gti_const.cuh
gti_general_kernels.o:  gti_simulationConst.h gti_kernelAttribute.cuh gti_general_kernels.cuh gti_const.cuh
gti_nonBond_kernels.o:  gti_simulationConst.h gti_kernelAttribute.cuh gti_def.h gti_utils.cuh ptxmacros.h gti_nonBond_kernels.cuh gti_const.cuh
gti_NBList_kernels.o: gti_simulationConst.h gti_kernelAttribute.cuh gti_utils.cuh gti_NBList_kernels.cuh
gti_gpuContext.o: gti_gpuContext.h gti_def.h
gti_simulationConst.o: gti_simulationConst.h
gti_simulationConst.h: gti_def.h
mdl_xray.o:	base_xrayDevConstants.h base_xrayHostContext.h
kXray.o:	mdl_xray.h
