set(PMEMD_CUDA_CUDA_SOURCES
	kForcesUpdate.cu kCalculateLocalForces.cu kCalculateGBBornRadii.cu
	kCalculatePMENonbondEnergy.cu kCalculateGBNonbondEnergy1.cu kNLRadixSort.cu
	kCalculateGBNonbondEnergy2.cu kShake.cu kNeighborList.cu kPMEInterpolation.cu
	kCalculateEFieldEnergy.cu kDataTransfer.cu kCalculateNEBForces.cu
    kXray.cu kCalculateCoarseGridEnergy.cu )

set(PMEMD_CUDA_GTI_SOURCES
	gti_energies_kernels.cu gti_general_kernels.cu
	gti_NBList_kernels.cu gti_nonBond_kernels.cu gti_PME_kernels.cu
	gti_cuda.cu)

set(PMEMD_CUDA_FORTRAN_SOURCES cuda_info.F90)

set(PMEMD_CUDA_CXX_SOURCES gpu.cpp
	base_gpuContext.cpp base_simulationConst.cpp gpuMemoryInfo.cpp
	gti_gpuContext.cpp gti_simulationConst.cpp
	gputypes.cpp gpuContext.cpp
	gti_gpu.cpp gti_f95.cpp
    base_xrayHostContext.cpp base_xrayDevConstants.cpp mdl_xray.cpp
	matrix.cpp bondRemap.cpp gti_schdule_functions.cpp)


# build settings
# --------------------------------------------------------------------
#NOTE: we don't (yet) have PMEMD_CXXFLAGS, so we just use the OPT_CXXFLAGS
set_property(SOURCE ${PMEMD_CUDA_CXX_SOURCES} PROPERTY COMPILE_FLAGS "${OPT_CXXFLAGS_SPC}")
set_property(SOURCE ${PMEMD_CUDA_CXX_SOURCES} PROPERTY COMPILE_DEFINITIONS ${PMEMD_C_DEFINITIONS})
set_property(SOURCE ${PMEMD_CUDA_CXX_SOURCES} PROPERTY INCLUDE_DIRECTORIES ${CUDA_INCLUDE_DIRS})

set_property(SOURCE ${PMEMD_CUDA_FORTRAN_SOURCES} PROPERTY COMPILE_FLAGS "${PMEMD_FFLAGS_SPC}")
set_property(SOURCE ${PMEMD_CUDA_FORTRAN_SOURCES} PROPERTY COMPILE_DEFINITIONS ${PMEMD_F_DEFINITIONS})

add_definitions(${PMEMD_CUDA_DEFINES})

#cuda libraries
#------------------------------------------------------------------------------------------
foreach(PRECISION ${PMEMD_CUDA_PRECISIONS})
	set(ADDITIONAL_FLAGS "")
	# HIP-TODO: Enable -ffast-math for DPFP when it's proven that the compiler enables
	# correct optimizations for DP (set it as PMEMD_NVCC_FLAGS in PMEMDCompilerFlags.cmake)
	if(HIP AND ${PRECISION} STREQUAL "SPFP")
		list(APPEND ADDITIONAL_FLAGS -ffast-math)
	endif()
	set(LIBRARY_TYPE STATIC)
	if(HIP AND HIP_RDC)
		# HIP-TODO: hipcc does not define something like __CUDACC_RDC__, so we need to do this
		# manually. Remove when a builtin alternative is added to the compiler.
		list(APPEND ADDITIONAL_FLAGS -fgpu-rdc -D __HIPCC_RDC__)
		set(LIBRARY_TYPE SHARED)
	endif()
	cuda_add_library(pmemd_cuda_${PRECISION} ${LIBRARY_TYPE} ${PMEMD_CUDA_CUDA_SOURCES} ${PMEMD_CUDA_GTI_SOURCES} ${PMEMD_CUDA_FORTRAN_SOURCES} ${PMEMD_CUDA_CXX_SOURCES}
		OPTIONS ${PMEMD_NVCC_FLAGS} ${ADDITIONAL_FLAGS} -Duse_${PRECISION})
	if(HIP AND HIP_RDC)
		set_property(TARGET pmemd_cuda_${PRECISION} PROPERTY LINKER_LANGUAGE HIP)
		install(TARGETS pmemd_cuda_${PRECISION} DESTINATION ${LIBDIR} COMPONENT pmemd_CUDA)
	endif()

	#pass the precision define to the non-CUDA sources
	target_compile_definitions(pmemd_cuda_${PRECISION} PRIVATE -Duse_${PRECISION})
	target_link_libraries(pmemd_cuda_${PRECISION} curand cublas)

	if(NOT HIP)
		target_link_libraries(pmemd_cuda_${PRECISION} cusolver)
	endif()

	if(VKFFT)
		if(NOT HIP)
			target_link_libraries(pmemd_cuda_${PRECISION} cuda cudadevrt nvrtc)
		endif()
	endif()

	if(NCCL)
		target_link_libraries(pmemd_cuda_${PRECISION} nccl)
	endif()

	if(MPI)
		cuda_include_directories(${MPI_CXX_INCLUDE_PATH})
		cuda_add_library(pmemd_cuda_${PRECISION}_mpi ${LIBRARY_TYPE} ${PMEMD_CUDA_CUDA_SOURCES} ${PMEMD_CUDA_GTI_SOURCES} ${PMEMD_CUDA_FORTRAN_SOURCES} ${PMEMD_CUDA_CXX_SOURCES}
			OPTIONS ${PMEMD_NVCC_FLAGS} ${ADDITIONAL_FLAGS} -Duse_${PRECISION} -DMPI)
		if(HIP AND HIP_RDC)
			set_property(TARGET pmemd_cuda_${PRECISION}_mpi PROPERTY LINKER_LANGUAGE HIP)
			install(TARGETS pmemd_cuda_${PRECISION}_mpi DESTINATION ${LIBDIR} COMPONENT pmemd_CUDA)
            target_link_libraries(pmemd_cuda_${PRECISION}_mpi open-pal)
		endif()

		#pass the precision define to the non-CUDA sources
		target_compile_definitions(pmemd_cuda_${PRECISION}_mpi PRIVATE -Duse_${PRECISION})
		target_link_libraries(pmemd_cuda_${PRECISION}_mpi curand cublas cufft mpi_c mpi_cxx mpi_fortran)

		if(NOT HIP)
			target_link_libraries(pmemd_cuda_${PRECISION}_mpi cusolver)
		endif()

		if(NCCL)
			target_link_libraries(pmemd_cuda_${PRECISION}_mpi nccl)
		endif()
	endif()

endforeach()
