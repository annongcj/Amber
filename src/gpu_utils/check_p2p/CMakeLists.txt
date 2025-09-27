#CMAKE build file for p2p check
set(P2P_CHECK_SOURCES gpuP2PCheck.cu)

#------------------------------------------------------------------------------------------

include_directories(.) #needed for CUDA compile
cuda_add_executable(gpuP2PCheck gpuP2PCheck.cu)

install(TARGETS gpuP2PCheck DESTINATION ${BINDIR} COMPONENT pmemd_CUDA)