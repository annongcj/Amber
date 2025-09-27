project(emil CXX)

set(EMIL_SOURCE amberInterface.C 
		amberNonbonds.C 
		defineTopo.C 
		forces.C 
        initialise.C 
		inout.C 
        jonkerVolgenant.C 
		lap.C 
		linearAssignment.C 
		mcMoves.C 
		mtrand.C 
		relocMove.C 
		reportFreeEnergy.C 
		setup.C 
		softForce.C 
		stepLambda.C 
		swapMove.C 
		tools.C 
		wellMove.C)

include_directories(.)
add_compile_options(${OPT_CXXFLAGS})
add_library(emil ${EMIL_SOURCE})
install_libraries(emil)

if(MPI)
	make_mpi_version(emil emil_mpi LANGUAGES CXX INSTALL)
	target_compile_definitions(emil_mpi PRIVATE USE_MPI)
endif()