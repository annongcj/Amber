# lists Amber's components and their descriptions.
# Must be included immediately before Packaging.cmake

include(CPackComponent)

# component groups: "AmberTools" and "Amber"
# --------------------------------------------------------------------
cpack_add_component_group(AmberTools DESCRIPTION "The open-source AmberTools molecular dynamics package" EXPANDED)

if(NOT AMBERTOOLS_ONLY)
	cpack_add_component_group(Amber DESCRIPTION "The closed-source, paid portion of Amber.  Essentially just the PMEMD simulator program and its variants." EXPANDED)
endif()

# --------------------------------------------------------------------
# component list


# the following causes install() commands with no COMPONENT argument to get added to the "Serial" component
# NOTE: we use the default component for the Amber Serial programs and libraries.
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME Serial)

cpack_add_component(Serial DISPLAY_NAME Serial DESCRIPTION "Non-parallelized AmberTools programs and libraries" GROUP AmberTools DEPENDS Data)
cpack_add_component(Data DISPLAY_NAME Data DESCRIPTION "Force field files and various other datafiles used by Amber" GROUP AmberTools)

# global list of components to install. 
set(AMBER_INSTALL_COMPONENTS Data Serial)

if(BUILD_PYTHON)
	if(DOWNLOAD_MINICONDA)
		set(PYTHON_COMPONENT_DESC "Python components of Amber.  Includes a self-contained Miniconda runtime.")
	else()
		set(PYTHON_COMPONENT_DESC "Python components of Amber.  You will need a compatible system python to use them.")
	endif()
	
	cpack_add_component(Python DISPLAY_NAME Python DESCRIPTION ${PYTHON_COMPONENT_DESC} GROUP AmberTools DEPENDS Serial)
	list(APPEND AMBER_INSTALL_COMPONENTS Python)
endif()


# the tests and examples don't work on Windows, and they're too large anyway.
if(TARGET_WINDOWS AND INSTALL_TESTS)
	message(WARNING "Tests cannot be packaged on Windows because they exceed CMake's hard limit of 1GB uncompressed package size")
endif()

cpack_add_component(Tests DISPLAY_NAME Tests DESCRIPTION "Test scripts for Amber and AmberTools" GROUP AmberTools)
cpack_add_component(Examples DISPLAY_NAME Examples DESCRIPTION "Example usage scripts for Amber and AmberTools" GROUP AmberTools)
cpack_add_component(Benchmarks DISPLAY_NAME Benchmarks DESCRIPTION "Benchmarking scripts for Amber and AmberTools" GROUP AmberTools)
list(APPEND AMBER_INSTALL_COMPONENTS Tests Examples Benchmarks)

if(MPI)
	cpack_add_component(MPI DISPLAY_NAME MPI DESCRIPTION "Components of AmberTools that are parallelized with MPI." GROUP AmberTools DEPENDS Data)
	list(APPEND AMBER_INSTALL_COMPONENTS MPI)
endif()

if(OPENMP)
	cpack_add_component(OpenMP DISPLAY_NAME OpenMP DESCRIPTION "Components of AmberTools that are parallelized with OpenMP" GROUP AmberTools DEPENDS Data)
	list(APPEND AMBER_INSTALL_COMPONENTS OpenMP)
endif()

if(CUDA)
	cpack_add_component(CUDA DISPLAY_NAME CUDA DESCRIPTION "Components of AmberTools that are parallelized with NVIDIA CUDA" GROUP AmberTools DEPENDS Data)
	list(APPEND AMBER_INSTALL_COMPONENTS CUDA)
endif()

# pmemd components
if(NOT AMBERTOOLS_ONLY)

	cpack_add_component(pmemd DISPLAY_NAME pmemd DESCRIPTION "Serial version of pmemd, the leaner and meaner molecular simulator" GROUP Amber DEPENDS Data)
	list(APPEND AMBER_INSTALL_COMPONENTS pmemd)
	
	if(MIC)
		cpack_add_component(pmemd_MIC DISPLAY_NAME pmemd.MIC DESCRIPTION "MIC-parallelized version of pmemd." GROUP Amber DEPENDS Data)
		list(APPEND AMBER_INSTALL_COMPONENTS pmemd_MIC)
	endif()
	
	if(MPI)
		cpack_add_component(pmemd_MPI DISPLAY_NAME pmemd.MPI DESCRIPTION "MPI version of pmemd" GROUP Amber DEPENDS Data)
		list(APPEND AMBER_INSTALL_COMPONENTS pmemd_MPI)
	endif()
	
	if(CUDA)
		# the cuda.MPI versions of pmemd are lumped with the regular CUDA versions.  I figure that makes sense because you need special hardware for CUDA, but MPI can be run anywhere.
	
		if(MPI)
			set(PMEMD_CUDA_DESCRIPTION "CUDA and CUDA.MPI versions of pmemd")
		else()
			set(PMEMD_CUDA_DESCRIPTION "CUDA versions of pmemd")
		endif()
		
		cpack_add_component(pmemd_CUDA DISPLAY_NAME pmemd.CUDA DESCRIPTION ${PMEMD_CUDA_DESCRIPTION} GROUP Amber DEPENDS Data)
		list(APPEND AMBER_INSTALL_COMPONENTS pmemd_CUDA)
	endif()
		
endif()

# add "make install_<component>" targets for all components
# --------------------------------------------------------------------

foreach(COMPONENT ${AMBER_INSTALL_COMPONENTS})

	string(TOLOWER ${COMPONENT} COMPONENT_LCASE)
	
	# add "make install" target	
	ADD_CUSTOM_TARGET(install_${COMPONENT_LCASE}
	  ${CMAKE_COMMAND}
	  -D "CMAKE_INSTALL_COMPONENT=${COMPONENT}"
	  -P "${CMAKE_BINARY_DIR}/cmake_install.cmake"
	  )
 endforeach()
 
 set(CPACK_COMPONENTS_ALL ${AMBER_INSTALL_COMPONENTS})