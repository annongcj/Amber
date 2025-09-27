set(AMBERTOOLS_INC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

# To add a new subdirectory, see cmake/WhichTools.cmake and cmake/3rdPartyTools.cmake

#3rd party programs
foreach(DIR ${3RDPARTY_SUBDIRS})
	if(${DIR} STREQUAL libtorch)
		include(${DIR}/libtorch.cmake)
	else()
		add_subdirectory(${DIR})
	endif()
endforeach()

# indicate when we switch from 3rd to 1st party programs
message(STATUS "Reticulating splines...")

#1st party programs
foreach(TOOL ${AMBER_TOOLS})
	if(NOT "${TOOL}" IN_LIST TOOLS_IN_SRC)
		add_subdirectory(${TOOL})
	endif()
endforeach()
