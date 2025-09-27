# the old Makefile build system wrote a Makefile include file called "config.h", which all the Makefiles read
# Many tests and other things also read this file to get their configuration parameters,
# so this script generates a compatibility "config.h" with the variables that are needed for the tests


# if we are installing tests, then we can safely set AMBER_SOURCE=$(AMBERHOME) since all the relevant test files can be read from there.
# Otherwise, we should set it to the source directory, so that test files can be read from there.
if(INSTALL_TESTS)
	set(CONFH_AMBER_SOURCE "$(AMBERHOME)") # This gets evaluated by Make, not CMake 
else()
	set(CONFH_AMBER_SOURCE "${CMAKE_SOURCE_DIR}")
endif()

# --------------------------------------------------------------------
# Convert logical variables

if(BUILD_PYTHON)
	set(CONFH_SKIP_PYTHON "no")
else()
	set(CONFH_SKIP_PYTHON "yes")
endif()

# We can't fully emulate "installtype" since CMake does multiple "installtypes" at once, but we at least try
if(MPI)
	set(CONFH_INSTALLTYPE "parallel")
elseif(OPENMP)
	set(CONFH_INSTALLTYPE "openmp")
else()
	set(CONFH_INSTALLTYPE "serial")
endif()

if(TARGET_OSX)
	# on OS X, anything that links with AMBER libraries, including the tests, must use this flag
	set(CONFH_LDFLAGS "-Wl,-rpath,\$(AMBERHOME)/lib")
else()
	set(CONFH_LDFLAGS "")
endif()

if("${AMBER_TOOLS}" MATCHES "rism")
	set(CONFH_TESTRISMSFF "testrism")
else()
	set(CONFH_TESTRISMSFF "")
endif()

if(GTI)
	set(CONFH_FEP_MODE "gti")
else()
	set(CONFH_FEP_MODE "afe")
endif()

if("${AMBER_TOOLS}" MATCHES "emil")
	set(CONFH_EMIL "EMIL")
else()
	set(CONFH_EMIL "")
endif()

if(BUILD_QUICK)
	set(CONFH_QUICK "yes")
else()
	set(CONFH_QUICK "no")
endif()

if(BUILD_TCPB)
	set(CONFH_TCPB "yes")
else()
	set(CONFH_TCPB "no")
endif()

if(BUILD_REAXFF_PUREMD)
	set(CONFH_REAXFF_PUREMD "yes")
else()
	set(CONFH_REAXFF_PUREMD "no")
endif()

if(AMBERTOOLS_ONLY)
	set(CONFH_HAS_AMBER "")
else()
	set(CONFH_HAS_AMBER "amber")
endif()

if(DOWNLOAD_MINICONDA)
	# use path to miniconda in install directory
	set(CONFH_PYTHON_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_POSTFIX}miniconda/bin/python)
else()
	# use system python executable
	set(CONFH_PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})
endif()

if(mkl_ENABLED)
	set(CONFH_MKL_HOME ${MKL_HOME})
else()
	set(CONFH_MKL_HOME "")
endif()


if(USE_DFTBPLUS)
	set(CONFH_DFTBPLUS "yes")
else()
	set(CONFH_DFTBPLUS "no")
endif()

if(USE_XTB)
	set(CONFH_XTB "yes")
else()
	set(CONFH_XTB "no")
endif()


if(USE_DEEPMDKIT)
	set(CONFH_DEEPMDKIT "yes")
else()
	set(CONFH_DEEPMDKIT "no")
endif()


# --------------------------------------------------------------------
# Calculate FLIBS (the link line used by AMBER to link to netlib)

# figure out which libraries to even put in them
set(FLIBS_LIBRARIES sff pbsa rism fftw3 netlib netcdf)
set(FLIBSF_LIBRARIES netlib xblas)

if(OPENMP)
	list(APPEND FLIBS_LIBRARIES openmp_c)
	list(APPEND FLIBSF_LIBRARIES openmp_fortran)
endif()

if(MPI)
	set(FLIBS_MPI_LIBRARIES sff_mpi pbsa_mpi rism_mpi fftw3_mpi fftw3 netlib netcdf)
	
	if(OPENMP)
		list(APPEND FLIBS_MPI_LIBRARIES openmp_c)
	endif()
endif()

set(ALL_FLIBS_LISTS FLIBS FLIBSF)
if(MPI)
	list(APPEND ALL_FLIBS_LISTS FLIBS_MPI)
endif()

# now, convert target lists to link lines
foreach(LIST ${ALL_FLIBS_LISTS})
	resolve_cmake_library_list(${LIST}_LIBS_PATHS ${${LIST}_LIBRARIES})
	resolved_lib_list_to_link_line(${LIST}_LINK_LINE ${LIST}_LINK_DIRECTORIES ${${LIST}_LIBS_PATHS})
	
	# build up text of variable
	list_to_space_separated(CONFH_${LIST} ${${LIST}_LINK_LINE})
	foreach(LINK_DIR ${${LIST}_LINK_DIRECTORIES})
		set(CONFH_${LIST} "-L${LINK_DIR} ${CONFH_${LIST}}")
	endforeach()
endforeach()

# --------------------------------------------------------------------


configure_file(${CMAKE_CURRENT_LIST_DIR}/config.h.cmake.in ${CMAKE_BINARY_DIR}/config.h)

# The configure script uses symlinks to bridge between these two locations, but symlinks aren't portable and I don't want the hassle of dealing with them
# so we just install the file to both places
install(FILES ${CMAKE_BINARY_DIR}/config.h DESTINATION "${CMAKE_INSTALL_POSTFIX}.")
install(FILES ${CMAKE_BINARY_DIR}/config.h DESTINATION "${CMAKE_INSTALL_POSTFIX}AmberTools/src")

# install a simple Makefile to drive the tests from the install dir
install(FILES ${CMAKE_CURRENT_LIST_DIR}/top-level-test-makefile-pmemd DESTINATION . RENAME Makefile COMPONENT Tests ${TESTS_EXCLUDE_OPTION})
