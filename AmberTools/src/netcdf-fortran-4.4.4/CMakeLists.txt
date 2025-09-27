## This is a CMake file, part of Unidata's netCDF Fortran package.

# It has been heavily modified to slave it to the Amber build system

# Copyright 2012-2014, see the COPYRIGHT file for more information.

##################################
# Set Project Properties
##################################

#Project Name for netCDF-Fortran
PROJECT (NC4F Fortran C)
SET(NC4F_CTEST_PROJECT_NAME "netcdf-fortran")

set(PACKAGE "${NC4F_CTEST_PROJECT_NAME}")# CACHE STRING "")

#Project Version
SET(NC4F_VERSION_MAJOR 4)
SET(NC4F_VERSION_MINOR 4)
SET(NC4F_VERSION_PATCH 4)
SET(NC4F_VERSION_NOTE "")
SET(NC4F_VERSION ${NC4F_VERSION_MAJOR}.${NC4F_VERSION_MINOR}.${NC4F_VERSION_PATCH}${NC4F_VERSION_NOTE})
SET(VERSION ${NC4F_VERSION})
SET(NC4F_LIB_VERSION 6.1.1)
SET(NC4F_SO_VERSION 6)
SET(PACKAGE_VERSION ${VERSION})

# Set known dependencies.
SET(EXTRA_DEPS "")
# Get system configuration,
# Use it to determine osname,
# os release, cpu. These will be used
# when committing to CDash.
find_program(UNAME NAMES uname)
IF(UNAME)
  macro(getuname name flag)
    exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
  endmacro(getuname)
  getuname(osname -s)
  getuname(osrel  -r)
  getuname(cpu    -m)
  set(TMP_BUILDNAME        "${osname}-${osrel}-${cpu}")
ENDIF()

###
# Allow for some customization of the buildname.
# This will make it easier to identify different builds,
# based on values passed from command line/shell scripts.
#
# For ctest scripts, we can use CTEST_BUILD_NAME.
###

SET(BUILDNAME_PREFIX "")#CACHE STRING "")
SET(BUILDNAME_SUFFIX "")#CACHE STRING "")

IF(BUILDNAME_PREFIX)
  SET(TMP_BUILDNAME "${BUILDNAME_PREFIX}-${TMP_BUILDNAME}")
ENDIF()

IF(BUILDNAME_SUFFIX)
  SET(TMP_BUILDNAME "${TMP_BUILDNAME}-${BUILDNAME_SUFFIX}")
ENDIF()

IF(NOT BUILDNAME)
  SET(BUILDNAME "${TMP_BUILDNAME}")# CACHE STRING "Build name variable for CDash")
ENDIF()
###
# End BUILDNAME customization.
###

# For CMAKE_INSTALL_LIBDIR
INCLUDE(GNUInstallDirs)

IF(MSVC)
  SET(GLOBAL PROPERTY USE_FOLDERS ON)
ENDIF()

# auto-configure style checks, other CMake modules.
INCLUDE (TestBigEndian)
INCLUDE (GetPrerequisites)
FIND_PACKAGE(PkgConfig QUIET)

# Enable 'dist and distcheck'.
# File adapted from http://ensc.de/cmake/FindMakeDist.cmake
# TODO: this can't currently be found, need
# FIND_PACKAGE(MakeDist)
# End 'enable dist and distcheck'

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

################################
# End Project Properties
################################


################################
# Set CTest Properties
################################

#ENABLE_TESTING()
#INCLUDE(CTest)

# Copy the CTest customization file into binary directory, as required.
# TODO: need to provide suitable CTestCustom.cmake



# Set variable to define the build type.
INCLUDE(GenerateExportHeader)

################################
# End CTest Properties
################################


################################
# Compiler Configuration
################################


SET (LIB_TYPE STATIC)
IF (SHARED)
  SET(LIB_TYPE SHARED)
ENDIF()


ADD_DEFINITIONS()


################################
# End Compiler Configuration
################################


################################
# Define Utility Macros
################################

# Macro to append files to the EXTRA_DIST files.
SET(EXTRA_DIST "")
MACRO(ADD_EXTRA_DIST files)
  FOREACH(F ${files})
    SET(EXTRA_DIST ${EXTRA_DIST} ${CMAKE_CURRENT_SOURCE_DIR}/${F})
    SET(EXTRA_DIST ${EXTRA_DIST} PARENT_SCOPE)
  ENDFOREACH()
ENDMACRO()


# [LEFT OFF HERE] replacing macros with ones from netcdf-c macros using NC4F_prefix

# Binary tests, but ones which depend on value of 'TEMP_LARGE' being defined.
MACRO(add_bin_env_temp_large_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()


  ADD_TEST(${prefix}_${F} bash "-c" "TEMP_LARGE=${NC4F_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
  ENDIF()
ENDMACRO()


# Tests which are binary, but depend on a particular environmental variable.
MACRO(add_bin_env_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()


  ADD_TEST(${prefix}_${F} bash "-c" "TOPSRCDIR=${NC4F_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
  ENDIF()
ENDMACRO()

# Build a binary used by a script, but don't make a test out of it.
MACRO(build_bin_test F)
  ADD_EXECUTABLE(${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()
ENDMACRO()

# Binary tests which are used by a script looking for a specific name.
MACRO(add_bin_test_no_prefix F)
  build_bin_test(${F})
  ADD_TEST(${F} ${EXECUTABLE_OUTPUT_PATH}/${F})
  IF(MSVC)
    SET_PROPERTY(TEST ${F} PROPERTY FOLDER "tests/")
  ENDIF()
ENDMACRO()

MACRO(add_bin_test prefix F)

  IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${F}.F)
    set(ext ".F")
  ELSE()
    set(ext ".f90")
  ENDIF()

  ADD_EXECUTABLE(${prefix}_${F} ${F}${ext})
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()
  ADD_TEST(${prefix}_${F} ${EXECUTABLE_OUTPUT_PATH}/${prefix}_${F})
  IF(MSVC)
    SET_PROPERTY(TEST ${prefix}_${F} PROPERTY FOLDER "tests/")
  ENDIF()
ENDMACRO()

# A script to print out information at the end of the configuration step.
MACRO(print_conf_summary)
  MESSAGE("")
  MESSAGE("")
  MESSAGE(STATUS "CMake Summary:")
  MESSAGE("")
  MESSAGE(STATUS "\tBuild Type:\t\t\t" ${CMAKE_BUILD_TYPE})
  MESSAGE(STATUS "\tBuilding Shared Libraries:\t" ${BUILD_SHARED_LIBS})
  MESSAGE(STATUS "\tBuilding netCDF-4:\t\t" ${ENABLE_NETCDF4})
  MESSAGE(STATUS "\tBuilding DAP Support:\t\t" ${ENABLE_DAP})
  MESSAGE(STATUS "\tLogging Support:\t\t" ${STATUS_LOGGING})
  IF(STATUS_PNETCDF)
    MESSAGE(STATUS "\tUsing pnetcdf:\t\t\t" ${STATUS_PNETCDF})
  ENDIF(STATUS_PNETCDF)

  MESSAGE(STATUS "\tUsing Parallel IO:\t\t" ${STATUS_PARALLEL})

  IF(CMAKE_PREFIX_PATH)
    MESSAGE(STATUS "\tCMake Prefix Path:\t\t" "${CMAKE_PREFIX_PATH}")
  ENDIF()

  SET(ALL_DEP_LIBS "${NETCDF_C_LIBRARY};${EXTRA_DEPS}")
  MESSAGE(STATUS "\tLinking against:\t\t" "${ALL_DEP_LIBS}")
  FOREACH(_LIB ${EXTRA_DEPS})
    GET_FILENAME_COMPONENT(_LIB_NAME ${_LIB} NAME_WE)
    STRING(REGEX REPLACE "^lib" "" _NAME ${_LIB_NAME})
    LIST(APPEND NCF_LIBS "-l${_NAME}")
  ENDFOREACH()
  MESSAGE("")

ENDMACRO()
# Shell script Macro
MACRO(add_sh_test prefix F)
  IF(NOT MSVC)
    ADD_TEST(${prefix}_${F} bash "-c" "export srcdir=${CMAKE_CURRENT_SOURCE_DIR};export TOPSRCDIR=${NC4F_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}/${F}.sh")
  ENDIF()

ENDMACRO()


# Macro for replacing '/MD' with '/MT'.
# Used only on Windows, /MD tells VS to use the shared
# CRT libs, MT tells VS to use the static CRT libs.
#
# Taken From:
# 	http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_build_my_MSVC_application_with_a_static_runtime.3F
#
MACRO(specify_static_crt_flag)

  SET(vars
    CMAKE_C_FLAGS
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_MINSIZEREL
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELWITHDEBINFO)

  FOREACH(flag_var ${vars})
    IF(${flag_var} MATCHES "/MD")
      STRING(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
    ENDIF()
  ENDFOREACH()

  FOREACH(flag_var ${vars})
    MESSAGE(STATUS " '${flag_var}': ${${flag_var}}")
  ENDFOREACH()
  MESSAGE(STATUS "")

ENDMACRO()

################################
# End Macro Definitions
################################


#####
# Some utility macros/scripts.

# Tests which are binary, but depend on a particular environmental variable.
MACRO(add_bin_env_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()


  ADD_TEST(${prefix}_${F} bash "-c" "TOPSRCDIR=${NC4F_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
  ENDIF()
ENDMACRO()

# Build a binary used by a script, but don't make a test out of it.
MACRO(build_bin_test F ext)
  ADD_EXECUTABLE(${F} ${F}${ext})
  TARGET_LINK_LIBRARIES(${F} ncfortran netcdff ${NETCDF_C_LIBRARY})
  IF(MSVC)
    SET_TARGET_PROPERTIES(${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()
ENDMACRO()

set(LARGE_FILE_TESTS OFF)
set(BUILD_BENCHMARKS OFF)
set(TEST_WITH_VALGRIND OFF)
set(TEST_PARALLEL OFF)
IF(TEST_PARALLEL)
  SET(STATUS_PARALLEL ON)
ENDIF()

# Amber: hardcoded these to match the settings in the netcdf part of the configure2 script
set(USE_NETCDF4 FALSE)
set(USE_NETCDF_V2 TRUE)
set(USE_LOGGING FALSE)
set(BUILD_PARALLEL FALSE)

##
# Set some status messages for the summary.
##

SET(STATUS_LOGGING OFF)
IF(USE_LOGGING)
  SET(STATUS_LOGGING ON)
ENDIF()

SET(STATUS_PARALLEL OFF)
IF(BUILD_PARALLEL)
  SET(STATUS_PARALLEL ON)
ENDIF()

IF(BUILD_DAP)
  SET(ENABLE_DAP ON)
ENDIF()

##
# End status messages/variables.
##

IF(USE_NETCDF4)
  SET(ENABLE_NETCDF4 ON)
  ADD_DEFINITIONS(-DUSE_NETCDF4)
  SET(ENABLE_NETCDF_4 ON CACHE BOOL "")
  SET(ENABLE_NETCDF4 ON CACHE BOOL "")
ENDIF()

IF (UNIX AND ${CMAKE_SIZEOF_VOID_P} MATCHES "8")
  ADD_DEFINITIONS(-DLONGLONG_IS_LONG)
  # for easier debugging of cfortran.h
  
  IF (APPLE)
  	# NOTE: (some versions of) OS X gcc do not support this flag
  	check_c_compiler_flag("-fmacro-backtrace-limit=0" COMPILER_SUPPORTS_MACRO_BT)
    
    if(COMPILER_SUPPORTS_MACRO_BT)
    	SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmacro-backtrace-limit=0")
  	endif()
  endif ()
ENDIF()

# Determine C/Fortran pointer compatibility.
try_compile(HAVE_TS29113_SUPPORT ${CMAKE_CURRENT_BINARY_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/CMakeExtras/test_c_ptrdiff_t.f90")
IF (HAVE_TS29113_SUPPORT)
  message(STATUS "netcdf-fortran: Detected TS29113 support in Fortran compiler")
ELSE ()
  check_type_size(ptrdiff_t SIZEOF_PTRDIFF_T)
	
  if(NOT HAVE_SIZEOF_PTRDIFF_T)
  	message(FATAL_ERROR "Could not determine size of ptrdiff_t")
  endif()
  
  ADD_DEFINITIONS(-DSIZEOF_PTRDIFF_T=${SIZEOF_PTRDIFF_T})
  message(STATUS "netcdf-fortran: No TS29113 support in Fortran compiler")
endif()

IF(USE_NETCDF_V2)
  SET(BUILD_V2 ON)
ENDIF()

set(ENABLE_FORTRAN_TYPE_CHECKS TRUE)

# Toggle whether or not to run tests.
set(ENABLE_TESTS FALSE)

# Set the default fortran builds; default is to build f03
SET(BUILD_F90 "ON")
SET(BUILD_V2 "ON")
SET(BUILD_F03 "ON")

IF(ENABLE_FORTRAN_TYPE_CHECKS)
  include(${CMAKE_CURRENT_SOURCE_DIR}/MatchNetcdfFortranTypes.cmake)

ELSE()
  # Set Fortran types to default.
  SET(NCBYTE_T "byte")
  SET(NCSHORT_T "integer*2")
  SET(NF_INT1_T "byte")
  SET(NF_INT2_T "integer*2")
  # Default is for following to be true.  TODO: test instead
  SET(NF_INT1_IS_C_SIGNED_CHAR "1")
  SET(NF_INT2_IS_C_SHORT "1")
  SET(NF_INT_IS_C_INT "1")
  SET(NF_REAL_IS_C_FLOAT "1")
  SET(NF_DOUBLEPRECISION_IS_C_DOUBLE "1")
ENDIF(ENABLE_FORTRAN_TYPE_CHECKS)

# Create libsrc/nfconfig.inc (don't need nfconfig1.inc temporary file used with autoconf)
configure_file("${NC4F_SOURCE_DIR}/libsrc/nfconfig.in.cmake"
  "${NC4F_BINARY_DIR}/libsrc/nfconfig.inc")
# Create config.h file
configure_file("${NC4F_SOURCE_DIR}/config.h.cmake.in"
  "${NC4F_BINARY_DIR}/config.h")

# For now, just copy a stub file
FILE(COPY "${NC4F_SOURCE_DIR}/nf-config.cmake.in"
  DESTINATION "${NC4F_BINARY_DIR}"
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
FILE(RENAME "${NC4F_BINARY_DIR}/nf-config.cmake.in" "${NC4F_BINARY_DIR}/nf-config")

# make sure previous two files get cleaned up...
SET_DIRECTORY_PROPERTIES (DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${netcdf-fortran_BINARY_DIR}/config.h ${netcdf-fortran_BINARY_DIR}/nf-config)


###
# Allow the user to specify libraries
# to link against, similar to automakes 'LIBS' variable.
###
SET(NC_EXTRA_DEPS "")# CACHE STRING "Additional libraries to link against.")
IF(NC_EXTRA_DEPS)
  STRING(REPLACE " " ";" DEPS_LIST ${NC_EXTRA_DEPS})
  FOREACH(_DEP ${DEPS_LIST})
    STRING(REGEX REPLACE "^-l" "" _LIB ${_DEP})
    FIND_LIBRARY("${_LIB}_DEP" NAMES "${_LIB}" "lib${_LIB}")
    MESSAGE(${${_LIB}_DEP})
    IF(NOT "${_LIB}_DEP")
      MESSAGE(FATAL_ERROR "Error finding ${_LIB}.")
    ELSE()
      MESSAGE(STATUS "Found ${_LIB}: ${${_LIB}_DEP}")
    ENDIF()
    SET(EXTRA_DEPS ${EXTRA_DEPS} "${${_LIB}_DEP}")
  ENDFOREACH()
  MESSAGE("Extra deps: ${EXTRA_DEPS}")
  LIST(REMOVE_DUPLICATES EXTRA_DEPS)
ENDIF()
###
# End user-specified dependent libraries.
###

# Determine whether or not to generate documentation.
set(ENABLE_DOXYGEN OFF)
IF(ENABLE_DOXYGEN)
  IF(DOXYGEN_FOUND)
    # Offer the option to build internal documentation.
    set(ENABLE_INTERNAL_DOCS OFF)
    IF(ENABLE_INTERNAL_DOCS)
      SET(BUILD_INTERNAL_DOCS YES)
    ELSE()
      SET(BUILD_INTERNAL_DOCS NO)
    ENDIF()

    # Option to turn on the TODO list in the doxygen-generated documentation.
    set(ENABLE_DOXYGEN_TASKS OFF)
    IF(ENABLE_DOXYGEN_TASKS)
      SET(SHOW_DOXYGEN_TAG_LIST YES)
    ELSE(ENABLE_DOXYGEN_TASKS)
      SET(SHOW_DOXYGEN_TODO_LIST NO)
    ENDIF(ENABLE_DOXYGEN_TASKS)

    # Specify whether or not 'dot' was found on the system path.
    IF(NC_DOT)
      SET(HAVE_DOT YES)
    ELSE(NC_DOT)
      SET(HAVE_DOT NO)
    ENDIF(NC_DOT)

  ELSE(DOXYGEN_FOUND)
    MESSAGE(STATUS "Unable to build internal documentation.  Doxygen does not appear to be on your executable path. Install doxygen and configure the project again.")
    SET(ENABLE_DOXYGEN OFF)
  ENDIF(DOXYGEN_FOUND)
ENDIF()


INSTALL(PROGRAMS ${NC4F_BINARY_DIR}/nf-config
        DESTINATION "${BINDIR}")

# End autotools-style checks for config.h

# Recurse into other subdirectories.
ADD_SUBDIRECTORY(fortran)

#ADD_SUBDIRECTORY(libsrc)

IF(ENABLE_TESTS)
  ADD_SUBDIRECTORY(nf_test)
  ADD_SUBDIRECTORY(nf03_test)
ENDIF()

# Build the Fortran examples
#ADD_SUBDIRECTORY(examples)
#ADD_SUBDIRECTORY(docs)
#print_conf_summary()
