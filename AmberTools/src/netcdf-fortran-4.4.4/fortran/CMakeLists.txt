SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES("." "${NC4F_BINARY_DIR}/libsrc")

SET(netcdff_SOURCES
  typeSizes.f90 module_netcdf_nc_data.F90 module_netcdf_nc_interfaces.f90
  module_netcdf_nf_data.F90 module_netcdf_nf_interfaces.F90
  module_netcdf_f03.f90
  nf_attio.F90 nf_control.F90 nf_dim.f90 nf_misc.f90
  nf_genatt.f90 nf_geninq.f90 nf_genvar.f90
  nf_vario.F90 nf_var1io.F90 nf_varaio.F90
  nf_varmio.F90 nf_varsio.F90
)

IF (USE_NETCDF4)
  list(APPEND netcdff_SOURCES
    module_netcdf4_nc_interfaces.f90
    module_netcdf4_nf_interfaces.F90 module_netcdf4_f03.f90
    nf_lib.c nf_nc4.f90
    )
ENDIF(USE_NETCDF4)

IF (${USE_LOGGING})
  list(APPEND netcdff_SOURCES nf_logging.F90)
ENDIF(${USE_LOGGING})

# SET(netcdff_DEPENDENCIES
#   netcdf_attributes.f90 netcdf_constants.f90 netcdf_dims.f90
#   netcdf_expanded.f90 netcdf_externals.f90 netcdf_file.f90
#   netcdf_overloads.f90 netcdf_text_variables.f90 netcdf_variables.f90
#   netcdf_visibility.f90 netcdf_eightbyte.f90
#   )

IF (${USE_NETCDF4})
  SET(netcdff_SOURCES ${netcdff_SOURCES}
    netcdf4.f90
    )
  # SET(netcdff_DEPENDENCIES ${netcdff_DEPENDENCIES}
  #   netcdf4_func.f90 netcdf4_externals.f90 netcdf4_visibility.f90
  #   netcdf4.f90 netcdf4_eightbyte.f90
  #   netcdf4_variables.f90 netcdf4_file.f90
  #   )
ELSE ()
  SET(netcdff_SOURCES ${netcdff_SOURCES}
    netcdf.f90 netcdf3_file.f90
    )
ENDIF(${USE_NETCDF4})

if (${BUILD_V2})
  SET(netcdff_SOURCES ${netcdff_SOURCES}
    nf_v2compat.c module_netcdf_fortv2_c_interfaces.f90 nf_fortv2.f90
    )
endif (${BUILD_V2})

if (${BUILD_F03})
  SET(LCHEADERS
    ${NETCDF_FORTRAN_MOD_DIR}/netcdf_nc_data.mod 
    ${NETCDF_FORTRAN_MOD_DIR}/netcdf_nc_interfaces.mod
    ${NETCDF_FORTRAN_MOD_DIR}/netcdf_nf_data.mod
    ${NETCDF_FORTRAN_MOD_DIR}/netcdf_nf_interfaces.mod 
    ${NETCDF_FORTRAN_MOD_DIR}/netcdf_f03.mod
    )
  if (${BUILD_V2})
    SET(LCHEADERS ${LCHEADERS}
      ${NETCDF_FORTRAN_MOD_DIR}/netcdf_fortv2_c_interfaces.mod
      )
  endif (${BUILD_V2})
  if (${USE_NETCDF4})
    SET(LCHEADERS ${LCHEADERS}
      ${NETCDF_FORTRAN_MOD_DIR}/netcdf4_nc_interfaces.mod
      ${NETCDF_FORTRAN_MOD_DIR}/netcdf4_nf_interfaces.mod 
      ${NETCDF_FORTRAN_MOD_DIR}/netcdf4_f03.mod
      )
  endif (${USE_NETCDF4})
endif (${BUILD_F03})

# These are always used
SET(LCCOMMON ${NETCDF_FORTRAN_MOD_DIR}/typesizes.mod ${NETCDF_FORTRAN_MOD_DIR}/netcdf.mod)
SET(COMMONHEADERS ${LCHEADERS} ${LCCOMMON})

# SET(netcdff_DEPENDENCIES ${netcdff_DEPENDENCIES} ${COMMONHEADERS} )
SET(BUILT_SOURCES ${LCHEADERS} ${LCCOMMON})

# A macro for appending file2 to a string
MACRO(APPENDFILE string file2)
  file(READ ${file2} file2_bytes)
  set(${string} "${${string}} 
${file2_bytes}")
ENDMACRO(APPENDFILE file1 file2)

SET(NETCDF_INC ${CMAKE_CURRENT_BINARY_DIR}/netcdf.inc)

# Build netcdf.inc file from selected netcdf2, netcdf3 and netcdf4 apis
set(NETCDF_INC_CONTENTS "!     NetCDF-3.")
APPENDFILE(NETCDF_INC_CONTENTS "netcdf3.inc")

if (${USE_NETCDF4})
  set(NETCDF_INC_CONTENTS "${NETCDF_INC_CONTENTS}
!     NetCDF-4.")
  file(APPEND ${NETCDF_INC} "!     NetCDF-4.")
  APPENDFILE(NETCDF_INC_CONTENTS "netcdf4.inc")
endif (${USE_NETCDF4})

if (${BUILD_V2})
  set(NETCDF_INC_CONTENTS "${NETCDF_INC_CONTENTS}
!     NetCDF-2.")
  APPENDFILE(NETCDF_INC_CONTENTS "netcdf2.inc")
endif (${BUILD_V2})

if (${USE_LOGGING})
  set(NETCDF_INC_CONTENTS "${NETCDF_INC_CONTENTS}
!     This is to turn on netCDF internal logging.
      integer nf_set_log_level
      external nf_set_log_level")
endif (${USE_LOGGING})

file(GENERATE OUTPUT ${NETCDF_INC} CONTENT "${NETCDF_INC_CONTENTS}")

SET(BUILT_SOURCES ${BUILT_SOURCES} netcdf.inc)

include_directories(${NC4F_SOURCE_DIR} ${NC4F_SOURCE_DIR}/libsrc)

# This is what we are building: a convenience library of Fortran 2003 functions.
# This just builds SHARED, even though STATIC is also specified
# ADD_LIBRARY(netcdff STATIC SHARED ${netcdff_SOURCES})
# Builds only static, not shared
ADD_LIBRARY(netcdff ${netcdff_SOURCES})
TARGET_LINK_LIBRARIES(netcdff netcdf)
TARGET_LINK_LIBRARIES(netcdff ${EXTRA_DEPS})
SET_TARGET_PROPERTIES(netcdff PROPERTIES
                      VERSION ${NC4F_LIB_VERSION}
                      SOVERSION ${NC4F_SO_VERSION}
                      )
set_property(TARGET netcdff PROPERTY COMPILE_OPTIONS ${OPT_FFLAGS})
   
# set up module directories (has to be cached so that things in other dirs can see it)
set(NETCDF_FORTRAN_MOD_DIR "${CMAKE_CURRENT_BINARY_DIR}/modules" CACHE INTERNAL "" FORCE)
file(MAKE_DIRECTORY ${NETCDF_FORTRAN_MOD_DIR}) 
                      
config_module_dirs(netcdff ${NETCDF_FORTRAN_MOD_DIR} ${NC4F_BINARY_DIR} ${NC4F_BINARY_DIR}/libsrc ${NC4F_SOURCE_DIR} ${NC4F_SOURCE_DIR}/libsrc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Installation of the program
install_libraries(netcdff)

INSTALL(DIRECTORY ${NETCDF_FORTRAN_MOD_DIR}/
        DESTINATION ${INCDIR}
        FILES_MATCHING PATTERN "*.mod"
        )
        
install(FILES ${NETCDF_INC} DESTINATION ${INCDIR})
