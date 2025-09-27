## This is a CMake file, part of Unidata's netCDF package.
# Copyright 2012-2013, see the COPYRIGHT file for more information.

# Process these files with m4.
SET(m4_SOURCES attr ncx putget t_ncxx)
foreach (f ${m4_SOURCES})
  IF(HAVE_M4)
    GEN_m4(${f})
  ENDIF()
endforeach(f)

# use the freshly m4'ed versions of these sources if possible
if(HAVE_M4)
  set(M4_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR})
else()
  set(M4_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif()

SET(libsrc_SOURCES v1hpg.c nc3dispatch.c
  nc3internal.c var.c dim.c  lookup3.c ncio.c nc_hashmap.c
  ${M4_SRC_DIR}/attr.c
  ${M4_SRC_DIR}/ncx.c
  ${M4_SRC_DIR}/putget.c)

SET(libsrc_SOURCES ${libsrc_SOURCES} pstdint.h ncio.h ncx.h)

IF (BUILD_DISKLESS)
  SET(libsrc_SOURCES ${libsrc_SOURCES} memio.c)

  IF (BUILD_MMAP)
    SET(libsrc_SOURCES ${libsrc_SOURCES} mmapio.c)
  ENDIF( BUILD_MMAP)
ENDIF (BUILD_DISKLESS)

IF (USE_FFIO)
  SET(libsrc_SOURCES ${libsrc_SOURCES} ffio.c)
ELSEIF (USE_STDIO)
    SET(libsrc_SORUCES ${libsrc_SOURCES} ncstdio.c)
ELSE (USE_FFIO)
  SET(libsrc_SOURCES ${libsrc_SOURCES} posixio.c)
ENDIF (USE_FFIO)

add_library(netcdf3 OBJECT ${libsrc_SOURCES})

# The C API man page.

SET(MAN_FILES netcdf.3)

set(ARGS_MANPAGE -DAPI=C)

IF (USE_NETCDF4)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DNETCDF4=TRUE)
ENDIF(USE_NETCDF4)

IF (BUILD_DAP)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DDAP=TRUE)
ENDIF(BUILD_DAP)

IF (BUILD_PARALLEL)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DPARALLEL_IO=TRUE)
ENDIF(BUILD_PARALLEL)

