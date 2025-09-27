# Script to install the Amber license files

# install top-level license file (this is the modified version for the installation layout)
install(FILES ${CMAKE_CURRENT_LIST_DIR}/LICENSE_cmake_version DESTINATION ${CMAKE_INSTALL_POSTFIX}. RENAME LICENSE)

if(NOT PMEMD_ONLY)

# grab AmberTools subsidary license files, and put them in the correct subdirectories so that the top-level license file is still correct
install(FILES ${CMAKE_SOURCE_DIR}/AmberTools/GNU_LGPL_v3 DESTINATION .)
install(FILES ${CMAKE_SOURCE_DIR}/AmberTools/src/netcdf-4.6.1/COPYRIGHT DESTINATION ${LICENSEDIR}/netcdf-4.3.0)
install(FILES ${CMAKE_SOURCE_DIR}/AmberTools/src/reduce/LICENSE.txt DESTINATION ${LICENSEDIR}/reduce)
install(FILES ${CMAKE_SOURCE_DIR}/AmberTools/src/arpack/LICENSE DESTINATION ${LICENSEDIR}/arpack)
install(FILES ${CMAKE_SOURCE_DIR}/AmberTools/src/ucpp-1.3/README DESTINATION ${LICENSEDIR}/ucpp-1.3)
install(FILES ${CMAKE_SOURCE_DIR}/AmberTools/src/cusplibrary-cuda9/LICENSE DESTINATION ${LICENSEDIR}/cusplibrary-cuda9)
install(FILES ${CMAKE_SOURCE_DIR}/AmberTools/src/cifparse/README DESTINATION ${LICENSEDIR}/cifparse)

endif()
