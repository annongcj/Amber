function(add_fortran_library LIBRARY_NAME)

    add_library("${LIBRARY_NAME}" ${ARGN})

    get_target_property(BINARY_DIR "${LIBRARY_NAME}" BINARY_DIR)
    set(LIBRARY_MODULE_DIR "${BINARY_DIR}/modules/${LIBRARY_NAME}")
    set_target_properties(${LIBRARY_NAME} PROPERTIES Fortran_MODULE_DIRECTORY "${LIBRARY_MODULE_DIR}")

    target_include_directories(${LIBRARY_NAME} PUBLIC "${LIBRARY_MODULE_DIR}")
endfunction(add_fortran_library)