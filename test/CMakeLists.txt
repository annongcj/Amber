
# Some unit tests find ndiff and dacdif etc. relative to AMBERHOME,
# so always install these files.
install(PROGRAMS ndiff.awk dacdif numprocs.awk program_error.sh DESTINATION "${CMAKE_INSTALL_POSTFIX}test")

install(DIRECTORY . USE_SOURCE_PERMISSIONS DESTINATION "${CMAKE_INSTALL_POSTFIX}test/" COMPONENT Tests ${TESTS_EXCLUDE_OPTION})
