project(kmmd CXX)

set(KMMD_SOURCE kmmd_forcecalc.C
        hash.c
        pdb_io.C
        scanDB.C
        DB_vs_forcefield.C
        kmmd_context.C
        kmmd_externFortran.C
        kmmd_parseJSON.C)


include_directories(.)
add_compile_options(${OPT_CXXFLAGS})
add_library(kmmd ${KMMD_SOURCE})
install_libraries(kmmd)

