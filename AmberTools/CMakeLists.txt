add_subdirectory(src)
add_subdirectory(test)

if( NOT PMEMD_ONLY )
   add_subdirectory(examples)
endif()
