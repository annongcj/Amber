project(gpu_utils)

if(NOT HIP)
    add_subdirectory(check_p2p)
endif()
