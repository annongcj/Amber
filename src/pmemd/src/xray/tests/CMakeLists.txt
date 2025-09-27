find_package(PFUNIT 4.0 REQUIRED CONFIG)
find_package(GTest REQUIRED CONFIG)
include(GoogleTest)

macro(add_pmemd_xray_test arg)
    add_pfunit_ctest (
        "${arg}"
    TEST_SOURCES
        "${arg}.pf"
    LINK_LIBRARIES
        PMEMD::xray_cpu
    )
endmacro()

add_pmemd_xray_test(test_xray_atomic_scatter_factor_impl_cpu)
add_pmemd_xray_test(test_xray_bulk_solvent_mask_based_impl_cpu)
add_pmemd_xray_test(test_xray_debug_dump)
add_pmemd_xray_test(test_xray_interface2)
add_pmemd_xray_test(test_xray_interface2_data)
add_pmemd_xray_test(test_xray_max_likelihood)
add_pmemd_xray_test(test_xray_non_bulk_impl_cpu)
add_pmemd_xray_test(test_xray_pure_utils)
add_pmemd_xray_test(test_xray_scaling_impl_cpu)
add_pmemd_xray_test(test_xray_scaling_impl_cpu_vs_cctbx)
add_pmemd_xray_test(test_xray_unit_cell)


#add_executable(
#    test_xray_cuda
#    test_xray_cuda.cpp
#)
#
#target_link_libraries(
#    test_xray_cuda
#PUBLIC
#    pmemd_xray_cuda
#    GTest::gtest_main
#)
#
#gtest_discover_tests(
#    test_xray_cuda
#)

if(${CMAKE_VERSION} VERSION_LESS "3.14.0")
    message(WARNING "Some test data might be missing. Consider updating CMake version or copy data manually.")
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/data" "${CMAKE_CURRENT_BINARY_DIR}/data")
else()
    # Symbolic link creation requires cmake 3.14
    file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/data" "${CMAKE_CURRENT_BINARY_DIR}/data" SYMBOLIC)
endif()

if (ENABLE_PMEMD_XRAY_COVERAGE)
    include(CodeCoverage)
    setup_target_for_coverage_lcov(
    NAME
        xray_tests_coverage
    EXECUTABLE
        ctest -j 2
    DEPENDENCIES
        test_xray_atomic_scatter_factor_impl_cpu
        test_xray_bulk_solvent_mask_based_impl_cpu
        test_xray_debug_dump
        test_xray_interface2
        test_xray_interface2_data
#        test_xray_cuda
        test_xray_max_likelihood
        test_xray_non_bulk_impl_cpu
        test_xray_pure_utils
        test_xray_scaling_impl_cpu
        test_xray_scaling_impl_cpu_vs_cctbx
        test_xray_unit_cell
    )
endif()


add_executable(
    calc_xray_force_cpu
    calc_xray_force.f90
)

target_link_libraries(
    calc_xray_force_cpu
PUBLIC
    PMEMD::xray_cpu
)

if (CUDA)
    add_executable(
        calc_xray_force_gpu
        calc_xray_force.f90
    )

    target_link_libraries(
        calc_xray_force_gpu
    PUBLIC
        PMEMD::xray_gpu
    )
endif ()