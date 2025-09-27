list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if (ENABLE_PMEMD_XRAY_UNIT_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

include(AddFortranLibrary)

add_fortran_library(xray_fft_interface STATIC src/xray_fft_interface.F90)

if ("${PMEMD_XRAY_CPU_FFT_BACKEND}" STREQUAL "FFTW")
    message(FATAL_ERROR "FFTW is not compatible with AMBER license")
    target_compile_definitions(xray_fft_interface PRIVATE -DPMEMD_XRAY_FFT_BACKEND_IS_FFTW)
    target_sources(xray_fft_interface PRIVATE src/xray_fft_interface_impl_fftw.F90)
    target_link_libraries(xray_fft_interface PUBLIC fftw)
elseif("${PMEMD_XRAY_CPU_FFT_BACKEND}" STREQUAL "MKL")
    find_package(MKL REQUIRED)
    target_compile_definitions(xray_fft_interface PRIVATE -DPMEMD_XRAY_FFT_BACKEND_IS_MKL)
    target_sources(xray_fft_interface PRIVATE src/xray_fft_interface_impl_mkl.F90)
    target_link_libraries(xray_fft_interface PUBLIC mkl::cdft)
elseif("${PMEMD_XRAY_CPU_FFT_BACKEND}" STREQUAL "NONE")
    target_compile_definitions(xray_fft_interface PRIVATE -DPMEMD_XRAY_FFT_BACKEND_IS_NONE)
    target_sources(xray_fft_interface PRIVATE src/xray_fft_interface_impl_none.F90)
    message(WARNING "PMEMD_XRAY_CPU_FFT_BACKEND=NONE disables xray functionality of `pmemd` executable")
#elseif("${PMEMD_XRAY_CPU_FFT_BACKEND}" STREQUAL "FFTPACK")
#    target_compile_definitions(xray_fft_interface PRIVATE -DPMEMD_XRAY_FFT_BACKEND_IS_FFTPACK)
#    target_sources(xray_fft_interface PRIVATE src/xray_fft_interface_impl_fftpack.F90)
else()
    message(FATAL_ERROR "PMEMD_XRAY_CPU_FFT_BACKEND must be one of MKL, NONE, FFTW" )
endif()


add_fortran_library(
        pmemd_xray_contracts
    STATIC
    EXCLUDE_FROM_ALL
        src/xray_contracts.F90
)

add_library(PMEMD::xray_contracts ALIAS pmemd_xray_contracts)
target_compile_definitions(
    pmemd_xray_contracts
    PRIVATE
        $<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>,$<CONFIG:MinSizeRel>>:NDEBUG>
)

set(FORTRAN_SOURCES
    src/xray_atomic_scatter_factor.F90
    src/xray_atomic_scatter_factor_impl_cpu.F90
    src/xray_bulk_mask.F90
    src/xray_bulk_mask_data.F90
    src/xray_bulk_mask_impl_cpu.F90
    src/xray_bulk_model.F90
    src/xray_bulk_model_afonine_2013.F90
    src/xray_bulk_model_none.F90
    src/xray_bulk_model_simple.F90
    src/xray_debug_dump.F90
    src/xray_dpartial.F90
    src/xray_dpartial_data.F90
    src/xray_dpartial_impl_cpu.F90
    src/xray_interface2.F90
    src/xray_interface2_data.F90
    src/xray_non_bulk.F90
    src/xray_non_bulk_data.F90
    src/xray_non_bulk_impl_cpu.F90
    src/xray_pure_utils.F90
    src/xray_scaling.F90
    src/xray_scaling_data.F90
    src/xray_scaling_impl_cpu.F90
    src/xray_target.F90
    src/xray_target_least_squares.F90
    src/xray_target_max_likelihood.F90
    src/xray_target_max_likelihood_data.F90
    src/xray_target_max_likelihood_impl_cpu.F90
    src/xray_target_vector_least_squares.F90
    src/xray_target_vector_least_squares_data.F90
    src/xray_target_vector_least_squares_impl_cpu.F90
    src/xray_unit_cell.F90
)

set(FORTRAN_GPU_SOURCES
    src/xray_atomic_scatter_factor_impl_gpu.F90
    src/xray_bulk_mask_impl_gpu.F90
    src/xray_dpartial_impl_gpu.F90
    src/xray_non_bulk_impl_gpu.F90
    src/xray_scaling_impl_gpu.F90
    src/xray_target_vector_least_squares_impl_gpu.F90
)

add_fortran_library(
        pmemd_xray_cpu
    STATIC
    EXCLUDE_FROM_ALL
        "${FORTRAN_SOURCES}"
)

add_library(PMEMD::xray_cpu ALIAS pmemd_xray_cpu)

target_link_libraries(
        pmemd_xray_cpu
    PUBLIC
        xray_fft_interface
        PMEMD::xray_contracts
)

if (CUDA)
    add_subdirectory(cuda)

    add_fortran_library(
            pmemd_xray_gpu
            EXCLUDE_FROM_ALL
        STATIC
            "${FORTRAN_SOURCES}"
            "${FORTRAN_GPU_SOURCES}"
    )

    target_link_libraries(
            pmemd_xray_gpu
        PUBLIC
            xray_fft_interface
            PMEMD::xray_contracts
            PMEMD::xray_cuda
    )

    add_library(PMEMD::xray_gpu ALIAS pmemd_xray_gpu)
endif()

# Coverage

if (ENABLE_PMEMD_XRAY_COVERAGE)
    target_compile_options(
            pmemd_xray_cpu
        PRIVATE
            -fprofile-arcs
            -ftest-coverage
            -pg
    )

    target_link_libraries(
            pmemd_xray_cpu
        PUBLIC
            gcov
    )

    target_link_options(
            pmemd_xray_cpu
        PRIVATE
            --coverage
    )
    if (CUDA)
        target_compile_options(
                pmemd_xray_gpu
            PRIVATE
                -fprofile-arcs
                -ftest-coverage
                -pg
        )

        target_link_libraries(
                pmemd_xray_gpu
            PUBLIC
                gcov
        )

        target_link_options(
                pmemd_xray_gpu
            PRIVATE
                --coverage
        )
    endif()
endif()