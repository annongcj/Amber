project(pmemd C Fortran)

option(ENABLE_PMEMD_XRAY_UNIT_TESTS "Enable pmemd xray tests" OFF)
option(ENABLE_PMEMD_XRAY_COVERAGE "Enable pmemd xray coverage" OFF)
set(PMEMD_XRAY_CPU_FFT_BACKEND "NONE" CACHE STRING "Backend for FFT module. Must be one of MKL, NONE, FFTW")

add_subdirectory(src)
