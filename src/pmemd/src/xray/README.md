# Amber XRAY module

This folder contains Amber xray module. 

### Fortran interface

Public interface of the module is defined in [`src/xray_interface2.F90`](src/xray_interface2.F90), it contains only 6 public functions:

```Fortran
init            ! Initialize the module
finalize        ! Finalize module
calc_force      ! Calculate force on given set of coordinates
get_r_factors   ! Calculate current R-work and R-free factors
get_f_calc      ! Get current F-calc (in original order)
```

See [`src/xray_interface2.F90`](src/xray_interface2.F90) for more details.

### CMake

The module defines two public libraries: 
 - `PMEMD::xray_cpu`, CPU-only version of the interface
 - `PMEMD::xray_gpu`, GPU-only version of the interface (defined only if `CUDA` CMake variable is set)

Consumer target must be linked exactly to one of `PMEMD::xray_*`.

Configuration variables:
  - `ENABLE_PMEMD_XRAY_UNIT_TESTS:BOOL` (defaults to `OFF`).
     Enables pmemd xray tests, requires pFUnit and GoogleTest packages to be pre-installed
  
  - `ENABLE_PMEMD_XRAY_COVERAGE:BOOL` (defaults to `OFF`).
     Adds coverage instrumentation to xray library (works only with gcc and probably clang)

  - `PMEMD_XRAY_CPU_FFT_BACKEND:STRING` (defaults to `NONE`). 
    Defines CPU implementation of FFT used by the module. Possible values:
    - `NONE`, disables xray module for CPU-only versions of pmemd (a runtime error will be thrown)
    - `MKL`, requires MKL package
    - `FFTW`, permanently disabled due to license incompatibility

### Known issues and limitations (January 2022)

#### Optimize CPU-GPU communication

Currently, most of xray submodules has two interchangeable implementations: CPU and GPU. 
Implementations of a module share same interface. 
This enables one to build a version with mixed CPU/GPU modules. 
For example, one can implement new target function/bulk solvent model only on CPU and seamlessly utilize GPU implementations of other modules.
On the other hand, shared interface requires GPU implementation to be in sync with CPU.
This involves redundant CPU/GPU data transfer.

#### Update unit cell parameters 

Unit cell parameters are NOT updated during the simulation. 
This means that xray forces are correctly calculated/applied only in NVT/NVE ensembles. 
Although simulations with xray constraint in NPT ensembles seems to be a strange choice, nothing stops user to get somewhat unexpected results.
I think at least a runtime warning (if not an error) should be raised when xray module gets initialized within NPT ensemble.

#### Sync compiler flags with Amber 
This module is (almost) agnostic to compiler flags/definitions set by Amber CMake toolchain.
Probably needs to be fixed by a couple of `target_compile_definition()`/`target_compile_flags()` calls.

#### Configure precision 
Currently, precision is hardcoded as `double` on CPU and `single` on GPU. 
Probably needs to controlled by CMake flags synced with rest of Amber.
