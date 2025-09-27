# PMEMD UNIT TESTS

We use `pFUnit` and `GoogleTest` frameworks to unit test fortran and C++ code. Both are expected to be found by CMake,
it's easier to have them installed system-wide.

### Install prerequisites

1. pFUnit

```bash
cd /tmp
wget https://github.com/Goddard-Fortran-Ecosystem/pFUnit/releases/download/v4.2.1/pFUnit-4.2.1.tar
tar xf ./pFUnit-4.2.1.tar
cd pFUnit
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make tests
sudo make install
``` 


2. Google Test

```bash
cd /tmp
wget https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
unzip release-1.11.0.zip
cd googletest-release-1.11.0
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
sudo make install
```

3. Patched cctbx

To compile [patched](https://github.com/sizmailov/cctbx_project/commit/7c674714681f012568d45066c2d041d3f9d2c3f0) cctbx
run special `bootstrap.py`.

```bash
wget https://raw.githubusercontent.com/sizmailov/cctbx_project/fix-scaling/libtbx/auto_build/bootstrap.py
python bootstrap.py --use-conda
```

### Generate reference data via modified cctbx

```bash
# Use patched cctbx to generate reference data
source path-to-patched-cctbx/build/setpaths.sh

cd path-to-amber-sources/src/pmemd/tests/data
cctbx.python gen_scaling_atomic_test.py
cctbx.python get_ml_alpha_beta_in_zones_test.py
cctbx.python get_ml_alpha_beta_individual_test.py
cctbx.python get_scaling_log_binning_test.py

```

### Run tests and coverage

Now we can build `pmemd` unit tests

```bash
cd path-to-amber-sources
mkdir build && cd build 

# We don't need python for xray modules tests
cmake .. \
    -DCOMPILER=GNU \
    -DDOWNLOAD_MINICONDA=FALSE \
    -DBUILD_PYTHON=FALSE \
    -DCUDA=TRUE \
    -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda/ \
    -DCHECK_UPDATES=FALSE

# Build & run all tests with coverage (`xray_tests_coverage` triggers build of other targets) 
cmake --build ./ \
      --target xray_tests_coverage
      
      
# Alternatively, one can build & run individual tests as usual binaries:  
# Build 
cmake --build ./ \
      --target test_xray_max_likelihood
# Run 
./src/pmemd/tests/test_xray_max_likelihood
```
 