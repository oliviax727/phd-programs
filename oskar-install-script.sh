#!/bin/bash
# RUN VIA SOURCE!!

# shellcheck disable=SC1091,SC2154
source "${software}"/install-scripts/oskar-install.sh
    
source ./venv-install-script.sh
    
cd ..
    
wget https://github.com/OxfordSKA/OSKAR/archive/refs/tags/2.13.0.tar.gz
tar -xzf 2.13.0.tar.gz
cd OSKAR-2.13.0/ || return
    
mkdir -p build
cd build || return
    
cmake -DCMAKE_INSTALL_PREFIX=../../oskareor.data/ -DFIND_OPENCL=ON -DOpenCL_INCLUDE_DIR=/software/setonix/rocm/rocm-6.4.1/include/ -DOpenCL_LIBRARY=/usr/lib64/libOpenCL.so ../
    
make -j8
make install
    
cd ../../phd-programs || return
