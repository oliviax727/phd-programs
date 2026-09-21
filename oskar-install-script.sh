#!/bin/bash
# RUN VIA SOURCE!!

salloc --account=mwaeor-gpu --partition=gpu-dev --time=01:00:00 --gres=gpu:4

source $software/install-scripts/oskar-install.sh

cd ..

wget https://github.com/OxfordSKA/OSKAR/archive/refs/tags/2.13.0.tar.gz
tar -xzf 2.13.0.tar.gz
cd OSKAR-2.13.0/

mkdir -p build
cd build

cmake -DCMAKE_INSTALL_PREFIX=/scratch/mwaeor/ohrw/oskareor.data/ -DFIND_OPENCL=ON -DOpenCL_INCLUDE_DIR=/software/setonix/rocm/rocm-6.4.1/include/ -DOpenCL_LIBRARY=/usr/lib64/libOpenCL.so ../

make -j8
make install

cd ../../phd-programs
