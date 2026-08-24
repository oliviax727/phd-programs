#!/bin/bash

set -euo pipefail

rm -rf .venv

python3 -m venv .venv

# shellcheck disable=SC1091
source .venv/bin/activate
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m ipykernel install --user --name=.venv

# cmake -DCMAKE_INSTALL_PREFIX=/scratch/mwaeor/ohrw/oskareor.data/ -DFIND_OPENCL=ON -DOpenCL_INCLUDE_DIR=/software/setonix/rocm/rocm-6.4.1/include/ -DOpenCL_LIBRARY=/usr/lib64/libOpenCL.so ../