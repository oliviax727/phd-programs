OSKAR_REPO_DIR=$1
OSKAR_OUT_DIR=$2

cd $OSKAR_REPO_DIR
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$OSKAR_OUT_DIR -DFIND_OPENCL=ON -DBUILD_INFO=ON -DOpenCL_INCLUDE_DIR=/software/setonix/rocm/5.7.3/include/ -DOpenCL_LIBRARY=/software/setonix/rocm/5.7.3/opencl/lib/libOpenCL.so ../
make -j8
make install