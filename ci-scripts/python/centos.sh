#!/bin/bash

set -e -x
yum clean all
yum makecache fast
yum install -y cmake bzip2-devel zlib-devel xz-devel arpack-devel suitesparse-devel lapack-devel blas-devel matio-devel hdf5-devel SuperLU-devel
yum install -y scipy numpy
mkdir -p build
cd build
cmake ../ -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=ON -DOPENMP=OFF
make
python  python/setup.2.6.py build && python  python/setup.2.6.py install
cd tests/basic
make test

