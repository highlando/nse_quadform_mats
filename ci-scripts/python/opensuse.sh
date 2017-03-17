#!/bin/bash

set -e 

zypper -n install cmake libbz2-devel zlib-devel xz-devel arpack-devel suitesparse-devel lapack-devel blas-devel libmatio-devel 
zypper -n install python-devel python3-devel python-numpy python3-numpy python-scipy python3-scipy hdf5-devel python-numpy-devel python3-numpy-devel superlu-devel
mkdir -p build
cd build
cmake ../ -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=ON -DARPACK=ON -DOPENMP=OFF
make
python  python/setup.2.7.py build && python  python/setup.2.7.py install
python3 python/setup.3.4.py build && python3 python/setup.3.4.py install
cd tests/basic
make test
  
