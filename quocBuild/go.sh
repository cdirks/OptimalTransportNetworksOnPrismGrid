#!/bin/bash
cmake -G"Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DUSE_PNG=1 -DUSE_OPENMP=1 -DUSE_LAPACK=0 -DUSE_BLAS=0 -DUSE_MOSEK=0 -DUSE_C++11=1 -DUSE_TIFF=0 ../quocmesh

