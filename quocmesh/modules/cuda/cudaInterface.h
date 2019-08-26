/*
 * Copyright (c) 2001-2014 AG Rumpf, INS, Universitaet Bonn                      *
 *                                                                               *
 * The contents of this file are subject to the terms of the Common Development  *
 * and Distribution License Version 1.0 (the "License"); you may not use       *
 * this file except in compliance with the License. You may obtain a copy of     *
 * the License at http://www.opensource.org/licenses/CDDL-1.0                    *
 *                                                                               *
 * Software distributed under the License is distributed on an "AS IS" basis,  *
 * WITHOUT WARRANTY OF ANY KIND, either expressed or implied.                    *
 */

#ifndef __CUDAINTERFACE_H
#define __CUDAINTERFACE_H

// System includes
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// CUDA includes
#include <cuda.h>
#include <cublas.h>

// CUDA 5+ doesn't have cutil.h anymore. Thus, we need to supply supplements for the few things we use from cutil.h.
#if ( CUDA_VERSION >= 5000 )
#include <helper_cuda.h>
#ifndef CUDA_SAFE_CALL
#define CUDA_SAFE_CALL checkCudaErrors
#endif
#ifndef CUT_CHECK_ERROR
#define CUT_CHECK_ERROR getLastCudaError
#endif
#else
#include <cutil.h>
#endif

namespace cuda {

template <typename RealType>
RealType* allocate ( const int Length ) {
  RealType *pointer;
  CUDA_SAFE_CALL ( cudaMalloc(reinterpret_cast<void**>(&pointer), Length*sizeof(RealType)) );
  return pointer;
}

template <typename RealType>
void deallocate ( RealType* Ptr ) {
  CUDA_SAFE_CALL ( cudaFree(Ptr) );
}

template <typename RealType>
void memset ( RealType *Ptr, int Value, const int Length) {
  CUDA_SAFE_CALL ( cudaMemset(reinterpret_cast<void*>(Ptr), Value, Length*sizeof(RealType)) );
}

template <typename RealType>
void memcpyHostToGPU ( RealType *PGPUMem, const RealType * const PHostMem, const int Length ) {
  CUDA_SAFE_CALL ( cudaMemcpy(reinterpret_cast<void*>(PGPUMem), PHostMem, Length*sizeof(RealType),cudaMemcpyHostToDevice) );
}

template <typename RealType>
void memcpyGPUToGPU ( RealType *PGPUMemDest, const RealType * const PGPUMemArg, const int Length ) {
  CUDA_SAFE_CALL ( cudaMemcpy(reinterpret_cast<void*>(PGPUMemDest), reinterpret_cast<const void*>(PGPUMemArg), Length*sizeof(RealType),cudaMemcpyDeviceToDevice) );
}

template <typename RealType>
void memcpyGPUToHost ( RealType *PHostMem, const RealType * const PGPUMem, const int Length ) {
  CUDA_SAFE_CALL ( cudaMemcpy(reinterpret_cast<void*>(PHostMem), PGPUMem, Length*sizeof(RealType),cudaMemcpyDeviceToHost) );
}

void printCurrentDevice ( );

/**
 * Creates a dim3 for the grid size and one for the block size (necessary kernel arguments) given the
 * desired block size and the size of the 2D computational domain.
 *
 * \author Berkels
 */
struct GridBlockConfig2D {
  const dim3 dimBlock;
  const dim3 dimGrid;

  GridBlockConfig2D ( const int blockSizeX, const int blockSizeY, const int NumX, const int NumY )
    : dimBlock ( blockSizeX, blockSizeY ),
      // Make sure the grid covers the full image domain.
      dimGrid ( NumX / blockSizeX + static_cast<int> ( NumX % blockSizeX > 0 ),
                NumY / blockSizeY + static_cast<int> ( NumY % blockSizeY > 0 ) ) {}
};

struct GridBlockConfig3D {
  const dim3 dimBlock;
  const dim3 numBlocks;
  const dim3 dimGrid;

  GridBlockConfig3D ( const int blockSizeX, const int blockSizeY, const int blockSizeZ,
                  const int NumX, const int NumY, const int NumZ )
    : dimBlock ( blockSizeX, blockSizeY, blockSizeZ ),
      // Make sure the grid covers the full image domain.
      numBlocks ( NumX / blockSizeX + static_cast<int> ( NumX % blockSizeX > 0 ),
                  NumY / blockSizeY + static_cast<int> ( NumY % blockSizeY > 0 ),
                  NumZ / blockSizeZ + static_cast<int> ( NumZ % blockSizeZ > 0 )),
      dimGrid ( numBlocks.x, numBlocks.y * numBlocks.z ){}
};

inline __device__ int getGlobalIndex ( int x, int y, int nx ) {
  return ( x + y*nx ); 
}

inline __device__ int getGlobalIndex ( int x, int y, int z, int nx, int ny ) {
  return ( x + (z*ny + y)*nx );
}

/**
 * \author Berkels
 */
template <typename RealType>
void cublasScale ( const int Size, const RealType Value, RealType *PData );

/**
 * \author Berkels
 */
template <typename RealType>
RealType cublasScalarProduct ( const int Size, const RealType *PDataA, const RealType *PDataB );

/**
 * \author Berkels
 */
template <typename RealType>
RealType cublasNormSqr ( const int Size, const RealType *PData );

/**
 * \author Berkels
 */
template <typename RealType>
void cublasAddMultiple ( const int Size, const RealType Factor, const RealType *PDataArg, RealType *PDataDest );

}

#endif // __CUDAINTERFACE_H
