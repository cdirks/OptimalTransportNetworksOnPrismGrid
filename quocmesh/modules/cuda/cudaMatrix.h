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

#ifndef __CUDAMATRIX_H
#define __CUDAMATRIX_H

#include <cudaVector.h>

namespace cuda {

/**
 * \brief A full matrix residing in GPU memory.
 *
 * Automatically handles memory allocation and deallocation on the GPU by storing the entries in
 * a GPUVector. The entries are stored in column-major format because CUBLAS needs this kind of
 * storage.
 * 
 * \warning The CUBLAS library needs to be initialized once using cublasInit() before this class can be used.
 *
 * \author Berkels
 */
template <typename RealType>
class GPUFullMatrix {
  GPUVector<RealType> _data;
  const int _numRows, _numCols;
public:
  GPUFullMatrix ( const int NumRows, const int NumCols );

  ~GPUFullMatrix () {}

  void setZero ();

  int getNumRows() const {
    return _numRows;
  }

  int getNumCols() const {
    return _numCols;
  }

  void apply ( const GPUVector<RealType> &Arg, GPUVector<RealType> &Dest ) const;

  void applyAdd ( const GPUVector<RealType> &Arg, GPUVector<RealType> &Dest ) const;

  void copyFromHost ( const RealType * const PHostMem, const bool Flip );
};

/**
 * \brief A diagonal matrix residing in GPU memory.
 *
 * \warning The CUBLAS library needs to be initialized once using cublasInit() before this class can be used.
 *
 * \author Berkels
 */
template <typename RealType>
class GPUDiagonalMatrix : protected GPUVector<RealType> {
public:
  GPUDiagonalMatrix ( const int NumRowAndCols )
    : GPUVector<RealType> ( NumRowAndCols ) {}

  int getNumRows() const {
    return this->size();
  }

  int getNumCols() const {
    return this->size();
  }

  void apply ( const GPUVector<RealType> &Arg, GPUVector<RealType> &Dest ) const;

  using GPUVector<RealType>::copyFromHost;
};


/**
 * \brief A sparse matrix class corresponding to 2D quadatric grids residing in GPU memory.
 *
 * \note The implementation is not optimized yet.
 *
 * \author Berkels
 */
template <typename RealType>
class GPU2DQuadraticGridMatrix : protected GPUMemoryBlock<RealType> {
  const int _numRowAndCols;
  const int _gridWidth;
public:
  GPU2DQuadraticGridMatrix ( const int NumRowAndCols )
    : GPUMemoryBlock<RealType> ( 9*NumRowAndCols ),
      _numRowAndCols ( NumRowAndCols ),
      _gridWidth ( static_cast<int> ( sqrt ( static_cast<double> ( _numRowAndCols ) ) ) ) {}

  int getNumRows() const {
    return _numRowAndCols;
  }

  int getNumCols() const {
    return _numRowAndCols;
  }

  void apply ( const GPUVector<RealType> &Arg, GPUVector<RealType> &Dest ) const;

  using GPUMemoryBlock<RealType>::copyFromHost;
};

/**
 * Helper class that encapsulates the creation and destruction of a CUSPARSE handle.
 * In particular, its declaration doesn't need the inclusion of any CUDA header file.
 *
 * \author Berkels
 */
class CUDASparseHandle {
  void *_pHandle;
public:
  CUDASparseHandle ( );
  ~CUDASparseHandle ( );

  const void *getHandlePointer ( ) const {
    return _pHandle;
  }
};

/**
 * \brief A Compressed Sparse Row Format (CSR) matrix residing in GPU memory.
 *
 * \note Requires CUDA 3.2 or later.
 *
 * \author Berkels
 */
template <typename RealType>
class GPUCSRMatrix {
  GPUVector<RealType> _values;
  GPUMemoryBlock<int> _rowStartingIndices;
  GPUMemoryBlock<int> _columnIndices;
  const int _numRows, _numCols;
  void *_pMatDescr;
  const CUDASparseHandle &_sparseHandle;
public:
  GPUCSRMatrix ( const int NumRows,
                 const int NumCols,
                 const int NumEntries,
                 const RealType *Values,
                 const int *RowStartingIndices,
                 const int *ColumnIndices, 
                 const CUDASparseHandle &SparseHandle );

  ~GPUCSRMatrix ();

  void apply ( const GPUVector<RealType> &Arg, GPUVector<RealType> &Dest ) const;
};

}

#endif // __CUDAMATRIX_H
