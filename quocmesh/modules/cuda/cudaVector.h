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

#ifndef __CUDAVECTOR_H
#define __CUDAVECTOR_H

namespace cuda {

/**
 * \brief A memory block residing in GPU memory.
 *
 * Automatically handles memory allocation and deallocation on the GPU with CUDA and allows to copy
 * data to and from the host
 *
 * \author Berkels
 */
template <typename DataType>
class GPUMemoryBlock {
protected:
  DataType *_pData;
  const int _size;
public:
  GPUMemoryBlock ( int Size );

  ~GPUMemoryBlock ();

  void setZero ();

  int size ( ) const {
    return _size;
  }

  void copyFromHost ( const DataType * const PHostMem );

  void copyToHost ( DataType *PHostMem ) const;

  DataType *getDataPointer ( ) {
    return _pData;
  }

  const DataType *getDataPointer ( ) const {
    return _pData;
  }
};

/**
 * \brief A vector residing in GPU memory.
 *
 * Endows GPUMemoryBlock with some arithmetic operations (implemented using the CUBLAS library)
 * and allows to copy data to and from the host
 *
 * \warning The CUBLAS library needs to be initialized once using cublasInit() in a cu file or
 *          cuda::initLibCublas ( ) anywhere before any of the arithmetic operations can be
 *          used.
 *
 * \author Berkels
 */
template <typename RealType>
class GPUVector : public GPUMemoryBlock<RealType> {
public:
  typedef RealType DataType;

  GPUVector ( int Size ) : GPUMemoryBlock<RealType> ( Size ) { }

  //! CopyFlag == 1 -> deep copy, CopyFlag == 2 -> only copy size and initialize data with zero.
  GPUVector ( const GPUVector<RealType> &Vec, const int CopyFlag = 1 );

  GPUVector<RealType>& operator= ( const GPUVector<RealType> &Vec );

  GPUVector<RealType>& operator+= ( const GPUVector<RealType> &Vec );

  GPUVector<RealType>& operator-= ( const GPUVector<RealType> &Vec );

  GPUVector<RealType>& operator*= ( const RealType Value );

  RealType operator* ( const GPUVector<RealType> &Vec ) const;

  RealType norm ( ) const;

  RealType normSqr() const {
    return ( ( *this ) * ( *this ) );
  }

  GPUVector<RealType>& addMultiple ( const GPUVector<RealType>& Vec, RealType Factor );
};

/**
 * \brief A vector of vectors residing in GPU memory.
 *
 * \warning Arithmetic operations require CUBLAS to be initialized.
 *
 * \author Berkels
 */
template <typename RealType>
class GPUMultiVector {
protected:
  const int _numComponents;
  GPUVector<RealType> **_pVectors;
public:
  GPUMultiVector ( const int Components, const int Size )
    : _numComponents ( Components ) {
    _pVectors = new GPUVector<RealType>*[Components];
    for ( int i = 0; i < _numComponents; ++i )
      _pVectors[i] = new GPUVector<RealType> ( Size );
  }

  GPUMultiVector ( const int Components, const int *SizeVec )
    : _numComponents ( Components ) {
    _pVectors = new GPUVector<RealType>*[Components];
    for ( int i = 0; i < _numComponents; ++i )
      _pVectors[i] = new GPUVector<RealType> ( SizeVec[i] );
  }

  ~GPUMultiVector ( ) {
    for ( int i = 0; i < _numComponents; ++i )
      delete _pVectors[i];
    delete[] _pVectors;
  }

  GPUVector<RealType> &operator[] ( const int Index ) {
    return *(_pVectors[Index]);
  }

  const GPUVector<RealType> &operator[] ( const int Index ) const {
    return *(_pVectors[Index]);
  }

  GPUMultiVector<RealType>& operator= ( const GPUMultiVector<RealType> &MVec ) {
    for ( int i = 0; i < _numComponents; ++i )
      (*this)[i] = MVec[i];
    return *this;
  }

  GPUMultiVector<RealType>& operator-= ( const GPUMultiVector<RealType> &MVec ) {
    for ( int i = 0; i < _numComponents; ++i )
      (*this)[i] -= MVec[i];
    return *this;
  }

  void setZero () {
    for ( int i = 0; i < _numComponents; ++i )
      (*this)[i].setZero();
  }
};

/**
 * Initialize the CUBLAS library.
 *
 * \author Berkels
 */
void initLibCublas ( );

}

#endif // __CUDAVECTOR_H
