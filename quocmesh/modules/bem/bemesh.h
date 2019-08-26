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

#ifndef __BEMESH_H
#define __BEMESH_H

#include <smallVec.h>
#include <smallMat.h>
#include <vec.h>
#include <scalarArray.h>

namespace aol {

//! Debugging flags
struct debugging {

  // Status messages
  static const bool particle = true;      //< Particle validity
  static const bool progress = true;      //< Progress marker when computing timestep
  static const bool statistics = true;    //< Statistical data when writing timestep
  static const bool timestep = true;      //< Adapt time discretization
  static const bool arrayverbose = false; //< Verbose mode of ScalarArray
  static const bool refine = true;        //< Adapt geometry

  // Error messages
  static const bool finiteness = true;

  // Detailed debugging messages, used only temporarily
  static const bool lowlevel = false;

  // Show traversal of particle quadtree
  static const bool treetraverse = false;

  // Show discrete data
  static const bool matrix = false;
  static const bool vector = false;

  // Compute solutions in whole domain
  static const bool potential = false;
  static const bool elasticdeform = false;
  static const bool elasticstrain = false;
  static const bool elasticenergy = false;

  // Flow control
  static const bool exitearly = false;
  static const bool waiteach = false;
};

//! Dyadic logarithm, rounded arithmetically
int ld ( int x );

//! Rotation matrix, for counter-clockwise rotation by angle d (given as fraction of 1)
template <class DataType>
aol::Matrix22<DataType> rotationMatrix ( DataType d );

//! Compute Q * mat * Q^T, where Q is a counter-clockwise rotation by angle d (given as fraction of 1)
void rotate ( aol::Matrix22<double>& mat, double d );

//! Rotate vec counter-clockwise by d as a fraction of a full rotation by angle d (given as fraction of 1)
void rotate ( aol::Vec2<double>& vec, double d );

//! Assume that the range [p0, p1) within v contains point coordinates as x0 x1 x2 x3 ... y0 y1 y2 y3 ...,
//! and rotate all points counter-clockwise by angle d (given as fraction of 1)
void rotate ( aol::Vector<double>& v, int p0, int p1, double d );

//! Assume that v consists of two ranges [0, p) and [p, end) that each contain point coordinates as x0 x1 x2 x3 ... y0 y1 y2 y3 ...,
//! and rotate all points counter-clockwise by angle d (given as fraction of 1)
void rotate ( aol::Vector<double>& v, int p, double d );

//! Angle (as fraction of 1) formed by two segments of a line, normal order for left side, reverse order for right side
template <class DataType>
DataType angle ( aol::Vec2<DataType> v, aol::Vec2<DataType> w );

//! Angle (as fraction of 1) from x-axis to v, counter-clockwise. v need not be normalized.
template <class DataType>
DataType angle ( aol::Vec2<DataType> v );

//! Differentiate a scalar function in one direction
template <class DataType>
void differentiate ( const qc::ScalarArray<DataType, qc::QC_2D>& f, qc::ScalarArray<DataType, qc::QC_2D>& df, int dir );

/*******************/
//! Standard iterator
template <class DataType>
class Iterator
      : public iterator<bidirectional_iterator_tag, DataType> {

public:

  //! Default constructor
  Iterator ();

  //! Copy constructor
  Iterator ( const Iterator<DataType>& it );

  //! Assignment Operator
  Iterator<DataType>& operator = ( const Iterator<DataType>& it );

  //! Destructor
  ~Iterator ();
};

/*************************/
//! Standard const iterator
template <class DataType>
class ConstIterator
      : public iterator<input_iterator_tag, DataType> {

public:

  //! Default constructor
  ConstIterator ();

  //! Copy constructor
  ConstIterator ( const ConstIterator<DataType>& it );

  //! Assignment Operator
  ConstIterator<DataType>& operator = ( const ConstIterator<DataType>& it );

  //! Destructor
  ~ConstIterator ();
};



// Implementations

template <class DataType>
aol::Matrix22<DataType> rotationMatrix ( DataType d ) {
  Vec2<DataType> dir ( cos ( 2 * aol::NumberTrait<DataType>::getPi() * d ), sin ( 2 * aol::NumberTrait<DataType>::getPi() * d ) );
  Vec2<DataType> ndir ( dir ); ndir.rotateLeft ();
  return Matrix22<DataType> ( dir [0], ndir [0], dir [1], ndir [1] );
}

template <class DataType>
DataType angle ( aol::Vec2<DataType> v, aol::Vec2<DataType> w ) {
  DataType a = angle ( v ) - angle ( w ) - aol::NumberTrait<long double>::pi;
  while ( a < 0 ) a += 2 * aol::NumberTrait<long double>::pi;
  return a / ( 2*aol::NumberTrait<long double>::pi );
}

template <class DataType>
DataType angle ( aol::Vec2<DataType> v ) {
  DataType a = atan2 ( v [1], v[0] );
  if ( a < 0 ) a += 2 * aol::NumberTrait<long double>::pi;
  return a;
}

template <class DataType>
void differentiate ( const qc::ScalarArray<DataType, qc::QC_2D>& f, qc::ScalarArray<DataType, qc::QC_2D>& df, int dir ) {
  int j, k;
  for ( j = 0; j < f.getNumX (); ++j )
    for ( k = 0; k < f.getNumY (); ++k ) {
      DataType diff = ( dir ? f.dyFD ( j, k ) : f.dxFD ( j, k ) ) * ( dir ? f.getNumY () : f.getNumX () ) / 2;
      df.set ( j, k, diff );
    }
}

// Iterator<DataType>

template <class DataType>
Iterator<DataType>::Iterator ()
    : iterator<bidirectional_iterator_tag, DataType> () {}

template <class DataType>
Iterator<DataType>::Iterator ( const Iterator<DataType>& it )
    : iterator<bidirectional_iterator_tag, DataType> ( it ) {}

template <class DataType>
Iterator<DataType>& Iterator<DataType>::operator = ( const Iterator<DataType>& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  iterator<bidirectional_iterator_tag, DataType>::operator = ( it );

  return *this;
}

template <class DataType>
Iterator<DataType>::~Iterator () {}

// ConstIterator<DataType>

template <class DataType>
ConstIterator<DataType>::ConstIterator ()
    : iterator<bidirectional_iterator_tag, DataType> () {}

template <class DataType>
ConstIterator<DataType>::ConstIterator ( const ConstIterator<DataType>& it )
    : iterator<bidirectional_iterator_tag, DataType> ( it ) {}

template <class DataType>
ConstIterator<DataType>& ConstIterator<DataType>::operator = ( const ConstIterator<DataType>& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  iterator<bidirectional_iterator_tag, DataType>::operator = ( it );

  return *this;
}

template <class DataType>
ConstIterator<DataType>::~ConstIterator () {}

}

#endif
