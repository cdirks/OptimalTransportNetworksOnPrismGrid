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

#ifndef __CURVATUREOPS_H
#define __CURVATUREOPS_H

#include <boundary.h>

namespace bm {

/**********************************************************************/
//! Abstract integral operator
//! That is constructed from a full operator plus a bidiagonal structure
template <class GlobalOperatorType, class LocalOperatorType, class CyclicBandType>
class IntegralOperatorPlusCyclicBand
      : public GlobalOperatorType {

protected:

  //! Default constructor is forbidden
  IntegralOperatorPlusCyclicBand ();

public:

  //! Copy constructor
  IntegralOperatorPlusCyclicBand ( const IntegralOperatorPlusCyclicBand& op );

  //! Assignment operator
  IntegralOperatorPlusCyclicBand& operator = ( const IntegralOperatorPlusCyclicBand& op );

  //! Destructor
  ~IntegralOperatorPlusCyclicBand ();

  //! Construct from boundaries and local evaluation plus step times bidiagonal
  IntegralOperatorPlusCyclicBand ( const LocalOperatorType& localOp,
                                   const Boundary<typename LocalOperatorType::ParticleType>& boundary,
                                   typename LocalOperatorType::ParticleType::DataType step );

private:

  //! Adds a cyclic band subblock
  void addCyclicBand ( int start, const aol::MultiVector<typename LocalOperatorType::ParticleType::DataType>& vec );
};

/***************************************************/
//! Local diagonal matrix containing angles of pi
template <class ParticleType>
class ConstAngleLocalDiagonal {

public:

  //! Number of diagonals besides the main one
  static const int diags = 0;

  //! Computes a multivector containing the diagonals
  void computeLocally ( const ParticleType& particle, aol::MultiVector<typename ParticleType::DataType>& entries, typename ParticleType::DataType step );
};

/***********************************************/
//! Local bidiagonal matrix containing curvatures
template <class ParticleType>
class CurvatureLocalBidiagonal {

public:

  //! Number of diagonals besides the main one
  static const int diags = 1;

  //! Computes a multivector containing the diagonals
  void computeLocally ( const ParticleType& particle, aol::MultiVector<typename ParticleType::DataType>& entries, typename ParticleType::DataType step );
};

/***********************************************/
//! Local tridiagonal matrix containing curvatures
template <class ParticleType>
class CurvatureLocalTridiagonal {

public:

  //! Number of diagonals besides the main one
  static const int diags = 1;

  //! Computes a multivector containing the diagonals
  void computeLocally ( const ParticleType& particle, aol::MultiVector<typename ParticleType::DataType>& entries, typename ParticleType::DataType step );
};

/**************************************************/
//! Local pentadiagonal matrix containing curvatures
template <class ParticleType>
class CurvatureLocalPentadiagonal {

public:

  //! Number of diagonals besides the main one
  static const int diags = 2;

  //! Computes a multivector containing the diagonals
  void computeLocally ( const ParticleType& particle, aol::MultiVector<typename ParticleType::DataType>& entries, typename ParticleType::DataType step );
};




// Class IntegralOperatorPlusCyclicBand<GlobalOperatorType,LocalOperatorType,CyclicBandType>

template <class GlobalOperatorType, class LocalOperatorType, class CyclicBandType>
IntegralOperatorPlusCyclicBand<GlobalOperatorType, LocalOperatorType, CyclicBandType>::IntegralOperatorPlusCyclicBand ()
  : GlobalOperatorType () {}

template <class GlobalOperatorType, class LocalOperatorType, class CyclicBandType>
IntegralOperatorPlusCyclicBand<GlobalOperatorType, LocalOperatorType, CyclicBandType>::IntegralOperatorPlusCyclicBand
( const IntegralOperatorPlusCyclicBand& op )
    : GlobalOperatorType ( op ) {}

template <class GlobalOperatorType, class LocalOperatorType, class CyclicBandType>
IntegralOperatorPlusCyclicBand<GlobalOperatorType, LocalOperatorType, CyclicBandType>&
IntegralOperatorPlusCyclicBand<GlobalOperatorType, LocalOperatorType, CyclicBandType>::operator =
( const IntegralOperatorPlusCyclicBand& op ) {
  // Beware of self-assignment
  if ( this == &op ) return *this;

  GlobalOperatorType::operator = ( op );

  return *this;
}

template <class GlobalOperatorType, class LocalOperatorType, class CyclicBandType>
IntegralOperatorPlusCyclicBand<GlobalOperatorType, LocalOperatorType, CyclicBandType>::~IntegralOperatorPlusCyclicBand () {}

template <class GlobalOperatorType, class LocalOperatorType, class CyclicBandType>
IntegralOperatorPlusCyclicBand<GlobalOperatorType, LocalOperatorType, CyclicBandType>::IntegralOperatorPlusCyclicBand
( const LocalOperatorType& localOp, const Boundary<typename LocalOperatorType::ParticleType>& boundary,
  typename LocalOperatorType::ParticleType::DataType step )
    : GlobalOperatorType ( localOp, boundary ) {
  CyclicBandType band;

  int pos = 0;

  if ( this->getNumRows () != boundary.getNumberOfSegments ()
       || this->getNumCols () != boundary.getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "IntegralOperatorPlusCyclicBand::IntegralOperatorPlusCyclicBand, matrix size does not equal number of segments", __FILE__, __LINE__ );

  for ( typename std::list<typename LocalOperatorType::ParticleType>::const_iterator it = boundary.begin (); it != boundary.end (); ++it ) {

    int size = it->getNumberOfSegments ();
    aol::MultiVector<typename LocalOperatorType::ParticleType::DataType> temp ( 2 * band.diags + 1, size );
    band.computeLocally ( *it, temp, step );

    addCyclicBand ( pos, temp );
    pos += size;
  }
}

template <class GlobalOperatorType, class LocalOperatorType, class CyclicBandType>
void IntegralOperatorPlusCyclicBand<GlobalOperatorType, LocalOperatorType, CyclicBandType>::addCyclicBand
( int start, const aol::MultiVector<typename LocalOperatorType::ParticleType::DataType>& vec ) {
  int size = vec [0].size ();
  int offset = ( vec.numComponents () - 1 ) / 2;

  for ( int i = 0; i < size; ++i )
    for ( int j = - offset; j < offset + 1; ++j ) {

      int r = start + i;
      int c = start + ( i + j + size ) % size;
      this->set ( r, c, this->get ( r, c ) + vec [j+offset][i] );
    }
}

// Class ConstAngleLocalDiagonal<ParticleType>

template <class ParticleType>
void ConstAngleLocalDiagonal<ParticleType>::computeLocally
( const ParticleType& particle, aol::MultiVector<typename ParticleType::DataType>& entries, typename ParticleType::DataType step ) {
  int size = particle.getNumberOfSegments ();

  if ( entries.numComponents () != 2 * diags + 1 )
    throw aol::DimensionMismatchException ( "ConstAngleLocalDiagonal::computeLocally, vector size does not equal number of segments", __FILE__, __LINE__ );
  for ( int i = 0; i < diags; ++i )
    if ( entries [i].size () != size )
      throw aol::DimensionMismatchException ( "ConstAngleLocalDiagonal::computeLocally, vector size does not equal number of segments", __FILE__, __LINE__ );

  for ( int i = 0; i < size; ++i )
    entries [0][i] = - 0.5 * step;
}

// Class CurvatureLocalBidiagonal<ParticleType>

template <class ParticleType>
void CurvatureLocalBidiagonal<ParticleType>::computeLocally
( const ParticleType& particle, aol::MultiVector<typename ParticleType::DataType>& entries, typename ParticleType::DataType step ) {
  int size = particle.getNumberOfSegments ();

  if ( entries.numComponents () != 2 * diags + 1 )
    throw aol::DimensionMismatchException ( "CurvatureLocalBidiagonal::computeLocally, vector size does not equal number of segments", __FILE__, __LINE__ );

  for ( int i = 0; i < diags; ++i )
    if ( entries [i].size () != size )
      throw aol::DimensionMismatchException ( "CurvatureLocalBidiagonal::computeLocally, vector size does not equal number of segments", __FILE__, __LINE__ );

  aol::Vector<typename ParticleType::DataType> lengths ( size );

  particle.getLengths ( lengths );

  for ( int i = 0; i < size; ++i ) {
    typename ParticleType::DataType len;
    len = lengths.get ( i );
    entries [0][i] = entries [2][i] = 2 * step / ( len * len );
    entries [1][i] = 0;
  }
}

// Class CurvatureLocalTridiagonal<ParticleType>

template <class ParticleType>
void CurvatureLocalTridiagonal<ParticleType>::computeLocally
( const ParticleType& particle, aol::MultiVector<typename ParticleType::DataType>& entries, typename ParticleType::DataType step ) {
  int size = particle.getNumberOfSegments ();

  if ( entries.numComponents () != 2 * diags + 1 )
    throw aol::DimensionMismatchException ( "CurvatureLocalTridiagonal::computeLocally, vector size does not equal number of segments", __FILE__, __LINE__ );
  for ( int i = 0; i < 2 * diags + 1; ++i )
    if ( entries [i].size () != size )
    throw aol::DimensionMismatchException ( "CurvatureLocalTridiagonal::computeLocally, vector size does not equal number of segments", __FILE__, __LINE__ );

  aol::Vector<typename ParticleType::DataType> lengths ( size );
  aol::MultiVector<typename ParticleType::DataType> normals ( 2, size );

  particle.getLengths ( lengths );

  typename ParticleType::ConstSegmentIteratorType it = particle.beginSegment ();
  const typename ParticleType::ConstSegmentIteratorType end = particle.endSegment ();
  typename ParticleType::ConstSegmentIteratorType pit = end; --pit;

  for ( int i = 0; i < size; ++i, pit = it, ++it ) {

    aol::Vec2<typename ParticleType::DataType> n = cornerNormal ( *pit, *it );
    normals [0][i] = n[0];
    normals [1][i] = n[1];
  }

  for ( int i = 0; i < size; ++i, pit = it, ++it ) {

    aol::Vec2<typename ParticleType::DataType> pn, cn, nn;

    for ( int j = 0; j < 2; ++j ) {
      pn [j] = normals [j] [ ( i+size-1 ) %size];
      cn [j] = normals [j] [i];
      nn [j] = normals [j] [ ( i+1 ) %size];
    }

    entries [0][i] = - 2 * step * cn * pn / ( lengths [ ( i+size-1 ) %size] * ( lengths [i] + lengths [ ( i+size-1 ) %size] ) );
    entries [1][i] = 2 * step / ( lengths [i] * lengths [ ( i+size-1 ) %size] );
    entries [2][i] = - 2 * step * cn * nn / ( lengths [i] * ( lengths [i] + lengths [ ( i+size-1 ) %size] ) );
  }

  if ( aol::debugging::matrix ) cerr << "Matrix Entries for implicit curvature" << endl << entries << endl;
}

// Class CurvatureLocalPentadiagonal<ParticleType>

template <class ParticleType>
void neighborhood ( const ParticleType& particle,
                    const typename ParticleType::ConstSegmentIteratorType& center,
                    int before, int after, std::vector<typename ParticleType::ConstSegmentType>& neighborhood ) {
  if ( before < 0 || after < 0 ) throw aol::ParameterException ( "neighborhood, neighborhood size must be positive", __FILE__, __LINE__ );
  neighborhood.resize ( before + after + 1 );

  const typename ParticleType::ConstSegmentIteratorType begin = particle.beginSegment ();
  const typename ParticleType::ConstSegmentIteratorType end = particle.endSegment ();
  typename ParticleType::ConstSegmentIteratorType it = center;

  for ( int i = 0; i < before; ++i, --it ) if ( it == begin ) it = end;

  for ( int i = 0; i < before + after + 1; ++i, ++it ) {
    if ( it == end ) it = begin;
    neighborhood [i] = *it;
  }
}

// Generate implicit coefficents for curvature at a point, given velocities attack at faces, offset o within both vectors
template <class SegmentType, class DataType>
void curvcoeffpt ( std::vector<SegmentType>& neighborhood, aol::Vector<DataType>& coeff, int o ) {
  aol::Vec2<DataType> normal [4];
  for ( int i = 0; i < 4; ++i )
    normal [i] = neighborhood [i+o].getNormal ();

  DataType h1 = neighborhood [1+o].getLength (), h2 = neighborhood [2+o].getLength ();
  aol::Vec2<DataType> factor = ( normal [1] + normal [2] ) / ( h1 * h2 * ( h1 + h2 ) * 8 );

  coeff [3+o] += h1 * normal [3] * factor;
  coeff [2+o] -= h2 * normal [2] * factor;
  coeff [1+o] -= h1 * normal [1] * factor;
  coeff [0+o] += h2 * normal [0] * factor;
}

// Generate implicit coefficents for curvature at a face, given velocities attack at faces
template <class ParticleType>
void curvcoeff ( const ParticleType& particle,
                 const typename ParticleType::ConstSegmentIteratorType& center,
                 aol::Vector<typename ParticleType::DataType>& coeff ) {
  coeff.resize ( 5 ); coeff.setZero ();
  std::vector<typename ParticleType::ConstSegmentType> nb ( 5 );
  neighborhood ( particle, center, 2, 2, nb );
  curvcoeffpt ( nb, coeff, 0 );
  curvcoeffpt ( nb, coeff, 1 ); // More efficient - simply shift?
}

template <class ParticleType>
void CurvatureLocalPentadiagonal<ParticleType>::computeLocally
( const ParticleType& particle, aol::MultiVector<typename ParticleType::DataType>& entries, typename ParticleType::DataType step ) {
  int size = particle.getNumberOfSegments ();

  if ( entries.numComponents () != 2 * diags + 1 )
    throw aol::DimensionMismatchException ( "CurvatureLocalPentadiagonal::computeLocall, vector size does not equal number of segments", __FILE__, __LINE__ );
  for ( int i = 0; i < 2 * diags + 1; ++i )
    if ( entries [i].size () != size )
      throw aol::DimensionMismatchException ( "CurvatureLocalPentadiagonal::computeLocall, vector size does not equal number of segments", __FILE__, __LINE__ );

  aol::Vector<typename ParticleType::DataType> lengths ( size );

  particle.getLengths ( lengths );

  typename ParticleType::ConstSegmentIteratorType it = particle.beginSegment ();
  const typename ParticleType::ConstSegmentIteratorType end = particle.endSegment ();

  for ( int i = 0; i < size; ++i, ++it ) {

    aol::Vector<typename ParticleType::DataType> coeffs ( 5 );
    curvcoeff ( particle, it, coeffs );
    for ( int j = 0; j < 2 * diags + 1; ++j )
      entries [j][i] = - coeffs [j] * step;
  }

  if ( aol::debugging::matrix ) cerr << "Matrix Entries for implicit curvature" << endl << entries << endl;
}

}

#endif
