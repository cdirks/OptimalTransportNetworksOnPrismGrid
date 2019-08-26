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

#ifndef __BEMOPS_H
#define __BEMOPS_H

#include <boundary.h>

namespace bm {

/****************************/
//! Abstract integral operator
template <class LocalOperatorType>
class IntegralOperator
      : public aol::FullMatrix<typename LocalOperatorType::ParticleType::DataType> {

protected:

  //! Default constructor is forbidden
  IntegralOperator ();

public:

  //! Copy constructor
  IntegralOperator ( const IntegralOperator& op );

  //! Assignment operator
  IntegralOperator<LocalOperatorType>& operator = ( const IntegralOperator& op );

  //! Destructor
  ~IntegralOperator ();

  //! Construct from boundaries and local evaluation
  IntegralOperator ( const LocalOperatorType& localOp, const Boundary<typename LocalOperatorType::ParticleType>& boundary );

  //! Normalize diagonal so that A 1 = 0
  void normalize ();
};

/********************************************************/
//! Abstract integral operator for linear ansatz functions
template <class LocalOperatorType, bool midpointEvaluation = false>
class LinIntegralOperator
      : public IntegralOperator<LocalOperatorType> {

protected:

  //! Default constructor is forbidden
  LinIntegralOperator ();

public:

  //! Copy constructor
  LinIntegralOperator ( const LinIntegralOperator& op );

  //! Assignment operator
  LinIntegralOperator<LocalOperatorType, midpointEvaluation>& operator = ( const LinIntegralOperator& op );

  //! Destructor
  ~LinIntegralOperator ();

  //! Construct from boundaries and local evaluation
  LinIntegralOperator ( const LocalOperatorType& localOp, const Boundary<typename LocalOperatorType::ParticleType>& boundary );
};

/********************************************************/
//! Abstract integral operator for linear ansatz functions
/**
 *  @todo ML,BG: make nicer / review
 */
template <class LocalOperatorType>
class LinJumpIntegralOperator
      : public IntegralOperator<LocalOperatorType> {

protected:

  //! Default constructor is forbidden
  LinJumpIntegralOperator ();

public:

  //! Copy constructor
  LinJumpIntegralOperator ( const LinJumpIntegralOperator& op );

  //! Assignment operator
  LinJumpIntegralOperator<LocalOperatorType>& operator = ( const LinJumpIntegralOperator& op );

  //! Destructor
  ~LinJumpIntegralOperator ();

  //! Construct from boundaries and local evaluation
  LinJumpIntegralOperator ( const LocalOperatorType& localOp, const Boundary<typename LocalOperatorType::ParticleType>& boundary );
};

/****************************/
//! Abstract integral operator
template <class Implementation, class ParticleType>
class LocalOperator {

public:

  //! Types
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  //! Evaluates the kernel
  DataType evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center ) const;

  //! Computes the potential on [-1;1]^2
  void evaluate ( qc::ScalarArray<DataType, qc::QC_2D>& array, const Boundary<ParticleType>& boundary,
                  const aol::Vector<DataType>& bndvalues, DataType a = - 1, DataType b = 1 ) const;

  //! Computes the potential on [-1;1]^2
  void evaluateRot ( qc::ScalarArray<DataType, qc::QC_2D>& array, const Boundary<ParticleType>& boundary,
                     const aol::Vector<DataType>& bndvalues, DataType a = - 1, DataType b = 1, DataType angle = 0 ) const;

  //! Evaluates the potential
  DataType evaluate ( const Boundary<ParticleType>& boundary, aol::Vec2<DataType> point,
                      const aol::Vector<DataType>& bndvalues, double bndpointnumber = -1, double eps = 1e-8  ) const;

  //! Evaluates the potential
  DataType evaluateRot ( const Boundary<ParticleType>& boundary, aol::Vec2<DataType> point,
                         const aol::Vector<DataType>& bndvalues, DataType angle ) const;

  //! Evaluates operator locally for given element pair
  DataType evaluateLocally ( const ConstSegmentType& integrate, const ConstSegmentType& evaluate ) const;

  //! Redefine on subclass as evalution of concrete operator
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool isinterior, bool isend ) const;

private:

  Implementation& implementation ();
  const Implementation& implementation () const;
};

/*********************************************************/
//! Abstract integral operator with linear ansatz functions
template <class Implementation, class ParticleType>
class LinLocalOperator : public LocalOperator<Implementation, ParticleType> {

public:

  //! Types
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  //! Computes the potential on [-1;1]^2
  void evaluate ( qc::ScalarArray<DataType, qc::QC_2D>& array, const Boundary<ParticleType>& boundary,
                  const aol::Vector<DataType>& bndvalues, DataType a = - 1, DataType b = 1 ) const;

  //! Evaluates the potential
  DataType evaluate ( const Boundary<ParticleType>& boundary, aol::Vec2<DataType> point,
                      const aol::Vector<DataType>& bndvalues, double bndpointnumber = -1, double eps = 1e-8 ) const;

  //! Evaluates operator locally for given element pair
  //! bool firstHalf decides which part of a base function to integrate
  DataType evaluateLocally ( const ConstSegmentType& integrate, const ConstSegmentType& evaluate, bool firstHalf ) const;

  //! Redefine on subclass as evalution of concrete operator
  //! bool firstHalf decides which part of a base function to integrate
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool isinterior, bool isend, bool firstHalf ) const;

private:

  Implementation& implementation ();
  const Implementation& implementation () const;
};

/*********************************************************/
//! Abstract integral operator with linear ansatz functions
/**
 *  @todo ML,BG: make nicer / review
 */
template <class Implementation, class ParticleType>
class LinJumpLocalOperator : public LocalOperator<Implementation, ParticleType> {

public:

  //! Types
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  //! Computes the potential on [-1;1]^2
  void evaluate ( qc::ScalarArray<DataType, qc::QC_2D>& array, const Boundary<ParticleType>& boundary,
                  const aol::Vector<DataType>& bndvalues, DataType a = - 1, DataType b = 1 ) const;

  //! Evaluates the potential
  DataType evaluate ( const Boundary<ParticleType>& boundary, aol::Vec2<DataType> point,
                      const aol::Vector<DataType>& bndvalues ) const;

  //! Evaluates operator locally for given element pair
  //! bool firstHalf decides which part of a base function to integrate
  DataType evaluateLocally ( const ConstSegmentType& integrate, const ConstSegmentType& evaluate, bool firstHalf ) const;

  //! Redefine on subclass as evalution of concrete operator
  //! bool firstHalf decides which part of a base function to integrate
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool isinterior, bool isend, bool firstHalf ) const;

private:

  Implementation& implementation ();
  const Implementation& implementation () const;
};


// IntegralOperator<LocalOperatorType>

template <class LocalOperatorType>
IntegralOperator<LocalOperatorType>::IntegralOperator ()
    : aol::FullMatrix<typename LocalOperatorType::ParticleType::DataType> ( 1, 1 ) // Realloc later
{}

template <class LocalOperatorType>
IntegralOperator<LocalOperatorType>::IntegralOperator ( const IntegralOperator& op )
    : aol::FullMatrix<typename LocalOperatorType::ParticleType::DataType> ( op ) {}

template <class LocalOperatorType>
IntegralOperator<LocalOperatorType>& IntegralOperator<LocalOperatorType>::operator = ( const IntegralOperator& op ) {
  // Beware of self-assignment
  if ( this == &op ) return *this;

  aol::FullMatrix<typename LocalOperatorType::ParticleType::DataType>::operator = ( op );

  return *this;
}

template <class LocalOperatorType>
IntegralOperator<LocalOperatorType>::~IntegralOperator () {}

template <class LocalOperatorType>
IntegralOperator<LocalOperatorType>::IntegralOperator ( const LocalOperatorType& localOp,
                                                        const Boundary<typename LocalOperatorType::ParticleType>& boundary )
    : aol::FullMatrix<typename LocalOperatorType::DataType> ( boundary.getNumberOfSegments (), boundary.getNumberOfSegments () ) {
  int i, j;
  ConstAllSegmentIterator<typename LocalOperatorType::ParticleType> iIt = boundary.beginSegment ();
  ConstAllSegmentIterator<typename LocalOperatorType::ParticleType> jIt = boundary.beginSegment ();
  ConstAllSegmentIterator<typename LocalOperatorType::ParticleType> endIt = boundary.endSegment ();

  for ( i = 0, iIt = boundary.beginSegment (); iIt != endIt; ++i, ++iIt )
    for ( j = 0, jIt = boundary.beginSegment (); jIt != endIt; ++j, ++jIt ) {
      typename LocalOperatorType::DataType val = localOp.evaluateLocally ( *jIt, *iIt );
      this->set ( i, j, val );
    }
}

template <class LocalOperatorType>
void IntegralOperator<LocalOperatorType>::normalize ()
{
  aol::Vector<typename LocalOperatorType::DataType> one (this->getNumCols ()), zero (this->getNumRows ());
  one.setAll (1);
  this->apply (one, zero);
  for (int i = 0; i < zero.size (); ++i)
    this->add (i, i, -zero [i]);
}

// Class LinIntegralOperator<LocalOperatorType>

template <class LocalOperatorType, bool midpointEvaluation>
LinIntegralOperator<LocalOperatorType, midpointEvaluation>::LinIntegralOperator ()
  : IntegralOperator<LocalOperatorType> () {}

template <class LocalOperatorType, bool midpointEvaluation>
LinIntegralOperator<LocalOperatorType, midpointEvaluation>::LinIntegralOperator ( const LinIntegralOperator& op )
    : IntegralOperator<LocalOperatorType> ( op ) {}

template <class LocalOperatorType, bool midpointEvaluation>
LinIntegralOperator<LocalOperatorType, midpointEvaluation>& LinIntegralOperator<LocalOperatorType, midpointEvaluation>::operator = ( const LinIntegralOperator& op ) {
  // Beware of self-assignment
  if ( this == &op ) return *this;

  IntegralOperator<LocalOperatorType>::operator = ( op );

  return *this;
}

template <class LocalOperatorType, bool midpointEvaluation>
LinIntegralOperator<LocalOperatorType, midpointEvaluation>::~LinIntegralOperator () {}

template <class LocalOperatorType, bool midpointEvaluation>
LinIntegralOperator<LocalOperatorType, midpointEvaluation>::LinIntegralOperator ( const LocalOperatorType& localOp,
                                                                                  const Boundary<typename LocalOperatorType::ParticleType>& boundary )
    : IntegralOperator<LocalOperatorType> () {
  int ii=0, jj=0, i=0, j=0, is=0, js=0;
  typename std::list<typename LocalOperatorType::ParticleType>::const_iterator piIt = boundary.begin (),
                                                                               pjIt = boundary.begin (),
                                                                               pEnd = boundary.end ();
  typename LocalOperatorType::ParticleType::ConstSegmentIteratorType iIt = piIt->beginSegment (),
                                                                     jIt = pjIt->beginSegment (),
                                                                     iEnd = piIt->endSegment (),
                                                                     jEnd = pjIt->endSegment ();
  int siz = boundary.getNumberOfSegments ();
  this->reallocate ( siz, siz );
  this->setZero ();

  // Store ParticleIterator statuses and Particle sizes to make them accessible in a non sequential manner by OpenMP
  vector< typename std::list<typename LocalOperatorType::ParticleType>::const_iterator > iters;
  vector< int > sizes;
  for ( typename std::list<typename LocalOperatorType::ParticleType>::const_iterator pIt = boundary.begin(); pIt != pEnd; ii+= pIt->size(), ++pIt )
  {
    iters.push_back( pIt );
    sizes.push_back( ii );
  }

#ifdef _OPENMP
#pragma omp parallel for firstprivate(ii,jj,i,j,is,js,piIt,pjIt,pEnd,iIt,jIt,iEnd,jEnd)
#endif
  for ( int pseudo = 0; pseudo < static_cast<int>(boundary.size()); ++pseudo )  {
    piIt = iters[pseudo]; ii = sizes[pseudo];
    for ( i = 0, is = piIt->size (), iIt = piIt->beginSegment (), iEnd = piIt->endSegment (); iIt != iEnd; ++iIt, ++i )  {
      for ( jj = 0, pjIt = boundary.begin (); pjIt != pEnd; ++pjIt, jj += j )  {
        for ( j = 0, js = pjIt->size (), jIt = pjIt->beginSegment (), jEnd = pjIt->endSegment (); jIt != jEnd; ++jIt, ++j ) {
          typename LocalOperatorType::DataType val1 = localOp.evaluateLocally ( *jIt, *iIt, false );
          this->add ( ii + i, jj + j,  val1 );
          #ifndef _OPENMP  // possibly conflicting access, do this only if OpenMP is NOT used
          typename LocalOperatorType::DataType val2 = localOp.evaluateLocally ( *jIt, *iIt, true );
          this->add ( ii + ( midpointEvaluation ? i : ( i + 1 ) % is ), jj + ( j + 1 ) % js, val2 );
          #endif
        }
      }
    }
  }
#ifdef _OPENMP  // if OpenMP IS used do it separately
#pragma omp parallel for firstprivate(ii,jj,i,j,is,js,piIt,pjIt,pEnd,iIt,jIt,iEnd,jEnd)
  for ( int pseudo = 0; pseudo < static_cast<int>(boundary.size()); ++pseudo )  {
    piIt = iters[pseudo]; ii = sizes[pseudo];
    for ( i = 0, is = piIt->size (), iIt = piIt->beginSegment (), iEnd = piIt->endSegment (); iIt != iEnd; ++iIt, ++i )  {
      for ( jj = 0, pjIt = boundary.begin (); pjIt != pEnd; ++pjIt, jj += j )  {
        for ( j = 0, js = pjIt->size (), jIt = pjIt->beginSegment (), jEnd = pjIt->endSegment (); jIt != jEnd; ++jIt, ++j ) {
          typename LocalOperatorType::DataType val2 = localOp.evaluateLocally ( *jIt, *iIt, true );
          this->add ( ii + ( midpointEvaluation ? i : ( i + 1 ) % is ), jj + ( j + 1 ) % js, val2 );
        }
      }
    }
  }
#endif
}

// Class LinJumpIntegralOperator<LocalOperatorType>

template <class LocalOperatorType>
LinJumpIntegralOperator<LocalOperatorType>::LinJumpIntegralOperator ()
  : IntegralOperator<LocalOperatorType> () {}

template <class LocalOperatorType>
LinJumpIntegralOperator<LocalOperatorType>::LinJumpIntegralOperator ( const LinJumpIntegralOperator& op )
    : IntegralOperator<LocalOperatorType> ( op ) {}

template <class LocalOperatorType>
LinJumpIntegralOperator<LocalOperatorType>& LinJumpIntegralOperator<LocalOperatorType>::operator = ( const LinJumpIntegralOperator& op ) {
  // Beware of self-assignment
  if ( this == &op ) return *this;

  IntegralOperator<LocalOperatorType>::operator = ( op );

  return *this;
}

template <class LocalOperatorType>
LinJumpIntegralOperator<LocalOperatorType>::~LinJumpIntegralOperator () {}

template <class LocalOperatorType>
LinJumpIntegralOperator<LocalOperatorType>::LinJumpIntegralOperator ( const LocalOperatorType& localOp,
                                                                      const Boundary<typename LocalOperatorType::ParticleType>& boundary )
    : IntegralOperator<LocalOperatorType> () {
  int ii, jj, i, j, is, js;
  typename std::list<typename LocalOperatorType::ParticleType>::const_iterator piIt = boundary.begin (), pjIt = boundary.begin (),
      pEnd = boundary.end ();
  typename LocalOperatorType::ParticleType::ConstSegmentIteratorType iIt = piIt->beginSegment (), jIt = pjIt->beginSegment ();

  int siz = boundary.getNumberOfSegments (), jumps = boundary.size (); // Additional endpoint per "particle" (is no collocation point, but dof)
  this->reallocate ( siz - jumps, siz );
  this->setZero ();

  for ( ii = 0; piIt != pEnd; ++piIt, ii += i )
    for ( i = 0, is = piIt->size () - 1, iIt = piIt->beginSegment (); i < is; ++iIt, ++i )
      for ( jj = 0, pjIt = boundary.begin (); pjIt != pEnd; ++pjIt, jj += j + 1 )
        for ( j = 0, js = pjIt->size () - 1, jIt = pjIt->beginSegment (); j < js; ++jIt, ++j ) {
          typename LocalOperatorType::DataType val = localOp.evaluateLocally ( *jIt, *iIt, false );
          this->add ( ii + i, jj + j, val );
          val = localOp.evaluateLocally ( *jIt, *iIt, true );
          this->add ( ii + i, jj + j + 1, val ); // !!!!!!!!!!!!!!!!!!!!!
        }
}

// Class LocalOperator<Implementation,ParticleType>

template <class Implementation, class ParticleType>
void LocalOperator<Implementation, ParticleType>::evaluate ( qc::ScalarArray<DataType, qc::QC_2D>& array, const Boundary<ParticleType>& boundary,
                                                             const aol::Vector<DataType>& bndvalues, DataType a, DataType b ) const {
  int sizeX = array.getNumX ();
  int sizeY = array.getNumY ();
  int i, j;
  DataType l = b - a;

  DataType val;

  for ( i = 0; i < sizeX; ++i ) {
    for ( j = 0; j < sizeY; ++j ) {

      DataType x = a + ( l * i ) / (sizeX-1);
      DataType y = a + ( l * j ) / (sizeY-1);

      val = implementation ().evaluate ( boundary, aol::Vec2<DataType> ( x, y ), bndvalues );
      array.set ( i, j, val );
    }
  }
}

template <class Implementation, class ParticleType>
void LocalOperator<Implementation, ParticleType>::evaluateRot ( qc::ScalarArray<DataType, qc::QC_2D>& array, const Boundary<ParticleType>& boundary,
                                                                const aol::Vector<DataType>& bndvalues, DataType a, DataType b, DataType angle ) const {
  int sizeX = array.getNumX ();
  int sizeY = array.getNumY ();
  int i, j;
  DataType l = b - a;

  DataType val;

  for ( i = 0; i < sizeX; ++i ) {
    for ( j = 0; j < sizeY; ++j ) {

      DataType x = a + ( l * i ) / (sizeX-1);
      DataType y = a + ( l * j ) / (sizeY-1);
      aol::Vec2<DataType> p ( x, y );
      rotate ( p, angle );
      val = implementation ().evaluate ( boundary, p, bndvalues );
      array.set ( i, j, val );
    }
  }
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LocalOperator<Implementation, ParticleType>::evaluateRot ( const Boundary<ParticleType>& boundary,
    aol::Vec2<DataType> point,
    const aol::Vector<DataType>& bndvalues,
    DataType angle ) const {
  rotate ( point, angle );
  return evaluate ( boundary, point, bndvalues );
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LocalOperator<Implementation, ParticleType>::evaluate ( const Boundary<ParticleType>& boundary,
    aol::Vec2<DataType> point,
    const aol::Vector<DataType>& bndvalues,
    double bndpointnumber, double eps ) const {
  int i, j, siz = bndvalues.size ();
  DataType result = 0;
  ConstAllSegmentIterator<ParticleType> iIt = boundary.beginSegment ();
  ConstAllSegmentIterator<ParticleType> endIt = boundary.endSegment ();

  typename std::list<ParticleType>::const_iterator partit = boundary.begin (), partend = boundary.end ();
  for ( i = 0; partit != partend; ++partit, i += j ) {
    typename ParticleType::ConstSegmentIteratorType it = partit->beginSegment (), end = partit->endSegment ();
    siz = partit->size ();
    for ( j = 0; it != end; ++it, ++j ) {
      bool isint = false, isstart = false, isend = false;
      if ( bndpointnumber >= i && bndpointnumber < i + siz ) { // bndpoint actually in this particle
        isint = bndpointnumber > i + j && bndpointnumber < i + j + 1;
        isstart = abs (bndpointnumber - (i+j))   < eps || abs (bndpointnumber - siz - (i+j))   < eps;
        isend   = abs (bndpointnumber - (i+j+1)) < eps || abs (bndpointnumber + siz - (i+j+1)) < eps;
      }
      result += implementation ().evaluateLocally ( *it, point, isstart, isstart || isend || isint, isend )
        * bndvalues [i + j];
    }
  }
  return result;
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LocalOperator<Implementation, ParticleType>::evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center ) const {
  return implementation ().evaluate ( point, center );
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LocalOperator<Implementation, ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
    const ConstSegmentType& evaluate ) const {
  aol::Vec2<DataType> point = evaluate.getStart () + 0.5 * evaluate.getDirection (); // Segment center collocation

  return implementation ().evaluateLocally ( integrate, point, false, integrate == evaluate, false );
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LocalOperator<Implementation, ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
    const aol::Vec2<DataType>& evaluate,
    bool isstart, bool isinterior, bool isend ) const {
  return implementation ().evaluateLocally ( integrate, evaluate, isstart, isinterior, isend );
}

template <class Implementation, class ParticleType>
Implementation& LocalOperator<Implementation, ParticleType>::implementation () {
  return static_cast<Implementation&> ( *this );
}

template <class Implementation, class ParticleType>
const Implementation& LocalOperator<Implementation, ParticleType>::implementation () const {
  return static_cast<const Implementation&> ( *this );
}

// Class LinLocalOperator<Implementation,ParticleType>

template <class Implementation, class ParticleType>
void LinLocalOperator<Implementation, ParticleType>::evaluate ( qc::ScalarArray<DataType, qc::QC_2D>& array, const Boundary<ParticleType>& boundary,
                                                                const aol::Vector<DataType>& bndvalues, DataType a, DataType b ) const {
  int sizeX = array.getNumX ();
  int sizeY = array.getNumY ();
  int i, j;

  DataType l = b - a;

  DataType val;

  for ( i = 0; i < sizeX; ++i ) {
    for ( j = 0; j < sizeY; ++j ) {

      DataType x = a + ( l * i ) / (sizeX-1);
      DataType y = a + ( l * j ) / (sizeY-1);
      val = implementation ().evaluate ( boundary, aol::Vec2<DataType> ( x, y ), bndvalues );
      array.set ( i, j, val );
    }
  }
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LinLocalOperator<Implementation, ParticleType>::evaluate ( const Boundary<ParticleType>& boundary,
    aol::Vec2<DataType> point,
    const aol::Vector<DataType>& bndvalues,
    double bndpointnumber, double eps ) const {
  int i, j, siz = bndvalues.size ();
  DataType result = 0;
  ConstAllSegmentIterator<ParticleType> iIt = boundary.beginSegment ();
  ConstAllSegmentIterator<ParticleType> endIt = boundary.endSegment ();

  typename std::list<ParticleType>::const_iterator partit = boundary.begin (), partend = boundary.end ();
  for ( i = 0; partit != partend; ++partit, i += j ) {
    typename ParticleType::ConstSegmentIteratorType it = partit->beginSegment (), end = partit->endSegment ();
    siz = partit->size ();
    for ( j = 0; it != end; ++it, ++j ) {
      bool isint = false, isstart = false, isend = false;
      if ( bndpointnumber >= i && bndpointnumber < i + siz ) { // bndpoint actually in this particle
        isint = bndpointnumber > i + j && bndpointnumber < i + j + 1;
        isstart = abs (bndpointnumber - (i+j))   < eps || abs (bndpointnumber - siz - (i+j))   < eps;
        isend   = abs (bndpointnumber - (i+j+1)) < eps || abs (bndpointnumber + siz - (i+j+1)) < eps;
      }
      result += implementation ().evaluateLocally ( *it, point, isstart, isstart || isend || isint, isend, false )
        * bndvalues [i + j];
    }
  }

  partit = boundary.begin ();
  for ( i = 0; partit != partend; ++partit, i += j ) {
    typename ParticleType::ConstSegmentIteratorType it = partit->beginSegment (), end = partit->endSegment ();
    siz = partit->size ();
    for ( j = 0; it != end; ++it, ++j ) {
      bool isint = false, isstart = false, isend = false;
      if ( bndpointnumber >= i && bndpointnumber < i + siz ) { // bndpoint actually in this particle
        isint = bndpointnumber > i + j && bndpointnumber < i + j + 1;
        isstart = abs (bndpointnumber - (i+j))   < eps || abs (bndpointnumber - siz - (i+j))   < eps;
        isend   = abs (bndpointnumber - (i+j+1)) < eps || abs (bndpointnumber + siz - (i+j+1)) < eps;
      }
      result += implementation ().evaluateLocally ( *it, point, isstart, isstart || isend || isint, isend, true )
                * bndvalues [i + ( j+1 ) % siz];
    }
  }
  return result;
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LinLocalOperator<Implementation, ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
                          const ConstSegmentType& evaluate, bool firstHalf ) const {
  aol::Vec2<DataType> point;

  // Node collocation
  if ( firstHalf ) point = evaluate.getEnd ();
  else point = evaluate.getStart ();

  return implementation ().evaluateLocally ( integrate, point, point == integrate.getStart (),
                                             point == integrate.getStart () || point == integrate.getEnd (), point == integrate.getEnd (), firstHalf );
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LinLocalOperator<Implementation, ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
                          const aol::Vec2<DataType>& evaluate,
                          bool isstart, bool isinterior, bool isend, bool firstHalf ) const {
  return implementation ().evaluateLocally ( integrate, evaluate, isstart, isinterior, isend, firstHalf );
}

template <class Implementation, class ParticleType>
Implementation& LinLocalOperator<Implementation, ParticleType>::implementation () {
  return static_cast<Implementation&> ( *this );
}

template <class Implementation, class ParticleType>
const Implementation& LinLocalOperator<Implementation, ParticleType>::implementation () const {
  return static_cast<const Implementation&> ( *this );
}

// Class LinJumpLocalOperator<Implementation,ParticleType>

template <class Implementation, class ParticleType>
void LinJumpLocalOperator<Implementation, ParticleType>::evaluate ( qc::ScalarArray<DataType, qc::QC_2D>& array, const Boundary<ParticleType>& boundary,
                                                                    const aol::Vector<DataType>& bndvalues, DataType a, DataType b ) const {
  int sizeX = array.getNumX ();
  int sizeY = array.getNumY ();
  int i, j;

  DataType l = b - a;

  DataType val;

  for ( i = 0; i < sizeX; ++i ) {
    for ( j = 0; j < sizeY; ++j ) {

      DataType x = a + ( l * i ) / (sizeX-1);
      DataType y = a + ( l * j ) / (sizeY-1);
      val = implementation ().evaluate ( boundary, aol::Vec2<DataType> ( x, y ), bndvalues );
      array.set ( i, j, val );
    }
  }
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LinJumpLocalOperator<Implementation, ParticleType>::evaluate ( const Boundary<ParticleType>& boundary,
    aol::Vec2<DataType> point,
    const aol::Vector<DataType>& bndvalues ) const {
  int i, j;
  DataType result = 0;

  typename std::list<ParticleType>::const_iterator partit = boundary.begin (), partend = boundary.end ();
  for ( i = 0; partit != partend; ++partit, i += j + 1 ) {
    typename ParticleType::ConstSegmentIteratorType it = partit->beginSegment (), end = partit->endSegment ();
    for ( j = 0; j < partit->size () - 1; ++it, ++j ) {
      result += implementation ().evaluateLocally ( *it, point, false, false, false, false ) * bndvalues [i + j];
    }
  }

  partit = boundary.begin ();
  for ( i = 0; partit != partend; ++partit, i += j + 1 ) {
    typename ParticleType::ConstSegmentIteratorType it = partit->beginSegment (), end = partit->endSegment ();
    for ( j = 0; j < partit->size () - 1; ++it, ++j ) {
      result += implementation ().evaluateLocally ( *it, point, false, false, false, true ) * bndvalues [i + j+1];
    }
  }

  return result;
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LinJumpLocalOperator<Implementation, ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
                              const ConstSegmentType& evaluate, bool firstHalf ) const {
  aol::Vec2<DataType> point = evaluate.getStart ();

  return implementation ().evaluateLocally ( integrate, point, point == integrate.getStart (),
                                             point == integrate.getStart () || point == integrate.getEnd (), point == integrate.getEnd (), firstHalf );
}

template <class Implementation, class ParticleType>
typename ParticleType::DataType LinJumpLocalOperator<Implementation, ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
                              const aol::Vec2<DataType>& evaluate,
                              bool isstart, bool isinterior, bool isend, bool firstHalf ) const {
  return implementation ().evaluateLocally ( integrate, evaluate, isstart, isinterior, isend, firstHalf );
}

template <class Implementation, class ParticleType>
Implementation& LinJumpLocalOperator<Implementation, ParticleType>::implementation () {
  return static_cast<Implementation&> ( *this );
}

template <class Implementation, class ParticleType>
const Implementation& LinJumpLocalOperator<Implementation, ParticleType>::implementation () const {
  return static_cast<const Implementation&> ( *this );
}

}

#endif
