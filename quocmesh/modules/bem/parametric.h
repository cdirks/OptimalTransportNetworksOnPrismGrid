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

#ifndef __PARAMETRIC_H
#define __PARAMETRIC_H

#include <aol.h>

#include <segment.h>
#include <particle.h>
#include <qmException.h>

namespace bm {

// See below
template <class DataType, class IndexType> class ParaParticle;
template <class DataType, class IndexType> class ConstParaSegment;

/********************************************************************/
//! Concrete representation as a side of an axis-aligned parametric
//! No default constructor, since reference member must be initialized
//! Can be moved
template <class _DataType = double, class IndexType = int>
class ParaSegment
      : public ReferenceSegment<ParaSegment<_DataType, IndexType>, ParaParticle<_DataType, IndexType> > {

public:

  // Types
  typedef _DataType DataType;

  //! Default Constructor
  ParaSegment ();

  //! Constructs a side of an parametric
  ParaSegment ( ParaParticle<DataType, IndexType>& parametric, IndexType side );

  //! Copy constructor
  ParaSegment ( const ParaSegment& segment );

  //! Assignment operator
  ParaSegment<DataType, IndexType>& operator = ( const ParaSegment& segment );

  //! Destructor
  ~ParaSegment ();

  //! Move the segment a short step in the diparaion of the normal
  void move ( DataType step );

protected:

  //! Compute the cached values of the points
  void updatePoints ();
};

/********************************************************************/
//! Concrete representation as a side of an axis-aligned parametric
//! No default constructor, since reference member must be initialized
//! Cannot be moved
template <class DataType = double, class IndexType = int>
class ConstParaSegment
      : public ReferenceSegment<ConstParaSegment<DataType, IndexType>, const ParaParticle<DataType, IndexType> > {

public:

  //! Default Constructor
  ConstParaSegment ();

  //! Constructs a side of an parametric
  ConstParaSegment ( const ParaParticle<DataType, IndexType>& parametric, IndexType side );

  //! Copy constructor
  ConstParaSegment ( const ConstParaSegment& segment );

  //! Assignment operator
  ConstParaSegment<DataType, IndexType>& operator = ( const ConstParaSegment& segment );

  //! Destructor
  ~ConstParaSegment ();

protected:

  //! Compute the cached values of the points
  void updatePoints ();
};

// See below
template <class DataType, class IndexType> class ConstParaSegmentIterator;

/*******************************************************************/
//! Enumerates the four sides of the parametric, allows to change them
template <class DataType = double, class IndexType = int>
class ParaSegmentIterator
      : public ReferenceSegmentIterator<ParaParticle<DataType, IndexType>, ParaSegment<DataType, IndexType> > {

public:

  //! Default constructor
  ParaSegmentIterator ();

  //! Constructor for both sentinels
  ParaSegmentIterator ( ParaParticle<DataType, IndexType>& particle, bool end = false );

  //! Copy constructor
  ParaSegmentIterator ( const ParaSegmentIterator& it );

  //! Assignment Operator
  ParaSegmentIterator& operator = ( const ParaSegmentIterator& it );

  //! Destructor
  ~ParaSegmentIterator ();

  //! Dereference operators
  ParaSegment<DataType, IndexType>* operator -> ();
  ParaSegment<DataType, IndexType>& operator * ();
};

/********************************************/
//! Enumerates the four sides of the parametric
template <class DataType = double, class IndexType = int>
class ConstParaSegmentIterator
      : public ReferenceSegmentIterator<const ParaParticle<DataType, IndexType>, ConstParaSegment<DataType, IndexType> > {

public:

  //! Default constructor
  ConstParaSegmentIterator ();

  //! Constructor for both sentinels
  ConstParaSegmentIterator ( const ParaParticle<DataType, IndexType>& particle, bool end = false );

  //! Copy constructor
  ConstParaSegmentIterator ( const ConstParaSegmentIterator& it );

  //! Assignment Operator
  ConstParaSegmentIterator& operator = ( const ConstParaSegmentIterator& it );

  //! Destructor
  ~ConstParaSegmentIterator ();

  //! Dereference operators
  const ConstParaSegment<DataType, IndexType>* operator -> () const;
  const ConstParaSegment<DataType, IndexType>& operator * () const;
};

/**********************/
//! Parametric particle
template <class _DataType = double, class _IndexType = int>
class ParaParticle
      : public Particle < ParaParticle<_DataType, _IndexType>, ParaSegmentIterator<_DataType, _IndexType>, ConstParaSegmentIterator<_DataType, _IndexType>,
      ParaSegment<_DataType, _IndexType>, ConstParaSegment<_DataType, _IndexType>, _DataType > {

public:

  // Configurator
  typedef _IndexType IndexType;
  typedef _DataType DataType;
  typedef ParaSegment<DataType, IndexType> SegmentType;
  typedef ConstParaSegment<DataType, IndexType> ConstSegmentType;
  typedef ParaSegmentIterator<DataType, IndexType> SegmentIteratorType;
  typedef ConstParaSegmentIterator<DataType, IndexType> ConstSegmentIteratorType;

public:

  //! Default constructor
  ParaParticle ();

  //! Creates an ellipsis from center and radii
  ParaParticle ( aol::Vec2<DataType> center, DataType a, DataType b = -1,
                 IndexType numOfPoints = 32, bool clockwise = false, DataType rotation = 0, DataType norm = 2 );

  //! Creates something like a rectangle
  ParaParticle ( aol::Vec2<DataType> center, aol::Vec2<DataType> radii, IndexType numOfPointsEach = 8 );

  //! Creates  a line (not closed!)
ParaParticle ( IndexType numOfSegments, aol::Vec2<DataType> start, aol::Vec2<DataType> end, bool oldStyle = true );

  //! Creates a refined polygon (closed by default)
  ParaParticle ( const std::vector<aol::Vec2<DataType> > & points, IndexType numOfPointsEach = 8, bool startAtEnd = true, bool closed = true, bool oldStyle = true );

  //! Creates a refined polygon
  ParaParticle ( const std::vector<aol::Vec2<DataType> > & points, const aol::Vector<IndexType>& numOfPoints );

  //! Copy constructor
  ParaParticle ( const ParaParticle<DataType, IndexType>& particle );

  //! Assignment operator
  ParaParticle<DataType, IndexType>& operator = ( const ParaParticle& particle );

  //! Destructor
  ~ParaParticle ();

  //! Comparison operators
  bool operator == ( const ParaParticle& particle ) const;
  bool operator != ( const ParaParticle& particle ) const;

  //! Number of elements
  IndexType getNumberOfSegments () const;

  //! Number of Points
  IndexType getNumberOfPoints () const;

  //! Rotate clockwise
  void rotate ( DataType angle );

  //! Move the particle
  void move ( const aol::Vector<DataType>& velocity );

  //! Move the particle
  void move ( const aol::MultiVector<DataType>& displacement, bool facecentered );

  //! Area of the Particle
  DataType getVolume () const;

  //! Surface of the Particle
  DataType getSurface () const;

  //! Average radius of the Particle
  DataType getAverageRadius () const;

  //! Center of the particle
  aol::Vec2<DataType> getCenter () const;

  //! Start SegmentIterator
  ParaSegmentIterator<DataType, IndexType> beginSegment ();

  //! End sentinel
  ParaSegmentIterator<DataType, IndexType> endSegment ();

  //! Start Constant SegmentIterator
  ConstParaSegmentIterator<DataType, IndexType> beginSegment () const;

  //! Constant end sentinel
  ConstParaSegmentIterator<DataType, IndexType> endSegment () const;

  //! First index value
  IndexType beginIndex () const;

  //! End sentinel index value
  IndexType endIndex () const;

  //! Add copy of start-point to end for jump bem operators
  void makeJump ();

  //! Compute bounding box
  void getBoundingBox ( aol::Vec2<DataType>& ll, aol::Vec2<DataType>& ur ) const;

  //! Collision detection, with a given buffer distance
  bool collides ( const ParaParticle& particle, DataType buffer = 0 ) const;

  //! Find if particle exceeds given bounding box
  bool exceeds ( aol::Vec2<DataType> ll, aol::Vec2<DataType> ur ) const;

  //! Find if point is inside the particle
  bool inside ( aol::Vec2<DataType> point ) const;

  //! Return iterator to intersection segment, end if no intersection
  ConstSegmentIteratorType intersects ( const ConstSegmentType& segment ) const;

  //! Unify two particles
  //! Must work for colliding particles
  //! Resulting particle need not be subset of the set theoretical union
  void absorb ( const ParaParticle& particle );

  //! Check by which fraction of the velocity the particle can be moved without removing particles
  DataType getMovementFactor ( const aol::Vector<DataType>& velocity ) const;

  //! Reparametrize particle, angles in degrees, maximum factor between h_min and h_max, additionally removes very small segments
  void reparametrize ( DataType minAngle, DataType maxAngle, DataType hFactor, DataType minLength );

  //! Reparametrize particle
  void reparametrize ( DataType hFactor, DataType minLength );

  //! Remove every n-th point
  void coarsen ( IndexType delEach = 2 );

  //! Insert additional points
  void refine ( IndexType newPointsPerSegment = 8 );

  //! Insert additional points into one segment
  void refineSegment ( IndexType segment, IndexType newPoints, bool morePointsAtEnds = false );

  //! Insert one point between floor(index) and ceil(index), either linearly interpolated or (if on_arc = true) on the circle defined by the three closest points, returns index of new point
  IndexType insertPoint ( double index, bool on_arc = false );

  //! Calculate point index for given arc length, including fractional part (i.e. index = 2.5 means point lies in the middle between 2 and 3)
  double getIndexForArcLength ( double arclength ) const;

  //! Calculate arc length for point at given index, including fractional part (i.e. index = 2.5 means point lies in the middle between 2 and 3)
  double getArcLengthForIndex ( double index ) const;

  //! Insert one point at a given arc length, do not insert a point if distance to next point is less than mindist (measured in arclength). Returns true if a point was inserted.
  bool insertPointArcLength ( double arclength, double mindist = 0 );

  //! Compute angle at a given point, 0 means straight line
  DataType getAngle (IndexType i);

  //! Smooth with one Jacobi step
  void smooth ();

  //! Checks for validity
  bool valid () const;

  //! Output function
  ostream& print ( ostream& os ) const;

  //! Input function
  istream& read ( istream& is );

  //! Number of segments
  IndexType size () const;

  //! Get the coordinates of a point on the particle
  const aol::Vec2<DataType>& getPoint ( IndexType index ) const;

  //! Get a modifialbe reference to a point
  aol::Vec2<DataType>& getPoint ( IndexType index );

private:

  //! List of points, counter-clockwise
  std::vector<aol::Vec2<DataType> > _points;

  // Closed polygon?
  bool _closed;
};

// class ParaSegment<DataType,IndexType>

template <class DataType, class IndexType>
inline ParaSegment<DataType, IndexType>::ParaSegment ()
    : ReferenceSegment<ParaSegment<DataType, IndexType>, ParaParticle<DataType, IndexType> > () {}

template <class DataType, class IndexType>
inline ParaSegment<DataType, IndexType>::ParaSegment ( ParaParticle<DataType, IndexType>& particle, IndexType index )
    : ReferenceSegment<ParaSegment<DataType, IndexType>, ParaParticle<DataType, IndexType> > ( particle, index ) {
  this->updatePoints ();
}

template <class DataType, class IndexType>
inline ParaSegment<DataType, IndexType>::ParaSegment ( const ParaSegment& segment )
    : ReferenceSegment<ParaSegment<DataType, IndexType>, ParaParticle<DataType, IndexType> > ( segment ) {}

template <class DataType, class IndexType>
inline ParaSegment<DataType, IndexType>& ParaSegment<DataType, IndexType>::operator = ( const ParaSegment& segment ) {
  // Beware of self-assignment
  if ( this == &segment ) return *this;

  ReferenceSegment<ParaSegment<DataType, IndexType>, ParaParticle<DataType, IndexType> >::operator = ( segment );

  return *this;
}

template <class DataType, class IndexType>
inline ParaSegment<DataType, IndexType>::~ParaSegment () {}

template <class DataType, class IndexType>
  inline void ParaSegment<DataType, IndexType>::move ( DataType /*step*/ ) {
  throw ( aol::OperationNotPermittedException ( "ParaSegment::move, you should not move a segment of a parametric particle alone since it changes the normals of the neighbors",
                                       __FILE__, __LINE__ ) );
}

template <class DataType, class IndexType>
void ParaSegment<DataType, IndexType>::updatePoints () {
  this->_start = this->_particle->getPoint ( this->_index );
  this->_end = this->_particle->getPoint ( this->_index + 1 );
}

// class ConstParaSegment<DataType,IndexType>

template <class DataType, class IndexType>
inline ConstParaSegment<DataType, IndexType>::ConstParaSegment ()
    : ReferenceSegment<ConstParaSegment<DataType, IndexType>, const ParaParticle<DataType, IndexType> > () {}

template <class DataType, class IndexType>
inline ConstParaSegment<DataType, IndexType>::ConstParaSegment ( const ParaParticle<DataType, IndexType>& particle, IndexType index )
    : ReferenceSegment<ConstParaSegment<DataType, IndexType>, const ParaParticle<DataType, IndexType> > ( particle, index ) {
  this->updatePoints ();
}

template <class DataType, class IndexType>
inline ConstParaSegment<DataType, IndexType>::ConstParaSegment ( const ConstParaSegment& segment )
    : ReferenceSegment<ConstParaSegment<DataType, IndexType>, const ParaParticle<DataType, IndexType> > ( segment ) {}

template <class DataType, class IndexType>
inline ConstParaSegment<DataType, IndexType>& ConstParaSegment<DataType, IndexType>::operator = ( const ConstParaSegment& segment ) {
  // Beware of self-assignment
  if ( this == &segment ) return *this;

  ReferenceSegment<ConstParaSegment<DataType, IndexType>, const ParaParticle<DataType, IndexType> >::operator = ( segment );

  return *this;
}

template <class DataType, class IndexType>
inline ConstParaSegment<DataType, IndexType>::~ConstParaSegment () {}

template <class DataType, class IndexType>
void ConstParaSegment<DataType, IndexType>::updatePoints () {
  this->_start = this->_particle->getPoint ( this->_index );
  this->_end = this->_particle->getPoint ( this->_index + 1 );
}

// Class ParaSegmentIterator<DataType,IndexType>

template <class DataType, class IndexType>
inline ParaSegmentIterator<DataType, IndexType>::ParaSegmentIterator ()
    : ReferenceSegmentIterator<ParaParticle<DataType, IndexType>, ParaSegment<DataType, IndexType> > () {}

template <class DataType, class IndexType>
inline ParaSegmentIterator<DataType, IndexType>::ParaSegmentIterator ( ParaParticle<DataType, IndexType>& particle, bool end )
    : ReferenceSegmentIterator<ParaParticle<DataType, IndexType>, ParaSegment<DataType, IndexType> > ( particle, end ) {}

template <class DataType, class IndexType>
inline ParaSegmentIterator<DataType, IndexType>::ParaSegmentIterator ( const ParaSegmentIterator& it )
    : ReferenceSegmentIterator<ParaParticle<DataType, IndexType>, ParaSegment<DataType, IndexType> > ( it ) {}

template <class DataType, class IndexType>
inline ParaSegmentIterator<DataType, IndexType>& ParaSegmentIterator<DataType, IndexType>::operator = ( const ParaSegmentIterator& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  ReferenceSegmentIterator<ParaParticle<DataType, IndexType>, ParaSegment<DataType, IndexType> >::operator = ( it );

  this->updateSegment ();

  return *this;
}

template <class DataType, class IndexType>
inline ParaSegmentIterator<DataType, IndexType>::~ParaSegmentIterator () {}

template <class DataType, class IndexType>
inline ParaSegment<DataType, IndexType>* ParaSegmentIterator<DataType, IndexType>::operator -> () {
  return &this->_segment;
}

template <class DataType, class IndexType>
inline ParaSegment<DataType, IndexType>& ParaSegmentIterator<DataType, IndexType>::operator * () {
  return this->_segment;
}

// Class ConstParaSegmentIterator<DataType,IndexType>

template <class DataType, class IndexType>
inline ConstParaSegmentIterator<DataType, IndexType>::ConstParaSegmentIterator ()
    : ReferenceSegmentIterator<const ParaParticle<DataType, IndexType>, ConstParaSegment<DataType, IndexType> > () {}

template <class DataType, class IndexType>
inline ConstParaSegmentIterator<DataType, IndexType>::ConstParaSegmentIterator ( const ParaParticle<DataType, IndexType>& particle, bool end )
    : ReferenceSegmentIterator<const ParaParticle<DataType, IndexType>, ConstParaSegment<DataType, IndexType> > ( particle, end ) {}

template <class DataType, class IndexType>
inline ConstParaSegmentIterator<DataType, IndexType>::ConstParaSegmentIterator ( const ConstParaSegmentIterator& it )
    : ReferenceSegmentIterator<const ParaParticle<DataType, IndexType>, ConstParaSegment<DataType, IndexType> > ( it ) {}

template <class DataType, class IndexType>
inline ConstParaSegmentIterator<DataType, IndexType>& ConstParaSegmentIterator<DataType, IndexType>::operator = ( const ConstParaSegmentIterator& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  ReferenceSegmentIterator<const ParaParticle<DataType, IndexType>, ConstParaSegment<DataType, IndexType> >::operator = ( it );

  this->updateSegment ();

  return *this;
}

template <class DataType, class IndexType>
inline ConstParaSegmentIterator<DataType, IndexType>::~ConstParaSegmentIterator () {}

template <class DataType, class IndexType>
inline const ConstParaSegment<DataType, IndexType>* ConstParaSegmentIterator<DataType, IndexType>::operator -> () const {
  return &this->_segment;
}

template <class DataType, class IndexType>
inline const ConstParaSegment<DataType, IndexType>& ConstParaSegmentIterator<DataType, IndexType>::operator * () const {
  return this->_segment;
}

// Class ParaParticle<DataType,IndexType>

template <class DataType, class IndexType>
inline ParaParticle<DataType, IndexType>::ParaParticle ()
    : Particle < ParaParticle<DataType, IndexType>, ParaSegmentIterator<DataType, IndexType>, ConstParaSegmentIterator<DataType>,
  ParaSegment<DataType, IndexType>, ConstParaSegment<DataType, IndexType>, DataType > (), _points (), _closed ( true ) {}

template <class DataType, class IndexType>
inline ParaParticle<DataType, IndexType>::ParaParticle ( aol::Vec2<DataType> center, DataType a, DataType b,
                                                         IndexType numOfPoints, bool clockwise, DataType rotation, DataType norm )
    : Particle < ParaParticle<DataType, IndexType>, ParaSegmentIterator<DataType, IndexType>, ConstParaSegmentIterator<DataType>,
    ParaSegment<DataType, IndexType>, ConstParaSegment<DataType, IndexType>, DataType > (), _points ( numOfPoints ), _closed ( true )  {
  if ( b <= 0 ) b = a;
  if ( a <= 0 )
    throw aol::ParameterException ( "ParaParticle::ParaParticle, parametric particle must have positive radii", __FILE__, __LINE__ );

  DataType phi = rotation * 2 * aol::NumberTrait<long double>::pi;
  aol::Matrix22<DataType> rot ( cos ( phi ), - sin ( phi ), sin ( phi ), cos ( phi ) );
  DataType tau = 2 * aol::NumberTrait<long double>::pi / numOfPoints * ( clockwise ? -1 : 1 );
  DataType t = 0; //- 0.5 * tau;
  for ( IndexType i = 0; i < numOfPoints; ++i, t += tau ) {
    aol::Vec2<DataType> r ( cos ( t ), sin ( t ) );
    r /= r.norm ( norm );
    r [0] *= a; r [1] *= b;
    _points [i] = center + rot * r;
  }
}

template <class DataType, class IndexType>
inline ParaParticle<DataType, IndexType>::ParaParticle ( aol::Vec2<DataType> center, aol::Vec2<DataType> radii, IndexType numOfPointsEach )
    : Particle < ParaParticle<DataType, IndexType>, ParaSegmentIterator<DataType, IndexType>, ConstParaSegmentIterator<DataType>,
    ParaSegment<DataType, IndexType>, ConstParaSegment<DataType, IndexType>, DataType > (), _points ( 4 * numOfPointsEach ), _closed ( true )  {
  // Allow different orientation?
  // if (radii [0] <= 0 || radii [1] <= 0)
  //   throw aol::ParameterException ( "ParaParticle::ParaParticle, parametric particle must have positive radii", __FILE__, __LINE__ );

  aol::Vec2<DataType> p [5] = {aol::Vec2<DataType> ( -radii [0], -radii [1] ),
                               aol::Vec2<DataType> ( radii [0], -radii [1] ),
                               aol::Vec2<DataType> ( radii [0], radii [1] ),
                               aol::Vec2<DataType> ( -radii [0], radii [1] ),
                               aol::Vec2<DataType> ( -radii [0], -radii [1] ) };

  //DataType tau = 1.0 / numOfPointsEach, close = 1 - 0.1 * tau, far = 0.1 * tau;

  for ( IndexType i = 0; i < 4; ++i ) {

    //_points [i * numOfPointsEach + 0] = center + p [i]; // !!!!!!!!!!!!!!!!!!!!!!!
    //_points [i * numOfPointsEach + 1] = center + close * p [i] + far * p [i+1];
    //for (IndexType j = 2; j < numOfPointsEach-1; ++j) {
    //DataType alpha = static_cast<DataType> (j-1) / (numOfPointsEach-2), beta = 1 - alpha;
    for ( IndexType j = 0; j < numOfPointsEach; ++j ) {
      DataType alpha = static_cast<DataType> ( j ) / ( numOfPointsEach ), beta = 1 - alpha;
      _points [i * numOfPointsEach + j] = center + beta * p [i] + alpha * p [i+1];
    }
    //_points [(i+1) * numOfPointsEach - 1] = center + far * p [i] + close * p [i+1];
  }
  return;

  // Smooth rectangle
  /*DataType tau = 2 * aol::NumberTrait<long double>::pi / (numOfPointsEach * 4);
  DataType t = 0.25 * aol::NumberTrait<long double>::pi;
  for (IndexType i = 0; i < numOfPointsEach * 4; ++i, t += tau) {
  DataType aniso = 1.0;
  DataType tt = static_cast<DataType> (i % numOfPointsEach) / numOfPointsEach;
  tt = (tt - 0.5) * (tt - 0.5) * aniso + 1 - 0.25 * aniso;
  _points [i] = center + aol::Vec2<DataType> (tt * radii [0] * cos (t), tt * radii [1] * sin (t));
  }*/
}

template <class DataType, class IndexType>
  inline ParaParticle<DataType, IndexType>::ParaParticle ( IndexType numOfSegments, aol::Vec2<DataType> start, aol::Vec2<DataType> end, bool oldStyle )
    : Particle < ParaParticle<DataType, IndexType>, ParaSegmentIterator<DataType, IndexType>, ConstParaSegmentIterator<DataType>,
                 ParaSegment<DataType, IndexType>, ConstParaSegment<DataType, IndexType>, DataType > (), _points ( numOfSegments + 1 ), _closed ( oldStyle ) {
  for ( IndexType j = 0; j <= numOfSegments; ++j ) {
    DataType alpha = static_cast<DataType> ( j ) / ( numOfSegments ), beta = 1 - alpha;
    _points [j] = beta * start + alpha * end;
  }
}

// TODO: closed ?
template <class DataType, class IndexType>
inline ParaParticle<DataType, IndexType>::ParaParticle ( const std::vector<aol::Vec2<DataType> > & points,
                                                         IndexType numOfPointsEach, bool startAtEnd, bool closed, bool oldStyle )
    : Particle < ParaParticle<DataType, IndexType>, ParaSegmentIterator<DataType, IndexType>, ConstParaSegmentIterator<DataType>,
    ParaSegment<DataType, IndexType>, ConstParaSegment<DataType, IndexType>, DataType > (), _points ( ( points.size () - ( closed ? 0 : 1 ) ) * numOfPointsEach + ( closed ? 0 : 1 ) ), _closed ( oldStyle )  {
  if ( startAtEnd ) {
    typename std::vector<aol::Vec2<DataType> >::const_iterator it = points.begin (), pit = points.end (), end = points.end (); --pit;
    if ( !closed ) --end;
    for ( int i = 0; it != end; pit = it, ++it, ++i ) {

      for ( int j = 0; j < numOfPointsEach; ++j ) {
        double l1 = static_cast<double> ( j ) / numOfPointsEach, l0 = 1 - l1;
        _points [i*numOfPointsEach + j] = l0 * ( *pit ) + l1 * ( *it );
      }
    }
    if ( !closed ) _points [ _points.size()-1 ] = *pit;
  } else {
    typename std::vector<aol::Vec2<DataType> >::const_iterator it = points.begin (), pit = points.begin (), end = points.end (); ++it;
    if ( !closed ) --end;
    for ( int i = 0; pit != end; ++pit, ++it, it = ( it == points.end () ? points.begin () : it ), ++i ) {

      for ( int j = 0; j < numOfPointsEach; ++j ) {
        double l1 = static_cast<double> ( j ) / numOfPointsEach, l0 = 1 - l1;
        _points [i*numOfPointsEach + j] = l0 * ( *pit ) + l1 * ( *it );
      }
    }
    if ( !closed ) _points [ _points.size()-1 ] = *pit;
  }
}

template <class DataType, class IndexType>
inline ParaParticle<DataType, IndexType>::ParaParticle ( const std::vector<aol::Vec2<DataType> > & points,
                                                         const aol::Vector<IndexType>& numOfPoints )
    : Particle < ParaParticle<DataType, IndexType>, ParaSegmentIterator<DataType, IndexType>, ConstParaSegmentIterator<DataType>,
    ParaSegment<DataType, IndexType>, ConstParaSegment<DataType, IndexType>, DataType > (), _points ( numOfPoints.sum () ), _closed ( true )  {
  typename std::vector<aol::Vec2<DataType> >::const_iterator it = points.begin (), pit = points.end (); --pit;
  for ( int i = 0, p = 0; it != points.end (); pit = it, ++it, ++i ) {

    for ( int j = 0; j < numOfPoints [i]; ++j, ++p ) {
      double l1 = static_cast<double> ( j ) / numOfPoints [i], l0 = 1 - l1;
      _points [p] = l0 * ( *pit ) + l1 * ( *it );
    }
  }
}

template <class DataType, class IndexType>
inline ParaParticle<DataType, IndexType>::ParaParticle ( const ParaParticle<DataType, IndexType>& particle )
    : Particle < ParaParticle<DataType, IndexType>, ParaSegmentIterator<DataType, IndexType>, ConstParaSegmentIterator<DataType>,
  ParaSegment<DataType, IndexType>, ConstParaSegment<DataType, IndexType>, DataType > ( particle ), _points ( particle._points ), _closed ( particle._closed ) {}

template <class DataType, class IndexType>
inline ParaParticle<DataType, IndexType>& ParaParticle<DataType, IndexType>::operator = ( const ParaParticle& particle ) {
  // Beware of self-assignment
  if ( this == &particle ) return *this;

  Particle < ParaParticle<DataType, IndexType>, ParaSegmentIterator<DataType, IndexType>, ConstParaSegmentIterator<DataType>,
  ParaSegment<DataType, IndexType>, ConstParaSegment<DataType, IndexType>, DataType >::operator = ( particle );

  _points = particle._points;
  _closed = particle._closed;

  return *this;
}

template <class DataType, class IndexType>
inline ParaParticle<DataType, IndexType>::~ParaParticle () {}

template <class DataType, class IndexType>
inline void ParaParticle<DataType, IndexType>::makeJump () {
  _points.push_back ( _points [0] );
}

template <class DataType, class IndexType>
inline bool ParaParticle<DataType, IndexType>::operator == ( const ParaParticle& particle ) const {
  // Slow
  return _points == particle._points && _closed == particle._closed;
}

template <class DataType, class IndexType>
inline bool ParaParticle<DataType, IndexType>::operator != ( const ParaParticle& particle ) const {
  return ! ( *this == particle );
}

template <class DataType, class IndexType>
inline IndexType ParaParticle<DataType, IndexType>::getNumberOfSegments () const {
  return _points.size () - (_closed ? 0 : 1);
}

template <class DataType, class IndexType>
inline IndexType ParaParticle<DataType, IndexType>::getNumberOfPoints () const {
  return _points.size ();
}

template <class DataType, class IndexType>
void ParaParticle<DataType, IndexType>::rotate ( DataType angle ) {
  aol::Matrix22<DataType> rot ( cos ( angle ), sin ( angle ), - sin ( angle ), cos ( angle ) );

  const ParaParticle p = *this;
  ConstParaSegmentIterator<DataType, IndexType> it = p.beginSegment ();
  ConstParaSegmentIterator<DataType, IndexType> end = p.endSegment ();

  for ( int i = 0; i < static_cast<int> ( _points.size () ); ++i ) {

    aol::Vec2<DataType> p = _points [i];
    rot.mult ( p, _points [i] );
  }
}

template <class DataType, class IndexType>
void ParaParticle<DataType, IndexType>::move ( const aol::Vector<DataType>& velocity ) {
  if ( velocity.size () != size () )
    throw aol::DimensionMismatchException ( "ParaParticle::move, velocity vector length must match number of segments", __FILE__, __LINE__ );

  const ParaParticle p = *this;
  ConstParaSegmentIterator<DataType, IndexType> it = p.beginSegment ();
  ConstParaSegmentIterator<DataType, IndexType> end = p.endSegment ();
  ConstParaSegmentIterator<DataType, IndexType> pit = end; --pit;

  for ( int i = 0; it != end; pit = it, ++it, ++i ) {

    aol::Vec2<DataType> n = cornerNormal ( *pit, *it );
    _points [i] += n * velocity [i];
  }
}

template <class DataType, class IndexType>
void ParaParticle<DataType, IndexType>::move ( const aol::MultiVector<DataType>& displacement, bool facecentered ) {
  int n = size ();
  if ( displacement.numComponents () != 2 || displacement [0].size () != n || displacement [1].size () != n )
    throw aol::DimensionMismatchException ( "ParaParticle::move, velocity vector length must match number of segments", __FILE__, __LINE__ );

  for ( int i = 0; i != n; ++i ) {

    aol::Vec2<DataType> d ( displacement [0][i], displacement [1][i] );
    if ( facecentered ) {
      d += aol::Vec2<DataType> ( displacement [0][ ( i+n-1 ) %n], displacement [1][ ( i+n-1 ) %n] ); d *= 0.5;
    }

    _points [i] += d;
  }
}

template <class DataType, class IndexType>
DataType ParaParticle<DataType, IndexType>::getVolume () const {

  if ( !_closed ) return 0;

  ConstParaSegmentIterator<DataType, IndexType> it = beginSegment ();
  ConstParaSegmentIterator<DataType, IndexType> end = endSegment ();
  DataType volume;

  // With Gauss, approximative
  for ( volume = 0; it != end; ++it ) volume += it->getLength () * it->getNormal () * ( it->getStart () + it->getEnd () );

  // @todo ML: Is it ok to return absolute volume so that inner part of donut is accepted?
  return volume / 4;
}

template <class DataType, class IndexType>
DataType ParaParticle<DataType, IndexType>::getSurface () const {
  ConstParaSegmentIterator<DataType, IndexType> it = beginSegment ();
  ConstParaSegmentIterator<DataType, IndexType> end = endSegment ();
  DataType surface;

  for ( surface = 0; it != end; ++it ) surface += it->getLength ();

  return surface;
}

template <class DataType, class IndexType>
DataType ParaParticle<DataType, IndexType>::getAverageRadius () const {
  ConstParaSegmentIterator<DataType, IndexType> it = beginSegment ();
  ConstParaSegmentIterator<DataType, IndexType> end = endSegment ();
  DataType radius;
  aol::Vec2<DataType> center = getCenter ();

  for ( radius = 0; it != end; ++it ) radius += (it->getStart () - center).norm ();

  return radius / getNumberOfSegments ();
}

template <class DataType, class IndexType>
aol::Vec2<DataType> ParaParticle<DataType, IndexType>::getCenter () const {
  ConstParaSegmentIterator<DataType, IndexType> it = beginSegment ();
  ConstParaSegmentIterator<DataType, IndexType> end = endSegment ();
  aol::Vec2<DataType> center;

  // Crude approximation
  for ( ; it != end; ++it ) center += it->getLength () * ( it->getStart () + it->getEnd () );

  return center / ( getSurface () * 2 );
}

template <class DataType, class IndexType>
void ParaParticle<DataType, IndexType>::getBoundingBox ( aol::Vec2<DataType>& ll, aol::Vec2<DataType>& ur ) const {
  ll = ur = _points [0];

  for ( IndexType i = 1; i < size (); ++i ) {
    aol::Vec2<DataType> pt = _points [i];
    if ( pt [0] < ll [0] ) ll [0] = pt [0];
    if ( pt [1] < ll [1] ) ll [1] = pt [1];
    if ( pt [0] > ur [0] ) ur [0] = pt [0];
    if ( pt [1] > ur [1] ) ur [1] = pt [1];
  }
}

template <class DataType, class IndexType>
bool ParaParticle<DataType, IndexType>::collides ( const ParaParticle& particle, DataType buffer ) const {

  // Grow particle because of buffer
  ParaParticle grownpart = particle;
  aol::Vector<DataType> vel ( grownpart.size () );
  vel.setAll ( buffer );
  grownpart.move ( vel );

  // Some strange collision (e.g. very acute angles intersecting) are not detected
  for ( int i = 0; i < grownpart.size (); ++i )
    if ( inside ( grownpart._points [i] ) ) return true;

  for ( int i = 0; i < size (); ++i )
    if ( grownpart.inside ( _points [i] ) ) return true;

  return false;
}

template <class DataType, class IndexType>
bool ParaParticle<DataType, IndexType>::exceeds ( aol::Vec2<DataType> ll, aol::Vec2<DataType> ur ) const {
  for ( IndexType i = 0; i < size (); ++i ) {
    aol::Vec2<DataType> pt = _points [i];
    if ( pt [0] < ll [0] || pt [1] < ll [1] || pt [0] > ur [0] || pt [1] > ur [1] ) return true;
  }
  return false;
}

template <class DataType, class IndexType>
bool ParaParticle<DataType, IndexType>::inside ( aol::Vec2<DataType> point ) const {
  IndexType count = 0;
  ConstSegmentIteratorType it = beginSegment ();
  ConstSegmentIteratorType end = endSegment ();
  // Count all segments crossing the line y = point [1] left of point [0]
  for ( ; it != end; ++it ) {
    aol::Vec2<DataType> p0 = it->getStart ();
    aol::Vec2<DataType> p1 = it->getEnd ();
    if ( p0 [1] > p1 [1] ) {
      aol::Vec2<DataType> tmp = p0; p0 = p1; p1 = tmp;
    }
    if ( p0 [1] <= point [1] && point [1] < p1 [1] ) { // Attn: Semi-open interval to account for points or segments on the line
      aol::Vec2<DataType> dir = p1 - p0;
      DataType lambda = ( point [1] - p0 [1] ) / dir [1];
      aol::Vec2<DataType> tmp = p0 + lambda * dir;
      if ( tmp [0] < point [0] ) ++count;
    }
  }
  return ( count % 2 ) != 0;
}

template <class DataType, class IndexType>
  ConstParaSegmentIterator<DataType, IndexType> ParaParticle<DataType, IndexType>::intersects ( const ConstSegmentType& segment ) const {
  ConstSegmentIteratorType end = endSegment (), it = end;
  for ( it = beginSegment (); it != end; ++it ) {
    try { it->intersect ( segment ); }
    catch ( aol::ParameterException& e ) { e.consume (); continue; }
    // If no exception, collision was found
    // Exit only on intersection pointing outside
    if ( it->getDirection () * segment.getNormal () > 0 ) break;
  }

  return it; // it == end if no intersection found
}

template <class DataType, class IndexType>
void ParaParticle<DataType, IndexType>::absorb ( const ParaParticle& particle2 ) {

  // Copy oneself, so that we can overwrite
  const ParaParticle particle1 ( *this );
  _points.clear ();

  // Find outside point
  ConstSegmentIteratorType begin1 = particle1.beginSegment (), end1 = particle1.endSegment (),
    begin2 = particle2.beginSegment (), end2 = particle2.endSegment (), it1 = begin1;
  for ( it1 = begin1; it1 != end1; ++it1 )
    if ( ! particle2.inside ( it1->getStart () ) )
      break;

  // Continue until we find this point again
  ConstSegmentIteratorType newbegin = it1;
  do {
    _points.push_back ( it1->getStart () );
    ConstSegmentIteratorType it2 = particle2.intersects ( *it1 );
    // It intersection, continue in particle2 until next intersection
    if ( it2 != end2 ) {
      _points.push_back ( it1->intersect ( *it2 ) );
      do {
  ++it2; if ( it2 == end2 ) it2 = begin2;
  _points.push_back ( it2->getStart () );
  it1 = particle1.intersects ( *it2 );
      } while ( it1 == end1 );
      _points.push_back ( it1->intersect ( *it2 ) );
    }
    ++it1; if ( it1 == end1 ) it1 = begin1;
  } while ( it1 != newbegin );
}

template <class DataType, class IndexType>
DataType ParaParticle<DataType, IndexType>::getMovementFactor ( const aol::Vector<DataType>& velocity ) const {
  if ( velocity.size () != static_cast<IndexType> ( _points.size () ) )
    throw aol::DimensionMismatchException ( "ParaParticle::getMovementFactor, parametric particle needs correct number of velocity components", __FILE__, __LINE__ );

  const DataType CFL = 1;

  aol::Vector<DataType> lengths ( velocity.size () );
  this->getLengths ( lengths );

  int siz = velocity.size ();
  DataType factor = 1000;

  for ( int i = 0; i < siz; ++i )
    factor = std::min ( factor, abs ( CFL * std::min ( lengths [i], lengths [ ( i+siz-1 ) % siz] ) / velocity [i] ) );

  return factor;
}

template <class DataType, class IndexType>
bool ParaParticle<DataType, IndexType>::valid () const {
  //! @todo ML: Implement
  return true;
}

template <class DataType, class IndexType>
inline ParaSegmentIterator<DataType, IndexType> ParaParticle<DataType, IndexType>::beginSegment () {
  ParaSegmentIterator<DataType, IndexType> it ( *this );
  return it;
}

template <class DataType, class IndexType>
inline ParaSegmentIterator<DataType, IndexType> ParaParticle<DataType, IndexType>::endSegment () {
  ParaSegmentIterator<DataType, IndexType> it ( *this, true );
  if ( !_closed ) --it;
  return it;
}

template <class DataType, class IndexType>
inline ConstParaSegmentIterator<DataType, IndexType> ParaParticle<DataType, IndexType>::beginSegment () const {
  ConstParaSegmentIterator<DataType, IndexType> it ( *this );
  return it;
}

template <class DataType, class IndexType>
inline ConstParaSegmentIterator<DataType, IndexType> ParaParticle<DataType, IndexType>::endSegment () const {
  ConstParaSegmentIterator<DataType, IndexType> it ( *this, true );
  if ( !_closed ) --it;
  return it;
}

template <class DataType, class IndexType>
inline IndexType ParaParticle<DataType, IndexType>::beginIndex () const {
  return 0;
}

template <class DataType, class IndexType>
inline IndexType ParaParticle<DataType, IndexType>::endIndex () const {
  return size ();
}

template <class DataType, class IndexType>
ostream& ParaParticle<DataType, IndexType>::print ( ostream& os ) const {
  // Only packed and output is overloaded
  if ( this->getFormat () == Particle<ParaParticle, SegmentIteratorType, ConstSegmentIteratorType, SegmentType, ConstSegmentType, DataType>::PACKED ) {
    os << "@begin packed parametric2" << endl;
    os << size () << " " << _closed << endl;
    for ( IndexType i = 0; i < size (); ++i )
      os << _points [i] << endl;
    os << "@end packed parametric2" << endl;
  } else if ( this->getFormat () == Particle<ParaParticle, SegmentIteratorType, ConstSegmentIteratorType, SegmentType, ConstSegmentType, DataType>::BINARY ) {
    os << "@begin binary parametric2" << endl;
    IndexType s = size ();
    aol::writebinary ( os, s );
    aol::writebinary ( os, _closed );
    for ( IndexType i = 0; i < s; ++i ) {
      aol::writebinary ( os, _points [i] [0] );
      aol::writebinary ( os, _points [i] [1] );
    }
    os << "@end binary parametric2" << endl;
  } else return Particle < ParaParticle<DataType, IndexType>, ParaSegmentIterator<DataType, IndexType>, ConstParaSegmentIterator<DataType>,
                  ParaSegment<DataType, IndexType>, ConstParaSegment<DataType, IndexType>, DataType >::print ( os );

  return os;
}

template <class DataType, class IndexType>
istream& ParaParticle<DataType, IndexType>::read ( istream& is ) {
  string temp;
  // Only packed and binary input is possible
  getline ( is, temp );
  if ( temp == "@begin packed parametric" || temp == "@begin packed parametric2" ) {
    IndexType s;
    is >> s;
    if ( temp == "@begin packed parametric2" ) { is >> _closed; } else { _closed = true; }
    _points.resize ( s );
    for ( IndexType i = 0; i < s; ++i )
      is >> _points [i];
    getline ( is, temp );
    if ( temp != "@end packed parametric" && temp != "@end packed parametric2" )
      throw ( aol::FileFormatException ( "ParaParticle::read, expected end of packed parametric", __FILE__, __LINE__ ) );
  } else if ( temp == "@begin binary parametric" || temp == "@begin binary parametric2" ) {
    IndexType s;
    aol::readbinary ( is, s );
    _points.resize ( s );
    if ( temp == "@begin binary parametric2" ) { aol::readbinary ( is, _closed ); } else { _closed = true; }
    for ( IndexType i = 0; i < s; ++i ) {
      aol::readbinary ( is, _points [i] [0] );
      aol::readbinary ( is, _points [i] [1] );
    }
    getline ( is, temp );
    if ( temp != "@end binary parametric" && temp != "@end binary parametric2" )
      throw ( aol::FileFormatException ( "ParaParticle::read, expected end of binary parametric, found" + temp, __FILE__, __LINE__ ) );
  } else throw ( aol::FileFormatException ( "ParaParticle::read, can only read packed and binary format", __FILE__, __LINE__ ) );

  return is;
}

template <class DataType, class IndexType>
inline IndexType ParaParticle<DataType, IndexType>::size () const {
  return _points.size ();
}

template <class DataType, class IndexType>
inline void ParaParticle<DataType, IndexType>::smooth () {

  std::vector<aol::Vec2<DataType> > np;

  int n = size ();
  for ( int i = 0; i < n; ++i )
    np.push_back ( 0.5 * ( _points [ ( i+n-1 ) % n ] + _points [ ( i+1 ) % n ] ) );

  _points = np;
}

template <class DataType, class IndexType>
void ParaParticle<DataType, IndexType>::reparametrize ( DataType hFactor, DataType minLength ) {

  IndexType n = getNumberOfSegments ();

  // Find h_min, h_max
  DataType h_max = 0;
  DataType h_min = numeric_limits<DataType>::infinity ();
  DataType len = 0;

  for ( IndexType i = 0; i < n; ++i ) {

    aol::Vec2<DataType> dir = _points [i] - _points [ ( i + 1 ) % n];
    DataType h = dir.norm ();
    h_max = std::max ( h, h_max );
    h_min = std::min ( h, h_min );
    len += h;
  }

  if ( h_max / h_min < hFactor ) return; // No reparametrization

  int newn = std::min ( n, static_cast<int> ( len / minLength ) + 1 );
  if ( newn < 4 ) newn = 4;

  std::vector<aol::Vec2<DataType> > newpoints ( newn );

  // Reparametrize with equal arclength
  newpoints [0] = _points [0];
  DataType tau = 0, step = 0;
  for ( int i = 1, j = 1; i < newn; ++i ) {
    tau += len / newn;
    //if ( aol::debugging::refine ) cerr << "Next reparametrization step size " << tau << endl;
    while ( j <= n ) {
      step = (_points [j%n] - _points [j-1]).norm ();
      if ( step >= tau ) break;
      //if ( aol::debugging::refine ) cerr << "Skipped segment of length " << step << endl;
      tau -= step;
      ++j;
    }
    newpoints [i] = ( tau / step ) * _points [j%n] + ( 1 - tau / step ) * _points [j-1];
    //if ( aol::debugging::refine ) cerr << "Generate point " << i << " at length " << tau << " of " << step << " between points " << j-1 << " and " << j%n << endl;
  }

  _points = newpoints;

  if ( aol::debugging::refine ) cerr << "Reparametrized with respect to arclength with " << getNumberOfSegments () << " segments." << endl;
}

template <class DataType, class IndexType>
  DataType ParaParticle<DataType, IndexType>::getAngle (IndexType i) {

  IndexType n = size ();

  aol::Vec2<DataType> pp = _points [ ( i + n - 1 ) % n], cp = _points [i], np = _points [ ( i + 1 ) % n];
  aol::Vec2<DataType> pd = cp - pp, nd = np - cp;
  pd.normalize (); nd.normalize ();

  DataType angle = pd * nd;
  if ( angle >= 1 ) angle = 0;
  else angle = 360 * abs ( acos ( angle ) ) / ( 2 * aol::NumberTrait<long double>::pi );
  return angle;
}

template <class DataType, class IndexType>
void ParaParticle<DataType, IndexType>::reparametrize ( DataType minAngle, DataType maxAngle, DataType hFactor, DataType minLength ) {
  // First step: Refine and coarsen due to angles
  IndexType n = size ();
  IndexType nn = n;
  std::vector<IndexType> refine ( n );
  bool dosomething = false;

  // Find h_max
  DataType h_max = 0;
  DataType h_min = numeric_limits<DataType>::infinity ();
  for ( IndexType i = 0; i < n; ++i ) {

    aol::Vec2<DataType> dir = _points [i] - _points [ ( i + 1 ) % n];
    h_max = std::max ( dir.norm (), h_max );
    h_min = std::min ( dir.norm (), h_min );
  }

  for ( IndexType i = 0; i < n; ++i ) {

    aol::Vec2<DataType> pp = _points [ ( i + n - 1 ) % n], cp = _points [i], np = _points [ ( i + 1 ) % n];
    aol::Vec2<DataType> pd = cp - pp, nd = np - cp;
    DataType h = std::min ( pd.norm (), nd.norm () );
    pd.normalize (); nd.normalize ();

    DataType angle = pd * nd;
    if ( angle >= 1 ) angle = 0;
    else angle = 360 * abs ( acos ( angle ) ) / ( 2 * aol::NumberTrait<long double>::pi );

    // Do not refine if h < h_max / hFactor, do not coarsen in h > h_min * hFactor
    if ( angle < minAngle && h < h_min * hFactor ) {
      refine [i] = -1; --nn; dosomething = true;
    } else if ( angle > maxAngle && h > h_max / hFactor ) {
      refine [i] = + 1; nn += 2; dosomething = true;
    } else {
      refine [i] = 0;
    }
  }

  if ( dosomething ) {

    if ( aol::debugging::refine ) cerr << "Refining particle from " << n << " to " << nn << " points." << endl;

    std::vector<aol::Vec2<DataType> > newpoints ( nn );

    for ( IndexType oi = 0, ni = 0; oi < n; ++oi ) {
      if ( refine [oi] == 0 ) {
        newpoints [ni++] = _points [oi];
      }
      if ( refine [oi] == 1 ) {
        // Weight 2 : 1 so that refinement is close to sharp node (and one segment can be refined from both sides)
        newpoints [ni++] = ( _points [oi] * 2 + _points [ ( oi + n - 1 ) % n] ) / 3;
        newpoints [ni++] = _points [oi];
        newpoints [ni++] = ( _points [oi] * 2 + _points [ ( oi + 1 ) % n] ) / 3;
      }
    }

    _points = newpoints;
    n = size ();
  }

  // Second step: Collapse ends of very small segments
  IndexType src, dst;
  for ( src = 0, dst = 0; src < n; ++src, ++dst ) {

    if ( ( _points [src] - _points [ ( src+1 ) % n] ).norm () < minLength ) {
      _points [dst] = 0.5 * ( _points [src] + _points [ ( src+1 ) % n] );
      ++src;
    } else {
      _points [dst] = _points [src];
    }
  }

  // Removed last point, correct count
  if ( src > n ) { _points [0] = _points [dst - 1]; --src; --dst; }

  if ( src > dst ) {
    if ( aol::debugging::refine ) cerr << "Removed " << src - dst << " small of " << src << " total segments." << endl;
    _points.resize ( dst );
  }
}

template <class DataType, class IndexType>
void ParaParticle<DataType, IndexType>::coarsen ( IndexType delEach ) {
  if ( delEach < 2 ) delEach = 2;

  if ( ( delEach - 1 ) * size () / delEach < 3 ) return;

  std::vector<aol::Vec2<DataType> > newpoints;

  for ( IndexType i = 0; i < size (); ++i )
    if ( i % delEach != delEach - 1 ) newpoints.push_back ( _points [i] );

  _points = newpoints;
}

template <class DataType, class IndexType>
void ParaParticle<DataType, IndexType>::refine ( IndexType newPointsPerSegment ) {
  std::vector<aol::Vec2<DataType> > newpoints;
  IndexType numNewSegments = newPointsPerSegment + 1;
  for ( IndexType i = 0; i < size (); ++i ) {
    newpoints.push_back ( _points [i] );
    for ( IndexType j = 1; j < numNewSegments; ++j ) {
      DataType beta = static_cast<DataType> (j) / numNewSegments, alpha = 1 - beta;
      aol::Vec2<DataType> newpoint = alpha * _points [i] + beta * _points [(i+1)%size()];
      newpoints.push_back ( newpoint );
    }
  }

  _points = newpoints;
}

template <class DataType, class IndexType>
  void ParaParticle<DataType, IndexType>::refineSegment ( IndexType segment, IndexType newPoints, bool morePointsAtEnds ) {
  std::vector<aol::Vec2<DataType> > newpoints;
  IndexType numNewSegments = newPoints + 1;
  for ( IndexType i = 0; i < size (); ++i ) {
    newpoints.push_back ( _points [i] );
    if ( i == segment ) {
      for ( IndexType j = 1; j < numNewSegments; ++j ) {
  DataType beta = static_cast<DataType> (j) / numNewSegments;
  if ( morePointsAtEnds ) beta = 0.5 - 0.5 * cos ( aol::NumberTrait<long double>::pi * ( beta ) ); // cf. Tschebyscheff-quadrature
  DataType alpha = 1 - beta;
  aol::Vec2<DataType> newpoint = alpha * _points [i] + beta * _points [(i+1)%size()];
  newpoints.push_back ( newpoint );
      }
    }
  }

  _points = newpoints;
}

namespace {
template <class DataType> void findCenterFrom3Points (const aol::Vec2<DataType>& a, const aol::Vec2<DataType>& b, const aol::Vec2<DataType>& c, aol::Vec2<DataType>& center) {
  aol::Matrix33<DataType> M11 ( a[0], a[1], 1, b[0], b[1], 1, c[0], c[1], 1 );
  aol::Matrix33<DataType> M12 ( a.normSqr(), a[1], 1, b.normSqr(), b[1], 1, c.normSqr(), c[1], 1 );
  aol::Matrix33<DataType> M13 ( a[0], a.normSqr(), 1, b[0], b.normSqr(), 1, c[0], c.normSqr(), 1 );
  DataType d11 = M11.det(), d12 = M12.det(), d13 = M13.det ();
  center = aol::Vec2<DataType> ( 0.5 *d12 / d11, -0.5 * d13 / d11 );
}
}

template <class DataType, class IndexType>
IndexType ParaParticle<DataType, IndexType>::insertPoint ( double index, bool on_arc ) {
  IndexType N = getNumberOfSegments ();
  while (index <  0) index += N;
  while (index >= N) index -= N;
  if (index == floor (index)) index += 1e-6;
  IndexType ind_before = static_cast<int> (floor(index)), nearest_ind = static_cast<int> (index + 0.5);
  aol::Vec2<DataType> newpt;
  if (on_arc) {
    aol::Vec2<DataType> center;
    findCenterFrom3Points (getPoint ((nearest_ind+N-1)%N), getPoint (nearest_ind), getPoint ((nearest_ind+1)%N), center);
    aol::Vec2<DataType> dira = getPoint ( ind_before ) - center, dirb = getPoint ( (ind_before + 1) % N ) - center;
    DataType arca = atan2 (dira[1], dira[0]), arcb = atan2 (dirb[1], dirb[0]);
    if (arcb > arca + aol::NumberTrait<long double>::pi) arca += 2 * aol::NumberTrait<long double>::pi;
    if (arca > arcb + aol::NumberTrait<long double>::pi) arcb += 2 * aol::NumberTrait<long double>::pi;
    DataType beta = index - ind_before, alpha = 1 - beta, arc = alpha * arca + beta * arcb, radius = 0.5*(dira.norm () + dirb.norm());
    newpt = center + radius * aol::Vec2<DataType> ( cos (arc), sin (arc) );
  } else {
    double beta = index - ind_before, alpha = 1 - beta;
    newpt = alpha * getPoint ( ind_before ) + beta * getPoint ( (ind_before + 1) % N);
    // cout << "  " << index << " : " << alpha  << getPoint ( ind_before ) << " " << beta << getPoint ( (ind_before + 1) % N) << " -> " << newpt;

  }
  
  typename std::vector<aol::Vec2<DataType> >::iterator it = _points.begin ();
  it += ind_before + 1;
  _points.insert ( it, newpt );
  return ind_before + 1;
}

template <class DataType, class IndexType>
  double ParaParticle<DataType, IndexType>::getIndexForArcLength ( double arclength ) const {

  aol::Vector<DataType> len (getNumberOfSegments ());
  this->getLengths (len);

  // cerr << len << endl << arclength << " -> ";

  int seg = 0;
  while (arclength >= len [seg]) {
    arclength -= len [seg];
    ++seg;
    if (seg == getNumberOfSegments ()) seg = 0;
  }
  while (arclength < 0) {
    --seg;
    if (seg < 0) seg = getNumberOfSegments () - 1;
    arclength += len [seg];
  }

  // cerr << seg << " + " << arclength << " = " << seg + arclength / len [seg] << endl;

  return seg + arclength / len [seg];
}

template <class DataType, class IndexType>
  bool ParaParticle<DataType, IndexType>::insertPointArcLength ( double arclength, double mindist ) {

  double index = getIndexForArcLength (arclength);
  double dist = std::min (index - static_cast<int> (index), static_cast<int> (index + 1) - index);
  if (dist > mindist) {
    insertPoint (index);
    return true;
  }
  return false;
}

template <class DataType, class IndexType>
  double ParaParticle<DataType, IndexType>::getArcLengthForIndex ( double index ) const {

  aol::Vector<DataType> len (getNumberOfSegments ());
  this->getLengths (len);

  while (index < 0) index += this->getNumberOfSegments ();
  while (index >= this->getNumberOfSegments ()) index -= this->getNumberOfSegments ();

  int seg = static_cast<int> ( index );

  double arclength = 0;
  for (int i = 0; i < seg; ++i) arclength += len [i];

  arclength += len [seg] * (index - seg);
  return arclength;
}

template <class DataType, class IndexType>
inline const aol::Vec2<DataType>& ParaParticle<DataType, IndexType>::getPoint ( IndexType index ) const {
  return _points [index % size () ];
}

template <class DataType, class IndexType>
inline aol::Vec2<DataType>& ParaParticle<DataType, IndexType>::getPoint ( IndexType index ) {
  return _points [index % size () ];
}

} // namespace

#endif
