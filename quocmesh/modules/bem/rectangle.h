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

#ifndef __RECTANGLE_H
#define __RECTANGLE_H

#include <aol.h>

#include <segment.h>
#include <particle.h>
#include <qmException.h>

namespace bm {

/********************************************************/
//! Labels for the four sides of an axis-aligned rectangle
enum Side {TOP, LEFT, BOTTOM, RIGHT, UNDEFINED};

//! Goto previous or next side
Side& operator ++ ( Side& side );
Side& operator -- ( Side& side );

// See below
template <class DataType> class RectParticle;
template <class DataType> class ConstRectSegment;

/********************************************************************/
//! Concrete representation as a side of an axis-aligned rectangle
//! No default constructor, since reference member must be initialized
//! Can be moved
template <class DataType = double>
class RectSegment
      : public ReferenceSegment<RectSegment<DataType>, RectParticle<DataType> > {

public:

  //! Default Constructor
  RectSegment ();

  //! Constructs a side of an rectangle
  RectSegment ( RectParticle<DataType>& rectangle, Side side );

  //! Copy constructor
  RectSegment ( const RectSegment& segment );

  //! Assignment operator
  RectSegment<DataType>& operator = ( const RectSegment& segment );

  //! Destructor
  ~RectSegment ();

  //! Move the segment a short step in the direction of the normal
  void move ( DataType step );

  //! Wraparound of segment numbers, call after side++ or side-- to make valid side number
  void wrap ( Side& side ) const;

protected:

  //! Compute the cached values of the points
  void updatePoints ();
};

/********************************************************************/
//! Concrete representation as a side of an axis-aligned rectangle
//! No default constructor, since reference member must be initialized
//! Cannot be moved
template <class DataType = double>
class ConstRectSegment
      : public ReferenceSegment<ConstRectSegment<DataType>, const RectParticle<DataType> > {

public:

  //! Default Constructor
  ConstRectSegment ();

  //! Constructs a side of an rectangle
  ConstRectSegment ( const RectParticle<DataType>& rectangle, Side side );

  //! Copy constructor
  ConstRectSegment ( const ConstRectSegment<DataType>& segment );

  //! Assignment operator
  ConstRectSegment<DataType>& operator = ( const ConstRectSegment<DataType>& segment );

  //! Destructor
  ~ConstRectSegment ();

protected:

  //! Compute the cached values of the points
  void updatePoints ();
};

// See below
template <class DataType> class ConstRectSegmentIterator;

/*******************************************************************/
//! Enumerates the four sides of the rectangle, allows to change them
template <class DataType = double>
class RectSegmentIterator
      : public ReferenceSegmentIterator<RectParticle<DataType>, RectSegment<DataType> > {

public:

  //! Default constructor
  RectSegmentIterator ();

  //! Constructor for both sentinels
  RectSegmentIterator ( RectParticle<DataType>& particle, bool end = false );

  //! Copy constructor
  RectSegmentIterator ( const RectSegmentIterator& it );

  //! Assignment Operator
  RectSegmentIterator& operator = ( const RectSegmentIterator& it );

  //! Destructor
  ~RectSegmentIterator ();

  //! Dereference operators
  RectSegment<DataType>* operator -> ();
  RectSegment<DataType>& operator * ();
};

/********************************************/
//! Enumerates the four sides of the rectangle
template <class DataType = double>
class ConstRectSegmentIterator
      : public ReferenceSegmentIterator<const RectParticle<DataType>, ConstRectSegment<DataType> > {

public:

  //! Default constructor
  ConstRectSegmentIterator ();

  //! Constructor for both sentinels
  ConstRectSegmentIterator ( const RectParticle<DataType>& particle, bool end = false );

  //! Copy constructor
  ConstRectSegmentIterator ( const ConstRectSegmentIterator& it );

  //! Assignment Operator
  ConstRectSegmentIterator& operator = ( const ConstRectSegmentIterator& it );

  //! Destructor
  ~ConstRectSegmentIterator ();

  //! Dereference operators
  const ConstRectSegment<DataType>* operator -> () const;
  const ConstRectSegment<DataType>& operator * () const;
};

/**********************/
//! Rectangular particle
template <class _DataType = double>
class RectParticle
      : public Particle < RectParticle<_DataType>, RectSegmentIterator<_DataType>, ConstRectSegmentIterator<_DataType>,
      RectSegment<_DataType>, ConstRectSegment<_DataType>, _DataType > {

public:

  // Configurator
  typedef Side IndexType;
  typedef _DataType DataType;
  typedef RectSegment<DataType> SegmentType;
  typedef ConstRectSegment<DataType> ConstSegmentType;
  typedef RectSegmentIterator<DataType> SegmentIteratorType;
  typedef ConstRectSegmentIterator<DataType> ConstSegmentIteratorType;

public:

  //! Default constructor
  RectParticle ();

  //! Creates a particle from two points
  RectParticle ( aol::Vec2<DataType> ll, aol::Vec2<DataType> ur );

  //! Creates a particle from center and radii
  RectParticle ( aol::Vec2<DataType> center, DataType a, DataType b );

  //! Copy constructor
  RectParticle ( const RectParticle& particle );

  //! Assignment operator
  RectParticle<DataType>& operator = ( const RectParticle& particle );

  //! Destructor
  ~RectParticle ();

  //! Comparison operators
  bool operator == ( const RectParticle& particle ) const;
  bool operator != ( const RectParticle& particle ) const;

  //! Read corner points of rectangle
  const aol::Vec2<DataType>& getLowerLeft () const;
  const aol::Vec2<DataType>& getUpperRight () const;

  //! Number of elements
  int getNumberOfSegments () const;

  //! Move one side of the particle
  void moveBoundary ( Side side, DataType step );

  //! Move the particle
  void move ( const aol::Vector<DataType>& velocity );

  //! Area of the Particle
  DataType getVolume () const;

  //! Surface of the Particle
  DataType getSurface () const;

  //! Average Radius of the Particle
  DataType getAverageRadius () const;

  //! Center of the particle
  aol::Vec2<DataType> getCenter () const;

  //! Compute the bounding box
  void getBoundingBox ( aol::Vec2<DataType>& ll, aol::Vec2<DataType>& ur ) const;

  //! Start SegmentIterator
  RectSegmentIterator<DataType> beginSegment ();

  //! End sentinel
  RectSegmentIterator<DataType> endSegment ();

  //! Start Constant SegmentIterator
  ConstRectSegmentIterator<DataType> beginSegment () const;

  //! Constant end sentinel
  ConstRectSegmentIterator<DataType> endSegment () const;

  //! First index value
  Side beginIndex () const;

  //! End sentinel index value
  Side endIndex () const;

  //! Collision detection, with a given buffer distance
  bool collides ( const RectParticle& particle, DataType buffer = 0 ) const;

  //! Find if particle exceeds given bounding box
  bool exceeds ( aol::Vec2<DataType> ll, aol::Vec2<DataType> ur ) const;

  //! Find if point is inside the particle
  bool inside ( aol::Vec2<DataType> point ) const;

  //! Unify two particles
  //! Must work for colliding particles
  //! Resulting particle need not be subset of the set theoretical union
  void absorb ( const RectParticle& particle );

  //! Check by which fraction of the velocity the particle can be moved without removing particles
  DataType getMovementFactor ( const aol::Vector<DataType>& velocity ) const;

  //! Checks for validity
  bool valid () const;

  //! Output function
  ostream& print ( ostream& os ) const;

  //! Input function
  istream& read ( istream& is );

  //! Number of segments
  int size () const;

  //! Reparametrize particle, does nothing for rectagnle
  void reparametrize ( DataType minAngle, DataType maxAngle, DataType hFactor );

private:

  //! Two points defining the rectangle
  aol::Vec2<DataType> _lowerLeft, _upperRight;
};

// enum Side

inline Side& operator ++ ( Side& side ) {
  // Important: Union avoids aliasing problem!
  union { int i; Side s; } uni;
  uni.s = side;
  ++ uni.i;
  side = uni.s;
  return side;
}

inline Side& operator -- ( Side& side ) {
  // Important: Union avoids aliasing problem!
  union { int i; Side s; } uni;
  uni.s = side;
  -- uni.i;
  side = uni.s;
  return side;
}

// class RectSegment<DataType>

template <class DataType>
inline RectSegment<DataType>::RectSegment ()
    : ReferenceSegment<RectSegment<DataType>, RectParticle<DataType> > () {}

template <class DataType>
inline RectSegment<DataType>::RectSegment ( RectParticle<DataType>& particle, Side index )
    : ReferenceSegment<RectSegment<DataType>, RectParticle<DataType> > ( particle, index ) {
  this->updatePoints ();
}

template <class DataType>
inline RectSegment<DataType>::RectSegment ( const RectSegment& segment )
    : ReferenceSegment<RectSegment<DataType>, RectParticle<DataType> > ( segment ) {}

template <class DataType>
inline RectSegment<DataType>& RectSegment<DataType>::operator = ( const RectSegment& segment ) {
  // Beware of self-assignment
  if ( this == &segment ) return *this;

  ReferenceSegment<RectSegment<DataType>, RectParticle<DataType> >::operator = ( segment );

  return *this;
}

template <class DataType>
inline RectSegment<DataType>::~RectSegment () {}

template <class DataType>
inline void RectSegment<DataType>::move ( DataType step ) {
  this->_particle->moveBoundary ( this->_index, step );

  this->updatePoints ();
}

template <class DataType>
inline void RectSegment<DataType>::wrap ( Side& side ) const {
  union { int i; Side s; } curr, skip;
  curr.s = side;
  skip.s = UNDEFINED;

  while ( curr.s < TOP ) curr.i += skip.i;
  while ( curr.s > RIGHT ) curr.i -= skip.i;

  side = curr.s;
}

template <class DataType>
void RectSegment<DataType>::updatePoints () {
  aol::Vec2<DataType> ll = this->_particle->getLowerLeft ();
  aol::Vec2<DataType> ur = this->_particle->getUpperRight ();

  DataType top = ur.y ();
  DataType left = ll.x ();
  DataType bottom = ll.y ();
  DataType right = ur.x ();

  switch ( this->_index ) {

  case TOP:
    this->_start = aol::Vec2<DataType> ( right, top );
    this->_end = aol::Vec2<DataType> ( left, top );
    break;
  case LEFT:
    this->_start = aol::Vec2<DataType> ( left, top );
    this->_end = aol::Vec2<DataType> ( left, bottom );
    break;
  case BOTTOM:
    this->_start = aol::Vec2<DataType> ( left, bottom );
    this->_end = aol::Vec2<DataType> ( right, bottom );
    break;
  case RIGHT:
    this->_start = aol::Vec2<DataType> ( right, bottom );
    this->_end = aol::Vec2<DataType> ( right, top );
    break;
  case UNDEFINED:
    break;
  }
}

// class ConstRectSegment<DataType>

template <class DataType>
inline ConstRectSegment<DataType>::ConstRectSegment ()
  : ReferenceSegment<ConstRectSegment<DataType>, const RectParticle<DataType> > () {}

template <class DataType>
inline ConstRectSegment<DataType>::ConstRectSegment ( const RectParticle<DataType>& particle, Side index )
    : ReferenceSegment<ConstRectSegment<DataType>, const RectParticle<DataType> > ( particle, index ) {
  updatePoints ();
}

template <class DataType>
inline ConstRectSegment<DataType>::ConstRectSegment ( const ConstRectSegment& segment )
    : ReferenceSegment<ConstRectSegment<DataType>, const RectParticle<DataType> > ( segment ) {}

template <class DataType>
inline ConstRectSegment<DataType>& ConstRectSegment<DataType>::operator = ( const ConstRectSegment& segment ) {
  // Beware of self-assignment
  if ( this == &segment ) return *this;

  ReferenceSegment<ConstRectSegment<DataType>, const RectParticle<DataType> >::operator = ( segment );

  return *this;
}

template <class DataType>
inline ConstRectSegment<DataType>::~ConstRectSegment () {}

template <class DataType>
void ConstRectSegment<DataType>::updatePoints () {
  aol::Vec2<DataType> ll = this->_particle->getLowerLeft ();
  aol::Vec2<DataType> ur = this->_particle->getUpperRight ();

  DataType top = ur.y ();
  DataType left = ll.x ();
  DataType bottom = ll.y ();
  DataType right = ur.x ();

  switch ( this->_index ) {

  case TOP:
    this->_start = aol::Vec2<DataType> ( right, top );
    this->_end = aol::Vec2<DataType> ( left, top );
    break;
  case LEFT:
    this->_start = aol::Vec2<DataType> ( left, top );
    this->_end = aol::Vec2<DataType> ( left, bottom );
    break;
  case BOTTOM:
    this->_start = aol::Vec2<DataType> ( left, bottom );
    this->_end = aol::Vec2<DataType> ( right, bottom );
    break;
  case RIGHT:
    this->_start = aol::Vec2<DataType> ( right, bottom );
    this->_end = aol::Vec2<DataType> ( right, top );
    break;
  case UNDEFINED:
    break;
  default:
    throw aol::Exception ( "bem::ConstRectSegment::updatePoints: unknown SideType", __FILE__, __LINE__ );
  }
}

// Class RectSegmentIterator<DataType>

template <class DataType>
inline RectSegmentIterator<DataType>::RectSegmentIterator ()
  : ReferenceSegmentIterator<RectParticle<DataType>, RectSegment<DataType> > () {}

template <class DataType>
inline RectSegmentIterator<DataType>::RectSegmentIterator ( RectParticle<DataType>& particle, bool end )
    : ReferenceSegmentIterator<RectParticle<DataType>, RectSegment<DataType> > ( particle, end ) {}

template <class DataType>
inline RectSegmentIterator<DataType>::RectSegmentIterator ( const RectSegmentIterator& it )
    : ReferenceSegmentIterator<RectParticle<DataType>, RectSegment<DataType> > ( it ) {}

template <class DataType>
inline RectSegmentIterator<DataType>& RectSegmentIterator<DataType>::operator = ( const RectSegmentIterator& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  ReferenceSegmentIterator<RectParticle<DataType>, RectSegment<DataType> >::operator = ( it );

  this->updateSegment ();

  return *this;
}

template <class DataType>
inline RectSegmentIterator<DataType>::~RectSegmentIterator () {}

template <class DataType>
inline RectSegment<DataType>* RectSegmentIterator<DataType>::operator -> () {
  return &this->_segment;
}

template <class DataType>
inline RectSegment<DataType>& RectSegmentIterator<DataType>::operator * () {
  return this->_segment;
}

// Class ConstRectSegmentIterator<DataType>

template <class DataType>
inline ConstRectSegmentIterator<DataType>::ConstRectSegmentIterator ()
    : ReferenceSegmentIterator<const RectParticle<DataType>, ConstRectSegment<DataType> > () {}

template <class DataType>
inline ConstRectSegmentIterator<DataType>::ConstRectSegmentIterator ( const RectParticle<DataType>& particle, bool end )
    : ReferenceSegmentIterator<const RectParticle<DataType>, ConstRectSegment<DataType> > ( particle, end ) {}

template <class DataType>
inline ConstRectSegmentIterator<DataType>::ConstRectSegmentIterator ( const ConstRectSegmentIterator& it )
    : ReferenceSegmentIterator<const RectParticle<DataType>, ConstRectSegment<DataType> > ( it ) {}

template <class DataType>
inline ConstRectSegmentIterator<DataType>& ConstRectSegmentIterator<DataType>::operator = ( const ConstRectSegmentIterator& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  ReferenceSegmentIterator<const RectParticle<DataType>, ConstRectSegment<DataType> >::operator = ( it );

  this->updateSegment ();

  return *this;
}

template <class DataType>
inline ConstRectSegmentIterator<DataType>::~ConstRectSegmentIterator () {}

template <class DataType>
inline const ConstRectSegment<DataType>* ConstRectSegmentIterator<DataType>::operator -> () const {
  return &this->_segment;
}

template <class DataType>
inline const ConstRectSegment<DataType>& ConstRectSegmentIterator<DataType>::operator * () const {
  return this->_segment;
}

// Class RectParticle<DataType>

template <class DataType>
inline RectParticle<DataType>::RectParticle ()
    : Particle < RectParticle<DataType>, RectSegmentIterator<DataType>, ConstRectSegmentIterator<DataType>,
    RectSegment<DataType>, ConstRectSegment<DataType>, DataType > (),
    _lowerLeft (), _upperRight () {}

template <class DataType>
inline RectParticle<DataType>::RectParticle ( aol::Vec2<DataType> ll, aol::Vec2<DataType> ur )
    : Particle < RectParticle<DataType>, RectSegmentIterator<DataType>, ConstRectSegmentIterator<DataType>,
    RectSegment<DataType>, ConstRectSegment<DataType>, DataType > (),
    _lowerLeft ( ll ), _upperRight ( ur ) {
  if ( ll.x () >= ur.x () || ll.y () >= ur.y () )
    throw aol::ParameterException ( "RectParticle::RectParticle, lower left point of rectangle should be first", __FILE__, __LINE__ );
}

template <class DataType>
inline RectParticle<DataType>::RectParticle ( aol::Vec2<DataType> center, DataType a, DataType b )
    : Particle < RectParticle<DataType>, RectSegmentIterator<DataType>, ConstRectSegmentIterator<DataType>,
    RectSegment<DataType>, ConstRectSegment<DataType>, DataType > (),
    _lowerLeft ( center.x () - a, center.y () - b ), _upperRight ( center.x () + a, center.y () + b ) {
  if ( a <= 0 || b <= 0 )
    throw aol::ParameterException ( "RectParticle::RectParticle, rectangle must have positive radii", __FILE__, __LINE__ );
}

template <class DataType>
inline RectParticle<DataType>::RectParticle ( const RectParticle& particle )
    : Particle < RectParticle<DataType>, RectSegmentIterator<DataType>, ConstRectSegmentIterator<DataType>,
    RectSegment<DataType>, ConstRectSegment<DataType>, DataType > ( particle ),
  _lowerLeft ( particle._lowerLeft ), _upperRight ( particle._upperRight ) {}

template <class DataType>
inline RectParticle<DataType>& RectParticle<DataType>::operator = ( const RectParticle& particle ) {
  // Beware of self-assignment
  if ( this == &particle ) return *this;

  Particle < RectParticle<DataType>, RectSegmentIterator<DataType>, ConstRectSegmentIterator<DataType>,
  RectSegment<DataType>, ConstRectSegment<DataType>, DataType >::operator =
  ( particle );

  _lowerLeft = particle._lowerLeft;
  _upperRight = particle._upperRight;

  return *this;
}

template <class DataType>
inline RectParticle<DataType>::~RectParticle () {}

template <class DataType>
inline bool RectParticle<DataType>::operator == ( const RectParticle& particle ) const {
  return _lowerLeft == particle._lowerLeft && _upperRight == particle._upperRight;
}

template <class DataType>
inline bool RectParticle<DataType>::operator != ( const RectParticle& particle ) const {
  return ! ( *this == particle );
}

template <class DataType>
inline const aol::Vec2<DataType>& RectParticle<DataType>::getLowerLeft () const {
  return _lowerLeft;
}

template <class DataType>
inline const aol::Vec2<DataType>& RectParticle<DataType>::getUpperRight () const {
  return _upperRight;
}

template <class DataType>
inline int RectParticle<DataType>::getNumberOfSegments () const {
  return 4;
}

template <class DataType>
void RectParticle<DataType>::moveBoundary ( Side side, DataType step ) {
  switch ( side ) {
  case TOP:
    _upperRight.yref () += step;
    break;
  case LEFT:
    _lowerLeft.xref () -= step;
    break;
  case BOTTOM:
    _lowerLeft.yref () -= step;
    break;
  case RIGHT:
    _upperRight.xref () += step;
    break;
  case UNDEFINED:
    break;
  }
}

template <class DataType>
void RectParticle<DataType>::move ( const aol::Vector<DataType>& velocity ) {
  if ( velocity.size () != 4 ) throw aol::DimensionMismatchException ( "RectParticle::move, Rectangular particle needs four velocity components", __FILE__, __LINE__ );
  _upperRight.yref () += velocity [0];
  _lowerLeft.xref () -= velocity [1];
  _lowerLeft.yref () -= velocity [2];
  _upperRight.xref () += velocity [3];
}

template <class DataType>
DataType RectParticle<DataType>::getVolume () const {
  DataType dx = _upperRight.x () - _lowerLeft.x ();
  DataType dy = _upperRight.y () - _lowerLeft.y ();

  if ( dx <= 0 || dy <= 0 ) return 0;

  return dx * dy;
}

template <class DataType>
DataType RectParticle<DataType>::getSurface () const {
  DataType dx = _upperRight.x () - _lowerLeft.x ();
  DataType dy = _upperRight.y () - _lowerLeft.y ();

  if ( dx <= 0 || dy <= 0 ) return 0;

  return 2 * ( dx + dy );
}

template <class DataType>
DataType RectParticle<DataType>::getAverageRadius () const {
  return getSurface () / 8;
}

template <class DataType>
aol::Vec2<DataType> RectParticle<DataType>::getCenter () const {
  return aol::Vec2<DataType> ( ( _lowerLeft + _upperRight ) / 2 );
}

template <class DataType>
void RectParticle<DataType>::getBoundingBox ( aol::Vec2<DataType>& ll, aol::Vec2<DataType>& ur ) const {
  ll = _lowerLeft;
  ur = _upperRight;
}

template <class DataType>
bool RectParticle<DataType>::collides ( const RectParticle& particle, DataType buffer ) const {
  return particle._lowerLeft [0] - buffer <= _upperRight [0]
         && particle._lowerLeft [1] - buffer <= _upperRight [1]
         && particle._upperRight [0] + buffer >= _lowerLeft [0]
         && particle._upperRight [1] + buffer >= _lowerLeft [1];
}

template <class DataType>
bool RectParticle<DataType>::exceeds ( aol::Vec2<DataType> ll, aol::Vec2<DataType> ur ) const {
  return _lowerLeft [0] < ll [0]
         || _lowerLeft [1] < ll [1]
         || _upperRight [0] > ur [0]
         || _upperRight [1] > ur [1];
}

template <class DataType>
bool RectParticle<DataType>::inside ( aol::Vec2<DataType> point ) const {
  return _lowerLeft [0] < point [0]
         && _lowerLeft [1] < point [1]
         && _upperRight [0] > point [0]
         && _upperRight [1] > point [1];
}

template <class DataType>
void RectParticle<DataType>::absorb ( const RectParticle& particle ) {

  DataType volume = this->getVolume () + particle.getVolume ();

  aol::Vec2<DataType> newll ( std::min ( this->_lowerLeft [0], particle._lowerLeft [0] ),
            std::min ( this->_lowerLeft [1], particle._lowerLeft [1] ) ),
                      newur ( std::max ( this->_upperRight [0], particle._upperRight [0] ),
            std::max ( this->_upperRight [1], particle._upperRight [1] ) ),
                      center = 0.5 * ( newll + newur );

  DataType newvol = ( newur [0] - newll [0] ) * ( newur [1] - newll [1] ), fac = sqrt ( volume / newvol );

  _lowerLeft  = center + fac * ( newll - center );
  _upperRight = center + fac * ( newur - center );

  if ( abs ( ( volume - this->getVolume () ) / volume ) > 1E-8 ) cerr << "Unification Error!" << endl;
}

template <class DataType>
DataType RectParticle<DataType>::getMovementFactor ( const aol::Vector<DataType>& velocity ) const {
  if ( velocity.size () != 4 )
    throw aol::DimensionMismatchException ( "RectParticle<DataType>::getMovementFactor, rectangular particle needs four velocity components", __FILE__, __LINE__ );

  // Maximum relaxation allowed is to vanishing particle

  // Pointing outward
  DataType hvel = velocity [1] + velocity [3];
  DataType vvel = velocity [0] + velocity [2];

  DataType hsize = _upperRight [0] - _lowerLeft [0];
  DataType vsize = _upperRight [1] - _lowerLeft [1];

  DataType hf = abs ( hsize / hvel );
  DataType vf = abs ( vsize / vvel );

  return std::min ( hf, vf );
}

template <class DataType>
bool RectParticle<DataType>::valid () const {
  return _lowerLeft [0] < _upperRight [0] && _lowerLeft [1] < _upperRight [1];
}

template <class DataType>
inline RectSegmentIterator<DataType> RectParticle<DataType>::beginSegment () {
  RectSegmentIterator<DataType> it ( *this );
  return it;
}

template <class DataType>
inline RectSegmentIterator<DataType> RectParticle<DataType>::endSegment () {
  RectSegmentIterator<DataType> it ( *this, true );
  return it;
}

template <class DataType>
inline ConstRectSegmentIterator<DataType> RectParticle<DataType>::beginSegment () const {
  ConstRectSegmentIterator<DataType> it ( *this );
  return it;
}

template <class DataType>
inline ConstRectSegmentIterator<DataType> RectParticle<DataType>::endSegment () const {
  ConstRectSegmentIterator<DataType> it ( *this, true );
  return it;
}

template <class DataType>
inline Side RectParticle<DataType>::beginIndex () const {
  return TOP;
}

template <class DataType>
inline Side RectParticle<DataType>::endIndex () const {
  return UNDEFINED;
}

template <class DataType>
ostream& RectParticle<DataType>::print ( ostream& os ) const {
  // Only packed and output is overloaded
  if ( this->getFormat () == Particle<RectParticle, SegmentIteratorType, ConstSegmentIteratorType, SegmentType, ConstSegmentType, DataType>::PACKED ) {
    os << "@begin packed rectangle" << endl;
    os << _lowerLeft << endl;
    os << _upperRight << endl;
    os << "@end packed rectangle" << endl;
  } else if ( this->getFormat () == Particle<RectParticle, SegmentIteratorType, ConstSegmentIteratorType, SegmentType, ConstSegmentType, DataType>::BINARY ) {
    os << "@begin binary rectangle" << endl;
    writebinary ( os, _lowerLeft );
    writebinary ( os, _upperRight );
    os << "@end binary rectangle" << endl;
  } else return Particle < RectParticle<DataType>, RectSegmentIterator<DataType>,
                  ConstRectSegmentIterator<DataType>, RectSegment<DataType>, ConstRectSegment<DataType>, DataType >::print ( os );

  return os;
}

template <class DataType>
istream& RectParticle<DataType>::read ( istream& is ) {
  string temp;
  // Only packed and binary input is possible
  getline ( is, temp );
  if ( temp == "@begin packed rectangle" ) {
    is >> _lowerLeft;
    is >> _upperRight;
    getline ( is, temp );
    if ( temp != "@end packed rectangle" )
      throw ( aol::FileFormatException ( "RectParticle::read, expected end of packed rectangle", __FILE__, __LINE__ ) );
  } else if ( temp == "@begin binary rectangle" ) {
    readbinary ( is, _lowerLeft );
    readbinary ( is, _upperRight );
    getline ( is, temp );
    if ( temp != "@end binary rectangle" )
      throw ( aol::FileFormatException ( "RectParticle::read, expected end of binary rectangle", __FILE__, __LINE__ ) );
  } else throw ( aol::FileFormatException ( "RectParticle::read, can only read packed and binary format", __FILE__, __LINE__ ) );

  return is;
}

template <class DataType>
inline int RectParticle<DataType>::size () const {
  return 4;
}

template <class DataType>
inline void RectParticle<DataType>::reparametrize ( DataType /*minAngle*/, DataType /*maxAngle*/, DataType /*hFactor*/ ) {
  // No reparametriziation necessary
}
}

#endif
