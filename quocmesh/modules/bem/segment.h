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

#ifndef __SEGMENT_H
#define __SEGMENT_H

#include <vec.h>

#include <bemesh.h>
#include <qmException.h>

namespace bm {

/******************************************************************************/
//! Segment of a polygonal particle boundary, i.e. a line segment
//! Abstract class, concrete implementation depends on the type of particle used
template <class Implementation, class DataType> class Segment {

public:

  //! Default constructor
  Segment ();

  //! Copy constructor
  explicit Segment ( const Segment& segment );

  //! Assignment operator
  Segment& operator = ( const Segment& segment );

  //! Destructor
  ~Segment ();

  //! Starting point of the line
  //! In counter-clockwise numeration as seen from the particle
  aol::Vec2<DataType> getStart () const;

  //! Ending point of the line
  //! In counter-clockwise numeration as seen from the particle
  aol::Vec2<DataType> getEnd () const;

  //! Direction vector of the line
  //! Points from start to end and its length equals the segments length
  aol::Vec2<DataType> getDirection () const;

  //! Normalized normal vector of the line
  //! Points to the right, i.e. to the outside as seen from the particle
  aol::Vec2<DataType> getNormal () const;

  //! Starting point of the line
  DataType getLength () const;

  //! Output function
  ostream& print ( ostream& os ) const;

  //! Output format
  typedef enum {
    BINARY, PACKED, PRETTY, GNUPLOT
  }
  FormatType;

  //! Change output format
  void setFormat ( FormatType format ) const;

  //! Read old output format
  FormatType getFormat () const;

  //! Compute intersection of two Segments
  //! No guarantee on output if the Segments dont intersect
  aol::Vec2<DataType> intersect ( const Segment& segment ) const;

private:

  //! Store current format type
  mutable FormatType _format;

protected:

  //! Compute the cached values of the points
  void updatePoints ();

  //! Interpret self as implementation instance
  Implementation& implementation ();

  //! Interpret self as implementation instance
  const Implementation& implementation () const;

  //! Cache start and end point
  aol::Vec2<DataType> _start, _end;
};

/***********************************************************************/
//! A segment that is described by a reference to a particle and an index
template <class Implementation, class ParticleType>
class ReferenceSegment
      : public Segment<Implementation, typename ParticleType::DataType> {

public:

  //! Types
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::IndexType IndexType;

  //! Default Constructor
  ReferenceSegment ();

  //! Constructs a side of an rectangle
  ReferenceSegment ( ParticleType& particle, IndexType index );

  //! Copy constructor
  explicit ReferenceSegment ( const ReferenceSegment& segment );

  //! Assignment operator
  ReferenceSegment& operator = ( const ReferenceSegment& segment );

  //! Destructor
  ~ReferenceSegment ();

  //! Compares two segments for equality
  //! not geometrically, but logically by its position in a particles boundary
  bool operator == ( const ReferenceSegment& segment ) const;

  //! Compares two segments for equality
  //! not geometrically, but logically by its position in a particles boundary
  bool operator != ( const ReferenceSegment& segment ) const;

  //! Is this the next segment of the same particle?
  bool isFollowedBy ( const ReferenceSegment& segment ) const;

  //! Is this another segment of the same particle?
  bool sameParticle ( const ReferenceSegment& segment ) const;

  //! Wraparound of segment numbers, call after side++ or side-- to make valid side number
  void wrap ( IndexType& side ) const;

protected:

  //! Reference to the particle
  ParticleType* _particle;
  //! Segment index
  IndexType _index;
};

//! Dicrete curvature between two segments
template <class SegmentType>
aol::Vec2<typename SegmentType::DataType> curvature ( const SegmentType& a, const SegmentType& b ) {
  aol::Vec2<typename SegmentType::DataType> adir = a.getDirection (), bdir = b.getDirection ();
  typename SegmentType::DataType alen = a.getLength (), blen = b.getLength ();

  return - 2.0 * ( bdir / blen - adir / alen ) / ( alen + blen );
}

//! Normal vector at the corner between two segments
template <class SegmentType>
aol::Vec2<typename SegmentType::DataType> cornerNormal ( const SegmentType& a, const SegmentType& b ) {

  aol::Vec2<typename SegmentType::DataType> n = - curvature ( a, b ); // "curvature" points inwards!

  // Avoid instability
  if ( n.norm () < 1E-8 ) {
    n = a.getNormal () * a.getLength() + b.getNormal () * b.getLength();
  }
    
  n.normalize ();

  if ( n * a.getNormal () < 0 && n * b.getNormal () < 0 ) n *= -1;

  return n;
}

// class Segment<Implementation,DataType>

template <class Implementation, class DataType>
inline Segment<Implementation, DataType>::Segment ()
    : _format ( PRETTY ), _start (), _end () {}

template <class Implementation, class DataType>
inline Segment<Implementation, DataType>::Segment ( const Segment& segment )
    : _format ( PRETTY ), _start ( segment._start ), _end ( segment._end ) {}

template <class Implementation, class DataType>
inline Segment<Implementation, DataType>& Segment<Implementation, DataType>::operator = ( const Segment& segment ) {
  // Beware of self-assignment
  if ( this == &segment ) return *this;

  _start = segment._start;
  _end = segment._end;
  _format = segment._format;

  return *this;
}

template <class Implementation, class DataType>
inline Segment<Implementation, DataType>::~Segment () {}

template <class Implementation, class DataType>
inline aol::Vec2<DataType> Segment<Implementation, DataType>::getStart () const {
  return _start;
}

template <class Implementation, class DataType>
inline aol::Vec2<DataType> Segment<Implementation, DataType>::getEnd () const {
  return _end;
}

template <class Implementation, class DataType>
inline aol::Vec2<DataType> Segment<Implementation, DataType>::getDirection () const {
  return getEnd () - getStart ();
}

template <class Implementation, class DataType>
inline aol::Vec2<DataType> Segment<Implementation, DataType>::getNormal () const {
  aol::Vec2<DataType> temp = getDirection ();
  temp.rotateRight ();
  temp.normalize ();
  return temp;
}

template <class Implementation, class DataType>
inline DataType Segment<Implementation, DataType>::getLength () const {
  return getDirection ().norm ();
}

template <class Implementation, class DataType>
inline ostream& Segment<Implementation, DataType>::print ( ostream& os ) const {
  switch ( getFormat () ) {

  case BINARY:
    writebinary ( os, _start );
    writebinary ( os, _end );
    break;
  case PACKED:
    os << _start << " " << _end << endl;
    break;
  case PRETTY:
    os << "| | [ " << _start << " | " << _end << " ] | |"  << endl;
    break;
  case GNUPLOT:
    os << _start << endl << _end << endl;
    break;
  default:
    throw aol::UnimplementedCodeException ( "Segment::print: unknown FormatType", __FILE__, __LINE__ );
  }

  return os;
}

template <class Implementation, class DataType>
inline ostream& operator << ( ostream& os, const Segment<Implementation, DataType>& seg ) {
  return seg.print ( os );
}

template <class Implementation, class DataType>
inline void Segment<Implementation, DataType>::setFormat ( FormatType format ) const {
  _format = format;
}

template <class Implementation, class DataType>
inline typename Segment<Implementation, DataType>::FormatType Segment<Implementation, DataType>::getFormat () const {
  return _format;
}

template <class Implementation, class DataType>
inline aol::Vec2<DataType> Segment<Implementation, DataType>::intersect ( const Segment& segment2 ) const {
  // p1 + a1 d1 = p2 + a2 d2, d1 a1 - d2 a2 = p2 - p1
  aol::Vec2<DataType> p1 = getStart (), d1 = getDirection (), p2 = segment2.getStart (), d2 = segment2.getDirection ();
  aol::Matrix22<DataType> m;
  aol::Vec2<DataType> b = p2 - p1;
  for ( int i = 0; i < 2; ++i ) {
    m [i][0] = d1[i];
    m [i][1] = - d2[i];
  }

  if ( abs ( m.det () ) < 1E-8 ) throw aol::ParameterException ( "Segment::intersect, segments collinear" , __FILE__, __LINE__ );

  aol::Matrix22<DataType> mi = m.inverse ();
  aol::Vec2<DataType> a = mi * b;

  if ( a[0] < 0 || a[0] >= 1 ) throw aol::ParameterException ( "Segment::intersect, segments do not intersect" , __FILE__, __LINE__ );
  if ( a[1] < 0 || a[1] >= 1 ) throw aol::ParameterException ( "Segment::intersect, segments do not intersect" , __FILE__, __LINE__ );

  return aol::Vec2<DataType> ( p1 + a[0] * d1 );
}

template <class Implementation, class DataType>
inline void Segment<Implementation, DataType>::updatePoints () {
  this->implementation ().updatePoints ();
}

template <class Implementation, class DataType>
inline Implementation& Segment<Implementation, DataType>::implementation () {
  return static_cast<Implementation&> ( *this );
}

template <class Implementation, class DataType>
inline const Implementation& Segment<Implementation, DataType>::implementation () const {
  return static_cast<const Implementation&> ( *this );
}

// class ReferenceSegment<Implementation,ParticleType>

template <class Implementation, class ParticleType>
inline bool ReferenceSegment<Implementation, ParticleType>::operator != ( const ReferenceSegment& segment ) const {
  return ! ( *this == segment );
}

template <class Implementation, class ParticleType>
inline bool ReferenceSegment<Implementation, ParticleType>::isFollowedBy ( const ReferenceSegment& segment ) const {

  if ( ! this->sameParticle ( segment ) ) return false;

  IndexType i = _index;

  ++i; this->implementation().wrap ( i );

  return i == segment._index;
}

template <class Implementation, class ParticleType>
inline bool ReferenceSegment<Implementation, ParticleType>::sameParticle ( const ReferenceSegment& segment ) const {

  return *_particle == *segment._particle;
}

template <class Implementation, class ParticleType>
inline void ReferenceSegment<Implementation, ParticleType>::wrap ( IndexType& side ) const {
  int skip = _particle->getNumberOfSegments ();
  union { int i; IndexType s; } curr;
  curr.s = side;

  while ( curr.i < 0 ) curr.i += skip;
  while ( curr.i >= skip ) curr.i -= skip;

  side = curr.s;
}

template <class Implementation, class ParticleType>
inline ReferenceSegment<Implementation, ParticleType>::ReferenceSegment ()
  : Segment<Implementation, DataType> (), _particle ( 0 ), _index () {}

template <class Implementation, class ParticleType>
inline ReferenceSegment<Implementation, ParticleType>::ReferenceSegment ( ParticleType& particle, IndexType index )
    : Segment<Implementation, DataType> (), _particle ( &particle ), _index ( index ) {}

template <class Implementation, class ParticleType>
inline ReferenceSegment<Implementation, ParticleType>::ReferenceSegment ( const ReferenceSegment& segment )
    : Segment<Implementation, DataType> ( segment ), _particle ( segment._particle ), _index ( segment._index ) {}

template <class Implementation, class ParticleType>
inline ReferenceSegment<Implementation, ParticleType>&
ReferenceSegment<Implementation, ParticleType>::operator = ( const ReferenceSegment& segment ) {
  // Beware of self-assignment
  if ( this == &segment ) return *this;

  Segment<Implementation, DataType>::operator = ( segment );

  _particle = segment._particle;
  _index = segment._index;

  return *this;
}

template <class Implementation, class ParticleType>
inline ReferenceSegment<Implementation, ParticleType>::~ReferenceSegment () {}

template <class Implementation, class ParticleType>
inline bool ReferenceSegment<Implementation, ParticleType>::operator == ( const ReferenceSegment& segment ) const {
  return _particle == segment._particle && _index == segment._index;
}
}
#endif

