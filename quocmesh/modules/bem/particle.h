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

#ifndef __PARTICLE_H
#define __PARTICLE_H

#include <bemesh.h>
#include <qmException.h>

namespace bm {

/*******************************************************************************/
//! Particle boundary described by a closed polygon, i.e. an iterator of segments
//! Abstract with respect to the real storage of data
template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
class Particle {

public:

  enum DOFPOSITION { START, CENTER };

  //! Default constructor
  Particle ();

  //! Copy constructor
  explicit Particle ( const Particle& particle );

  //! Assignment operator
  Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>& operator = ( const Particle& particle );

  //! Destructor
  virtual ~Particle ();

  //! Start SegmentIterator
  IteratorType beginSegment ();

  //! End sentinel
  IteratorType endSegment ();

  //! Start Constant SegmentIterator
  ConstIteratorType beginSegment () const;

  //! Constant end sentinel
  ConstIteratorType endSegment () const;

  //! Number of elements
  int getNumberOfSegments () const;

  //! Lengths of all Segments
  void getLengths ( aol::Vector<DataType>& lengths, bool betweencenters = false, bool jumping = false ) const;

  //! Normals on all Segments
  void getNormals ( aol::MultiVector<DataType>& normals, bool onFaces = false ) const;

  //! Positions of all start points
  void getPositions ( aol::MultiVector<DataType>& positions, DOFPOSITION locpos = START ) const;

  //! Curvatures on all Segments
  void getCurvatures ( aol::Vector<DataType>& curvatures, bool onFaces = false ) const;

  //! Area of the Particle
  DataType getVolume () const;

  //! Surface of the Particle
  DataType getLength () const;

  //! Output function
  virtual ostream& print ( ostream& os ) const;

  //! Output format
  typedef enum {
    BINARY, PACKED, PRETTY, GNUPLOT
  }
  FormatType;

  //! Input function, only available on concrete classes
  istream& read ( istream& is );

  //! Change output format
  void setFormat ( FormatType format ) const;

  //! Read old output format
  FormatType getFormat () const;

protected:

  //! Interpret self as implementation instance
  Implementation& implementation ();

  //! Interpret self as implementation instance
  const Implementation& implementation () const;

private:

  //! Store current format type
  mutable FormatType _format;
};

/**************************************************************************************/
//! Iterator that references a particle and a segment index
//! Const-ness of dereference operators depends on the const-ness of template parameters
template <class ParticleType, class SegmentType>
class ReferenceSegmentIterator
      : public aol::Iterator<SegmentType> {

public:

  //! Types
  typedef typename ParticleType::IndexType IndexType;

  //! Default constructor
  ReferenceSegmentIterator ();

  //! Constructor for both sentinels
  ReferenceSegmentIterator ( ParticleType& particle, bool end = false );

  //! Copy constructor
  ReferenceSegmentIterator ( const ReferenceSegmentIterator& it );

  //! Assignment operator
  ReferenceSegmentIterator<ParticleType, SegmentType>& operator = ( const ReferenceSegmentIterator& it );

  //! Destructor
  ~ReferenceSegmentIterator ();

  //! Comparison operator
  bool operator == ( const ReferenceSegmentIterator& it ) const;
  bool operator != ( const ReferenceSegmentIterator& it ) const;

  //! Prefix operator ++
  ReferenceSegmentIterator& operator ++ ();

  //! Prefix operator --
  ReferenceSegmentIterator& operator -- ();

protected:

  //! Compute cached value of current segment
  void updateSegment ();

  //! Reference to particle
  ParticleType* _particle;
  //! Reference to side number
  IndexType _index;

  //! Cache current segment
  SegmentType _segment;
};

// Class Particle<Implementation,IteratorType,ConstIteratorType,SegmentType,ConstSegmentType,DataType>

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::Particle ()
    : _format ( PRETTY ) {}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::Particle ( const Particle& )
    : _format ( PRETTY ) {}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>&
Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::operator = ( const Particle& particle ) {
  // Beware of self-assignment
  if ( this == &particle ) return *this;

  _format = particle._format;

  return *this;
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::~Particle () {}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline void Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::getLengths ( aol::Vector<DataType>& lengths,
    bool betweencenters, bool jumping ) const {
  int i = 0, size = lengths.size ();

  if ( size != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::getLengths, vector length does not equal number of segments", __FILE__, __LINE__ );

  lengths.setZero ();

  for ( ConstIteratorType it = beginSegment (); it != endSegment (); ++it, ++i ) {
    DataType l = it->getLength ();
    if ( jumping && i == size - 1 ) l = 0;
    if ( betweencenters ) {
      lengths [i] += l / 2; lengths [ ( i+1 ) % size] += l / 2;
    } else lengths [i] = l;
  }
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline void Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::getPositions ( aol::MultiVector<DataType>& positions, DOFPOSITION locpos ) const {
  int i = 0;

  if ( positions.numComponents () != 2
       || positions [0].size () != getNumberOfSegments ()
       || positions [1].size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::getPositions, vector length does not equal number of segments", __FILE__, __LINE__ );

  for ( ConstIteratorType it = beginSegment (); it != endSegment (); ++it, ++i ) {
    aol::Vec2<DataType> x = it->getStart ();
    if ( locpos == CENTER ) { x += it->getEnd (); x *= 0.5; }
    positions [0][i] = x [0];
    positions [1][i] = x [1];
  }
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline void Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::getCurvatures ( aol::Vector<DataType>& curvatures, bool onFaces ) const {
  int i;

  if ( curvatures.size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::getCurvatures, vector length does not equal number of segments", __FILE__, __LINE__ );
  ConstIteratorType last = endSegment ();
  ConstIteratorType curr = beginSegment ();
  ConstIteratorType next = curr;
  ConstIteratorType end = last;
  --last; ++next;

  for ( i = 0; i < getNumberOfSegments (); last = curr, curr = next, ++next, ++i ) {
    if ( next == end ) next = beginSegment ();

    aol::Vec2<DataType> kappa, normal, normaldir;
    if ( onFaces ) {
      kappa = ( curvature ( *last, *curr ) + curvature ( *curr, *next ) ) / 2;
      normal = curr->getNormal ();
    } else {
      kappa = curvature ( *last, *curr );
      normal = cornerNormal ( *last, *curr );
      normal.normalize ();
    }
    curvatures [i] = kappa * normal;
  }
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline void Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::getNormals ( aol::MultiVector<DataType>& normals, bool onFaces ) const {
  int i = 0;

  if ( normals.numComponents () != 2
       || normals [0].size () != getNumberOfSegments ()
       || normals [1].size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::getNormals, vector length does not equal number of segments", __FILE__, __LINE__ );

  ConstIteratorType last = endSegment ();
  ConstIteratorType curr = beginSegment ();
  ConstIteratorType next = curr;
  ConstIteratorType end = last;
  --last; ++next;

  for ( i = 0; i < getNumberOfSegments (); last = curr, curr = next, ++next, ++i ) {
    if ( next == end ) next = beginSegment ();

    aol::Vec2<DataType> normal;
    if ( onFaces ) {
      normal = curr->getNormal ();
    } else {
      normal = cornerNormal ( *last, *curr );
    }
    normals [0][i] = normal [0];
    normals [1][i] = normal [1];
  }
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
ostream& Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::print ( ostream& os ) const {
  switch ( getFormat () ) {

  case BINARY:
    os << "@begin binary particle" << endl;
    for ( ConstIteratorType it = beginSegment (); it != endSegment (); ++it ) {
      it->setFormat ( ConstSegmentType::BINARY );
      os << *it;
    }
    os << "@end binary particle" << endl;
    break;

  case PACKED:
    os << "@begin packed particle" << endl;
    for ( ConstIteratorType it = beginSegment (); it != endSegment (); ++it ) {
      it->setFormat ( ConstSegmentType::PACKED );
      os << *it;
    }
    os << "@end packed particle" << endl;
    break;

  case PRETTY:
    os << "| /====================== Particle =====================\\ |" << endl;
    for ( ConstIteratorType it = beginSegment (); it != endSegment (); ++it ) {
      it->setFormat ( ConstSegmentType::PRETTY );
      os << *it;
    }
    os << "| \\=====================================================/ |" << endl;
    os << "|                                                         |" << endl;
    break;

  case GNUPLOT:
    for ( ConstIteratorType it = beginSegment (); it != endSegment (); ++it ) {
      it->setFormat ( ConstSegmentType::GNUPLOT );
      os << *it;
    }
    break;
    
  default:
      throw aol::UnimplementedCodeException ( "Particle::print: unknown FormatType", __FILE__, __LINE__ );
  }

  return os;
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
ostream& operator << ( ostream& os, const Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>& part ) {
  return part.print ( os );
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
istream& operator >> ( istream& is, Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>& part ) {
  return part.read ( is );
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline void Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::setFormat ( FormatType format ) const {
  _format = format;
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline typename Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::FormatType
Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::getFormat () const {
  return _format;
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline IteratorType Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::beginSegment () {
  return implementation ().beginSegment ();
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline IteratorType Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::endSegment () {
  return implementation ().endSegment ();
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline ConstIteratorType Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::beginSegment () const {
  return implementation ().beginSegment ();
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline ConstIteratorType Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::endSegment () const {
  return implementation ().endSegment ();
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline int Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::getNumberOfSegments () const {
  return implementation ().getNumberOfSegments ();
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline DataType Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::getVolume () const {
  return implementation ().getVolume ();
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline DataType Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::getLength () const {
  aol::Vector<DataType> temp ( this->getNumberOfSegments () );
  getLengths ( temp );
  return temp.sum ();
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline istream& Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::read ( istream& is ) {
  return implementation ().read ( is );
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline Implementation& Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::implementation () {
  return static_cast <Implementation&> ( *this );
}

template <class Implementation, class IteratorType, class ConstIteratorType, class SegmentType, class ConstSegmentType, class DataType>
inline const Implementation& Particle<Implementation, IteratorType, ConstIteratorType, SegmentType, ConstSegmentType, DataType>::implementation () const {
  return static_cast <const Implementation&> ( *this );
}

// Class ReferenceSegmentIterator<ParticleType,SegmentType>

template <class ParticleType, class SegmentType>
inline ReferenceSegmentIterator<ParticleType, SegmentType>::ReferenceSegmentIterator ()
    : aol::Iterator<SegmentType> (), _particle ( 0 ), _index (), _segment ( *this->_particle, this->_index ) {}

template <class ParticleType, class SegmentType>
inline ReferenceSegmentIterator<ParticleType, SegmentType>::ReferenceSegmentIterator ( ParticleType& particle, bool end )
    : aol::Iterator<SegmentType> (), _particle ( &particle ), _index ( end ? particle.endIndex () : particle.beginIndex () ), _segment ( *this->_particle, this->_index ) {}

template <class ParticleType, class SegmentType>
inline ReferenceSegmentIterator<ParticleType, SegmentType>::ReferenceSegmentIterator ( const ReferenceSegmentIterator& it )
: aol::Iterator<SegmentType> (), _particle ( it._particle ), _index ( it._index ), _segment ( *this->_particle, this->_index ) {}

template <class ParticleType, class SegmentType>
inline ReferenceSegmentIterator<ParticleType, SegmentType>&
ReferenceSegmentIterator<ParticleType, SegmentType>::operator = ( const ReferenceSegmentIterator& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  aol::Iterator<SegmentType>::operator = ( it );

  this->_particle = it._particle;
  this->_index = it._index;
  this->_segment = it._segment;

  return *this;
}

template <class ParticleType, class SegmentType>
inline ReferenceSegmentIterator<ParticleType, SegmentType>::~ReferenceSegmentIterator () {}

template <class ParticleType, class SegmentType>
inline bool ReferenceSegmentIterator<ParticleType, SegmentType>::operator == ( const ReferenceSegmentIterator& it ) const {
  return this->_segment == it._segment;
}

template <class ParticleType, class SegmentType>
inline bool ReferenceSegmentIterator<ParticleType, SegmentType>::operator != ( const ReferenceSegmentIterator& it ) const {
  return ! ( it == *this );
}

template <class ParticleType, class SegmentType>
inline ReferenceSegmentIterator<ParticleType, SegmentType>&
ReferenceSegmentIterator<ParticleType, SegmentType>::operator ++ () {
  ++this->_index;
  this->updateSegment ();
  return *this;
}

template <class ParticleType, class SegmentType>
inline ReferenceSegmentIterator<ParticleType, SegmentType>&
ReferenceSegmentIterator<ParticleType, SegmentType>::operator -- () {
  --this->_index;
  this->updateSegment ();
  return *this;
}

template <class ParticleType, class SegmentType>
inline void ReferenceSegmentIterator<ParticleType, SegmentType>::updateSegment () {
  this->_segment = SegmentType ( *this->_particle, this->_index );
}
}

#endif
