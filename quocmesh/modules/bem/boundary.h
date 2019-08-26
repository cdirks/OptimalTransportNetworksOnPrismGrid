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

#ifndef __BOUNDARY_H
#define __BOUNDARY_H

#include <aol.h>
#include <typedParser.h>
#include <bemesh.h>
#include <particle.h>

namespace bm {

template <class ParticleType> class Boundary;

/**********************************************************/
//! Iterator that enumerates segments of particles in a list
//! Abstract w.r. to const-ness
template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
class PipingIterator
      : public aol::Iterator<SegmentType> {

protected:

  //! Default constructor is not allowed
  PipingIterator ();

public:

  //! Copy constructor
  PipingIterator ( const PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>& it );

  //! Assignment operator
  PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>& operator =
  ( const PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>& it );

  //! Destructor
  ~PipingIterator ();

  //! Constructor for start and end sentinel
  PipingIterator ( BoundaryType& bnd, bool end = false );

  //! Comparison operator
  bool operator == ( const PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>& it ) const;

  //! Comparison operator
  bool operator != ( const PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>& it ) const;

  //! Prefix operator ++
  PipingIterator& operator ++ ();

  //! Prefix operator --
  PipingIterator& operator -- ();

  //! Number of segments in current particle
  int getCurrentNumberOfSegments () const;

  //! Segment index in current particle
  int getCurrentIndex () const;

protected:

  //! Updates current segment after iterator changes
  void updateSegment ( bool end = false );

  ParticleIteratorType _particleBegin;
  ParticleIteratorType _particleEnd;
  ParticleIteratorType _particleCurrent;

  SegmentIteratorType _segmentBegin;
  SegmentIteratorType _segmentEnd;
  SegmentIteratorType _segmentCurrent;

  int _index;
};

/***********************************/
//! Segment iterator for all segments
template <class ParticleType>
class AllSegmentIterator
      : public PipingIterator < Boundary<ParticleType>, typename Boundary<ParticleType>::iterator,
      typename ParticleType::SegmentIteratorType, typename ParticleType::SegmentType > {

protected:

  //! Default constructor is not allowed
  AllSegmentIterator ();

public:

  //! Copy constructor
  AllSegmentIterator ( const AllSegmentIterator<ParticleType>& it );

  //! Assignment operator
  AllSegmentIterator<ParticleType>& operator = ( const AllSegmentIterator<ParticleType>& it );

  //! Destructor
  ~AllSegmentIterator ();

  //! Constructor for start and end sentinel
  AllSegmentIterator ( Boundary<ParticleType>& bnd, bool end = false );

  //! Dereference operator
  typename ParticleType::SegmentType* operator -> ();

  //! Dereference operator
  typename ParticleType::SegmentType& operator * ();
};

/********************************************/
//! Constant segment iterator for all segments
template <class ParticleType>
class ConstAllSegmentIterator
      : public PipingIterator < const Boundary<ParticleType>, typename Boundary<ParticleType>::const_iterator,
      typename ParticleType::ConstSegmentIteratorType, typename ParticleType::ConstSegmentType > {

protected:

  //! Default constructor is not allowed
  ConstAllSegmentIterator ();

public:

  //! Copy constructor
  ConstAllSegmentIterator ( const ConstAllSegmentIterator<ParticleType>& it );

  //! Assignment operator
  ConstAllSegmentIterator<ParticleType>& operator = ( const ConstAllSegmentIterator<ParticleType>& it );

  //! Destructor
  ~ConstAllSegmentIterator ();

  //! Constructor for start and end sentinel
  ConstAllSegmentIterator ( const Boundary<ParticleType>& bnd, bool end = false );

  //! Dereference operator
  const typename ParticleType::ConstSegmentType* operator -> () const;

  //! Dereference operator
  const typename ParticleType::ConstSegmentType& operator * () const;
};

/********************************************************/
//! A list of all particle boundaries in a problem setting
template <class ParticleType>
class Boundary
      : public list<ParticleType> {

public:
  using list<ParticleType>::insert;

  //! Default constructor
  Boundary ();

  //! Copy constructor
  Boundary ( const Boundary<ParticleType>& bnd );

  //! Construct from startName or deformedName in parameter file
  Boundary ( const aol::TypedParameterParser& parameter, bool deformed = false );

  //! Assignment operator
  Boundary<ParticleType>& operator = ( const Boundary<ParticleType>& bnd );

  //! Destructor
  ~Boundary ();

  //! Start SegmentIterator
  AllSegmentIterator<ParticleType> beginSegment ();

  //! End sentinel
  AllSegmentIterator<ParticleType> endSegment ();

  //! Start constant SegmentIterator
  ConstAllSegmentIterator<ParticleType> beginSegment () const;

  //! Constant end sentinel
  ConstAllSegmentIterator<ParticleType> endSegment () const;

  //! Number of segments
  int getNumberOfSegments () const;

  //! Number of points
  int getNumberOfPoints () const;

  //! Lengths of all Segments
  void getLengths ( aol::Vector<typename ParticleType::DataType>& lengths, bool betweencenters = false, bool jumping = false ) const;

  //! Curvatures on all Segments
  void getCurvatures ( aol::Vector<typename ParticleType::DataType>& curvatures, bool onFaces = false ) const;

  //! Positions of all Segment Start Points
  void getPositions ( aol::MultiVector<typename ParticleType::DataType>& positions, typename ParticleType::DOFPOSITION locpos = ParticleType::START ) const;

  //! Normals on all segments
  void getNormals ( aol::MultiVector<typename ParticleType::DataType>& normals, bool onFaces = false ) const;

  //! Output function
  ostream& print ( ostream& os ) const;

  //! Input function
  istream& read ( istream& is );

  //! Move by computed velocity
  void move ( const aol::Vector<typename ParticleType::DataType>& velocity );

  //! Move by displacement
  void move ( const aol::MultiVector<typename ParticleType::DataType>& displacement, bool facecentered );

  //! Rotate counter-clockwise
  void rotate ( typename ParticleType::DataType angle );

  //! Check by which fraction of the velocity the boundary can be moved without removing particles
  typename ParticleType::DataType getMovementFactor ( const aol::Vector<typename ParticleType::DataType>& velocity ) const;

  //! Remove all particles that would vanish with this velocity
  int countVanishingParticles ( const aol::Vector<typename ParticleType::DataType>& velocity, typename ParticleType::DataType factor ) const;

  //! Remove all particles that would vanish with this velocity
  void removeVanishingParticles ( const aol::Vector<typename ParticleType::DataType>& velocity,
          Boundary<ParticleType>& removedparticles,
          aol::Vector<int>& mapcoarsetofull, aol::Vector<int>& mapfinetofull,
          typename ParticleType::DataType minfactor = 1 );

  //! Compute total volume of particles
  typename ParticleType::DataType getVolume () const;

  //! Compute total surface of particles
  typename ParticleType::DataType getSurface () const;

  //! Compute average radius of particles
  typename ParticleType::DataType getAverageRadius () const;

  //! Compute average aspect ratio of bounding box, weighted by volume
  typename ParticleType::DataType getAverageBBAR () const;

  //! Compute minimal segment length
  typename ParticleType::DataType getMinLength () const;

  //! Remove small Particles
  bool cleanup ( typename ParticleType::DataType minVolume = 1E-12 );

  //! Remove small Particles
  bool cleanup ( typename ParticleType::DataType minVolume, aol::Vector<int>& mapnewtoold );

  //! Unify colliding particles
  void unifyCollisions ();

  //! Delete every n-th points
  void coarsen ( typename ParticleType::IndexType delEach = 2 );

  //! Remove small Particles
  void reparametrize ( typename ParticleType::DataType minAngle, typename ParticleType::DataType maxAngle,
                       typename ParticleType::DataType hFactor, typename ParticleType::DataType minLength );

  //! Remove small Particles
  void reparametrize ( typename ParticleType::DataType hFactor, typename ParticleType::DataType minLength );

  //! Find whether a given point lies within the particle phase
  bool inside ( aol::Vec2<typename ParticleType::DataType> point ) const;

  //! Output format
  typedef enum {
    BINARY, PACKED, PRETTY, GNUPLOT_AS_ONE_CURVE, GNUPLOT_AS_SEPARATE_CURVES, BITMAP, FILLEDBITMAP, FILLEDPDF, POSTSCRIPT 
  }
  FormatType;

  //! Change output format
  void setFormat ( FormatType format ) const;

  //! Read old output format
  FormatType getFormat () const;

  //! Change resolution for gnuplot output
  void setResolution ( int resolution ) const;

  //! Read old resolution
  int getResolution () const;

  //! Modify plotting of points
  void setPlotPoints ( bool plotThem ) const;

  //! Are points plotted?
  bool getPlotPoints () const;

  //! Set rectangular plotting area
  void setPlotBox ( aol::Vec2<typename ParticleType::DataType> ll, aol::Vec2<typename ParticleType::DataType> ur ) const;

  //! Set square plotting area
  void setPlotBox ( typename ParticleType::DataType a, typename ParticleType::DataType b ) const;

  //! Get plotting area
  void getPlotBox ( aol::Vec2<typename ParticleType::DataType> &ll, aol::Vec2<typename ParticleType::DataType> &ur ) const;

  //! Set colors for plotting, color channels range from 0 to one, particle colors may repeat
  void setColors ( const aol::Vec3<double>& background, const std::vector<aol::Vec3<double> >& particles ) const;

  //! Add additional point to all particles
  void makeJump ();

  //! Insert a new particle.
  //! Particle is absorbed (possibly multiple times recursively) or ignored if it hits an existing one
  //! You can also prescribe a minimum distance between particles
  bool insert ( const ParticleType& particle, typename ParticleType::DataType buffer = 0, bool unifyCollisions = true );

  //! Remove all particles that exceed a given bounding box
  void clipBoundingBox ( aol::Vec2<typename ParticleType::DataType> ll, aol::Vec2<typename ParticleType::DataType> ur );

  //! Compute a pointwise deformation based on a deformed version of the boundary
  void computeDeformation ( const Boundary& deformed, aol::MultiVector<typename ParticleType::DataType>& deformation,
                            bool facecentered = false ) const;

private:

  //! Store current format type
  mutable FormatType _format;
  //! Store gnuplot resolution
  mutable int _resolution;
  //! Plot points also in graphical ouptut
  mutable bool _plotPoints;
  //! Lower left and upper right edges of plotting area
  mutable aol::Vec2<double> _ll;
  mutable aol::Vec2<double> _ur;
  //! Plotting colors
  mutable aol::Vec3<double> _bgcolor;
  mutable std::vector<aol::Vec3<double> > _partcolors;
};

// Class PipingIterator<BoundaryType,ParticleIteratorType,SegmentIteratorType,SegmentType>

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
inline PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::PipingIterator ()
    : aol::Iterator<SegmentType> (), _particleBegin (), _particleEnd (), _particleCurrent (),
    _segmentBegin (), _segmentEnd (), _segmentCurrent (), _index ( 0 ) {}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
inline PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::PipingIterator
( const PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>& it )
    : aol::Iterator<SegmentType> ( it ), _particleBegin ( it._particleBegin ), _particleEnd ( it._particleEnd ), _particleCurrent ( it._particleCurrent ),
    _segmentBegin ( it._segmentBegin ), _segmentEnd ( it._segmentEnd ), _segmentCurrent ( it._segmentCurrent ), _index ( it._index ) {}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
inline PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>&
PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::operator =
( const PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  aol::Iterator<SegmentType>::operator = ( it );

  _particleBegin = it._particleBegin;
  _particleEnd = it._particleEnd;
  _particleCurrent = it._particleCurrent;
  _segmentBegin = it._segmentBegin;
  _segmentEnd = it._segmentEnd;
  _segmentCurrent = it._segmentCurrent;
  _index = it._index;

  return *this;
}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
inline PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::~PipingIterator () {}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
inline PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::PipingIterator ( BoundaryType& bnd, bool end )
    : aol::Iterator<SegmentType> (),
    _particleBegin ( bnd.begin () ), _particleEnd ( bnd.end () ), _particleCurrent ( end ? _particleEnd : _particleBegin ),
    _segmentBegin ( _particleBegin->beginSegment () ), _segmentEnd ( _particleBegin->endSegment () ), _segmentCurrent ( _segmentBegin ), _index ( 0 ) {}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
inline bool PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::operator ==
( const PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>& it ) const {
  return ( _particleCurrent == it._particleCurrent && _segmentCurrent == it._segmentCurrent ) ||
         ( _particleCurrent == _particleEnd && it._particleCurrent == it._particleEnd );
}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
inline bool PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::operator !=
( const PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>& it ) const {
  return ! ( *this == it );
}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
inline PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>&
PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::operator ++ () {
  ++_segmentCurrent; ++_index;
  if ( _segmentCurrent == _segmentEnd ) {
    ++_particleCurrent;
    updateSegment ();
  }

  return *this;
}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
inline PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>&
PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::operator -- () {
  if ( _segmentCurrent == _segmentBegin ) {
    --_particleCurrent;
    updateSegment ( true );
  }
  --_segmentCurrent; --_index;

  return *this;
}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
inline void PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::updateSegment ( bool end ) {
  if ( _particleCurrent == _particleEnd ) return;

  _segmentBegin = _particleCurrent->beginSegment ();
  _segmentEnd = _particleCurrent->endSegment ();

  _segmentCurrent = end ? _segmentEnd : _segmentBegin;
  _index = end ? _particleCurrent->getNumberOfSegments () : 0;
}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
int PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::getCurrentNumberOfSegments () const {
  return _particleCurrent->getNumberOfSegments ();
}

template <class BoundaryType, class ParticleIteratorType, class SegmentIteratorType, class SegmentType>
int PipingIterator<BoundaryType, ParticleIteratorType, SegmentIteratorType, SegmentType>::getCurrentIndex () const {
  return _index;
}

// Class AllSegmentIterator<ParticleType>

template <class ParticleType>
inline AllSegmentIterator<ParticleType>::AllSegmentIterator ()
    : PipingIterator < Boundary<ParticleType>, typename Boundary<ParticleType>::iterator,
    typename ParticleType::SegmentIteratorType, typename ParticleType::SegmentType > () {}

template <class ParticleType>
inline AllSegmentIterator<ParticleType>::AllSegmentIterator ( const AllSegmentIterator<ParticleType>& it )
    : PipingIterator < Boundary<ParticleType>, typename Boundary<ParticleType>::iterator,
    typename ParticleType::SegmentIteratorType, typename ParticleType::SegmentType > ( it ) {}

template <class ParticleType>
inline AllSegmentIterator<ParticleType>& AllSegmentIterator<ParticleType>::operator = ( const AllSegmentIterator<ParticleType>& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  PipingIterator < Boundary<ParticleType>, typename Boundary<ParticleType>::iterator,
  typename ParticleType::SegmentIteratorType, typename ParticleType::SegmentType >::operator = ( it );

  return *this;
}

template <class ParticleType>
inline AllSegmentIterator<ParticleType>::~AllSegmentIterator () {}

template <class ParticleType>
inline AllSegmentIterator<ParticleType>::AllSegmentIterator ( Boundary<ParticleType>& bnd, bool end )
    : PipingIterator < Boundary<ParticleType>, typename Boundary<ParticleType>::iterator,
    typename ParticleType::SegmentIteratorType, typename ParticleType::SegmentType > ( bnd, end ) {}

template <class ParticleType>
inline typename ParticleType::SegmentType* AllSegmentIterator<ParticleType>::operator -> () {
  return & ( *this->_segmentCurrent );
}

template <class ParticleType>
inline typename ParticleType::SegmentType& AllSegmentIterator<ParticleType>::operator * () {
  return *this->_segmentCurrent;
}

// Class ConstAllSegmentIterator<ParticleType>

template <class ParticleType>
inline ConstAllSegmentIterator<ParticleType>::ConstAllSegmentIterator ()
    : PipingIterator < const Boundary<ParticleType>, typename Boundary<ParticleType>::const_iterator,
    typename ParticleType::ConstSegmentIteratorType, typename ParticleType::ConstSegmentType > () {}

template <class ParticleType>
inline ConstAllSegmentIterator<ParticleType>::ConstAllSegmentIterator ( const ConstAllSegmentIterator<ParticleType>& it )
    : PipingIterator < const Boundary<ParticleType>, typename Boundary<ParticleType>::const_iterator,
    typename ParticleType::ConstSegmentIteratorType, typename ParticleType::ConstSegmentType > ( it ) {}

template <class ParticleType>
inline ConstAllSegmentIterator<ParticleType>& ConstAllSegmentIterator<ParticleType>::operator = ( const ConstAllSegmentIterator<ParticleType>& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  PipingIterator < const Boundary<ParticleType>, typename Boundary<ParticleType>::const_iterator,
  typename ParticleType::ConstSegmentIteratorType, typename ParticleType::ConstSegmentType >::operator = ( it );

  return *this;
}

template <class ParticleType>
inline ConstAllSegmentIterator<ParticleType>::~ConstAllSegmentIterator () {}

template <class ParticleType>
inline ConstAllSegmentIterator<ParticleType>::ConstAllSegmentIterator ( const Boundary<ParticleType>& bnd, bool end )
    : PipingIterator < const Boundary<ParticleType>, typename Boundary<ParticleType>::const_iterator,
    typename ParticleType::ConstSegmentIteratorType, typename ParticleType::ConstSegmentType > ( bnd, end ) {}

template <class ParticleType>
inline const typename ParticleType::ConstSegmentType* ConstAllSegmentIterator<ParticleType>::operator -> () const {
  return & ( *this->_segmentCurrent );
}

template <class ParticleType>
inline const typename ParticleType::ConstSegmentType& ConstAllSegmentIterator<ParticleType>::operator * () const {
  return *this->_segmentCurrent;
}

// Class Boundary<ParticleType>

template <class ParticleType>
inline Boundary<ParticleType>::Boundary ()
  : list<ParticleType> (), _format ( PRETTY ), _resolution ( 1024 ), _plotPoints ( false ), _ll ( -0.1, -0.1 ), _ur ( 1.1, 1.1 ), _bgcolor (1, 0.94, 0.78), _partcolors () {
  _partcolors.push_back (aol::Vec3<double> (0,0,1));
 }

template <class ParticleType>
inline Boundary<ParticleType>::Boundary ( const Boundary<ParticleType>& bnd )
  : list<ParticleType> ( bnd ), _format ( bnd._format ), _resolution ( bnd._resolution ), _plotPoints ( bnd._plotPoints ), _ll ( bnd._ll ), _ur ( bnd._ur ),
  _bgcolor ( bnd._bgcolor ), _partcolors ( bnd._partcolors ) {}

template <class ParticleType>
inline Boundary<ParticleType>::Boundary ( const aol::TypedParameterParser& parameter, bool deformed )
  : list<ParticleType> (), _format ( PACKED ), _resolution ( 1024 ), _plotPoints ( false ), _ll ( -0.1, -0.1 ), _ur ( 1.1, 1.1 ), _bgcolor (1, 1, 0), _partcolors () {
  std::string filename;
  if ( deformed ) parameter.get ( "deformedName", filename );
  else parameter.get ( "startName", filename );
  filename += ".dat";
  aol::ipfstream file ( filename.c_str () );
  read ( file );

  bool binaryData;
  parameter.get ( "binaryData", binaryData );
  if ( binaryData ) setFormat ( BINARY );
  _partcolors.push_back (aol::Vec3<double> (0,0,1));
}

template <class ParticleType>
inline Boundary<ParticleType>& Boundary<ParticleType>::operator = ( const Boundary<ParticleType>& bnd ) {
  // Beware of self-assignment
  if ( this == &bnd ) return *this;

  list<ParticleType>::operator = ( bnd );
  _format = bnd._format;
  _resolution = bnd._resolution;
  _plotPoints = bnd._plotPoints;
  _ll = bnd._ll;
  _ur = bnd._ur;
  _bgcolor = bnd._bgcolor;
  _partcolors = bnd._partcolors;

  return *this;
}

template <class ParticleType>
inline Boundary<ParticleType>::~Boundary () {}

template <class ParticleType>
inline AllSegmentIterator<ParticleType> Boundary<ParticleType>::beginSegment () {
  AllSegmentIterator<ParticleType> it ( *this );
  return it;
}

template <class ParticleType>
inline AllSegmentIterator<ParticleType> Boundary<ParticleType>::endSegment () {
  AllSegmentIterator<ParticleType> it ( *this, true );
  return it;
}

template <class ParticleType>
inline ConstAllSegmentIterator<ParticleType> Boundary<ParticleType>::beginSegment () const {
  ConstAllSegmentIterator<ParticleType> it ( *this );
  return it;
}

template <class ParticleType>
inline ConstAllSegmentIterator<ParticleType> Boundary<ParticleType>::endSegment () const {
  ConstAllSegmentIterator<ParticleType> it ( *this, true );
  return it;
}

template <class ParticleType>
inline int Boundary<ParticleType>::getNumberOfSegments () const {
  int size = 0;
  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it )
    size += it->getNumberOfSegments ();
  return size;
}

template <class ParticleType>
inline int Boundary<ParticleType>::getNumberOfPoints () const {
  int size = 0;
  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it )
    size += it->getNumberOfPoints ();
  return size;
}

template <class ParticleType>
inline void Boundary<ParticleType>::makeJump () {
  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ++it )
    it->makeJump ();
}

template <class ParticleType>
void Boundary<ParticleType>::getLengths ( aol::Vector<typename ParticleType::DataType>& lengths, bool betweencenters, bool jumping ) const {
  int pos = 0;

  if ( lengths.size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::getLengths, vector length does not equal number of segments", __FILE__, __LINE__ );

  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it ) {

    int size = it->getNumberOfSegments ();
    aol::Vector<typename ParticleType::DataType> temp ( size );
    it->getLengths ( temp, betweencenters, jumping );

    lengths.setBlock ( pos, temp );
    pos += size;
  }
}

template <class ParticleType>
void Boundary<ParticleType>::getCurvatures ( aol::Vector<typename ParticleType::DataType>& curvatures, bool onFaces ) const {
  int pos = 0;

  if ( curvatures.size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::getCurvatures, vector length does not equal number of segments", __FILE__, __LINE__ );

  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it ) {

    int size = it->getNumberOfSegments ();
    aol::Vector<typename ParticleType::DataType> temp ( size );
    it->getCurvatures ( temp, onFaces );

    curvatures.setBlock ( pos, temp );
    pos += size;
  }
}

template <class ParticleType>
void Boundary<ParticleType>::getPositions ( aol::MultiVector<typename ParticleType::DataType>& positions, typename ParticleType::DOFPOSITION locpos ) const {
  int pos = 0;

  if ( positions.numComponents () != 2 ||
       positions [0].size () != getNumberOfSegments () || positions [1].size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::getPositions, vector length does not equal number of segments", __FILE__, __LINE__ );

  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it ) {

    int size = it->getNumberOfSegments ();
    aol::MultiVector<typename ParticleType::DataType> temp ( 2, size );
    it->getPositions ( temp, locpos );

    for ( int i = 0; i < 2; ++i )
      positions [i].setBlock ( pos, temp [i] );
    pos += size;
  }
}

template <class ParticleType>
void Boundary<ParticleType>::getNormals ( aol::MultiVector<typename ParticleType::DataType>& positions, bool onFaces ) const {
  int pos = 0;

  if ( positions.numComponents () != 2 ||
       positions [0].size () != getNumberOfSegments () || positions [1].size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::getNormals, vector length does not equal number of segments", __FILE__, __LINE__ );

  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it ) {

    int size = it->getNumberOfSegments ();
    aol::MultiVector<typename ParticleType::DataType> temp ( 2, size );
    it->getNormals ( temp, onFaces );

    for ( int i = 0; i < 2; ++i )
      positions [i].setBlock ( pos, temp [i] );
    pos += size;
  }
}

namespace {
  string gnuplotColor ( aol::Vec3<double>& col ) {
    std::stringstream result ( "x" );
    for ( int i = 0; i < 3; ++i ) {
      double v = col [i] * 256;
      if ( v < 0 ) v = 0;
      if ( v > 255 ) v = 255;
      int c = static_cast<unsigned char> ( v );
      result.flags ( ios::right | ios::hex );
      result.width ( 2 );
      result.fill ( '0' );
      result << c;
    }
    return result.str ();
  }
}

template <class ParticleType>
inline ostream& Boundary<ParticleType>::print ( ostream& os ) const {

  bool plotbnd = false;
  string intx = string ( "[" ) + aol::detailedFormat ( _ll[ 0 ] ) + ":" + aol::detailedFormat ( _ur[ 0 ] ) + "]";
  string inty = string ( "[" ) + aol::detailedFormat ( _ll[ 1 ] ) + ":" + aol::detailedFormat ( _ur[ 1 ] ) + "]";
  double ratio = ( _ur[1] - _ll[1] ) / ( _ur[0] - _ll[0] );

  switch ( getFormat () ) {

  case BINARY:
    os << "@begin binary boundary" << endl;
    for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it ) {
      it->setFormat ( ParticleType::BINARY );
      os << *it;
    }
    os << "@end binary boundary" << endl << endl;
    break;

  case PACKED:
    os << "@begin packed boundary" << endl;
    for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it ) {
      it->setFormat ( ParticleType::PACKED );
      os << *it;
    }
    os << "@end packed boundary" << endl << endl;
    break;

  case PRETTY:
    os << "/************************ Boundary ***********************\\" << endl;
    os << "|                                                         |" << endl;
    for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it ) {
      it->setFormat ( ParticleType::PRETTY );
      os << *it;
    }
    os << "\\*********************************************************/" << endl << endl;
    break;

  case GNUPLOT_AS_ONE_CURVE:
    for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it ) {
      it->setFormat ( ParticleType::GNUPLOT );
      os << *it << endl;
    }
    break;

  case GNUPLOT_AS_SEPARATE_CURVES:
    for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it ) {
      it->setFormat ( ParticleType::GNUPLOT );
      os << *it << "e" << endl;
    }
    break;

  case BITMAP:
    os << "set terminal pbm color" << endl;
    os << "set size " << static_cast<double> ( _resolution ) / 640 << ", " << static_cast<double> ( _resolution ) / 480 << endl;
    if ( !plotbnd ) os << "set noborder" << endl;
    os << "set noxtics" << endl;
    os << "set noytics" << endl;
    os << "plot " << intx << " " << inty << " '-' notitle with " << ( _plotPoints ? "linespoints" : "lines" ) << endl;
    setFormat ( GNUPLOT_AS_ONE_CURVE );
    os << *this;
    setFormat ( BITMAP );
    os << "e" << endl;
    break;

  case FILLEDBITMAP:
    os << "set terminal png size " << _resolution << ", " << _resolution << " enhanced x" << gnuplotColor ( _bgcolor ) << " x000000 x00000";
    for ( unsigned int i = 0; i < this->size (); ++i ) {
      os << " x" << gnuplotColor (_partcolors [i % _partcolors.size ()]);
    }
    os << endl;
    os << "set noborder" << endl;
    os << "set noxtics" << endl;
    os << "set noytics" << endl;
    os << "plot " << intx << " " << inty << " '-' notitle with filledcurves";
    for ( unsigned int i = 1; i < this->size (); ++i ) { os << ", '-' notitle with filledcurves"; }
    os << endl;
    setFormat ( GNUPLOT_AS_SEPARATE_CURVES );
    os << *this;
    setFormat ( FILLEDBITMAP );
    break;

  case FILLEDPDF:
    os << "set terminal pdfcairo color rounded size 6cm,6cm" << endl;
    os << "set noborder" << endl;
    os << "set noxtics" << endl;
    os << "set noytics" << endl;
    os << "plot " << intx << " " << inty << " '-' notitle with filledcurves linecolor rgb \"#" << gnuplotColor (_partcolors [0]) << "\"";
    for ( unsigned int i = 1; i < this->size (); ++i ) { os << ", '-' notitle with filledcurves linecolor rgb \"#" << gnuplotColor (_partcolors [i % _partcolors.size ()]) << "\""; }
    os << endl;
    setFormat ( GNUPLOT_AS_SEPARATE_CURVES );
    os << *this;
    setFormat ( FILLEDPDF );
    break;

  case POSTSCRIPT:
    os << "set terminal postscript enhanced color eps" << endl;
    os << "set size ratio " << ratio << " 1, " << 10.0/7.0 * ratio << endl; // Correct scaling
    if ( !plotbnd ) os << "set noborder" << endl;
    os << "set noxtics" << endl;
    os << "set noytics" << endl;
    os << "plot " << intx << " " << inty << " '-' notitle with " << ( _plotPoints ? "linespoints" : "lines" ) << endl;
    setFormat ( GNUPLOT_AS_ONE_CURVE );
    os << *this;
    setFormat ( POSTSCRIPT );
    os << "e" << endl;
    break;

  default:
    throw aol::Exception ( "bem::Boundary::print: unknown FormatType", __FILE__, __LINE__ );
  }

  return os;
}

template <class ParticleType>
istream& Boundary<ParticleType>::read ( istream& is ) {
  // Delete old particles
  this->clear ();

  string temp;

  // Only packed input is possible
  getline ( is, temp );
  if ( temp != "@begin packed boundary" && temp != "@begin binary boundary" )
    throw ( aol::FileFormatException ( "Boundary::read, can only read packed and binary format, found line [" + temp + "]", __FILE__, __LINE__ ) );

  getline ( is, temp );
  while ( temp != "@end packed boundary" && temp != "@end binary boundary" ) {
    is.seekg ( - static_cast<int> ( temp.size () ) - 1, ios_base::cur ); // Go back a line, static_cast necessary for correct negation
    ParticleType particle;
    is >> particle;
    this->push_back ( particle );
    getline ( is, temp );
  }

  return is;
}


template <class ParticleType>
inline ostream& operator << ( ostream& os, const Boundary<ParticleType>& bnd ) {
  return bnd.print ( os );
}

template <class ParticleType>
inline istream& operator >> ( istream& is, Boundary<ParticleType>& bnd ) {
  return bnd.read ( is );
}

template <class ParticleType>
void Boundary<ParticleType>::move ( const aol::Vector<typename ParticleType::DataType>& velocity ) {
  int pos = 0;
  if ( velocity.size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::move, vector length does not equal number of segments", __FILE__, __LINE__ );

  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ++it ) {

    int size = it->getNumberOfSegments ();
    aol::Vector<typename ParticleType::DataType> temp ( size );
    velocity.getBlock ( pos, temp );
    pos += size;
    it->move ( temp );
  }
}

template <class ParticleType>
void Boundary<ParticleType>::move ( const aol::MultiVector<typename ParticleType::DataType>& displacement, bool facecentered ) {
  int pos = 0;
  if ( displacement.numComponents () != 2 ||
       displacement [0].size () != getNumberOfSegments () || displacement [1].size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::move, vector length does not equal number of segments", __FILE__, __LINE__ );

  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ++it ) {

    int size = it->getNumberOfSegments ();
    aol::MultiVector<typename ParticleType::DataType> temp ( 2, size );
    for ( int i = 0; i < 2; ++i )
      displacement [i].getBlock ( pos, temp [i] );
    pos += size;
    it->move ( temp, facecentered );
  }
}

template <class ParticleType>
typename ParticleType::DataType Boundary<ParticleType>::getMovementFactor ( const aol::Vector<typename ParticleType::DataType>& velocity ) const {
  int i = 0;
  typename ParticleType::DataType factor = 1000; // Do not overrelaxate larger than this

  if ( velocity.size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::move, vector length does not equal number of segments", __FILE__, __LINE__ );

  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); i += it->size (), ++it ) {
    aol::Vector<typename ParticleType::DataType> localVelocity ( it->size () );
    velocity.getBlock ( i, localVelocity );
    typename ParticleType::DataType newFactor = it->getMovementFactor ( localVelocity );
    factor = std::min ( factor, newFactor );
  }

  return factor;
}

template <class ParticleType>
int Boundary<ParticleType>::countVanishingParticles ( const aol::Vector<typename ParticleType::DataType>& velocity,
                  typename ParticleType::DataType factor ) const {

  int i = 0, n = 0;

  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); i += it->size (), ++it ) {
    aol::Vector<typename ParticleType::DataType> localVelocity ( it->size () );
    velocity.getBlock ( i, localVelocity );
    if ( it->getMovementFactor ( localVelocity ) <= factor ) ++n;
  }

  return n;
}

template <class ParticleType>
void Boundary<ParticleType>::removeVanishingParticles ( const aol::Vector<typename ParticleType::DataType>& velocity,
              Boundary<ParticleType>& removedparticles,
              aol::Vector<int>& mapcoarsetofull, aol::Vector<int>& mapfinetofull,
              typename ParticleType::DataType minfactor ) {
  int i = 0, cj = 0, fj = 0;

  if ( velocity.size () != getNumberOfSegments () )
    throw aol::DimensionMismatchException ( "Particle::removeVanishingParticles, vector length does not equal number of segments", __FILE__, __LINE__ );

  mapcoarsetofull.resize ( velocity.size () );
  mapcoarsetofull.setAll ( -1 );

  mapfinetofull.resize ( velocity.size () );
  mapfinetofull.setAll ( -1 );

  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ) {
    aol::Vector<typename ParticleType::DataType> localVelocity ( it->size () );
    velocity.getBlock ( i, localVelocity );
    if ( it->getMovementFactor ( localVelocity ) <= minfactor ) {

      removedparticles.push_back ( *it );
      for ( int k = 0; k < it->size (); ++k )
  mapfinetofull [ fj + k ] = i + k;
      fj += it->size ();
      i += it->size ();
      it = this->erase ( it ); // includes increment
    }
    else {
      for ( int k = 0; k < it->size (); ++k )
  mapcoarsetofull [ cj + k ] = i + k;
      cj += it->size ();
      i += it->size ();
      ++it;
    }
  }

  mapcoarsetofull.resize ( getNumberOfSegments () );
  mapfinetofull.resize ( removedparticles.getNumberOfSegments () );
}

template <class ParticleType>
typename ParticleType::DataType Boundary<ParticleType>::getVolume () const {
  typename ParticleType::DataType volume ( 0 );
  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it )
    volume += it->getVolume ();

  return volume;
}

template <class ParticleType>
typename ParticleType::DataType Boundary<ParticleType>::getSurface () const {
  typename ParticleType::DataType surface = 0;
  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it )
    surface += it->getSurface ();

  return surface;
}

template <class ParticleType>
typename ParticleType::DataType Boundary<ParticleType>::getAverageRadius () const {
  typename ParticleType::DataType radius = 0;
  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it )
    radius += it->getAverageRadius ();

  return radius / this->size ();
}

template <class ParticleType>
typename ParticleType::DataType Boundary<ParticleType>::getAverageBBAR () const {
  typename ParticleType::DataType ar = 0, vol = 0;
  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it ) {
    aol::Vec2<typename ParticleType::DataType> ll, ur;
    it->getBoundingBox ( ll, ur );
    typename ParticleType::DataType w = ur [0] - ll [0], h = ur [1] - ll [1], v = it->getVolume ();
    ar += (w < h ? h / w : w / h) * v;
    vol += v;
  }
  return ar / vol;
}

template <class ParticleType>
typename ParticleType::DataType Boundary<ParticleType>::getMinLength () const {
  typename ParticleType::DataType minlen = std::numeric_limits<typename ParticleType::DataType>::infinity ();
  for ( ConstAllSegmentIterator<ParticleType> it = this->beginSegment (); it != this->endSegment (); ++it )
    minlen = std::min ( minlen, it->getLength () );

  return minlen;
}

template <class ParticleType>
bool Boundary<ParticleType>::cleanup ( typename ParticleType::DataType minVolume ) {
  bool dunnit = false;
  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ++it )
    if ( it->getVolume () < minVolume ) {
      typename std::list<ParticleType>::iterator target = it;
      ++it;
      this->erase ( target );
      --it;
      if ( aol::debugging::particle ) cerr << "Particle small, removed." << endl;
      dunnit = true;
    }
  return dunnit;
}

template <class ParticleType>
bool Boundary<ParticleType>::cleanup ( typename ParticleType::DataType minVolume, aol::Vector<int>& mapnewtoold ) {
  bool dunnit = false;
  int oldi = 0, newi = 0;
  mapnewtoold.resize ( getNumberOfSegments () );
  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ++it ) {
    if ( it->getVolume () < minVolume ) {
      oldi += it->size ();
      typename std::list<ParticleType>::iterator target = it;
      ++it;
      this->erase ( target );
      --it;
      if ( aol::debugging::particle ) cerr << "Particle small, removed." << endl;
      dunnit = true;
    } else {
      for ( int j = 0; j < it->size (); ++j, ++newi, ++oldi ) mapnewtoold [ newi ] = oldi;
    }
  }
  mapnewtoold.resize ( newi );
  return dunnit;
}

template <class ParticleType>
  void Boundary<ParticleType>::unifyCollisions () {

  cerr << "Testing for collisions ... ";

  Boundary bnd = *this;
  this->clear ();
  for ( typename std::list<ParticleType>::iterator it = bnd.begin (); it != bnd.end (); ++it )
    insert ( *it, 0, true );

  cerr << "done." << endl;
}

template <class ParticleType>
void Boundary<ParticleType>::rotate ( typename ParticleType::DataType angle ) {
  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ++it )
    it->rotate ( angle );
}

template <class ParticleType>
void Boundary<ParticleType>::coarsen ( typename ParticleType::IndexType delEach ) {
  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ++it )
    it->coarsen ( delEach );
}

template <class ParticleType>
void Boundary<ParticleType>::reparametrize ( typename ParticleType::DataType minAngle, typename ParticleType::DataType maxAngle,
                                             typename ParticleType::DataType hFactor, typename ParticleType::DataType minLength ) {
  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ++it )
    it->reparametrize ( minAngle, maxAngle, hFactor, minLength );
}

template <class ParticleType>
void Boundary<ParticleType>::reparametrize ( typename ParticleType::DataType hFactor, typename ParticleType::DataType minLength ) {
  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ++it )
    it->reparametrize ( hFactor, minLength );
}

template <class ParticleType>
bool Boundary<ParticleType>::inside ( aol::Vec2<typename ParticleType::DataType> point ) const {
  for ( typename std::list<ParticleType>::const_iterator it = this->begin (); it != this->end (); ++it )
    if ( it->inside ( point ) )
      return true;

  return false;
}

template <class ParticleType>
bool Boundary<ParticleType>::insert ( const ParticleType& particle, typename ParticleType::DataType buffer, bool unifyCollisions ) {
  ParticleType part ( particle );
  typename std::list<ParticleType>::iterator it = this->begin ();
  while ( it != this->end () ) {
    if ( it->collides ( particle, buffer ) ) {

      // Collision found
      if ( unifyCollisions ) {
        part.absorb ( *it );
        typename std::list<ParticleType>::iterator target = it;
        this->erase ( target );
        if ( aol::debugging::particle ) cerr << "Particle collides, united." << endl;

        // Try again with larger particle
        it = this->begin ();
      } else {
        // Ignore
        if ( aol::debugging::particle ) cerr << "Particle collides, ignored." << endl;
        return false;
      }
    } else {
      ++it;
    }
  }

  this->push_back ( part );
  return true;
}

template <class ParticleType>
void Boundary<ParticleType>::clipBoundingBox ( aol::Vec2<typename ParticleType::DataType> ll, aol::Vec2<typename ParticleType::DataType> ur ) {
  for ( typename std::list<ParticleType>::iterator it = this->begin (); it != this->end (); ++it )
    if ( it->exceeds ( ll, ur ) ) {
      typename std::list<ParticleType>::iterator target = it;
      ++it;
      this->erase ( target );
      --it;
      if ( aol::debugging::particle ) cerr << "Particle outside, removed." << endl;
    }
}

template <class ParticleType>
void Boundary<ParticleType>::computeDeformation ( const Boundary& deformed,
                                                  aol::MultiVector<typename ParticleType::DataType>& deformation,
                                                  bool facecentered ) const {
  int i = 0;
  for ( ConstAllSegmentIterator<ParticleType> uit = beginSegment(), uend = endSegment (),
        dit = deformed.beginSegment (), dend = deformed.endSegment (); uit != uend; ++uit, ++dit, ++i ) {
    aol::Vec2<typename ParticleType::DataType> upt = facecentered ? 0.5 * ( uit->getStart () + uit->getEnd () ) : uit->getStart (),
                                                     dpt = facecentered ? 0.5 * ( dit->getStart () + dit->getEnd () ) : dit->getStart (), def = dpt - upt;

    deformation [0][i] = def [0];
    deformation [1][i] = def [1];
  }
}

template <class ParticleType>
inline void Boundary<ParticleType>::setFormat ( FormatType format ) const {
  _format = format;
}

template <class ParticleType>
inline typename Boundary<ParticleType>::FormatType Boundary<ParticleType>::getFormat () const {
  return _format;
}

template <class ParticleType>
inline void Boundary<ParticleType>::setResolution ( int resolution ) const {
  _resolution = resolution;
}

template <class ParticleType>
inline int Boundary<ParticleType>::getResolution () const {
  return _resolution;
}

template <class ParticleType>
inline void Boundary<ParticleType>::setPlotPoints ( bool plotThem ) const {
  _plotPoints = plotThem;
}

template <class ParticleType>
inline bool Boundary<ParticleType>::getPlotPoints () const {
  return _plotPoints;
}

template <class ParticleType>
inline void Boundary<ParticleType>::setPlotBox ( typename ParticleType::DataType a, typename ParticleType::DataType b ) const {
  _ll[ 0 ] = a; _ll[ 1 ] = a;
  _ur[ 0 ] = b; _ur[ 1 ] = b;
}

template <class ParticleType>
inline void Boundary<ParticleType>::setPlotBox ( aol::Vec2<typename ParticleType::DataType> ll, aol::Vec2<typename ParticleType::DataType> ur ) const {
  _ll = ll;
  _ur = ur;
}

template <class ParticleType>
inline void Boundary<ParticleType>::getPlotBox ( aol::Vec2<typename ParticleType::DataType> &ll, aol::Vec2<typename ParticleType::DataType> &ur ) const {
  ll = _ll;
  ur = _ur;
}

template <class ParticleType>
inline void Boundary<ParticleType>::setColors ( const aol::Vec3<double>& background, const std::vector<aol::Vec3<double> >& particles ) const {
  _bgcolor = background;
  _partcolors = particles;
}
} // namespace bm

#endif
