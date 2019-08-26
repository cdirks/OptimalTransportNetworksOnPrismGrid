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

#ifndef __NARROWBANDGRID_H
#define __NARROWBANDGRID_H

#include <quoc.h>
#include <scalarArray.h>
#include <eikonalNA3d.h>
#include <qmException.h>
#include <configurators.h>
#include <hashSet.h>

namespace nb {

//-----------------------------------------------------------------------------------------------------------------------------

template <typename _FullGridConfType, typename _InitType>
class FullGridConfiguratorHull;

//-----------------------------------------------------------------------------------------------------------------------------

template <typename _GridType>
class ElementIterator {
  typedef typename aol::HashSet<typename _GridType::ElementType>::const_iterator Base;
public:
  typedef _GridType                                 GridType;
  typedef typename GridType::ElementType            IteratedType;
  typedef GridType                                  BeginType;
  typedef qc::EndElement                            EndType;
  typedef ElementIterator<GridType>                 Self;

  ElementIterator ( const BeginType & Grid ) {
    _cur = Grid.getHashSet().begin();
    _endElement = Grid.getHashSet().end();
  }

  //! Initialize iterator from grid
  const Self & operator= ( const BeginType & Grid ) {
    // set current element to first element
    _cur =  Grid.getHashSet().begin();

    // set termination element
    _endElement = Grid.getHashSet().end();
    return *this;
  }

  //! Initialize iterator from configurator
  template <typename _FullGridConfType>
  const Self & operator= ( const FullGridConfiguratorHull<_FullGridConfType, GridType> & Configurator ) {
    operator= ( Configurator.getInitializer() );
  }

  //! test for finish
  bool operator== ( const EndType & ) const {
    return _cur == _endElement;
  }

  //! test for continuing
  bool operator!= ( const EndType & ) const {
    return _cur != _endElement;
  }

  const Self & operator ++ () {
    _cur++;
    return *this;
  }
  const IteratedType & operator* () {
    return _cur.operator*();
  }

  IteratedType const * operator-> () {
    return _cur.operator->();
  }

protected:
  Base _endElement;
  Base _cur;
};

//-----------------------------------------------------------------------------------------------------------------------------

template <typename _FullGridType>
class NarrowBandGridBase : protected _FullGridType {

public:
  typedef _FullGridType                             FullGridType;
  typedef FullGridType                              Base;
  typedef typename FullGridType::ElementType        ElementType;
  typedef aol::HashSet<ElementType>                 HashType;

  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

protected:
  HashType _hash;

public:
  explicit NarrowBandGridBase ( const FullGridType &grid )
  : _FullGridType ( grid )
  , _hash () {}

  void insert ( const ElementType & el ) {
    _hash.insert ( el );
  }

  bool exists ( const ElementType & el ) const {
    return ( _hash.find ( el ) != _hash.end() );
  }

  void clear( ) {
    _hash.clear( );
  }

  //! \todo: Name misleading, rename to getFullGridRef
  const FullGridType & getFullGrid () const {
    return static_cast<const FullGridType & >(*this);
  }

  //! \todo: Name misleading, rename to getHashSetRef
  virtual const HashType & getHashSet( ) const { return _hash; }

  using Base::H;
  using Base::getNumberOfNodes;
  using Base::getGridDepth;
  using Base::getNumX;
  using Base::getNumY;
  using Base::getNumZ;
  using Base::getSize;
  using Base::getDimOfWorld;
  using Base::getElementIndex;

  using Base::isAdaptive;
  using Base::checkForHangingNode;
};

//-----------------------------------------------------------------------------------------------------------------------------

template <typename _FullGridType, qc::Dimension Dim> class NarrowBandGrid {};

template <typename _FullGridType>
class NarrowBandGrid<_FullGridType, qc::QC_2D>
  : public NarrowBandGridBase<_FullGridType> {

public:
  typedef _FullGridType                             FullGridType;
  typedef NarrowBandGridBase<FullGridType>          Base;
  typedef NarrowBandGrid<FullGridType, qc::QC_2D>   Self;
  typedef typename FullGridType::CubicGridType      CubicGridType;
  typedef typename Base::ElementType                ElementType;

  typedef ElementIterator<Self>                     OldFullElementIterator;
  //typedef NodeIterator<Self, qc::QC_2D>             OldFullNodeIterator;
  class OldFullNodeIterator {};

  typedef typename OldFullElementIterator::BeginType   BeginIterType;
  typedef typename OldFullElementIterator::EndType     EndIterType;

  static const qc::Dimension DimOfWorld = qc::QC_2D;

  // we provide an empty type Face which is requested by the configurator
  class Face {};

  class Edge {
  public:
    Edge ( const qc::Element &el, int loc1, int loc2 )
        : _el ( el ) {
      _locNode[0] = loc1;
      _locNode[1] = loc2;
    }

    qc::Element _el;
    int _locNode[2];
  };

  typedef typename vector<Edge>::const_iterator eiterator;
  typedef typename vector<qc::CoordType>::const_iterator biterator;

protected:
  vector<Edge> _edges;
  vector<qc::CoordType> _boundaryNodes;

public:
  explicit NarrowBandGrid ( const FullGridType & grid )
  : NarrowBandGridBase<_FullGridType> ( grid ) {}

  const BeginIterType & begin() const {
    return *this;
  }

  EndIterType end() const {
    return EndIterType();
  }

  void extractEdges( );

  eiterator ebegin( ) const { return _edges.begin();  }
  eiterator eend( ) const { return _edges.end(); }

  biterator bbegin( ) const { return _boundaryNodes.begin(); }
  biterator bend( )  const { return _boundaryNodes.end(); }
};

//-----------------------------------------------------------------------------------------------------------------------------

template <typename _FullGridType>
class NarrowBandGrid<_FullGridType, qc::QC_3D>
  : public NarrowBandGridBase<_FullGridType> {

public:
  typedef _FullGridType                             FullGridType;
  typedef NarrowBandGridBase<FullGridType>          Base;
  typedef NarrowBandGrid<FullGridType, qc::QC_3D>   Self;
  typedef typename FullGridType::CubicGridType      CubicGridType;
  typedef typename Base::ElementType                ElementType;

  typedef ElementIterator<Self>                     OldFullElementIterator;
  //typedef NodeIterator<Self, qc::QC_3D>             OldFullNodeIterator;
  class OldFullNodeIterator {};

  typedef typename OldFullElementIterator::BeginType   BeginIterType;
  typedef typename OldFullElementIterator::EndType     EndIterType;

  static const qc::Dimension DimOfWorld = qc::QC_3D;

  class Face {
  public:
    Face ( const qc::Element & el, int loc1, int loc2, int loc3, int loc4 )
        : _el ( el ) {
      _locNode[0] = loc1;
      _locNode[1] = loc2;
      _locNode[2] = loc3;
      _locNode[3] = loc4;
    }

    qc::Element _el;
    int _locNode[4];
  };

  typedef typename vector<Face>::const_iterator eiterator;
  typedef typename vector<qc::CoordType>::const_iterator biterator;

protected:
  vector<Face>    _edges;
  vector<qc::CoordType> _boundaryNodes;

public:
  explicit NarrowBandGrid ( const FullGridType & grid )
  : NarrowBandGridBase<_FullGridType> ( grid ) {}

  const BeginIterType & begin() const {
    return *this;
  }

  EndIterType end() const {
    return EndIterType();
  }

  void extractEdges( );

  eiterator ebegin( ) const {
    return _edges.begin();
  }

  eiterator eend( ) const {
    return _edges.end();
  }

  biterator bbegin( ) const { return _boundaryNodes.begin(); }
  biterator bend( ) const { return _boundaryNodes.end(); }
};

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, qc::Dimension Dim> class NarrowBandMaskedGrid { };

//-----------------------------------------------------------------------------------------------------------------------------

template <typename _RealType>
class NarrowBandMaskedGrid<_RealType, qc::QC_3D>
  : public NarrowBandGrid <qc::GridDefinition, qc::QC_3D> {
public:
  typedef _RealType Real;
  typedef OldFullElementIterator iterator;
  enum ITER_MODE { INT, EXT };
protected:
  // use a pointer to the member in Willmore class here instead
  const qc::ScalarArray<_RealType, qc::QC_3D> &_mask_img;
  // const qc::ScalarArray<_RealType, qc::QC_3D> &_levelset;
  vector<Face> _edgesD;
  vector<qc::CoordType> _boundaryNodesD;
  NarrowBandGrid<qc::GridDefinition, qc::QC_3D>::HashType _hashInterior;

  vector<qc::CoordType> _outerBoundaryNodes;

  ITER_MODE _iterMode;

public:
  NarrowBandMaskedGrid ( const qc::GridDefinition &grid, const qc::ScalarArray<Real, qc::QC_3D> &mask_img )
    : NarrowBandGrid<qc::GridDefinition, qc::QC_3D> ( grid ),
      _mask_img ( mask_img ),
      /* _levelset( levelset ), */
      _hashInterior (),
      _iterMode ( EXT ) {}

  void setIteratorMode ( ITER_MODE mode ) { _iterMode = mode; }

  ITER_MODE getIteratorMode( ) const { return _iterMode; }

  bool existsInterior ( const qc::Element &el ) const {
    return ( _hashInterior.find ( el ) != _hashInterior.end() );
  }

  bool checkMask ( const qc::Element &el ) {
    qc::ElementNodeIterator<qc::QC_3D> nit = this->elementNodeBegin<qc::QC_3D> ( el );
    for ( ; nit != this->elementNodeEnd<qc::QC_3D> ( el ); ++nit ) {
      if ( _mask_img.get ( *nit ) >= 0 ) return true;
    }
    return false;
  }

  bool checkBoundaryCell ( const qc::Element &el ) {
    qc::ElementNodeIterator<qc::QC_3D> nit = this->elementNodeBegin<qc::QC_3D> ( el );
    for ( ; nit != this->elementNodeEnd<qc::QC_3D> ( el ); ++nit ) {
      if ( _mask_img.get ( *nit ) < .0 ) return true;
    }
    return false;
  }

  void insert ( const qc::Element &el ) {
    if ( checkMask ( el ) ) _hash.insert ( el );
  }

  void extractEdges( );
  eiterator eDbegin( ) const { return _edgesD.begin();  }
  eiterator eDend( ) const { return _edgesD.end(); }
  biterator bDbegin( ) const { return _boundaryNodesD.begin(); }
  biterator bDend( ) const { return _boundaryNodesD.end(); }

  biterator oBbegin( ) const { return _outerBoundaryNodes.begin(); }
  biterator oBend( ) const { return _outerBoundaryNodes.end(); }

  virtual const HashType & getHashSet( ) const {
    if ( _iterMode == INT )
      return _hashInterior;
    else
      return _hash;
  }
};

//-----------------------------------------------------------------------------------------------------------------------------

} // end of namespace nb.

#endif
