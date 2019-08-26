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

#ifndef __NARROWBANDCONFIGURATORS_H
#define __NARROWBANDCONFIGURATORS_H

#include <aol.h>
#include <quoc.h>
#include <gridSize.h>
#include <fastUniformGridMatrix.h>
#include <FEOpInterface.h>
#include <smallVec.h>
#include <narrow.h>
#include <indexMapper.h>

#include <subGridSparseMatrix.h>
#include <tfeBaseFunctionSet.h>

namespace nb {

//-----------------------------------------------------------------------------------------------------------------------------

template <typename _FullGridConfType, typename _InitType>
class FullGridConfiguratorHull;

//-----------------------------------------------------------------------------------------------------------------------------

template <typename _FullGridConfType, typename _InitType>
class NodeIterator {
public:
  typedef _InitType                                             GridType;
  typedef NodeIterator<_FullGridConfType, _InitType>            _Self;

  typedef qc::CoordType                                         IteratedType;
  typedef FullGridConfiguratorHull<_FullGridConfType, GridType> BeginType;
  typedef qc::EndElement                                        EndType;

  static const qc::Dimension Dim = _FullGridConfType::Dim;

  NodeIterator ( const BeginType & Configurator )
  : _mask ( Configurator.getInitializer().getSize() )
  , _size ( Configurator.getInitializer().getSize() )
  , _done ( false ) {
    Configurator.writeNodeExistMaskTo ( _mask );
  }

  const _Self & operator= ( const BeginType & Configurator ) {
    _done = false;
    _mask.resize ( Configurator.getFullGrid() );
    _size = Configurator.getInitializer().getSize();
    Configurator.writeNodeExistMaskTo ( _mask );
    return *this;
  }

  bool operator== ( const EndType & ) const {
    return _done;
  }

  bool operator!= ( const EndType & ) const {
    return !_done;
  }

  const _Self & operator++ () {
    do {
      _cur.xref() ++;
      if ( _cur.x() == _size.x() ) {
        _cur.xref() = 0;
        _cur.yref() ++;
        if ( _cur.y() == _size.y() ) {
          _cur.yref() = 0;
          _cur.zref() ++;
          if ( _cur.z() == _size.z() )
            _done = true;
        }
      }
    } while ( !(_mask.get(_cur) || _done) );
    return *this;
  }

  IteratedType& operator*() {
    return _cur;
  }

  IteratedType* operator->() {
    return &_cur;
  }

protected:
  typename qc::BitArray<Dim> _mask;
  IteratedType _cur;
  qc::CoordType _size;
  bool _done;
};

//-----------------------------------------------------------------------------------------------------------------------------

template <typename _FullGridConfType, typename _InitType>
class FullGridConfiguratorHullWithFullGridMatrix : protected _FullGridConfType {
public:
  typedef _FullGridConfType                                       Base;
  typedef FullGridConfiguratorHullWithFullGridMatrix<Base, _InitType> Self;
  typedef _InitType                                               InitType;

  typedef Base                                                    FullGridConfiguratorType;
  typedef typename Base::RealType                                 RealType;
  typedef typename Base::ElementType                              ElementType;
  typedef typename qc::BitArray<Base::DimOfWorld>                 MaskType;
  typedef typename Base::BaseFuncSetType                          BaseFuncSetType;
  typedef typename Base::QuadType                                 QuadType;

  typedef typename Base::VectorType                               VectorType;
  typedef typename Base::ArrayType                                ArrayType;
  typedef typename Base::MatrixType                               MatrixType;
  typedef typename Base::FullMatrixType                           FullMatrixType;

  typedef typename Base::VecType                                  VecType;
  typedef typename Base::DomVecType                               DomVecType;
  typedef typename Base::MatType                                  MatType;

  typedef typename InitType::OldFullElementIterator                  ElementIteratorType;
  typedef NodeIterator<_FullGridConfType, InitType>               NodeIteratorType;

  typedef typename InitType::BeginIterType                        BeginIterType;
  typedef typename InitType::EndIterType                          EndIterType;

  typedef typename InitType::Face                                 FaceType;

  // make inherited static member variables public:
  using Base::Dim;
  using Base::DomDim;
  using Base::DimOfWorld;
  using Base::IndexMode;
  using Base::maxNumLocalDofs;

  //! usual constructor
  FullGridConfiguratorHullWithFullGridMatrix ( const InitType & NBGrid )
      : _FullGridConfType ( NBGrid.getFullGrid() )
      , _NBGrid ( NBGrid )
  {}

  //! constructor only for TFE configurator subclass
  FullGridConfiguratorHullWithFullGridMatrix ( const InitType & NBGrid,
                                               const VectorType & LevelValues,
                                               RealType BandRadius )
      : _FullGridConfType ( NBGrid.getFullGrid(), LevelValues, BandRadius )
      , _NBGrid ( NBGrid )
  {}

  //! copy constructor
  //!
  //! copies everything like the implicit copy constructor,
  //! only the Mask pointer is not copied and will be newly
  //! created when accessed on the copy for the first time.
  FullGridConfiguratorHullWithFullGridMatrix ( const Self & other )
      : _FullGridConfType ( other )
      , _NBGrid ( other._NBGrid )
  {}

  const InitType & getInitializer () const {
    return _NBGrid;
  }

  const typename InitType::FullGridType & getFullGrid () const {
    return _NBGrid.getFullGrid ();
  }

  const BeginIterType & begin () const {
    return _NBGrid.begin ();
  }

  const EndIterType end () const {
    return _NBGrid.end ();
  }

  void writeNodeExistMaskTo ( MaskType & existMask ) const {
    existMask.setAll ( false );
    for ( ElementIteratorType iter = this->begin(); iter != this->end(); ++iter ) {
      for ( int i = 0; i < getNumLocalDofs ( *iter ); ++i )
        existMask.set ( localToGlobal ( *iter, i ), true );
    }
  }

  //! computes and returns node exist mask
  const MaskType & getNodeExistMask () const {
    _nodeExistMask.reset ( new MaskType ( qc::GridSize<Dim> ( this->getInitializer() ) ) );
    writeNodeExistMaskTo ( *_nodeExistMask );
    return *_nodeExistMask;
  }

  // make inherited functions public:
  using Base::localToGlobal;
  using Base::vol;
  using Base::H;
  using Base::getNumLocalDofs;
  using Base::getNumGlobalDofs;
  using Base::getBaseFunctionSet;
  using Base::createNewMatrix;
  using Base::maxNumQuadPoints;

protected:
  const InitType &           _NBGrid;
  mutable auto_ptr<MaskType> _nodeExistMask;
};

//-----------------------------------------------------------------------------------------------------------------------------

template <typename _FullGridConfType, typename _InitType>
class FullGridConfiguratorHull
  : public FullGridConfiguratorHullWithFullGridMatrix<_FullGridConfType, _InitType> {

  typedef FullGridConfiguratorHullWithFullGridMatrix<_FullGridConfType, _InitType> Base;

public:
  // unfortunately, no types are inherited.
  typedef typename Base::RealType                 RealType;
  typedef typename Base::InitType                 InitType;
  typedef typename Base::QuadType                 QuadType;
  typedef typename Base::VecType                  VecType;
  typedef typename Base::DomVecType               DomVecType;
  typedef typename Base::MatType                  MatType;
  typedef typename Base::ArrayType                ArrayType;
  typedef typename Base::MaskType                 MaskType;
  typedef typename Base::BaseFuncSetType          BaseFuncSetType;
  typedef typename Base::FullGridConfiguratorType FullGridConfiguratorType;
  typedef typename Base::ElementType              ElementType;
  typedef typename Base::ElementIteratorType      ElementIteratorType;
  typedef typename Base::NodeIteratorType         NodeIteratorType;
  typedef typename Base::BeginIterType            BeginIterType;
  typedef typename Base::EndIterType              EndIterType;
  typedef typename Base::FullMatrixType           FullMatrixType;
  typedef typename Base::VectorType               VectorType;
  typedef typename Base::FaceType                 FaceType;

  // this is the only one we really want:
  typedef SubGridSparseMatrix<RealType, Base >    MatrixType;

  //! usual constructor
  FullGridConfiguratorHull ( const _InitType &NarrowGrid ) :
      FullGridConfiguratorHullWithFullGridMatrix<_FullGridConfType, _InitType> ( NarrowGrid ) {}

  //! constructor only for TFE configurator subclass
  FullGridConfiguratorHull ( const InitType & NBGrid,
                             const VectorType & LevelValues,
                             RealType BandRadius ) :
      FullGridConfiguratorHullWithFullGridMatrix<_FullGridConfType, _InitType> ( NBGrid, LevelValues, BandRadius ) {}

  //! create a new, clean matrix
  MatrixType * createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( *this );
    mat->setZero( );
    return mat;
  }

  int getConsecutiveElementNumber ( const ElementType &/*El*/ ) const {
    throw aol::UnimplementedCodeException( "FullGridConfiguratorHull::getConsecutiveElementNumber", __FILE__, __LINE__);
    return -1;
  }
};

//-----------------------------------------------------------------------------------------------------------------------------

/** This class is only provided for backward compatibility and supplies no relevant
 *  functionality over its base class
 */
template < typename _RealType,
qc::Dimension Dim,
typename _QuadType,
typename _InitType >
class NarrowBandConfiguratorTraitMultiLin
  : public FullGridConfiguratorHull<qc::QuocConfiguratorTraitMultiLin<_RealType, Dim, _QuadType>,
                                    _InitType> {

public:
  // **** typedefs ****
  typedef FullGridConfiguratorHull<qc::QuocConfiguratorTraitMultiLin<_RealType, Dim, _QuadType>, _InitType> Base;

  typedef typename Base::RealType                 RealType;
  typedef typename Base::InitType                 InitType;
  typedef typename Base::QuadType                 QuadType;
  typedef typename Base::VectorType               VectorType;
  typedef typename Base::MatrixType               MatrixType;
  typedef typename Base::VecType                  VecType;
  typedef typename Base::DomVecType               DomVecType;
  typedef typename Base::MatType                  MatType;
  typedef typename Base::ArrayType                ArrayType;
  typedef typename Base::MaskType                 MaskType;
  typedef typename Base::BaseFuncSetType          BaseFuncSetType;
  typedef typename Base::FullGridConfiguratorType FullGridConfiguratorType;
  typedef typename Base::ElementType              ElementType;
  typedef typename Base::ElementIteratorType      ElementIteratorType;
  typedef typename Base::FullMatrixType           FullMatrixType;
  typedef typename Base::FaceType                 FaceType;

  static const qc::Dimension DimOfWorld = Dim;

  // **** constructor ****
  NarrowBandConfiguratorTraitMultiLin ( const _InitType &NarrowGrid )
      : FullGridConfiguratorHull<qc::QuocConfiguratorTraitMultiLin<_RealType, Dim, _QuadType>, _InitType> ( NarrowGrid )
  {}
};

//-----------------------------------------------------------------------------------------------------------------------------

/** This class is only provided for backward compatibility and supplies no relevant
 *  functionality over its base class
 */
template < typename _RealType,
qc::Dimension Dim,
typename _QuadType,
typename _InitType >
class NarrowBandConfiguratorTraitMultiQuad
  : public FullGridConfiguratorHull<qc::QuocConfiguratorTraitMultiQuad<_RealType, Dim, _QuadType>,
                                    _InitType> {


public:
  // **** typedefs ****
  typedef FullGridConfiguratorHull<qc::QuocConfiguratorTraitMultiQuad<_RealType, Dim, _QuadType>, _InitType> Base;

  typedef typename Base::RealType                 RealType;
  typedef typename Base::InitType                 InitType;
  typedef typename Base::QuadType                 QuadType;
  typedef typename Base::VectorType               VectorType;
  typedef typename Base::MatrixType               MatrixType;
  typedef typename Base::VecType                  VecType;
  typedef typename Base::DomVecType               DomVecType;
  typedef typename Base::MatType                  MatType;
  typedef typename Base::ArrayType                ArrayType;
  typedef typename Base::MaskType                 MaskType;
  typedef typename Base::BaseFuncSetType          BaseFuncSetType;
  typedef typename Base::FullGridConfiguratorType FullGridConfiguratorType;
  typedef typename Base::ElementType              ElementType;
  typedef typename Base::ElementIteratorType      ElementIteratorType;
  typedef typename Base::FullMatrixType           FullMatrixType;
  typedef typename Base::FaceType                 FaceType;

  // **** constructor ****
  NarrowBandConfiguratorTraitMultiQuad ( const _InitType &NarrowGrid )
      : FullGridConfiguratorHull<qc::QuocConfiguratorTraitMultiQuad<_RealType, Dim, _QuadType>, _InitType> ( NarrowGrid )
  {}
};

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType,
          qc::Dimension Dim,
          typename QuadType,
          typename InitType>
class QuocConfiguratorTraitTFE {};

//-----------------------------------------------------------------------------------------------------------------------------

//! configurator for cut-off linear basis functions on tetrahedra
//! \author von Deylen (feb 2009)
template <typename _RealType, typename _QuadType, typename _InitType>
class QuocConfiguratorTraitTFE<_RealType, qc::QC_3D, _QuadType, _InitType>
  : public qc::QuocConfiguratorTraitBase<_RealType, qc::QC_3D> {

public:
  typedef QuocConfiguratorTraitTFE<_RealType, qc::QC_3D,
                                   _InitType, _QuadType>        Self;
  typedef qc::QuocConfiguratorTraitBase<_RealType, qc::QC_3D>   Base;
  typedef _QuadType                                             QuadType;
  typedef _RealType                                             RealType;
  typedef _InitType                                             InitType;
  typedef typename aol::Vector<RealType>                        VectorType;
  typedef typename Base::ElementType                            ElementType;
  typedef typename InitType::OldFullElementIterator                ElementIteratorType;
  typedef typename aol::FullMatrix<RealType>                    FullMatrixType;
  typedef typename Base::MaskType                               MaskType;
  typedef aol::Vec3<RealType>                                   VecType;
  typedef aol::Vec3<RealType>                                   DomVecType;
  typedef aol::Matrix33<RealType>                               MatType;
  typedef qc::MultilinFEBandMatrix<RealType, qc::QC_3D>         MatrixType;
  typedef typename qc::ScalarArray<RealType, qc::QC_3D>         ArrayType;
  typedef qc::BaseFunctionSetTFE<RealType, QuadType>            BaseFuncSetType;
  typedef qc::FastILexMapper<qc::QC_3D>                         IndexMapperType;
  typedef Self                                                  FullGridConfiguratorType;

  static const int maxNumLocalDofs      = 8;
  static const qc::Dimension Dim        = qc::QC_3D;
  static const qc::Dimension DimOfWorld = qc::QC_3D;
  static const qc::Dimension DomDim     = qc::QC_3D;
  static const aol::GridGlobalIndexMode IndexMode = InitType::IndexMode;

  QuocConfiguratorTraitTFE ( const InitType & Grid,
                             const VectorType & LevelValues,
                             RealType BandRadius )
  : qc::QuocConfiguratorTraitBase<_RealType, qc::QC_3D> ( Grid )
  , _indexMapper ( qc::GridSize<qc::QC_3D> ( Grid ) )
  , _baseFuncSet ( LevelValues, BandRadius, Grid.H() ) {

    ElementIteratorType elIter = Grid.begin();
    computeElementToNodeOffsets ( *elIter, _elementToNodeOffsets );
  }

  QuocConfiguratorTraitTFE ( const InitType & Grid )
  : qc::QuocConfiguratorTraitBase<_RealType, qc::QC_3D> ( Grid )
  , _indexMapper ( qc::GridSize<qc::QC_3D> ( Grid ) ) {

    ElementIteratorType elIter = Grid.begin();
    computeElementToNodeOffsets ( *elIter, _elementToNodeOffsets );
  }

  void destroy () {
    _baseFuncSet.destroy();
  }

  int maxNumQuadPoints () const {
    return _baseFuncSet.maxNumQuadPoints();
  }

  inline int getNumLocalDofs ( const ElementType & ) const {
    // actually, at each point inside the element, only four
    // basis functions may be nonzero. But as over the whole
    // element there are 8 basis functions to be taken into
    // account, we
    return 8;
  }

  int getNumGlobalDofs () const {
    return this->getInitializer().getNumberOfNodes();
  }
/*
  void splitGlobalIndex ( int globalIndex, int & x, int & y, int & z ) const {
    _indexMapper.splitGlobalIndex ( globalIndex, x, y, z );
  }
*/
  const BaseFuncSetType & getBaseFunctionSet ( const ElementType & El ) const {
    /* if a non-uniform grid is used, the array elementToNodeOffsets
       has to be computed before every initialization. In this case,
       one should have an additional initializeElement() routine to
       be called only when necessary. Up to now, the array may
       stay the same for all elements.
    */
    // computeElementToNodeOffsets ( el, _elementToNodeOffsets );

    // this method computes the virtual tetrahedra, the weights etc. pp.
    _baseFuncSet.setToElement ( localToGlobal ( El, 0 ), _elementToNodeOffsets );
    return _baseFuncSet;
  }

  void computeElementToNodeOffsets ( const ElementType & El, int * ElementToNodeOffsets ) {
    for (int i = 0; i < getNumLocalDofs ( El ); ++i)
      ElementToNodeOffsets[i] = localToGlobal ( El, i ) - localToGlobal ( El, 0 );
  }

  int localToGlobal ( const ElementType & El, int LocalIndex ) const {
    return _indexMapper.localToGlobal ( El, LocalIndex );
  }

  inline void localToGlobal ( const ElementType &El,
                              int localIndex0, int localIndex1,
                              aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  RealType vol ( const ElementType & El ) const {
    return getBaseFunctionSet( El ).getVolume() * aol::Cub ( this->getInitializer().H() );
  }

    //! create a new, clean matrix
  MatrixType * createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( this->getInitializer() );
    mat->setZero( );
    return mat;
  }

protected:
  IndexMapperType         _indexMapper;
  mutable BaseFuncSetType _baseFuncSet;
  int                     _elementToNodeOffsets[8];
};

//-----------------------------------------------------------------------------------------------------------------------------

/** This class is only provided for backward compatibility and supplies no relevant
 *  functionality over its base class
 */
template < typename _RealType,
qc::Dimension Dim,
typename _QuadType,
typename _InitType >
class NarrowBandConfiguratorTraitTFE
  : public FullGridConfiguratorHull<QuocConfiguratorTraitTFE<_RealType, Dim, _QuadType, typename _InitType::FullGridType>,
                                    _InitType> {

public:
  // **** typedefs ****
  typedef FullGridConfiguratorHull<QuocConfiguratorTraitTFE<_RealType, Dim, _QuadType, typename _InitType::FullGridType>, _InitType> Base;

  typedef typename Base::RealType                 RealType;
  typedef typename Base::InitType                 InitType;
  typedef typename Base::QuadType                 QuadType;
  typedef typename Base::VectorType               VectorType;
  typedef typename Base::MatrixType               MatrixType;
  typedef typename Base::VecType                  VecType;
  typedef typename Base::DomVecType               DomVecType;
  typedef typename Base::MatType                  MatType;
  typedef typename Base::ArrayType                ArrayType;
  typedef typename Base::MaskType                 MaskType;
  typedef typename Base::BaseFuncSetType          BaseFuncSetType;
  typedef typename Base::FullGridConfiguratorType FullGridConfiguratorType;
  typedef typename Base::ElementType              ElementType;
  typedef typename Base::ElementIteratorType      ElementIteratorType;
  typedef typename Base::FullMatrixType           FullMatrixType;
  typedef typename Base::FaceType                 FaceType;

  // **** constructor ****
  NarrowBandConfiguratorTraitTFE ( const InitType & NarrowGrid,
                                   const VectorType & LevelValues,
                                   RealType BandRadius )
  : FullGridConfiguratorHull
      < QuocConfiguratorTraitTFE<_RealType, Dim, _QuadType, typename _InitType::FullGridType>,
        _InitType> ( NarrowGrid, LevelValues, BandRadius )
  {}

  NarrowBandConfiguratorTraitTFE ( const InitType & NarrowGrid )
  : FullGridConfiguratorHull
      < QuocConfiguratorTraitTFE<_RealType, Dim, _QuadType, typename _InitType::FullGridType>,
        _InitType> ( NarrowGrid )
  {}

  using Base::destroy;
};

//-----------------------------------------------------------------------------------------------------------------------------

} // end namespace nb

#endif
