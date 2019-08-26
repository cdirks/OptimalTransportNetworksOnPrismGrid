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

#ifndef __TPCFEELEMENT_H
#define __TPCFEELEMENT_H

#include <tpCFEBasics.h>
#include <tpCFELookup.h>
#include <tpCFETetrahedron.h>
#include <quoc.h>

namespace tpcfe {

/** A class which stores the minimum information needed to identify a cfe element.
 */

class CFETopoElement : public qc::Element {
protected:
  qc::GridSize<qc::QC_3D> _gridSize;
  CFEType _cfeType;

public:
  //! Default Constructor
  CFETopoElement ( ) : qc::Element(), _gridSize ( 0, 0, 0 ), _cfeType ( ) {
  }

  CFETopoElement ( const CFEType CfeType ) : qc::Element(), _gridSize ( 0, 0, 0 ), _cfeType ( CfeType ) {
  }

  // \todo GSize or ref? constructor used?
  CFETopoElement ( const int X, const int Y, const int Z, const qc::GridSize<qc::QC_3D> &GSize, const CFEType CfeType )
    : qc::Element ( X, Y, Z, 0, 0 ), _gridSize ( GSize ), _cfeType ( CfeType ) {
  }

  CFETopoElement ( const qc::Element &qcEl, const qc::GridSize<qc::QC_3D> &GSize, const CFEType &CfeType )
    : qc::Element ( qcEl ), _gridSize ( GSize ), _cfeType ( CfeType ) {
  }

  //! Copy constructor
  CFETopoElement ( const CFETopoElement& other ) : qc::Element ( other ), _gridSize ( other._gridSize ), _cfeType ( other._cfeType ) {
  }

  //! Assignment operator
  CFETopoElement& operator= ( const CFETopoElement& other ) {
    // self-assignment should be no problem.
    qc::Element::operator= ( other );
    _gridSize = other._gridSize;
    _cfeType = other._cfeType;

    return ( *this );
  }


  /** Return the cfe type including structure information (in upper bits)
   */
  CFEType cfeType () const  {
    return _cfeType;
  }

  unsigned char pureCFEType ( ) const {
    return ( _cfeType._pureType );
  }

  /** Return the global index of this element
   */
  inline int globalIndex () const  {
    return qc::ILexCombine3 ( coords[0], coords[1], coords[2], _gridSize[0], _gridSize[1] );
  }

  void set ( const qc::Element &qcEl, const qc::GridSize<qc::QC_3D> &GSize, const CFEType &CfeType ) {
    qc::Element::operator= ( qcEl );
    _gridSize = GSize;
    _cfeType = CfeType;
  }

  void setCfeType ( const CFEType CfeType ) {
    _cfeType = CfeType;
  }

  inline signed char getNonInterfacedSign ( ) const {
    return ( _cfeType.signForNonInterfacedType() );
  }
};


/** Cubic Elements for CFE computations containing geometric data (and some more topological data we need less often).
 *  Information can be recomputed by calling
 *  - computeAssembleData
 */
template < typename RealType  >
class CFEElement : public CFETopoElement {
  friend class CFETetraInElementIterator<RealType>;

protected:
  RealType _structureValue[8];
  RealType _cutRelation[8][8];
  std::vector< tpcfe::CFETetra<RealType> > _tetraVec;
  bool _assembleDataComputed;

public:
  /** Default constructor
   */
  CFEElement ( CFEType CfeType = CFEType() ) : CFETopoElement ( CfeType ), _tetraVec(), _assembleDataComputed ( false ) {
    memset ( _structureValue, 0, sizeof ( _structureValue ) );
    for ( short i = 0; i < 8; ++i )
      for ( short j = 0; j < 8; ++j )
        _cutRelation[i][j] = aol::NumberTrait<RealType>::NaN;
  }

  /** Conversion constructor
   */
  explicit CFEElement ( const CFETopoElement &sparseOther ) : CFETopoElement ( sparseOther ), _tetraVec(), _assembleDataComputed ( false ) {
    memset ( _structureValue, 0, sizeof ( _structureValue ) );
    for ( short i = 0; i < 8; ++i )
      for ( short j = 0; j < 8; ++j )
        _cutRelation[i][j] = aol::NumberTrait<RealType>::NaN;
  }

  /** Initialize from a vec3
   */
  explicit CFEElement ( const aol::Vec3<RealType> &vec ) : CFETopoElement ( static_cast<int> ( vec[0] ), static_cast<int> ( vec[1] ), static_cast<int> ( vec[2] ), qc::GridSize<qc::QC_3D>( 0, 0, 0 ), CFEType() ), _tetraVec(), _assembleDataComputed ( false ) {
    memset ( _structureValue, 0, sizeof ( _structureValue ) );
    for ( short i = 0; i < 8; ++i )
      for ( short j = 0; j < 8; ++j )
        _cutRelation[i][j] = aol::NumberTrait<RealType>::NaN;
  }

  /** Copy constructor
   */
  explicit CFEElement ( const CFEElement & other ) : CFETopoElement ( other ), _tetraVec( other._tetraVec ), _assembleDataComputed ( other._assembleDataComputed ) {
    memcpy ( _cutRelation, other._cutRelation, sizeof ( _cutRelation ) );
    memcpy ( _structureValue, other._structureValue, sizeof ( _structureValue ) );
  }

  /** Assignment operator
   */
  CFEElement & operator= ( const CFEElement & other ) {
    // self-assignment should be no problem.
    CFETopoElement::operator= ( other );
    memcpy ( _cutRelation, other._cutRelation, sizeof ( _cutRelation ) );
    memcpy ( _structureValue, other._structureValue, sizeof ( _structureValue ) );
    _tetraVec = other._tetraVec;
    _assembleDataComputed = other._assembleDataComputed;
    return *this;
  }

  //! \todo think about &GSize - is this constructor used?
  CFEElement ( const qc::GridSize<qc::QC_3D> &GSize, const int X, const int Y, const int Z = 0,
               const CFEType CfeType = CFEType() ) :
      CFETopoElement ( X, Y, Z, GSize, CfeType ), _tetraVec(), _assembleDataComputed ( false ) {
    memset ( _structureValue, 0, sizeof ( _structureValue ) );

    for ( short i = 0; i < 8; ++i )
      for ( short j = 0; j < 8; ++j )
        _cutRelation[i][j] = aol::NumberTrait<RealType>::NaN;
  }

  /** \todo document
   */
  CFEElement ( const qc::Element &qcEl, const qc::GridSize<qc::QC_3D> &GSize, const CFEType &CfeType )
    : CFETopoElement ( qcEl, GSize, CfeType ), _tetraVec(), _assembleDataComputed ( false ) {
    memset ( _structureValue, 0, sizeof ( _structureValue ) );

    for ( short i = 0; i < 8; ++i )
      for ( short j = 0; j < 8; ++j )
        _cutRelation[i][j] = aol::NumberTrait<RealType>::NaN;
  }

  void set ( const qc::Element &qcEl, const qc::GridSize<qc::QC_3D> &GSize, const CFEType &CfeType ) {
    CFETopoElement::set ( qcEl, GSize, CfeType );
    memset ( _structureValue, 0, sizeof ( _structureValue ) );

    for ( short i = 0; i < 8; ++i )
      for ( short j = 0; j < 8; ++j )
        _cutRelation[i][j] = aol::NumberTrait<RealType>::NaN;

    _tetraVec.clear();
    _assembleDataComputed = false;
  }

  inline void setStructureValue ( const int i, const RealType val ) {
#ifdef BOUNDS_CHECK
    if ( i < 0 || i > 7 )
      throw aol::OutOfBoundsException ( "index out of bounds", __FILE__, __LINE__ );
#endif
    _structureValue[i] = val;
  }

  inline void setCutRelation ( const int i, const int j, const RealType val ) {
    _cutRelation[i][j] = val;
  }

  inline qc::CoordType getGlobalPos ( const int i ) const {
    return ( *this + CFETopoLookup::_hexNbOffsets[i] );
  }

  inline RealType getStructureValue ( const int i ) const {
#ifdef BOUNDS_CHECK
    if ( i < 0 || i > 7 )
      throw aol::OutOfBoundsException ( "index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _structureValue[i] );
  }

  inline RealType getCutRelation ( const int i, const int j ) const {
#ifdef BOUNDS_CHECK
    if ( i < 0 || i > 7 || j < 0 || j > 7 )
      throw aol::OutOfBoundsException ( "index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _cutRelation[i][j] );
  }

  /** Set new grid width
   */
  void setGridSize ( const qc::GridSize<qc::QC_3D> &GSize ) {
    _gridSize = GSize;
  }

  /** Return the local index of node i of this element.
   */
  inline int globalIndex ( const int i ) const {
#ifdef BOUNDS_CHECK
    if ( i < 0 || i > 7 )
      throw aol::OutOfBoundsException ( "index out of bounds", __FILE__, __LINE__ );
#endif
    return ( this->globalIndex() + qc::ILexCombine3 ( tpcfe::CFETopoLookup::_hexNbOffsets[i][0], tpcfe::CFETopoLookup::_hexNbOffsets[i][1], tpcfe::CFETopoLookup::_hexNbOffsets[i][2], _gridSize[0], _gridSize[1] ) );
  }

  using tpcfe::CFETopoElement::globalIndex;

  /** Print information about this CFE Element
   */
  void print () const  {
    cerr << "GlobalIndex: ";
    for ( short i = 0; i < 8; ++i )
      cerr << globalIndex(i) << " ";
    qc::Element::print ();
    dumpCutRelations();
    dumpStructureValues();
  }


  /** Print the cut relations to a stream
   */
  void dumpCutRelations ( ostream &out = cerr ) const {
    for ( int i = 0; i < 8; ++i ) {
      for ( int j = 0; j < 8; ++j ) {
        out << _cutRelation[i][j] << " ";
      }
      out << endl;
    }
  }

  void dumpStructureValues ( ostream &out = cerr ) const {
    for ( int i = 0; i < 8; ++i ) {
      out << _structureValue[i] << " ";
    }
    out << endl;
  }


  /** determine the vertices \f$ v_0,\dots, v_3 \f$ which span the given tetrahedron
   *  IMPORTANT: The cut relations must have been computed before use of this method
   */
  void getTetraVertices ( CFETetra<RealType> &t ) const;

  /** Compute the normal of the interface approximated in the ith regular tetrahedron
   */
  aol::Vec3<RealType> interfaceNormalOnRegTetra ( const int i ) const;

  /** Compute all data needed for assembly of matrices. This is
   *  - cut relations on the interfaces
   *  - volumes of the sub-tetrahedra
   *  - vertices coordinates of the sub-tetrahedra
   *  This method makes all computations for grid-width 1. I.e. a
   *  standard tetrahedron has size 1./6. (as set in CFEGrid.cpp) and
   *  the volume of a tetrahedron which is computed here has to be seen
   *  relative to the unit cube. Thus, for assembly of matrices a
   *  scaling with the grid-width in the correct power has to be implemented.
   *
   */
  template < tpcfe::ConstraintType CT, typename NodalCoeffType >
  void computeAssembleData ( const CFEGridBase< RealType, CT, NodalCoeffType > &grid ) {

    if (  _cfeType.representsInterfaced() ) {
      grid.getStructureValues ( *this );
      grid.getCutRelations ( *this );
    }

    _tetraVec.clear();
    for ( CFETopoTetraIterator it ( _cfeType ); it.notAtEnd(); ++it ) {
      _tetraVec.push_back ( CFETetra<RealType> ( *it ) );
    }

    for ( typename std::vector< tpcfe::CFETetra<RealType> >::iterator it = _tetraVec.begin(); it != _tetraVec.end(); ++it ) {
      getTetraVertices ( *it );
      it->setVolume ( it->determinant() / 6 );
    }

    _assembleDataComputed = true;
  }
};


/** Iterator iterating over (regular or virtual) tetrahedra inside one CFE element, depending on type, can iterate over all, positive and negative tetrahedra.
 *  \author Preusser, Schwen
 */
template< typename RealType >
class CFETetraInElementIterator {
protected:
  typename std::vector< tpcfe::CFETetra<RealType> >::const_iterator _startIt, _curIt, _endIt;
  const tpcfe::CFEElement<RealType> &_el;
  signed char _sign;

public:
  CFETetraInElementIterator ( const tpcfe::CFEElement<RealType> &el, const signed char sign = 0 ) : _el ( el ), _sign ( sign ) {
    if ( el._assembleDataComputed == false ) {
      throw aol::Exception ( " tpcfe::CFEElement::CFETetraInElementIterator cannot be initialized without prior computation of assemble data", __FILE__, __LINE__ );
    }
    restart ( sign );
  }

  // copy constructor does the correct thing; assignment operator?

  /** Check whether the end has been reached
   */
  bool notAtEnd() const {
    return ( ! ( atEnd() ) );
  }

  /** Check whether the end has been reached
   */
  bool atEnd() const {
    return ( _curIt == _endIt );
  }

  /** prefix-increment operator
   */
  void operator++ () {
    ++_curIt;
    findNextSignValidTet();
  }

  /** Return reference to current tetrahedron of the iteration
   */
  const CFETetra<RealType>& operator* () const {
    if ( _curIt == _endIt ) {
      throw aol::Exception ( "tpcfe::CFETopoTetraInElementIterator::operator*: cannot dereference *endIt", __FILE__, __LINE__ );
    }
    return ( *_curIt );
  }

  const CFETetra<RealType>* operator-> () const {
    if ( _curIt == _endIt ) {
      throw aol::Exception ( "tpcfe::CFETopoTetraInElementIterator::operator->: cannot dereference *endIt", __FILE__, __LINE__ );
    }
    return ( &(*_curIt) );
  }

  /** restart existing iterator
   */
  void restart ( const signed char sign = 0 ) {
    _sign = sign;
    _startIt = _el._tetraVec.begin();
    _curIt   = _startIt;
    _endIt   = _el._tetraVec.end();
    findNextSignValidTet();
  }

 protected:
  inline void findNextSignValidTet ( ) {
    if ( _sign == 0 ) {
      return;
    }
    while ( ( _curIt != _endIt ) && ( _curIt->getSign() != _sign ) ) {
      ++_curIt;
    }
  }

  CFETetraInElementIterator();
  // copy constructor and assignment operator do the correct thing, but are probably slow
  CFETetraInElementIterator<RealType>( const CFETetraInElementIterator<RealType>& );
  CFETetraInElementIterator<RealType>& operator= ( const CFETetraInElementIterator<RealType>& );

  // end class CFETetraInElementIterator
};


}
#endif
