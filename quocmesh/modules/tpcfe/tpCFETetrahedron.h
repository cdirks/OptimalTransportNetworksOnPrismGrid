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

#ifndef __TPCFETETRAHEDRON_H
#define __TPCFETETRAHEDRON_H

#include <quoc.h>
#include <smallMat.h>

#include <tpCFEBasics.h>
#include <tpCFELookup.h>

namespace tpcfe {

//! (i,j)=1 <=> local node j is contained in stdTetra i

const bool __CFE_LOCAL_IN_STDTETRA[8][6] = { {1, 1, 1, 0, 0, 0},
                                             {1, 1, 0, 1, 1, 1},
                                             {1, 0, 0, 1, 1, 0},
                                             {0, 0, 0, 1, 0, 0},
                                             {0, 0, 1, 0, 0, 0},
                                             {0, 1, 1, 0, 0, 1},
                                             {1, 1, 1, 0, 1, 1},
                                             {0, 0, 0, 1, 1, 1}
                                           };

/** \brief Basis class for CFE tetrahedra containing topological data that can usefully be stored in the lookup tables.
 *
 *  In the composite finite element (CFE) context, a cube-element
 *  disaggregates into tetrahedra.
 *  If a tetrahedron-node lies on an edge between two cube-nodes a and b,
 *  it is described be the index pair (a,b) and referred to as a <b>virtual node<\b>.
 *  -> the 4 tetrahedral nodes are stored in array[4][2]
 *
 *  There is a standard-partition of a cube in 6 tetrahedra of equal Volume.
 *  Those tetrahedra are referred to as std-tetrahedra and are numbered from 0 to 5
 *
 *  If a std-tetrahedron disaggregates into further tetrahedra ( sub-tetrahedra )
 *  it is referred to as the parent.
 *
 *  \author Preusser, Schwen
 */
class CFETopoTetra {
public:
  enum VolType {  P = 0,                // <- does not split up and is positive

                  PPPN  = 1,             // <- splits up in 3 positve and 1 negative
                  PPPN1 = 2,
                  PPPN2 = 3,
                  PPPN3 = 4,

                  NNNP  = 5,             // <-  splits up in 1 positve and 3 negative
                  NNNP1 = 6,
                  NNNP2 = 7,
                  NNNP3 = 8,

                  NNNPPP  = 9,           // <-  splits up in 3 positve and 3 negative
                  NNNPPP1 = 10,
                  NNNPPP2 = 11,
                  NNNPPP3 = 12,
                  NNNPPP4 = 13,
                  NNNPPP5 = 14,

                  PPPNNN  = 15,          // <-  splits up in 3 negative and 3 positve
                  PPPNNN1 = 16,
                  PPPNNN2 = 17,
                  PPPNNN3 = 18,
                  PPPNNN4 = 19,
                  PPPNNN5 = 20,

                  N = 21,               // <- does not split up and is negative

                  NO_TYPE = 22          // <- no type assigned
               };
protected:
  VolType _volType;             //!< inidicates type of tetra for computation of volume
  short _virtualNodeNum;        //!< Number of virtual nodes in this tetrahedron
  signed char _sign;            //!< On which side of the interface is this tetra
  short _parent;                //!< references corresponding standard-tetra
  short _node[4][2];            //!< stores local indices (second entry != NODE_NOT_VIRTUAL indicates a virtual node)

  static const char volumeTypeName[23][8];

public:
  CFETopoTetra ( ) : _volType ( NO_TYPE ), _virtualNodeNum ( 0 ), _sign ( 0 ), _parent ( 11 ) {
    _node[0][0] = _node[0][1] = NODE_NOT_VIRTUAL;
    _node[1][0] = _node[1][1] = NODE_NOT_VIRTUAL;
    _node[2][0] = _node[2][1] = NODE_NOT_VIRTUAL;
    _node[3][0] = _node[3][1] = NODE_NOT_VIRTUAL;
  }

  //! Create a tetrahedron with the given data
  CFETopoTetra ( const short x0, const short y0,
                 const short x1, const short y1,
                 const short x2, const short y2,
                 const short x3, const short y3,
                 const short Parent, const VolType volType, const short virtualNodeNum, const signed char sign ) : _volType ( volType ), _virtualNodeNum ( virtualNodeNum ), _sign ( sign ),  _parent ( Parent ) {
    _node[0][0] = x0;    _node[0][1] = y0;
    _node[1][0] = x1;    _node[1][1] = y1;
    _node[2][0] = x2;    _node[2][1] = y2;
    _node[3][0] = x3;    _node[3][1] = y3;
  }

  //! copy constructor
  CFETopoTetra ( const CFETopoTetra &other ) : _volType ( other._volType ), _virtualNodeNum ( other._virtualNodeNum ), _sign ( other._sign ), _parent ( other._parent ) {
    memcpy ( _node, other._node, sizeof ( _node ) );
  }

  //! Assignment operator
  CFETopoTetra & operator= ( const CFETopoTetra &other ) {
    // Beware of self-assignment
    if ( this == &other )
      return *this;

    _volType = other._volType;
    _virtualNodeNum = other._virtualNodeNum;
    _sign = other._sign;
    _parent = other._parent;
    memcpy ( _node, other._node, sizeof ( _node ) );

    return *this;
  }

  //! Set parent tetra of this one
  inline void setParent ( const short p ) { _parent = p; }

  //! Get parent tetra of this one
  inline short getParent () const  { return _parent; }

  inline void setVolType ( const VolType t ) { _volType = t;}

  inline VolType getVolType () const  { return _volType;}

  //! Returns the sign of this tetrahedron, i.e. the side on which this tetra lies
  inline signed char getSign() const { return _sign; }

  //! Set the sign of this tetrahedron, i.e. the side on which this tetra lies
  inline void setSign ( const signed char sign ) { _sign = sign; }

  //! Checks whether the node with index k is a virtual node
  inline bool isVirtualNode ( const short k ) const { return _node[k][1] != NODE_NOT_VIRTUAL; }

  inline void setNodes ( const short x0, const short y0, const short x1, const short y1, const short x2, const short y2, const short x3, const short y3 ) {
    _node[0][0] = x0;    _node[0][1] = y0;
    _node[1][0] = x1;    _node[1][1] = y1;
    _node[2][0] = x2;    _node[2][1] = y2;
    _node[3][0] = x3;    _node[3][1] = y3;
  }

  //! Sets  k-th node to x or to virtual node x,y if argument y is supplied and not equal to NODE_NOT_VIRTUAL
  inline void setNode ( const short k, const short x, const short y = NODE_NOT_VIRTUAL ) {
    _node[k][0] = x;
    _node[k][1] = y;
  }

  //! returns first component of k-th node
  inline short getNode ( const short k ) const {
    return _node[k][0];
  }

  //! returns j-th component of i-th node
  inline short operator () ( const short i, const short j ) const {
#ifdef BOUNDS_CHECK
    if ( i < 0 || i > 3 || j < 0 || j > 1 )
      throw aol::OutOfBoundsException ( "tpcfe::CFETetrahedron::operator( const short i, const short j ) must have i = 0,1,2,3 and j = 0,1", __FILE__, __LINE__ );
#endif
    return _node[i][j];
  }

  //! returns const pointer to i-th node
  inline const short *operator () ( const short i ) const { return _node[i]; }

  //! true if tetrahedron is one of 6 stdTetrahedra
  inline bool isStandard () const { return ( _volType == P ) || ( _volType == N ); }

  //! returns true if LocIndex appears in description or tetrahedron
  bool containsLocalNode ( const short LocIndex ) const {
    if ( isStandard () )
      return __CFE_LOCAL_IN_STDTETRA[LocIndex][ getParent() ];
    else
      return ( _node[0][0] == LocIndex || _node[0][1] == LocIndex ||
               _node[1][0] == LocIndex || _node[1][1] == LocIndex ||
               _node[2][0] == LocIndex || _node[2][1] == LocIndex ||
               _node[3][0] == LocIndex || _node[3][1] == LocIndex   );
  }

  //! returns position of LocIndex in tetrahedron desciption
  short localToIndex ( const short LocIndex, const short LocIndex2 = NODE_NOT_VIRTUAL ) const {
    for ( short i = 0; i < 4; ++i )
      if ( _node[i][0] == LocIndex && _node[i][1] == LocIndex2 )
        return i;
    throw aol::Exception ( "CFETetra::localToIndex: no occurance of LocIndex in tetra", __FILE__, __LINE__ );
    return -1;
  }

  //! returns a string representation of the volume type
  const char* volTypeToString() const {
    return volumeTypeName[_volType];
  }

  short getVirtualNodeNum() const {
    return ( _virtualNodeNum );
  }

  // swaps indices 0 and 1 in permutation of indices and  in assignment of tetra-generating-vectors
  //( needed to set orientation to RIGHT and track node-mapping e[k(i)]--f-->x[i] between T(e0,e1,e2,e3),T(x[k0],..., x[k3]))
  inline void swap () {
    const short tmp0 = _node[0][0];
    const short tmp1 = _node[0][1];

    _node[0][0] = _node[1][0];
    _node[0][1] = _node[1][1];

    _node[1][0] = tmp0;
    _node[1][1] = tmp1;
  }

  //! Return the local coordinates of a node in a grid with width h=1
  template< typename RealType >
  void computeLocalCoordinate ( aol::Vec3<RealType> &tn0x, const tpcfe::CFEElement<RealType> &el, const int index ) const;

  //! Return the local coordinates of a node in a grid with width h=1
  template< typename RealType >
  void computeLocalCoordinate ( aol::Vec3<RealType> &tn0x, RealType _cutRelation[8][8], const int index ) const;

  /** Return the global coordinates of a node in a grid with width h=1
   */
  template< typename RealType >
  void computeGlobalCoordinate ( aol::Vec3<RealType> &tn0x, const tpcfe::CFEElement<RealType> &el, const int index ) const;

  //! Return the global coordinates of the barycenter of this tetrahedron in a grid with width h=1
  template< typename RealType >
  void computeGlobalCoordinateOfBarycenter ( aol::Vec3<RealType> &bary, const tpcfe::CFEElement<RealType> &el ) const {
    bary.setZero();

    aol::Vec3<RealType> baryc;
    for ( short i = 0; i < 4; ++i ) {
      computeGlobalCoordinate ( baryc, el, i );
      baryc *= static_cast<RealType> ( 0.25 );
      bary += baryc;
    }
  }

};


/** CFE Tetrahedra including geometric information
 *
 *  \author Preusser, Schwen
 */
template <typename RealType>
class CFETetra  : public CFETopoTetra {

public:
  aol::Vec3<RealType>              _edge[3];                  //!< three vectors spanning the tetrahedron
  mutable aol::Mat<4, 4, RealType> _lsm;                      //!< the local standard stiffness matrix produced by this tetrahedron
  mutable aol::Matrix33<RealType>  _barycentricCoordinates;   //!< barycentric Coords for the computation of local matrices
  // mutable is necessary here due to "lazy evaluation"

protected:
  RealType  _volume;          //!< stores volume

public:
  //! Standard constructor which creates a dummy tetrahedron. The values set here must not be used since they do not make sense.
  CFETetra () : CFETopoTetra(), _volume ( tpcfe::CFELookup<RealType>::_stdTetraVolume ) {
    // _edge, _lsm, _barycentricCoordinates initialized automatically
  }

  //! Create a tetrahedron with the given data
  CFETetra ( const short x0, const short y0,
             const short x1, const short y1,
             const short x2, const short y2,
             const short x3, const short y3,
             const short Parent, const VolType volType, const short virtualNodeNum, const signed char sign ) : CFETopoTetra ( x0, y0, x1, y1, x2, y2, x3, y3, Parent, volType, virtualNodeNum, sign ), _volume ( tpcfe::CFELookup<RealType>::_stdTetraVolume ) {
  }

  //! conversion constructor from basis class copying its data
  explicit CFETetra ( const CFETopoTetra &topoTet ) : CFETopoTetra ( topoTet ), _volume ( tpcfe::CFELookup<RealType>::_stdTetraVolume ) {
  }

  //! Copy constructor
  CFETetra ( const CFETetra<RealType> & other ) : CFETopoTetra ( other ), _lsm ( other._lsm ), _barycentricCoordinates ( other._barycentricCoordinates ), _volume ( other._volume ) {
    for ( short i = 0; i < 3; ++i ) {
      _edge[i] = other._edge[i];
    }
  }

  //! Assignment operator
  CFETetra<RealType> & operator= ( const CFETetra<RealType> & other ) {
    // Beware of self-assignment
    if ( this == &other )
      return *this;

    CFETopoTetra::operator= ( other );

    for ( short i = 0; i < 3; ++i )
      _edge[i] = other._edge[i];
    _lsm = other._lsm;
    _barycentricCoordinates = other._barycentricCoordinates;
    _volume = other._volume;

    return *this;
  }

  //! Return the volume of this tetra
  inline RealType getVolume() const { return _volume; }

  //! Set the volume of this tetra
  inline void setVolume ( const RealType Volume ) { _volume = Volume; }

  //! returns this tetrahedron as a string
  void print ( ostream &out ) const {
    for ( int i = 0; i < 4; ++i ) {
      out << "(" << _node[i][0] << ", " << _node[i][1] << ") ";
    }
    out << "sign " << static_cast<int>( _sign ) <<  " VT " << _volType << " par: " << _parent;
  }

  /** returns the determinant of the matrix spanned by the vectors \f$ v_1-v_0,\dots, v_3-v_0\f$
   *  N.B. The determinant equals 6 times the volume of this tetrahedron
   */
  RealType determinant() const;

  /** Invert the transformation matrix in order to obtain the barycentric coordinates
   *  of the given tetrahedron. Given the matrix
   *  \f$ K = [ v_1 - v_0; v_2 - v_0; v_3 - v_0 ] \f$ this method computes \f$ det(K)K^{-1} \f$,
   *  where the vectors \f$ v_0,\dots, v_3\f$ span this tetrahedron.
   */
  void computeInverseTransformation() const;

  /** compute the scaled stiffness matrix out of the inverted barycentric coordinates
   *  The stiffness-matrix is given by \f$ (D\lambda)(D\lambda)^T \f$, but this method computes
   *  is scaled with a factor \f$ det(K)^2 \f$. This must be remembered when the local matrix
   *  is assembled into the global one. The 4x3 matrix of barycentric coordinates is given by
   *  \f$ D\lambda = \left(\begin{array}{c}-\sum_i (K^{-1})_{ij}\\[.5ex]\hline K^{-1}\end{array}\right)\f$.
   *  The transformation matrix \f$ K^{-1} \f$ is determined with computeInverseTransformation
   */
  void computeStiffnessMatrix() const;

  //! Compute the center of mass of this tetrahedron, relative to 0th vertex
  void computeCenterOfMass ( aol::Vec3<RealType> &com ) const {
    cerr << "computeCenterOfMass called" << endl;
    com = ( _edge[0] + _edge[1] + _edge[2] ) / 3;
  }


  //! destructor has nothing to do
  ~CFETetra () { }

  void dump( ostream& out = cout ) const;

};

template <typename RealType>
ostream & operator<< ( ostream &out, const CFETetra<RealType> &t ) {
  t.print ( out );
  return out;
}


/** Iterator iterating over (regular or virtual) tetrahedra for one topological splitting case, depending on type, can iterate over all, positive and negative tetrahedra.
 *  \author Preusser, Schwen
 */
class CFETopoTetraIterator {
protected:
  std::vector < CFETopoTetra >::const_iterator _startIt, _curIt, _endIt;
  signed char _sign;

public:
  // standard constructor not useful but necessary
  CFETopoTetraIterator ( ) : _sign ( 2 ) {
  }


  explicit CFETopoTetraIterator ( const CFEType type, const signed char sign = 0 ) : _sign ( sign ) {
    restart ( type, sign );
  }

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
  const CFETopoTetra& operator* () const {
    if ( _curIt == _endIt ) {
      throw aol::Exception ( "tpcfe::CFETopoTetraIterator::operator*: cannot dereference *endIt", __FILE__, __LINE__ );
    }
    return ( *_curIt );
  }

  const CFETopoTetra* operator-> () const {
    if ( _curIt == _endIt ) {
      throw aol::Exception ( "tpcfe::CFETopoTetraIterator::operator->: cannot dereference *endIt", __FILE__, __LINE__ );
    }
    return ( &(*_curIt) );
  }

  /** restart existing iterator
   */
  void restart ( const CFEType type, const signed char sign = 0 ) {
    _sign = sign;
    _startIt = tpcfe::CFETopoLookup::_topo[type._pureType]->getTopoTetraVec().begin();
    _curIt   = _startIt;
    _endIt   = tpcfe::CFETopoLookup::_topo[type._pureType]->getTopoTetraVec().end();
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

  // copy constructor and assignment operator do the correct thing, but are probably slow
  CFETopoTetraIterator( const CFETopoTetraIterator& );
  CFETopoTetraIterator& operator= ( const CFETopoTetraIterator& );
};

}

#endif
