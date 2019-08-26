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

#ifndef __TFEBASEFUNCTIONSET_H
#define __TFEBASEFUNCTIONSET_H

#include <aol.h>
#include <quoc.h>

#include <tpCFELookup.h>
#include <tpCFETetrahedron.h>

namespace qc {

// --------------------------------------------------------------------------

template <typename RealType, qc::Dimension Dim>
class MidpointQuadrature {};

// --------------------------------------------------------------------------

template <typename RealType>
class MidpointQuadrature<RealType, qc::QC_3D> {
  static const RealType _barycentricPoint[4];
public:
  enum { numQuadPoints = 1 };
  inline static const RealType* getRefCoord ( int /*QuadPoint*/ );
  inline static RealType getWeight ( int /*QuadPoint*/ );
};

// --------------------------------------------------------------------------

// definition of static members
template <typename RealType>
const RealType MidpointQuadrature<RealType, qc::QC_3D>::_barycentricPoint[4] = { 0.25, 0.25, 0.25, 0.25 };

// class functions implementation
template <typename RealType>
inline const RealType* MidpointQuadrature<RealType, qc::QC_3D>::getRefCoord ( int /*QuadPoint*/ ) {
  return _barycentricPoint;
}

template <typename RealType>
inline RealType MidpointQuadrature<RealType, qc::QC_3D>::getWeight ( int /*QuadPoint*/ ) {
  return 1.;
}

// ==========================================================================

// forward declaration
template <typename RealType, class QuadRuleType>
class BaseFunctionSetTFE;

// --------------------------------------------------------------------------

/**
 *  The base function set for affine linear elements in 3d intersected
 *  by the boundary of a narrow band (TFE = truncated finite elements).
 *
 *  This class is for internal use via BaseFunctionSetTFE only. Do
 *  never instantiate it explicitely.
 */
template <typename RealType, class QuadRuleType>
class BaseFunctionSetTFE_Static {
  friend class BaseFunctionSetTFE<RealType, QuadRuleType>;
  friend class auto_ptr<BaseFunctionSetTFE_Static<RealType, QuadRuleType> >;

private:
  BaseFunctionSetTFE_Static ( const aol::Vector<RealType> & levelSetValues,
                              RealType radius,
                              RealType h );
  ~BaseFunctionSetTFE_Static ();

public:
  void printStats();

  int maxNumQuadPoints();
  int numQuadPoints ( int element ) const;
  inline RealType getWeight ( int element, int QuadPoint ) const;
  void getRefCoord ( int element, int QuadPoint, aol::Vec3<RealType>&RefCoord ) const;

  void evaluateGradient ( int BaseFuncNum, const aol::Vec3<RealType> &RefCoord, aol::Vec3<RealType> &Gradient ) const;
  inline const aol::Vec3<RealType>& evaluateGradient ( int element, int BaseFuncNum, int QuadPoint ) const;
  RealType evaluate ( int BaseFuncNum, const aol::Vec3<RealType> &RefCoord ) const;
  inline const RealType evaluate ( int element, int BaseFuncNum, int QuadPoint ) const;
  RealType getVolume ( int element ) const;

  bool nodeIsUsed ( int Element, int NodeNumber ) const;
  void initializeElement ( int Element, const int ElementToNodeOffsets[] );

private:
  typedef tpcfe::CFELookup<RealType> Lookup;

  //! Find out which nodes are actually used for computation.
  //! This algorithm iterates over all tetrahedra.
  unsigned char findUsedNodes ( const vector<tpcfe::CFETopoTetra> &tv );

  //! Find out which nodes are actually used for computation.
  //! This algorithm uses the cube element's signature and iteration
  //! over all edges and should give the same result as
  //! findUsedNodes(). It is only used in debug mode for
  //! testing coherence of tetrahedra vector and edge list.
  unsigned char findUsedNodesViaSignature ( unsigned short cfeType );
  unsigned short findUsedBit ( unsigned short n ) const;

  //! compute the cut relations between the interface and the tetrahedra
  void computeCutRelations( );

  //! compute vertices and volumes of sub-tetrahedra
  void computeVerticesAndVolumes ( vector<tpcfe::CFETetra<RealType> > &tv ) const;
  //! compute the weights of tetrahedra
  void computeWeights ( int element, const vector<tpcfe::CFETetra<RealType> > &tv );
  void computeQuadPointToTetraMap ( int element, const vector<tpcfe::CFETetra<RealType> > &tv );
  //! Compute the barycentric coordinates of all the quadpoints in the basis of the regular
  //! parent tetrahedron
  void computeRegTetraQuadPoints ( int element, const vector<tpcfe::CFETetra<RealType> > &tv );
  void initializeTetraCache();
  void initializeQuadCache();
  void initializeTetraVertex ();

  // upper bound (which is pessimistic) on maxNumQuadPoints determined by
  // #(regular tetrahedra) * #(max virtual tetrahedra inside
  // a regular tetrahedron and inside the narrow band) * #quadpoints
  // 3 comes from splitting a prism into 3 virtual tetrahedra
  static const int _maxNumQuadPoints = 6 * 3 * QuadRuleType::numQuadPoints;

  //! contains the tables that have to be cached + a bool that specifies
  //! whether the information for an element is already computed.
  struct tableCache {
    //! are the informations for the element already computed?
    bool isAlreadyComputed;

    //! the fractions of the volumina of the virtual tetrahedra divided by
    //! the volume of the regular tetrahedra
    RealType _weightTable[ BaseFunctionSetTFE_Static::_maxNumQuadPoints ];

    //! mapping from quadrature point to regular tetrahedron
    //! (used in evaluateGradient() and evaluate() methods)
    int _quadPointToTetra[ BaseFunctionSetTFE_Static::_maxNumQuadPoints ];

    //! mapping from a quadpoint number to the quadpoint's
    //! barycentric coordinates in the regular tetrahedron
    RealType _regTetraQuadpoint[ BaseFunctionSetTFE_Static::_maxNumQuadPoints ][ 4 ];

    //! volume of the used part of this element
    //! (sum of tetrahedra's parts' volume that lie inside the band)
    //! \note volume of an uncut element is 1, not \f$ h^3 \f$!
    RealType _volume;

    //! bit array (of length 8), where the i-th bit says if
    //! element node number i is a computational node.
    unsigned char nodeIsUsed;

    //! number of quadpoints = overall number of virtual tetrahedra inside nb
    //! times QuadRuleType::NumQuadPoints
    int _numQuadPoints;

    // constructor
    tableCache ()
        : isAlreadyComputed  ( false ) { }
  };

  //! level set values on the vertices.
  //! Only needed while initializeElement() is being performed.
  RealType _structureValue[8];

  //! the cut relations of the interface with the tetrahedra.
  //! Only needed while initializeElement() is being performed.
  RealType _cutRelation[8][8];

  //! distance from which an interface will be moved to a node
  RealType _adjustLevelset;
  RealType _adjustVN;

  //! the level function (surprise)
  RealType * _levelSetValues;

  // needs hopefully less memory than a stl::vector
  tableCache * _tableCacheVector;

  //! mapping that provides the vertex number of a certain regular
  //! tetrahedron to a cube vertex number.
  //! -1 means that the cube vertex doesn't belong to the tetrahedron.
  //! first arg: cube vertex, second: tetrahedron index
  int _stdCubeVertexTetraVertex[8][6];

  //! The 6 regular parent tetrahedra and their transformations cached.
  tpcfe::CFETetra<RealType> regularTetra[6];

  //! gradients of the 8 basis functions living on a real,
  //! sub-triangulated by 6 tetrahedra. The table is initialized by the
  //! constructor by scaling stdGradTable with _h.
  aol::Vec3<RealType> basisTetraGradients[8][6];

  RealType    _h;
  RealType    _narrowBandRadius;

  mutable int _numNodes;
  mutable int _numCalls;
  mutable int _numCallsWithLocalCache;
  mutable int _numCallsWithGlobalCache;
};

// --------------------------------------------------------------------------

/**
 *  Base function set for FE functions that are truncated to a
 *  narrow band around the zero level set of a given level
 *  function.
 *
 *  Assumptions on input data:
 *  - problem dimension is 3D
 *  - QuadRuleType is a quadrature rule on simplices
 *
 *  The base functions are just as usual on points where
 *  levelSetFct's absolute value is smaller than narrowBandRadius,
 *  gut zero anywhere else. During the setToElement() execution,
 *  enough quadrature points are computed such that the
 *  quadrature rule is just as exact as it would be on nun-truncated
 *  functions.
 *
 *  You may use this class also for full grid problems, but of
 *  course there is no use in doing so. Mostly you will use
 *  it in grids where only elements are stored which cut the
 *  narrow band, e. g. the nb::NarrowBandGrid or dt::DTGrid classes.
 *
 *  This class has nearly no own functionality, it is only a
 *  wrapper for the BaseFunctionSetTFE_Static. It redirects
 *  all work to be done to the static BaseFunctionSetTFE_Static variable
 *  which is (thanks to auto_ptr) released at the end of the program.
 *
 *
 *  BaseFunctionSetTFE_Static has many disadvantages -- perhaps the
 *  the largest one is that the bfs of any element has to be
 *  manually initialized by calling  initializeElement(). This is
 *  done in this class when you call setToElement(). Any configurator
 *  that uses this class should do so when his getBaseFuncSet() method
 *  is called.
 *
 *  \author Nielsen, Nemitz, von Deylen
 */
template <typename RealType, class QuadRuleType>
class BaseFunctionSetTFE {

public:
  BaseFunctionSetTFE ();

  BaseFunctionSetTFE ( const aol::Vector<RealType> & levelSetFct,
                       RealType radius, RealType h );

  void destroy();
  void setToElement ( int Element, const int ElementToNodeOffsets[] );
  int maxNumQuadPoints() const;
  int numQuadPoints( ) const;
  RealType getWeight ( int QuadPoint ) const;
  const aol::Vec3<RealType> getRefCoord ( int QuadPoint ) const;
  void evaluateGradient ( int /* BaseFuncNum */,
                          const aol::Vec3<RealType> & /* refCoord */,
                          aol::Vec3<RealType> & /* gradient */ ) const;
  const aol::Vec3<RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const;
  const aol::Vec3<RealType>& evaluateGradient ( int BaseFuncNum, int TetraNum,
                                                const aol::BarCoord<3, RealType> & /* refCoord */ ) const;
  RealType evaluate ( int, const aol::Vec3<RealType> & ) const;
  RealType evaluate ( int BaseFuncNum, int QuadPoint ) const;
  RealType evaluate ( int BaseFuncNum, int TetraNum, const aol::BarCoord<3, RealType> & refCoord ) const;
  bool nodeIsUsed ( int NodeNumber ) const;
  RealType getVolume () const;

protected:
  int _element;

private:
  static auto_ptr<BaseFunctionSetTFE_Static<RealType, QuadRuleType> > bfs_static;
};

// --------------------------------------------------------------------------

template <typename RealType, class QuadRuleType>
auto_ptr<BaseFunctionSetTFE_Static<RealType, QuadRuleType> >
BaseFunctionSetTFE<RealType, QuadRuleType>::bfs_static;

/* ==========================================================================
 *
 *      IMPLEMENTATION BaseFunctionSetTFE_Static
 *
 */

template <typename RealType, typename QuadRuleType>
BaseFunctionSetTFE_Static<RealType, QuadRuleType>::BaseFunctionSetTFE_Static ( const aol::Vector<RealType> & levelSetValues,
                                                                               RealType radius,
                                                                               RealType h )
    : _adjustLevelset ( 0.05 * h )
    , _adjustVN ( 0.05 )
    , _levelSetValues ( levelSetValues.getData() )
    , _tableCacheVector ( new tableCache[levelSetValues.size() ] )
    , _h ( h )
    , _narrowBandRadius ( radius )
    , _numNodes ( levelSetValues.size() )
    , _numCalls ( 0 )
    , _numCallsWithLocalCache ( 0 )
    , _numCallsWithGlobalCache ( 0 ) {

  Lookup::construct();
  initializeQuadCache();
  initializeTetraCache();
  initializeTetraVertex();
}

// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
BaseFunctionSetTFE_Static<RealType, QuadRuleType>::~BaseFunctionSetTFE_Static ( ) {
  /* printStats(); */
  delete[] ( _tableCacheVector );
  Lookup::destruct();
 }

// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::printStats() {
  cout << "==================================" << endl
       << "CFENBBaseFunctionSet Static Stats:" << endl
       << "numNodes = " << _numNodes << endl
       << "numCalls = " << _numCalls << endl
       << "numCalls (with local cache) = " << _numCallsWithLocalCache << endl
       << "number of calls removed by local cache = "
       << _numCalls - _numCallsWithLocalCache << endl
       << "numCalls (with global cache) = "
       << _numCallsWithGlobalCache << endl
       << "average calls pr node = "
       << 1. * ( _numCalls ) / _numNodes << endl
       << "average calls pr node (with local cache) = "
       << 1. * ( _numCallsWithLocalCache ) / _numNodes << endl
       << "average calls pr node (with global cache) = "
       << 1. * ( _numCallsWithGlobalCache ) / _numNodes << endl
       << "SizeOf( Cache per element ) = " << sizeof ( tableCache ) << endl
       << "==================================" << endl;
}

// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
int BaseFunctionSetTFE_Static<RealType, QuadRuleType>::maxNumQuadPoints()                          {   return _maxNumQuadPoints;   }
// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
int BaseFunctionSetTFE_Static<RealType, QuadRuleType>::numQuadPoints ( int element ) const {   return _tableCacheVector[element]._numQuadPoints;   }
// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
inline RealType BaseFunctionSetTFE_Static<RealType, QuadRuleType>::getWeight ( int element, int QuadPoint ) const {   return _tableCacheVector[element]._weightTable[QuadPoint];   }
// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::getRefCoord ( int element, int QuadPoint,
                                                                      aol::Vec3<RealType>& RefCoord ) const {

  // t is the number of the regular tetrahedron which contains the
  // given quad point:
  short t = _tableCacheVector[element]._quadPointToTetra[QuadPoint];

  // Lookup::_stdTetraVertex is a short*, containing the (local)
  // indices of t's nodes.
  aol::Vec<4, short> v ( Lookup::_stdTetraVertex[t] );

  // _tableCache::_regTetraQuadpoint gives the barycentric coords
  // of the quad point wrt the surrounding regular tetrahedron.
  aol::Vec<4, RealType> barCoord ( _tableCacheVector[element]._regTetraQuadpoint[QuadPoint] );
  for ( int i = 0; i < 3; i++ ) {
    RefCoord[i] = 0;
    for ( int j = 0; j < 4; ++j )
      // Lookup::_map stores the reference coordinates of
      // the cube nodes wrt this cube.
      RefCoord[i] += Lookup::_hexNbOffs[v[j]][i] * barCoord[j];
  }
}
// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::evaluateGradient ( int /*BaseFuncNum*/,
                                                                           const aol::Vec3<RealType> &/*RefCoord*/,
                                                                           aol::Vec3<RealType> &/*Gradient*/ ) const {
  throw aol::UnimplementedCodeException ( "evaluateGradient with RefCoord: "
                                          "Not implemented yet!", __FILE__, __LINE__ );
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
inline const aol::Vec3<RealType>&
BaseFunctionSetTFE_Static<RealType, QuadRuleType>::evaluateGradient ( int element, int BaseFuncNum, int QuadPoint ) const {
  return basisTetraGradients[BaseFuncNum][_tableCacheVector[element]
                                          ._quadPointToTetra[QuadPoint]];
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
RealType BaseFunctionSetTFE_Static<RealType, QuadRuleType>::evaluate ( int /*BaseFuncNum*/, const aol::Vec3<RealType> &/*RefCoord*/ ) const {
  throw aol::Exception ( "evaluate with RefCoord: Not implemented yet!",
                         __FILE__, __LINE__ );
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
inline const RealType BaseFunctionSetTFE_Static<RealType, QuadRuleType>::evaluate ( int element, int BaseFuncNum, int QuadPoint ) const {

  // Find out in which regular tetrahedron our quadpoint lies
  int regularTetrahedron
  = _tableCacheVector[element]._quadPointToTetra[QuadPoint];

  // Given a regularTetrahedron, and an index of the cube vertex
  // (baseFuncNum), we wish to find out which vertex of the tetrahedron
  // this corresponds to (regularTetrahedronVertex).
  int regularTetrahedronVertex
  = _stdCubeVertexTetraVertex[BaseFuncNum][regularTetrahedron];

  // If regularTetrahedronVertex is -1, this cube vertex does not exis
  // in this tetrahedron, and hence the baseFunction of this vertex is
  // zero at the quadpoint, otherwise the baseFunction at the quadpoint
  // is the regularTetrahedronVertex'th value of the quadpoint in local
  // coordinates of the regular tetrahedron
  return ( regularTetrahedronVertex == -1 )
         ? 0.0
         : _tableCacheVector[element].
         _regTetraQuadpoint[QuadPoint][regularTetrahedronVertex];
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
RealType BaseFunctionSetTFE_Static<RealType, QuadRuleType>::getVolume ( int element ) const {
  return _tableCacheVector[element]._volume;
}
// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
inline bool BaseFunctionSetTFE_Static<RealType, QuadRuleType>::nodeIsUsed ( int Element, int NodeNumber ) const {   return _tableCacheVector[Element].nodeIsUsed & ( 1 << NodeNumber );       }
// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::initializeElement ( int Element, const int ElementToNodeOffsets[] ) {
  _numCallsWithLocalCache++;

  if ( _tableCacheVector[Element].isAlreadyComputed ) return;

  _numCallsWithGlobalCache++;

  // Get the level set values at the vertices
  for ( int i = 0; i < 8; i++ ) {
    RealType p = fabs ( _levelSetValues[Element + ElementToNodeOffsets[i] ] )
                 - _narrowBandRadius;
    // TODO: Would this adjustment be necessary?
    /*
            if ( fabs ( p ) < _adjustLevelset ) {
              // Level-set function is too close to zero at grid point
              if ( p < 0 ) p -= _adjustLevelset;
              else p += _adjustLevelset;
            }
    */
    _structureValue[i] = p;
  }

  // compute the signature of the element (in [0,255]) stored in cfeType
  unsigned short cfeType = 0x0;
  unsigned char  ONE     = 0x1;

  // Run through all nodes of element
  for ( int i = 0; i < 8; i++ ) {
    if ( _structureValue[i] < 0 ) cfeType |= ONE;
    ONE <<= 1;
  }
  /*
    if ( !cfeType )
      throw aol::Exception ( "Unused element inserted into narrow band!",
                             __FILE__, __LINE__ );
  */
  // this gives us a vector of tetrahedra depending on the topology
  // in case of no intersection, this vector contains the regular tetrahedra,
  // if there's intersection, it contains the regular tetrahedra
  // which are completely contained in the narrow band
  // and the virtual tetrahedra inside.
  // To get only inner tetrahedra, we use -1 in the iterator
  vector < tpcfe::CFETopoTetra > topotv;

  for ( tpcfe::CFETopoTetraIterator tit ( tpcfe::CFEType ( -1/*domain*/, cfeType ), -1 ); tit.notAtEnd(); ++tit ) {
    topotv.push_back ( *tit );
  }

  // find out which nodes of the element are computational nodes.
  _tableCacheVector[Element].nodeIsUsed = findUsedNodes ( topotv );
#ifdef DEBUG
  if ( findUsedNodesViaSignature ( cfeType )
       != _tableCacheVector[Element].nodeIsUsed )
    throw aol::Exception ( "findUsedNodes() and findUsedNodesViaSignature() "
                           "do not return the same used masks.", __FILE__, __LINE__ );
#endif


  // set number of quadpoints. This is equal to the number
  // of tetrahedra in the tetrahedra list
  _tableCacheVector[Element]._numQuadPoints = static_cast<int> ( topotv.size() );

  std::vector< tpcfe::CFETetra<RealType> > tv;
  tv.reserve( topotv.size() );
  for ( typename::vector< tpcfe::CFETopoTetra >::const_iterator it = topotv.begin(); it != topotv.end(); ++it ) {
    tv.push_back ( tpcfe::CFETetra<RealType> ( *it ) );
  }

  // Compute cut relations on edges
  computeCutRelations();

  // Compute vertices and volumes of sub-tetrahedra
  computeVerticesAndVolumes ( tv );

  // compute weights
  computeWeights ( Element, tv );

  // compute quadpoint to regular tetrahedron mapping
  computeQuadPointToTetraMap ( Element, tv );

  // transform the quadpoints in barycentric coordinates
  // of the cut tetrahedra to barycentric coordinates
  // of the regular tetrahedra
  computeRegTetraQuadPoints ( Element, tv );

  _tableCacheVector[Element]._volume = 0.;
  typename vector<tpcfe::CFETetra<RealType> >::const_iterator it;
  for ( it = tv.begin(); it != tv.end(); ++it )
    _tableCacheVector[Element]._volume += it->getVolume();

  for ( int i = 0; i < _tableCacheVector[Element]._numQuadPoints; ++i )
    _tableCacheVector[Element]._weightTable[i] /= _tableCacheVector[Element]._volume;

  _tableCacheVector[Element].isAlreadyComputed = true;
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
unsigned char BaseFunctionSetTFE_Static<RealType, QuadRuleType>::findUsedNodes ( const vector< tpcfe::CFETopoTetra > &tv ) {
  // initalize nodeIsInBand vector as false
  unsigned char nodeUsed = 0;

  // iterate over all standard and virtual tetrahedra
  // entirely contained in the computation domain
  for ( unsigned short i_tetra = 0; i_tetra < tv.size(); ++i_tetra ) {
    // iterate over all nodes of one tetrahedron
    for ( unsigned short i_node = 0; i_node < 4; ++i_node ) {
      nodeUsed |= ( 1 << ( tv[i_tetra] ) ( i_node, 0 ) );
      // If node[i_node][1] is not NODE_NOT_VIRTUAL, the stored
      // tetrahedron node is actually not the cube node with number
      // node[i_node][0], but a virtual node between node[i_node][0]
      // and node[i_node][1]. Thus, also node[i_node][1] is
      // a cube node which is used for computation.
      if ( ( tv[i_tetra] ) ( i_node, 1 ) != tpcfe::NODE_NOT_VIRTUAL )
        nodeUsed |= ( 1 << ( tv[i_tetra] ) ( i_node, 1 ) );
    }
  }
  return nodeUsed;
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
unsigned char BaseFunctionSetTFE_Static<RealType, QuadRuleType>::findUsedNodesViaSignature ( unsigned short cfeType ) {
  // initalize nodeIsInBand vector as false
  unsigned char nodeUsed = 0;

  // iterate over all nodes of the cube element
  for ( unsigned short i_node = 0; i_node < 8; ++i_node )
    // if i-th bit of the signature is set, the node lies
    // inside the band. (1 << i) is a one shifted i bits
    // to the left.
    if ( cfeType & ( 1 << i_node ) ) {
      nodeUsed |= ( 1 << i_node );
      for ( unsigned short i_edge = 0; i_edge < 19; ++i_edge ) {
        // edge[i] has two bits set, one for
        // each end point. If one of them is
        // the i-th node, this node is entirely
        // surrounded by regular and virtual tetrahedra,
        // and so every edge-connected neighbour will also
        // be a computational node.
        if ( Lookup::_edge[i_edge] & ( 1 << i_node ) ) {
          size_t second_node = findUsedBit (
                                 static_cast<unsigned short> (
                                   Lookup::_edge[i_edge] - ( 1 << i_node ) ) );
          nodeUsed |= ( 1 << second_node );
        }
      }
    }
  return nodeUsed;
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
unsigned short BaseFunctionSetTFE_Static<RealType, QuadRuleType>::findUsedBit ( unsigned short n ) const {
  for ( unsigned short i = 0; i < 8 * sizeof ( unsigned short ); ++i )
    if ( n & ( 1 << i ) )
      return i;
  return static_cast<unsigned short> ( -1 );
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::computeCutRelations( ) {
  // traverse possible vertex-pairings
  for ( int i = 0; i < 8; i++ ) {
    for ( int j = i; j < 8; j++ ) {
      const RealType a = _structureValue[i];
      const RealType b = _structureValue[j];

      const RealType denom = a - b;
      if ( denom ) {
        RealType tmp = a / denom;
        /*
        // TODO: Would this adjustment be necessary?
        // _cutRelations are computed on all edges of interfaced elements,
        // even if the edge is not interfaced, even if the edge
        // is not a tetrahedron edge at all.

        // We only adjust the _cutRelations on interfaced edges

        if( tmp >= 0.0 && tmp < 0.0 + _adjustVN )
          tmp = 0.0 + _adjustVN;
        if( tmp <= 1.0 && tmp > 1.0 - _adjustVN )
          tmp = 1.0 - _adjustVN;
        */
        _cutRelation[i][j] = tmp;
        _cutRelation[j][i] = 1.0 - tmp;
      } else
        _cutRelation[i][j] = _cutRelation[j][i] = 0.0;
    }
  }
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::computeVerticesAndVolumes ( vector<tpcfe::CFETetra<RealType> > &tv ) const {
  // iterate through all sub-tetrahedra
  // (ie. both those regular that are left and the virtual ones too)
  for ( typename vector<tpcfe::CFETetra<RealType> >::iterator it = tv.begin(); it != tv.end(); ++it ) {
    tpcfe::CFETetra<RealType> &t = *it;

    aol::Vec3<RealType> vertex[4];
    for ( int l = 0; l < 4; l++ ) {
      const short a = t ( l, 0 );
      const short b = t ( l, 1 );
      if ( b == tpcfe::NODE_NOT_VIRTUAL ) {
        vertex[l][0] = Lookup::_hexNbOffs[a][0];
        vertex[l][1] = Lookup::_hexNbOffs[a][1];
        vertex[l][2] = Lookup::_hexNbOffs[a][2];
      } else {        // Is a virtual node and thus an interpolation between two real nodes
        RealType w1 = _cutRelation[a][b];
        RealType w2 = _cutRelation[b][a];

        vertex[l][0] = w2 * Lookup::_hexNbOffs[a][0] + w1 * Lookup::_hexNbOffs[b][0];
        vertex[l][1] = w2 * Lookup::_hexNbOffs[a][1] + w1 * Lookup::_hexNbOffs[b][1];
        vertex[l][2] = w2 * Lookup::_hexNbOffs[a][2] + w1 * Lookup::_hexNbOffs[b][2];
      }
    }
    t._edge[0] = vertex[1] - vertex[0];
    t._edge[1] = vertex[2] - vertex[0];
    t._edge[2] = vertex[3] - vertex[0];

    t.setVolume ( t.determinant() / 6.0 );
  }
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::computeWeights ( int element, const vector<tpcfe::CFETetra<RealType> > &tv ) {
  typename vector<tpcfe::CFETetra<RealType> >::const_iterator it;
  int j, i;
  for ( i = 0, it = tv.begin(); it != tv.end(); ++it )
    for ( j = 0; j < QuadRuleType::numQuadPoints; j++, i++ )
      // The weight of a quad point is computed as the volume fraction times the weight
      // of the quadpoint as specified from the QuadRuleType
      _tableCacheVector[element]._weightTable[i] =
        it->getVolume() * QuadRuleType::getWeight ( j );
}


// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::computeQuadPointToTetraMap ( int element, const vector<tpcfe::CFETetra<RealType> > &tv ) {
  typename vector<tpcfe::CFETetra<RealType> >::const_iterator it;
  int i, j;
  for ( it = tv.begin(), i = 0; it != tv.end(); ++it )
    for ( j = 0; j < QuadRuleType::numQuadPoints; j++, i++ )
      _tableCacheVector[element]._quadPointToTetra[i] = it->getParent();
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::computeRegTetraQuadPoints ( int element, const vector<tpcfe::CFETetra<RealType> > &tv ) {
  // transform the quadpoints in barycentric coordinates of the cut tetrahedra
  // to barycentric coordinates of the regular tetrahedra
  int i, j;

  // iterate through the tetrahedra and the quadpoints in each tetrahedron
  typename vector<tpcfe::CFETetra<RealType> >::const_iterator it;
  for ( it = tv.begin(), i = 0; it != tv.end(); ++it ) {
    // Now determine in which regular tetrahedron this quad-point lies.
    // This is equal to the parent tetrahedron
    short regularTetrahedron = it->getParent();

    for ( j = 0; j < QuadRuleType::numQuadPoints; j++, i++ ) {
      // Get the coordinates of the quadpoint in barycentric coordinates of the
      // virtual tetrahedron
      const RealType *virtualRefCoord = QuadRuleType::getRefCoord ( j );

      // Now transform the local ref coord to ref coords in the cube (local element)
      // p is the coordinates of the quad point in the coordinates of the cube (local element)
      // and p can be computed as the sum of the local barycentric coordinates in the virtual tetrahedron times the
      // coordinates of the virtual tetrahedron in the coordinates of the cube (local element)
      aol::Vec3<RealType> p;
      // tnx0 is temporary storage for a vertex of the tetrahedron in the coordinates of the local cube
      aol::Vec3<RealType> tn0x;

      // TODO: OPT: If there are several quadpoints, this only have to be computed once
      // per virtual tetrahedron!
      it->computeLocalCoordinate ( tn0x, _cutRelation, 0 );
      p  = virtualRefCoord[0] * tn0x;
      it->computeLocalCoordinate ( tn0x, _cutRelation, 1 );
      p += virtualRefCoord[1] * tn0x;
      it->computeLocalCoordinate ( tn0x, _cutRelation, 2 );
      p += virtualRefCoord[2] * tn0x;
      it->computeLocalCoordinate ( tn0x, _cutRelation, 3 );
      p += virtualRefCoord[3] * tn0x;

      // Given the ref coord in (x,y,z) coordinates of the cube, that is p which we just computed, compute
      // the barycentric coordinates of the ref coord in the regular tetrahedron
      // First create the regular tetrahedron
      // Here we need the nodes of the regular tetrahedron in the ordering of the cube
      // As parent we give the regular tetrahedron itself (given by its global ordering)
      aol::Vec3<RealType> vertex0;
      vertex0[0] = Lookup::_hexNbOffs[ Lookup::_stdTetraVertex[regularTetrahedron][0] ][0];
      vertex0[1] = Lookup::_hexNbOffs[ Lookup::_stdTetraVertex[regularTetrahedron][0] ][1];
      vertex0[2] = Lookup::_hexNbOffs[ Lookup::_stdTetraVertex[regularTetrahedron][0] ][2];
      // we need to subtract from p the 0'th vertex of this regular tetrahedron
      aol::Vec3<RealType> regularRefCoord = regularTetra[regularTetrahedron]._barycentricCoordinates * ( p - vertex0 );
      _tableCacheVector[element]._regTetraQuadpoint[ i ][ 0 ] = 1.0 - regularRefCoord[0] - regularRefCoord[1] - regularRefCoord[2];
      _tableCacheVector[element]._regTetraQuadpoint[ i ][ 1 ] = regularRefCoord[0];
      _tableCacheVector[element]._regTetraQuadpoint[ i ][ 2 ] = regularRefCoord[1];
      _tableCacheVector[element]._regTetraQuadpoint[ i ][ 3 ] = regularRefCoord[2];
    }
  }

}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::initializeTetraCache() {
  int regularTetrahedron;
  for ( regularTetrahedron = 0; regularTetrahedron < 6; regularTetrahedron++ ) {
    short x0 = Lookup::_stdTetraVertex[regularTetrahedron][0];
    short x1 = Lookup::_stdTetraVertex[regularTetrahedron][1];
    short x2 = Lookup::_stdTetraVertex[regularTetrahedron][2];
    short x3 = Lookup::_stdTetraVertex[regularTetrahedron][3];
    regularTetra[regularTetrahedron].setNodes ( x0, tpcfe::NODE_NOT_VIRTUAL,
                                                x1, tpcfe::NODE_NOT_VIRTUAL,
                                                x2, tpcfe::NODE_NOT_VIRTUAL,
                                                x3, tpcfe::NODE_NOT_VIRTUAL );
    regularTetra[regularTetrahedron].setParent ( regularTetrahedron );

    // *** Initialize edges of the regular tetrahedra ***
    aol::Vec3<RealType> vertex[4];
    for ( int l = 0; l < 4; l++ ) {
      const short a = regularTetra[regularTetrahedron] ( l, 0 );
      vertex[l][0] = Lookup::_hexNbOffs[a][0];
      vertex[l][1] = Lookup::_hexNbOffs[a][1];
      vertex[l][2] = Lookup::_hexNbOffs[a][2];
    }
    regularTetra[regularTetrahedron]._edge[0] = vertex[1] - vertex[0];
    regularTetra[regularTetrahedron]._edge[1] = vertex[2] - vertex[0];
    regularTetra[regularTetrahedron]._edge[2] = vertex[3] - vertex[0];
    // ***************************************************

    // Compute the inverse transformation on the tetrahedron
    // enabling us to transform a vertex in (x,y,z) coordinates
    // of the cube to barycentric coordinates of this regular tetrahedron
    regularTetra[regularTetrahedron].computeInverseTransformation();
  }
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::initializeQuadCache() {
  // gradients of the 8 basis functions living on the std cube,
  // sub-triangulated in the standard way by 6 tetrahedra. The vertex
  // mapping is as defined by the "stdCubeVertices" and
  // "stdTetraVertexMap" tables below ("#if 0" part after private:
  // tag).
  //
  // The indexing is: grad_table[#tetra][#bfct][d], so it is a 6x8
  // lookup-table storing the gradients as vectors.
  //
  // QuadPoint serves as tetra number as long as we only use
  // p.w. constant gradients.
  //
  static const RealType _stdGradTable[6][8][3] = {
                                                   { // 0th tetra, vertices 0, 1, 2, 6
                                                     { -1.0, -1.0,  0.0 },
                                                     {  1.0,  0.0,  0.0 },
                                                     {  0.0,  1.0, -1.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  1.0 },
                                                     {  0.0,  0.0,  0.0 }
                                                   },
                                                   { // 1st tetra, vertices 1, 0, 5, 6
                                                     { -1.0, -1.0,  0.0 },
                                                     {  1.0,  1.0, -1.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0, -1.0,  1.0 },
                                                     {  0.0,  1.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 }
                                                   },
                                                   { // 2nd tetra, vertices 0, 4, 5, 6
                                                     {  0.0,  0.0, -1.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     { -1.0, -1.0,  1.0 },
                                                     {  1.0,  0.0,  0.0 },
                                                     {  0.0,  1.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 }
                                                   },
                                                   { // 3rd tetra, vertices 2, 1, 3, 7
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0, -1.0,  0.0 },
                                                     { -1.0,  0.0,  0.0 },
                                                     {  1.0,  1.0, -1.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  1.0 }
                                                   },
                                                   { // 4th tetra, vertices 1, 2, 6, 7
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0, -1.0,  0.0 },
                                                     {  0.0,  1.0, -1.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     { -1.0, -1.0,  1.0 },
                                                     {  1.0,  1.0,  0.0 }
                                                   },
                                                   { // 5th tetra, vertices 5, 1, 6, 7
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0, -1.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0,  0.0,  0.0 },
                                                     {  0.0, -1.0,  1.0 },
                                                     { -1.0,  0.0,  0.0 },
                                                     {  1.0,  1.0,  0.0 }
                                                   }
                                                 };
  for ( int b = 0; b < 8; b++ )
    for ( int q = 0; q < 6; q++ )
      for ( int d = 0; d < 3; d++ )
        basisTetraGradients[b][q][d] = _stdGradTable[q][b][d] / _h;
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE_Static<RealType, QuadRuleType>::initializeTetraVertex () {
  static const int stdCubeVertexTetraVertex[8][6] = {
                                                      {  0,  1,  0, -1, -1, -1 },    // 0
                                                      {  1,  0, -1,  1,  0,  1 },    // 1
                                                      {  2, -1, -1,  0,  1, -1 },    // 2
                                                      { -1, -1, -1,  2, -1, -1 },    // 3
                                                      { -1, -1,  1, -1, -1, -1 },    // 4
                                                      { -1,  2,  2, -1, -1,  0 },    // 5
                                                      {  3,  3,  3, -1,  2,  2 },    // 6
                                                      { -1, -1, -1,  3,  3,  3 }     // 7
                                                    };

  for ( int i = 0; i < 8; ++i )
    for ( int j = 0; j < 6; ++j )
      _stdCubeVertexTetraVertex[i][j] = stdCubeVertexTetraVertex[i][j];
}

/* ==========================================================================
 *
 *      IMPLEMENTATION BaseFunctionSetTFE
 *
 */

//! standard constructor. Works only if the static base
//! function set has already been constructed via a call
//! to the other constructor.
template <typename RealType, typename QuadRuleType>
BaseFunctionSetTFE<RealType, QuadRuleType>::BaseFunctionSetTFE ()
    : _element ( -1 ) {
  if ( !bfs_static.get() )
    throw aol::Exception ( "BaseFunctionSetTFE: First instance has to be "
                           "constructed with level set, band radius and "
                           "element width!", __FILE__, __LINE__ );
}

// --------------------------------------------------------------------------

//! first-time constructor. This one has to be called before
//! any call to the standard constructor in order that the
//! static base function set can be constructed.
template <typename RealType, typename QuadRuleType>
BaseFunctionSetTFE<RealType, QuadRuleType>::BaseFunctionSetTFE ( const aol::Vector<RealType> & levelSetFct,
                                                                 RealType radius, RealType h )
    : _element ( -1 ) {
  bfs_static.reset ( new BaseFunctionSetTFE_Static<RealType, QuadRuleType>
                     ( levelSetFct, radius, h ) );
}

// --------------------------------------------------------------------------

//! destroys the static base function set. Can be called
//! before terminating the program (to free memory before
//! deconstruction of the vector manager) or to enforce a
//! re-construction with the next constructor call (e. g.
//! when the band width has changed).
template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE<RealType, QuadRuleType>::destroy()                                         {   bfs_static.reset();   }
// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE<RealType, QuadRuleType>::setToElement ( int Element, const int ElementToNodeOffsets[] ) {
  bfs_static->_numCalls++;
  if ( Element == _element )
    return;

  bfs_static->initializeElement ( Element, ElementToNodeOffsets );
  _element = Element;
}

// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
int BaseFunctionSetTFE<RealType, QuadRuleType>::maxNumQuadPoints() const       {   return bfs_static->maxNumQuadPoints();   }
// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
int BaseFunctionSetTFE<RealType, QuadRuleType>::numQuadPoints( ) const {   return bfs_static->numQuadPoints ( _element );   }
// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
RealType BaseFunctionSetTFE<RealType, QuadRuleType>::getWeight ( int QuadPoint ) const {   return bfs_static->getWeight ( _element, QuadPoint );   }
// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
RealType BaseFunctionSetTFE<RealType, QuadRuleType>::getVolume () const         {   return bfs_static->getVolume ( _element );   }
// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
const aol::Vec3<RealType> BaseFunctionSetTFE<RealType, QuadRuleType>::getRefCoord ( int QuadPoint ) const {
  aol::Vec3<RealType> ret;
  bfs_static->getRefCoord ( _element, QuadPoint, ret );
  return ret;
}
// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFunctionSetTFE<RealType, QuadRuleType>::evaluateGradient ( int, const aol::Vec3<RealType> &,
                                                                    aol::Vec3<RealType> & ) const {
  throw aol::UnimplementedCodeException ( "evaluateGradient with RefCoord: "
                                          "Not implemented yet!", __FILE__, __LINE__ );
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
const aol::Vec3<RealType>& BaseFunctionSetTFE<RealType, QuadRuleType>::evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
  return bfs_static->evaluateGradient ( _element, BaseFuncNum, QuadPoint );
}

// --------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
RealType BaseFunctionSetTFE<RealType, QuadRuleType>::evaluate ( int, const aol::Vec3<RealType> & ) const {
  throw aol::UnimplementedCodeException ( "evaluate with RefCoord: "
                                          "Not implemented yet!", __FILE__, __LINE__ );
}

// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
RealType BaseFunctionSetTFE<RealType, QuadRuleType>::evaluate ( int BaseFuncNum, int QuadPoint ) const {   return bfs_static->evaluate ( _element, BaseFuncNum, QuadPoint );   }
// --------------------------------------------------------------------------
template <typename RealType, typename QuadRuleType>
bool BaseFunctionSetTFE<RealType, QuadRuleType>::nodeIsUsed ( int NodeNumber ) const {   return bfs_static->nodeIsUsed ( _element, NodeNumber );   }
// --------------------------------------------------------------------------

} // end of namespace qc.

#endif
