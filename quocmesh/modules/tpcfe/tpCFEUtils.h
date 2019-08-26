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

#ifndef __TPCFEUTILS_H
#define __TPCFEUTILS_H

#include <bzipiostream.h>
#include <scalarArray.h>
#include <triangMesh.h>
#include <tpCFEGrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFETetrahedron.h>

namespace tpcfe {

//! compute Lame-Navier constant lambda from Young's modulus and Poisson's ratio
template < typename RealType >
RealType computeLambda ( const RealType E, const RealType nu ) {
  return ( nu * E / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) ) );
}

//! compute Lame-Navier constant mu from Young's modulus and Poisson's ratio
template < typename RealType >
RealType computeMu ( const RealType E, const RealType nu ) {
  return ( E / ( 2.0 * ( 1.0 + nu ) ) );
}


template< typename RealType >
void setLameNavierTensorFromLambdaMu ( const RealType lambda, const RealType mu, RealType LNTensor[3][3][3][3] ) {
  for ( short i = 0; i < 3; ++i )
    for ( short j = 0; j < 3; ++j )
      for ( short k = 0; k < 3; ++k )
        for ( short l = 0; l < 3; ++l )
          LNTensor[i][j][k][l] = lambda * ( i == j ) * ( k == l ) + mu * ( ( i == k ) * ( j == l ) + ( i == l ) * ( j == k ) );
}


template< typename RealType >
void setLameNavierTensorFromENu ( const RealType E, const RealType nu, RealType LNTensor[3][3][3][3] /* This is automatically call-by-reference. */ ) {
  setLameNavierTensorFromLambdaMu ( computeLambda ( E, nu ), computeMu ( E, nu ), LNTensor );
}


template< typename GridType, typename MatrixType >
void setNotInnerDomainNodesToIdentity ( const GridType &grid, MatrixType &matrix, typename GridType::RealType diag_value = 1.0 ) {
  aol::ProgressBar<> pb ( "Setting nearBoundary entries to 1" );
  pb.start ( grid.getNumberOfNodes() );

  for ( int i = 0; i < grid.getNumberOfNodes(); ++i ) {
    if (  grid.isDomainNode ( i ) && !grid.isInnerNode ( i ) ) {
      matrix.setRowColToZero ( i );
      matrix.set ( i, i, diag_value );
    }
#ifdef VERBOSE
    pb++;
#endif
  }
#ifdef VERBOSE
  pb.finish();
#endif
}

template< typename GridType, typename MatrixType >
void setNonDomainEntriesToIdentity ( const GridType &grid, MatrixType &matrix, typename GridType::RealType diag_value = 1.0 ) {
  aol::ProgressBar<> pb ( "Setting nonDomain entries to 1" );
  pb.start ( grid.getNumberOfNodes() );

  for ( int i = 0; i < grid.getNumberOfNodes(); ++i ) {
    if ( ! ( grid.isDomainNode ( i ) ) ) {
      matrix.set ( i, i, diag_value );
    }
#ifdef VERBOSE
    pb++;
#endif
  }
#ifdef VERBOSE
  pb.finish();
#endif
}


template< typename GridType, typename MatrixType >
void restrictDirichletEntries ( const GridType &grid, MatrixType &matrix, typename GridType::RealType diag_value = 1.0 ) {
  aol::ProgressBar<> pb ( "Restricting Dirichlet entries" );
  pb.start ( grid.getNumberOfNodes() );

  for ( int i = 0; i < grid.getNumberOfNodes(); ++i ) {
    if ( grid.isDirichletNode ( i ) ) {
      matrix.setRowColToZero ( i );
      matrix.set ( i, i, diag_value );
    }
#ifdef VERBOSE
    pb++;
#endif
  }
#ifdef VERBOSE
  pb.finish();
#endif

}

template< typename GridType, typename MatrixType >
void restrictNonDomainEntries ( const GridType &grid, MatrixType &matrix, typename GridType::RealType diag_value = 1.0 ) {
  aol::ProgressBar<> pb ( "Restricting nonDomain entries" );
  pb.start ( grid.getNumberOfNodes() );

  for ( int i = 0; i < grid.getNumberOfNodes(); ++i ) {
    if ( ! ( grid.isDomainNode ( i ) ) ) {
      matrix.setRowColToZero ( i );
      matrix.set ( i, i, diag_value );
    }
#ifdef VERBOSE
    pb++;
#endif
  }
#ifdef VERBOSE
  pb.finish();
#endif
}

template< typename GridType, typename MatrixType >
void restrictNonDOFEntries ( const GridType &grid, MatrixType &matrix, typename GridType::RealType diag_value = 1.0 ) {
  aol::ProgressBar<> pb ( "Restricting nonDOF entries" );
  pb.start ( grid.getNumberOfNodes() );

  for ( int i = 0; i < grid.getNumberOfNodes(); ++i ) {
    if ( ! ( grid.isDOFNode ( i ) ) ) {
      matrix.setRowColToZero ( i );
      matrix.set ( i, i, diag_value );
    }
#ifdef VERBOSE
    pb++;
#endif
  }
#ifdef VERBOSE
  pb.finish();
#endif
}

template< typename GridType, typename MatrixType >
void restrictOuterDOFEntries ( const GridType &grid, MatrixType &matrix, typename GridType::RealType diag_value = 1.0 ) {
  aol::ProgressBar<> pb ( "Setting outerDOF entries to 1" );
  pb.start ( grid.getNumberOfNodes() );
  for ( int i = 0; i < grid.getNumberOfNodes(); ++i ) {
    if (  grid.isDomainNode ( i ) && ! ( grid.isInnerNode ( i ) ) ) {
      matrix.setRowColToZero ( i );
      matrix.set ( i, i, diag_value );
    }
#ifdef VERBOSE
    pb++;
#endif
  }
#ifdef VERBOSE
  pb.finish();
#endif
}

template< typename GridType >
void determineElementAndTetraForWC ( const GridType &grid, const aol::Vec3< typename GridType::RealType > &positionWC,
                                     CFEElement< typename GridType::RealType > &element, CFETopoTetra &tetra, aol::Vec< 4, typename GridType::RealType > &baryCO, bool &foundInside,
                                     const typename GridType::RealType tolerance = 1.0e-14 ) {
  typedef typename GridType::RealType RealType;

  aol::Vec3<RealType> position ( positionWC );           // global = image coordinates
  position *= static_cast<RealType> ( 1.0 / grid.H() );

  short x = static_cast<short> ( floor ( position[0] ) ), y = static_cast<short> ( floor ( position[1] ) ), z = static_cast<short> ( floor ( position[2] ) );
  const short n = static_cast<short> ( grid.getNumXYZ() );

  // width - 1 is last node but not an element
  if ( x == ( n - 1 ) )      x -= 1;
  if ( y == ( n - 1 ) )      y -= 1;
  if ( z == ( n - 1 ) )      z -= 1;

  const CFEType elType = grid.getElType ( x, y, z );

  element = CFEElement<RealType> ( grid.getSize(), x, y, z, elType ); // this is the element in which the point lies.
  if ( element.cfeType().representsInterfaced() ) // allow computeCutRelations always? ... why does this not work always??
    grid.getCutRelations ( element );

  short whichTetra = 0;

  switch ( GridType::CT ) {
    case CFE_CD:
      whichTetra = -1;
      break;
    case CFE_TPOS:
    case CFE_TPOSELAST:
    case CFE_LIEHR:
      whichTetra = 0;
      break;
    default:
      throw aol::Exception ( "tpcfe::determineElementAndTetraForWC: only implemented for CFE_CD, CFE_LIEHR, CFE_TPOS, CFE_TPOSELAST", __FILE__, __LINE__ );
  }

  for ( CFETopoTetraIterator tit ( elType, whichTetra ); tit.notAtEnd(); ++tit ) {
    aol::Vec3<RealType> vertices[4], shifted_position, bary3; // in global coordinates

    const CFETopoTetra &tet = *tit;

    for ( int i = 0; i < 4; ++i )
      tet.computeGlobalCoordinate ( vertices[i], element, i );

    aol::Matrix33<RealType> matinv;

    aol::Matrix33<RealType> mat;
    mat.setCol ( 0, vertices[0] - vertices[3] );  mat.setCol ( 1, vertices[1] - vertices[3] );  mat.setCol ( 2, vertices[2] - vertices[3] );

    if ( fabs ( mat.det() ) < tolerance ) {
      cerr << "Matrix not invertible for " << endl
           << "CFE type " << static_cast<int> ( element.pureCFEType() ) << endl
           << "position in image coordinates: " << position[0] << " " << position[1] << " " << position[2] << endl
           << "element: " << x << " " << y << " " << z << endl
           << "tetrahedra vertices: " << vertices[0] << endl << vertices[1] << endl << vertices[2] << endl << vertices[3] << endl
           << "matrix: " << endl << mat << endl;
    }

    matinv = mat.inverse();

    shifted_position = position - vertices[3];

    bary3 = matinv * shifted_position;

    // if point lies inside tetrahedron
    if ( ( bary3[0] >= -tolerance  && bary3[1] >= -tolerance && bary3[2] >= -tolerance && ( bary3[0] + bary3[1] + bary3[2] <= 1 + tolerance ) ) ) {
      for ( short i = 0; i < 3; ++i )
        baryCO[i] = bary3[i];
      baryCO[3] = 1.0 - bary3.sum();
      foundInside = true;
      tetra = tet;
      return;
    }
  }

  foundInside = false;
  tetra = CFETopoTetra();
}


/** Based on CFEGrid, evaluate given data at given position (in world
 *  coordinates) The method first finds out on which (inner)
 *  tetrahedron to interpolate, using which barycentric coordinates,
 *  then does the interpolation.
 *  \author Schwen
 */
template< typename GridType >
typename GridType::RealType interpolateDataAtPositionWC ( const GridType &grid,
                                                          const qc::ScalarArray< typename GridType::RealType, qc::QC_3D> &data,
                                                          const aol::Vec3< typename GridType::RealType > &positionWC,
                                                          const typename GridType::RealType tolerance = 1.0e-14 ) {
  typedef typename GridType::RealType RealType;

  { // first check whether we wish to evaluate at a grid node and just return the value there -- in this case, we need not search for tetrahedron.
    aol::Vec3<RealType> position ( positionWC ); // global = image coordinates
    position *= static_cast<RealType> ( 1.0 / grid.H() );

    const short x = static_cast<short> ( floor ( position[0] ) ), y = static_cast<short> ( floor ( position[1] ) ), z = static_cast<short> ( floor ( position[2] ) );

    if ( aol::Sqr ( position[0] - x ) + aol::Sqr ( position[1] - y ) + aol::Sqr ( position[2] - z ) < tolerance ) {
      return ( data.get ( x, y, z ) );
    }
  }

  tpcfe::CFEElement<RealType> el;
  tpcfe::CFETopoTetra tetra;
  aol::Vec<4, RealType> baryCO;
  bool foundInside = false;

  determineElementAndTetraForWC<GridType> ( grid, positionWC, el, tetra, baryCO, foundInside, tolerance );

  if ( foundInside ) {
    return ( interpolateDataAtPositionBary ( grid, data, el, tetra, baryCO ) );
  }

  // else:
  cerr << endl << "tpCFEGrid::interpolateDataAtPositionWC: no suitable tetrahedron found, returning NaN for " << positionWC << endl;
  return ( aol::NumberTrait<RealType>::NaN );

}


template< typename GridType >
aol::Vec3<typename GridType::RealType> interpolateDataAtPositionWC ( const GridType &grid,
                                                                     const qc::MultiArray< typename GridType::RealType, qc::QC_3D > &data,
                                                                     const aol::Vec3< typename GridType::RealType > &positionWC,
                                                                     const typename GridType::RealType tolerance = 1.0e-14 ) {
  typedef typename GridType::RealType RealType;

  { // first check whether we wish to evaluate at a grid node and just return the value there -- in this case, we need not search for tetrahedron.
    aol::Vec3<RealType> position ( positionWC ); // global = image coordinates
    position *= static_cast<RealType> ( 1.0 / grid.H() );

    const short x = static_cast<short> ( floor ( position[0] ) ), y = static_cast<short> ( floor ( position[1] ) ), z = static_cast<short> ( floor ( position[2] ) );

    if ( aol::Sqr ( position[0] - x ) + aol::Sqr ( position[1] - y ) + aol::Sqr ( position[2] - z ) < tolerance ) {
      return ( data.get ( qc::CoordType ( x, y, z ) ) );
    }
  }

  tpcfe::CFEElement<RealType> el;
  tpcfe::CFETopoTetra tetra;
  aol::Vec<4, RealType> baryCO;
  bool foundInside = false;

  tpcfe::determineElementAndTetraForWC<GridType> ( grid, positionWC, el, tetra, baryCO, foundInside, tolerance );

  if ( foundInside ) {
    return ( interpolateDataAtPositionBary ( grid, data, el, tetra, baryCO ) );
  }

  // else:
  cerr << positionWC << endl;
  throw aol::Exception ( "tpcfe::interpolateDataAtPositionWC: no suitable tetrahedron found", __FILE__, __LINE__ );
}


/** CFE interpolation of discrete data at barycentric coordinates
 *  inside a given tetrahedron in an element.
 *  \author Schwen
 */
template< typename GridType >
typename GridType::RealType interpolateDataAtPositionBary ( const GridType &grid,
                                                            const qc::ScalarArray< typename GridType::RealType, qc::QC_3D> &data,
                                                            const tpcfe::CFEElement<typename GridType::RealType> &el,
                                                            const tpcfe::CFETopoTetra &tetra,
                                                            const aol::Vec<4, typename GridType::RealType > &baryCO ) {
  typedef typename GridType::RealType RealType;

  aol::Vec<4, RealType> valueAtVertex;

  for ( short i = 0; i < 4; ++i ) { // loop over all vertices
    const int gIdx0 = el.globalIndex ( tetra ( i, 0 ) );
    const int gIdx1 = ( ( tetra ( i, 1 ) == NODE_NOT_VIRTUAL ) ?  -1  :  el.globalIndex ( tetra ( i, 1 ) ) );

    if ( gIdx1 == -1 ) { // non-interpolated node
      valueAtVertex[i] = data.get ( gIdx0 );
    } else {           // interpolated node
      const typename GridType::VNType &vn = grid.getVirtualNodeRef ( gIdx0, gIdx1 );
      for ( unsigned char j = 0; j < vn._numOfConstraints; ++j ) {
        valueAtVertex[i] += vn.weight ( j ) * data.get ( vn.constraint ( j ) );
      }
    }
  }

  // linear interpolation based on the barycentric coordinates computed above.
  return ( baryCO * valueAtVertex );
}


template< typename GridType >
aol::Vec3<typename GridType::RealType> interpolateDataAtPositionBary ( const GridType &grid,
                                                                       const qc::MultiArray< typename GridType::RealType, qc::QC_3D > &data,
                                                                       const tpcfe::CFEElement<typename GridType::RealType> &el,
                                                                       const tpcfe::CFETopoTetra &tetra,
                                                                       const aol::Vec< 4, typename GridType::RealType > &baryCO ) {
  typedef typename GridType::RealType RealType;

  aol::Vec3<RealType> valueAtVertex[4];

  for ( short i = 0; i < 4; ++i ) { // loop over all vertices
    const int gIdx0 = el.globalIndex ( tetra ( i, 0 ) );
    const int gIdx1 = ( ( tetra ( i, 1 ) == tpcfe::NODE_NOT_VIRTUAL ) ?  -1  :  el.globalIndex ( tetra ( i, 1 ) ) );

    if ( gIdx1 == -1 ) { // non-interpolated node
      data.getTo ( gIdx0, valueAtVertex[i] );
    } else {           // interpolated node
      const typename GridType::VNType &vn = grid.getVirtualNodeRef ( gIdx0, gIdx1 );
      for ( unsigned char j = 0; j < vn._numOfConstraints; ++j ) {
        aol::Vec3<RealType> dwarf;
        vn.weight ( j ).mult ( data.get ( vn.constraint ( j ) ), dwarf );
        valueAtVertex[i] += dwarf ;
      }
    }
  }

  // linear interpolation based on the barycentric coordinates computed above.
  aol::Vec3<RealType> ret;

  for ( short i = 0; i < 4; ++i )
    for ( short j = 0; j < 3; ++j )
      ret[j] += baryCO[i] * valueAtVertex[i][j];

  return ( ret );
}

//! set coefficients in an AArray depending on sign of the level set function
template< typename RealType, typename NodalCoeffType >
void setCoeffForLevelset ( qc::AArray< NodalCoeffType, qc::QC_3D > &coeff, const qc::ScalarArray< RealType, qc::QC_3D > &levelset, const NodalCoeffType coeffMinus, const NodalCoeffType coeffPlus ) {
  for ( int i = 0; i < levelset.size(); ++i ) {
    if ( levelset[i] < 0 )
      coeff[i] = coeffMinus;
    else if ( levelset[i] > 0 )
      coeff[i] = coeffPlus;
    else {
      cerr << i << " " << levelset[i] << endl;
      throw aol::Exception ( "tpcfe::setCoeffForLevelset: illegal levelset value (0 or nan?)", __FILE__, __LINE__ );
    }
  }
}


/** Class for interface triangles
 *  \author Schwen
 */
template <typename RealType>
class CFEInterfaceTriang {
protected:
  typedef CFEVirtualNode<RealType, CFE_NONE, RealType>    VNType;
  typedef typename VNType::IndexType   IndexType;
  CFEElement<RealType>        _element;   //!< element in which the interface triang lies
  CFETopoTetra                _tetra;     //!< tetrahedron which the interface triang is part of (should be inner tetra)
  aol::Vec3< unsigned char >  _tNodes;    //!< node indices of the tetra that form the triang (in 0..3)

public:
  //! Standard constructor that constructs a dummy
  CFEInterfaceTriang() : _element (), _tetra(), _tNodes() {
    // do nothing else
  }

public:
  //! Copy constructor
  CFEInterfaceTriang ( const CFEInterfaceTriang<RealType> &other ) : _element ( other._element ), _tetra ( other._tetra ), _tNodes ( other._tNodes ) {
    // do nothing else
  }

public:
  //! This constructor should be used
  CFEInterfaceTriang ( const CFEElement<RealType> &el, const CFETopoTetra &tetra, const aol::Vec3< unsigned char > &tns ) : _element ( el ), _tetra ( tetra ), _tNodes ( tns ) {
    // do nothing else
  }

public:
  //! Assignment operator
  CFEInterfaceTriang<RealType>& operator= ( const CFEInterfaceTriang<RealType> &other ) {
    set ( other._element, other._tetra, other._tNodes );
  }

  void set ( const CFEElement<RealType> &el, const CFETopoTetra &tetra, const aol::Vec3< unsigned char > &tns ) {
    _element = el;
    _tetra   = tetra;
    _tNodes  = tns;
  }

  // access functions:
  const qc::Element& getElement() const {
    return _element;
  }

  //! lIdx0 and lIdx1 become the local indices (wrt element) of the two nodes forming the edge
  //! containing vertex i (being a virtual node)
  void getLocalEdgeIndicesWrtElement ( const unsigned char i, unsigned char &lIdx0, unsigned char &lIdx1 ) const {
    lIdx0 = _tetra ( _tNodes[i], 0 );
    lIdx1 = _tetra ( _tNodes[i], 1 );
  }

  //! lIdx0 and lIdx1 become the local indices (wrt regular tetrahedron) of the two nodes forming the edge
  //! containing vertex i (being a virtual node)
  void getLocalEdgeIndicesWrtTetra ( const unsigned char i, unsigned char &lIdx0, unsigned char &lIdx1 ) const {
    const aol::Vec2<unsigned char> triangVertexLocIndWrtElem ( _tetra ( _tNodes[i], 0 ), _tetra ( _tNodes[i], 1 ) );
    const unsigned char regTetraIndWrtElem = _tetra.getParent();

    lIdx0 = lIdx1 = 5;
    for ( unsigned char regTetVI = 0; regTetVI < 4; ++regTetVI ) {
      if ( tpcfe::CFETopoLookup::_stdTetraVertex[ regTetraIndWrtElem ][ regTetVI ] == triangVertexLocIndWrtElem[0] ) {
        lIdx0 = regTetVI;
      }
      if ( tpcfe::CFETopoLookup::_stdTetraVertex[ regTetraIndWrtElem ][ regTetVI ] == triangVertexLocIndWrtElem[1] ) {
        lIdx1 = regTetVI;
      }
    }
    if ( ( lIdx0 == 5 ) || ( lIdx1 == 5 ) || ( lIdx0 == lIdx1 ) ) {
      throw aol::Exception ( "no appropriate index found", __FILE__, __LINE__ );
    }
  }

  //! Return local coordinates of vertices (not scaled by h) in the element
  //! note that this is inefficient
  std::vector< aol::Vec3< RealType > > getLocalCoordinates ( ) const {
    std::vector< aol::Vec3< RealType > > ret ( 3 );
    for ( int i = 0; i < 3; ++i )
      _tetra.computeLocalCoordinate ( ret[i], _element, _tNodes[i] );
    return ( ret );
  }

  //! Return global coordinates of vertices (not scaled by h)
  //! note that this is inefficient
  std::vector< aol::Vec3< RealType > > getGlobalCoordinates ( ) const {
    std::vector< aol::Vec3< RealType > > ret ( 3 );
    for ( int i = 0; i < 3; ++i )
      _tetra.computeGlobalCoordinate ( ret[i], _element, _tNodes[i] );
    return ( ret );
  }

  //! Return global virtual node indices of vertices
  aol::Vec3< IndexType > getVNIndices ( ) const {
    aol::Vec3< IndexType > ret;
    for ( int i = 0; i < 3; ++i ) {
      const int lIdx0 = _tetra ( _tNodes[i], 0 ), lIdx1 = _tetra ( _tNodes[i], 1 );
      const int gIdx0 = _element.globalIndex ( lIdx0 ), gIdx1 = _element.globalIndex ( lIdx1 );

#ifdef DEBUG
      if ( lIdx0 == NODE_NOT_VIRTUAL )
        throw aol::Exception ( "CFEInterfaceTriang: Something is wrong, first node index must not indicate regular node ...", __FILE__, __LINE__ );
#endif

      ret[i] = VNType::mapIndex ( gIdx0, gIdx1 );
    }
    return ( ret );
  }

  //! Return world coordinates of vertices (scaled by h)
  std::vector< aol::Vec3< RealType > > getWorldCoordinates ( const RealType h ) const {
    std::vector< aol::Vec3< RealType > > ret ( 3 );
    for ( int i = 0; i < 3; ++i ) {
      _tetra.computeGlobalCoordinate ( ret[i], _element, _tNodes[i] );
      ret[i] *= h;
    }
    return ( ret );
  }

  //! Compute Area of this triangle.
  RealType getArea ( const RealType h ) const {
    std::vector< aol::Vec3< RealType > > localCoords = getLocalCoordinates();
    aol::Vec3<RealType> temp ( localCoords[2] - localCoords[0] );
    return ( ( h * h / 2 ) * ( temp.crossProduct ( localCoords[2] - localCoords[1] ) ).norm() );
  }

  void dump ( ostream &os = cerr ) const {
    std::vector< aol::Vec3<RealType> > globalCoord = this->getGlobalCoordinates();
    os << globalCoord[0] << " " << globalCoord[1] << " " << globalCoord[2] << endl;
  }

  // Additional functions:

  // RealType getArea(){
  // }
  // getNormal?
  // global indices etc?
  // global / world coordinates?

};


// \todo maybe implement similar virtual tetrahedra iterator to avoid nested element/tetra loop?

template <typename GridType>
class CFEInterfaceTriangleIterator {
protected:
  typedef typename GridType::RealType RealType;

  const GridType &_grid;
  const qc::GridSize<qc::QC_3D> _gridSize;
  qc::RectangularGrid<qc::QC_3D>::FullElementIterator _elit;
  CFETopoTetraIterator _tit;
  CFEInterfaceTriang<RealType> _triang;
  CFEElement< typename GridType::RealType> _currentElement;

public:
  explicit CFEInterfaceTriangleIterator ( const GridType &grid ) : _grid ( grid ), _gridSize ( grid ), _elit ( _grid ), _triang(), _currentElement () {
    restart();
  }

public:
  void restart ( ) {
    _elit.restart ( _grid );
    _currentElement.set ( *_elit, _gridSize, _grid.getElType ( *_elit ) );
    while ( _elit.notAtEnd() && ! ( _currentElement.cfeType().representsInterfaced() ) ) {
      ++_elit;
      if ( _elit.notAtEnd() ) {
        _currentElement.set ( *_elit, _gridSize, _grid.getElType ( *_elit ) );
      } else { // else _currentElement is unusable
        _currentElement.set ( *_elit, _gridSize, CFEType ( ) );
      }
    }

    if ( !_elit.atEnd() ) {
      _grid.getCutRelations ( _currentElement );
      _tit.restart ( _currentElement.cfeType(), -1 );
      next_interface_tet();
      setTriang();
    } else {
      setDummyTriang(); // if no interfaced element exists, loop will terminate right away - still we need to set something (do we?)
    }

  }

  bool atEnd ( ) const {
    return ( ! ( notAtEnd() ) );
  }

  bool notAtEnd ( ) const {
    return ( _elit.notAtEnd() );
  }

  // prefix increment operator
  CFEInterfaceTriang<RealType>& operator++ () {
    increment_tit();
    next_interface_tet();

    setTriang();

    return ( _triang );
  }

  const CFEInterfaceTriang< RealType >& getTriangRef () const {
    return _triang;
  }

  const CFETopoTetra& getTetraRef ( ) const {
    return *_tit;
  }

  const CFEElement<RealType>& getElementRef ( ) const {
    return ( _currentElement );
  }

protected:
  void setTriang ( ) {
    if ( _tit.notAtEnd() ) { // which should only happen at the very end.
      if ( ! ( ( *_tit ).isVirtualNode ( 0 ) ) ) {
        _triang.set ( _currentElement, *_tit, aol::Vec3< unsigned char > ( 1, 2, 3 ) );
      } else if ( ! ( ( *_tit ).isVirtualNode ( 1 ) ) ) {
        _triang.set ( _currentElement, *_tit, aol::Vec3< unsigned char > ( 0, 2, 3 ) );
      } else if ( ! ( ( *_tit ).isVirtualNode ( 2 ) ) ) {
        _triang.set ( _currentElement, *_tit, aol::Vec3< unsigned char > ( 0, 1, 3 ) );
      } else if ( ! ( ( *_tit ).isVirtualNode ( 3 ) ) ) {
        _triang.set ( _currentElement, *_tit, aol::Vec3< unsigned char > ( 0, 1, 2 ) );
#ifdef DEBUG
      } else {
        throw aol::Exception ( "CFEInterfaceTriangleIterator: Something is wrong ...", __FILE__, __LINE__ );
#endif
      }
    } else {
      setDummyTriang();
    }
  }

  void setDummyTriang ( ) {
    // at the end
    CFEElement<RealType> dummyel;
    CFETopoTetra dummytet;
    aol::Vec3<unsigned char> dummyvec;
    _triang.set ( dummyel, dummytet, dummyvec ); // not to be used as we're at the end ...
  }

  //! go forward to next (or stay at this) interface tetrahedron, characterized by having exactly three virtual nodes
  void next_interface_tet ( ) {
    while ( !_elit.atEnd() && ( *_tit ).getVirtualNodeNum() != 3 ) {
      increment_tit();
    }
  }

  void increment_tit ( ) {
    ++_tit;
    if ( _tit.atEnd() ) {
      increment_elit();
    }
  }

  void increment_elit ( ) {
    ++_elit;
    while ( _elit.notAtEnd() && ! ( _grid.getElType ( *_elit ).representsInterfaced() ) ) {
      ++_elit;
    }

    if ( !_elit.atEnd() ) {
      _currentElement.set ( *_elit, _gridSize, _grid.getElType ( *_elit ) );
      _grid.getCutRelations ( _currentElement );
      _tit.restart ( _currentElement.cfeType(), -1 );
    }

  }

private:
  CFEInterfaceTriangleIterator ( );
  CFEInterfaceTriangleIterator ( const CFEInterfaceTriangleIterator<GridType> & );
  CFEInterfaceTriangleIterator<GridType>& operator= ( const CFEInterfaceTriangleIterator<GridType> & );
};



template < class GridType, class Float_Type = float >
class CFEInterfaceTriangulationGenerator {
public:
  typedef typename GridType::RealType                          RealType;
  typedef Float_Type                                           FloatType;

protected:
  typedef typename GridType::VNType                            VNodeType;
  typedef typename VNodeType::IndexType                        VNodeIndType;
  typedef aol::Vec3< VNodeIndType >                            TriangType;

  typedef std::map< VNodeIndType, int >                        VNodeIndexMapperType;
  typedef std::map< TriangType, int >                          TriangIndexMapperType;

public:
  typedef std::vector< aol::Vec3< FloatType > >                VertexVectorType;
  typedef std::vector< aol::Vec3< int > >                      TriangVectorType;
  typedef std::vector< FloatType >                             VertexDataVectorType;
  typedef std::vector< FloatType >                             TriangDataVectorType;

protected:
  const GridType &_grid;

  VertexVectorType       _interfaceVertexVector,     _boundaryVertexVector;
  TriangVectorType       _interfaceTriangVector,     _boundaryTriangVector;
  VertexDataVectorType   _interfaceVertexDataVector, _boundaryVertexDataVector;
  TriangDataVectorType   _interfaceTriangDataVector, _boundaryTriangDataVector;

  VNodeIndexMapperType   _interfaceVertexIndexmapper, _boundaryVertexIndexmapper;
  TriangIndexMapperType  _interfaceTriangIndexmapper, _boundaryTriangIndexmapper;
  bool
  _haveVertexData,
  _haveTriangData,
  _interfaceTriangulationComputed;

  int
  _numInterfaceVertices, _numBoundaryVertices, _numInterfaceTriangs, _numBoundaryTriangs;

public:
  explicit CFEInterfaceTriangulationGenerator ( const GridType &grid ) : _grid ( grid ), _haveVertexData ( false ), _haveTriangData ( false ), _interfaceTriangulationComputed ( false ),
      _numInterfaceVertices ( 0 ), _numBoundaryVertices ( 0 ), _numInterfaceTriangs ( 0 ), _numBoundaryTriangs ( 0 ) {
  }

public:
  virtual ~CFEInterfaceTriangulationGenerator() { }

public:
  void determineInterfaceTriangulation ( ) ;

  void determineInterfaceAndBoundaryTriangulation ( ) {
    determineInterfaceAndBoundaryTriangulation ( aol::Vec3<int> ( 0, 0, 0 ), _grid.getSize() - aol::Vec3<int> ( 1, 1, 1 ) );
  }

  void determineInterfaceAndBoundaryTriangulation ( const aol::Vec3<int> &clip_lower, const aol::Vec3<int> &clip_upper ) {
    determineInterfaceTriangulation();
    determineClippedBoundaryTriangulation ( clip_lower, clip_upper );
  }

  void determineBoundaryTriangulation ( const signed char which ) {
    determineClippedBoundaryTriangulation ( aol::Vec3<int> ( 0, 0, 0 ), _grid.getSize() - aol::Vec3<int> ( 1, 1, 1 ), 1.0e-9, which );
  }

  void determineClippedBoundaryTriangulation ( const aol::Vec3<int> &clip_lower, const aol::Vec3<int> &clip_upper, const RealType geomTolerance  = 1.0e-9, const signed char which = -1 );

  void deformVertices ( const qc::MultiArray<RealType, 3> &deformations, const FloatType scaling ) {
    deformVertices ( _interfaceVertexIndexmapper, _interfaceVertexVector, deformations, scaling );
    deformVertices ( _boundaryVertexIndexmapper,  _boundaryVertexVector,  deformations, scaling );
  }

  void deformVertices ( const qc::MultiArray<RealType, 3> &deformations, const aol::Vec3<FloatType> scaling ) {
    deformVertices ( _interfaceVertexIndexmapper, _interfaceVertexVector, deformations, scaling );
    deformVertices ( _boundaryVertexIndexmapper,  _boundaryVertexVector,  deformations, scaling );
  }

  void saveToPLYFile ( const char* filename ) const;

  void writeToTriangMesh ( aol::TriangMesh< FloatType > &TriangMesh ) const;

  void getVertexVector ( VertexVectorType &vertexVector );
  void getTriangVector ( TriangVectorType &triangVector );
  void getVertexDataVector ( VertexDataVectorType &vertexDataVector );
  void getTriangDataVector ( TriangDataVectorType &triangDataVector );


  //! returns reference to internal data -- use with caution
  const VertexVectorType& getInterfaceVertexVectorRef ( ) const {
    return ( _interfaceVertexVector );
  }
  //! returns reference to internal data -- use with caution
  const VertexVectorType& getBoundaryVertexVectorRef ( ) const {
    return ( _boundaryVertexVector );
  }

  //! returns reference to internal data -- use with caution
  const TriangVectorType& getInterfaceTriangVectorRef ( ) const {
    return ( _interfaceTriangVector );
  }
  //! returns reference to internal data -- use with caution
  const TriangVectorType& getBoundaryTriangVectorRef ( ) const {
    return ( _boundaryTriangVector );
  }

  //! returns reference to internal data -- use with caution
  const VertexDataVectorType& getInterfaceVertexDataVectorRef ( ) const {
    return ( _interfaceVertexDataVector );
  }
  //! returns reference to internal data -- use with caution
  const VertexDataVectorType& getBoundaryVertexDataVectorRef ( ) const {
    return ( _boundaryVertexDataVector );
  }


  //! UNTESTED
  void determineInterfaceTriangulation ( VertexVectorType &interfaceVertices, TriangVectorType &interfaceTriangles ) {
    determineInterfaceTriangulation();
    getVertexVector ( interfaceVertices );
    getTriangVector ( interfaceTriangles );
  }

  //! UNTESTED
  void determineInterfaceAndBoundaryTriangulation ( VertexVectorType &allVertices, TriangVectorType &allTriangles ) {
    determineInterfaceAndBoundaryTriangulation();
    getVertexVector ( allVertices );
    getTriangVector ( allTriangles );
  }

  //! UNTESTED
  void determineClippedBoundaryTriangulation ( VertexVectorType &boundaryVertices, TriangVectorType &boundaryTriangles, const aol::Vec3<int> &clip_lower, const aol::Vec3<int> &clip_upper, const RealType geomTolerance = 1.0e-9 ) {
    determineClippedBoundaryTriangulation ( clip_lower, clip_upper, geomTolerance );
    getVertexVector ( boundaryVertices );
    getTriangVector ( boundaryTriangles );
  }


  void determineSliceTriangulation ( const qc::Comp direction, const short slice, const RealType geomTolerance = 1.0e-9, const signed char sign = 0 );

  int getNumberOfVertices ( ) const {
    return ( static_cast<int> ( _interfaceVertexVector.size() + _boundaryVertexVector.size() ) );
  }

  int getNumberOfTriangles ( ) const {
    return ( static_cast<int> ( _interfaceTriangVector.size() + _boundaryTriangVector.size() ) );
  }

protected:

  virtual void computeTetraData ( const tpcfe::CFEElement<RealType>& /*el*/, const tpcfe::CFETopoTetra& /*tetra*/,
                                  const aol::Vec3<FloatType>& /*v0*/, const aol::Vec3<FloatType>& /*v1*/, const aol::Vec3<FloatType>& /*v2*/, const aol::Vec3<FloatType>& /*v3*/,
                                  const aol::Vec<4, int>& /*m*/, const bool /*onBoundary*/ ) {
    // overload this function on derived classes e g for displaying von Mises stress
  }

  virtual void postprocessInterfaceVertexData ( ) { }  // overload these functions on derived classes where necessary
  virtual void postprocessBoundaryVertexData ( ) { }
  virtual void postprocessInterfaceTriangData ( ) { }
  virtual void postprocessBoundaryTriangData ( ) { }

  virtual void computeVertexData ( const VNodeIndType /*VNodeIndex*/, VertexDataVectorType& /*vertexDataVector*/ ) const { }  // overload this function on derived classes e g for displaying scalar data

  void deformVertices ( VNodeIndexMapperType &IM, VertexVectorType &vertices, const qc::MultiArray<RealType, 3> &deformations, const aol::Vec3<FloatType> scaling );

  void deformVertices ( VNodeIndexMapperType &IM, VertexVectorType &vertices, const qc::MultiArray<RealType, 3> &deformations, const FloatType scaling ) {
    deformVertices ( IM, vertices, deformations, aol::Vec3<FloatType> ( scaling, scaling, scaling ) );
  }

};
// end of class CFEInterfaceTriangulationGenerator



template< class GridType>
class CFEInterfaceTriangulationWithScalarDataGenerator : public CFEInterfaceTriangulationGenerator<GridType> {

  // inherit typedefs

public:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::RealType              RealType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::FloatType             FloatType;

protected:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VNodeType             VNodeType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VNodeIndType          VNodeIndType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::TriangType            TriangType;

  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VNodeIndexMapperType  VNodeIndexMapperType;

public:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VertexVectorType      VertexVectorType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::TriangVectorType      TriangVectorType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VertexDataVectorType  VertexDataVectorType;

protected:
  const qc::ScalarArray<RealType, qc::QC_3D> &_scalarData;

public:
  CFEInterfaceTriangulationWithScalarDataGenerator ( const GridType &grid, const qc::ScalarArray<RealType, qc::QC_3D> &scalarData ) : CFEInterfaceTriangulationGenerator<GridType> ( grid ), _scalarData ( scalarData ) {
    this->_haveVertexData = true;
  }

protected:
  virtual void computeVertexData ( const VNodeIndType VNodeIndex, VertexDataVectorType &vertexDataVector ) const {
    int gIdx0 = -1, gIdx1 = -1;
    VNodeType::splitIndex ( VNodeIndex, gIdx0, gIdx1 );
    if ( gIdx1 == 0 ) {
      vertexDataVector.push_back ( _scalarData[gIdx0] );
    } else {
      const VNodeType &vn = this->_grid.getVirtualNodeRef ( gIdx0, gIdx1 );
      vertexDataVector.push_back ( vn.extrapolate ( _scalarData ) );
    }
  }

};


template< class GridType>
class CFEInterfaceTriangulationWithComponentDataGenerator : public CFEInterfaceTriangulationGenerator<GridType> {

  // inherit typedefs

public:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::RealType              RealType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::FloatType             FloatType;

protected:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VNodeType             VNodeType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VNodeIndType          VNodeIndType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::TriangType            TriangType;

  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VNodeIndexMapperType  VNodeIndexMapperType;

public:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VertexVectorType      VertexVectorType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::TriangVectorType      TriangVectorType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VertexDataVectorType  VertexDataVectorType;

protected:
  const qc::MultiArray<RealType, qc::QC_3D> &_displacement;
  const qc::Comp _component;

public:
  CFEInterfaceTriangulationWithComponentDataGenerator ( const GridType &Grid, const qc::MultiArray<RealType, qc::QC_3D> &Displacement, qc::Comp Component )
      : CFEInterfaceTriangulationGenerator<GridType> ( Grid ),
      _displacement ( Displacement ),
      _component ( Component ) {
    this->_haveVertexData = true;
  }

protected:
  virtual void computeVertexData ( const VNodeIndType VNodeIndex, VertexDataVectorType &vertexDataVector ) const {
    int gIdx0 = -1, gIdx1 = -1;
    VNodeType::splitIndex ( VNodeIndex, gIdx0, gIdx1 );
    if ( gIdx1 == 0 ) {
      vertexDataVector.push_back ( _displacement[ _component ][ gIdx0 ] );
    } else {
      const VNodeType &vn = this->_grid.getVirtualNodeRef ( gIdx0, gIdx1 );
      aol::Vec3<RealType> value;
      vn.extrapolate ( _displacement, value );
      vertexDataVector.push_back ( value[ _component ] );
    }
  }

};


template< class GridType >
class CFEInterfaceTriangulationWithVonMisesStressGeneratorBase : public CFEInterfaceTriangulationGenerator<GridType> {

  // inherit typedefs

public:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::RealType              RealType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::FloatType             FloatType;

protected:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VNodeType             VNodeType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VNodeIndType          VNodeIndType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::TriangType            TriangType;

  typedef typename CFEInterfaceTriangulationGenerator<GridType>::TriangIndexMapperType TriangIndexMapperType;
  typedef std::map< TriangType, float >                                                TriangDataMapType;

public:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::TriangVectorType      TriangVectorType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::TriangDataVectorType  TriangDataVectorType;

protected:
  TriangDataMapType _interfaceTriangDataCollector, _boundaryTriangDataCollector;
  const qc::MultiArray<RealType, 3> &_deformation;
  const RealType _aleph;

public:
  CFEInterfaceTriangulationWithVonMisesStressGeneratorBase ( const GridType &grid, const qc::MultiArray<RealType, 3> &deformation, const RealType aleph ) :
      CFEInterfaceTriangulationGenerator<GridType> ( grid ), _deformation ( deformation ), _aleph ( aleph ) {
    this->_haveTriangData = true;
  }

public:
  using CFEInterfaceTriangulationGenerator<GridType>::deformVertices;
  void deformVertices ( const FloatType scaling ) {
    this->deformVertices ( _deformation, scaling );
  }

  void deformVertices ( const aol::Vec3<FloatType> scaling ) {
    this->deformVertices ( _deformation, scaling );
  }

protected:
  virtual void postprocessInterfaceTriangData ( ) {
    this->_interfaceTriangDataVector.resize ( this->_interfaceTriangVector.size() );
    for ( typename TriangDataMapType::iterator it = _interfaceTriangDataCollector.begin(); it != _interfaceTriangDataCollector.end(); ++it ) {
      this->_interfaceTriangDataVector[ this->_interfaceTriangIndexmapper[ it->first ] ] = it->second;
    }
    _interfaceTriangDataCollector.clear();
  }

  virtual void postprocessBoundaryTriangData ( ) {
    this->_boundaryTriangDataVector.resize ( this->_boundaryTriangVector.size() );
    for ( typename TriangDataMapType::iterator it = _boundaryTriangDataCollector.begin(); it != _boundaryTriangDataCollector.end(); ++it ) {
      this->_boundaryTriangDataVector[ this->_boundaryTriangIndexmapper[ it->first ] ] = it->second;
    }
    _boundaryTriangDataCollector.clear();
  }

  virtual void computeTetraData ( const tpcfe::CFEElement<RealType> &el, const tpcfe::CFETopoTetra &tetra,
                                  const aol::Vec3<FloatType> &v0, const aol::Vec3<FloatType> &v1, const aol::Vec3<FloatType> &v2, const aol::Vec3<FloatType> &v3,
                                  const aol::Vec<4, int> &m, const bool onBoundary ) ;

  virtual void computeStress ( const aol::Matrix33<RealType> &gradient_deformations, const tpcfe::CFEElement<RealType> &el, const tpcfe::CFETopoTetra &tetra, aol::Matrix33<RealType> &sigma ) = 0;
};


//! basis class for template specialization
template< typename GridType >
class CFEInterfaceTriangulationWithVonMisesStressGenerator : public CFEInterfaceTriangulationWithVonMisesStressGeneratorBase < GridType > {
};


//! complicated domain, microscopic isotropy given by Lame-Navier constants
template< class _RealType >
class CFEInterfaceTriangulationWithVonMisesStressGenerator< tpcfe::CFEGrid< _RealType, tpcfe::CFE_CD, _RealType > > : public CFEInterfaceTriangulationWithVonMisesStressGeneratorBase < tpcfe::CFEGrid< _RealType, tpcfe::CFE_CD, _RealType > > {
public:
  typedef tpcfe::CFEGrid< _RealType, tpcfe::CFE_CD, _RealType >             GridType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::RealType   RealType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::FloatType  FloatType;

protected:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VNodeType             VNodeType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::TriangType            TriangType;
  typedef std::map< TriangType, float >                                                TriangDataMapType;

  const RealType _lambda, _mu;

public:
  CFEInterfaceTriangulationWithVonMisesStressGenerator ( const GridType &grid, const qc::MultiArray<RealType, 3> &deformation, const RealType lambda, const RealType mu, const RealType aleph ) :
      CFEInterfaceTriangulationWithVonMisesStressGeneratorBase<GridType> ( grid, deformation, aleph ), _lambda ( lambda ) , _mu ( mu ) {
  }

protected:
  virtual void computeStress ( const aol::Matrix33<RealType> &gradientDeformations, const tpcfe::CFEElement<RealType> &, const tpcfe::CFETopoTetra &, aol::Matrix33<RealType> &sigma ) {
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        sigma.set ( i, j, _mu * ( gradientDeformations[i][j] + gradientDeformations[j][i] ) );
      }
      sigma.add ( i, i, _lambda * ( gradientDeformations[0][0] + gradientDeformations[1][1] + gradientDeformations[2][2] ) );
    }
  }

};


//! jumping coefficient elasticity,
template< class _RealType >
class CFEInterfaceTriangulationWithVonMisesStressGenerator< tpcfe::CFEGrid< _RealType, tpcfe::CFE_TPOSELAST, tpcfe::IsotropicElasticityCoefficient<_RealType> > > : public CFEInterfaceTriangulationWithVonMisesStressGeneratorBase < tpcfe::CFEGrid< _RealType, tpcfe::CFE_TPOSELAST, tpcfe::IsotropicElasticityCoefficient<_RealType> > > {
public:
  typedef tpcfe::CFEGrid< _RealType, tpcfe::CFE_TPOSELAST, tpcfe::IsotropicElasticityCoefficient<_RealType> >  GridType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::RealType                                      RealType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::FloatType                                     FloatType;

protected:
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::VNodeType                                     VNodeType;
  typedef typename CFEInterfaceTriangulationGenerator<GridType>::TriangType                                    TriangType;
  typedef std::map< TriangType, float >                                                                        TriangDataMapType;


  typedef typename GridType::NodalCoeffType                                                                   NodalCoeffType;

  const qc::AArray< NodalCoeffType, qc::QC_3D > &_nodalCoeff;

public:
  CFEInterfaceTriangulationWithVonMisesStressGenerator ( const GridType &grid, const qc::AArray< NodalCoeffType, qc::QC_3D > &NodalCoeff, const qc::MultiArray<RealType, 3> &deformation, const RealType aleph ) :
      CFEInterfaceTriangulationWithVonMisesStressGeneratorBase<GridType> ( grid, deformation, aleph ), _nodalCoeff ( NodalCoeff ) {
  }

protected:
  virtual void computeStress ( const aol::Matrix33<RealType> &gradientDeformations, const tpcfe::CFEElement<RealType> &el, const tpcfe::CFETopoTetra &tetra, aol::Matrix33<RealType> &sigma )  {
    aol::Matrix33<RealType> epsilon;

    for ( short i = 0; i < 3; ++i )
      for ( short j = 0; j < 3; ++j )
        epsilon[i][j] += 0.5 * ( gradientDeformations[i][j] + gradientDeformations[j][i] );

    RealType localTensor[3][3][3][3];

    CFEWeightProvider<RealType, NodalCoeffType> weightProvider ( _nodalCoeff, el );
    weightProvider.meanWeight ( tetra.getSign() ).getAnisotropicTensor ( localTensor );

    sigma.setZero();

    for ( short i = 0; i < 3; ++i )
      for ( short j = 0; j < 3; ++j )
        for ( short k = 0; k < 3; ++k )
          for ( short l = 0; l < 3; ++l )
            sigma[i][j] += localTensor[i][j][k][l] * epsilon [k][l];
  }

};

} // end of namespace tpcfe

#endif

