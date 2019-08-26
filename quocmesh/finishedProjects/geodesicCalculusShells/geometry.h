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

#ifndef __GEOMETRY_H
#define __GEOMETRY_H

#include <aol.h>

using namespace aol;
using namespace qc;

//!=============================================================================================================
//!=============================================================================================================

//! \brief Class which only stores all important topology information of a given mesh.
//! \author Heeren
template< typename MeshType >
class MeshTopologySaver {

  static const int IndexNotSet = -1;  
  const MeshType& _mesh;

  // saves global indices of nodes belonging to each triangle
  RandomAccessContainer< aol::Vec3<int> > _nodesOfTriangles;
  // saves global indices of neighbourin triangles for each triangle (these are -1 at he boundary )
  RandomAccessContainer< aol::Vec3<int> > _neighboursOfTriangles;
  // saves global indices of node in neighboring triangle which is not contained in this
  RandomAccessContainer< aol::Vec3<int> > _oppositeNodesOfNeighboringTriangles;

  // saves global indices of neighbourin triangles for each edge (these are -1 at he boundary )
  RandomAccessContainer< aol::Vec2<int> > _neighboursOfEdges;
  
  // saves global indices of nodes belonging to each edge
  RandomAccessContainer< aol::Vec2<int> > _nodesOfEdges;
  // saves global indices of nodes opposite of each edge (one of them is -1 for boundary edges)
  RandomAccessContainer< aol::Vec2<int> > _oppositeNodesOfEdges;

  const int _numOfElems;
  const int _numOfNodes;
  const int _numOfEdges;

public:
  MeshTopologySaver( const MeshType& mesh )
      : _mesh(mesh),
        _nodesOfTriangles( mesh.getNumFaces(), aol::Vec3<int>(IndexNotSet) ),
        _neighboursOfTriangles( mesh.getNumFaces(), aol::Vec3<int>(IndexNotSet) ),
        _oppositeNodesOfNeighboringTriangles( mesh.getNumFaces(), aol::Vec3<int>(IndexNotSet) ),
        _neighboursOfEdges( mesh.getNumEdges(), aol::Vec2<int>(IndexNotSet) ),
        _nodesOfEdges( mesh.getNumEdges(), aol::Vec2<int>(IndexNotSet) ),
        _oppositeNodesOfEdges( mesh.getNumEdges(), aol::Vec2<int>(IndexNotSet) ),
        _numOfElems( mesh.getNumFaces() ),
        _numOfNodes( mesh.getNumVertices() ),
        _numOfEdges( mesh.getNumEdges() )
  {
     //
     for ( typename MeshType::ElementIteratorType iter = mesh; iter.notAtEnd(); ++iter ){
       int idx = iter.getIndex();
       _nodesOfTriangles[idx] = mesh.getNodeIndices( idx );
       _neighboursOfTriangles[idx] =  mesh.getTriangleIndicesOfNeighbours( idx );
       _oppositeNodesOfNeighboringTriangles[idx] = mesh.getOppositeNodesOfNeighboringTriangles( idx );
     }

     //
     for ( typename MeshType::EdgeIteratorType iter = mesh; iter.notAtEnd(); ++iter )
     {
       int idx = iter.getIndex();
       iter.getAdjacentNodes( _nodesOfEdges[idx][0], _nodesOfEdges[idx][1] );
       iter.getOppositeNodes( _oppositeNodesOfEdges[idx][0], _oppositeNodesOfEdges[idx][1] );
       iter.getAdjacentFaces( _neighboursOfEdges[idx][0], _neighboursOfEdges[idx][1] );       
     }

  }
  
  const MeshType& getGrid() const {
    return _mesh;
  }

  int getNodeOfTriangle( int globalFaceIndex, int localNodeIndex ) const {
    return _nodesOfTriangles[globalFaceIndex][localNodeIndex];
  }

  int getNeighbourOfTriangle( int globalFaceIndex, int localNodeIndex ) const {
    return _neighboursOfTriangles[globalFaceIndex][localNodeIndex];
  }

  int getOppositeNodeOfNeighbouringTriangle( int globalFaceIndex, int localNodeIndex ) const {
    return _oppositeNodesOfNeighboringTriangles[globalFaceIndex][localNodeIndex];
  }

  int getAdjacentNodeOfEdge( int globalEdgeIndex, int localIndex ) const{
    return _nodesOfEdges[globalEdgeIndex][localIndex];
  }
  
  int getAdjacentTriangleOfEdge( int globalEdgeIndex, int localIndex ) const{
    return _neighboursOfEdges[globalEdgeIndex][localIndex];
  }

  int getOppositeNodeOfEdge( int globalEdgeIndex, int localIndex ) const{
    return _oppositeNodesOfEdges[globalEdgeIndex][localIndex];
  }

  int getNumEdges() const {
    return _numOfEdges;
  }

  int getNumFaces() const {
    return _numOfElems;
  }

  int getNumVertices() const {
    return _numOfNodes;
  }
};


//! Reconstructs a mesh from its geometry (i.e. its nodal positions given as a MultiVector)
//! and its topology ( i.e. its connectivity given as a MeshTopologySaver)
//! \author Heeren
//! The mesh has to be empty when calling this function!
template< typename MeshType >
void createMeshFromGeometryAndTopology( const MeshTopologySaver<MeshType>& topology,
                                        const aol::MultiVector<typename MeshType::RealType>& geometry,
                                        MeshType& mesh ){
  if( (mesh.getNumVertices() > 0) || (mesh.getNumFaces() > 0) )
    mesh.clearMesh();  
  if( geometry.numComponents() != 3 )
    throw aol::Exception ( "createMeshFromGeometryAndTopology(): geometry has wrong size!", __FILE__, __LINE__ );  
  if( geometry[0].size() != topology.getNumVertices() )
    throw aol::Exception ( "createMeshFromGeometryAndTopology(): wrong sizes!", __FILE__, __LINE__ );

  for ( int i = 0; i < topology.getNumVertices(); i++ ){
    aol::Vec3<typename MeshType::RealType> coords;
    geometry.getTo(i, coords);
    mesh.addVertex( coords );
  }

  for ( int i = 0; i < topology.getNumFaces(); i++ ){
    aol::Vec3<int> idx;
    for( int j = 0; j < 3; j++ )
      idx[j] = topology.getNodeOfTriangle( i, j );
    mesh.addFace( idx );
  }
}

//!=============================================================================================================
//!==============================================================================================================
// Average = ( A + B ) / 2
template <typename VectorType >
void averageVector( const VectorType& A, const VectorType& B, VectorType& Average ){
  Average.reallocate( A );
  Average = A;
  Average += B;
  Average /= 2.;
}

// B = ( A + Res ) / 2 => Res = 2B - A
template <typename VectorType >
void extrapolateVector( const VectorType& A, const VectorType& B, VectorType& Result ){
  Result.reallocate( A );
  Result = B;
  Result *= 2.;
  Result -= A;  
}
  
// dot product
template<typename RealType>
RealType dotProduct( const aol::Vec3<RealType>& a, const aol::Vec3<RealType>& b ) {
     return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// squared edge length of edge connecting the two given node positions
template<typename RealType>
RealType lengthSqr( const aol::Vec3<RealType>& a, const aol::Vec3<RealType>& b ) {
    return dotProduct( a-b, a-b );
}

//
template<typename RealType, int Dim>
void getWeightedMatrixSum( RealType a, const aol::Mat<Dim,Dim,RealType>& mat1, RealType b, const aol::Mat<Dim,Dim,RealType>& mat2, aol::Mat<Dim,Dim,RealType>& res ) {
    for( int i=0; i<Dim; i++ )
      for( int j=0; j<Dim; j++ )
        res.set( i, j, a*mat1.get(i,j) + b*mat2.get(i,j) );
}

//
template<typename RealType>
void getWeightedVectorSum( RealType a, const aol::Vec3<RealType>& vec1, RealType b, const aol::Vec3<RealType>& vec2, aol::Vec3<RealType>& res ) {
    for( int i=0; i<3; i++ )
      res[i] = a*vec1[i]+b*vec2[i];
}

template<typename RealType>
void getDx( const aol::RandomAccessContainer< aol::Vec3<RealType> >& P, aol::Mat<3,2,RealType>& Dx ) {
  Dx.setCol( 0, P[0]-P[1] );
  Dx.setCol( 1, P[0]-P[2] );
}

template<typename RealType>
void getDx( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk,  aol::Mat<3,2,RealType>& Dx ) {
  Dx.setCol( 0, Pi-Pj );
  Dx.setCol( 1, Pi-Pk );
}

template<typename RealType>
void getMetric( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Matrix22<RealType>& g ) {
  aol::Mat<3,2,RealType> Dx;
  getDx( Pi, Pj, Pk, Dx);
  g.makeProductAtransposedB( Dx, Dx );
}

template<typename RealType>
void getMetric( const aol::RandomAccessContainer< aol::Vec3<RealType> >& P, aol::Matrix22<RealType>& g ) {
  getMetric( P[0], P[1], P[2], g );
}


// squared area of triangle (Pi, Pj, Pk)
template<typename RealType>
RealType getAreaSqr( const aol::Vec3<RealType>& Pi,
                     const aol::Vec3<RealType>& Pj,
                     const aol::Vec3<RealType>& Pk ) {
  aol::Vec3<RealType> normal;
  normal.makeCrossProduct( Pk-Pj, Pi-Pk );
  return normal.normSqr() / 4. ;
}

// squared area of triangle (Pi, Pj, Pk) (here stored in a container)
template<typename RealType>
RealType getAreaSqr( const aol::RandomAccessContainer< aol::Vec3<RealType> >& P ) {
     assert( P.size() == 3 );
    return getAreaSqr( P[0], P[1], P[2] );
}

// returns dihedral angle between n1 and n2 at an edge e
template<typename RealType>
RealType getDihedralAngle( const aol::Vec3<RealType>& nk,
                           const aol::Vec3<RealType>& nl,
                           const aol::Vec3<RealType>& e ) {
    aol::Vec3<RealType> crossprod;
    crossprod.makeCrossProduct( nk, nl);
    //return std::asin( temp*e / e.norm() );
    return crossprod*e < 0. ? -std::acos( min( nk*nl, 1. ) ) : std::acos( min( nk*nl, 1. ) );
}

//
template<typename RealType>
RealType getDihedralAngle( const aol::Vec3<RealType>& Pi,
                           const aol::Vec3<RealType>& Pj,
                           const aol::Vec3<RealType>& Pk,
                           const aol::Vec3<RealType>& Pl ) {
    aol::Vec3<RealType> nk, nl;
    nk.makeCrossProduct( Pk-Pj, Pi-Pk );
    nk.normalize();
    nl.makeCrossProduct( Pl-Pi, Pj-Pl );
    nl.normalize();
    return getDihedralAngle( nk, nl, Pj-Pi );
}
 
// cotang at vertex i in a triangle (Pi, Pj, Pk)
template<typename RealType>
RealType getCotang ( const aol::Vec3<RealType>& Pi,
                  const aol::Vec3<RealType>& Pj,
                  const aol::Vec3<RealType>& Pk ) {
    return dotProduct( Pj-Pi, Pk-Pi ) / ( 2.*sqrt( getAreaSqr(Pi, Pj, Pk) ) );
}

// derivative of the area of a triangle (Pi, Pj, Pk) w.r.t. node Pk
template<typename RealType>
void nablaAreak( const aol::Vec3<RealType>& Pi,
                 const aol::Vec3<RealType>& Pj,
                 const aol::Vec3<RealType>& Pk,
                 aol::Vec3<RealType>& areak ) {
    aol::Vec3<RealType> temp;
    temp.makeCrossProduct( Pj-Pi, Pk-Pj );
    temp /= 2.*temp.norm();
    areak.makeCrossProduct( temp, Pj-Pi);
}

// returns negative(!!!) gradient of dihedral angle
template<typename RealType>
void nablaThetak( const aol::Vec3<RealType>& Pi,
                  const aol::Vec3<RealType>& Pj,
                  const aol::Vec3<RealType>& Pk,
                  aol::Vec3<RealType>& thetak ) {
    thetak.makeCrossProduct( Pk-Pj, Pi-Pk );
    thetak *= aol::Vec3<RealType>(Pj-Pi).norm() / thetak.normSqr();
}

//!====================================================================================================
//! TAKEN FROM PETERS MATHEMATICA CODE

// forwar declaration
template<typename RealType>
void getHijDiffQuotient( const aol::Vec3<RealType>&, const aol::Vec3<RealType>&, const aol::Vec3<RealType>&, const aol::Vec3<RealType>&, aol::Matrix33<RealType>&);
 

// 
template<typename RealType>
void getCrossOp( const aol::Vec3<RealType>& a, aol::Matrix33<RealType>& matrix ) {
  matrix.setZero();
  matrix.set( 0, 1, -a[2]); matrix.set( 0, 2,  a[1]);
  matrix.set( 1, 0,  a[2]); matrix.set( 1, 2, -a[0]);
  matrix.set( 2, 0, -a[1]); matrix.set( 2, 1,  a[0]);
}

template<typename RealType>
void getProjection( const aol::Vec3<RealType>& x, aol::Matrix33<RealType>& m ) {
  m.setIdentity();
  aol::Matrix33<RealType> temp;
  temp.makeTensorProduct( x, x );
  m.addMultiple( temp, -1.0);
}

template<typename RealType>  
void getReflection( const aol::Vec3<RealType>& x, aol::Matrix33<RealType>& m ) {
  m.setIdentity();
  aol::Matrix33<RealType> temp;
  temp.makeTensorProduct( x, x );
  m.addMultiple( temp, -2.0);
}

template<typename RealType>
RealType getArea( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk ) {
  aol::Vec3<RealType> temp;   
  temp.makeCrossProduct( Pk-Pj, Pi-Pk );
  return 0.5 * temp.norm();
}

template<typename RealType>
RealType getWeightedNormal( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Vec3<RealType>& normal  ) {
  normal.makeCrossProduct( Pk-Pj, Pi-Pk);
  return normal.norm();
}

template<typename RealType>
void getNormal( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Vec3<RealType>& normal  ) {
  normal.makeCrossProduct( Pk-Pj, Pi-Pk);
  normal.normalize();
}

template<typename RealType>
void getAreaGradK( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Vec3<RealType>& grad  ) {
  aol::Vec3<RealType> a(Pi-Pk), d(Pk-Pj), e(Pj-Pi);
  RealType area = getArea( Pi, Pj, Pk );
  RealType temp1( -0.25 * dotProduct(e,a) / area ), temp2( 0.25 * dotProduct(e,d) / area );
  getWeightedVectorSum( temp1, d, temp2, a, grad );
}

template<typename RealType>
void getLengthGradk( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, aol::Vec3<RealType>& grad  ) {
  grad = Pi - Pj;
  grad.normalize();
}

template<typename RealType>
void getThetaGradK( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Vec3<RealType>& grad  ) {
  aol::Vec3<RealType> e(Pj-Pi);
  getNormal( Pi, Pj, Pk, grad );
  grad *= -0.5 * e.norm() / getArea( Pi, Pj, Pk );
}

template<typename RealType>
void getThetaGradILeftPart( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Vec3<RealType>& grad  ) {
  aol::Vec3<RealType> e(Pj-Pi), d(Pk-Pj);
  getThetaGradK( Pi, Pj, Pk, grad );
  grad *= dotProduct(d,e) / dotProduct(e,e);
}

template<typename RealType>
void getThetaGradJLeftPart( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Vec3<RealType>& grad  ) {
  aol::Vec3<RealType> e(Pj-Pi), a(Pi-Pk);
  getThetaGradK( Pi, Pj, Pk, grad );
  grad *= dotProduct(a,e) / dotProduct(e,e);
}

template<typename RealType>
void getThetaGradI( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, const aol::Vec3<RealType>& Pl, aol::Vec3<RealType>& grad  ) {
  aol::Vec3<RealType> temp;
  getThetaGradILeftPart( Pi, Pj, Pk, grad );
  getThetaGradILeftPart( Pi, Pj, Pl, temp );
  grad -= temp;
}

template<typename RealType>
void getThetaGradJ( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, const aol::Vec3<RealType>& Pl, aol::Vec3<RealType>& grad  ) {
  aol::Vec3<RealType> temp;
  getThetaGradJLeftPart( Pi, Pj, Pk, grad );
  getThetaGradJLeftPart( Pi, Pj, Pl, temp );
  grad -= temp;
}


template<typename RealType>
void getHessAreaKK( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Matrix33<RealType>& Hess ) {
  aol::Vec3<RealType> eNormalized(Pj-Pi), gradAreaK;
  aol::Matrix33<RealType> proj;
  eNormalized.normalize();
  getAreaGradK( Pi, Pj, Pk, gradAreaK );
  Hess.makeTensorProduct( gradAreaK, gradAreaK );
  getProjection( eNormalized, proj );
  Hess.addMultiple( proj, -0.25 * dotProduct(Pj-Pi,Pj-Pi) );
  Hess *= -1. / getArea( Pi, Pj, Pk );
}

template<typename RealType>
void getHessAreaIK( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Matrix33<RealType>& Hess ) {
  aol::Vec3<RealType> e(Pj-Pi), d(Pk-Pj), temp1, temp2;
  getAreaGradK( Pj, Pk, Pi, temp1 );
  getAreaGradK( Pi, Pj, Pk, temp2 );
  Hess.makeTensorProduct( temp1, temp2 );
  aol::Matrix33<RealType> auxMat;  
  auxMat.makeTensorProduct( e, d );
  Hess.addMultiple( auxMat, 0.25 );  
  Hess.addToDiagonal( -0.25 * dotProduct(d,e) );     
  Hess *= -1. / getArea( Pi, Pj, Pk );  
  getNormal( Pi, Pj, Pk, temp1);
  getCrossOp( temp1, auxMat );
  Hess.addMultiple( auxMat, 0.5 );
}

template<typename RealType>
void getAreaIKDiffQuotient( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Matrix33<RealType>& Hik ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    aol::Vec3<RealType> PiPlusH( Pi ), grad1, grad2;    
    PiPlusH[i] += H;
    getAreaGradK( PiPlusH, Pj, Pk, grad1 );
    getAreaGradK( Pi, Pj, Pk, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hik.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHessThetaKK( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Matrix33<RealType>& Hkk ) {
  
  RealType areaSqr = getArea( Pi, Pj, Pk ) * getArea( Pi, Pj, Pk );
   
  aol::Vec3<RealType> e(Pj-Pi), gradArea, normal;
  getAreaGradK( Pi, Pj, Pk, gradArea );
  getNormal( Pi, Pj, Pk, normal );
        
  aol::Matrix33<RealType> mat1, mat2;    
  getCrossOp( e, mat1 );
  mat2.makeTensorProduct( gradArea, normal );
    
  getWeightedMatrixSum( e.norm() / (4. * areaSqr), mat1,  e.norm() / areaSqr, mat2, Hkk );
}

template<typename RealType>
void getHessThetaIK( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Matrix33<RealType>& Hik ) {
        
  RealType area = getArea( Pi, Pj, Pk );
  RealType areaSqr = area * area;
   
  aol::Vec3<RealType> e(Pj-Pi), d(Pk-Pj), gradArea, normal;    
  getAreaGradK( Pj, Pk, Pi, gradArea );
  getNormal( Pi, Pj, Pk, normal );
    
  aol::Matrix33<RealType> mat1, mat2, mat3;    
  mat1.makeTensorProduct( e, normal );
  getCrossOp( d, mat2 );    
  getWeightedMatrixSum( 1. / (2.*area*e.norm()), mat1,  e.norm() / (4.*areaSqr), mat2, mat3 );
    
  mat1.makeTensorProduct( gradArea, normal );
  getWeightedMatrixSum( 1., mat3, e.norm() / areaSqr, mat1, Hik );
}

template<typename RealType>
void getHessThetaJK( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Matrix33<RealType>& Hjk ) {
    
  RealType area = getArea( Pi, Pj, Pk );
  RealType areaSqr = area * area;
    
  aol::Vec3<RealType> e(Pi-Pj), a(Pi-Pk), gradArea, normal;      
  getAreaGradK( Pk, Pi, Pj, gradArea );
  getNormal( Pi, Pj, Pk, normal );
    
  aol::Matrix33<RealType> mat1, mat2, mat3;    
  mat1.makeTensorProduct( e, normal );
  getCrossOp( a, mat2 );    
  getWeightedMatrixSum( 1. / (2.*area*e.norm()), mat1,  e.norm() / (4.*areaSqr), mat2, mat3 );
    
  mat1.makeTensorProduct( gradArea, normal );
  getWeightedMatrixSum( 1., mat3, e.norm() / areaSqr, mat1, Hjk );
}

template<typename RealType>
void getHessThetaILeftPartI( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Matrix33<RealType>& HiLeft ) {
  aol::Vec3<RealType> e(Pj-Pi), d(Pk-Pj), eNormalized(Pj-Pi), gradThetaK, temp;
  eNormalized.normalize();
  aol::Matrix33<RealType> mat1, mat2, Refl;
  getThetaGradK( Pi, Pj, Pk, gradThetaK );
  getReflection( eNormalized, Refl );
  Refl.mult( d, temp );
  mat1.makeTensorProduct( temp, gradThetaK );
  getHessThetaIK( Pi, Pj, Pk, mat2 );
  getWeightedMatrixSum( -1. / dotProduct(e,e), mat1, dotProduct(d,e) / dotProduct(e,e), mat2, HiLeft );
}

template<typename RealType>
void getHessThetaJLeftPartI( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, aol::Matrix33<RealType>& HjLeft ) {
  aol::Vec3<RealType> e(Pj-Pi), d(Pk-Pj), eNormalized(Pj-Pi), gradThetaK, temp, thetak;
  eNormalized.normalize();
  aol::Matrix33<RealType> mat1, mat2, Refl;
  getThetaGradK( Pi, Pj, Pk, gradThetaK );
  getReflection( eNormalized, Refl );
  Refl.mult( d-e, temp );
  mat1.makeTensorProduct( temp, gradThetaK );
  getHessThetaJK( Pi, Pj, Pk, mat2 );
  getWeightedMatrixSum( 1. / dotProduct(e,e), mat1, dotProduct(d,e) / dotProduct(e,e), mat2, HjLeft );
 }

template<typename RealType>
void getHessThetaII( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, const aol::Vec3<RealType>& Pl, aol::Matrix33<RealType>& Hii ) {    
  aol::Matrix33<RealType> temp;
  getHessThetaILeftPartI(Pi, Pj, Pk, Hii);
  getHessThetaILeftPartI(Pi, Pj, Pl, temp);
  Hii -= temp;
}

template<typename RealType>
void getHessThetaJI( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, const aol::Vec3<RealType>& Pl, aol::Matrix33<RealType>& Hji ) {
  //!TODO replace difference quotient by derivative which is for some reason buggy!
  getHijDiffQuotient( Pi, Pj, Pk, Pl, Hji );
  return;
  
  aol::Matrix33<RealType> temp;
  getHessThetaJLeftPartI(Pi, Pj, Pk, Hji);
  getHessThetaJLeftPartI(Pi ,Pj, Pl, temp);
  Hji -= temp;
}

// DIFFERENCE QUOTIENTS

template<typename RealType>
void getHiiDiffQuotient( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, const aol::Vec3<RealType>& Pl, aol::Matrix33<RealType>& Hii ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    aol::Vec3<RealType> PiPlusH( Pi ), grad1, grad2;    
    PiPlusH[i] += H;
    getThetaGradI( PiPlusH, Pj, Pk, Pl, grad1 );
    getThetaGradI( Pi, Pj, Pk, Pl, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hii.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHjjDiffQuotient( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, const aol::Vec3<RealType>& Pl, aol::Matrix33<RealType>& Hjj ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    aol::Vec3<RealType> PjPlusH( Pj ), grad1, grad2;    
    PjPlusH[i] += H;
    getThetaGradJ( Pi, PjPlusH, Pk, Pl, grad1 );
    getThetaGradJ( Pi, Pj, Pk, Pl, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hjj.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHijDiffQuotient( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, const aol::Vec3<RealType>& Pl, aol::Matrix33<RealType>& Hij ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    aol::Vec3<RealType> PjPlusH( Pj ), grad1, grad2;    
    PjPlusH[i] += H;
    getThetaGradI( Pi, PjPlusH, Pk, Pl, grad1 );
    getThetaGradI( Pi, Pj, Pk, Pl, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hij.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHjiDiffQuotient( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, const aol::Vec3<RealType>& Pl, aol::Matrix33<RealType>& Hji ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    aol::Vec3<RealType> PiPlusH( Pi ), grad1, grad2;    
    PiPlusH[i] += H;
    getThetaGradJ( PiPlusH, Pj, Pk, Pl, grad1 );
    getThetaGradJ( Pi, Pj, Pk, Pl, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hji.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHikDiffQuotient( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, const aol::Vec3<RealType>& Pl, aol::Matrix33<RealType>& Hik ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    aol::Vec3<RealType> PiPlusH( Pi ), grad1, grad2;    
    PiPlusH[i] += H;
    getThetaGradK( PiPlusH, Pj, Pk, grad1 );
    getThetaGradK( Pi, Pj, Pk, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hik.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHkiDiffQuotient( const aol::Vec3<RealType>& Pi, const aol::Vec3<RealType>& Pj, const aol::Vec3<RealType>& Pk, const aol::Vec3<RealType>& Pl, aol::Matrix33<RealType>& Hki ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    aol::Vec3<RealType> PkPlusH( Pk ), grad1, grad2;    
    PkPlusH[i] += H;
    getThetaGradI( Pi, Pj, PkPlusH, Pl, grad1 );
    getThetaGradI( Pi, Pj, Pk, Pl, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hki.set( j, i, grad1[j] / H );
  }
}

#endif