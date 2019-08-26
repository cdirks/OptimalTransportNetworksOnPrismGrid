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

#ifndef USE_EXTERNAL_OPENMESH
#error OpenMesh external is needed to compile module openmesh
#endif

#include <triMesh.h>
#include <openMeshOps.h>
#include <preconditioner.h>

namespace om {

template <typename RealType>
void TriMesh<RealType>::setVertexToAverageOfNeighbours( const VertexHandle &vh ) {
  std::vector<int> indices;
  get1RingVertices( vh.idx(), indices );
  aol::Vec3<RealType> average, vertex;
  average.setZero();
  for( uint i = 0; i < indices.size(); i++ ){
    getVertex( indices[i], vertex );
    average += vertex;
  }
  average *= 1./indices.size();
  setVertex( vh, average );
}
  
template <typename RealType>
RealType TriMesh<RealType>::enclosedVolume( ) const {
  RealType v = 0.;
  int num = 0;

  for ( ElementIteratorType it = *this; it.notAtEnd(); ++it ) {
    aol::Vec3<RealType> n, c;
    it->barycenter ( c );
    it->normalizedNormal ( n );
    v += ( n * c ) * it->area();
    num++;
  }
  return v / 3.;
}

template <typename RealType>
RealType TriMesh<RealType>::area( ) const {
  RealType a = 0.;
  int num = 0;
  for ( ElementIteratorType it = *this; it.notAtEnd(); ++it ) {
    a += it->area();
    num++;
  }
#ifdef VERBOSE
  std::cerr << "anum = " << num << std::endl;
#endif
  return a;
}

template <typename RealType>
RealType TriMesh<RealType>::getAspectRatio( const FaceHandle &fh ) const{

  // aspect ratio = h / r, where r = radius of biggest inner circle = 2 * vol(fh) / (sum of edge lengths)
  aol::Vec3<RealType> v;
  int i = 0;
  typename OpenMeshType::ConstFaceHalfedgeIter cfh_it = cfh_iter( fh );
  for(; cfh_it.is_valid(); ++cfh_it)
    v[i++]  = getOpenMeshObject().calc_edge_length( *cfh_it );
  RealType ar = max( v[0], max(v[1], v[2]) );
  ar *= (v[0] + v[1] + v[2] );
  getFaceNormal ( fh, v, false );
  return ar / v.norm() ;

}

template <typename RealType>
void TriMesh<RealType>::translate( const aol::Vec3<RealType>& offset ) {
  for ( int i = 0; i < getNumVertices(); i++ ) {
    Point &p = getOpenMeshObject().point ( getOpenMeshObject().vertex_handle ( i ) );
    for( int j = 0; j < 3; j++ )
      p[j] += offset[j];
  }
}

template <typename RealType>
void TriMesh<RealType>::rotate( const aol::Matrix33<RealType>& RotationMatrix ) {
  for ( int i = 0; i < getNumVertices(); i++ ) {
    const Point &p = getOpenMeshObject().point ( getOpenMeshObject().vertex_handle ( i ) );
    aol::Vec3<RealType> vertex ( p[0], p[1], p[2] ), rotatedVertex;
    RotationMatrix.mult( vertex, rotatedVertex );
    setVertex (i, rotatedVertex);
  }
}

template <typename RealType>
void TriMesh<RealType>::scaleSizeByFactor( RealType factor ) {
  for ( int i = 0; i < getNumVertices(); i++ ) {
    Point &p = getOpenMeshObject().point ( getOpenMeshObject().vertex_handle ( i ) );
    for( int j = 0; j < 3; j++ )
      p[j] *= factor;
  }
}

template <typename RealType>
void TriMesh<RealType>::normalizeTo01( bool translate ) {
  const int size = getNumVertices();
  aol::Vec3<RealType> min, max;

  RealType factor = getMeshBoundingBox ( min, max );
  Point offset;
  for ( int i = 0; i < 3; ++i)
    offset[i] = min[i];
#ifdef VERBOSE
  std::cerr << "maximum diam = " << factor << endl;
#endif
  if ( (factor > 1.0) || translate )
    for ( int i = 0; i < size; i++ ) {
      Point &p = getOpenMeshObject().point ( getOpenMeshObject().vertex_handle ( i ) );
      if (translate)
        p -= offset;
      if (factor > 1.0)
        for ( int j = 0; j < 3; j++ )
          p[j] /= factor;
    }
}


template <typename RealType>
void TriMesh<RealType>::getFaceNormal ( const FaceHandle &fh, aol::Vec3<RealType> &wn, bool unit) const  {

  // Dx = [x1-x0 | x2-x0] = [ e_1 | -e_0], e_1 x (-e_0) = e_0 x e_1
  typename OpenMeshType::ConstFaceHalfedgeIter feIter = this->getOpenMeshObject().cfh_iter( fh );
  Point a, b;
  this->getOpenMeshObject().calc_edge_vector( *feIter, a );
  ++feIter;
  this->getOpenMeshObject().calc_edge_vector( *feIter, b );

  for( int i = 0; i < 3; i++ )
    wn[i] = a[(i+1)%3] *  b[(i+2)%3]  -  a[(i+2)%3]  *  b[(i+1)%3];

  if(unit){
    if( wn.normSqr() == 0.){
      cerr<<"Face "<<fh.idx()<<" is degenerated! Edges ";
      ConstFaceHalfedgeIter fh_it = this->getOpenMeshObject().cfh_iter(fh);
      for(; fh_it.is_valid(); ++fh_it)
       cerr<<static_cast<int>( floor( fh_it->idx() / 2. ) )<<", ";
      cerr<<"\n";
      ConstFaceVertexIter fv_it = this->getOpenMeshObject().cfv_iter(fh);
      for(; fv_it.is_valid(); ++fv_it) {
        aol::Vec3<RealType> v;
        getVertex( *fv_it, v);
        cerr<<" v["<<fv_it->idx()<<"] = ("<<v[0]<<", "<<v[1]<<", "<<v[2]<<").\n";
      }
      throw aol::Exception ( "TriMesh<RealType>::getFaceNormal: degenerated face!", __FILE__, __LINE__ );
    }
    wn.normalize();
  }
}


template <typename RealType>
void TriMesh<RealType>::getDx(const FaceHandle &fh,  aol::Mat<3,2,RealType> &Dx ) const {

  // Dx = [x1-x0 | x2-x0] = [ e_1 | -e_0]
  Dx.setZero();
  typename OpenMeshType::ConstFaceHalfedgeIter feIter = this->getOpenMeshObject().cfh_iter( fh );
  for( int i = 0; i < 2; i++ ){
    Point edge;
    this->getOpenMeshObject().calc_edge_vector( *feIter, edge );
    for ( int j = 0; j < 3; j++ )
      Dx.set( j, (i+1)%2, (2.*i - 1.) * edge[j] );
    ++feIter;
  }

}

template <typename RealType>
void TriMesh<RealType>::get1RingVertices( const int idx, std::vector<int>& indices ) const {
  indices.clear();
  for ( ConstVertexVertexIter vv_it = getOpenMeshObject().cvv_iter ( vertexHandle(idx) ); vv_it.is_valid(); ++vv_it )
    indices.push_back(  static_cast<int>( vv_it->idx() ) );
}
  
template <typename RealType>
void TriMesh<RealType>::get2RingVertices( const int idx, std::vector<int>& indices ) const {
    std::set<int> tempSet;
  std::vector<int> oneRingIndices, tempVector;
  get1RingVertices( idx, oneRingIndices );
  
  // get 2-ring indices (without counting them twice, hence using std::set<>)
  for( uint i = 0; i < oneRingIndices.size(); i++ ){
    tempSet.insert( oneRingIndices[i] );
    get1RingVertices( oneRingIndices[i], tempVector );    
    for( uint j = 0; j < tempVector.size(); j++ )
      tempSet.insert( tempVector[j] );
  }
  
  // write back
  indices.clear();
  for ( std::set<int>::iterator it = tempSet.begin(); it != tempSet.end(); ++it)
    indices.push_back( *it );
}
  
  

template <typename RealType>
void TriMesh<RealType>::getFirstFundForm(const FaceHandle &fh, aol::Matrix22<RealType> &g ) const{
  aol::Mat<3,2,RealType> Dx;
  getDx(fh,Dx);
  g.makeProductAtransposedB(Dx,Dx);
}


template <typename RealType>
void TriMesh<RealType>::toVector ( aol::MultiVector<RealType> &dst ) const {
  const int size = getNumVertices();
  if( (dst.numComponents() != 3) || (dst[0].size() != size) )
    dst.reallocate( 3, size );
  for ( int i = 0; i < size; i++ ) {
    const Point &p = this->point ( getOpenMeshObject().vertex_handle ( i ) );
    for ( int j = 0; j < 3; j++ ) {
      dst[j][i] = p[j];
    }
  }
}

template <typename RealType>
void TriMesh<RealType>::toVector ( aol::Vector<RealType> &dst ) const {
  dst.reallocate( 3 * getNumVertices());
  int idx = 0;
  for ( NodeIteratorType vIter = *this; vIter.notAtEnd(); ++vIter ) {
    const Point &p = this->point ( vIter.vertexHandle ( ) );
    for ( int j = 0; j < 3; j++ )
      dst[idx++] = p[j];
  }
}

template <typename RealType>
void TriMesh<RealType>::fromVector ( const aol::MultiVector<RealType> &src ) {

  if (  src.numComponents() != 3 )
    throw aol::Exception ( "source must have 3 components", __FILE__, __LINE__ );

  const int size = getNumVertices();
  if( src[0].size() != size )
    throw aol::Exception ( "sizes do not match", __FILE__, __LINE__ );

  for ( int i = 0; i < size; i++ ) {
    Point &p = this->point ( getOpenMeshObject().vertex_handle ( i ) );
    for ( int j = 0; j < 3; j++ )
      p[j] = src[j][i];
  }
}

template <typename RealType>
void TriMesh<RealType>::fromVector ( const aol::Vector<RealType> &src ) {

  if( src.size() != 3 * getNumVertices() )
    throw aol::Exception ( "sizes do not match", __FILE__, __LINE__ );

  int idx = 0;
  for ( NodeIteratorType vIter = *this; vIter.notAtEnd(); ++vIter  ) {
      Point p;
      for(int i = 0; i < 3; ++i)
        p[i] = src[idx++];
      this->set_point( vIter.vertexHandle() , p );
  }
}


template <typename RealType>
void TriMesh<RealType>::regularizeMesh ( ) {
  for ( NodeIteratorType vIter = *this; vIter.notAtEnd(); ++vIter ) {
    int faces = 0;
    aol::Vec3<RealType> wbary, wn;
    RealType totalvol = 0;

    wbary.setZero();
    wn.setZero();
    for ( VertexFaceIter vf_it = getOpenMeshObject().vf_iter ( vIter.vertexHandle() ); vf_it.is_valid(); ++vf_it ) {

      // compute the barycenter
      aol::Vec3<RealType> bary, coords[3], v[2], n;
      bary.setZero();
      int i = 0;
      for ( FaceVertexIter iv_it = getOpenMeshObject().fv_iter ( *vf_it ); iv_it.is_valid(); ++iv_it ) {
        const Point &p = this->point ( *iv_it );
        for ( int j = 0; j < 3; j++ ) { coords[i][j] = p[j]; }
        bary += coords[i];
        ++i;
      }
      bary /= 3.;
      v[0] = coords[1]; v[0] -= coords[0];
      v[1] = coords[2]; v[1] -= coords[0];
      n = v[0].crossProduct ( v[1] );
      n *= 0.5;
      totalvol += n.norm();
      wn += n;

      bary *= n.norm();
      wbary += bary;
      ++faces;
    }
    wn /= totalvol;
    //wbary /= (RealType)faces;
    wbary /= totalvol;

    // second cycle: compute t
    RealType nom = 0, denom = 0;
    for ( VertexFaceIter vf_it = getOpenMeshObject().vf_iter ( vIter.vertexHandle() ); vf_it.is_valid(); ++vf_it ) {

      // compute the barycenter
      aol::Vec3<RealType> coords[3];
      int center = -1;
      int i = 0;
      for ( FaceVertexIter iv_it = getOpenMeshObject().fv_iter ( *vf_it ); iv_it.is_valid(); ++iv_it ) {
        if ( vIter.vertexHandle() == *iv_it ) {
          center = i;
        }
        ++i;
      }
#ifdef VERBOSE
      std::cerr << "center = " << center << endl;
#endif
      i = 0;
      for ( FaceVertexIter iv_it = getOpenMeshObject().fv_iter ( *vf_it ); iv_it.is_valid(); ++iv_it ) {
        const Point &p = this->point ( *iv_it );
        for ( int j = 0; j < 3; j++ ) { coords[ ( i+3-center ) %3][j] = p[j]; }
        ++i;
      }
      for ( int j = 0; j < 3; j++ ) {
        coords[j] -= wbary;
      }
      aol::Vec3<RealType> cross;
      cross = coords[0].crossProduct ( coords[1] );
      nom += ( cross * coords[2] );
      cross = wn.crossProduct ( coords[1] );
      denom += cross * coords[2];
    }

    RealType t = nom / denom;
    Point &p = this->point ( vIter.vertexHandle() );
    // wn *=  t / (1 + aol::Abs(t*100));
    wn *= 0.8 * t;
    //std::cerr << " t = " << t ;
    wbary += wn;

    for ( int j = 0; j < 3; j++ ) {
      p[j] = wbary[j];
    }
  }
}


template <typename RealType>
RealType TriMesh<RealType>::getMeshBoundingBox ( aol::Vec3<RealType> &min,
                                                    aol::Vec3<RealType> &max ) const {
  const int size = getNumVertices();
  const Point &p = this->point ( getOpenMeshObject().vertex_handle ( 0 ) );
  min[0] = p[0]; min[1] = p[1]; min[2] = p[2];
  max = min;

  for ( int i = 1; i < size; i++ ) {
    const Point &p = getOpenMeshObject().point ( getOpenMeshObject().vertex_handle ( i ) );
    for ( int j = 0; j < 3; j++ ) {
      max[j] = aol::Max ( ( RealType ) p[j], max[j] );
      min[j] = aol::Min ( ( RealType ) p[j], min[j] );
    }
  }
  return aol::Max ( aol::Max ( max[0] - min[0], max[1] - min[1] ),  max[2] - min[2] );
}


template <typename RealType>
void TriMesh<RealType>::loadFromUDPLY ( const std::string filename ) {

  aol::TriangMesh<RealType> triangMesh;
  triangMesh.loadFromUDPLY ( filename.c_str() );
  importFromTriangMesh(triangMesh);
}

template <typename RealType>
void TriMesh<RealType>::loadFromPLY ( const std::string filename ) {

  aol::TriangMesh<RealType> triangMesh;
  triangMesh.loadFromPLY ( filename.c_str() );
  importFromTriangMesh(triangMesh);
}

template <typename RealType>
void TriMesh<RealType>::loadFromOBJ ( const std::string filename, bool mode2D ) {

  aol::TriangMesh<RealType> triangMesh;
  triangMesh.loadFromOBJ ( filename.c_str(), mode2D );
  importFromTriangMesh(triangMesh);
}

// convert this to TriangMesh
template <typename RealType>
void TriMesh<RealType>::convertToTriangMesh ( aol::TriangMesh<RealType> & triangMesh ) const {

  const int numVertex = getNumVertices();
  const int numTriang = getNumFaces();
  triangMesh.resize ( numVertex, numTriang );

  int nIndex = 0;
  for ( NodeIteratorType nIter = *this; nIter.notAtEnd(); ++nIter, ++nIndex ) {
    const Point &p = getOpenMeshObject().point ( nIter.vertexHandle() );

    aol::Vec3<RealType> vertex ( p[0], p[1], p[2] );
    triangMesh.setVertex (nIndex, vertex);
  }

  int nFace = 0;
  for ( ElementIteratorType fIter = *this; fIter.notAtEnd(); ++fIter, ++nFace ) {

    int n = 0;
    aol::Vec3<int> triangle;
    for ( ConstFaceVertexIter v_it = getOpenMeshObject().cfv_iter ( fIter.faceHandle() ); v_it.is_valid(); ++v_it ) {
      triangle[n] = v_it->idx();
      n++;
    }
    triangMesh.setTriang ( nFace, triangle );
  }
}

// import from aol::TriangMesh
template <typename RealType>
TriMesh<RealType>& TriMesh<RealType>::importFromTriangMesh(const aol::TriangMesh<RealType>& triangMesh){
  this->clear();

  const int numVertex = triangMesh.getNumVertices();
  const int numTriang = triangMesh.getNumTriangs();

  for ( int i = 0; i < numVertex; i++ )
    addVertex ( triangMesh.getVertex ( i ) );

  for ( int i = 0; i < numTriang; i++ )
    addFace ( triangMesh.getTriang ( i ) );

  return (*this);

}


// save as OBJ
template <typename RealType>
void TriMesh<RealType>::saveAsOBJ ( const std::string filename, int precision, bool planar ) const {
  
    aol::Bzipofstream out ( filename.c_str() );
    out.precision( precision );

    if ( !out )
      throw ( aol::Exception ( "TriMesh::saveAsOBJ(): could not open file \""
                               + filename + "\" for writing", __FILE__, __LINE__ ) );

    out << "# " << this->getNumVertices () << " vertices, "<<this->getNumFaces () << " faces" <<endl;

    // vertex coordinates and data
    for ( NodeIteratorType it = *this; it.notAtEnd(); ++it ) {
      aol::Vec3<RealType> coords = it.getCoords ();
      out<<"v "<< coords[0] <<" "<<coords[1];
      if( !planar ) out <<" "<<coords[2];
      out<<endl;
    }
    // face date
    for ( ElementIteratorType it = *this; it.notAtEnd(); ++it ) {
      aol::Vec3<int> indices = it.getNodeIndices();
      out<<"f "<< indices[0]+1 <<" "<<indices[1]+1<<" "<<indices[2]+1<<endl;
    }
}

//=======================================================================================================
//  TriMeshWithEdgeNormals
//=======================================================================================================

template <typename RealType>
void TriMeshWithEdgeNormals<RealType>::getEdgeNormal ( const FaceHandle &fh, int localIdx, aol::Vec3<RealType> &wn, bool normalized ) const{
  if( localIdx > 2 )
    throw aol::Exception ( "unvalid local index (must be 0,1,2)!", __FILE__, __LINE__ );

  this->getFaceNormal( fh, wn );

  // get right position of HalfedgeIterator
  typename TriMesh<RealType>::ConstFaceHalfedgeIter feIter = this->getOpenMeshObject().cfh_iter( fh );
  for(int i = 0; i < localIdx; i++)
    ++feIter;

  // if the oppsoite Halfedge is not a boundary Halfedge, add its face normal
  typename TriMeshWithEdgeNormals<RealType>::HalfedgeHandle oh( this->getOpenMeshObject().opposite_halfedge_handle( *feIter ));
  if( !this->getOpenMeshObject().is_boundary( oh ) ){
    aol::Vec3<RealType> temp;
    this->getFaceNormal( this->getOpenMeshObject().face_handle( oh ), temp );
    wn += temp;
  }

  if(normalized){
    if( wn.normSqr() == 0.){
      cerr<<"Face "<<fh.idx()<<" is degenerated:\n";
      typename TriMeshWithEdgeNormals<RealType>::ConstFaceVertexIter fv_it = this->getOpenMeshObject().cfv_iter(fh);
      for(; fv_it.is_valid(); ++fv_it) {
        aol::Vec3<RealType> v;
        this->getVertex( *fv_it, v);
        cerr<<" v["<<fv_it->idx()<<"] = ("<<v[0]<<", "<<v[1]<<", "<<v[2]<<").\n";
      }
      throw aol::Exception ( "TriMeshWithEdgeNormals<RealType>::getEdgeNormal : degenerated face!", __FILE__, __LINE__ );
    }
    wn.normalize();
  }
}


template <typename RealType>
void TriMeshWithEdgeNormals<RealType>::getEdgeNormal ( const EdgeHandle &eh, aol::Vec3<RealType> &wn, bool normalized ) const{
  if(!_edgeNormalProp)
    throw aol::Exception ( "edge normal property has not been set yet!", __FILE__, __LINE__ );
  wn = this->getOpenMeshObject().property( _edgeNormals, eh );
  if( normalized )
    wn.normalize();
}

template <typename RealType>
void TriMeshWithEdgeNormals<RealType>::releaseEdgeNormals (){
  if( _edgeNormalCount > 0 ){
    _edgeNormalCount--;
    if ( _edgeNormalCount == 0 ) {
      this->getOpenMeshObject().remove_property ( _edgeNormals );
      _edgeNormalProp  = false;
    }
  }
}

template <typename RealType>
bool TriMeshWithEdgeNormals<RealType>::hasEdgeNormals() const {
  return _edgeNormalProp;
}

template <typename RealType>
int  TriMeshWithEdgeNormals<RealType>::updateEdgeNormals (){
  requestEdgeNormals ( );
  return _edgeNormalCount;
}


template <typename RealType>
void TriMeshWithEdgeNormals<RealType>::requestEdgeNormals(  ) {

  if(_edgeNormalProp==false)
    this->add_property(_edgeNormals);

  // iterate over all faces
  for ( typename TriMesh<RealType>::ElementIteratorType fIter = *this; fIter.notAtEnd(); ++fIter ) {

    // get face normal
    aol::Vec3<RealType> faceNormal;
    this->getFaceNormal ( fIter.faceHandle(), faceNormal );

    // iterate over edges
    typename TriMesh<RealType>::ConstFaceEdgeIter feIter = this->getOpenMeshObject().cfe_iter( fIter.faceHandle() );
    for(; feIter.is_valid(); ++feIter)
      this->getOpenMeshObject().property( _edgeNormals, *feIter ) += faceNormal;

  }

  _edgeNormalProp = true;
  _edgeNormalCount++;
}

template <typename RealType>
void TriMeshWithEdgeNormals<RealType>::getSecondFundForm( const FaceHandle &fh, aol::Matrix22<RealType> &h ) const{

  // c = 1/(2|\omega|^2) where \omega represents the unit reference element in \R^2 having vertices at (0,0), (1,0), (0,1)
  RealType c = 2.;

  // get edges e_i
  aol::Vec3<RealType> edges[3];
  typename TriMesh<RealType>::Point points[3];
  int i = 0;
  for( typename TriMesh<RealType>::ConstFaceVertexIter fv_it = this->getOpenMeshObject().cfv_iter(fh); fv_it.is_valid(); ++fv_it )
    points[i++] = this->getOpenMeshObject().point ( *fv_it );

  // in Openmesh vertex i has e_i as incoming edge and e_{(i+1)%3} as outgoing edge
  for( int j = 0; j<3; j++ )
   for( int k = 0; k<3; k++ )
     edges[j][k] = points[j][k] - points[ (j+2)%3 ][k];

  // edge normals
  aol::Vec3<RealType> normals[3];
  for( int i = 0; i < 3; i++ )
    getEdgeNormal ( fh, i, normals[i] );

  // h = c \sum_{i=0}^2 ( N_i, e_{i+2) B_i, where N_i is the normal on edge e_i, i = i%3 is the local coordinate
  // rotation of e_i by 90 degrees leads to t_i; set B_i := t_i t_i^T
  // t_0 = (-1, 0) => B_0 = (1 0; 0 0), t_1 = (0 , -1) =>  B_1 = (0 0; 0 1), t_2 = (1, 1) => B_2 = (1 1; 1 1)
  RealType aux = normals[2]*edges[1];
  h.set (0, 0, aux + normals[0]*edges[2]);
  h.set (0, 1, aux);
  h.set (1, 0, aux);
  h.set (1, 1, aux + normals[1]*edges[0]);
  h *= c;

}

template <typename RealType>
void TriMeshWithEdgeNormals<RealType>::getSecondFundFormFast( const FaceHandle &fh, aol::Matrix22<RealType> &h ) const{
  if(!_edgeNormalProp)
    throw aol::Exception ( "edge normal property has not been set yet!", __FILE__, __LINE__ );

  // c = 1/(2|\omega|^2) where \omega represents the unit reference element in \R^2 having vertices at (0,0), (1,0), (0,1)
  RealType c = 2.;

  // get edges e_i
  aol::Vec3<RealType> edges[3];
  typename TriMesh<RealType>::Point points[3];
  int i = 0;
  for( typename TriMesh<RealType>::ConstFaceVertexIter fv_it = this->getOpenMeshObject().cfv_iter(fh); fv_it.is_valid(); ++fv_it )
    points[i++] = this->getOpenMeshObject().point ( *fv_it );

  // in Openmesh vertex i has e_i as incoming edge and e_{(i+1)%3} as outgoing edge
  for( int j = 0; j<3; j++ )
   for( int k = 0; k<3; k++ )
     edges[j][k] = points[j][k] - points[ (j+2)%3 ][k];

  // edge normals
  aol::Vec3<RealType> normals[3];
  i = 0;
  for( typename TriMesh<RealType>::ConstFaceEdgeIter fe_it = this->getOpenMeshObject().cfe_iter(fh); fe_it.is_valid(); ++fe_it )
    getEdgeNormal ( *fe_it, normals[i++] );

  // h = c \sum_{i=0}^2 ( N_i, e_{\pi(i)} ) B_i, where N_i is the normal on edge e_i, \pi(i):= (i+2)%3, i is the local coordinate
  // rotation of e_i by 90 degrees leads to t_i; set B_i := t_i t_i^T
  // t_0 = (0, -1) => B_0 = (1 0; 0 0), t_1 = (1 , 1) =>  B_1 = (1 1; 1 1), t_2 = (-1, 0) => B_2 = (0 0; 0 1)
  RealType aux = normals[1]*edges[0];
  h.set (0, 0, aux + normals[0]*edges[2]);
  h.set (0, 1, aux);
  h.set (1, 0, aux);
  h.set (1, 1, aux + normals[2]*edges[1]);
  h *= c;
}

template <typename RealType>
void TriMeshWithEdgeNormals<RealType>::getShapeOperatorFast( const FaceHandle &fh, aol::Matrix22<RealType> &S ) const {
    aol::Matrix22<RealType> temp;
    this->getFirstFundForm(fh,temp);
    S.makeInverse(temp);
    getSecondFundFormFast(fh,temp);
    S*=temp;
}

//=======================================================================================================
//=======================================================================================================

template <typename RealType>
void TriMeshWithVertexNormals<RealType>::getVertexNormal ( const VertexHandle &vh, aol::Vec3<RealType> &wn, bool unit) const {

  wn.setZero();
  int faces = 0;
  for ( typename TriMesh<RealType>::ConstVertexFaceIter vf_it = this->getOpenMeshObject().cvf_iter ( vh ); vf_it.is_valid(); ++vf_it ) {
    // compute the barycenter
    aol::Vec3<RealType> coords[3], v[2], n;
    int i = 0;
    for ( typename TriMesh<RealType>::ConstFaceVertexIter iv_it = this->getOpenMeshObject().cfv_iter ( *vf_it ); iv_it.is_valid(); ++iv_it ) {
      const typename TriMesh<RealType>::Point &p = this->getOpenMeshObject().point ( *iv_it );
      for ( int j = 0; j < 3; j++ ) { coords[i][j] = p[j]; }
      ++i;
    }
    v[0] = coords[1]; v[0] -= coords[0];
    v[1] = coords[2]; v[1] -= coords[0];
    n = v[0].crossProduct ( v[1] );
    n.normalize();
    wn += n;

    ++faces;
  }

  if(unit)
    wn.normalize();
  else
    wn /= static_cast<RealType>(faces);
}

template <typename RealType>
void TriMeshWithVertexNormals<RealType>::getDn(const FaceHandle &fh, aol::Mat<3,2,RealType> &Dn ) const {

  Dn.setZero();
  aol::Vec3<RealType> wn[3];
  int i =0;
  for( typename TriMesh<RealType>::ConstFaceVertexIter fv_it = this->getOpenMeshObject().cfv_iter(fh); fv_it.is_valid(); ++fv_it )
    getVertexNormal(*fv_it,wn[i++]);

  // Dn = [n1-n0 | n2-n0]
  aol::Vec<2,RealType> row;
  for(int i =0; i<3; ++i){
    row[0] = wn[1][i] - wn[0][i];
    row[1] = wn[2][i] - wn[0][i];
    Dn[i]  = row;
  }
}

template <typename RealType>
void TriMeshWithVertexNormals<RealType>::getSecondFundForm( const FaceHandle &fh, aol::Matrix22<RealType> &h ) const{
  aol::Mat<3,2,RealType> Dx, Dn;
  this->getDx(fh,Dx);
  getDn(fh,Dn);
  h.makeSymmetricProduct(Dn,Dx);
}

template <typename RealType>
void TriMeshWithVertexNormals<RealType>::noise ( RealType factor, bool fixBoundary ) {
  for ( typename TriMesh<RealType>::NodeIteratorType vIter = *this; vIter.notAtEnd(); ++vIter ) {
    if( fixBoundary & vIter.isAtBoundary() )
      continue;

    aol::Vec3<RealType> wn;
    getVertexNormal ( vIter.vertexHandle(), wn );

    typename TriMesh<RealType>::Point &p = this->point ( vIter.vertexHandle() );
    RealType r = ( static_cast<RealType> ( rand() ) / static_cast<RealType> ( RAND_MAX ) * 2. - 1. ) * factor;
    for ( int j = 0;j < 3;j++ ) {
      p[j] += r * wn[j];
    }
  }
}
//=======================================================================================================
//=======================================================================================================

template <typename RealType, typename Imp>
void TriMeshWithGeomProps<RealType, Imp>::getCurvature( const FaceHandle &fh, aol::Vec2<RealType> &curv ) const{
  aol::Matrix22<RealType> shapeOp;
  getShapeOperator(fh, shapeOp );
  shapeOp.eigenValues(curv);
}

template <typename RealType, typename Imp>
void TriMeshWithGeomProps<RealType, Imp>::writePrincipalCurvaturesTo ( aol::MultiVector<RealType> & pc ) {

  if ( pc.numComponents() != 2 || ! pc.allDimsEqual ( this->getNumFaces() ) )
    throw aol::DimensionMismatchException ( "om::TrimMesh::writePrincipalCurvaturesTo(): "
                                            "passed MultiVector has not the right size.",
                                            __FILE__, __LINE__ );

  aol::Vec2<RealType> temp;
  for ( typename TriMesh<RealType>::ElementIteratorType fIter = *this; fIter.notAtEnd(); ++fIter ) {

    getCurvature( fIter.faceHandle(), temp );

    pc[0][fIter.getIndex()] = temp[0];
    pc[1][fIter.getIndex()] = temp[1];
  }
}

//! exports mesh with curvatures as vtk-file
template <typename RealType, typename Imp>
void TriMeshWithGeomProps<RealType, Imp>::exportCurvature( const std::string constfilename ){

  int n = this->getNumFaces();
  aol::MultiVector<RealType> curv ( 2, n );
  writePrincipalCurvaturesTo ( curv );

  aol::Vector<RealType> meanCurv  ( n ),
                        gaussCurv ( n ),
                        big20     ( n ),
                        big50     ( n ),
                        big100    ( n );

  for ( int i = 0; i < n; ++i ) {
    meanCurv[i]  = curv[0][i] + curv[1][i];
    gaussCurv[i] = curv[0][i] * curv[1][i];
    big20[i]     = std::max(curv[0][i], curv[1][i])> 20  ? 1 : 0;
    big50[i]     = std::max(curv[0][i], curv[1][i])> 50  ? 1 : 0;
    big100[i]    = std::max(curv[0][i], curv[1][i])> 100 ? 1 : 0;
  }

  // Does file end with ".vtk"?
  string filename = constfilename.substr( 0, constfilename.rfind(".") ) + ".vtk";

  aol::MeshWithData<TriMesh<RealType> > ( *this )
    .addData ( curv[0],   "k1", aol::FACE_DATA )
    .addData ( curv[1],   "k2", aol::FACE_DATA )
    .addData ( gaussCurv, "K",  aol::FACE_DATA )
    .addData ( meanCurv,  "H",  aol::FACE_DATA )
    .addData ( big20,     "big20",   aol::FACE_DATA )
    .addData ( big50,     "big50",   aol::FACE_DATA )
    .addData ( big100,    "big100",  aol::FACE_DATA )
    .saveAsLegacyVTK ( filename );
}

//=======================================================================================================
//=======================================================================================================

template class TriMesh<float>;
template class TriMesh<double>;
template class TriMesh<long double>;

template class TriMeshWithVertexNormals<float>;
template class TriMeshWithVertexNormals<double>;
template class TriMeshWithVertexNormals<long double>;

template class TriMeshWithEdgeNormals<float>;
template class TriMeshWithEdgeNormals<double>;
template class TriMeshWithEdgeNormals<long double>;

template class TriMeshWithGeomProps< float,  TriMeshWithEdgeNormals<float> >;
template class TriMeshWithGeomProps< float,  TriMeshWithVertexNormals<float> >;
template class TriMeshWithGeomProps< double, TriMeshWithEdgeNormals<double> >;
template class TriMeshWithGeomProps< double, TriMeshWithVertexNormals<double> >;

} // end of namespace om.
