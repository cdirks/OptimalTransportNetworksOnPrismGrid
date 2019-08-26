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

#ifndef __PROGRESSIVEMESH_H
#define __PROGRESSIVEMESH_H

#include <openMeshIncludes.h>

#include <aol.h>
#include <configurators.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <triangMesh.h>
#include<progressBar.h>


namespace om {

/**
 * This class describes the collapse of a single edge e = ( vt, vs ) to one vertex vs_new.
 * The collapse is performed in all of the N meshes that are decimated simultaneously ( cf. describtion of SimulProgMesh<> ).
 * The index of the vertex vs in the original mesh(es) is stored in vs_orig.
 * The displacements delta_vs = vs_new - vs (and delta_vt = vs_new - vt) are stored as well for every mesh to perform the prolongation afterwards.
 *
 * \author Heeren
 */
template <typename RealType >
class Collapse{

public:
Collapse( int N ) : delta_vs( N ), delta_vt( N ){}

Collapse( const Collapse& other ) : vsOrig( other.vsOrig ) {
  delta_vs.clear();
  delta_vt.clear();
  for( int i = 0; i < other.delta_vs.size(); i++ ){
    delta_vs.pushBack( other.delta_vs[i] );
    delta_vt.pushBack( other.delta_vt[i] );
  }
}

// index of vs in original meshes (ie. the meshes that are commited in SimulProgMesh<>::calculateDecimation () )
int vsOrig;
// displacements (needed by prolongation)
aol::RandomAccessContainer< aol::Vec3<RealType> > delta_vs;
aol::RandomAccessContainer< aol::Vec3<RealType> > delta_vt;

};


/**
 * The class ProlongInfo::RelativeCoords stores the relative coordiantes of a point x \in \R^3
 * with respect to the closest triangle T = (t_0, t_1, t_2), t_i \in \R^3, in a triangular mesh.
 * The index of T is stored in 'face', the barycentric coordinates of the projection P(x) are given by bar0 and bar1,
 * projN = ( x - P(x), N ), where N is the weighted normal on T.
 * Note: P is the orthogonal projection on the plane that is uniquely defined by T. Hence bar0 and bar1 are not barycentric coordinates
 * in the sense that 0 <= bar0, bar1 <= 1, bar0 + bar1 <= 1, but the local coord. of P(x) on the plane:
 * P(x) = t_0 + bar0*( t_1 - t_0 ) + bar1*( t_2 - t_0 )
 *
 * ProlongInfos stores the relative coordinates of all vertices vt and vs that has been deleted/reset after the collabation of e = (vt,vs)
 * Let M_0 be an initial mesh (with n vertices) and M_k a decimated mesh that has been created out of M_0 by k edge collabations e -> vs.
 * Each vertex 0 <= i <= n-k-1 in M_k has a corresponding index in M_0, which is stored in vertexMap[i]
 * The index of the vertex vt that was deleted in the jth decimation, 1<=j<=k, is stored in vertexMap[ n-k-1+j ].
 * The index of the vertex vs with respect to the initial mesh M_0 is stored in the std::vector vs.
 *
 * \author Heeren
 */
template<typename RealType>
class ProlongInfo{

public:
ProlongInfo(){}
ProlongInfo( const ProlongInfo<RealType>& other )
                       : relCoordsCont( other.relCoordsCont.size() ),
                         vertexMap( other.vertexMap.size() ),
                         vs( other.vs.size() ) {
  for( int i = 0; i < relCoordsCont.size(); i++ )
    relCoordsCont[i] = other.relCoordsCont[i];
  for( int i = 0; i < vertexMap.size(); i++ )
    vertexMap[i] = other.vertexMap[i];
  for( int i = 0; i < vs.size(); i++ )
    vs[i] = other.vs[i];
}

class RelativeCoords{
public:
  RelativeCoords( ) : face(-1), bar0(0.), bar1(0.), projN(0.) {}
  RelativeCoords( const RelativeCoords& other ) : face(other.face), bar0(other.bar0), bar1(other.bar1),  projN(other.projN) {}
  RelativeCoords( int f, RealType b0, RealType b1, RealType pN ) : face(f), bar0(b0), bar1(b1), projN(pN) {}

  int face;
  RealType bar0, bar1;
  RealType projN;

};

void clear(){
  relCoordsCont.clear();
  vertexMap.clear();
  vs.clear();
}

// store all relative coordinates
aol::RandomAccessContainer< RelativeCoords > relCoordsCont;
// store vertexmap from DecimationInfo
std::vector<int> vertexMap;
// store original indices of deleted vertices (vsOrig in Collapse)
// Mind the order: the first entry in vs is the index of the vertex, that has been deleted first!
std::vector<int> vs;
};


/**
 * The class DecimationInfo stores the mesh decimation of an initial mesh M_0 to a coarser mesh M_k.
 * The decimation has been performed by k edge collapses e=(vt,vs) to vs
 * The indices of the collapsed edges are stores in colIndices, detailed describtion of the collapse is stored in collapses.
 * vertexMap stores the indices of the deleted vertices vt as well as the indices of all vertices in M_k w.r.t. M_0 (see above)
 *
 * \author Heeren
 */
template <typename RealType >
class DecimationInfo{

public:
DecimationInfo(){}
DecimationInfo( int n ) : vertexMap(n), vertexProjection(n) {
  for( int i = 0; i < n; i++ ){
    vertexMap[i] = i;
    vertexProjection[i] = i;
  }
}

void reallocate( int n ){
  vertexMap.clear();
  vertexProjection.clear();
  colIndices.clear();
  collapses.clear();
  for( int i = 0; i < n; i++ ){
    vertexMap.push_back( i );
    vertexProjection.push_back( i );
  }
}

aol::RandomAccessContainer< Collapse<RealType> > collapses;
// After k < n decimations vertexMap[i] is  a) the index of vertex i in the original mesh for i < n-k
//                                          b) the index of the vertex vt, that was removed in the jth decimation, for i = n-k+j, j=0,...,k-1
std::vector<int> vertexMap;
// vertexProjection[i] is the vertex index in the coarse mesh, whose position is (likely to be the) closest to vertex i in the original mesh
std::vector<int> vertexProjection;
// stores all indices of the removed edges
std::vector<int> colIndices;
};


/**
 * The class EdgeWithValue stores the index of the edge e =(vt,vs), the energy value that describes the cost of deleting this edge
 * and the optimal position of vs = v_bar after collapsing e -> vs
 *
 * \author Heeren
 */
template < typename RealType >
class EdgeWithValue{

public:
EdgeWithValue( ) : index( -1 ), legal( false ), energy( -1. ) {}
EdgeWithValue( int idx, RealType e, const aol::MultiVector<RealType>& newpos ) : index( idx ), legal( true ), energy( e ), v_bar( newpos ) {}
EdgeWithValue( const EdgeWithValue& other ) : index( other.index ), legal( other.legal ), energy( other.energy ), v_bar( other.v_bar ) {}

EdgeWithValue<RealType> &operator= ( const EdgeWithValue<RealType> &other) {
  index = other.index;
  legal = other.legal;
  energy = other.energy;
  v_bar.reallocate( other.v_bar );
  v_bar = other.v_bar;
  return *this;
}

int index;
bool legal;
RealType energy;
aol::MultiVector<RealType> v_bar;

};

//! class that only exists of a pointer to an EdgeWithValue
template < typename RealType >
class EdgeWithValuePointer{

public:
EdgeWithValuePointer( const EdgeWithValuePointer<RealType>& other ) : pointer( other.pointer ){}
EdgeWithValuePointer( const EdgeWithValue<RealType>& other ) : pointer( &other ){}

const EdgeWithValue<RealType>* pointer;
};

//! comparison functions in order to define std::list< EdgeWithValuePointer<> > --
template <typename RealType>
bool operator<( const EdgeWithValuePointer < RealType > & e1,
                const EdgeWithValuePointer < RealType > & e2 ){ return e1.pointer->energy < e2.pointer->energy; }

template <typename RealType>
bool operator>( const EdgeWithValuePointer < RealType > & e1,
                const EdgeWithValuePointer < RealType > & e2 ){ return e1.pointer->energy > e2.pointer->energy; }

template <typename RealType>
bool operator==( const EdgeWithValuePointer < RealType > & e1,
                 const EdgeWithValuePointer < RealType > & e2 ){ return e1.pointer->energy == e2.pointer->energy; }

template <typename RealType>
bool operator!=( const EdgeWithValuePointer < RealType > & e1,
                 const EdgeWithValuePointer < RealType > & e2 ){ return e1.pointer->energy != e2.pointer->energy; }



/**
 * Projection of a point x to the triangle T = (t_0, t_1, t_2), t_i \in \R^3
 * Partly taken from "Distance Between Point and Triangle in 3D" by David Eberly, Geometric Tools, LLC,
 * http://www.geometrictools.com/
 *
 * If barCoordWRTwholeSubspace = false, then the barycentric coord. 0 <= s,t <= 1 with s + t <= 1 of the closest point to x in T are returned.
 * Othwerwise, the local coordiantes (s,t) \in \R^2 of the orthogonal projection P(x) are returned,
 * where x is projected onto the plane that is uniquely defined by T.
 *
 * \author Heeren
 */
template<typename RealType, typename Point>
RealType computeProjectionOnTriangle( const Point& t0,
                            const Point& t1,
                            const Point& t2,
                            const aol::Vec3<RealType>& x,
                            aol::Vec2<RealType>& barCoord,
                            aol::Vec3<RealType>& projection,
                            bool barCoordWRTwholeSubspace = false ){

    // The triangle T = (t0, t1, t2) is given in barycentric coordinates by:
    // T(s,t) = t0 + s*(t1-t0) + t*(t2-t0),
    // The distance function to X is given by Q(s,t) = | T(s,t) - X |^2
    // Q(s,t) = a*s^2 + 2*b*s*t + c*t^2 + 2*d*s + 2*e*t+ f
    aol::Vec3<RealType> D( t0[0] - x[0], t0[1] - x[1], t0[2] - x[2] );
    Point e0 = t1 - t0, e1 = t2 - t0;
    RealType a ( e0[0]*e0[0] + e0[1]*e0[1] + e0[2]*e0[2] ),
       b ( e0[0]*e1[0] + e0[1]*e1[1] + e0[2]*e1[2] ),
           c ( e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2] ),
       d ( D[0]*e0[0] + D[1]*e0[1] + D[2]*e0[2] ),
       e ( D[0]*e1[0] + D[1]*e1[1] + D[2]*e1[2] );

    // Q'(s,t) = 2*( as+bt+d, bs+ct+e ) = ( 0, 0 ) <==> s = (b*e - c*d) / det , t = (b*d - a*e) / det
    // Note det = a*c-b*b = |e0 x e1|^2 = |normal(T)|^2 > 0
    RealType det( a*c - b*b ), s( b*e - c*d ), t( b*d - a*e );

/*
    // Numbering of the 7 regions, where the argmin could be
    //       ^t
.   //   \ 2 \
    //    \  \
.   //     \ \
.   //      \\
.   //       \                                \
.   //       \ \      1                       \ \
    //       \  \                             \  \  line 0
    //   3   \   \                            \   \
    //       \    \                   line 2  \    \
    //       \  0  \                          \     \
    //       \      \                         \      \
    // ----------------------------> s        --------
    //   4   \    5   \    6                   line 1
    //       \         \
*/

    // determine region and line
    int region( -1 ), line( -1 );
    if ( s + t <= det ){
      if ( s < 0. ){
  if ( t < 0. ) { region = 4; }
  else { line = 2; } // region = 3
      }
      else{
        if ( t < 0. ) { line = 1; } // region = 5
  else { barCoord[0] = s/det; barCoord[1] = t/det; }  // region = 0
      }
    }
    else{
      if ( s < 0. ) { region = 2; }
      else{
  if ( t < 0. ) { region = 6; }
  else { line = 0; }  // region = 1
      }
    }

    if( (line == -1) & (region > 0) ){
      switch( region ){
  // Grad(Q) = 2(as+bt+d,bs+ct+e)
  // (1,-1)*Grad(Q(0,1)) = (1,-1)*(b+d,c+e) = (b+d)-(c+e)
  // min on edge s+t=1 if (1,-1)*Grad(Q(0,1)) < 0 )
  // min on edge s=0 otherwise
  case 2:
      if ( c + e > b + d) // minimum on edge s+t=1
        line = 0;
    else // minimum on edge s=0
            line = 2;
    break;

  // Grad(Q) = 2(as+bt+d,bs+ct+e), Grad(Q(0,0)) = 2(d,e)
  // (1,0)*Grad(Q(0,0)) = (1,0)*(d,e) = d
  // min on edge t=0 if (1,0)*Grad(Q(0,0)) < 0 )
  // min on edge s=0 otherwise
  case 4:
    if ( d < 0. ) // minimum on edge t=0 (cf case 5)
      line = 1;
    else // minimum on edge s=0 (cf case 3)
      line = 2;
    break;

  // Grad(Q) = 2(as+bt+d,bs+ct+e), Grad( Q(1,0)) = (a+d, b+e)
  // (-1,1)*Grad(Q(1,0)) = (-1,1)*(a+d,b+e) = -(a+d)+(b+e)
  // min on edge s+t=1 if (1,-1)*Grad(Q(0,1)) < 0 )
  // min on edge s=0 otherwise
  case 6:
      if ( a + d > b + e ) // minimum on edge s+t=1
        line = 0;
    else// minimum on edge t=0
            line = 1;
    break;

  default:
    cerr<<"Error while calculating projection.\n";
    break;
      }
    }

    if( line != -1 ){
      switch( line ){

  // min on edge s+t=1
  // Intersect graph of Q with plane { (s,t) : t + s = 1 }:
  // F(s) = Q(s,1-s) = (a-2b+c)s^2 + 2(b-c+d-e)s + (c+2e+f)
  // F’(s)/2 = (a-2b+c)s + (b-c+d-e)
  // F’(s) = 0 when s = (c+e-b-d)/(a-2b+c)
  // a-2b+c = |edge1 - edge2|^2 > 0, so only sign of c+e-b-d need to be considered
  case 0: {
    RealType numer = c+e-b-d;
    if ( numer <= 0. )
      barCoord[0] = 0.;
    else{
      RealType denom = a - 2.*b + c; // positive quantity
      barCoord[0] = ( numer >= denom ? 1. : numer/denom );
    }
    barCoord[1] = 1. - barCoord[0];}
          break;

  // min on edge t=0
  // Intersect graph of Q with plane { (s,t) : t = 0 }:
  // F(s) = Q(s,0) = as^2 + 2ds + f
  // F’(s)/2 = as + d
  // F’(s) = 0 when s = -d/a where a = |edge0|^2 > 0
  case 1:
    barCoord[1] = 0.;
          barCoord[0] = d > 0. ? 0. : ( -d > a ? 1. : -d/a );
    break;

  // min on edge s=0
  // Intersect graph of Q with plane { (s,t) : s = 0 }:
  // F(t) = Q(0,t) = ct^2 + 2et + f
  // F’(t)/2 = ct + e
  // F’(t) = 0 when t = -e/c where c = |edge1|^2 > 0
  case 2:
    barCoord[0] = 0.;
          barCoord[1] = e > 0. ? 0. : ( -e > c ? 1. : -e/c );
    break;

  default:
    cerr<<"Error while calculating projection.\n";
    break;
      }
    }

    // return euclidean distance
    RealType dist = 0.;
    for( int i = 0; i < 3; i++ ){
      projection[i] = t0[i] + barCoord[0]*e0[i] + barCoord[1] * e1[i];
      dist += (projection[i] - x[i]) * (projection[i] - x[i]);
    }

    // compute "barycentric" coordinates of the projection of x onto the plane that is defined by T
    if( barCoordWRTwholeSubspace ){
      barCoord[0] = s/det;
      barCoord[1] = t/det;
      for( int i = 0; i < 3; i++ )
        projection[i] = t0[i] + barCoord[0]*e0[i] + barCoord[1]*e1[i];
    }

    return sqrt( dist );
}


/**
 *   Implementation of Simultaneous Progressive Meshes (based on edge decimation).
 *
 *   PROGESSIVE MESHES:
 *   Given: an initial mesh M_0 and an integer k (which represents the number of vertices to be removed).
 *   Aim: Reduce the mesh iterively by collapsing k edges e = (vt,vs) to a single vertex vs_new
 *   Proceeding:
 *   i)   For each edge we calculate a specific cost (which depends on the definition of the error functional) and an optimal position for vs_new.
 *        The triple ( edge index, cost, vs_new) is inserted into a heap
 *   ii)  For i = 1, ..., k:
 *          iia) Sort the heap w.r.t. the costs and  get the cheapest edge e = (vt,vs)
 *          iib) Collapse e to vs_new and update the costs of all edges in a neighbourhood of e
 *
 *   For further information see Hoppe, "Progressive meshes", in SIGGRAPH 96
 *   The energy/error functional is taken from Garland and Heckbert, "Surface Simplification Using Quadric Error Metrics"
 *
 *   SIMULTANEOUS DECIMATION:
 *   Here n meshes (which are topologically equivalent) can be decimated simultaneously.
 *   This means there the energy functional is defined for each mesh and the cost of removing a
 *   particular edge in all meshes simultaneously is given by the sum of all individual energy functionals.
 *
 *   The decimation is performed by calcDecimation(), the (inverse) prolongation is computed by  calcProlongation().
 *
 *   Note: At the moment, boundary edges cannot be removed. ( which does not matter for meshes without boundary, of course. )
 *
 * \author Heeren
 */

template < typename TriMeshType, typename RealType >
class SimulProgMesh{

  typedef typename TriMeshType::Edge      Edge;
  typedef typename TriMeshType::Halfedge  Halfedge;
  typedef typename TriMeshType::EdgeHandle  EdgeHandle;
  typedef typename TriMeshType::HalfedgeHandle  HalfedgeHandle;
  typedef typename TriMeshType::Point Point;
  typedef typename TriMeshType::FaceHandle  FaceHandle;

  typedef typename TriMeshType::VertexVertexIter VViter;
  typedef typename TriMeshType::ConstVertexVertexIter CVViter;
  typedef typename TriMeshType::VertexOHalfedgeIter VOHiter;
  typedef typename TriMeshType::ConstVertexOHalfedgeIter CVOHiter;
  typedef typename TriMeshType::ConstEdgeIter CEiter;
  typedef typename TriMeshType::ConstVertexFaceIter CVFiter;


  typedef EdgeWithValue< RealType > EdgeType;
  typedef EdgeWithValuePointer< RealType > EdgePointerType;
  typedef std::list< EdgePointerType > HeapType;
  typedef std::map< int, EdgeType > EdgeEnergyMap;
  typedef aol::RandomAccessContainer< aol::MultiVector<RealType> > MatrixContainer;
  typedef aol::RandomAccessContainer< TriMeshType > TriMeshContainer;

  typedef typename ProlongInfo<RealType>::RelativeCoords RelativeCoords;

  const RealType Inf;

public:
  SimulProgMesh( ) : Inf( 1e+10 ) {}

  //! decimation
  void calcDecimation ( TriMeshContainer& Trimeshes, RealType theta, DecimationInfo<RealType>& decInfo ) {
    if( !(theta < 1.) )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::calcDecimation(): invalid theta", __FILE__, __LINE__ );
    calcDecimation ( Trimeshes, static_cast<int>( floor( Trimeshes[0].getNumVertices() * (1. - theta)) ), decInfo );
  }

  //! decimation
  void calcDecimation ( TriMeshContainer& Trimeshes, int n, DecimationInfo<RealType>& decInfo ) {

    if( Trimeshes.size() == 0 )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::calcDecimation(): TriMeshContainer is empty!", __FILE__, __LINE__ );
    int numVertices = Trimeshes[0].getNumVertices();
    int currentNumEdges = Trimeshes[0].n_edges();

    decInfo.reallocate( numVertices );

    HeapType heap;
    EdgeEnergyMap emap;
    MatrixContainer matCont;

    // store indices of all old vertices (i.e. preimages) that are merged into the new vertices in the coarse mesh in a map
    // This eases the calculation of the projections in the prolongation
    std::map<int, std::vector<int> > preImages;
    for( int i = 0; i < numVertices; i++ ){
      std::vector<int> preImage;
      preImage.push_back( i );
      preImages[i] = preImage;
    }

    // fill heap and sort
    fillHeap( Trimeshes, matCont, heap, emap );
    heap.sort();
    bool decimationStillPossible = (heap.size() > 0);

    int counter = 0;
    // start progressbar
    aol::ProgressBar<> pb ( "Decimation " );
    pb.start ( n );
    // start decimation
    while( (counter < n) & decimationStillPossible ){

      // get edge with lowest cost
      const EdgeType* edge( heap.front().pointer );
      if( !(edge->energy < Inf) ){
         cerr<<"Minimum energy exceeds "<<Inf<<" after "<<counter<<" decimations!\n\n";
         decimationStillPossible = false;
         break;
      }

      // define a Halfedge on the edge, that is not on the boundary
      int auxI = Trimeshes[0].is_boundary( HalfedgeHandle( 2 * edge->index ) )? 1 : 0;
      HalfedgeHandle heh( 2 * edge->index + auxI );
      
      if( Trimeshes[0].is_boundary(heh) )
	throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::calcDecimation(): handle is on boundary!", __FILE__, __LINE__ );
      for( int j = 0; j < Trimeshes.size(); j++ )
        if( Trimeshes[j].face_handle(heh) == TriMeshType::InvalidFaceHandle )
	  throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::calcDecimation(): invalid face handle!", __FILE__, __LINE__ );
	
      // check whether collapse is ok
      if( !Trimeshes[0].is_collapse_ok( heh ) ){
	//cerr<<"Collapse of edge " << edge->index << " not ok! Proceed...\n";
	// pop critical edge at front and insert at back with energy = INF
	EdgeType newEdge( *edge );
	newEdge.energy = Inf;
	heap.pop_front();
	emap[ newEdge.index ] = newEdge;
        heap.push_back( EdgePointerType( emap[ newEdge.index ] ) );
	continue;
      }
      
       // store index of collapsed edge
      decInfo.colIndices.push_back( edge->index );

      // define relevant halfedges in the neighbourhood
      //HalfedgeHandle nh  = Trimeshes[0].next_halfedge_handle(heh);
      HalfedgeHandle oh  = Trimeshes[0].opposite_halfedge_handle(heh);
      HalfedgeHandle noh = Trimeshes[0].next_halfedge_handle( Trimeshes[0].opposite_halfedge_handle(heh) );

      // these are the edges which are going to be removed
      int colIdx  = edge->index;
      int prevIdx = static_cast<int>( floor( Trimeshes[0].prev_halfedge_handle(heh).idx() / 2 ) );
      int opNextIdx = Trimeshes[0].is_boundary( oh ) ? -1 : static_cast<int>( floor( noh.idx() / 2 ) );

      // these are the faces that are going to be removed (wegen dem OM index swapping müssen diese neuberechnet werden!)
      int faceIdx[2];
      faceIdx[0] = Trimeshes[0].face_handle( heh ).idx();
      faceIdx[1] = Trimeshes[0].is_boundary( oh ) ? -1 : Trimeshes[0].face_handle( oh ).idx();

      // define collapses
      Collapse<RealType> col( Trimeshes.size() );
      int vs = Trimeshes[0].to_vertex_handle( heh ).idx();
      int vt = Trimeshes[0].from_vertex_handle( heh ).idx();
      col.vsOrig = decInfo.vertexMap[ vs ];
      
      // If it happens that vs is the last vertex in the current mesh, then vs_new will recieve the index of vt,
      // since vt is removed after decimation and OpenMesh performs a corresponding index swapping
      vs = ( vs == Trimeshes[0].getNumVertices() - 1 ) ? vt : vs;

      // get the new position for vs_new and store the corresponding displacements vt -> vs_new, vs -> vs_new
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for( int j = 0; j < Trimeshes.size(); j++ ){
        aol::Vec3<RealType> v_new ( edge->v_bar[j][0], edge->v_bar[j][1], edge->v_bar[j][2] );
        Point vs_coord ( Trimeshes[j].point( Trimeshes[j].to_vertex_handle(heh) ) );
        Point vt_coord ( Trimeshes[j].point( Trimeshes[j].from_vertex_handle(heh) ) );
        // delta_vt = vs_new - vt, delta_vs = vs_new - vs
        // => vs = vs_new - delta_vs, vt = vs_new + delta_vs
        aol::Vec3<RealType> delta_vs, delta_vt;
        for( int i = 0; i < 3; i++ ){
          delta_vs[i] = v_new[i] - vs_coord[i];
          delta_vt[i] = v_new[i] - vt_coord[i];
        }
        col.delta_vs[j] = delta_vs;
        col.delta_vt[j] = delta_vt;

        // remove edge and swap indices
        Trimeshes[j].collapse( heh );
        Trimeshes[j].garbage_collection();
        // set vs to its new position
        Trimeshes[j].setVertex( vs, v_new );
      }

      // update number of edges
      currentNumEdges -= 2;
      if( opNextIdx != -1 )
        currentNumEdges -= 1;

      counter++;
      // update preimages: vt gets all preimages of the last vertex which was deleted after decimation
      if( vt != numVertices - counter ){
        for( unsigned int i = 0; i < preImages[numVertices - counter].size(); i++ )
          preImages[vt].push_back( preImages[numVertices - counter][i] );
        preImages[numVertices - counter].clear();
      }
      // update preimages: if vs was not the last index, it gets all preimages of vt
      if( vt != vs ){
        for( unsigned int i = 0; i < preImages[vt].size(); i++ )
          preImages[vs].push_back( preImages[vt][i] );
        preImages[vt].clear();
      }

      // update vertexmap
      swap( decInfo.vertexMap[ vt ], decInfo.vertexMap[ numVertices - counter ] );

      // add the collapse to decimationInfo
      decInfo.collapses.pushBack( col );

      // update the heap
      if( (counter < n) & decimationStillPossible ){
        // Note: For performance reasons, all edges that have to be re-computed are first removed
        // from the heap and then inserted again. To be removed, their costs are set to Inf.

        // write all edge indices whose energy values have to be re-computed to newEdgesCandidates
        // Note: these might be counted twice!
        std::list<int> newEdgesCandidates;
        // Openmesh swaps the edge indices of the deleted edges with the 3 (or 2) last edge indices.
        // Then the last indices are removed. Hence the edges which carry now the "deleted" indices have to be re-computed.
        newEdgesCandidates.push_back( colIdx );
        emap[ colIdx ].energy = Inf;
        emap[ currentNumEdges ].energy = Inf;

        newEdgesCandidates.push_back( prevIdx );
        emap[ prevIdx ].energy = Inf;
        emap[ currentNumEdges + 1 ].energy = Inf;

        if( opNextIdx != -1 ){
          newEdgesCandidates.push_back( opNextIdx );
          emap[ opNextIdx ].energy = Inf;
          emap[ currentNumEdges + 2 ].energy = Inf;
        }

        // Furthermore all edges in a local neighbourhood around the new vertex vs have to be re-computed.
        // However, the particular size of that local neighbourhood depends on the respective energy method
        // All edges that own a vertex which is in the vertex 1-ring of vs as well as their successer have to be updated
        for ( VViter vv_it = Trimeshes[0].vv_iter( Trimeshes[0].vertex_handle(vs) ); vv_it.is_valid(); ++vv_it )
          for ( VOHiter voh_it = Trimeshes[0].voh_iter( *vv_it ); voh_it.is_valid(); ++voh_it ){
            newEdgesCandidates.push_back(  static_cast<int>( floor( voh_it->idx() / 2 ) ) );
            emap[ static_cast<int>( floor( voh_it->idx() / 2 ) ) ].energy = Inf;
            HalfedgeHandle sucHandle = Trimeshes[0].next_halfedge_handle( *voh_it );
            newEdgesCandidates.push_back(  static_cast<int>( floor( sucHandle.idx() / 2 ) ) );
            emap[ static_cast<int>( floor( sucHandle.idx() / 2 ) ) ].energy = Inf;
          }


        // update face matrices
        for( int j = 0; j < Trimeshes.size(); j++ ){
          // update deleted faces (due to index swapping in OM)
          if( faceIdx[0] < Trimeshes[j].getNumFaces() )
            updateFaceMatrix( Trimeshes[j], faceIdx[0], matCont[j][ faceIdx[0] ] );
          if( (faceIdx[1] != -1) & ( faceIdx[1] < Trimeshes[j].getNumFaces() ) )
            updateFaceMatrix( Trimeshes[j], faceIdx[1], matCont[j][ faceIdx[1] ] );
          // update adjacent faces
          updateFaceMatrices( Trimeshes[j], vs, matCont[j] );
        }

        // remove edges that have to be re-computed (i.e. the edges whose costs were set to Inf previously)
        aol::MultiVector<RealType> mv;
        heap.remove( EdgePointerType( EdgeType( 0, Inf, mv ) ) );

        // built temporary heap which stores all "new" edges
        HeapType tmpHeap;
        newEdgesCandidates.sort();
        std::list<int>::iterator iter = newEdgesCandidates.begin();
        int oldIdx = *iter - 1;
        while( (iter != newEdgesCandidates.end()) & (*iter < currentNumEdges) ){
          // ignore indices, that have been counted twice
          if( *iter > oldIdx){
            RealType e (0.);
            aol::MultiVector<RealType> v_cont( Trimeshes.size(), 3 );
            for( int i = 0; i < Trimeshes.size(); i++ )
              e += calcEnergyAtHalfedge( Trimeshes[i], matCont[i], HalfedgeHandle( 2 * (*iter) ), v_cont[i] );
            emap[ *iter ] = EdgeType( *iter, e, v_cont );
            tmpHeap.push_back( EdgePointerType( emap[ *iter ] ) );
            oldIdx = *iter;
          }
          ++iter;
        }

        // sort temp. heap and merge into the other heap
        tmpHeap.sort();
        heap.merge( tmpHeap );

      }
      // increment progressbar
      pb++;
    }

    // finish progressbar
    pb.finish();

    // build projection map
    for( int i = 0; i < numVertices - counter; i++ )
      for( unsigned int j = 0; j < preImages[i].size(); j++ )
        decInfo.vertexProjection[ preImages[i][j] ] = i;

  }

  //! Prolongation by means of the stored displacements 
  //! Note: decInfo must correspond to the coarse meshes!
  void prolongate ( const TriMeshContainer& coarseMeshes, TriMeshContainer& fineMeshes, const DecimationInfo<RealType>& decInfo  ){

    if( coarseMeshes.size() != fineMeshes.size() )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::prolongate(): sizes of input data do not match!", __FILE__, __LINE__ );
    int N = decInfo.collapses.size();
    int n = coarseMeshes[0].getNumVertices();
    int numVertices = fineMeshes[0].getNumVertices();
    if( n + N != numVertices )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::prolongate(): invalid sizes!", __FILE__, __LINE__ );

    // for all meshes
    for( int j = 0; j < coarseMeshes.size(); j++ ){
      // first: reset the vertices that have not been removed
      for( int i = 0; i < n; i++ ){
        aol::Vec3<RealType> v;
        coarseMeshes[j].getVertex( i, v );
        fineMeshes[j].setVertex( decInfo.vertexMap[ i ], v );
      }

      // second: compute positions of additional vertices by applying the displacements
      for( int i = 0; i < N; i++ ){
        Collapse<RealType> col = decInfo.collapses[N-1-i];
        aol::Vec3<RealType> v;
        // position of vs
        fineMeshes[j].getVertex(  col.vsOrig, v );
        fineMeshes[j].setVertex( col.vsOrig, v - col.delta_vs[j]  );
        // position of vt
        fineMeshes[j].setVertex( decInfo.vertexMap[ n + i ], v - col.delta_vt[j]  );
      }
    }
  }

  //! Prolongation of all meshes by means of the stored displacements
  //! Note: decInfo must correspond to the coarse meshes meshes_c!
  //! The "relative" prolongation is stored in pInfoContainer and can then be applied to other meshes afterwards!
  void calcProlongation ( const TriMeshContainer& coarseMeshes,
                          TriMeshContainer& fineMeshes,
                          const DecimationInfo<RealType>& decInfo,
                          aol::RandomAccessContainer< ProlongInfo<RealType> >& pInfoContainer ){

    if( coarseMeshes.size() != fineMeshes.size() )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::calcProlongation(): sizes of input data do not match(1)!", __FILE__, __LINE__ );
    if( coarseMeshes.size() != pInfoContainer.size() )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::calcProlongation(): sizes of input data do not match(2)!", __FILE__, __LINE__ );

    for( int j = 0; j < coarseMeshes.size(); j++ )
      calcProlongation ( coarseMeshes[j], fineMeshes[j], decInfo,  pInfoContainer[j], j );

  }


  //! Prolongation of a single mesh by means of the stored displacements
  //! The "relative" prolongation is stored in pInfo and can then be applied to other meshes afterwards!
  //!
  //! The integer 'activeMesh' refers to the number of coarseMesh if it has been decimated with other meshes simultaneously.
  //! As decInfo stores information (e.g. of the displacements) which differs from mesh to mesh, 'activeMesh' must coincide with the
  //! number of coarseMesh in the corresponding decimation!
  void calcProlongation (  const TriMeshType& coarseMesh,
                           TriMeshType& fineMesh,
                           const DecimationInfo<RealType>& decInfo,
                           ProlongInfo<RealType>& pInfo,
                           int activeMesh = 0  ){

    int N = decInfo.collapses.size();
    int n = coarseMesh.getNumVertices();
    if( n + N != fineMesh.getNumVertices() )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::calcProlongation(): invalid sizes!", __FILE__, __LINE__ );
    RealType h = coarseMesh.H();

    pInfo.clear();

    // first: reset the vertices that have not been removed
    for( int i = 0; i < n; i++ ){
      aol::Vec3<RealType> v;
      coarseMesh.getVertex( i, v );
      fineMesh.setVertex( decInfo.vertexMap[ i ], v );
      pInfo.vertexMap.push_back( decInfo.vertexMap[ i ] );
    }

    // second: compute positions of additional vertices by applying the displacements
    //         and store the displacement in relative coordinates
    for( int i = 0; i < N; i++ ){
      Collapse<RealType> col = decInfo.collapses[N-1-i];
      aol::Vec3<RealType> v_new;
      pInfo.vs.push_back( col.vsOrig );

      // calculate projection of v_old = vs_old = vs_new - vsplit.delta_vs on the coarse mesh meshes_c
      fineMesh.getVertex(  col.vsOrig, v_new );
      aol::Vec3<RealType> v_old( v_new );
      v_old -= col.delta_vs[ activeMesh ];
      // look for vertex in neighbourhood of vs in the original mesh that is closest to vs_new
      RelativeCoords vsInRelCoord;
      RealType dist = computeProjection( coarseMesh, decInfo.vertexProjection[  col.vsOrig  ], v_old, vsInRelCoord );
      // if the distance is still too big (i.e. the prior was too bad ) look in the 1-ring neighbourhood
      if( dist > h/100. )
        dist = computeProjectionOn1Ring( coarseMesh, decInfo.vertexProjection[  col.vsOrig  ],v_old, vsInRelCoord );
      // if the distance is still too big try each and every face in the mesh ( = O(#faces) )
      if( dist > h/100. )
        dist = computeProjectionGreedy( coarseMesh, v_old, vsInRelCoord );
      // add relative coordinates to pInfo
      pInfo.relCoordsCont.pushBack( vsInRelCoord );

      // compute (relative) position of vs_old
      aol::Vec3<RealType> vProlongated;
      computePosition( coarseMesh, vsInRelCoord, vProlongated );
      fineMesh.setVertex( col.vsOrig, vProlongated );

      // calculate projection of v_old = vt_old = vs_new - vsplit.delta_vt on the coarse mesh meshes_c
      v_old = v_new;
      v_old -= col.delta_vt[ activeMesh ];
      RelativeCoords vtInRelCoord;
      dist = computeProjection( coarseMesh, decInfo.vertexProjection[  col.vsOrig  ], v_old, vtInRelCoord );
      if( dist > h/100. )
        dist = computeProjectionOn1Ring( coarseMesh, decInfo.vertexProjection[  col.vsOrig  ], v_old, vtInRelCoord );
      if( dist > h/100. )
        dist = computeProjectionGreedy( coarseMesh, v_old, vtInRelCoord );
      pInfo.relCoordsCont.pushBack( vtInRelCoord );

      // compute (relative) position of vt_old
      computePosition( coarseMesh, vtInRelCoord, vProlongated );
      fineMesh.setVertex( decInfo.vertexMap[ n + i ], vProlongated );

      // copy information of decInfo to pInfo ( Hence we do not need decInfo any more afterwards!)
      pInfo.vertexMap.push_back( decInfo.vertexMap[ n + i ] );
    }
  }


//! Prolongation of a single mesh only based on the relative coordinates
//! Note: This prolongation can be applied to an arbitrary mesh geometry (However, the topology must be appropriate!)
void prolongate ( const TriMeshType& coarseMesh, TriMeshType& fineMesh, const ProlongInfo<RealType>& pInfo ){

    int N = pInfo.vs.size();
    int n = coarseMesh.getNumVertices();
    int numVertices = fineMesh.getNumVertices();
    if( n + N != numVertices )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::prolongate(): invalid sizes(1)!", __FILE__, __LINE__ );
    if( pInfo.vertexMap.size() != numVertices )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::prolongate(): invalid sizes(2)!", __FILE__, __LINE__ );
    if( pInfo.relCoordsCont.size() != N + N )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::prolongate(): invalid sizes(3)!", __FILE__, __LINE__ );

    // first: reset the vertices that have not been removed
    for( int i = 0; i < n; i++ ){
      aol::Vec3<RealType> v;
      coarseMesh.getVertex( i, v );
      fineMesh.setVertex( pInfo.vertexMap[ i ], v );
    }

    // second: compute positions of additional vertices by the relative coordinates
    for( int i = 0; i < N; i++ ){
      // calc. position of vs_old
      aol::Vec3<RealType> vProlongated;
      RelativeCoords vsInRelCoord( pInfo.relCoordsCont[2*i] );
      computePosition( coarseMesh, vsInRelCoord, vProlongated );
      fineMesh.setVertex( pInfo.vs[i], vProlongated );

      // calc. position of vt_old
      RelativeCoords vtInRelCoord( pInfo.relCoordsCont[2*i + 1] );
      computePosition( coarseMesh, vsInRelCoord, vProlongated );
      fineMesh.setVertex( pInfo.vertexMap[ n + i ], vProlongated );
    }
  }


//! Prolongation of a single mesh only based on the relative coordinates
//! Note: This prolongation can be applied to an arbitrary mesh geometry (However, the topology must be appropriate!)
//! The refined mesh is here given by its corresponding aol::MultiVector<>
void prolongate ( const TriMeshType& coarseMesh, aol::MultiVector<RealType>& fineMeshMV, const ProlongInfo<RealType>& pInfo ){

    int N = pInfo.vs.size();
    int n = coarseMesh.getNumVertices();
    if( (fineMeshMV.numComponents() == 0) && fineMeshMV.isDataOwned() )
      fineMeshMV.reallocate( 3, n+N );
 
    if( pInfo.vertexMap.size() != n + N )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::prolongate(): invalid sizes(1)!", __FILE__, __LINE__ );
    if( pInfo.relCoordsCont.size() != N + N )
      throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::prolongate(): invalid sizes(2)!", __FILE__, __LINE__ );  

    // first: reset the vertices that have not been removed
    for( int i = 0; i < n; i++ ){
      aol::Vec3<RealType> v;
      coarseMesh.getVertex( i, v );
      for( int j = 0; j < 3; j++ )
        fineMeshMV[j][ pInfo.vertexMap[ i ] ] = v[j];
    }

    // second: compute positions of additional vertices by the relative coordinates
    for( int i = 0; i < N; i++ ){
      // calc. position of vs_old
      aol::Vec3<RealType> vProlongated;
      RelativeCoords vsInRelCoord( pInfo.relCoordsCont[2*i] );
      computePosition( coarseMesh, vsInRelCoord, vProlongated );
      for( int j = 0; j < 3; j++ )
        fineMeshMV[j][ pInfo.vs[i] ] = vProlongated[j];

      // calc. position of vt_old
      RelativeCoords vtInRelCoord( pInfo.relCoordsCont[2*i + 1] );
      computePosition( coarseMesh, vtInRelCoord, vProlongated );
      for( int j = 0; j < 3; j++ )
        fineMeshMV[j][ pInfo.vertexMap[ n + i ] ] = vProlongated[j];
    }

  }


protected:

  //! Compute position from relative coordinates in relCoords.
  //! Let T = (t0, t1, t2) be the triangle given by relCoords.face, N its weighted normal.
  //! v = t0 + bar01*(t1 - t0) + bar1*(t2 - t0) + c/|N|^2 * N
  //! Here c corresponds to the additional coordinate given by relCoords.projN
  void computePosition( const TriMeshType& mesh, const RelativeCoords& relCoords, aol::Vec3<RealType>& v ){
      // Note: The numbering of the edges in the face must correspond to the
      //       ordering that was used when the rel. coord. were computed! ( cf. computeProjectionOnTriangle() )
      HalfedgeHandle h =  mesh.halfedge_handle( mesh.face_handle( relCoords.face ) );
      HalfedgeHandle nh = mesh.next_halfedge_handle( h );
      Point p( mesh.point( mesh.from_vertex_handle( h ) ) );
      aol::Vec3<RealType> e0 ( mesh.point( mesh.to_vertex_handle( h ) ) - p );
      aol::Vec3<RealType> e1 ( mesh.point( mesh.to_vertex_handle( nh ) ) - p );
      v = p;
      v.addMultiple( e0, relCoords.bar0 );
      v.addMultiple( e1, relCoords.bar1 );
      aol::Vec3<RealType> wn;
      mesh.getFaceNormal ( mesh.face_handle( relCoords.face ), wn, false );
      v.addMultiple( wn, relCoords.projN / (wn*wn) );
  }


  //! Compute projection of x onto the face T specified by its handle.
  //! If x is closer to T than 'dist', then its relative coordinates are stored and dist is overwritten.
  RealType computeProjection( const TriMeshType& mesh, const FaceHandle& fh, const aol::Vec3<RealType>& x, RelativeCoords& relCoords, RealType& dist ){

    aol::Vec2<RealType> barCoord;
    aol::Vec3<RealType> projection;

    HalfedgeHandle h =  mesh.halfedge_handle( fh );
    HalfedgeHandle nh = mesh.next_halfedge_handle( h );
    RealType tmp = computeProjectionOnTriangle( mesh.point( mesh.from_vertex_handle( h ) ),
                                                mesh.point( mesh.to_vertex_handle( h ) ),
                                                mesh.point( mesh.to_vertex_handle( nh) ), x, barCoord, projection, true );

    if( tmp < dist ){
        dist = tmp;
        relCoords.face = fh.idx();
        relCoords.bar0 = barCoord[0];
        relCoords.bar1 = barCoord[1];
        aol::Vec3<RealType> wn;
        mesh.getFaceNormal ( fh, wn, false );
        aol::Vec3<RealType> diff( x - projection );
        relCoords.projN = diff*wn;
    }

    return dist;

  }

  //! Compute projection of x onto the neighbourhood (i.e. all adjacent faces) of the vertex v
  RealType computeProjection( const TriMeshType& mesh, int v, const aol::Vec3<RealType>& x, RelativeCoords& relCoords ){

    RealType dist = Inf;
    for ( CVFiter vf_it = mesh.cvf_iter( mesh.vertex_handle( v ) ); vf_it.is_valid(); ++vf_it )
      dist = computeProjection( mesh, *vf_it, x, relCoords, dist );
    return dist;

  }

  //! Compute projection of x onto the neighbourhood (i.e. all adjacent faces of the vertex 1-ring) of the vertex v
  RealType computeProjectionOn1Ring( const TriMeshType& mesh, int v, const aol::Vec3<RealType>& x, RelativeCoords& relCoords ){

    std::list<int> faces;
    RealType dist = Inf;

    for ( CVViter vv_it = mesh.cvv_iter( mesh.vertex_handle( v ) ); vv_it.is_valid(); ++vv_it )
      for ( CVFiter vf_it = mesh.cvf_iter( *vv_it ); vf_it.is_valid(); ++vf_it )
        faces.push_back( vf_it->idx() );

    // ignore faces, that have been counted twice
    faces.sort();
    std::list<int>::iterator iter = faces.begin();
    int old = *iter -1;
    for( ; iter != faces.end(); ++iter )
      if( *iter > old ){
        dist = computeProjection( mesh, mesh.face_handle(*iter), x, relCoords, dist );
        old = *iter;
    }

     return dist;
  }

  //! Compute projection of x onto the whole mesh ( linear time )
  RealType computeProjectionGreedy( const TriMeshType& mesh, const aol::Vec3<RealType>& x, RelativeCoords& relCoords ){

     RealType dist = Inf;
     for( int i = 0; i < mesh.getNumFaces(); i++ )
      dist = computeProjection( mesh, mesh.face_handle( i ), x, relCoords, dist );
     return dist;

  }

  //! Does the edge or at least one of its vertices touch the boundary of the mesh?
  bool edgeIsAtBoundary( const TriMeshType& mesh, const HalfedgeHandle& heh ){
    bool aux = mesh.is_boundary( heh ) || mesh.is_boundary( mesh.opposite_halfedge_handle(heh) );
    return aux || mesh.is_boundary( mesh.from_vertex_handle( heh ) ) || mesh.is_boundary( mesh.to_vertex_handle( heh ) );
  }

  //! Does the edge or at least one of its vertices touch the boundary of the mesh?
  bool edgeIsAtBoundary( const TriMeshType& mesh, int edge ){
    return edgeIsAtBoundary( mesh, HalfedgeHandle( 2 * edge ) ) || edgeIsAtBoundary( mesh, HalfedgeHandle( 2 * edge + 1 ) );
  }

  //! Calculate the cost of collapsing an edge for every edge in the mesh
  void fillHeap( const TriMeshContainer& Trimeshes, MatrixContainer& matCont, HeapType& heap , EdgeEnergyMap& emap  ){

    emap.clear();
    heap.clear();

    matCont.reallocate( Trimeshes.size() );
    for( int j = 0; j < Trimeshes.size(); j++ )
     calcErrorQuadricMatrices( Trimeshes[j], matCont[j] );

    for ( CEiter e_it = Trimeshes[0].edges_begin(); e_it != Trimeshes[0].edges_end(); ++e_it){
      RealType e = 0.;
      aol::MultiVector<RealType> v_cont( Trimeshes.size(), 3 );
      for( int i = 0; i < Trimeshes.size(); i++ )
        e += calcEnergyAtHalfedge( Trimeshes[i], matCont[i], HalfedgeHandle( 2 * e_it->idx() ), v_cont[i] );
      emap[ e_it->idx() ] = EdgeType( e_it->idx(), e, v_cont );
      heap.push_back( EdgePointerType( emap[ e_it->idx() ] ) );
    }

  }


  //! The all important energy functional defined by Garland and Heckbert, "Surface Simplification Using Quadric Error Metrics"
  //! Notation corresponds to that used in the article.
  RealType calcEnergyAtHalfedge( const TriMeshType& mesh, const aol::MultiVector<RealType>& Mats, const HalfedgeHandle& heh, aol::Vector<RealType>& v_bar ){
    
   if( v_bar.size() != 3 )
     throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::calcEnergyAtHalfedge(): invalid size!", __FILE__, __LINE__ );

   //!TODO Edges at boundary (ie having one or two boundary nodes) should be treated more carefully!
     
   // edge is not boundary but both nodes are boundary nodes => edge should not be removed!
   HalfedgeHandle oh = mesh.opposite_halfedge_handle(heh);
   v_bar.setZero();
   if( mesh.is_boundary(mesh.from_vertex_handle( heh ))  && mesh.is_boundary(mesh.to_vertex_handle( heh ) ) )
     if( !edgeIsAtBoundary( mesh, heh ) && !edgeIsAtBoundary( mesh, oh ) )     
       return Inf;

   // Q  = Q_1 + Q_2, where Q_i is the matrix corresponding to v_i, e = (v_1, v_2)
   // Q_i is the sum of all face matrices, whose face touches v_i
   // Note: Q is a symmetric 4x4-matrix and can therefore be stored in a vector of 10 entries
   aol::Vector<RealType> Q (10);
   Q.setZero();   
   for ( CVFiter vf_it = mesh.cvf_iter( mesh.from_vertex_handle( heh ) ); vf_it.is_valid(); ++vf_it )
     Q += Mats[ vf_it->idx() ];
   for ( CVFiter vf_it = mesh.cvf_iter( mesh.to_vertex_handle( heh ) ); vf_it.is_valid(); ++vf_it ){
     // Don't count overlapping faces twice!
     if( ( *vf_it != mesh.face_handle(heh) ) ){
       if( !mesh.is_boundary( oh ) & ( *vf_it == mesh.face_handle( oh ) ) )
         continue;
       else
         Q += Mats[ vf_it->idx() ];
     }
   }

    // Optimal position v_bar is the last row of A^{-1}:
    //         | Q0 Q1 Q2 Q3 | 0 |      | Q0 Q1 Q2 0 | -Q3 |
    // [A|b] = | Q1 Q4 Q5 Q6 | 0 | <->  | Q1 Q4 Q5 0 | -Q6 |
    //         | Q2 Q5 Q7 Q8 | 0 |      | Q2 Q5 Q7 0 | -Q8 |
    //         | 0  0  0  1  | 1 |      | 0  0  0  1 |  1  |

    // If A is not invertible or v_i is boundary vertex, we set v_bar = 0.5 * (v_1 + v_2)
    RealType detA = Q[0]*Q[4]*Q[7] + 2.*Q[1]*Q[2]*Q[5] - Q[2]*Q[2]*Q[4] - Q[5]*Q[5]*Q[0] - Q[1]*Q[1]*Q[7];
    if( (std::abs( detA ) > std::numeric_limits<RealType>::epsilon()) &&  !mesh.is_boundary(mesh.from_vertex_handle( heh ))  && !mesh.is_boundary(mesh.to_vertex_handle( heh )) ){
      aol::Matrix33Symm< RealType > A( Q[0], Q[1], Q[2], Q[4], Q[5], Q[7] );
      aol::Vec3< RealType > b( -Q[3], -Q[6], -Q[8] );
      aol::Vec3< RealType > old( v_bar );
      A.applyInverseTo( b, v_bar );
    }
    else{
      Point p = mesh.point( mesh.from_vertex_handle( heh ) );
      Point q = mesh.point( mesh.to_vertex_handle( heh ) );
      RealType lambda = mesh.is_boundary( mesh.from_vertex_handle( heh ) ) ? 1.0 : ( mesh.is_boundary( mesh.to_vertex_handle( heh ) ) ? 0.0 : 0.5 );
      for( int i = 0; i < 3; i++ )
        v_bar[i] = lambda * p[i] + (1. - lambda) * q[i];
    }

    // v = (v_bar1, v_bar2, v_bar3, 1), where v_bar would be the optimal position for the vertex e is going to be collapsed to
    aol::Vector<RealType> v (4), Qv (4);
    for( int i = 0; i < 3; i++ )
      v[i] = v_bar[i];
    v[3] = 1.;

    Qv[0] = Q[0]*v[0] + Q[1]*v[1] + Q[2]*v[2] + Q[3];
    Qv[1] = Q[1]*v[0] + Q[4]*v[1] + Q[5]*v[2] + Q[6];
    Qv[2] = Q[2]*v[0] + Q[5]*v[1] + Q[7]*v[2] + Q[8];
    Qv[3] = Q[3]*v[0] + Q[6]*v[1] + Q[8]*v[2] + Q[9];

    return v*Qv;
  }

  //! update face matrices of the faces that contain vertex vs (cf. Garland and Heckbert, "Surface Simplification Using Quadric Error Metrics" )
  void updateFaceMatrices( const TriMeshType& mesh, int vs, aol::MultiVector<RealType>& faceMatrices ){
    for ( CVFiter vf_it = mesh.cvf_iter( mesh.vertex_handle( vs ) ); vf_it.is_valid(); ++vf_it )
      calcErrorQuadricMatrix( mesh, faceMatrices[ vf_it->idx() ], *vf_it );
  }

  //! update one face matrix K_p (cf. Garland and Heckbert, "Surface Simplification Using Quadric Error Metrics" )
  void updateFaceMatrix( const TriMeshType& mesh, int idx, aol::Vector<RealType>& K ){
    calcErrorQuadricMatrix( mesh, K, mesh.face_handle( idx ) );
  }

  //! compute K_p for all faces (cf. Garland and Heckbert, "Surface Simplification Using Quadric Error Metrics" )
  void calcErrorQuadricMatrices( const TriMeshType& mesh, aol::MultiVector<RealType>& faceMatrices ) const {
    faceMatrices.reallocate( mesh.getNumFaces(), 10 );
    for ( typename TriMeshType::ElementIteratorType fIter = mesh; fIter.notAtEnd(); ++fIter )
      calcErrorQuadricMatrix( mesh, faceMatrices[ fIter.getIndex() ], fIter.faceHandle() );
  }

  //! compute one single K_p (cf. Garland and Heckbert, "Surface Simplification Using Quadric Error Metrics" )
  void calcErrorQuadricMatrix( const TriMeshType& mesh, aol::Vector<RealType>& K, const typename TriMeshType::FaceHandle& fhandle ) const {
      if( K.size() != 10 )
	throw aol::Exception ( "om::SimulProgMesh<TriMeshType,RealType>::calcErrorQuadricMatrix(): invalid size!", __FILE__, __LINE__ );
      aol::Vec3<RealType> n;
      mesh.getFaceNormal ( fhandle, n );
      Point point = mesh.point( mesh.to_vertex_handle( mesh.halfedge_handle( fhandle ) ) );
      RealType d = - ( n[0] * point[0] + n[1] * point[1] + n[2] * point[2] );
      K[0] = n[0] * n[0]; K[1] = n[0] * n[1]; K[2] = n[0] * n[2]; K[3] = n[0] * d;
      K[4] = n[1] * n[1]; K[5] = n[1] * n[2]; K[6] = n[1] * d;
      K[7] = n[2] * n[2]; K[8] = n[2] * d;
      K[9] = d * d;
   }
};


} // namespace om

#endif // __PROGRESSIVEMESH_H
