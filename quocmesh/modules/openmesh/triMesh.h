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

#ifndef __TRIMESH_H
#define __TRIMESH_H

#include <openMeshIncludes.h>
#include <aol.h>
#include <smallVec.h>
#include <smallMat.h>
#include <triangMesh.h>
#include <meshWithData.h>

#include <progressiveMesh.h>

namespace om {


//!TODO add comment
struct PMTraits : public OpenMesh::DefaultTraits
{
  VertexAttributes( OpenMesh::Attributes::Status );
  EdgeAttributes( OpenMesh::Attributes::Status );
  FaceAttributes( OpenMesh::Attributes::Status );
};



//=======================================================================================================
//
//       TriMesh
//
//-------------------------------------------------------------------------------------------------------
  //! Note 1: TriMesh does only provide a normal field on faces. To further get a notion of a normal related to vertices or edges,
  //!         please use the derived classes TriMeshWithVertexNormals or TriMeshWithEdgeNormals, respectively.
  //!         ( Especially if one needs to calculate with geometric quantities like curvatures. )
  //! Note 2: There is no notion of any arithmetric operations on this class. If one needs to perform such kind of operations,
  //!         please use the derived class CompTriMesh (which inherits from TriMeshWith{Vertex,Edge}Normals by a template argument).
template <typename _RealType>
class TriMesh : protected OpenMesh::TriMesh_ArrayKernelT<PMTraits> {

public:

  static const int IndexNotSet = -1;
  
  typedef OpenMesh::TriMesh_ArrayKernelT<PMTraits> OpenMeshType;
  typedef _RealType RealType;
  typedef OpenMesh::FaceHandle FaceHandle;
  typedef OpenMesh::VertexHandle VertexHandle;
  typedef OpenMeshType::Point Point;

  friend class SimulProgMesh< TriMesh<RealType>, RealType >;

  //----------- Geometric objects -----------------------------------------------------------------------

  // --- Edge ---
  struct Edge {

    typedef typename TriMesh<RealType>::EdgeHandle HandleType;

    Edge( ) : _handle ( 0 ), _length ( 0 ) { }

    Edge ( HandleType handle, RealType length )
        : _handle ( handle ), _length ( length ) {}


    HandleType _handle;
    RealType _length;

    HandleType get_handle() const { return _handle; }

    bool operator< ( const Edge &e ) const {
      return _length < e._length;
    }

    bool operator> ( const Edge &e ) const {
      return _length > e._length;
    }
  };

  //  --- Element ---
  //! Element is the type of element used in a TriMesh.
  //! Unlike in (regular) planar meshes ( where the single elements do not have many properties )
  //! the elements of hyperplanes should be able to provide more information, eg metric or normal.
  //!
  //! By dereferencing the corresponding element iterator in TriMesh, the metric of the iterated triangle
  //! is computed and the Element object is filled. Thus - when using an element iterator - it is
  //! reasonable to compute quantities like the metric or the gradient of the parametrization Dx
  //! on the element instead of using the corresponding routines on the mesh.
  //!
  //! Furthermore, Element contains a notion of neighbours.


  class Element : public aol::Triangle<RealType>{

    typedef TriMesh<RealType> MeshType;
    typedef typename TriMesh<RealType>::FaceHandle FaceHandle;

  protected:
    MeshType *_mesh;
    // spanning vectors of element
    aol::Vec3<RealType> v[2];
    // global indices of element and nodes
    int _globIdx;
    aol::Vec<3, int> _globNodeIdx;
    // first fund.form
    aol::Matrix22<RealType> g;
    // area
    RealType _area;

    // global indices of neighbours
    mutable FaceHandle globIdxOfNeighbours[3];
    mutable bool hasNeighbours;

    void setMesh ( const MeshType& newMesh ) {
      if ( _mesh )
        delete _mesh;
      _mesh = new MeshType( newMesh );
    }

  public:
    // ElementIterator is friend to set mesh pointer and fill the mesh
    friend class TriMesh::ElementIter;

  public:
    Element() : aol::Triangle<RealType>(), _mesh(0), hasNeighbours( false ){}

    Element( const MeshType& Mesh, int globalIdx ) : aol::Triangle<RealType>(), _globIdx(globalIdx), hasNeighbours( false ){
      _mesh = new MeshType( Mesh );
      fillElement( _mesh->face_handle(_globIdx) );
    }

    ~Element(){
      if ( _mesh )
       delete _mesh;
    }

    int globNodeIdx( int localIndex ) const {
      return _globNodeIdx[ localIndex ];
    }

    int globIdx(  ) const {
      return _globIdx;
    }

    void getDx( aol::Mat<3,2,RealType> &Dx ) const {
      // Dx = [ x1-x0 | x2-x0 ] = [ v0 | v1 ]
      for ( int i = 0; i < 2; i++ )
        Dx.setCol ( i, v[i] );
    }

    void getDxT( aol::Mat<2,3,RealType> &DxT ) const {
      // Dx = [ x1-x0 | x2-x0 ] = [ v0 | v1 ]
      for ( int i = 0; i < 2; i++ )
        DxT.setRow ( i, v[i] );
    }

    // compute first fundamental form
    void computeMetric( ) {
      v[0] = this->edge( 0, 1 );
      v[1] = this->edge( 0, 2 );
      // g = [v0|v1]^T[v0|v1]
      for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < 2; j++ )
          g[i][j] = v[i] * v[j];
      _area = .5 * sqrt( g.det() );
    }

    const aol::Matrix22<RealType>& metric( ) const {
      return g;
    }

    RealType detg( ) const {
      return g.det();
    }

    aol::Matrix22<RealType> ginv() const {
      return g.inverse();
    }

    RealType area( ) const {
      return aol::Max ( _area, std::numeric_limits<RealType>::epsilon () );
    }

    // returns FaceHandle of *this
    FaceHandle faceHandle() const {
      if ( !_mesh )
        throw aol::Exception ( "Trimesh::Element: mesh has not been set yet.", __FILE__, __LINE__ );
      return _mesh->faceHandle( this->_globIdx );
    }

    // fills vector of neighbours
    void computeNeighbours() const {
      if ( !_mesh )
        throw aol::Exception ( "Trimesh::Element: mesh has not been set yet.", __FILE__, __LINE__ );

      if( !hasNeighbours ){
        int localIdx = 0;
        // Openmesh FaceFaceIter does not detect boundary edges! Hence FaceEdgeIter is used here!
        // Note: If *this contains a boundary edge, FaceFaceIter is incremented automatically!
        ConstFaceEdgeIter feIter = _mesh->getOpenMeshObject().cfe_iter( this->faceHandle() );
        ConstFaceFaceIter ffIter = _mesh->getOpenMeshObject().cff_iter( this->faceHandle() );
        for(; feIter.is_valid(); ++feIter){
            // boundary edges get an invalid handle
            if( _mesh->getOpenMeshObject().is_boundary( *feIter ) ){
              globIdxOfNeighbours[localIdx++] = TriMesh<RealType>::InvalidFaceHandle;
            }
            else{
              globIdxOfNeighbours[localIdx++] = *ffIter;
              ++ffIter;
            }
        }
        hasNeighbours = true;
      }
    }

    // returns local index of the edge in the adjacent face (specified by its local index) that also belongs to *this
    int getLocalIndexOfThisFromAdjacentFace( int locIndexOfNeighbour ) const {
      if ( !_mesh )
        throw aol::Exception ( "Trimesh::Element: mesh has not been set yet.", __FILE__, __LINE__ );
      if( !hasNeighbours )
        computeNeighbours();

      // boundary?
      if( globIdxOfNeighbours[ locIndexOfNeighbour ] == TriMesh<RealType>::InvalidFaceHandle )
        return -1;

      // find this face from the neighbour's point of view
      int index = 0;
      // find edge handle of edge with local index locIndexOfNeighbour
      ConstFaceEdgeIter feIter = _mesh->getOpenMeshObject().cfe_iter( this->faceHandle() );
      for( int i = 0; i < locIndexOfNeighbour; i++ )
        ++feIter;
      // find the same edge from the neighbour`s point of view
      ConstFaceEdgeIter feIter_nb = _mesh->getOpenMeshObject().cfe_iter( globIdxOfNeighbours[ locIndexOfNeighbour ] );
      while( ( *feIter != *feIter_nb ) && (index < 3) ){
        ++feIter_nb;
        ++index;
      }

      return ( index < 3 )? index : -1;
    }

    // returns global index of the face adjacent to edge with local index locIndex
    FaceHandle getGlobalIndexOfNeighbour( int locIndex ) const {
      if ( !_mesh )
        throw aol::Exception ( "Trimesh::Element: mesh has not been set yet.", __FILE__, __LINE__ );
      if( !hasNeighbours )
        computeNeighbours();

      return  globIdxOfNeighbours[ locIndex ];
    }


    // For compatibility with qc::Element
    short level() const {
      throw aol::Exception ( "om::Element does not have a level.", __FILE__, __LINE__ );
      return -1;
    }
    unsigned char type() const {
      throw aol::Exception ( "om::Element does not have a type.", __FILE__, __LINE__ );
      return 23;
    }

  protected:
    void fillElement( const FaceHandle& fh ) {
      if ( !_mesh )
        throw aol::Exception ( "Trimesh::Element: mesh has not been set yet.", __FILE__, __LINE__ );

      typename MeshType::ConstFaceVertexIter v_it;
      int i = 0;
      for ( v_it = _mesh->cfv_iter ( fh );  v_it.is_valid(); ++v_it ) {
        typename OpenMeshType::Point p = _mesh->point ( *v_it );
        for( int j = 0; j < 3; ++j )
          this->getNode(i)[j] = p[j];

        _globNodeIdx[i++] = v_it->idx();
      }
      _globIdx = fh.idx();
      computeMetric();
    }

  };

  //----------- Iterators -------------------------------------------------------------------------------

  // --- ElementIterator ---
  //
  //! Iterates over all faces/triangles. If you only use iter.getIndex() or iter.faceHandle(),
  //! this iterator is as fast as the underlying OpenMesh::ConstFaceIter.
  //! By dereferencing, the first fundamental form of the iterated triangle
  //! will be computed and the Element object will be filled (if it hasnt been filled yet!)
  //!
  //! Note that all QuocMesh iterators are const iterators
  //! (unless explicitely called "nonConst").

class ElementIter {
    typedef TriMesh<RealType> MeshT;
  protected:
    mutable Element _cur;
    const MeshT &_mesh;
    typename OpenMeshType::ConstFaceIter _iter;
    mutable bool _filled;

  public:
    typedef ElementIter Self;
    typedef MeshT        BeginType;
    typedef MeshT        EndType;
    typedef Element     IteratedType;

    // constructor from mesh and OpenMesh face iterator
    ElementIter ( const MeshT &mesh, typename OpenMeshType::FaceIter iter )
        : _mesh ( mesh ),
          _iter ( iter ),
          _filled ( false )
    {_cur.setMesh( _mesh );}

    // constructor from mesh
    ElementIter ( const MeshT &mesh )
        : _mesh ( mesh ),
          _iter ( mesh.faces_begin() ),
          _filled ( false )
    {_cur.setMesh( _mesh );}

    bool operator== ( const EndType & ) const {
      return atEnd();
    }

    bool operator!= ( const EndType & ) const {
      return notAtEnd();
    }

    bool atEnd () const {
      typename MeshT::ConstFaceIter end_it = _mesh.faces_end();
      return _iter == end_it;
    }

    bool notAtEnd () const {
      return ! atEnd ();
    }

    Self& operator++( ) {
      ++_iter;
      _filled = false;
      _cur.hasNeighbours = false;
      return *this;
    }

    const IteratedType& operator*() const {
      if (!_filled)
        _cur.fillElement( *_iter );
      _filled = true;
      return _cur;
    }

    const IteratedType* operator->() const {
      if (!_filled)
        _cur.fillElement( *_iter );
      _filled = true;
      return &_cur;
    }

    typename MeshT::FaceHandle faceHandle( ) const {
      return *_iter;
    }

    int getIndex () const {
      return faceHandle().idx();
    }

    aol::Vec3<int> getNodeIndices() const {
      aol::Vec3<int> ret;
      typename MeshT::ConstFaceVertexIter v_it;
      int i = 0;
      for ( v_it = _mesh.cfv_iter ( *_iter ); v_it.is_valid(); ++v_it )
        ret[i++] = v_it->idx();
      return ret;
    }

  };

  // --- NodeIterator ---
  //
  //! Iterates over all nodes in the mesh. Virtually no overhead in regard
  //! to the underlying OpenMesh::ConstVertexIter.
  //!
  //! Note that all QuocMesh iterators are const iterators
  //! (unless explicitely called "nonConst").
  class FullNodeIter {
    typedef TriMesh<RealType>   MeshT;
  protected:
    const MeshT & _mesh;
    typename OpenMeshType::ConstVertexIter _iter;

  public:
    typedef FullNodeIter        Self;
    typedef MeshT               BeginType;
    typedef MeshT               EndType;
    typedef int                 IteratedType;

    // constructor from mesh and OpenMesh face iterator
    FullNodeIter ( const MeshT &mesh, typename OpenMeshType::VertexIter iter )
        : _mesh ( mesh ),
          _iter ( iter )
    {}

    // constructor from mesh
    FullNodeIter ( const MeshT &mesh )
        : _mesh ( mesh ),
          _iter ( mesh.vertices_begin() )
    {}

    bool operator== ( const EndType & ) const {
      return atEnd();
    }

    bool operator!= ( const EndType & ) const {
      return notAtEnd();
    }

    bool atEnd () const {
      return _iter == _mesh.vertices_end();
    }

    bool notAtEnd () const {
      return ! atEnd ();
    }

    Self& operator++( ) {
      ++_iter;
      return *this;
    }

    IteratedType operator*() const {
      return getIndex();
    }

    typename MeshT::VertexHandle vertexHandle( ) const {
      return *_iter;
    }

    int getIndex () const {
      return vertexHandle().idx();
    }

    bool isAtBoundary() const {
      return _mesh.is_boundary( *_iter );
    }

    aol::Vec3<RealType> getCoords () const {
      aol::Vec3<RealType> coords;
      for ( int i = 0; i < 3; ++i )
        coords[i] = const_cast<MeshT &> ( _mesh ).point ( *_iter ) [i];
      return coords;
    }
    
    int getValence() const {
      int ct = 0;
      ConstVertexFaceIter vfIter = _mesh.getOpenMeshObject().cvf_iter (  *_iter );
      for (; vfIter.is_valid(); ++vfIter)	ct++;
      return ct;	
    }
    
    void getAdjacentFaces ( aol::RandomAccessContainer<int>& adjFaces ) const {  
      adjFaces.clear();
      typename OpenMeshType::ConstVertexFaceIter vfIter = _mesh.getOpenMeshObject().cvf_iter (  *_iter );
      for (; vfIter.is_valid(); ++vfIter ) 
	adjFaces.pushBack( vfIter->idx() );
    }

    void getAdjacentNodes ( aol::RandomAccessContainer<int>& adjNodes ) const {  
      adjNodes.clear();
      typename OpenMeshType::ConstVertexVertexIter vvIter = _mesh.getOpenMeshObject().cvv_iter (  *_iter );
      for (; vvIter.is_valid(); ++vvIter ) 
	adjNodes.pushBack( vvIter->idx() );
    }
  };

  // --- BoundaryNodeIterator ---
  //
  //! Cave: This iterator does not iterate over the whole (possibly non-connected)
  //! boundary, but only over the part that belongs to the half-edge given via its
  //! constructor. If no half-edge is given, it will iterate over the boundary component
  //! that belongs to the first boundary edge.
  //!
  //! Note that all QuocMesh iterators are const iterators
  //! (unless explicitely called "nonConst").
  class BoundaryNodeIter {
    typedef TriMesh<RealType>       MeshT;
  protected:
    const MeshT & _mesh;
    typename OpenMeshType::HalfedgeHandle _iter;
    typename OpenMeshType::HalfedgeHandle _begin;
    bool _atStart;

  public:
    typedef BoundaryNodeIter        Self;
    typedef MeshT                   BeginType;
    typedef MeshT                   EndType;
    typedef int                     IteratedType;

    //! constructor from mesh
    BoundaryNodeIter ( const MeshT & mesh )
        : _mesh ( mesh ),
          _iter ( getFirstBoundaryHalfedge ( mesh ) ),
          _begin ( _iter ),
          _atStart ( true )
    {}

    BoundaryNodeIter ( const MeshT & mesh, const typename MeshT::HalfedgeHandle & halfedge )
        : _mesh ( mesh ),
          _iter ( halfedge ),
          _begin ( _iter ),
          _atStart ( true )
    {}

    bool operator== ( const EndType & ) const {
      return atEnd();
    }

    bool operator!= ( const EndType & ) const {
      return notAtEnd();
    }

    bool atEnd () const {
      return ! ( _iter.is_valid() ) || ( !_atStart && _iter == _begin );
    }

    bool notAtEnd () const {
      return ! atEnd ();
    }

    Self& operator++( ) {
      _iter = _mesh.next_halfedge_handle ( _iter );
      _atStart = false;
      return *this;
    }

    IteratedType operator*() const {
      return getIndex();
    }

    typename MeshT::VertexHandle vertexHandle( ) const {
      return const_cast<MeshT &> ( _mesh ).to_vertex_handle ( _iter );
    }

    int getIndex () const {
      return vertexHandle().idx ();
    }

    aol::Vec3<RealType> getCoords () const {
      aol::Vec3<RealType> coords;
      typename MeshT::VertexHandle vHandle = _mesh.to_vertex_handle ( _iter );
      for ( int i = 0; i < 3; ++i )
        coords[i] = const_cast<MeshT &> ( _mesh ).point ( vHandle ) [i];
      return coords;
    }
   

  protected:
    static typename MeshT::HalfedgeHandle getFirstBoundaryHalfedge (const MeshT & mesh ) {
      ConstVertexIter vIter = mesh.vertices_begin();
      for (; !mesh.is_boundary ( *vIter ); ++vIter) ;
      ConstVertexIHalfedgeIter hIter ( mesh, *vIter );
      for (; !mesh.is_boundary ( *hIter ); ++hIter ) ;
      return *hIter;
    }
  };
  
  // --- EdgeIterator ---
  //
  //! Iterates over all edges in the mesh. Virtually no overhead in regard
  //! to the underlying OpenMesh::ConstEdgeIter.
  //!
  //! Note that all QuocMesh iterators are const iterators
  //! (unless explicitely called "nonConst").
  class EdgeIter {
    typedef TriMesh<RealType>   MeshT;
  protected:
    const MeshT & _mesh;
    typename OpenMeshType::ConstEdgeIter _iter;
    typedef typename OpenMeshType::HalfedgeHandle  HalfEdgeHandle;
    
  public:
    typedef EdgeIter        Self;
    typedef MeshT           BeginType;
    typedef MeshT           EndType;
    typedef int             IteratedType;

    // constructor from mesh and OpenMesh face iterator
    EdgeIter ( const MeshT &mesh, typename OpenMeshType::EdgeIter iter )
        : _mesh ( mesh ),
          _iter ( iter )
    {}

    // constructor from mesh
    EdgeIter ( const MeshT &mesh )
        : _mesh ( mesh ),
          _iter ( mesh.edges_begin() )
    {}

    bool operator== ( const EndType & ) const {
      return atEnd();
    }

    bool operator!= ( const EndType & ) const {
      return notAtEnd();
    }

    bool atEnd () const {
      return _iter == _mesh.edges_end();
    }

    bool notAtEnd () const {
      return ! atEnd ();
    }

    Self& operator++( ) {
      ++_iter;
      return *this;
    }

    IteratedType operator*() const {
      return getIndex();
    }

    typename MeshT::EdgeHandle edgeHandle( ) const {
      return *_iter;
    }

    int getIndex () const {
      return edgeHandle().idx();
    }
    
    //!TODO check, how dihedral angle is actually computed in OpenMesh!
    RealType getDihedralAngle () const {
      return _mesh.calc_dihedral_angle( *_iter );
    }
    
    RealType getLength () const {
      return _mesh.calc_edge_length( *_iter );
    }
    
    void getAdjacentFaces ( int& i, int& j ) const {
      // index equals -1 if *iter is boundary edge
      i = _mesh.face_handle( hEdgeHandle(0) ).idx();
      j = _mesh.face_handle( hEdgeHandle(1) ).idx();
    }
    
    void getAdjacentNodes ( int& i, int& j ) const {
      i = _mesh.from_vertex_handle( hEdgeHandle(0) ).idx();
      j = _mesh.to_vertex_handle( hEdgeHandle(0) ).idx();
    }
    
    void getOppositeNodes ( int& i, int& j ) const {
      // index equals -1 if *iter is boundary edge
      i = _mesh.opposite_vh( hEdgeHandle(0) ).idx();
      j = _mesh.opposite_vh( hEdgeHandle(1) ).idx();
    }
  
    bool isBoundary(  ) const {
      return _mesh.is_boundary( *_iter );
    }
 
    protected:
    inline HalfEdgeHandle hEdgeHandle( int i ) const {
      return _mesh.halfedge_handle( *_iter, i);
    }


  };
  //---- end of iterators -------------------------------------------------------------------------------

  typedef Element            ElementType;
  typedef ElementIter        ElementIteratorType;
  typedef FullNodeIter       NodeIteratorType;
  typedef BoundaryNodeIter   BoundaryNodeIteratorType;  
  typedef EdgeIter           EdgeIteratorType;

  //! Remark on iterators and local coordinates in Openmesh and Trimesh:
  //! Fixing one particular triangle, the local numbering of the vertices, edges and adjacent faces
  //! corresponds to the running order of the Openmesh face-iterators Face{Vertex,Edge,Face}Iterator, respectively.
  //! Vertex i (in local coordinates) has edge i (in local coordinates) as incoming directed edge, and edge (i+2) mod 3 as outgoing directed edge.
  //! Edge i (in local coordinates) has face i (in local coordinates) as adjacent face.


  // constructor
  TriMesh( ) { }

  TriMesh( const std::string filename ) {
    loadFromPLY ( filename );
  }

  // get (protected, thus normally invisible) parent class
  OpenMeshType & getOpenMeshObject () {
    return *this;
  }

  // get (protected, thus normally invisible) parent class
  const OpenMeshType & getOpenMeshObject () const {
    return *this;
  }

  int getNumberOfNodes () const {
    return static_cast<int> ( n_vertices() );
  }

  int getNumVertices () const {
    return static_cast<int> ( n_vertices() );
  }

  int getNumFaces () const {
    return static_cast<int> ( n_faces() );
  }
  
  int getNumberOfEdges () const {
    return static_cast<int> ( n_edges() );
  }
  
  int getNumEdges () const {
    return static_cast<int> ( n_edges() );
  }

  VertexHandle vertexHandle( int idx ) const {
    if ( idx < 0 || idx >= getNumVertices()  ) {
      throw aol::Exception ( "wrong index", __FILE__, __LINE__ );
    }
    return getOpenMeshObject().vertex_handle ( idx );
  }

  FaceHandle faceHandle( int idx ) const {
    if ( idx < 0 || idx >= getNumFaces() ) {
      throw aol::Exception ( "wrong index", __FILE__, __LINE__ );
    }
    return getOpenMeshObject().face_handle ( idx );
  }

  // only considers one connected component of boundary!
  // cf. documentation of BoundaryNodeIterator.
  void fillBoundaryMask( aol::BitVector& mask ) const {
    if( mask.size() != static_cast<int> ( n_vertices() ) )
      throw aol::Exception ( "om::TriMesh::fillBoundaryMask(): mask has wrong size.", __FILE__, __LINE__ );
    mask.setZero();
    for ( BoundaryNodeIteratorType bvIter = *this; bvIter.notAtEnd(); ++bvIter )
      mask.set( bvIter.getIndex() , true );
  }
  
  void fillFullBoundaryMask( aol::BitVector& mask ) const {
    if( mask.size() != static_cast<int> ( n_vertices() ) )
      throw aol::Exception ( "om::TriMesh::fillFullBoundaryMask(): mask has wrong size.", __FILE__, __LINE__ );
    mask.setZero();
    for ( NodeIteratorType iter = *this; iter.notAtEnd(); ++iter )
      if( iter.isAtBoundary() )
        mask.set( iter.getIndex() , true );
  }

  //---- mesh modification ------------------------------------------------------------------------------
  int addVertex ( const aol::Vec3<RealType> & coords ) {
    VertexHandle vHandle = getOpenMeshObject().add_vertex ( Point( coords[0], coords[1], coords[2] ) );
    return vHandle.idx();
  }

  int setVertex( const VertexHandle &vh, const aol::Vec3<RealType> & v ){
      Point p;
      for(int i = 0; i < 3; ++i) p[i] = v[i];
      this->set_point( vh, p );
      return vh.idx();
  }

  int setVertex( const int &idx, const aol::Vec3<RealType> & v ){
     return setVertex( vertexHandle( idx ), v );
  }

  int getVertex( const VertexHandle &vh, aol::Vec3<RealType> & v ) const{
      Point p = getOpenMeshObject().point ( vh );
      for(int i = 0; i < 3; ++i) v[i] = p[i];
      return vh.idx();
  }

  int getVertex( const int &idx, aol::Vec3<RealType> & v ) const{
      return getVertex( vertexHandle( idx ), v );
  }

  int addFace ( const aol::Vec3<int> & vertexIdx ) {
    FaceHandle fHandle = getOpenMeshObject().add_face ( vertexHandle(vertexIdx[0]), vertexHandle(vertexIdx[1]), vertexHandle(vertexIdx[2]) );
    if( fHandle == InvalidFaceHandle )
       throw aol::Exception ( "probably orientation problem (OpenMesh can only deal with orientable surfaces!)", __FILE__, __LINE__ );
    return fHandle.idx();
  }

  // mind the comment on local numbering above!!!
  void getEdge( const FaceHandle &fh, int localCoord, aol::Vec3<RealType> & e ) const{
      if( !( localCoord < 3 ) )
       throw aol::Exception ( "wrong local coordinate!", __FILE__, __LINE__ );
      Normal omEdge;
      ConstFaceEdgeIter fe_it = getOpenMeshObject().cfe_iter( fh );
      for( int i = 0; i < localCoord; i++ )
        ++fe_it;
      getOpenMeshObject().calc_edge_vector( *fe_it, omEdge );
      for( int i = 0; i < 3; i++ )
        e[i] = omEdge[i];
  }
  
  void setVertexToAverageOfNeighbours( const VertexHandle &vh );
  
  void setVertexToAverageOfNeighbours( const int &idx ){
    return setVertexToAverageOfNeighbours( vertexHandle( idx ) );
  }

  //---- geometric properties ---------------------------------------------------------------------------

  RealType enclosedVolume( ) const;
  RealType area( ) const;

  RealType vol( FaceHandle fh ) const {
    aol::Vec3<RealType> n;
    getFaceNormal(fh, n, false);
    return n.norm()/2.0;
  }
  RealType vol( int i ) const {
    return vol ( faceHandle( i ) );
  }

  RealType getAspectRatio( const FaceHandle &fh ) const;

  double H () const {
    double h = 0.;
    // iterate over all edges
    for (OpenMeshType::ConstEdgeIter e_it=this->edges_begin(); e_it!=this->edges_end(); ++e_it)
      if ( OpenMeshType::calc_edge_length( *e_it ) > h)
        h = OpenMeshType::calc_edge_length( *e_it );
    return h;
  }
  
  aol::Vec3<int> getNodeIndices( int faceIdx ) const {
    aol::Vec3<int> idx( IndexNotSet );
    int i = 0;
    typename OpenMeshType::ConstFaceVertexIter fv_it = this->getOpenMeshObject().cfv_iter( this->getOpenMeshObject().face_handle(faceIdx) );
    for(; fv_it.is_valid(); ++fv_it)
      idx[i++] = fv_it->idx();
    return idx;
  }
  
  // ordering convention: i'th neighboring triangle is opposite of node with local index i
  aol::Vec3<int> getTriangleIndicesOfNeighbours( int faceIdx ) const {
    aol::Vec3<int> idx( IndexNotSet );
    typename OpenMeshType::HalfedgeHandle heh = this->getOpenMeshObject().halfedge_handle( this->getOpenMeshObject().face_handle(faceIdx) );
    for( int i = 0; i < 3; i++ ) {      
      idx[(i+1)%3] = this->getOpenMeshObject().opposite_face_handle( heh ).idx();
      heh = this->getOpenMeshObject().next_halfedge_handle(heh);
    }
    return idx;
  }
  
  // ordering convention: i'th entry is node in i'th neighboring triangle that does not belong to faceIdx
  aol::Vec3<int> getOppositeNodesOfNeighboringTriangles( int faceIdx ) const {
    aol::Vec3<int> idx( IndexNotSet );
    typename OpenMeshType::HalfedgeHandle heh = this->getOpenMeshObject().halfedge_handle( this->getOpenMeshObject().face_handle(faceIdx) );
    for( int i = 0; i < 3; i++ ) {      
      idx[(i+1)%3] = this->getOpenMeshObject().opposite_he_opposite_vh( heh ).idx();
      heh = this->getOpenMeshObject().next_halfedge_handle(heh);
    }
    return idx;
  }

  // returns L_infty diam
  RealType getMeshBoundingBox ( aol::Vec3<RealType> &min, aol::Vec3<RealType> &max ) const;
  
  RealType computeBoundingBox ( aol::Vec3<RealType> &min, aol::Vec3<RealType> &max ) const{
    return getMeshBoundingBox( min, max );
  }

  // global mesh deformations
  void translate( const aol::Vec3<RealType>& offset );
  void shiftByOffset( const aol::Vec3<RealType>& offset ){ translate(offset); }
  void rotate( const aol::Matrix33<RealType>& RotationMatrix );  
  void scaleSizeByFactor( RealType factor );
  
  
  // normalize: if diam > 1.0 divide each point by diam, if translate: shift min to (0,0,0)
  void normalizeTo01( bool translate = false );
  void regularizeMesh ( );

  // returns weighted face normal
  void getFaceNormal ( const FaceHandle &fh, aol::Vec3<RealType> &wn, bool normalized = true ) const;
  void getFaceNormal ( const int idx, aol::Vec3<RealType> &wn, bool normalized = true ) const {
    getFaceNormal( faceHandle(idx), wn, normalized );
  }

  // returns Jacobian of parametrization x on element fh
  void getDx  ( const FaceHandle &fh,  aol::Mat<3,2,RealType> &Dx ) const;
  // get first fundamental form g of triangle fh (depends on parametrization!)
  void getFirstFundForm ( const FaceHandle &fh, aol::Matrix22<RealType> &g ) const;
  
  void get1RingVertices( const int idx, std::vector<int>& indices ) const;
  void get2RingVertices( const int idx, std::vector<int>& indices ) const;

  //---- conversion and IO ------------------------------------------------------------------------------

  void toVector ( aol::MultiVector<RealType> &dst ) const;
  void fromVector ( const aol::MultiVector<RealType> &src );

  void toVector ( aol::Vector<RealType> &dst ) const;
  void fromVector ( const aol::Vector<RealType> &src );

  // load routines
  void loadFromUDPLY ( const std::string filename );
  void loadFromPLY ( const std::string filename );
  void loadFromOBJ ( const std::string filename, bool planar = false );

  // save routines
  void saveAsPLY ( const std::string filename, int precision = 8 ) const {
    aol::MeshWithData<TriMesh<RealType> > ( *this, precision ).saveAsPLY ( filename );
  }
  void saveAsUDPLY ( const std::string filename ) const {
    aol::MeshWithData<TriMesh<RealType> > ( *this ).saveAsUDPLY ( filename );
  }
  void saveAsVTK( const std::string constfilename ) const {
    string filename = constfilename.substr( 0, constfilename.rfind(".") ) + ".vtk";
    aol::MeshWithData<TriMesh<RealType> >  ( *this ).saveAsLegacyVTK( filename );
  }
  
  void saveAsOBJ ( const std::string filename, int precision = 8, bool planar = false ) const;

  // converts this to aol::TriangMesh
  void convertToTriangMesh ( aol::TriangMesh<RealType> & triangMesh ) const;

  // converts aol::TriangMesh to TriMesh
  om::TriMesh<RealType>& importFromTriangMesh ( const aol::TriangMesh<RealType>& );
  
  void clearMesh() {
    this->clear();
  }

//! TODO implement base class of a general grid (not necessarily regular or cubic) which can be herited by qc::GridStructure and om::TriMesh
//!      and remove the subsequent functions afterwards!!!
  //---- compatibility with qc::GridStructure----------------------------------------------------------

  bool isAdaptive() const {
    throw aol::Exception ( "om::TriMesh::isAdaptive() is not implemented!", __FILE__, __LINE__ );
    return false;
  }

  bool isAdmissibleNode (const qc::CoordType & ) const {
    throw aol::Exception ( "om::TriMesh::isAdmissibleNode() is not implemented!", __FILE__, __LINE__ );
    return false;
  }

  int getGridDepth() const {
    throw aol::Exception ( "om::TriMesh::getGridDepth() is not implemented!", __FILE__, __LINE__ );
    return -1;
  }

  int checkForHangingNode ( const ElementType &, int ) const {
    throw aol::Exception ( "om::TriMesh::checkForHangingNode() is not implemented!", __FILE__, __LINE__ );
    return -1;
  }

  int getNumX() const {
    throw aol::Exception ( "om::TriMesh::getNumX() is not implemented!", __FILE__, __LINE__ );
    return -1;
  }

  int getNumY() const {
    throw aol::Exception ( "om::TriMesh::getNumY() is not implemented!", __FILE__, __LINE__ );
    return -1;
  }

  int getNumZ() const {
    throw aol::Exception ( "om::TriMesh::getNumZ() is not implemented!", __FILE__, __LINE__ );
    return -1;
  }

  aol::Vec3< int > getSize () const {
    throw aol::Exception ( "om::TriMesh::getSize() is not implemented!", __FILE__, __LINE__ );
    return aol::Vec3< int >( -1, -1, -1 );
  }


  qc::Dimension getDimOfWorld () const {
    return qc::QC_3D;
  }

  int getNumberOfBoundaryNodes () const{
    throw aol::Exception ( "om::TriMesh::getNumberOfBoundaryNodes() is not implemented!", __FILE__, __LINE__ );
    return -1;
  }

};

//=======================================================================================================
//
//  TriMesh with geometric properties
//
//-------------------------------------------------------------------------------------------------------
//! Provides functions to calculate geometric quantities like the shape operator or curvature.
//! Please use one of the derived classes TriMeshWith{Vertex, Edge}Normals to specify a notion of a normal field.
template <typename RealType, typename Imp>
class TriMeshWithGeomProps : public TriMesh<RealType>{

public:
  typedef typename TriMesh<RealType>::FaceHandle   FaceHandle;
  typedef typename TriMesh<RealType>::Element      ElementType;
  typedef typename TriMesh<RealType>::ElementIter  ElementIteratorType;

public:
  TriMeshWithGeomProps( ) : TriMesh<RealType>( ){}
  
  TriMeshWithGeomProps( const std::string filename ) : TriMesh<RealType>( filename ){}

  void getShapeOperator( const FaceHandle &fh, aol::Matrix22<RealType> &S ) const {
    aol::Matrix22<RealType> temp;
    this->getFirstFundForm(fh,temp);
    S.makeInverse(temp);
    this->asImp().getSecondFundForm(fh,temp);
    S*=temp;
  }

  //! interface function, has to be provided in derived classes.
  void getSecondFundForm( const FaceHandle &fh, aol::Matrix22<RealType> &h ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getSecondFundForm( fh, h );
  }

  // computes principal curvatures as the eigenvalues of the shape operator (in matrix notation)
  void getCurvature ( const FaceHandle &fh, aol::Vec2<RealType> &curv ) const;

  // computes principal curvatures (if not previously stored) and writes
  // them into pc, a MultiVector whose size is equal to our number of faces
  void writePrincipalCurvaturesTo ( aol::MultiVector<RealType> & pc );

  // exports mesh with curvatures as vtk-file
  void exportCurvature( const std::string filename );

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};

//=======================================================================================================
//
//  TriMeshWithVertexNormals
//
//-------------------------------------------------------------------------------------------------------
template <typename RealType>
class TriMeshWithVertexNormals : public TriMeshWithGeomProps< RealType, TriMeshWithVertexNormals<RealType> >{

public:
  typedef typename TriMesh<RealType>::FaceHandle FaceHandle;
  typedef typename TriMesh<RealType>::VertexHandle VertexHandle;

  typedef TriMeshWithGeomProps<RealType, TriMeshWithVertexNormals<RealType> > MeshType;
  typedef typename  MeshType::Element       ElementType;
  typedef typename  MeshType::ElementIter   ElementIteratorType;

public:
  TriMeshWithVertexNormals( ) : TriMeshWithGeomProps< RealType, TriMeshWithVertexNormals<RealType> >( ){}
  TriMeshWithVertexNormals( const std::string filename ) : TriMeshWithGeomProps< RealType, TriMeshWithVertexNormals<RealType> >( filename ){}

  // returns weighted vertex normal
  void getVertexNormal( const int &idx, aol::Vec3<RealType> &wn, bool normalized = true ) const {
      getVertexNormal( this->vertexHandle( idx ), wn, normalized );
  }
  void getVertexNormal ( const VertexHandle &vh, aol::Vec3<RealType> &wn, bool normalized = true ) const;
  // returns Jacobian of normal (depends on parametrization x!) on element fh
  void getDn ( const FaceHandle &fh, aol::Mat<3,2,RealType> &Dn ) const;
  //! returns second fundamental form (needed by base class!)
  void getSecondFundForm( const FaceHandle &fh, aol::Matrix22<RealType> &h ) const;

  // creates noise in normal direction
  void noise ( RealType factor, bool fixBoundary = false );

};

//=======================================================================================================
//
//  TriMeshWithEdgeNormals
//
//-------------------------------------------------------------------------------------------------------

template <typename RealType>
class TriMeshWithEdgeNormals : public TriMeshWithGeomProps<RealType, TriMeshWithEdgeNormals<RealType> >{

public:
  typedef typename TriMesh<RealType>::FaceHandle  FaceHandle;
  typedef typename TriMesh<RealType>::EdgeHandle  EdgeHandle;
  typedef typename TriMesh<RealType>::HalfedgeHandle  HalfedgeHandle;

  typedef TriMeshWithGeomProps<RealType, TriMeshWithEdgeNormals<RealType> > MeshType;
  typedef typename  MeshType::Element       ElementType;
  typedef typename  MeshType::ElementIter   ElementIteratorType;

public:
  TriMeshWithEdgeNormals( ) : TriMeshWithGeomProps<RealType, TriMeshWithEdgeNormals<RealType> >( ){ _edgeNormalProp = false; _edgeNormalCount = 0;}
  TriMeshWithEdgeNormals( const std::string filename ): TriMeshWithGeomProps<RealType, TriMeshWithEdgeNormals<RealType> >( filename ){ _edgeNormalProp = false; _edgeNormalCount = 0;}

  //! returns second fundamental form (needed by base class!)
  void getSecondFundForm( const FaceHandle &fh, aol::Matrix22<RealType> &h ) const;

  // uses pre-computed edge normals (i.e. one has to call requestEdgeNormals() first)
  void getSecondFundFormFast( const FaceHandle &fh, aol::Matrix22<RealType> &h ) const;
  // uses pre-computed edge normals (i.e. one has to call requestEdgeNormals() first)
  void getShapeOperatorFast( const FaceHandle &fh, aol::Matrix22<RealType> &S ) const;

  // routines to handle edge Normals as standard porperty
  void requestEdgeNormals ();
  void releaseEdgeNormals ();
  int  updateEdgeNormals ();
  bool hasEdgeNormals() const;

  // returns weighted edge normal
  void getEdgeNormal ( const FaceHandle &fh, int localCoord, aol::Vec3<RealType> &wn, bool normalized = true ) const;
  // returns weighted edge normal
  void getEdgeNormal ( const EdgeHandle &eh, aol::Vec3<RealType> &wn, bool normalized = true ) const;

private:
  // variables to treat edge normals as properties in Openmesh style
  bool _edgeNormalProp;
  int _edgeNormalCount;
  OpenMesh::EPropHandleT< aol::Vec3<RealType> > _edgeNormals;
};

//=======================================================================================================
//
//  CompTriMesh: extension of TriMesh which provides arithmetics, an inner product and a norm
//
//-------------------------------------------------------------------------------------------------------
template <typename RealType, typename TriMeshWithNormalType = TriMeshWithVertexNormals<RealType> >
class CompTriMesh : public TriMeshWithNormalType{

public:
  typedef typename TriMesh<RealType>::Point                      Point;
  typedef typename TriMeshWithNormalType::ElementType            ElementType;
  typedef typename TriMeshWithNormalType::ElementIteratorType    ElementIteratorType;
  typedef typename TriMesh<RealType>::NodeIteratorType           NodeIteratorType;
  typedef typename TriMesh<RealType>::BoundaryNodeIteratorType   BoundaryNodeIteratorType;

  typedef CompTriMesh<RealType, TriMeshWithNormalType> CompTriMeshType;


public:
  CompTriMesh( ) : TriMeshWithNormalType( ){}
  CompTriMesh( const std::string filename ) : TriMeshWithNormalType( filename ){}

  CompTriMesh( const CompTriMeshType& other, aol::CopyFlag copyFlag = aol::DEEP_COPY) : TriMeshWithNormalType( ) {

    if (copyFlag == aol::DEEP_COPY) {
      NodeIteratorType v_iter ( other );
        for ( ;  v_iter.notAtEnd();  ++v_iter )
          this->addVertex ( v_iter.getCoords() );

      ElementIteratorType e_iter ( other );
        for ( ;  e_iter.notAtEnd();  ++e_iter )
          this->addFace ( e_iter.getNodeIndices() );
    }
    else
      throw aol::Exception ( "om::CompTriMesh<RealType, TriMeshWithNormalType>::CompTriMesh( const CompTriMesh& other, CopyFlag copyFlag ): illegal copy flag", __FILE__, __LINE__ );
  }

  CompTriMesh( const TriMesh<RealType>& other, aol::CopyFlag copyFlag = aol::DEEP_COPY) {

    if (copyFlag == aol::DEEP_COPY) {
      NodeIteratorType v_iter ( other );
        for ( ;  v_iter.notAtEnd();  ++v_iter )
          this->addVertex ( v_iter.getCoords() );

      ElementIteratorType e_iter ( other );
        for ( ;  e_iter.notAtEnd();  ++e_iter )
          this->addFace ( e_iter.getNodeIndices() );
    }
    else
      throw aol::Exception ( "om::CompTriMesh<RealType, TriMeshWithNormalType>::CompTriMesh( const TriMesh<RealType>& other, CopyFlag copyFlag ): illegal copy flag", __FILE__, __LINE__ );
  }

  //---- mesh arithmetic ------------------------------------------------------------------------------

  void setZero(){
    (*this) *= aol::ZOTrait<RealType>::zero;
  }

  // Note: only makes sense for meshes with same connectivity and a one-to-one node correspondence
  CompTriMeshType &  operator+=(  const CompTriMeshType& other ){

    return addMultiple( other, aol::ZOTrait<RealType>::one );
  }

  // Note: only makes sense for meshes with same connectivity and a one-to-one node correspondence
  CompTriMeshType &  operator+=(  const aol::MultiVector<RealType>& other ){

    return addMultiple( other, aol::ZOTrait<RealType>::one );
  }

  // Note: only makes sense for meshes with same connectivity and a one-to-one node correspondence
  CompTriMeshType &  operator-=(  const CompTriMeshType& other ){

    return addMultiple( other, -aol::ZOTrait<RealType>::one );
  }

  // Note: only makes sense for meshes with same connectivity and a one-to-one node correspondence
  CompTriMeshType &  operator-=(  const aol::MultiVector<RealType>& other ){

    return addMultiple( other, -aol::ZOTrait<RealType>::one );
  }

  CompTriMeshType &  operator*=(  const RealType value ){

    for ( int i = 0; i < this->getNumVertices(); i++ ) {
      Point &p = this->point ( this->getOpenMeshObject().vertex_handle ( i ) );
      for ( int j = 0; j < 3; j++ )
        p[j] *= value;
    }
    return *this;
  }

  // Note: only makes sense for meshes with same connectivity and a one-to-one node correspondence
  CompTriMeshType&  addMultiple( const CompTriMeshType& other, const RealType value ){

    if ( this->getNumVertices() != other.getNumVertices() )
      throw aol::Exception ( "sizes don't match", __FILE__, __LINE__ );

    for ( NodeIteratorType iter = other; iter.notAtEnd(); ++iter ) {
      Point &p = this->point ( iter.vertexHandle( ) );
      aol::Vec3<RealType> q = iter.getCoords ();
      for ( int j = 0; j < 3; j++ )
        p[j] += value * q[j];
    }
    return *this;
  }

  // Note: only makes sense for meshes with same connectivity and a one-to-one node correspondence
  CompTriMeshType&  addMultiple( const aol::MultiVector<RealType>& other, const RealType value ){

    if ( this->getNumVertices() != other[0].getTotalSize()  )
      throw aol::Exception ( "sizes don't match", __FILE__, __LINE__ );

    for ( int i = 0; i < this->getNumVertices(); i++ ) {
      Point &p = this->point ( this->getOpenMeshObject().vertex_handle ( i ) );
      for ( int j = 0; j < 3; j++ )
        p[j] += value * other[j][i];
    }
    return *this;
  }

  // Note: only makes sense for meshes with same connectivity and a one-to-one node correspondence
  RealType dotProduct ( const CompTriMeshType& other ) const {

    if ( this->getNumVertices() != other.getNumVertices()  )
      throw aol::Exception ( "sizes don't match", __FILE__, __LINE__ );

    RealType dot = aol::ZOTrait<RealType>::zero;
    for ( NodeIteratorType iter = other; iter.notAtEnd(); ++iter ) {
      Point &p = this->point ( iter.vertexHandle( ) );
      aol::Vec3<RealType> q = iter.getCoords ();
      for ( int j = 0; j < 3; j++ )
        dot += p[j] * q[j];
    }
    return dot;
  }

  // norm
  RealType normSqr ( ) const {
    RealType norm = aol::ZOTrait<RealType>::zero;
    for ( NodeIteratorType iter = *this; iter.notAtEnd(); ++iter ) {
      const Point &p = this->point ( iter.vertexHandle( ) );
      norm += p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
    }
    return norm;
  }

  RealType norm ( ) const {
    return sqrt( normSqr() );
  }


};


template <typename _RealType>
inline std::ostream &operator<< ( std::ostream &os, const typename TriMesh<_RealType>::Element &t ) {
  return os << t.coords[0] << "    " << t.coords[1] << "    " << t.coords[2];
}

} // namespace om

#endif // __TRIMESH_H
