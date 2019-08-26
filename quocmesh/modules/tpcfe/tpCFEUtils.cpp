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

#include <tpCFEUtils.h>

namespace tpcfe {
template< class GridType, class FloatType >
void CFEInterfaceTriangulationGenerator<GridType, FloatType>::determineInterfaceTriangulation ( ) {

  if ( _interfaceTriangulationComputed )
    return;

  //  grid.detectVirtualNodes(); // is this necessary? not if always done before ...
  const FloatType grid_spacing = static_cast<FloatType> ( _grid.H() );

  for ( tpcfe::CFEInterfaceTriangleIterator< GridType > itit ( _grid ); itit.notAtEnd(); ++itit ) {
    aol::Vec3<FloatType>  FVertexCoords[4]; // coordinates of vertices (world coordinates) in Float

    for ( int i = 0; i < 4; ++i ) { // loop over vertices of tetra
      aol::Vec3<RealType> tmp;
      itit.getTetraRef().computeGlobalCoordinate ( tmp, itit.getElementRef(), i );
      for ( int j = 0; j < 3; ++j ) {
        FVertexCoords[i][j] =  static_cast<FloatType> ( tmp[j] ) * grid_spacing;
      }
    }

    aol::Vec<4,int> m; // vertex reindexing such that m[3] is the non-interface index
    m[0] = -1; m[1] = -1; m[2] = -1; m[3] = -1;
    if ( ! ( itit.getTetraRef().isVirtualNode ( 0 ) ) ) {
      m[0] = 1;   m[1] = 3;   m[2] = 2;  m[3] = 0;
    } else if ( ! ( itit.getTetraRef().isVirtualNode ( 1 ) ) ) {
      m[0] = 0;   m[1] = 2;   m[2] = 3;  m[3] = 1;
    } else if ( ! ( itit.getTetraRef().isVirtualNode ( 2 ) ) ) {
      m[0] = 0;   m[1] = 3;   m[2] = 1;  m[3] = 2;
    } else if ( ! ( itit.getTetraRef().isVirtualNode ( 3 ) ) ) {
      m[0] = 0;   m[1] = 1;   m[2] = 2;  m[3] = 3;
#ifdef DEBUG
    } else {
      throw aol::Exception ( "This should not happen", __FILE__, __LINE__ );
#endif
    }

    computeTetraData ( itit.getElementRef(), itit.getTetraRef(), FVertexCoords[m[0]], FVertexCoords[m[1]], FVertexCoords[m[2]], FVertexCoords[m[3]], m, false );

    TriangType new_triang = itit.getTriangRef().getVNIndices();

    // \todo think about whether it is necessary to flip some of the triangles
    if ( ! ( itit.getTetraRef().isVirtualNode ( 0 ) ) || ! ( itit.getTetraRef().isVirtualNode ( 2 ) ) ) {
      VNodeIndType tmp = new_triang[1];
      new_triang[1] = new_triang[2];
      new_triang[2] = tmp;
    }

    // store vertices
    for ( int tri = 0; tri < 3; ++tri ) {
      if ( _interfaceVertexIndexmapper.find ( new_triang[tri] ) == _interfaceVertexIndexmapper.end() ) {
        _interfaceVertexVector.push_back ( FVertexCoords[ m[tri] ] );
        computeVertexData ( new_triang[tri], _interfaceVertexDataVector );
        _interfaceVertexIndexmapper[ new_triang[tri] ] = _numInterfaceVertices;
        ++_numInterfaceVertices;
      }
    }

    // store triangle
    _interfaceTriangVector.push_back ( aol::Vec3<int> ( _interfaceVertexIndexmapper[ new_triang[0] ],  _interfaceVertexIndexmapper[ new_triang[1] ],  _interfaceVertexIndexmapper[ new_triang[2] ] ) );
    _interfaceTriangIndexmapper [ new_triang ] = this->_numInterfaceTriangs;
    ++_numInterfaceTriangs;
  }

  postprocessInterfaceVertexData();
  postprocessInterfaceTriangData();

  _interfaceTriangulationComputed = true;

} // method determineInterfaceTriangulation


template< class GridType, class FloatType >
void CFEInterfaceTriangulationGenerator<GridType, FloatType>::determineClippedBoundaryTriangulation ( const aol::Vec3<int> &clip_lower, const aol::Vec3<int> &clip_upper, const RealType geomTolerance, const signed char which ) {

  const FloatType grid_spacing = static_cast<FloatType> ( _grid.H() );

#ifdef VERBOSE
  cerr << "geometric clipping (image coordinates): "
       << "  clip_x_min = " << clip_lower[0] << ", clip_x_max = " << clip_upper[0]
       << ", clip_y_min = " << clip_lower[1] << ", clip_y_max = " << clip_upper[1]
       << ", clip_z_min = " << clip_lower[2] << ", clip_z_max = " << clip_upper[2] << endl;
#endif

  const qc::GridSize<qc::QC_3D> gridSize ( _grid );
  for ( typename GridType::FullElementIterator it ( _grid ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> curEl ( *it, gridSize, _grid.getElType ( *it ) );

    // only need to consider elements near the clipping planes
    if ( aol::Sqr ( curEl.x() - clip_lower[0] ) < 2 || aol::Sqr ( curEl.x() - clip_upper[0] ) < 2 ||
         aol::Sqr ( curEl.y() - clip_lower[1] ) < 2 || aol::Sqr ( curEl.y() - clip_upper[1] ) < 2 ||
         aol::Sqr ( curEl.z() - clip_lower[2] ) < 2 || aol::Sqr ( curEl.z() - clip_upper[2] ) < 2    ) {

      const CFEType cfeType = curEl.cfeType();

      // now we are interested in all inner tetrahedra

      if ( cfeType.representsInterfaced() )
        _grid.getCutRelations ( curEl ); // this fails for types 0x00 and 0xFF, but in this case there is no interface, so there are no cut relations.

      for ( tpcfe::CFETopoTetraIterator tit ( cfeType, which ); tit.notAtEnd(); ++tit ) {
        const tpcfe::CFETopoTetra tetra = *tit;

        aol::Vec3<RealType>    VertexCoords[4]; // coordinates of vertices (global, i. e. image coordinates)
        aol::Vec3<FloatType>  FVertexCoords[4]; // coordinates of vertices (global, i. e. image coordinates) in Float

        int glob_ind_vertex[4][2];
        bool vertexOnBoundary[4][3][2];

        for ( short i = 0; i < 4; ++i ) { // loop over vertices of tetra

          tetra.computeGlobalCoordinate ( VertexCoords[i], curEl, i );
          const int lIdx0 = tetra ( i, 0 );
          const int lIdx1 = tetra ( i, 1 );

          glob_ind_vertex[i][0] = curEl.globalIndex ( lIdx0 );
          glob_ind_vertex[i][1] = ( lIdx1 == NODE_NOT_VIRTUAL ? -1 : curEl.globalIndex ( lIdx1 ) ); // mapIndex ( a, -1 ) = mIC ( -1, a ) = a * 2^32

          const aol::Vec3<RealType> VC = VertexCoords[i];

          for ( short c = 0; c < 3; ++c ) {
            vertexOnBoundary[i][c][0] = ( aol::Sqr ( VC[c] - clip_lower[c] ) < geomTolerance );
            vertexOnBoundary[i][c][1] = ( aol::Sqr ( VC[c] - clip_upper[c] ) < geomTolerance );
          }

          for ( int j = 0; j < 3; ++j ) {
            FVertexCoords[i][j] =  static_cast<FloatType> ( VertexCoords[i][j] ) * grid_spacing;
          }

        } // for vertices of tetra

        aol::Vec<4,int> m[4];
        m[0][0] = 0;        m[0][1] = 1;        m[0][2] = 2;        m[0][3] = 3;
        m[1][0] = 0;        m[1][1] = 3;        m[1][2] = 1;        m[1][3] = 2;
        m[2][0] = 2;        m[2][1] = 3;        m[2][2] = 0;        m[2][3] = 1;
        m[3][0] = 2;        m[3][1] = 1;        m[3][2] = 3;        m[3][3] = 0;

        bool draw_triang[4];
        for ( short i = 0; i < 4; ++i ) {
          draw_triang[i] = ( ( vertexOnBoundary[ m[i][0] ][ 0 ][ 0 ] && vertexOnBoundary[ m[i][1] ][ 0 ][ 0 ] && vertexOnBoundary[ m[i][2] ][ 0 ][ 0 ] ) ||
                             ( vertexOnBoundary[ m[i][0] ][ 1 ][ 0 ] && vertexOnBoundary[ m[i][1] ][ 1 ][ 0 ] && vertexOnBoundary[ m[i][2] ][ 1 ][ 0 ] ) ||
                             ( vertexOnBoundary[ m[i][0] ][ 2 ][ 0 ] && vertexOnBoundary[ m[i][1] ][ 2 ][ 0 ] && vertexOnBoundary[ m[i][2] ][ 2 ][ 0 ] ) ||
                             ( vertexOnBoundary[ m[i][0] ][ 0 ][ 1 ] && vertexOnBoundary[ m[i][1] ][ 0 ][ 1 ] && vertexOnBoundary[ m[i][2] ][ 0 ][ 1 ] ) ||
                             ( vertexOnBoundary[ m[i][0] ][ 1 ][ 1 ] && vertexOnBoundary[ m[i][1] ][ 1 ][ 1 ] && vertexOnBoundary[ m[i][2] ][ 1 ][ 1 ] ) ||
                             ( vertexOnBoundary[ m[i][0] ][ 2 ][ 1 ] && vertexOnBoundary[ m[i][1] ][ 2 ][ 1 ] && vertexOnBoundary[ m[i][2] ][ 2 ][ 1 ] ) );
        }

        for ( int t = 0; t < 4; ++t ) {
          if ( draw_triang[t] ) {

            TriangType new_triang ( VNodeType::mapIndex ( glob_ind_vertex[ m[t][0] ][0] , glob_ind_vertex[ m[t][0] ][1] ),
                                    VNodeType::mapIndex ( glob_ind_vertex[ m[t][1] ][0] , glob_ind_vertex[ m[t][1] ][1] ),
                                    VNodeType::mapIndex ( glob_ind_vertex[ m[t][2] ][0] , glob_ind_vertex[ m[t][2] ][1] ) );

            for ( int i = 0; i < 3; ++i ) {
              if ( _boundaryVertexIndexmapper.find ( new_triang[i] ) == _boundaryVertexIndexmapper.end() ) {
                _boundaryVertexVector.push_back ( FVertexCoords[ m[t][i] ] );
                computeVertexData ( new_triang[i], _boundaryVertexDataVector );
                _boundaryVertexIndexmapper[ new_triang[i] ] = _numBoundaryVertices;
                ++_numBoundaryVertices;
              }
            }

            _boundaryTriangVector.push_back ( aol::Vec3<int> ( _boundaryVertexIndexmapper[ new_triang[0] ],  _boundaryVertexIndexmapper[ new_triang[1] ],  _boundaryVertexIndexmapper[ new_triang[2] ] ) );
            _boundaryTriangIndexmapper [ new_triang ] = this->_numBoundaryTriangs;
            ++_numBoundaryTriangs;

            // TODO: maybe float accuracy is insufficient here (matrices for computing gradients not invertible) ... think about this.

            computeTetraData ( curEl, tetra, FVertexCoords[ m[t][0] ], FVertexCoords[ m[t][1] ], FVertexCoords[ m[t][2] ], FVertexCoords[ m[t][3] ], m[t], true );

          }
        }

      } // tetra iterator
    } // interesting element
  } // element iterator

  postprocessBoundaryVertexData();
  postprocessBoundaryTriangData();

} // end determineClippedBoundaryTriangulation


template< class GridType, class FloatType >
void CFEInterfaceTriangulationGenerator<GridType, FloatType>::deformVertices ( VNodeIndexMapperType &IM, VertexVectorType &vertices, const qc::MultiArray<RealType, 3> &deformations, const aol::Vec3<FloatType> scaling ) {
  for ( typename VNodeIndexMapperType::iterator it = IM.begin(); it != IM.end(); ++it ) {
    int gIdx0 = 0, gIdx1 = 0;
    VNodeType::splitIndex ( it->first, gIdx0, gIdx1 );
    if ( gIdx1 == 0 ) {
      vertices[ it->second ] [0] += scaling[0] * deformations[0][gIdx0];
      vertices[ it->second ] [1] += scaling[1] * deformations[1][gIdx0];
      vertices[ it->second ] [2] += scaling[2] * deformations[2][gIdx0];
    } else {
      const VNodeType &vn = _grid.getVirtualNodeRef ( gIdx0, gIdx1 );
      aol::Vec3<RealType> localDisplacement;
      vn.extrapolate ( deformations, localDisplacement );
      localDisplacement[0] *= scaling[0];
      localDisplacement[1] *= scaling[1];
      localDisplacement[2] *= scaling[2];
      vertices[ it->second ] [0] += localDisplacement[0];
      vertices[ it->second ] [1] += localDisplacement[1];
      vertices[ it->second ] [2] += localDisplacement[2];
    }

  }
} // end deformVertices


template< class GridType, class FloatType >
void CFEInterfaceTriangulationGenerator<GridType, FloatType>::getVertexVector ( VertexVectorType &vertexVector ) {

  vertexVector.clear();
  vertexVector.reserve ( _interfaceVertexVector.size() + _boundaryVertexVector.size() ); // reserve memory, do not resize

  for ( typename VertexVectorType::iterator it = _interfaceVertexVector.begin(); it != _interfaceVertexVector.end(); ++it ) {
    vertexVector.push_back ( *it );
  }

  for ( typename VertexVectorType::iterator it = _boundaryVertexVector.begin(); it != _boundaryVertexVector.end(); ++it ) {
    vertexVector.push_back ( *it );
  }

} // end getVertexVector


template< class GridType, class FloatType >
void CFEInterfaceTriangulationGenerator<GridType, FloatType>::getTriangVector ( TriangVectorType &triangVector ) {

  triangVector.clear();
  triangVector.reserve ( _interfaceTriangVector.size() + _boundaryTriangVector.size() );

  for ( typename TriangVectorType::iterator it = _interfaceTriangVector.begin(); it != _interfaceTriangVector.end(); ++it ) {
    triangVector.push_back ( *it );
  }

  aol::Vec3<int> offset_triang ( _interfaceVertexVector.size(), _interfaceVertexVector.size(), _interfaceVertexVector.size() );

  for ( typename TriangVectorType::iterator it = _boundaryTriangVector.begin(); it != _boundaryTriangVector.end(); ++it ) {
    triangVector.push_back ( ( *it ) + offset_triang );
  }

} // end getTriangVector


// This method does not modify the index mappers
template< class GridType, class FloatType >
void CFEInterfaceTriangulationGenerator<GridType, FloatType>::getVertexDataVector ( VertexDataVectorType &vertexDataVector ) {

  vertexDataVector.clear();
  vertexDataVector.reserve ( _interfaceVertexDataVector.size() + _boundaryVertexDataVector.size() ); // reserve memory, do not resize

  for ( typename VertexDataVectorType::iterator it = _interfaceVertexDataVector.begin(); it != _interfaceVertexDataVector.end(); ++it ) {
    vertexDataVector.push_back ( *it );
  }

  for ( typename VertexDataVectorType::iterator it = _boundaryVertexDataVector.begin(); it != _boundaryVertexDataVector.end(); ++it ) {
    vertexDataVector.push_back ( *it );
  }

} // end getVertexDataVector


// This method does not modify the index mappers
template< class GridType, class FloatType >
void CFEInterfaceTriangulationGenerator<GridType, FloatType>::getTriangDataVector ( TriangDataVectorType &triangDataVector ) {

  triangDataVector.clear();
  triangDataVector.reserve ( _interfaceTriangDataVector.size() + _boundaryTriangDataVector.size() ); // reserve memory, do not resize

  for ( typename TriangDataVectorType::iterator it = _interfaceTriangDataVector.begin(); it != _interfaceTriangDataVector.end(); ++it ) {
    triangDataVector.push_back ( *it );
  }

  for ( typename TriangDataVectorType::iterator it = _boundaryTriangDataVector.begin(); it != _boundaryTriangDataVector.end(); ++it ) {
    triangDataVector.push_back ( *it );
  }

} // end getTriangDataVector



template< class GridType, class FloatType >
void CFEInterfaceTriangulationGenerator<GridType, FloatType>::saveToPLYFile ( const char* filename ) const {

  aol::Bzipofstream file ( filename );

  const int numberOfVertices = this->getNumberOfVertices(), numberOfTriangles = this->getNumberOfTriangles();

  // write UD-ply header
  file << "ply" << endl
       << "format ascii 1.0" << endl
       << "comment written by tpcfe::CFEInterfaceTriangulationGenerator<GridType,FloatType>::saveToPLYFile" << endl
       << "element vertex " << numberOfVertices << endl
       << "property float x" << endl
       << "property float y" << endl
       << "property float z" << endl;

  if ( _haveVertexData )
    file << "property float generic_scalar_values" << endl;

  file << "element face " <<  numberOfTriangles << endl
       << "property list uchar int vertex_index" << endl;

  if ( _haveTriangData )
    file << "property float generic_scalar_values" << endl;

  file << "end_header" << endl;

  for ( unsigned int i = 0; i < _interfaceVertexVector.size(); ++i ) {
    file << _interfaceVertexVector[i][0] << " " << _interfaceVertexVector[i][1] << " " << _interfaceVertexVector[i][2] ;
    if ( _haveVertexData )
      file << " " << _interfaceVertexDataVector[i];
    file << endl;
  }

  for ( unsigned int i = 0; i < _boundaryVertexVector.size(); ++i ) {
    file << _boundaryVertexVector[i][0] << " " << _boundaryVertexVector[i][1] << " " << _boundaryVertexVector[i][2] ;
    if ( _haveVertexData )
      file << " " << _boundaryVertexDataVector[i];
    file << endl;
  }

  for ( unsigned int j = 0; j < _interfaceTriangVector.size(); ++j ) {
    file << "3 " << _interfaceTriangVector[j][0] << " " << _interfaceTriangVector[j][1] << " " << _interfaceTriangVector[j][2];
    if ( _haveTriangData )
      file << " " << _interfaceTriangDataVector[j];
    file << endl;
  }

  const int offset = _interfaceVertexVector.size();
  for ( unsigned int j = 0; j < _boundaryTriangVector.size(); ++j ) {
    file << "3 " << _boundaryTriangVector[j][0] + offset << " " << _boundaryTriangVector[j][1] + offset << " " << _boundaryTriangVector[j][2] + offset;
    if ( _haveTriangData )
      file << " " << _boundaryTriangDataVector[j];
    file << endl;
  }

}

template< class GridType, class FloatType >
void CFEInterfaceTriangulationGenerator<GridType, FloatType>::writeToTriangMesh ( aol::TriangMesh<FloatType> &TriangMesh ) const {
  TriangMesh.clear();
  if ( _haveVertexData )
    TriangMesh.createVertexData ( 1 );

  if ( _haveTriangData )
    TriangMesh.createTriangData ( 1 );

  TriangMesh.reserve ( this->getNumberOfVertices(), this->getNumberOfTriangles() );

  for ( unsigned int i = 0; i < _interfaceVertexVector.size(); ++i ) {
    TriangMesh.pushBackVertex ( _interfaceVertexVector[i] );
    if ( _haveVertexData )
      TriangMesh.getVertexData(0)[ TriangMesh.getNumVertices() - 1 ] = _interfaceVertexDataVector[i];
  }

  for ( unsigned int i = 0; i < _boundaryVertexVector.size(); ++i ) {
    TriangMesh.pushBackVertex ( _boundaryVertexVector[i] );
    if ( _haveVertexData )
      TriangMesh.getVertexData(0)[ TriangMesh.getNumVertices() - 1 ] = _boundaryVertexDataVector[i];
  }

  for ( unsigned int j = 0; j < _interfaceTriangVector.size(); ++j ) {
    TriangMesh.pushBackTriang ( _interfaceTriangVector[j] );
    if ( _haveTriangData )
      TriangMesh.getTriangData(0)[ TriangMesh.getNumTriangs() - 1 ] = _interfaceTriangDataVector[j];
  }

  const aol::Vec3<int> offset ( _interfaceVertexVector.size(), _interfaceVertexVector.size(), _interfaceVertexVector.size() );
  for ( unsigned int j = 0; j < _boundaryTriangVector.size(); ++j ) {
    TriangMesh.pushBackTriang ( _boundaryTriangVector[j] + offset );
    if ( _haveTriangData )
      TriangMesh.getTriangData(0)[ TriangMesh.getNumTriangs() - 1 ] = _boundaryTriangDataVector[j];
  }
}


template< class GridType, typename Float_Type >
void CFEInterfaceTriangulationGenerator< GridType, Float_Type >::determineSliceTriangulation ( const qc::Comp direction, const short slice, const RealType geomTolerance, const signed char sign ) {

  const FloatType gridSpacing = static_cast<FloatType> ( this->_grid.H() );

  const qc::GridSize<qc::QC_3D> gridSize ( _grid );
  for ( typename GridType::FullElementIterator it ( _grid ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> curEl ( *it, gridSize, _grid.getElType ( *it ) );

    // only need to consider elements on one side of the plane; this in particular means we cannot display the "last" slice.
    if ( curEl[direction] == slice ) {
      const CFEType cfeType = curEl.cfeType();

      if ( cfeType.representsInterfaced() )
        this->_grid.getCutRelations ( curEl );

      signed char usedSign = sign;
      if ( ( GridType::CT == tpcfe::CFE_CDWI_ELAST ) &&  this->_grid.hasDomain() ) {
        const int structureNo =  curEl.cfeType()._structureNo;
        if ( curEl.cfeType().representsInterfaced() ) {
          if ( structureNo == tpcfe::MAX_STRUCT_ID ) {
            // at domain boundary, only loop over inner tet
            usedSign = -1;
          } else {
            // at interface (between coefficient domains) boundary
            usedSign = 0;
          }
        } else {
          // not interfaced: pureType refers to domain, so all inner will use whichever coefficient domain
          usedSign = -1;
        }
      } else {
        usedSign = 0; // loop over all tet if there is no domain
      }

      for ( tpcfe::CFETopoTetraIterator tit ( cfeType, usedSign ); tit.notAtEnd(); ++tit ) {
        const tpcfe::CFETopoTetra tetra = *tit;

        aol::Vec3<RealType>    VertexCoords[4]; // coordinates of vertices (global, i. e. image coordinates)
        aol::Vec3<FloatType>  FVertexCoords[4]; // coordinates of vertices (global, i. e. image coordinates) in Float

        int glob_ind_vertex[4][2];
        bool vertexOnSlice[4];

        for ( short i = 0; i < 4; ++i ) { // loop over vertices of tetra

          tetra.computeGlobalCoordinate ( VertexCoords[i], curEl, i );
          const int lIdx0 = tetra ( i, 0 );
          const int lIdx1 = tetra ( i, 1 );

          glob_ind_vertex[i][0] = curEl.globalIndex ( lIdx0 );
          glob_ind_vertex[i][1] = ( lIdx1 == NODE_NOT_VIRTUAL ? -1 : curEl.globalIndex ( lIdx1 ) );

          for ( short c = 0; c < 3; ++c ) {
            FVertexCoords[i][c] =  static_cast<FloatType> ( VertexCoords[i][c] ) * gridSpacing;
          }

          vertexOnSlice[i] = ( aol::Sqr ( FVertexCoords[i][ direction ] - gridSpacing * slice ) < geomTolerance );

        } // for vertices of tetra


        aol::Vec<4,int> m[4];
        m[0][0] = 0;        m[0][1] = 1;        m[0][2] = 2;        m[0][3] = 3;
        m[1][0] = 0;        m[1][1] = 3;        m[1][2] = 1;        m[1][3] = 2;
        m[2][0] = 2;        m[2][1] = 3;        m[2][2] = 0;        m[2][3] = 1;
        m[3][0] = 2;        m[3][1] = 1;        m[3][2] = 3;        m[3][3] = 0;

        bool draw_triang[4];
        for ( short t = 0; t < 4; ++t ) {
          draw_triang[t] = ( vertexOnSlice[ m[t][0] ] && vertexOnSlice[ m[t][1] ] && vertexOnSlice[ m[t][2] ] );
        }

        for ( short t = 0; t < 4; ++t ) {
          if ( draw_triang[t] ) {
            TriangType new_triang ( VNodeType::mapIndex ( glob_ind_vertex[ m[t][0] ][0] , glob_ind_vertex[ m[t][0] ][1] ),
                                    VNodeType::mapIndex ( glob_ind_vertex[ m[t][1] ][0] , glob_ind_vertex[ m[t][1] ][1] ),
                                    VNodeType::mapIndex ( glob_ind_vertex[ m[t][2] ][0] , glob_ind_vertex[ m[t][2] ][1] ) );

            for ( short i = 0; i < 3; ++i ) {
              if ( _boundaryVertexIndexmapper.find ( new_triang[i] ) == _boundaryVertexIndexmapper.end() ) {
                _boundaryVertexVector.push_back ( FVertexCoords[ m[t][i] ] );
                computeVertexData ( new_triang[i], _boundaryVertexDataVector );
                _boundaryVertexIndexmapper[ new_triang[i] ] = _numBoundaryVertices;
                ++_numBoundaryVertices;
              }
            }
            _boundaryTriangVector.push_back ( aol::Vec3<int> ( _boundaryVertexIndexmapper[ new_triang[0] ],  _boundaryVertexIndexmapper[ new_triang[1] ],  _boundaryVertexIndexmapper[ new_triang[2] ] ) );
            _boundaryTriangIndexmapper [ new_triang ] = _numBoundaryTriangs;
            ++_numBoundaryTriangs;

            computeTetraData ( curEl, tetra, FVertexCoords[ m[t][0] ], FVertexCoords[ m[t][1] ], FVertexCoords[ m[t][2] ], FVertexCoords[ m[t][3] ], m[t], true );

          }
        }
      }
    }
  }

  this->postprocessBoundaryVertexData();
  this->postprocessBoundaryTriangData();

}

template< class GridType >
void CFEInterfaceTriangulationWithVonMisesStressGeneratorBase< GridType >::computeTetraData ( const tpcfe::CFEElement<RealType> &el, const tpcfe::CFETopoTetra &tetra,
                                                                                              const aol::Vec3<FloatType> &v0, const aol::Vec3<FloatType> &v1,
                                                                                              const aol::Vec3<FloatType> &v2, const aol::Vec3<FloatType> &v3,
                                                                                              const aol::Vec<4,int> &m, const bool onBoundary ) {
  aol::Matrix33<RealType> directionsMatrix;

  directionsMatrix.setRow ( 0, v0 - v3 );
  directionsMatrix.setRow ( 1, v1 - v3 );
  directionsMatrix.setRow ( 2, v2 - v3 );

  directionsMatrix *= this->_aleph; // conversion to SI units

  aol::Matrix33<RealType> tetraMatrix;

  if ( directionsMatrix.det() == 0.0 ) {

    // TODO: think about what to do in this case -- maybe need double accuracy?
    cerr << "Trying to computeTetraData on degenerate tetrahedron, directions are: " << directionsMatrix << endl << endl;
    tetraMatrix.setIdentity();

  } else {

    tetraMatrix = directionsMatrix.inverse();

  }

  aol::Vec3<RealType> deformation_at_vertex[4];
  uint64_t nodeIndex[4];

  for ( short i = 0; i < 4; ++i ) { // loop over vertices of tetra

    const int lIdx0 = tetra ( i, 0 ), lIdx1 = tetra ( i, 1 );
    const int gIdx0 = el.globalIndex ( lIdx0 ), gIdx1 = ( lIdx1 == NODE_NOT_VIRTUAL ? -1 : el.globalIndex ( lIdx1 ) );

    if ( lIdx1 == NODE_NOT_VIRTUAL ) {
      deformation_at_vertex[i] = this->_deformation.get(gIdx0);
    } else {
      this->_grid.getVirtualNodeRef ( gIdx0, gIdx1 ).extrapolate ( this->_deformation, deformation_at_vertex[i] );
    }
    nodeIndex[i] = VNodeType::mapIndex ( gIdx0, gIdx1 );
  }

  aol::Vec3<RealType> difference_deformations[3];
  aol::Matrix33<RealType> gradient_deformations, sigma;
  for ( int i = 0; i < 3; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      difference_deformations[j][i] = deformation_at_vertex[ m[i] ][j] - deformation_at_vertex[ m[3] ][j];
    }
  }

  for ( int i = 0; i < 3; ++i ) { // here, we compute the gradient
    gradient_deformations[i] = tetraMatrix * difference_deformations[i];
  }

  computeStress ( gradient_deformations, el, tetra, sigma );

  // this is the definition of von Mises stress
  const RealType vMStress = sqrt ( (  1.0 ) * ( aol::Sqr ( sigma[0][0] )  + aol::Sqr ( sigma[1][1] )  + aol::Sqr ( sigma[2][2] ) )  +
                                   ( -1.0 ) * ( sigma[0][0] * sigma[1][1] + sigma[0][0] * sigma[2][2] + sigma[1][1] * sigma[2][2] ) +
                                   (  3.0 ) * ( sigma[0][1] * sigma[0][1] + sigma[0][2] * sigma[0][2] + sigma[1][2] * sigma[1][2] )   );


  TriangDataMapType &dataCollector = ( onBoundary ? this->_boundaryTriangDataCollector : this->_interfaceTriangDataCollector );
  dataCollector[ TriangType ( nodeIndex[ m[0] ], nodeIndex[ m[1] ], nodeIndex[ m[2] ] ) ] = vMStress;
}



template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< float, tpcfe::CFE_NONE, float >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< double, tpcfe::CFE_NONE, double >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< double, tpcfe::CFE_NONE, double >, double >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< long double, tpcfe::CFE_NONE, long double >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< long double, tpcfe::CFE_NONE, long double >, long double >;

template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< float, tpcfe::CFE_CD, float >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< double, tpcfe::CFE_CD, double >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< long double, tpcfe::CFE_CD, long double >, float >;

template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< float, tpcfe::CFE_TPOS, float >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< double, tpcfe::CFE_TPOS, double >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< long double, tpcfe::CFE_TPOS, long double >, float >;

template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< float, tpcfe::CFE_CDWI_TPOS, float >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< double, tpcfe::CFE_CDWI_TPOS, double >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< double, tpcfe::CFE_CDWI_TPOS, double >, double >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< long double, tpcfe::CFE_CDWI_TPOS, long double >, float >;

template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< float, tpcfe::CFE_TPOSELAST, tpcfe::IsotropicElasticityCoefficient<float> >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< double, tpcfe::CFE_TPOSELAST, tpcfe::IsotropicElasticityCoefficient<double> >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< long double, tpcfe::CFE_TPOSELAST, tpcfe::IsotropicElasticityCoefficient<long double> >, float >;

template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< float, tpcfe::CFE_TPOSELAST, tpcfe::VoigtElasticityCoefficient<float> >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< double, tpcfe::CFE_TPOSELAST, tpcfe::VoigtElasticityCoefficient<double> >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< long double, tpcfe::CFE_TPOSELAST, tpcfe::VoigtElasticityCoefficient<long double> >, float >;

template class CFEInterfaceTriangulationWithVonMisesStressGeneratorBase< tpcfe::CFEGrid< float, tpcfe::CFE_CD, float > >;
template class CFEInterfaceTriangulationWithVonMisesStressGeneratorBase< tpcfe::CFEGrid< double, tpcfe::CFE_CD, double > >;
template class CFEInterfaceTriangulationWithVonMisesStressGeneratorBase< tpcfe::CFEGrid< long double, tpcfe::CFE_CD, long double > >;

template class CFEInterfaceTriangulationWithVonMisesStressGeneratorBase< tpcfe::CFEGrid< float, tpcfe::CFE_TPOSELAST, tpcfe::IsotropicElasticityCoefficient<float> > >;
template class CFEInterfaceTriangulationWithVonMisesStressGeneratorBase< tpcfe::CFEGrid< double, tpcfe::CFE_TPOSELAST, tpcfe::IsotropicElasticityCoefficient<double> > >;
template class CFEInterfaceTriangulationWithVonMisesStressGeneratorBase< tpcfe::CFEGrid< long double, tpcfe::CFE_TPOSELAST, tpcfe::IsotropicElasticityCoefficient<long double> > >;

  // this is only partially useful
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< double, tpcfe::CFE_CDWI_ELAST, tpcfe::IsotropicElasticityCoefficient<double> >, float >;
template class CFEInterfaceTriangulationGenerator < tpcfe::CFEGrid< double, tpcfe::CFE_CDWI_ELAST, tpcfe::VoigtElasticityCoefficient<double> >, float >;

  // end namespace
}
