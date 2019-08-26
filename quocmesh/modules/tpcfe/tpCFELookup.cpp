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

#include <tpCFELookup.h>
#include <tpCFEElement.h>

namespace tpcfe {

// all tetra have positive orientation.
const short CFETopoLookup::_stdTetraVertex[6][4] = { { 0, 1, 2, 6 },
                                                     { 1, 0, 5, 6 },
                                                     { 0, 4, 5, 6 },
                                                     { 2, 1, 3, 7 },
                                                     { 1, 2, 6, 7 },
                                                     { 5, 1, 6, 7 } };

const short CFETopoLookup::_stdTetraEdges[6][12] =  { { 0, 1, 0, 2, 0, 6, 1, 2, 1, 6, 2, 6 },
                                                      { 0, 1, 0, 5, 0, 6, 1, 5, 1, 6, 5, 6 },
                                                      { 0, 4, 0, 5, 0, 6, 4, 5, 4, 6, 5, 6 },
                                                      { 1, 2, 1, 3, 1, 7, 2, 3, 2, 7, 3, 7 },
                                                      { 1, 2, 1, 6, 1, 7, 2, 6, 2, 7, 6, 7 },
                                                      { 1, 5, 1, 6, 1, 7, 5, 6, 5, 7, 6, 7 } };

// The diagonal (triangular faces) faces, indices are ordered positive
const short CFETopoLookup::_elementFaces[11][4] = {  { 1, 6, 2, -1},
                                                     { 1, 7, 2, -1},
                                                     { 0, 5, 6, -1},
                                                     { 1, 5, 6, -1},
                                                     { 0, 1, 6, -1},

                                                     { 0, 1, 3, 2},   // The quadratic faces
                                                     { 4, 5, 7, 6},
                                                     { 0, 4, 6, 2},
                                                     { 1, 5, 7, 3},
                                                     { 0, 1, 5, 4},
                                                     { 2, 3, 7, 6} };

const aol::Vec3<short> CFETopoLookup::_tetNbOffsets[ NUM_NEIGHBORS ] = { aol::Vec3<short> (  0, -1, -1 ),   // 15 neighboring nodes connected by the edges of the six standard tetrahedra
                                                                         aol::Vec3<short> (  1, -1, -1 ),
                                                                         aol::Vec3<short> ( -1,  0, -1 ),
                                                                         aol::Vec3<short> (  0,  0, -1 ),
                                                                         aol::Vec3<short> (  0, -1,  0 ),
                                                                         aol::Vec3<short> (  1, -1,  0 ),
                                                                         aol::Vec3<short> ( -1,  0,  0 ),
                                                                         aol::Vec3<short> (  0,  0,  0 ),   // _tetNbOffset[ NEIGHBOR_SELF + i ] = - _tetNbOffset[ NEIGHBOR_SELF - i ]
                                                                         aol::Vec3<short> (  1,  0,  0 ),
                                                                         aol::Vec3<short> ( -1,  1,  0 ),
                                                                         aol::Vec3<short> (  0,  1,  0 ),
                                                                         aol::Vec3<short> (  0,  0,  1 ),
                                                                         aol::Vec3<short> (  1,  0,  1 ),
                                                                         aol::Vec3<short> ( -1,  1,  1 ),
                                                                         aol::Vec3<short> (  0,  1,  1 ) };



const aol::Vec3<short> CFETopoLookup::_hexNbOffsets[ 8 ] = { aol::Vec3<short> (  0, 0, 0 ),
                                                             aol::Vec3<short> (  1, 0, 0 ),
                                                             aol::Vec3<short> (  0, 1, 0 ),
                                                             aol::Vec3<short> (  1, 1, 0 ),
                                                             aol::Vec3<short> (  0, 0, 1 ),
                                                             aol::Vec3<short> (  1, 0, 1 ),
                                                             aol::Vec3<short> (  0, 1, 1 ),
                                                             aol::Vec3<short> (  1, 1, 1 ) };

const short CFETopoLookup::_cubeEdges[12][2] = { {0, 1}, {1, 3}, {3, 2}, {0, 2}, {4, 5}, {5, 7}, {7, 6}, {6, 4}, {0, 4}, {1, 5}, {3, 7}, {2, 6} };

const short CFETopoLookup::_tetraEdges[ NUM_TETRA_EDGES_IN_CUBE ][2] = { {0, 1}, {0, 2}, {1, 2}, {1, 3}, {2, 3}, {0, 4}, {0, 5}, {1, 5}, {4, 5}, {0, 6}, {1, 6}, {2, 6}, {4, 6}, {5, 6}, {1, 7}, {2, 7}, {3, 7}, {5, 7}, {6, 7} };

const short CFETopoLookup::_edge[ NUM_TETRA_EDGES_IN_CUBE ] =          {     3,       5,      6,     10,     12,     17,     33,     34,     48,     65,     66,     68,     80,     96,    130,    132,    136,    160,    192 };

const bool CFETopoLookup::_admissibleSignature[ 256 ] = { 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1,
                                                          1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1,
                                                          1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1,
                                                          1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1,
                                                          1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1,
                                                          1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1,
                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                          1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
                                                          1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1,
                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                          1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1,
                                                          1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1,
                                                          1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1,
                                                          1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1,
                                                          1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1,
                                                          1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1 };

const CFECube *CFETopoLookup::_topo[256];


template <typename RealType>
int CFELookup<RealType>::refCounter = 0;

template <typename RealType>
int CFELookup<RealType>::mem = 0;

template <typename RealType>
const RealType CFELookup<RealType>::_stdTetraVolume = 1. / 6.;

template <typename RealType>
const aol::Vec3<RealType> CFELookup<RealType>::_hexNbOffs[8] = { aol::Vec3<RealType> ( 0, 0, 0 ),    // 0
                                                                 aol::Vec3<RealType> ( 1, 0, 0 ),    // 1
                                                                 aol::Vec3<RealType> ( 0, 1, 0 ),    // 2
                                                                 aol::Vec3<RealType> ( 1, 1, 0 ),    // 3
                                                                 aol::Vec3<RealType> ( 0, 0, 1 ),    // 4
                                                                 aol::Vec3<RealType> ( 1, 0, 1 ),    // 5
                                                                 aol::Vec3<RealType> ( 0, 1, 1 ),    // 6
                                                                 aol::Vec3<RealType> ( 1, 1, 1 ) };  // 7


template <typename RealType>
const RealType CFELookup<RealType>::_localMassMatrixRefTetra[4][4] = { { 1. /  60., 1. / 120., 1. / 120., 1. / 120. },
                                                                       { 1. / 120., 1. /  60., 1. / 120., 1. / 120. },
                                                                       { 1. / 120., 1. / 120., 1. /  60., 1. / 120. },
                                                                       { 1. / 120., 1. / 120., 1. / 120., 1. /  60. } };

template <typename RealType>
const RealType CFELookup<RealType>::_localStiffMatrixRefTetra[4][4] = { { + 3. / 6., - 1. / 6., - 1. / 6., - 1. / 6. },
                                                                        { - 1. / 6., + 1. / 6.,        0.,        0. },
                                                                        { - 1. / 6.,        0., + 1. / 6.,        0. },
                                                                        { - 1. / 6.,        0.,        0., + 1. / 6. } };


//! gradients of the standard basis functions on the reference tetrahedron
template <typename RealType>
const RealType CFELookup<RealType>::_gradientStdTetraBasisFct[4][3] = { { -1., -1., -1. },
                                                                        {  1.,  0.,  0. },
                                                                        {  0.,  1.,  0. },
                                                                        {  0.,  0.,  1. } };

template <typename RealType>
aol::Vec3<RealType> CFELookup<RealType>::_stdTetraBasisGradients[6][4];

template <typename RealType>
RealType CFELookup<RealType>::_localStiffnessMatrixStdTetra[6][4][4];

template <typename RealType>
aol::Vec3<RealType> CFELookup<RealType>::_faceNormal[11];

template <typename RealType>
bool CFELookup<RealType>::_initialized = false;


template <typename RealType>
void CFELookup<RealType>::createStiffnessForStdTetra() {
  const vector <CFETopoTetra > &tv = CFETopoLookup::_topo[0]->getTopoTetraVec();
  const int size = tv.size();
  CFEElement<RealType>     element ( qc::GridSize<qc::QC_3D>( 1, 1, 1 ), 0, 0, 0, CFEType() );

  // Run over tetra of non-interfaced element
  for ( int i = 0; i < size; ++i ) {
    CFETetra<RealType> t ( tv[i] );

    element.getTetraVertices ( t );

    t.computeInverseTransformation();
    t.computeStiffnessMatrix();

    const RealType volume = static_cast< RealType > ( 1.0 / 6.0 );

    for ( int a = 0; a < 4; ++a )
      for ( int b = 0; b < 4; ++b ) {
        // Because for the standard tetrahedra the determinant is always 1 we
        // do not have to divide by it
        _localStiffnessMatrixStdTetra[i][a][b] = t._lsm.get ( a, b ) * volume; // / det2;
      }
  }
}

#if 0
template <typename RealType>
struct SortStruct {
  RealType value;
  int id;
};

template <typename RealType>
int compare ( const void *a, const void *b ) {
  const SortStruct<RealType> *A = static_cast<const SortStruct<RealType>* > ( a );
  const SortStruct<RealType> *B = static_cast<const SortStruct<RealType>* > ( b );

  const RealType diff = A->value - B->value;
  if ( diff > 0 ) return -1;
  if ( diff == 0 ) return 0;
  if ( diff < 0 ) return 1;
  return 0;
}
#endif


template <typename RealType>
void CFELookup<RealType>::computeNormals() {
  for ( int i = 0; i < 11; ++i ) {  // For each face
    RealType vertex[3][3];
    for ( int v = 0; v < 3; ++v ) {  // For each vertex
      const int index = _elementFaces[i][v];
      for ( int j = 0; j < 3; ++j ) {  // For each coordinate
        vertex[v][j] = _hexNbOffs[index][j];
        if ( v > 0 ) vertex[v][j] -= vertex[0][j];
      }
    }

    // Compute cross-product of edges
    _faceNormal[i][0] =     vertex[1][1] * vertex[2][2] - vertex[2][1] * vertex[1][2];
    _faceNormal[i][1] = - ( vertex[1][0] * vertex[2][2] - vertex[2][0] * vertex[1][2] );
    _faceNormal[i][2] =     vertex[1][0] * vertex[2][1] - vertex[2][0] * vertex[1][1];

    // normalize normals
    const RealType l = sqrt ( aol::Sqr ( _faceNormal[i][0] ) + aol::Sqr ( _faceNormal[i][1] ) + aol::Sqr ( _faceNormal[i][2] ) );
    for ( int k = 0; k < 3; ++k ) _faceNormal[i][k] /= l;
  }
}

template <typename RealType>
void CFELookup<RealType>::createBasisGradientsForStdTetra() {
  // Loop over all std tetra
  const vector < CFETopoTetra > &tv = CFETopoLookup::_topo[0]->getTopoTetraVec();
  const int size = tv.size();
  CFEElement<RealType>     element ( qc::GridSize<qc::QC_3D>( 1, 1, 1 ), 0, 0, 0, CFEType() );

  // Run over tetra of non-interfaced element
  for ( int i = 0; i < size; ++i ) {
    CFETetra<RealType> t ( tv[i] );

    element.getTetraVertices ( t );
    t.computeInverseTransformation();

#ifdef VERBOSE
    cerr << "---------------- type " << i << endl;
    cerr << "tetra is " << t << endl;
    // print vertices?

    cerr << "barycentric coordinates are" << endl;

    for ( int j = 0; j < 3; ++j ) {
      for ( int k = 0; k < 3; ++k ) {
        cerr << t._barycentricCoordinates[j][k] << ", ";
      }
      cerr << endl;
    }
#endif

    for ( int k = 0; k < 3; ++k ) { // For each coordinate
      RealType value = 0.0;
      for ( int j = 1; j < 4; ++j ) {  // For each basis function
        const RealType tmp = t._barycentricCoordinates[j-1][k];

        _stdTetraBasisGradients[i][j][k] = tmp;
        value -= tmp;
      }
      _stdTetraBasisGradients[i][0][k] = value;
    }

#ifdef VERBOSE
    cerr << "Tetra " << t << endl;
    cerr << "BasisFct Gradients = " << endl;
    for ( int j = 0; j < 4; ++j ) {
      cerr << _stdTetraBasisGradients[i][j] << endl;
    }
    cerr << "<-+ Key!" << endl;
#endif

#ifdef VERBOSE
    // check whether gradient is orthogonal to face opposite to node.
    /*
      for ( int j = 0; j < 4; ++j ) {
        // get indices of opposite face's vertices
        int i1 = t ( ( j + 1 ) % 4, 0 );
        int i2 = t ( ( j + 2 ) % 4, 0 );
        int i3 = t ( ( j + 3 ) % 4, 0 );

        cerr << "Face opposite to node " << t ( j, 0 ) << " has vertex indices " << i1 << ", " << i2 << ", " << i3 << endl;
        aol::Vec3<RealType> c  = t.vertex[ ( j+1 ) %4]; // FIXME: vertex no longer exists
        aol::Vec3<RealType> v1 = t.vertex[ ( j+2 ) %4];
        aol::Vec3<RealType> v2 = t.vertex[ ( j+3 ) %4];

        cerr << "Vertices are " << c << " and " << v1 << " and " << v2 << endl;

        v1 -= c;  v2 -= c;
        cerr << "Face is spanned by vectors " << v1 << " and " << v2 << endl;
        cerr << "Scalar products are : " << aol::color::red;
        cerr << _stdTetraBasisGradients[i][j] * v1 << ", ";
        cerr << _stdTetraBasisGradients[i][j] * v2 << aol::color::black << endl << endl;
      }
    */

    // Check whether sum of all gradients equals to zero
    aol::Vec3<RealType> sum;

    for ( int j = 0; j < 4; ++j ) {
      sum += _stdTetraBasisGradients[i][j];
    }
    cerr << "Sum of all gradients is " << sum << endl;

    getchar();
#endif
  }
}

template< typename RealType >
void CFELookup<RealType>::writeRefTetraVecTo ( std::vector< tpcfe::CFETetra<RealType> > &refTetraVec ) {
  refTetraVec.clear();
  refTetraVec.reserve( 6 );
  for ( CFETopoTetraIterator it ( CFEType ( -1, 0 ) ); it.notAtEnd(); ++it ) {
    tpcfe::CFETetra<RealType> tetra( *it );
    for ( short i = 0; i < 3; ++i ) {
      tetra._edge[i] = CFELookup<RealType>::_hexNbOffs[ (*it)( 1+i, 0 ) ] - CFELookup<RealType>::_hexNbOffs[ (*it)( 0, 0 ) ];
    }
    tetra.computeInverseTransformation();
    refTetraVec.push_back ( tetra );
  }
}

template <typename RealType>
int CFELookup<RealType>::construct () {
#ifdef _OPENMP
#pragma omp critical (tpcfe_CFELookup)
#endif
  {
    if ( refCounter == 0 ) {
      mem = 0;

      // Create all sample elements
      for ( int i = 0; i < 256; ++i ) {
        _topo[i] = new CFECube ( i );
        mem += _topo[i]->memoryUsed();
        mem += sizeof ( CFECube* );
      }

      // Create gradients of basis functions on std tetrahedra
      createBasisGradientsForStdTetra();

      // Create stiffness matrices for std tetrahedra
      createStiffnessForStdTetra();

      // Compute the normals of the outer and inner faces of an element
      computeNormals();

      _initialized = true;
    }

    ++refCounter;
  }

  return ( mem );
}


template <typename RealType>
void CFELookup<RealType>::destruct () {
#ifdef _OPENMP
#pragma omp critical (tpcfe_CFELookup)
#endif
  {
    if ( refCounter == 1 ) {
      for ( int i = 0; i < 256; ++i ) delete ( _topo[i] );
      mem = 0;
      _initialized = false;
    }
    --refCounter;
  }
}


CFECube::CFECube ( int id ) : _interfaced ( 256 ), _tetraVec ( 0 ) {
  for ( int i = 0; i < NUM_TETRA_EDGES_IN_CUBE; i++ ) {
    const short ei = tpcfe::CFETopoLookup::_edge[i];
    const short tmp = ei & id;
    if ( tmp == 0 || tmp == ei )  // Edge fully inside or fully outside structure
      _interfaced[ei] = false;
    else
      _interfaced[ei] = true;
  }
  for ( int i = 0; i < 6; i++ ) {
    splitInTetra ( i, id );
  }
  makeOrientation();
}

const vector < CFETopoTetra > & CFECube::getTopoTetraVec ( ) const {
  return ( _tetraVec );
}

int CFECube::memoryUsed() const {
  return  ( sizeof ( CFECube ) +
            sizeof ( CFETopoTetra ) * _tetraVec.size() );
}

template <typename RealType>
void CFECube::getTetraEdges ( const CFETopoTetra &t, aol::Vec3<RealType> edge[3] ) {
  aol::Vec3<RealType> vertex[4];
  for ( int l = 0; l < 4; ++l ) {
    const short a = t ( l, 0 );
    const short b = t ( l, 1 );
    if ( b == NODE_NOT_VIRTUAL ) { // Is not a virtual node
      vertex[l][0] = CFELookup<RealType>::_hexNbOffs[a][0];
      vertex[l][1] = CFELookup<RealType>::_hexNbOffs[a][1];
      vertex[l][2] = CFELookup<RealType>::_hexNbOffs[a][2];
    } else {       // Is a virtual node and thus an interpolation between two real nodes
      vertex[l][0] = 0.5 * ( CFELookup<RealType>::_hexNbOffs[a][0] + CFELookup<RealType>::_hexNbOffs[b][0] );
      vertex[l][1] = 0.5 * ( CFELookup<RealType>::_hexNbOffs[a][1] + CFELookup<RealType>::_hexNbOffs[b][1] );
      vertex[l][2] = 0.5 * ( CFELookup<RealType>::_hexNbOffs[a][2] + CFELookup<RealType>::_hexNbOffs[b][2] );
    }
  }

  edge[0] = vertex[1] - vertex[0];
  edge[1] = vertex[2] - vertex[0];
  edge[2] = vertex[3] - vertex[0];
}


void CFECube::makeOrientation() {
  for ( std::vector < CFETopoTetra >::iterator it = _tetraVec.begin(); it != _tetraVec.end(); ++it ) {
    aol::Vec3<double>   edge[3];
    getTetraEdges ( *it, edge );
    if ( determinant ( edge[0], edge[1], edge[2] ) < 0 ) {
      it->swap();
    }
  }
}


/** This is the method which does the splitting of the cubes into the set of tetrahedra
 *  based on the type of the cube (or element), i.e. the signs of the level-set function.
 *  Most of this method was directly adopted from the original code of F. Liehr.
 */
void CFECube::splitInTetra ( short par, int id ) {
  const short *knots;                      // handle for nodes of par-th standardtetra
  const short *x;                          // handle for edges of par-th standardtetra
  short a = 0, b = 0, c = 0, d = 0;        // dummies
  //  short negSub1, negSub2, posSub1, posSub2, neg1, neg2, pos1, pos2;                       // dummies for caching indices
  bool foundInterface = false;

  knots = CFETopoLookup::_stdTetraVertex[par];
  x = CFETopoLookup::_stdTetraEdges[par];

  d = knots[0] + knots[1] + knots[2] + knots[3];

  // if tetra t has at least one _interfaced edge:
  // t splits in tetra|prism  <=> there are two incident non_interfaced edges
  // t splits in tetra|tetra  <=> there are two nonincident non_interfaced edges

  short cnt = 0;                    // counts collected nodes of non_interfaced edges
  bool done = false;
  int i = 0;

  CFETopoTetra t;

  // Set the parent for all terahedra which are created in the following. This means that a regular tetrahedron (not split) is its own parent.
  t.setParent ( par );

  while ( !done ) {

    const short nnv = NODE_NOT_VIRTUAL; // abbreviation

    if ( !_interfaced[edgeBetween ( x[i], x[i + 1] ) ] ) {                  //  (x[i],x[i+1]) not interfaced?
      if ( cnt == 0 ) {                 //  first node to collect?
        t.setNode ( cnt++, a = x[i++] ); //  collect both  nodes
        t.setNode ( cnt++, b = x[i++] );
      } else {                   // if first node is already collected
        if ( t.getNode ( 0 ) == x[i] || t.getNode ( 1 ) == x[i] ) {  // insert x[i+1] and remaining node
          // (case: 2 incident, non_interfaced edges)
          t.setNode ( 2, c = x[i + 1] );
          t.setNode ( 3, d = d - a - b - c );

          while ( i < 12 ) {         // check trivial case (all edges non_interfaced)
            if ( _interfaced[edgeBetween ( x[i], x[i + 1] ) ] ) {
              foundInterface = true;
              i = 12;
            }
            i = i + 2;
          }

          if ( !foundInterface ) {
            //case: tetra does not split up
            if ( ( id & ( static_cast < short > ( ( 1 << a ) + ( 1 << b ) + ( 1 << c ) + ( 1 << d ) ) ) ) == 0 ) {
              // all nodes positive  => insert in positve list
              t.setVolType ( CFETopoTetra::P );
              t.setSign ( + 1 );
              _tetraVec.push_back ( t );
            } else {       // insert in negative list
              t.setVolType ( CFETopoTetra::N );
              t.setSign ( - 1 );
              _tetraVec.push_back ( t );
            }
            done = true;
          } else {         // <--- end <no splitting>
            // case: tetra-prism
            // now t consists of nodes (a,b,c,d) with a<b<c. a,b,c span non_interfaced edges of t
            if ( ( id & ( static_cast < short > ( 1 << d ) ) ) != 0 ) {       // is corner at d negative?
              // local indices reference levelset-values ( + + + - )
              _tetraVec.push_back ( CFETopoTetra ( a, nnv, b, nnv, d,   b, d, c, par, CFETopoTetra::PPPN , 2, + 1 ) ); //, 0, 4, 2));    // DO NOT CHANGE!!!
              _tetraVec.push_back ( CFETopoTetra ( b, nnv, c, nnv, a, nnv, d, c, par, CFETopoTetra::PPPN1, 1, + 1 ) ); //, 4, 5, 0));
              _tetraVec.push_back ( CFETopoTetra ( a, nnv, d,   a, d,   b, d, c, par, CFETopoTetra::PPPN2, 3, + 1 ) ); //, 0, 1, 2));
              _tetraVec.push_back ( CFETopoTetra ( d, nnv, d,   a, d,   b, d, c, par, CFETopoTetra::PPPN3, 3, - 1 ) ); //, 3, 1, 2));
              done = true;
            } else {
              // ( - - - + )
              _tetraVec.push_back ( CFETopoTetra ( a, nnv, b, nnv, b,   d, c, d, par, CFETopoTetra::NNNP , 2, - 1 ) ); //, 0, 4, 2));    // DO NOT CHANGE!!!
              _tetraVec.push_back ( CFETopoTetra ( b, nnv, c, nnv, a, nnv, c, d, par, CFETopoTetra::NNNP1, 1, - 1 ) ); //, 4, 5, 0));
              _tetraVec.push_back ( CFETopoTetra ( a, nnv, a,   d, b,   d, c, d, par, CFETopoTetra::NNNP2, 3, - 1 ) ); //, 0, 1, 2));
              _tetraVec.push_back ( CFETopoTetra ( d, nnv, a,   d, b,   d, c, d, par, CFETopoTetra::NNNP3, 3, + 1 ) ); //, 3, 1, 2));
              done = true;
            }
          }
        } else {
          // insert  x[i], x[i+1] (case: 2 nonincident, noninterfaced edges)
          t.setNode ( 2, c = x[i] );
          t.setNode ( 3, d = d - a - b - c );
          while ( i < 12 ) {         // checking trivial case (all edges noninterfaced)
            if ( _interfaced[edgeBetween ( x[i], x[i + 1] ) ] ) {
              foundInterface = true;
              i = 12;
            }
            i = i + 2;
          }
          if ( !foundInterface ) {
            //case: tetra does not split up
            if ( ( id & ( static_cast < short > ( ( 1 << a ) + ( 1 << b ) + ( 1 << c ) + ( 1 << d ) ) ) ) == 0 ) {
              t.setVolType ( CFETopoTetra::P );
              t.setSign ( + 1 );
              _tetraVec.push_back ( t );
            } else {
              t.setVolType ( CFETopoTetra::N );
              t.setSign ( - 1 );
              _tetraVec.push_back ( t );
            }
            done = true;
          } else {     // <--- end <noninterfaced>
            // case: prism-prism
            // now t contains nodes (a,b,c,d), with noninterfaced edges (a,b) and (c,d) only
            if ( ( id & ( static_cast < short > ( ( 1 << a ) + ( 1 << b ) ) ) ) != 0 ) {       //i.e. nodes a,b are negative and c,d positive
              // ( - - + + )
              _tetraVec.push_back ( CFETopoTetra ( a, nnv, b, nnv, b, d, b, c, par, CFETopoTetra::NNNPPP , 2, - 1 ) ); //, 0, 1, 2));    // DO NOT CHANGE
              _tetraVec.push_back ( CFETopoTetra ( a, nnv, a,   c, a, d, b, c, par, CFETopoTetra::NNNPPP1, 3, - 1 ) ); //, 0, 6, 7));    // IN ANY WAY!
              _tetraVec.push_back ( CFETopoTetra ( a, nnv, b,   d, a, d, b, c, par, CFETopoTetra::NNNPPP2, 3, - 1 ) ); //, 0, 2, 7));
              _tetraVec.push_back ( CFETopoTetra ( c, nnv, b,   d, a, d, b, c, par, CFETopoTetra::NNNPPP3, 3, + 1 ) ); //, 8, 2, 7));
              _tetraVec.push_back ( CFETopoTetra ( c, nnv, a,   c, a, d, b, c, par, CFETopoTetra::NNNPPP4, 3, + 1 ) ); //, 8, 6, 7));
              _tetraVec.push_back ( CFETopoTetra ( c, nnv, d, nnv, a, d, b, d, par, CFETopoTetra::NNNPPP5, 2, + 1 ) ); //, 3, 4, 5));
            } else {       // i.e. nodes a,b are positve and c,d negative
              // ( + + - - )
              _tetraVec.push_back ( CFETopoTetra ( a, nnv, b, nnv, d, b, c, b, par, CFETopoTetra::PPPNNN , 2, + 1 ) ); //, 0, 1, 2));
              _tetraVec.push_back ( CFETopoTetra ( a, nnv, c,   a, d, a, c, b, par, CFETopoTetra::PPPNNN1, 3, + 1 ) ); //, 0, 6, 7));
              _tetraVec.push_back ( CFETopoTetra ( a, nnv, d,   b, d, a, c, b, par, CFETopoTetra::PPPNNN2, 3, + 1 ) ); //, 0, 2, 7));
              _tetraVec.push_back ( CFETopoTetra ( c, nnv, d,   b, d, a, c, b, par, CFETopoTetra::PPPNNN3, 3, - 1 ) ); //, 8, 2, 7));
              _tetraVec.push_back ( CFETopoTetra ( c, nnv, c,   a, d, a, c, b, par, CFETopoTetra::PPPNNN4, 3, - 1 ) ); //, 8, 6, 7));
              _tetraVec.push_back ( CFETopoTetra ( c, nnv, d, nnv, d, a, d, b, par, CFETopoTetra::PPPNNN5, 2, - 1 ) ); //, 3, 4, 5));
            }
            done = true;
          }           // <--- end <prism-prism>
        }
      }
    } else {          // <--- end <if not _interfaced>
      i = i + 2;
      foundInterface = true;
    }
  }                   // <--- end <while>
}                     // <--- end <splitInTetra>

template class CFELookup < float >;
template class CFELookup < double >;
template class CFELookup < long double >;

}
