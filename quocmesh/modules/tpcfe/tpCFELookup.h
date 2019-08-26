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

#ifndef __TPCFELOOKUP_H
#define __TPCFELOOKUP_H

#include <vec.h>
#include <tpCFEBasics.h>

namespace tpcfe {

class CFETopoLookup {
public:
  /** The vertices of the six standard tetrahedra
   */
  static const short _stdTetraVertex[6][4];

  /** The edges of the six standard tetrahedra
   */
  static const short _stdTetraEdges[6][12];

  /** The faces of a cube and the inscribed tetrahedra
   *  The first 5 entries describe vertices of the triangular faces (which lie diagonally inside the cube)
   *  the last 6 entries gives the standard quadratic outer faces of the cube.
   */
  static const short _elementFaces[11][4];

  /** Offsets to the neighboring points in regular tetrahedral grid
   */
  static const aol::Vec3<short> _tetNbOffsets[ NUM_NEIGHBORS ];

  // /** Offsets to neighboring points (in the sense of overlapping support of basis functions) in case of jumping coefficients/TPOS and similar
  //  */
  // static const aol::Vec3<short> _jcNbOffsets[ NUM_NEIGHBORS_JC ];

  /** Offsets to the neighboring points in regular hexahedral grid, sorted according to local index.
   */
  static const aol::Vec3<short> _hexNbOffsets[8];

  /** Stores the edges of the unit cube
   */
  static const short _cubeEdges[12][2];

  /** Stores the edges of the standard divison of a cube into
   *  tetrahedra.
   */
  static const short _tetraEdges[ NUM_TETRA_EDGES_IN_CUBE ][2];

  /** Once again the edges of the standard division of a cube into
   *  tetrahedra, for an edge from node i to j, both bits i and j are
   *  set. */
  static const short _edge[ NUM_TETRA_EDGES_IN_CUBE ];

  /** A signature is admissible if the corresponding sign pattern on a cube represents at most one cut of the interface through the cube.
   */
  static const bool _admissibleSignature[256];

  /** Contains a cube for each topology, i.e. for each configuration of pos and
   *  neg. level-set vertices
   */
  static const CFECube *_topo[256];
};

template <typename RealType>
class CFELookup : public CFETopoLookup {
private:
  static bool _initialized;
  static int refCounter;
  static int mem;

  /** Create the gradients of basis functions on the std tetrahedra
   */
  static void createBasisGradientsForStdTetra();

  /** Create the stiffness matrices for standard tetrahedra
   */
  static void createStiffnessForStdTetra();

  /** Compute the normals of the faces of a cube
   */
  static void computeNormals();

public:
  /** The volume of a standard tetrahedron = 1./6.;
   */
  static const RealType _stdTetraVolume;

  /** Maps the edges of the reference cube to coordinates in R^3
   */
  static const aol::Vec3<RealType> _hexNbOffs[8];

  /** Stores the normal of the element faces
   */
  static aol::Vec3<RealType> _faceNormal[11];

  /** Contains the gradients of the basis functions on the standard tetraheda.
   *  Used e.g. for the computations of the normals of interfaces at virtual nodes.
   */
  static aol::Vec3<RealType> _stdTetraBasisGradients[6][4];

  /** Stores the mass and stiffness matrix for the reference tetrahedron spanned by (0,0,0), (0,0,1), (0,1,0), (1,0,0).
   *  Note that the mass matrices are the same for all six standard tetrahedra, the stiffness matrices are not.
   */
  static const RealType _localMassMatrixRefTetra[4][4];
  static const RealType _localStiffMatrixRefTetra[4][4];

  /** Stores the stiffness matrix for the standard tetrahedra. The standard tetrahedra
   *  are the ones that combine to the unit cube. I.e. they have volume 1./6.
   */
  static RealType _localStiffnessMatrixStdTetra[6][4][4];

  // see cpp file
  static const RealType _gradientStdTetraBasisFct[4][3];

  static void writeRefTetraVecTo ( std::vector< tpcfe::CFETetra<RealType> > &refTetraVec );

public:
  CFELookup () {}

  ~CFELookup( ) {}

  static bool isInitialized () {
    return _initialized;
  }

  /** Build up all the tables and increase the ref count by one
   */
  static int construct();

  /** Decrease the ref counter and clean up all used memory is ref count == 0
   */
  static void destruct();
};



/** Class for storage of topological data in lookup table.
 */
class CFECube {
protected:
  /** Tells whether an edge between two nodes represented by the bits
   *  of the above variable @see edge is interfaced
   */
  std::vector < bool > _interfaced;
  std::vector < CFETopoTetra > _tetraVec;

public:
  CFECube ( const int id );

  int memoryUsed() const;

  /** Return a vector of tetrahedra this cube splits into.
   *  sign determines which tetrahedra are iterated: 0=all, -1=only (-) tetra,
   *  +1= only (+) tetra
   */
  const vector < CFETopoTetra > &getTopoTetraVec ( ) const ;


protected:
  /** Split a standard tetrahedron into its sub-tetrahedra
   */
  void splitInTetra ( short par, int id );

  /** Manipulates index_Permutaion and tetra_generating_Vector's indices for all tetrahedra that turn
   *  out to be a leftSystem if one proceeds as in CFEElement.
   */
  void makeOrientation();

  //! auxiliary-function of makeOrientation() ( usage independend of dofs -> no excessive use )
  template< typename RealType >
  inline RealType determinant ( aol::Vec3 < RealType >&a, aol::Vec3 < RealType >&b, aol::Vec3 < RealType >&c ) const {
    return ( a.x() * b.y() * c.z() + a.z() * b.x() * c.y() +
             a.y() * b.z() * c.x() - a.x() * b.z() * c.y() -
             a.y() * b.x() * c.z() - a.z() * b.y() * c.x() );
  }

  inline int edgeBetween ( short x, short y ) { return static_cast < int > ( ( 1 << x ) + ( 1 << y ) ); }

  template< typename RealType >
  void getTetraEdges ( const CFETopoTetra &t, aol::Vec3<RealType> vertex[3] );

};




}

#endif
