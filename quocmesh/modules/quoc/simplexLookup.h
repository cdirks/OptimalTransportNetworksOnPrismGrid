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

#ifndef __SIMPLEXLOOKUP_H
#define __SIMPLEXLOOKUP_H

#include <quoc.h>
#include <smallVec.h>

namespace qc {

namespace simplex {

// --------------------------------------------------------------------------

template <Dimension Dim>
struct TopologyLookup {};

template <typename RealType, Dimension Dim>
struct Lookup {};

// --------------------------------------------------------------------------

template <>
struct TopologyLookup<QC_2D> {
  //! localIndicesSimplexToCube[i][j] is the local
  //! index wrt square element of the j-th vertex
  //! in the i-th triangle.
  static const short localIndicesSimplexToCube[2][3];
  static const short edges[3][2];
  static const short edgesAsBits[3];

  static const short numSimplexesPerCube = 2;
  static const short numEdges = 3;
  static const short maxEdgeBit = 6;
};

// --------------------------------------------------------------------------

template <typename RealType>
struct Lookup<RealType, QC_2D> :
  public TopologyLookup<QC_2D> {

  typedef aol::Vec2<RealType> VecType;
  //! gradients[i][j] is the cartesian gradient
  //! of the j-th basis function on the i-th triangle.
  static const VecType gradients[2][3];

  static const VecType cartesianCubeNodeCoords[4];
};

// --------------------------------------------------------------------------

template <>
struct TopologyLookup<QC_3D> {
  //! localIndicesSimplexToCube[i][j] is the local
  //! index wrt cube element of the j-th vertex
  //! in the i-th tetrahedron.
  static const short localIndicesSimplexToCube[6][4];
  static const short edges[6][2];
  static const short edgesAsBits[6];

  static const short numSimplexesPerCube = 6;
  static const short numEdges = 6;
  static const short maxEdgeBit = 12;
};

// --------------------------------------------------------------------------

template <typename RealType>
struct Lookup<RealType, QC_3D> :
  public TopologyLookup<QC_3D> {

  typedef aol::Vec3<RealType> VecType;
  //! gradients[i][j] is the cartesian gradient
  //! of the j-th basis function on the i-th tetrahedron.
  static const VecType gradients[6][4];

  static const VecType cartesianCubeNodeCoords[8];
};

// --------------------------------------------------------------------------

} // end of namespace simplex.

} // end of namespace qc.

#endif
