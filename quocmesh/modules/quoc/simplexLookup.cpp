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

#include <simplexLookup.h>

namespace qc {

namespace simplex {

// ==========================================================================

const short TopologyLookup<QC_2D>::
localIndicesSimplexToCube[2][3] = { { 0, 1, 2 },
                                    { 1, 2, 3 } };

// --------------------------------------------------------------------------

const short TopologyLookup<QC_2D>::
edges[3][2] = { { 0, 1 },
                { 0, 2 },
                { 1, 2 } };

// --------------------------------------------------------------------------

// each entry has 2 bits set, according to the indices
// of the edge.
const short TopologyLookup<QC_2D>::       // edges:
edgesAsBits[3] = { (1 << 0) + (1 << 1),   // 0 -> 1
                   (1 << 0) + (1 << 2),   // 0 -> 2
                   (1 << 1) + (1 << 2) }; // 1 -> 2

// --------------------------------------------------------------------------

template <typename RealType>
const typename Lookup<RealType, QC_2D>::VecType Lookup<RealType, QC_2D>::
gradients[2][3] = { { aol::Vec2<RealType> ( -1, -1 ),
                      aol::Vec2<RealType> (  1,  0 ),
                      aol::Vec2<RealType> (  0,  1 ) },
                    { aol::Vec2<RealType> (  0, -1 ),
                      aol::Vec2<RealType> ( -1,  0 ),
                      aol::Vec2<RealType> (  1,  1 ) } };

// --------------------------------------------------------------------------

template <typename RealType>
const typename Lookup<RealType, QC_2D>::VecType Lookup<RealType, QC_2D>::
cartesianCubeNodeCoords[4] = { aol::Vec2<RealType> ( 0, 0 ),   // 0
                               aol::Vec2<RealType> ( 1, 0 ),   // 1
                               aol::Vec2<RealType> ( 0, 1 ),   // 2
                               aol::Vec2<RealType> ( 1, 1 ) }; // 3

// ==========================================================================

// all tetra have positive orientation.
const short TopologyLookup<QC_3D>::
localIndicesSimplexToCube[6][4] = { { 0, 1, 2, 6 },
                                    { 1, 0, 5, 6 },
                                    { 0, 4, 5, 6 },
                                    { 2, 1, 3, 7 },
                                    { 1, 2, 6, 7 },
                                    { 5, 1, 6, 7 } };

// --------------------------------------------------------------------------

const short TopologyLookup<QC_3D>::
edges[6][2] = { { 0, 1 },
                { 0, 2 },
                { 0, 3 },
                { 1, 2 },
                { 1, 3 },
                { 2, 3 } };

// --------------------------------------------------------------------------

// each entry has 2 bits set, according to the indices
// of the edge.
const short TopologyLookup<QC_3D>::       // edges:
edgesAsBits[6] = { (1 << 0) + (1 << 1),   // 0 -> 1
                   (1 << 0) + (1 << 2),   // 0 -> 2
                   (1 << 0) + (1 << 3),   // 0 -> 3
                   (1 << 1) + (1 << 2),   // 1 -> 2
                   (1 << 1) + (1 << 3),   // 1 -> 3
                   (1 << 2) + (1 << 3) }; // 2 -> 3

// --------------------------------------------------------------------------

template <typename RealType>
const typename Lookup<RealType, QC_3D>::VecType Lookup<RealType, QC_3D>::
gradients[6][4] = { { aol::Vec3<RealType> ( -1, -1,  0 ),
                      aol::Vec3<RealType> (  1,  0,  0 ),
                      aol::Vec3<RealType> (  0,  1, -1 ),
                      aol::Vec3<RealType> (  0,  0,  1 ) },
           /* 1 */  { aol::Vec3<RealType> (  1,  1, -1 ),
                      aol::Vec3<RealType> ( -1, -1,  0 ),
                      aol::Vec3<RealType> (  0, -1,  1 ),
                      aol::Vec3<RealType> (  0,  1,  0 ) },
           /* 2 */  { aol::Vec3<RealType> (  0,  0, -1 ),
                      aol::Vec3<RealType> ( -1, -1,  1 ),
                      aol::Vec3<RealType> (  1,  0,  0 ),
                      aol::Vec3<RealType> (  0,  1,  0 ) },
           /* 3 */  { aol::Vec3<RealType> ( -1,  0,  0 ),
                      aol::Vec3<RealType> (  0, -1,  0 ),
                      aol::Vec3<RealType> (  1,  1, -1 ),
                      aol::Vec3<RealType> (  0,  0,  1 ) },
           /* 4 */  { aol::Vec3<RealType> (  0, -1,  0 ),
                      aol::Vec3<RealType> (  0,  1, -1 ),
                      aol::Vec3<RealType> ( -1, -1,  1 ),
                      aol::Vec3<RealType> (  1,  1,  0 ) },
           /* 5 */  { aol::Vec3<RealType> (  0, -1,  1 ),
                      aol::Vec3<RealType> (  0,  0, -1 ),
                      aol::Vec3<RealType> ( -1,  0,  0 ),
                      aol::Vec3<RealType> (  1,  1,  0 ) } };

// --------------------------------------------------------------------------

template <typename RealType>
const typename Lookup<RealType, QC_3D>::VecType Lookup<RealType, QC_3D>::
cartesianCubeNodeCoords[8] = { aol::Vec3<RealType> ( 0, 0, 0 ),   // 0
                               aol::Vec3<RealType> ( 1, 0, 0 ),   // 1
                               aol::Vec3<RealType> ( 0, 1, 0 ),   // 2
                               aol::Vec3<RealType> ( 1, 1, 0 ),   // 3
                               aol::Vec3<RealType> ( 0, 0, 1 ),   // 4
                               aol::Vec3<RealType> ( 1, 0, 1 ),   // 5
                               aol::Vec3<RealType> ( 0, 1, 1 ),   // 6
                               aol::Vec3<RealType> ( 1, 1, 1 ) }; // 7

// ==========================================================================

template struct Lookup<float, QC_2D>;
template struct Lookup<double, QC_2D>;
template struct Lookup<long double, QC_2D>;

template struct Lookup<float, QC_3D>;
template struct Lookup<double, QC_3D>;
template struct Lookup<long double, QC_3D>;

} // end of namespace simplex.

} // end of namespace qc.
