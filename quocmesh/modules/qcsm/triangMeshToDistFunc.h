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

#ifndef __TRIANGMESHTODISTFUNC_H
#define __TRIANGMESHTODISTFUNC_H


#include <scalarArray.h>
#include <triangMesh.h>
#include <triangle.h>

namespace qcsm {

/**
 * \brief Helper class to calculate the distance function of an aol::TriangMesh.
 *
 * \note Only calculates the distance on the grid nodes close to the aol::TriangMesh, but also
 *       supplies seed points that can be used by eik::EikonalNA3D to calculate the distance
 *       function on the whole quoc grid.
 *
 * \author Notthoff
 */
template<class RealType>
class TriangMeshToDistFunc {

public:
  TriangMeshToDistFunc ( const aol::TriangMesh<RealType> &Mesh,
                         qc::ScalarArray<RealType, qc::QC_3D> &Dist,
                         vector<qc::CoordType> &Seedpoints,
                         const bool SignedDist = false )
    : _mesh ( Mesh ),
      _dist ( Dist ),
      _distanceToPlaneThroughClosestTriangle ( Dist, aol::STRUCT_COPY ),
      _seedpoints ( Seedpoints ),
      _signedDist ( SignedDist ),
      _h ( 1. / ( Dist.getNumXYZ() - 1. ) ) {}

  /** Calculate the distance for each triangle of the SurfMesh
   */
  void makeDistFkt ( );

private:
  void distFkt ( RealType min_max[6], int SmCellId );

  const aol::TriangMesh<RealType> &_mesh;
  qc::ScalarArray<RealType, qc::QC_3D> &_dist;
  qc::ScalarArray<RealType, qc::QC_3D> _distanceToPlaneThroughClosestTriangle;
  vector<qc::CoordType> &_seedpoints;
  const bool _signedDist;
  //! edgelength of one cube
  const RealType _h;
};

}
#endif // __TRIANGMESHTODISTFUNC_H
