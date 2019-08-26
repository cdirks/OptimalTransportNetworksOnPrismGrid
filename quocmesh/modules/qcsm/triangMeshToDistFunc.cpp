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

#include <triangMeshToDistFunc.h>



namespace qcsm {

template<class RealType>
void TriangMeshToDistFunc<RealType>::makeDistFkt ( ) {
  _dist.setAll ( aol::NumberTrait<RealType>::Inf );
  const int smnots = _mesh.getNumTriangs();
  RealType min_max[6];//={0};

  for ( int i = 0; i < smnots; ++i ) {
    memset ( min_max, 0, sizeof ( min_max ) );//set to 0
    min_max[3] = min_max[4] = min_max[5] = 1;
    for ( int j = 0; j < 3; ++j ) {
      const RealType x = _mesh.getVertexCoord ( _mesh.getTriangNodeIdx ( i, j ), 0 );
      const RealType y = _mesh.getVertexCoord ( _mesh.getTriangNodeIdx ( i, j ), 1 );
      const RealType z = _mesh.getVertexCoord ( _mesh.getTriangNodeIdx ( i, j ), 2 );

      if ( x > min_max[0] ) min_max[0] = x;
      if ( y > min_max[1] ) min_max[1] = y;
      if ( z > min_max[2] ) min_max[2] = z;

      if ( x < min_max[3] ) min_max[3] = x;
      if ( y < min_max[4] ) min_max[4] = y;
      if ( z < min_max[5] ) min_max[5] = z;
    }

    distFkt ( min_max, i );

  }

}

template<class RealType>
void TriangMeshToDistFunc<RealType>::distFkt ( RealType min_max[6], int SmCellId ) {

  const RealType numericAccuracy = 1.e-16;

  const int gridWidthM1 = static_cast<int> ( 1. / _h ) - 1;

  const int startX = aol::Max ( 0, static_cast<int> ( min_max[3] / _h ) );
  const int startY = aol::Max ( 0, static_cast<int> ( min_max[4] / _h ) );
  const int startZ = aol::Max ( 0, static_cast<int> ( min_max[5] / _h ) );

  const int endX = aol::Min ( gridWidthM1, static_cast<int> ( min_max[0] / _h + 1 ) );
  const int endY = aol::Min ( gridWidthM1, static_cast<int> ( min_max[1] / _h + 1 ) );
  const int endZ = aol::Min ( gridWidthM1, static_cast<int> ( min_max[2] / _h + 1 ) );

  if ( endX < 0 || endY < 0 || endZ < 0 || startX > gridWidthM1 || startY > gridWidthM1 || startZ > gridWidthM1 )
    return;

  const aol::Triangle<RealType> triangle ( _mesh, SmCellId );
  aol::Vec3<RealType> normal( triangle.weightedNormal() );
  const RealType area = normal.norm();

  if( area == 0 )
    return;
  normal /= area;

  for ( int z = startZ ; z <= endZ ; ++z ) {
    for ( int y = startY ; y <= endY ; ++y ) {
      for ( int x = startX ; x <= endX ; ++x ) {

        const aol::Vec3<RealType> qcp ( x * _h, y * _h, z * _h );

        RealType signedDistToPlane = 0.;
        const RealType euclidianDistToPlane = triangle.calcDist ( qcp, normal, signedDistToPlane );
        if ( !_signedDist )
          signedDistToPlane = 1.;
        if ( aol::Abs ( aol::Abs ( _dist.get ( x, y, z ) ) - euclidianDistToPlane ) < numericAccuracy ) { // another triangle has same distance
          if ( _distanceToPlaneThroughClosestTriangle.get ( x, y, z ) < aol::Abs ( signedDistToPlane ) ) { // the new triangle distance sign is more reliable
            _dist.set ( x, y, z, euclidianDistToPlane * aol::signum ( signedDistToPlane ) );
            _distanceToPlaneThroughClosestTriangle.set ( x, y, z, aol::Abs ( signedDistToPlane ) );
          }
        } else if ( aol::Abs ( _dist.get ( x, y, z ) ) > euclidianDistToPlane ) { // the new triangle is closest
          _dist.set ( x, y, z, euclidianDistToPlane * aol::signum ( signedDistToPlane ) );
          _distanceToPlaneThroughClosestTriangle.set ( x, y, z, aol::Abs ( signedDistToPlane ) );
          qc::CoordType seedPoint ( x, y, z );
          _seedpoints.push_back ( seedPoint );
        }
      }
    }
  }
}

template class  TriangMeshToDistFunc<float>;
template class  TriangMeshToDistFunc<double>;

}
