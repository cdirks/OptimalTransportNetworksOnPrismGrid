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

#include <signedDistanceSweeping.h>

#include <bitArray.h>
#include <scalarArray.h>
#include <sweeping.h>
#include <tpCFEGrid.h>
#include <tpCFEUtils.h>


namespace qc {

template< typename RealType >
void SignedDistanceSweepingOp3D<RealType>::applyAdd ( const qc::ScalarArray<RealType, qc::QC_3D> &Arg, qc::ScalarArray<RealType, qc::QC_3D> &Dest, const qc::BitArray<qc::QC_3D>* const pIsRelevant ) const {
  qc::ScalarArray<RealType, qc::QC_3D> tmpArray ( Dest, aol::STRUCT_COPY );

  qc::RectangularGrid< qc::QC_3D > grid ( qc::GridSize<qc::QC_3D>::createFrom ( Arg ) );
  tmpArray.setAll ( aol::NumberTrait<double>::Inf );
  qc::BitArray<qc::QC_3D> seedPoints ( grid );

  tpcfe::CFEGrid< RealType, tpcfe::CFE_NONE > cfeGrid ( qc::GridSize<qc::QC_3D>::createFrom ( Arg ) );
  cfeGrid.setAdjustLevelset ( 1e-10 );

  cfeGrid.setDomainFrom ( Arg );
  const RealType h = grid.H();

  for ( tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > itit( cfeGrid ); itit.notAtEnd(); ++itit ) {
    std::vector< aol::Vec3< RealType > > interfaceVertex = itit.getTriangRef().getGlobalCoordinates();
    for ( int i = 0; i < 3; ++i ) {
      interfaceVertex[i] *= h;
    }
    const aol::Triangle<RealType> triangle ( interfaceVertex[0], interfaceVertex[1], interfaceVertex[2] );
    // interfaceVertex contains world coordinates of the vertices of the current interface triangle

    const qc::CoordType elCoord ( itit.getTriangRef().getElement().x(), itit.getTriangRef().getElement().y(), itit.getTriangRef().getElement().z() );
    for ( short i = 0; i < 2 ; ++i ) {
      for ( short j = 0 ; j < 2 ; ++j ) {
        for ( short k = 0 ; k < 2 ; ++k ) {
          const qc::CoordType qcp = elCoord + qc::CoordType ( k, j, i );           // one cube vertex in global coordinates
          const aol::Vec3<RealType> coord ( qcp[0] * h, qcp[1] * h, qcp[2] * h );  // the same cube vertex in world coordinates

          const RealType tmpdist = triangle.calcDist ( coord );
          if ( aol::isNaN ( tmpdist ) )
            throw aol::Exception( "tmpdist is NaN", __FILE__, __LINE__);

          // not trusting the orientation of the triangle, the sign needs to be determined from the levelset function
          if ( fabs ( tmpArray.get ( qcp ) ) > fabs ( tmpdist ) ) {
            tmpArray.set ( qcp, aol::signum ( Arg.get ( qcp ) ) * tmpdist );
          }
          seedPoints.set ( qcp, true );
        } // k
      } // j
    } // i
  }

  qc::DistanceSweeper3d< RealType, qc::BitArray<qc::QC_3D>, qc::ScalarArray< RealType, qc::QC_3D >, qc::RectangularGrid< qc::QC_3D > > distanceSweeper ( seedPoints, pIsRelevant, grid );
  distanceSweeper.setVerboseMode ( true );

  distanceSweeper.computeDistances ( tmpArray, _delta );

  Dest += tmpArray;
}


template class SignedDistanceSweepingOp3D<float>;
template class SignedDistanceSweepingOp3D<double>;
template class SignedDistanceSweepingOp3D<long double>;

}
