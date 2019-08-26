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

#include <tpCFEGrid.h>
#include <tpCFEUtils.h>
#include <tpCFETetrahedron.h>

typedef double RealType;
typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD > GridType;

// const char* FNMask = "out/rre_3c/Def%d/Diam_0.400_0.400_0.400_Rmv_0.000_0.000_0.000_Rnd_0_fix_2_shift_2%s";

static const short n_rods = 3;
static const int nPts = 1000;

int main ( int, char** argv ) {
  try {

    const int compLevel = atoi ( argv[1] ), gstdLevel = atoi ( argv[2] );
    const char* FNMask = argv[3];

    GridType gstdGrid ( gstdLevel ), compGrid ( compLevel );
    qc::ScalarArray<RealType, qc::QC_3D> gstdLevelset ( gstdGrid ), compLevelset ( compGrid );

    char filenamemask[1024];

    sprintf ( filenamemask, FNMask, gstdLevel, ".dat.bz2" );
    gstdLevelset.load ( filenamemask );

    sprintf ( filenamemask, FNMask, compLevel, ".dat.bz2" );
    compLevelset.load ( filenamemask );

    gstdGrid.setDomainFrom ( gstdLevelset );
    gstdGrid.detectAndInitVirtualNodes();

    compGrid.setDomainFrom ( compLevelset );
    compGrid.detectAndInitVirtualNodes();

    qc::MultiArray<RealType, qc::QC_3D> gstdDeformations ( gstdGrid ), compDeformations ( compGrid );

    sprintf ( filenamemask, FNMask, gstdLevel, "_Def%d.dat.bz2" );
    gstdDeformations.load ( filenamemask );

    sprintf ( filenamemask, FNMask, compLevel, "_Def%d.dat.bz2" );
    compDeformations.load ( filenamemask );

#if 1
    RealType L1sum = 0, L2sum = 0, Linf = 0, volSum = 0;

    const qc::GridSize<qc::QC_3D> gridSize ( compGrid );
    for ( GridType::FullElementIterator it ( compGrid ); it.notAtEnd(); ++it ) {
      tpcfe::CFEElement<RealType> el ( *it, gridSize, compGrid.getElType ( *it ) );

      el.computeAssembleData ( compGrid );

      for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el, -1 ); tit.notAtEnd(); ++tit ) {
        const tpcfe::CFETetra<RealType> &t = *tit;

        aol::Vec3<RealType> barycenter, gstdValue, compValue;
        t.computeGlobalCoordinateOfBarycenter ( barycenter, el );
        barycenter *= compGrid.H(); // global -> world coordinates; barycenter of current tetrahedron

        bool foundInside = false;
        {
          tpcfe::CFEElement<RealType> elem;
          tpcfe::CFETetra<RealType> tetra;
          aol::Vec<4,RealType> baryCO;
          tpcfe::determineElementAndTetraForWC<GridType> ( gstdGrid, barycenter, elem, tetra, baryCO, foundInside );
        }

        if ( foundInside ) {
          for ( short d = 0; d < 3; ++d ) {
            gstdValue[d] = tpcfe::interpolateDataAtPositionWC ( gstdGrid, gstdDeformations[d], barycenter );
            compValue[d] = tpcfe::interpolateDataAtPositionWC ( compGrid, compDeformations[d], barycenter );
          }

          const RealType
            vol = t.getVolume() * aol::Cub ( compGrid.H() ), // volume of current tetrahedron
            diffNorm = ( gstdValue - compValue ).norm();

          L1sum += vol * ( diffNorm );
          L2sum += vol * aol::Sqr ( diffNorm );
          Linf = aol::Max ( Linf, diffNorm );
          volSum += vol;
        }

      }
    }

    const RealType
      BCDeformation = 0.01,
      L1NormScaled = L1sum / BCDeformation / volSum,
      L2NormScaled = sqrt ( L2sum / volSum) / BCDeformation,
      LinfScaled = Linf / BCDeformation;

    cerr << "Vol: " << aol::detailedFormat ( volSum ) << " abs_avg: " << aol::detailedFormat ( L1NormScaled ) << " rms_avg: " << aol::detailedFormat ( L2NormScaled ) << " max: " << aol::detailedFormat ( LinfScaled ) << endl;
#else
    aol::MultiVector<RealType> differences ( 3, 3 * nPts * aol::Sqr ( n_rods ) );

    int mvi = 0;

    for ( int ir = 0; ir < n_rods; ++ir ) {
      for ( int jr = 0; jr < n_rods; ++jr ) {
        const RealType A = 0.0 + 1.0 * ( 2.0 * ir + 1.0 ) / ( 2.0 * n_rods );
        const RealType B = 0.0 + 1.0 * ( 2.0 * jr + 1.0 ) / ( 2.0 * n_rods );
        for ( int pr = 0; pr < nPts; ++pr ) {
          const aol::Vec3<RealType> xRpos ( pr / (  nPts - 1.0 ), A, B );
          const aol::Vec3<RealType> yRpos ( A, pr / (  nPts - 1.0 ), B );
          const aol::Vec3<RealType> zRpos ( A, B, pr / (  nPts - 1.0 ) );

          for ( short d = 0; d < 3; ++d )
            differences[d][mvi] = ( tpcfe::interpolateDataAtPositionWC ( gstdGrid, gstdDeformations[d], xRpos ) - tpcfe::interpolateDataAtPositionWC ( compGrid, compDeformations[d], xRpos ) );
          ++mvi;

          for ( short d = 0; d < 3; ++d )
            differences[d][mvi] = ( tpcfe::interpolateDataAtPositionWC ( gstdGrid, gstdDeformations[d], yRpos ) - tpcfe::interpolateDataAtPositionWC ( compGrid, compDeformations[d], yRpos ) );
          ++mvi;

          for ( short d = 0; d < 3; ++d )
            differences[d][mvi] = ( tpcfe::interpolateDataAtPositionWC ( gstdGrid, gstdDeformations[d], zRpos ) - tpcfe::interpolateDataAtPositionWC ( compGrid, compDeformations[d], zRpos ) );
          ++mvi;

        }
      }
    }

    const RealType L1factor = ( 1.0 / 0.01 ) *      ( 1.0 / nPts ) / ( 3 * aol::Sqr ( n_rods ) ); // relative to prescribed deformation; averaged over nPts along axes of length 1
    const RealType L2factor = ( 1.0 / 0.01 ) * sqrt ( 1.0 / nPts ) / ( 3 * aol::Sqr ( n_rods ) ); // same; one square root of nPts already hidden in l2-norm
    const RealType Lifactor = ( 1.0 / 0.01 ); // only relative to prescribed deformation;

    cerr << "  abs x " << aol::detailedFormat ( differences[0].lpNorm ( 1.0 ) * L1factor )
         << "  rms x " << aol::detailedFormat ( differences[0].norm()  * L2factor )
         << "  max x " << aol::detailedFormat ( differences[0].getMaxAbsValue() * Lifactor ) << endl
         << "  abs y " << aol::detailedFormat ( differences[1].lpNorm ( 1.0 ) * L1factor )
         << "  rms y " << aol::detailedFormat ( differences[1].norm() * L2factor )
         << "  max y " << aol::detailedFormat ( differences[1].getMaxAbsValue() * Lifactor ) << endl
         << "  abs z " << aol::detailedFormat ( differences[2].lpNorm ( 1.0 ) * L1factor )
         << "  rms z " << aol::detailedFormat ( differences[2].norm() * L2factor )
         << "  max z " << aol::detailedFormat ( differences[2].getMaxAbsValue() * Lifactor ) << endl;
#endif

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
