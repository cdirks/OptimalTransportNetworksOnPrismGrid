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

#ifndef __BONEELAST_THRESHOLDING_H
#define __BONEELAST_THRESHOLDING_H

#include <progressBar.h>
#include <scalarArray.h>
#include <tpCFELevelsets.h>

/** Determine segmentation threshold similar to Ridler, Calvard:
 *  Picture Thresholding Using an Interactive Selection Method, IEEE
 *  Transactions on Systems, Man, and Cybernetics 1978, SMC-8,
 *  pp. 630-632.
 *  \author Schwen
 */
template< typename RealType >
RealType RiCaThreshold ( const qc::ScalarArray< RealType, qc::QC_3D> &dataset, const RealType HistoAnalyzeMin, const RealType HistoAnalyzeMax, const RealType volScaleFactor, const char* histoFNRoot = NULL ) {

  qc::ScalarArray< RealType, qc::QC_3D> datasetCopy ( dataset, aol::DEEP_COPY );

  const RealType
    nX = dataset.getNumX(),
    nY = dataset.getNumY(),
    nZ = dataset.getNumZ(),
    nM = aol::Max ( nX, nY, nZ ),
    nL = aol::Max ( nX - 1, nY - 1, nZ - 1 ),
    scaleFac = ( (nX-1) * (nY-1) * (nZ-1) / aol::Cub(nL) ) / ( nX * nY * nZ / aol::Cub ( nM )  ),
    voxelVol = scaleFac / aol::Cub ( nM ),
    GVOffsetB = - dataset.getMinValue();

  datasetCopy.addToAll ( GVOffsetB );
  cerr << "GVOffsetB = " << aol::detailedFormat ( GVOffsetB ) << endl;

  const int histoSteps = static_cast< int > ( datasetCopy.getMaxValue() ) + 1;
  aol::Vector<int> histogram ( histoSteps );
  for ( int i = 0; i < datasetCopy.size(); ++i ) {
    if ( datasetCopy[i] - GVOffsetB > HistoAnalyzeMin && datasetCopy[i] - GVOffsetB < HistoAnalyzeMax ) {
      ++( histogram[ static_cast<int>( datasetCopy[i]  ) ] );
    }
  }

  if ( histoFNRoot != NULL ) {
    char histoFN[1024];
    sprintf ( histoFN, "%s_HistoVoxelNoisy.dat", histoFNRoot );
    ofstream histo ( histoFN );
    for ( int i = 0; i < histogram.size(); ++i ) {
      histo << aol::detailedFormat ( i - GVOffsetB ) << " " << aol::detailedFormat ( voxelVol * histogram.size() * histogram[i] ) << endl;
    }

    sprintf ( histoFN, "%s_HistoVoxelNoisyIntegrated.dat", histoFNRoot );
    ofstream histoIntegr ( histoFN );
    RealType histoSum = 0.0;
    for ( int i = histogram.size() - 1 ; i >= 0; --i ) {
      histoSum += voxelVol * histogram[i];
      histoIntegr << aol::detailedFormat ( i - GVOffsetB ) << " " << aol::detailedFormat ( histoSum ) << endl;
    }
  }

  RealType estimatedThreshold = datasetCopy.sum() / datasetCopy.size();
  for ( int iter = 0; iter < 20; ++iter ) {
    RealType lowerNum = 0, lowerDen = 0, upperNum = 0, upperDen = 0;
    for ( int i = 0; i <= static_cast<int>( estimatedThreshold ); ++i ) {
      lowerNum += i * histogram[i];
      lowerDen += 2 * histogram[i];
    }
    for ( int i = static_cast<int>( estimatedThreshold ) + 1; i < histoSteps ; ++i ) {
      upperNum += i * histogram[i];
      upperDen += 2 * histogram[i];
    }

    estimatedThreshold = lowerNum / lowerDen + upperNum / upperDen;

    cerr << aol::detailedFormat ( estimatedThreshold ) << endl;
  }

  cerr << "Estimated threshold (voxel-based) is " << aol::detailedFormat ( estimatedThreshold - GVOffsetB );

  RealType histoSum = 0.;
  for ( int i = histogram.size() - 1 ; i >= estimatedThreshold; --i ) {
    histoSum += voxelVol * histogram[i];
  }
  cerr << " with volume " << aol::detailedFormat ( volScaleFactor * histoSum ) << endl;

  return ( estimatedThreshold - GVOffsetB );

}


// resampleFactor is changed to the one actually used
template< typename RealType >
void resampleDatasetByFactor ( const qc::ScalarArray< RealType, qc::QC_3D> &datasetOrig, const RealType outResolution, RealType &resampleFactor, RealType &volScaleFactor, qc::ScalarArray< RealType, qc::QC_3D> &datasetSmall ) {
  const short
    roiX = datasetOrig.getNumX(),
    roiY = datasetOrig.getNumY(),
    roiZ = datasetOrig.getNumZ(),
    resX = static_cast<short> ( resampleFactor * roiX ),
    resY = static_cast<short> ( resampleFactor * roiY ),
    resZ = static_cast<short> ( resampleFactor * roiZ );

  volScaleFactor = aol::Cub ( outResolution * ( resZ - 1 ) );

  resampleFactor = aol::Min ( ( 1.0 * resX ) / roiX, ( 1.0 * resY ) / roiY, ( 1.0 * resZ ) / roiZ );
  datasetSmall.reallocate ( static_cast<short> ( roiX * resampleFactor ), static_cast<short> ( roiY * resampleFactor ), static_cast<short> ( roiZ * resampleFactor ) );

  cerr << "resampling by factor " << resampleFactor << " to size " << datasetSmall.getSize() << endl;

  for ( qc::RectangularIterator<qc::QC_3D> rit ( datasetSmall ); rit.notAtEnd(); ++rit ) {
    datasetSmall.set ( *rit, datasetOrig.interpolate ( aol::Vec3<RealType>( (*rit) ) / resampleFactor ) );
  }

}

// resampleFactor is changed to the one actually used
template< typename RealType >
void resampleDatasetByIntegerFactor ( const qc::ScalarArray< RealType, qc::QC_3D> &datasetOrig, qc::ScalarArray< RealType, qc::QC_3D> &datasetSmall, const short resampleIntFactor ) {
  const short
    roiX = datasetOrig.getNumX(),
    roiY = datasetOrig.getNumY(),
    roiZ = datasetOrig.getNumZ(),
    resX = roiX / resampleIntFactor, // integer division used on purpose
    resY = roiY / resampleIntFactor,
    resZ = roiZ / resampleIntFactor;

  datasetSmall.reallocate ( resX, resY, resZ );

  for ( qc::RectangularIterator<qc::QC_3D> rit ( datasetSmall ); rit.notAtEnd(); ++rit ) {
    RealType averageVal = 0;
    for ( qc::RectangularIterator<qc::QC_3D> bit ( qc::CoordType( 0, 0, 0 ), qc::CoordType( resampleIntFactor, resampleIntFactor, resampleIntFactor ) ); bit.notAtEnd(); ++bit ) {
      averageVal += datasetOrig.get ( resampleIntFactor * (*rit) + (*bit) );
    }
    datasetSmall.set ( *rit, averageVal / aol::Cub ( resampleIntFactor ) );
  }
}


template< typename RealType >
void thresholdAndDepoempelize ( const qc::ScalarArray< RealType, qc::QC_3D> &dataset, const RealType threshold, const RealType volScaleFactor, qc::ScalarArray< RealType, qc::QC_3D > &dest, const bool removePoempels = true ) {

  dest = dataset;
  dest.addToAll ( -threshold );
  dest *= -1.0;

  if ( removePoempels ) {
    tpcfe::TopBottomPoempelRemover<RealType> tbpr;
    tbpr.applySingle ( dest );
    tpcfe::SmallComponentsPoempelRemover<RealType> scpr;
    scpr.applySingle ( dest );
  }

  {
    tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > gridTmp ( qc::GridSize<qc::QC_3D>::createFrom ( dest ) );

    int innerVoxels = 0;
    for ( int i = 0; i < dest.size(); ++i ) {
      if ( dest[i] < 0 ) {
        ++innerVoxels;
      }
    }
    RealType totVol = aol::Cub ( 1.0 / aol::Max ( dest.getNumX(), dest.getNumY(), dest.getNumZ() ) ) * innerVoxels; // different interpretation of voxels!
    cerr << "Total voxel-based volume = " << aol::detailedFormat ( totVol )
         << " for voxel threshold " << aol::detailedFormat ( threshold )
         << " corresponding to " << aol::detailedFormat ( volScaleFactor * totVol ) << "m^3"<< endl;

    gridTmp.setDomainFrom ( dest );
    gridTmp.detectAndInitVirtualNodes();
    totVol = gridTmp.getTotalInnerVolume();
    cerr << "Total segmented volume = " << aol::detailedFormat ( totVol )
         << " for voxel threshold " << aol::detailedFormat ( threshold )
         << " corresponding to " << aol::detailedFormat ( volScaleFactor * totVol ) << " m^3" << endl;
  }
}


template< typename RealType >
RealType determineCFEVolumeSensitivity ( const qc::ScalarArray< RealType, qc::QC_3D> &dataset, const RealType thresholdDelta, const RealType volScaleFactor ) {
  RealType volMinus = -1, volPlus = -1;
  {
    qc::ScalarArray< RealType, qc::QC_3D> datasetTemp ( dataset, aol::DEEP_COPY );
    datasetTemp.addToAll ( -thresholdDelta );
    tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > gridTmp ( qc::GridSize<qc::QC_3D>::createFrom ( datasetTemp ) );
    gridTmp.setDomainFrom ( datasetTemp );
    gridTmp.detectAndInitVirtualNodes();
    volMinus = gridTmp.getTotalInnerVolume();
    cerr << datasetTemp.getMinValue() << " " << datasetTemp.getMaxValue() << " " << volScaleFactor * volMinus << endl;
  }
  {
    qc::ScalarArray< RealType, qc::QC_3D> datasetTemp ( dataset, aol::DEEP_COPY );
    datasetTemp.addToAll ( +thresholdDelta );
    tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > gridTmp ( qc::GridSize<qc::QC_3D>::createFrom ( datasetTemp ) );
    gridTmp.setDomainFrom ( datasetTemp );
    gridTmp.detectAndInitVirtualNodes();
    volPlus = gridTmp.getTotalInnerVolume();
    cerr << datasetTemp.getMinValue() << " " << datasetTemp.getMaxValue() << " " << volScaleFactor * volPlus << endl;
  }
  return ( volScaleFactor * ( volPlus - volMinus ) );
}


void voxelsToCubeVisualization ( const qc::ScalarArray<double, qc::QC_3D> & Levelset, const char* Filename ) {
  aol::Vec3<int> imgSize = Levelset.getSize();
  double ha = 1.0 / aol::Max ( imgSize[0] - 1, imgSize[1] - 1, imgSize[2] - 1 );

  aol::TriangMesh<double> tMesh;
  tMesh.reserve ( ( imgSize[0] + 1 ) * ( imgSize[1] + 1 ) * ( imgSize[2] + 1 ), Levelset.size() / 10 /* wild estimate */ );

  std::map< int, aol::Vec3<double> > vertices;

  // store all points, not only used ones (for now)
  aol::ProgressBar<> pbp ( "Storing Nodes" );
  pbp.start ( 2 * ( imgSize[0] + 1 ) * ( imgSize[1] + 1 ) * ( imgSize[2] + 1 ) );
  for ( qc::RectangularIterator<qc::QC_3D> rit ( aol::Vec3<int> ( 0, 0, 0 ), imgSize + aol::Vec3<int> ( 1, 1, 1 ) ); rit.notAtEnd(); ++rit, pbp++ ) {
    vertices[ qc::ILexCombine3 ( (*rit)[0], (*rit)[1], (*rit)[2], Levelset.getNumX() + 1, Levelset.getNumY() + 1 ) ]
      = ha * aol::Vec3<double> ( (*rit)[0] - 0.5, (*rit)[1] - 0.5, (*rit)[2] - 0.5 );
  }

  for ( std::map< int, aol::Vec3<double> >::const_iterator it = vertices.begin(); it != vertices.end(); ++it, pbp++ ) {
    tMesh.pushBackVertex ( it->second );
  }
  pbp.finish();

  aol::ProgressBar<> pb ( "Converting" );
  pb.start ( Levelset.size() );
  // draw all cubes that are not surrounded (face-nbhd) by interior cubes
  for ( qc::RectangularIterator<qc::QC_3D> rit ( Levelset ); rit.notAtEnd(); ++rit, pb++ ) {
    if ( Levelset.get ( *rit ) < 0 ) {
      // bool displayCube = false;
      for ( int ox = -1; ox < 2; ox += 2 ) {
        for ( int oy = -1; oy < 2; oy += 2 ) {
          for ( int oz = -1; oz < 2; oz += 2 ) {
            if ( ( (*rit)[0] + ox >= 0 ) &&  ( (*rit)[0] + ox < Levelset.getNumX() ) &&
                 ( (*rit)[1] + oy >= 0 ) &&  ( (*rit)[1] + oy < Levelset.getNumY() ) &&
                 ( (*rit)[2] + oz >= 0 ) &&  ( (*rit)[2] + oz < Levelset.getNumZ() ) &&
                 ( Levelset.get ( (*rit)[0] + ox, (*rit)[1] + oy, (*rit)[2] + oz )  ) ) {
              // displayCube = true;
              ox = oy = oz = 3;
            }
          }
        }
      }

      const int
        node0 = qc::ILexCombine3 ( (*rit)[0] + 0, (*rit)[1] + 0, (*rit)[2] + 0, Levelset.getNumX() + 1, Levelset.getNumY() + 1 ),
        node1 = qc::ILexCombine3 ( (*rit)[0] + 1, (*rit)[1] + 0, (*rit)[2] + 0, Levelset.getNumX() + 1, Levelset.getNumY() + 1 ),
        node2 = qc::ILexCombine3 ( (*rit)[0] + 0, (*rit)[1] + 1, (*rit)[2] + 0, Levelset.getNumX() + 1, Levelset.getNumY() + 1 ),
        node3 = qc::ILexCombine3 ( (*rit)[0] + 1, (*rit)[1] + 1, (*rit)[2] + 0, Levelset.getNumX() + 1, Levelset.getNumY() + 1 ),
        node4 = qc::ILexCombine3 ( (*rit)[0] + 0, (*rit)[1] + 0, (*rit)[2] + 1, Levelset.getNumX() + 1, Levelset.getNumY() + 1 ),
        node5 = qc::ILexCombine3 ( (*rit)[0] + 1, (*rit)[1] + 0, (*rit)[2] + 1, Levelset.getNumX() + 1, Levelset.getNumY() + 1 ),
        node6 = qc::ILexCombine3 ( (*rit)[0] + 0, (*rit)[1] + 1, (*rit)[2] + 1, Levelset.getNumX() + 1, Levelset.getNumY() + 1 ),
        node7 = qc::ILexCombine3 ( (*rit)[0] + 1, (*rit)[1] + 1, (*rit)[2] + 1, Levelset.getNumX() + 1, Levelset.getNumY() + 1 );
      tMesh.pushBackTriang ( aol::Vec3<int> ( node0, node1, node2 ) );   tMesh.pushBackTriang ( aol::Vec3<int> ( node1, node2, node3 ) );
      tMesh.pushBackTriang ( aol::Vec3<int> ( node4, node5, node6 ) );   tMesh.pushBackTriang ( aol::Vec3<int> ( node5, node6, node7 ) );
      tMesh.pushBackTriang ( aol::Vec3<int> ( node0, node2, node4 ) );   tMesh.pushBackTriang ( aol::Vec3<int> ( node2, node4, node6 ) );
      tMesh.pushBackTriang ( aol::Vec3<int> ( node1, node3, node5 ) );   tMesh.pushBackTriang ( aol::Vec3<int> ( node3, node5, node7 ) );
      tMesh.pushBackTriang ( aol::Vec3<int> ( node0, node1, node4 ) );   tMesh.pushBackTriang ( aol::Vec3<int> ( node1, node4, node5 ) );
      tMesh.pushBackTriang ( aol::Vec3<int> ( node2, node3, node6 ) );   tMesh.pushBackTriang ( aol::Vec3<int> ( node3, node6, node7 ) );
    }
  }
  pb.finish();

  tMesh.saveAsPLY ( Filename );
}

template< typename RealType >
void saveSlice ( const qc::ScalarArray< RealType, qc::QC_3D> &dataset, const short slice, const char* fnroot ) {
  char outfile[1024];
  qc::ScalarArray< unsigned char, qc::QC_2D > slic ( dataset.getNumX(), dataset.getNumY() );
  for ( qc::RectangularIterator<qc::QC_3D> rit ( qc::CoordType( 0, 0, slice ), qc::CoordType( dataset.getNumX(), dataset.getNumY(), slice + 1 ) ); rit.notAtEnd(); ++rit ) {
    slic.set ( (*rit)[0], (*rit)[1], ( dataset.get ( *rit ) < 0 ? 255 : 0 ) );
  }
  sprintf ( outfile, "%s_binarized.png", fnroot );
  slic.savePNG ( outfile );

  qc::ScalarArray< RealType, qc::QC_3D > slices ( dataset.getNumX(), dataset.getNumY(), 2 );
  for ( qc::RectangularIterator<qc::QC_3D> rit ( qc::CoordType( 0, 0, slice ), qc::CoordType( dataset.getNumX(), dataset.getNumY(), slice + 2 ) ); rit.notAtEnd(); ++rit ) {
    slices.set ( (*rit)[0], (*rit)[1], (*rit)[2] - slice, dataset.get ( *rit ) );
  }

  tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > gridTmp ( qc::GridSize<qc::QC_3D>::createFrom ( slices ) );
  gridTmp.setDomainFrom ( slices );
  gridTmp.detectAndInitVirtualNodes();

  tpcfe::CFEInterfaceTriangulationGenerator< tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > > itg ( gridTmp );
  itg.determineSliceTriangulation ( qc::QC_Z, 0, 1.0e-9, -1 );
  aol::TriangMesh<float> tMesh;
  itg.writeToTriangMesh ( tMesh );
  sprintf ( outfile, "%s_CFESlice.ply", fnroot );

  tMesh.saveAsUDPLY ( outfile );
}

#endif
