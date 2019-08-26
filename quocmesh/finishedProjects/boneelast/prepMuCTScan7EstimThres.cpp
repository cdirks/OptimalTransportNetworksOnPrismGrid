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

// estimate a segmentation threshold for an already determined ROI (Uwe's HR CT scans from Vienna)

#include <configurators.h>
#include <linearSmoothOp.h>
#include <parameterParser.h>
#include <scalarArray.h>

#include <tpCFEGrid.h>
#include <tpCFELevelsets.h>


// todo:  use code from thresholding.h"

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


template< typename RealType >
void thresholdAndDepoempelize ( const qc::ScalarArray< RealType, qc::QC_3D> &dataset, const RealType threshold, const RealType volScaleFactor, qc::ScalarArray< RealType, qc::QC_3D > &dest ) {

  dest = dataset;
  dest.addToAll ( -threshold );
  dest *= -1.0;

  tpcfe::TopBottomPoempelRemover<RealType> tbpr;
  tbpr.applySingle ( dest );
  tpcfe::SmallComponentsPoempelRemover<RealType> scpr;
  scpr.applySingle ( dest );

  {
    tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > gridTmp ( qc::GridSize<qc::QC_3D>::createFrom ( dest ) );
    gridTmp.setDomainFrom ( dest );
    gridTmp.detectAndInitVirtualNodes();
    RealType totVol = gridTmp.getTotalInnerVolume();
    cerr << "Total segmented volue = " << aol::detailedFormat ( totVol )
         << " for voxel threshold " << aol::detailedFormat ( threshold )
         << " corresponding to " << aol::detailedFormat ( volScaleFactor * totVol ) << " m^3" << endl;
  }
}



typedef double RealType;

int main ( int argc, char **argv ) {
  if ( argc != 2 ) {
    cerr << "Usage: " << argv[0] << " parameter_file" << endl;
    abort();
  }

  try {

    aol::ParameterParser params ( argv[1] );

    char loadFilename[1024];
    params.getString ( "loadFilename", loadFilename );
    qc::ScalarArray< RealType, qc::QC_3D > ROI ( loadFilename );
    cerr << "Loaded dataset with dimensions " << ROI.getSize() << endl;
    cerr << "Minimal and maximal values are: " << ROI.getMinValue() << " " << ROI.getMaxValue() << endl;

    char histoFNRoot[1024];
    params.getString ( "histoFNRoot", histoFNRoot );

    const RealType
      HistoAnalyzeMin = params.getDouble ( "HistoAnalyzeMin" ),
      HistoAnalyzeMax = params.getDouble ( "HistoAnalyzeMax" ),
      outResolution = params.getDouble("outResolution");

    RealType resampleFactor = params.getReal<RealType> ( "resampleFactor" ); // will be modified

    RealType volScaleFactor = 1.0;

    const RealType voxelThresholdBeforeDenoising = RiCaThreshold<RealType> ( ROI, HistoAnalyzeMin, HistoAnalyzeMax, outResolution, histoFNRoot );

    qc::ScalarArray< RealType, qc::QC_3D > ROIdown;

    resampleDatasetByFactor<RealType> ( ROI, outResolution, resampleFactor, volScaleFactor, ROIdown );

    if ( params.hasVariable ( "numLinsmoothSteps" ) ) {
      {
        typedef qc::RectangularGridConfigurator< RealType, qc::QC_3D, aol::GaussQuadrature< RealType, qc::QC_3D, 3 > > ConfigType;
        tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > grid ( qc::GridSize<qc::QC_3D>::createFrom ( ROIdown ) );
        qc::GeneralLinearSmoothOp< ConfigType > lso ( grid, params.getDouble ( "linsmoothSigmaFactor" ) * grid.H() );
        lso.setSolverParameters ( 1.0e-16, 10000 );
        cerr << "computing " <<  params.getInt ( "numLinsmoothSteps" ) << " steps of isotropic diffusion with sigma = " << params.getDouble ( "linsmoothSigmaFactor" ) * grid.H() << endl;
        for ( int iter = 0 ; iter < params.getInt ( "numLinsmoothSteps" ) ; ++iter ) {
          lso.applySingle ( ROIdown );
        }
      }

      RiCaThreshold<RealType> ( ROI, HistoAnalyzeMin, HistoAnalyzeMax, volScaleFactor, histoFNRoot );
    }

    qc::ScalarArray< RealType, qc::QC_3D > ROIdownSeg ( ROIdown, aol::STRUCT_COPY );
    thresholdAndDepoempelize< RealType > ( ROIdown, voxelThresholdBeforeDenoising, volScaleFactor, ROIdownSeg );

    cerr << "Save dataset with dimensions " << ROIdownSeg.getSize() << " ...";
    char geometryFilename[1024];
    params.getString ( "geometryFilename", geometryFilename );

    ROIdownSeg.save( geometryFilename, qc::PGM_FLOAT_BINARY );

#ifdef SAVE_HOMOG_CUBES
    // PROBABLY NEEDS TO BE ADAPTED!!
    {
      qc::ScalarArray<RealType, qc::QC_3D> cubeA ( 129, 129, 129 ), cubeB ( 129, 129, 129 );
      qc::CoordType
        offsetA ( ( ROI.getNumX() - 129 ) / 2, ( ROI.getNumY() - 129 ) / 2, 5 ),
        offsetB ( ( ROI.getNumX() - 129 ) / 2, ( ROI.getNumY() - 129 ) / 2, ROI.getNumZ() - 135 );
      for ( qc::RectangularIterator<qc::QC_3D> rit ( cubeA ); rit.notAtEnd(); ++rit ) {
        cubeA.set ( *rit, ROIdown.get ( *rit + offsetA ) );
        cubeB.set ( *rit, ROIdown.get ( *rit + offsetB ) );
      }

      tpcfe::AllFacePoempelRemover<RealType> afpr;
      afpr.applySingle ( cubeA );
      afpr.applySingle ( cubeB );

      char cubeFileName[1024];
      sprintf ( cubeFileName, "%s_cubeA.dat.bz2", histoFNRoot );      cubeA.save ( cubeFileName, qc::PGM_FLOAT_BINARY );
      sprintf ( cubeFileName, "%s_cubeB.dat.bz2", histoFNRoot );      cubeB.save ( cubeFileName, qc::PGM_FLOAT_BINARY );
    }
#endif

  } catch ( aol::Exception &e ) {
    e.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}



#ifdef DO_CFE_HISTOGRAM
    // this probably no longer works without adapting code

    { // determine threshold based on CFE-segmented volumes
      ROI.addToAll ( GVOffset );

      RealType estimatedThreshold = ROI.sum() / ROI.size();

      const int numHistoSteps = params.getInt ( "numCFEHistoSteps" );
      const RealType histoRange = ROI.getMaxValue(), histoStep = histoRange / ( numHistoSteps - 1 );
      {
        aol::Vector< RealType > currentVolumes ( numHistoSteps + 1 );
        for ( int i = 0; i < numHistoSteps+1; ++i ) {
          const RealType currentHistoValue = i * histoStep;
          if ( currentHistoValue - GVOffset > HistoAnalyzeMin && currentHistoValue  - GVOffset < HistoAnalyzeMax ) {
            qc::ScalarArray< RealType, qc::QC_3D > tmpImage ( ROI, aol::DEEP_COPY );
            tmpImage.addToAll ( - currentHistoValue );
            tmpImage *= -1;
            tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > gridTmp ( qc::GridSize<qc::QC_3D>::createFrom ( tmpImage ), qc::QC_3D );
            gridTmp.setDomainFrom ( tmpImage );
            gridTmp.detectAndInitVirtualNodes();
            const RealType currentVolume = gridTmp.getTotalInnerVolume();
            currentVolumes[i] = currentVolume;
            cerr << aol::detailedFormat ( currentHistoValue - GVOffset ) << " " << aol::detailedFormat ( currentVolumes[i] ) << endl;
          }
        }

        aol::Vector<RealType> histogram ( numHistoSteps );
        for ( int i = 0; i < numHistoSteps; ++i ) {
          if ( i * histoStep - GVOffset > HistoAnalyzeMin && (i+1) * histoStep - GVOffset < HistoAnalyzeMax ) {
            histogram[i] = currentVolumes[i] - currentVolumes[i+1];
          }
        }

        sprintf ( histoFN, "%s_HistoCFE.dat", histoFNRoot );
        ofstream histo ( histoFN );
        for ( int i = 0; i < histogram.size(); ++i )
          histo << aol::detailedFormat ( (i+0.5) * histoStep - GVOffset ) << " " << aol::detailedFormat ( histogram[i] ) << endl;

        sprintf ( histoFN, "%s_HistoCFEIntegrated.dat", histoFNRoot );
        ofstream histoIntegr ( histoFN );
        for ( int i = 0; i < histogram.size(); ++i ) {
          histoIntegr << aol::detailedFormat ( (i+0.5) * histoStep - GVOffset ) << " " << aol::detailedFormat ( currentVolumes[i] ) << endl;
        }

        for ( int iter = 0; iter < 10; ++iter ) {
          RealType lowerNum = 0, lowerDen = 0, upperNum = 0, upperDen = 0;
          for ( int i = 0; i < histogram.size(); ++i ) {
            if ( (i+0.5) * histoStep < estimatedThreshold ) {
              lowerNum += (i+0.5) * histoStep * histogram[i];
              lowerDen += 2 * histogram[i];
            } else {
              upperNum += (i+0.5) * histoStep * histogram[i];
              upperDen += 2 * histogram[i];
            }
          }

          estimatedThreshold = lowerNum / lowerDen + upperNum / upperDen;

          cerr << aol::detailedFormat ( estimatedThreshold ) << endl;
        }
      }
      estimatedThreshold -= GVOffset;
      ROI.addToAll ( - GVOffset );

      cerr << "Estimated threshold for segmentation is " << aol::detailedFormat ( estimatedThreshold ) << "; thresholding and removing poempels ... " << endl;
    }
#endif
