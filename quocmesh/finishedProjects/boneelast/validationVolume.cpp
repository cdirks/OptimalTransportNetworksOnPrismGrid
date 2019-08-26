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
#include <tpCFEUtils.h>

#include "thresholding.h"


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

    const RealType
      HistoAnalyzeMin = params.getDouble ( "HistoAnalyzeMin" ),
      HistoAnalyzeMax = params.getDouble ( "HistoAnalyzeMax" );
    // ThresholdDelta  = params.getDouble ( "ThresholdDelta" );

    const RealType origRes = params.getDouble ( "origRes" );

    const RealType voxelThresholdFinest = ( params.hasVariable ( "DeviceThreshold" ) ?
              params.getDouble ( "DeviceThreshold" ) :
              RiCaThreshold<RealType> ( ROI, HistoAnalyzeMin, HistoAnalyzeMax, aol::Cub ( ROI.getNumZ() * origRes ) ) );

    qc::ScalarArray< RealType, qc::QC_3D > ROIcopy ( ROI, aol::STRUCT_COPY );
    thresholdAndDepoempelize< RealType > ( ROI, voxelThresholdFinest, aol::Cub ( ROI.getNumZ() * origRes ), ROIcopy );

    // cerr << "CFE volume sensitivity for ThresholdDelta = " << aol::longScientificFormat ( ThresholdDelta )
    //      << " is " << aol::longScientificFormat ( determineCFEVolumeSensitivity ( ROIcopy, ThresholdDelta, aol::Cub ( ROI.getNumZ() * origRes ) ) ) << endl;

    for ( int resS = 2; resS < 10; ++resS ) {

      qc::ScalarArray< RealType, qc::QC_3D > ROIdown;

      resampleDatasetByIntegerFactor<RealType> ( ROI, ROIdown, resS  );

      cerr << "integer resample factor " << resS << " resulting in resolution " << origRes * resS << " m per voxel, size " << ROIdown.getSize() << endl;

      qc::ScalarArray< RealType, qc::QC_3D > ROIdownSeg ( ROIdown, aol::STRUCT_COPY );

      thresholdAndDepoempelize< RealType > ( ROIdown, voxelThresholdFinest, aol::Cub ( origRes * resS * ROIdown.getNumZ()  ), ROIdownSeg );

      // cerr << "CFE volume sensitivity for ThresholdDelta = " << aol::longScientificFormat ( ThresholdDelta )
      //      << " is " << aol::longScientificFormat ( determineCFEVolumeSensitivity ( ROIdownSeg, ThresholdDelta, aol::Cub ( origRes * resS * ROIdown.getNumZ() ) ) ) << " m^3" << endl;

#if 0
      const short slicNo = 120 / resS;
      char fnroot[1024];
      sprintf ( fnroot, "out/testSlice%02d", resS );
      saveSlice<RealType> ( ROIdownSeg, slicNo, fnroot );
#endif

    }

  } catch ( aol::Exception &e ) {
    e.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}

