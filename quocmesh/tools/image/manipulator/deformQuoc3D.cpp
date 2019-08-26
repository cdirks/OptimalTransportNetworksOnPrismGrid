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

/**
 * \file
 * \brief Deforms the 3D quoc input array as specified by the (possibly coarser) input displacement field and saves the result with double precision.
 *
 * \author Berkels, Iglesias
 */

#include <aol.h>
#include <op.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <multiArray.h>
#include <levelSetDrawer.h>
#include <timestepSaver.h>
#include <configurators.h>
#include <deformations.h>

typedef double RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_3D, aol::GaussQuadrature<RType, qc::QC_3D, 3> > ConfType;

int main ( int argc, char **argv ) {

  try {
    if ( ( argc < 4 ) || ( argc > 4 ) ) {
      cerr << "USAGE: " << argv[0] << " <input> <input_def_mask> <out_file>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputImageFileName = argv[1];
    const string inputDefFileMask = argv[2];
    const string outputFileName = argv[3];

    qc::ScalarArray<RType, qc::QC_3D> dataArray ( inputImageFileName.c_str() );
    qc::MultiArray<RType, qc::QC_3D> deformation;
    deformation.load( inputDefFileMask.c_str() );
    ConfType::InitType fineGrid ( dataArray.getSize() );
    ConfType::InitType coarseGrid ( deformation[0].getSize() );
    qc::ScalarArray<RType, qc::QC_3D> deformedArray ( fineGrid );

    qc::deformImageWithCoarseDeformation<ConfType> ( dataArray, fineGrid, coarseGrid, deformedArray, deformation );

    deformedArray.save ( outputFileName.c_str(), qc::PGM_DOUBLE_BINARY );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
