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
 * \brief Deforms the 2D quoc input array as specified by the input displacement field and saves the result with double precision.
 *
 * \author Berkels
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
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType, qc::QC_2D, 3> > ConfType;

int main ( int argc, char **argv ) {

  try {
    if ( ( argc < 5 ) || ( argc > 7 ) ) {
      cerr << "USAGE: " << argv[0] << " <input_png> <input_file def_x> <input_file def_y> <out_file> [invertDev] [NNinterp]" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputImageFileName = argv[1];
    const string inputFileNameDefX = argv[2];
    const string inputFileNameDefY = argv[3];
    const string outputFileName = argv[4];

    const bool useInverseDeformation = ( argc > 5 ) ? ( atoi ( argv[5] ) != 0 ) : false;
    const bool NearestNeighborInterpolation = ( argc > 6 ) ? ( atoi ( argv[6] ) != 0 ) : false;

    qc::ScalarArray<RType, qc::QC_2D> dataArray ( inputImageFileName.c_str() );
    qc::MultiArray<RType, qc::QC_2D> deformation ( inputFileNameDefX, inputFileNameDefY );
    ConfType::InitType grid ( dataArray.getSize() );
    qc::ScalarArray<RType, qc::QC_2D> deformedArray ( grid );

    if ( useInverseDeformation == false )
      qc::DeformImage<ConfType>( dataArray, grid, deformedArray, deformation, true, 0, NearestNeighborInterpolation );
    else
      qc::InvDeformImage<ConfType>( dataArray, grid, deformedArray, deformation );

    deformedArray.save ( outputFileName.c_str(), qc::PGM_DOUBLE_BINARY );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
