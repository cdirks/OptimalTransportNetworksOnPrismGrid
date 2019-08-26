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
 * \brief Deforms the colored input PNG as specified by the input displacement field.
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
    char inputImageFileName[1024];
    char inputFileNameDefX[1024];
    char inputFileNameDefY[1024];
    char outputFileName[1024];

    if ( ( argc < 5 ) || ( argc > 7 ) ) {
      cerr << "USAGE: " << argv[0] << " <input_png> <input_file def_x> <input_file def_y> <out_file> [invertDev] [NNinterp]" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    } else {
      sprintf ( inputImageFileName, "%s",  argv[1] );
      sprintf ( inputFileNameDefX, "%s",  argv[2] );
      sprintf ( inputFileNameDefY, "%s",  argv[3] );
      sprintf ( outputFileName, "%s",  argv[4] );
    }

    const bool useInverseDeformation = ( argc > 5 ) ? ( atoi ( argv[5] ) != 0 ) : false;
    const bool NearestNeighborInterpolation = ( argc > 6 ) ? ( atoi ( argv[6] ) != 0 ) : false;

    qc::deformAndSaveColoredImage<ConfType> ( inputImageFileName, inputFileNameDefX, inputFileNameDefY, outputFileName, useInverseDeformation, NearestNeighborInterpolation );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
