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

/** \file
 *  \brief Adds noise to your data (works both in 2d and 3d).
 *
 * Usage: addNoise inputImg outputImg input output noiseType [ratio=0.25 min=0 max=1]
 *
 * Possible noisetypes: SALT_AND_PEPPER, GAUSSIAN, EQUALLY_DISTRIBUTED
 *
 * Example-args for a distance-function: addNoise input.bz2 output.bz2 GAUSSIAN 0.005 0 0.01
 *
 * \attention The output-image will only be saved with double values if the suffix is '.bz2'!
 *
 * \author Nemitz
 */

#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>
#include <scalarArray.h>
#include <string.h>

#include <qmException.h>
#include <fstream>
#include <auxiliary.h>

#include <op.h>

typedef double RealType;

using namespace aol::color;

void finish ( const char *msg ) {
  cerr << red << msg << reset << endl;
  cerr << "usage: addNoise input output noiseType [ratio=0.25 min=0 max=1]\n" << reset;
  cerr << "with possible noisetypes: SALT_AND_PEPPER, GAUSSIAN, EQUALLY_DISTRIBUTED. \n";
  exit ( 0 );
}


int main ( int argc, char **argv ) {
  try {
    if ( argc < 4 ) {
      cerr << "Adds some noise to your nice clean image. \n";
      cerr << red << "Possible noisetypes: SALT_AND_PEPPER, GAUSSIAN, EQUALLY_DISTRIBUTED. \n";
      cerr << red << "usage: " << argv[0] << " input output noiseType [ratio=0.25 min=0 max=1]\n" << reset;
      cerr << "Example for a distance-function: addNoise input.bz2 output.bz2 GAUSSIAN 0.005 0 0.01" << endl;
      cerr << "Attention: If the suffix of the output is '.bz2', it will be saved with double values, otherwise not!" << endl;
      return EXIT_FAILURE;
    }

    // get the noise type
    aol::NoiseOperator<double>::NOISE_TYPE noiseType = aol::NoiseOperator<double>::SALT_AND_PEPPER;
    if ( !strcmp ( argv[3], "SALT_AND_PEPPER" ) ) {
      noiseType = aol::NoiseOperator<double>::SALT_AND_PEPPER;
    } else if ( !strcmp ( argv[3], "GAUSSIAN" ) ) {
      noiseType = aol::NoiseOperator<double>::GAUSSIAN;
    } else if ( !strcmp ( argv[3], "EQUALLY_DISTRIBUTED" ) ) {
      noiseType = aol::NoiseOperator<double>::EQUALLY_DISTRIBUTED;
    } else {
      finish ( "Error: No valid noise_type given!" );
    }

    // get possible further arguments
    double ratio = 0.25;
    double min   = 0.;
    double max   = 1.;

    if ( argc >= 5 ) ratio = atof ( argv[4] );
    if ( argc >= 6 ) min   = atof ( argv[5] );
    if ( argc == 7 ) max   = atof ( argv[6] );

    // declare the noise-operator and set the noise ratio
    aol::NoiseOperator<double> noiseOp ( noiseType, min, max );

    const qc::Dimension dim = qc::getDimensionFromArrayFile ( argv[1] );

    // ------------------------- --------------------------------------------
    // ------------------------- 3D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( dim == qc::QC_3D ) {
      qc::ScalarArray<RealType, qc::QC_3D> inputImg ( argv[1] );
      qc::ScalarArray<RealType, qc::QC_3D> outputImg ( inputImg );
      cerr << blue << "Adding " << red << argv[3] << "-noise" << blue << " to the 2d-image with " << red << "ratio " << red << ratio << blue;
      cerr << " and " << red << "min=" << min << blue << ", " << red << "max=" << red << max << blue << "...";

      noiseOp.applyAdd ( inputImg, outputImg );

      cerr << "done! Saving...\n" << aol::color::black;
      outputImg.save ( argv[2], qc::SaveTypeTrait<RealType>::BinarySaveType );
    }

    // ------------------------- --------------------------------------------
    // ------------------------- 2D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( dim == qc::QC_2D ) {
      qc::ScalarArray<RealType, qc::QC_2D> inputImg ( argv[1] );
      qc::ScalarArray<RealType, qc::QC_2D> outputImg ( inputImg );
      cerr << blue << "Adding " << red << argv[3] << "-noise" << blue << " to the 2d-image with " << red << "ratio " << red << ratio << blue;
      cerr << " and " << red << "min=" << min << blue << ", " << red << "max=" << red << max << blue << "...";

      noiseOp.applyAdd ( inputImg, outputImg );

      cerr << "done! Saving...\n" << aol::color::black;
      qc::recognizeEndingAndSave2d ( outputImg, argv[2] );
    }

    cerr << aol::color::blue << "done! Thanx for using addNoise.\n";
    cerr << aol::color::black;
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
}
