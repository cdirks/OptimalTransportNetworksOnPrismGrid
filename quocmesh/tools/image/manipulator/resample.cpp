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
 *  \brief Up- or downsampling of an image
 *
 *  Up- or downsampling of an image (i.e. doubling or halving the resolution in each direction). \\
 *  Usage: resample input output sampleMode \\
 *  Possible sampleModes: UP, DOWN.
 *
 *  \author Oliver Nemitz
 */

#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>
#include <scalarArray.h>
#include <string.h>

#include <qmException.h>
#include <fstream>
#include <auxiliary.h>
#include <restriction.h>
#include <prolongation.h>

typedef double RealType;
using namespace aol::color;

// get the first character of the header
char getFirstHeaderChar ( const char* fileName ) {
  aol::ipfstream file ( fileName );
  qc::ArrayHeader header;
  qc::ReadArrayHeader ( file, header );
  file.close();

  return header.magic[0];
}

void finish ( const char *msg ) {
  cerr << red << msg << reset << endl;
  cerr << "usage: resample input output sampleMode\n" << reset;
  cerr << "with possible sampleMode: UP, DOWN. \n";
  exit ( 0 );
}

enum sampleMode {
  UP,
  DOWN
};

int main ( int argc, char **argv ) {
  try {
    if ( argc < 4 ) {
      cerr << "Up- or downsampling of an image (i.e. doubling or halving the resolution in each direction). \n";
      cerr << red << "usage: " << argv[0] << " input output sampleMode\n" << reset;
      cerr << red << "Possible sampleModes: UP, DOWN. \n" << reset;
      return EXIT_FAILURE;
    }

    // get the sample mode
    sampleMode sMode = UP;
    if ( !strcmp ( argv[3], "DOWN" ) ) {
      sMode = DOWN;
    } else if ( !strcmp ( argv[3], "UP" ) ) {
      sMode = UP;
    } else {
      finish ( "Error: No valid sampleMode given!" );
    }

    // get possible further arguments
    int numSteps = 1;
    if ( argc >= 5 ) numSteps = atoi ( argv[4] );

    char c = getFirstHeaderChar ( argv[1] );

    // ------------------------- --------------------------------------------
    // ------------------------- 3D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( c == 'Q' ) {
      qc::ScalarArray<RealType, qc::QC_3D> inputImg ( argv[1] );

      // get a grid according to the image
      int N = inputImg.getNumX ();
      int d = qc::logBaseTwo ( N );
      qc::GridDefinition grid ( d, qc::QC_3D );

      // construct a grid according to the sampled image
      if ( sMode == UP ) d += numSteps;
      else d -= numSteps;
      if ( d <= 0 ) finish ( "Error: Grid is too coarse for desired downsampling!" );
      qc::GridDefinition resultGrid ( d, qc::QC_3D );

      qc::ScalarArray<RealType, qc::QC_3D> outputImg ( resultGrid );

      if ( sMode == UP ) {
        cerr << blue << "Upsampling " << red << argv[1] << blue << "...";

        qc::ProlongOp<double> Prolong ( grid, resultGrid );
        Prolong.apply ( inputImg, outputImg );

      } else {
        cerr << blue << "Downsampling " << red << argv[1] <<  blue << "...";

        qc::RestrictOp<double, qc::STD_QUOC_RESTRICT> Restrict ( resultGrid, grid );
        Restrict.apply ( inputImg, outputImg );
      }

      cerr << "done! Saving '" << red << argv[2] << blue << "'..." << endl << aol::color::black;
      outputImg.save ( argv[2], qc::SaveTypeTrait<RealType>::BinarySaveType );
    }

    // ------------------------- --------------------------------------------
    // ------------------------- 2D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( c == 'P' ) {
      qc::ScalarArray<RealType, qc::QC_2D> inputImg ( argv[1] );

      // get a grid according to the image
      int N = inputImg.getNumX ();
      int d = qc::logBaseTwo ( N );
      qc::GridDefinition grid ( d, qc::QC_2D );

      // construct a grid according to the sampled image
      if ( sMode == UP ) d += numSteps;
      else d -= numSteps;
      if ( d <= 0 ) finish ( "Error: Grid is too coarse for desired downsampling!" );
      qc::GridDefinition resultGrid ( d, qc::QC_2D );

      qc::ScalarArray<RealType, qc::QC_2D> outputImg ( resultGrid );
      outputImg.quietMode = true;

      if ( sMode == UP ) {
        cerr << blue << "Upsampling " << red << argv[1] << blue << "...";

        qc::ProlongOp<double> Prolong ( grid, resultGrid );
        Prolong.apply ( inputImg, outputImg );

      } else {
        cerr << blue << "Downsampling " << red << argv[1] <<  blue << "...";

        qc::RestrictOp<double, qc::STD_QUOC_RESTRICT> Restrict ( resultGrid, grid );
        Restrict.apply ( inputImg, outputImg );
      }

      cerr << "done! Saving '" << red << argv[2] << blue << "'..." << endl << aol::color::black;
      qc::recognizeEndingAndSave2d ( outputImg, argv[2] );
    }

    if ( c != 'P' && c != 'Q' )  finish ( "Can't determine whether your argument is a 2d or a 3d-image!" );

    cerr << blue << "done! Thanx for using resample.\n" << reset;
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
