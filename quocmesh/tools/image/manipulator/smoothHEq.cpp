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
 *  \brief Applies $n$ steps of the heat equation with $\tau=h$ to the input-image.
 *
 *  Applies $n$ steps of the heat equation with $\tau=h$ to the input-image
 *  (works both in 2d and 3d). Standard filter width is $w*h$ with $w=4$.
 *  Usage: smoothHEq inputImg outputImg n [w]
 *  {\bf Attention:} The output-image will only be saved with double values if the suffix is '.bz2'!
 *
 *  \author Oliver Nemitz
 */

#include <quoc.h>
#include <scalarArray.h>
#include <linearSmoothOp.h>

typedef double RealType;
using namespace aol::color;

// get the first character of the header
char getFirstHeaderChar ( char* fileName ) {
  aol::ipfstream file ( fileName );
  qc::ArrayHeader header;
  qc::ReadArrayHeader ( file, header );
  file.close();

  return header.magic[0];
}

int main ( int argc, char **argv ) {
  try {
    if ( argc < 4 ) {
      cerr << "Applies n steps of the heat equation with tau=h to the image. Standard filter width is w*h with w=4.\n";
      cerr << aol::color::red << "usage: " << argv[0] << " input output n <w>\n" << aol::color::reset;
      return EXIT_FAILURE;
    }

    int n = atoi ( argv[3] );         // number of steps
    double w = 4.;
    if ( argc == 5 ) {
      w = atof ( argv[4] );           // change the filter-width
      cerr << "Using" << green << " new filter width " << w << "*h." << endl << reset;
    }
    cerr << "Loading image '" << green << argv[1] << reset << "'...\n";

    char c = getFirstHeaderChar ( argv[1] );

    // ------------------------- --------------------------------------------
    // ------------------------- 3D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( c == 'Q' ) {
      qc::ScalarArray<RealType, qc::QC_3D> image ( argv[1] );
      image.quietMode = true;

      // get a grid according to the image
      int N = image.getNumX ();
      int d = qc::logBaseTwo ( N );
      qc::GridDefinition grid ( d, qc::QC_3D );

      // decleare the operator
      qc::LinearSmoothOp<RealType> smoothOp;
      smoothOp.setCurrentGrid ( grid );
      smoothOp.setSigma ( w * grid.H() );

      cerr << endl;
      for ( int i = 0; i < n; i++ ) {
        cerr << "\r3D: Smoothing the image, step " << red << i + 1 << reset << " of " << n << "...";
        smoothOp.apply ( image, image );
      }
      cerr << "done! \n\nSaving '" << green << argv[2] << reset << "'..." << reset;
      image.save ( argv[2], qc::SaveTypeTrait<RealType>::BinarySaveType );
    }


    // ------------------------- --------------------------------------------
    // ------------------------- 2D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( c == 'P' ) {
      qc::ScalarArray<RealType, qc::QC_2D> image ( argv[1] );
      image.quietMode = true;

      // get a grid according to the image
      int N = image.getNumX ();
      int d = qc::logBaseTwo ( N );
      qc::GridDefinition grid ( d, qc::QC_2D );

      // decleare the operator
      qc::LinearSmoothOp<RealType> smoothOp;
      smoothOp.setCurrentGrid ( grid );
      smoothOp.setSigma ( w * grid.H() );

      cerr << endl;
      for ( int i = 0; i < n; i++ ) {
        cerr << "\r2D: Smoothing the image, step " << red << i + 1 << reset << " of " << n << "...";
        smoothOp.apply ( image, image );
      }

      cerr << "done! Saving...\n" << aol::color::black;
      qc::recognizeEndingAndSave2d ( image, argv[2] );
    }





    cerr << "done! \nThanx for using smoothHEq (Non-registered version!!)\n";
    cerr << aol::color::black;
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
