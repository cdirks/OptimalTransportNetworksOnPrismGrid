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
 *  \brief combines two 3D images to a single one with alternating stripes of both inputs
 *
 *  This tool loads two 3D images (not necessarily of the same size) and produces
 *  a third 3D image consisting of alternating diagonal stripes from both original
 *  images (i.e. a mixture of both images).
 *
 *  Also, three 2D images are produced of cross-sections in all three directions at the middle of the 3D result.
 *
 *  Usage: generateSlantedStripesCombinedImage3D [GridLevelResult] [NumberOfStripes] [InputImage1] [InputImage2] [FileNameResult]
 *
 *  \author Wirth
 */

#include <scalarArray.h>
#include <scalarArray.h>

typedef double RealType;

int main ( int argc, char **argv ) {

  // in case of wrong number of parameters, notify user of correct usage
  if ( argc != 6 ) {
    cerr << "usage: " << argv[0] << " <level> <factor> <input1> <input2> <output>\n\n";
    return EXIT_FAILURE;
  }

  try {
    // define the grid of the result and its width
    qc::GridDefinition grid ( atoi ( argv[1] ) , qc::QC_3D );
    const int size = grid.getWidth();

    // load the input images and pad them into images of the correct size
    cerr << "load input images..." << endl;

    qc::ScalarArray<RealType, qc::QC_3D> in ( argv[3] );
    qc::ScalarArray<RealType, qc::QC_3D> in2 ( argv[4] );

    qc::ScalarArray<RealType, qc::QC_3D> input ( grid );
    qc::ScalarArray<RealType, qc::QC_3D> input2 ( grid );

    input.padFrom ( in );
    input2.padFrom ( in2 );

    cerr << "done" << endl;

    // produce combination of both images and save the result
    cerr << "produce stripe image..." << endl;

    qc::ScalarArray<RealType, qc::QC_3D> target ( grid );

    int stripeWidth = 2 * size / atoi ( argv[2] );
    for ( int X = 0; X < size; X++ ) {
      for ( int Y = 0; Y < size; Y++ ) {
        for ( int Z = 0; Z < size; Z++ ) {
          // put data from both images in alternating stripes
          if ( ( ( X + Y + Z ) / stripeWidth ) % 2  ) {
            target.set ( X, Y, Z, input.get ( X, Y, Z ) );
          } else {
            target.set ( X, Y, Z, input2.get ( X, Y, Z ) );
          }
        }
      }
    }

    target.save ( argv[5], qc::PGM_DOUBLE_BINARY );

    cerr << "done" << endl;

    // save cross-sections through the middle in x-, y- and z-direction
    cerr << "produce cross-sections..." << endl;

    qc::ScalarArray<RealType, qc::QC_2D> slice ( size, size );

    target.getSlice ( qc::QC_X, size >> 1, slice );
    slice.save ( "slice_x.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );

    target.getSlice ( qc::QC_Y, size >> 1, slice );
    slice.save ( "slice_y.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );

    target.getSlice ( qc::QC_Z, size >> 1, slice );
    slice.save ( "slice_z.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );

    cerr << "done" << endl;

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
  return EXIT_SUCCESS;
}
