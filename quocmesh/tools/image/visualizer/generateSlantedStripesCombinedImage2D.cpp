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
 *  \brief combines two 2D images to a single one with alternating stripes of both inputs
 *
 *  This tool loads two 2D images (not necessarily of the same size) and produces
 *  a third 2D image consisting of alternating diagonal stripes from both original
 *  images (i.e. a mixture of both images).
 *
 *  Also, three 2D images are produced of cross-sections in all three directions at the middle of the 3D result.
 *
 *  Usage: generateSlantedStripesCombinedImage2D [NumberOfStripes] [InputImage1] [InputImage2] [FileNameResult]
 *
 *  \author Wirth
 */

#include <scalarArray.h>

typedef double RealType;

int main ( int argc, char **argv ) {

  // in case of wrong number of parameters, notify user of correct usage
  if ( argc != 6 ) {
    cerr << "usage: " << argv[0] << " <factor> <input1> <input2> <output>\n\n";
    return EXIT_FAILURE;
  }

  try {
    // load the input images
    cerr << "load input images..." << endl;

    qc::ScalarArray<RealType, qc::QC_2D> input ( argv[2] );
    qc::ScalarArray<RealType, qc::QC_2D> input2 ( argv[3] );
    const int size = input.getNumX();

    // throw exception, if image dimensions mismatch
    if ( !input.equalDimensions ( input2 ) )
      throw aol::Exception ( "input images do not have the same dimension" );

    cerr << "done" << endl;

    // produce combination of both images and save the result
    cerr << "produce stripe image..." << endl;

    qc::ScalarArray<RealType, qc::QC_2D> target ( input );

    int stripeWidth = 2 * size / atoi ( argv[1] );
    for ( int X = 0; X < size; X++ ) {
      for ( int Y = 0; Y < size; Y++ ) {
        // put data from both images in alternating stripes
        if ( ( ( X + Y ) / stripeWidth ) % 2  ) {
          target.set ( X, Y, input.get ( X, Y ) );
        } else {
          target.set ( X, Y, input2.get ( X, Y ) );
        }
      }
    }

    target.save ( argv[4], qc::PGM_UNSIGNED_CHAR_BINARY );

    cerr << "done" << endl;

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
  return EXIT_SUCCESS;
}
