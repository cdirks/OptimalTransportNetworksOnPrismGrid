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
 *  \brief change contrast in an image
 *
 *  This tool loads a 3D image, scales its pixel values such that they lie
 *  within [0,1] and then takes each pixel value to the power of some constant 1/gamma.
 *
 *  For 2D just replace ScalarArray QC\_3D by ScalarArray QC\_2D.
 *
 *  Usage: rootContrastChange [3dInputImage] [gamma] [FileNameResult]
 *
 *  \author Wirth
 */

#include <scalarArray.h>
#include <scalarArray.h>

typedef double RealType;

int main ( int argc, char **argv ) {

  // in case of wrong number of parameters, notify user of correct usage
  if ( argc != 4 ) {
    cerr << "usage: " << argv[0] << " <input> <gamma> <output>\n\n";
    return EXIT_FAILURE;
  }

  try {

    // load the input image and scale pixel values to the interval [0,1]
    qc::ScalarArray<RealType, qc::QC_2D> input ( argv[1] );
    input.addToAll ( - input.getMinValue() );
    input /= input.getMaxValue();

    // take the gammath root of each pixel value
    RealType gamma = atof ( argv[2] );
    const int size = input.size();
    for ( int i = 0; i < size; i++ ) {
      input[i] = pow ( input[i], 1. / gamma );
    }

    // save the result
    input.save ( argv[3], qc::PGM_DOUBLE_BINARY );


  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
  return EXIT_SUCCESS;
}
