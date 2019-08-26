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
 *  \brief compare two PGM images
 *
 *  Compares two pgm (grey) images, computes number of different
 *  pixels and sum of the absolute of the difference and saves the
 *  difference image.
 *
 *  Usage: compareImages inimage1 inimage2 outdiffimage
 *
 *  \author Schwen
 */

#include <iostream>
#include "scalarArray.h"


// compare two pgm images, compute number of different pixels and sum of the absolute of the difference, save difference image.

int main ( int argc, char **argv ) {
  if ( argc == 4 ) {
    qc::ScalarArray<int, qc::QC_2D> image1 ( argv[1] ), image2 ( argv[2] ), diffim ( image1 );
    diffim -= image2;

    int different_pixels = 0, diffsum = 0;

    for ( int i = 0; i < diffim.size(); i++ ) {
      if ( diffim[i] != 0 ) {
        different_pixels++;
        diffim[i] = abs ( diffim[i] );
        diffsum += diffim[i];
      }
    }

    diffim.save ( argv[3], aol::fileNameEndsWith ( argv[3], ".png" ) ? qc::PNG_2D : ( aol::fileNameEndsWith ( argv[3], ".dat.bz2" ) ? qc::PGM_DOUBLE_BINARY : qc::PGM_UNSIGNED_CHAR_BINARY ) );

    cout << "In total, " << different_pixels << " pixels were different, that is "
    << ( 100.0 * different_pixels ) / diffim.size() << " percent." << endl;
    cout << "The total difference is " << diffsum << ", the average difference is "
    << ( 1.0 * diffsum ) / different_pixels << " among these, "
    << ( 1.0 * diffsum ) / diffim.size() << " for all pixels." << endl;

    return ( 0 );
  } else {
    cerr << "usage: " << endl << "compareImages image_1.pgm image_2.pgm output_difference_image.pgm" << endl;
    return ( 1 );
  }
}

