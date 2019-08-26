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
 *  \brief linear interpolation between two images
 *
 *  Merges two given images according to
 *  destinationImage = t*sourceImage1 + (1-t)*sourceImage2.
 *
 *  Usage: mergeImages sourceImage1 sourceImage2 destinationImage t
 *
 *  \author Wirth
 */

#include <iostream>
#include "scalarArray.h"

// get the first character of the header
char getSaveType ( const char* FileName ) {
  aol::ipfstream file ( FileName );
  qc::ArrayHeader header;
  qc::ReadArrayHeader ( file, header );
  file.close();

  return header.magic[1];
}

static const qc::Dimension Dim = qc::QC_2D;
typedef double RealType;

// compare two pgm images, compute number of different pixels and sum of the absolute of the difference, save difference image.

int main ( int argc, char **argv ) {
  if ( argc == 5 ) {
    qc::ScalarArray<RealType,Dim> image1 ( argv[1] ), image2 ( argv[2] ), mergeIm ( image1, aol::DEEP_COPY );
    char magicNum[1];
    magicNum[0] = getSaveType( argv[1] );
    qc::SaveType type = static_cast<qc::SaveType>( atoi( magicNum ) );
    RealType t = strtod( argv[4], NULL );

    mergeIm *= t;
    image2 *= 1. - t;
    mergeIm += image2;

    mergeIm.save( argv[3], type );

    return ( 0 );
  } else {
    cerr << aol::color::red << "USAGE: " << argv[0] << " sourceImage1 sourceImage2 destinationImage t" << endl << aol::color::reset;
    return EXIT_FAILURE;
  }
}
