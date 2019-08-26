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
 *  \brief stretching the contrast in a png image
 *
 *  Stretching the contrast in a png image (8 bit/channel) from given
 *  range to [0,255] in each color channel.
 *
 *  Usage: colorContrastStretch infile outfile minvalue maxvalue
 *
 *  \author Schwen
 */


#include <aol.h>
#include <FEOpInterface.h>
#include <multiArray.h>

typedef double RealType;
int main ( int argc, char **argv ) {
  try {

    char inputFileName[1024];
    char outputFileName[1024];

    if ( argc != 5 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <output_file> <min> <max>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    sprintf ( inputFileName, "%s",  argv[1] );
    sprintf ( outputFileName, "%s",  argv[2] );
    const RealType min = atof ( argv[3] );
    const RealType max = atof ( argv[4] );

    cerr << "Reading from " << inputFileName << endl;

    qc::MultiArray< unsigned char, qc::QC_2D, 3 > a;
    a.loadPNG ( inputFileName );


    for ( short comp = 0; comp < 3; ++comp ) {
      for ( int i = 0; i < a[comp].size(); ++i ) {
        const RealType dwarf = aol::Clamp<RealType> ( a[comp][i], min, max ) - min;
        a[comp][i] = static_cast<unsigned char> ( 255 * dwarf / static_cast<RealType> ( max - min ) );
      }
    }

    a.savePNG ( outputFileName );

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
