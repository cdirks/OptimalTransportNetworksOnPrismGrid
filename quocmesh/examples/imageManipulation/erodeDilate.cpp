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
 * \brief First applies an erosion and then an dilation to the input image. This can be used to denoise characteristic functions.
 *
 * \author Berkels
 */

#include <scalarArray.h>

typedef double RType;
int main ( int argc, char **argv ) {
  try {
    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <radius>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const string outputFileNameBase = aol::getBaseFileName ( inputFileName );
    qc::ScalarArray<double, qc::QC_2D> a ( inputFileName );
    qc::ScalarArray<double, qc::QC_2D> tmp ( a, aol::STRUCT_COPY_UNINIT );

    const int radius = atoi ( argv[2] );

    // Erode
    for ( int y = 0; y < a.getNumY(); ++y ) {
      for ( int x = 0; x < a.getNumX(); ++x ) {
        const int XMin = aol::Max ( 0, x - radius );
        const int YMin = aol::Max ( 0, y - radius );
        const int XMax = aol::Min ( a.getNumX() - 1, x + radius );
        const int YMax = aol::Min ( a.getNumY() - 1, y + radius );

        qc::ScalarArray<double, qc::QC_2D> block ( ( XMax - XMin + 1 ), ( YMax - YMin + 1 ) );
        a.copyBlockTo( XMin, YMin, block );
        tmp.set ( x, y, block.getMinValue() );
      }
    }

    // Dilate
    for ( int y = 0; y < a.getNumY(); ++y ) {
      for ( int x = 0; x < a.getNumX(); ++x ) {
        const int XMin = aol::Max ( 0, x - radius );
        const int YMin = aol::Max ( 0, y - radius );
        const int XMax = aol::Min ( a.getNumX() - 1, x + radius );
        const int YMax = aol::Min ( a.getNumY() - 1, y + radius );

        qc::ScalarArray<double, qc::QC_2D> block ( ( XMax - XMin + 1 ), ( YMax - YMin + 1 ) );
        tmp.copyBlockTo( XMin, YMin, block );
        a.set ( x, y, block.getMaxValue() );
      }
    }

    tmp.savePNG ( aol::strprintf( "%s_%derode.png", outputFileNameBase.c_str(), radius ).c_str() );
    a.savePNG ( aol::strprintf( "%s_%derodeDilate.png", outputFileNameBase.c_str(), radius ).c_str() );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
