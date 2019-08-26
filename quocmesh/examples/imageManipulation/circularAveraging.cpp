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
 * \brief Applies the circular averaging filter with different kernel sizes to a 2D ScalarArray.
 *        Also saves the corresponding kernels as PNGs.
 *
 * \author Berkels
 */

#include <scalarArray.h>

typedef double RType;
int main ( int argc, char **argv ) {
  try {
    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const string outputFileNameBase = aol::getBaseFileName ( inputFileName );
    qc::ScalarArray<double, qc::QC_2D> a ( inputFileName );

    for ( int i = 1; i < 10; ++++i ) {
      qc::CircleAverageKernel2d<RType> kernel ( i );
      qc::ScalarArray<double, qc::QC_2D> b ( a, aol::STRUCT_COPY_UNINIT );
      a.applyLinearFilterTo ( kernel, b );
      b.savePNG ( aol::strprintf( "%s_%d.png", outputFileNameBase.c_str(), i ).c_str() );

      const int offset = kernel.getOffset();
      qc::ScalarArray<double, qc::QC_2D> c ( kernel.getSize(), kernel.getSize() );
      for ( int y = -offset; y <= offset; ++y )
        for ( int x = -offset; x <= offset; ++x )
          c.set ( x + offset, y + offset, kernel.getValue( x, y ) );

      c.setOverflowHandlingToCurrentValueRange();
      c.savePNG ( aol::strprintf( "%s_%d_kernel.png", outputFileNameBase.c_str(), i ).c_str() );
    }
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
