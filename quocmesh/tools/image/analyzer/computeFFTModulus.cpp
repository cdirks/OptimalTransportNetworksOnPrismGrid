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
 * \brief Computes the modulus of the FFT of an input image. The modulus is scaled with log and stretched to the PNG range.
 *
 * If p is specified, the lowest p*numPixles values of the modulus are mapped to black in the PNG.
 * The remaining values are still stretched to PNG range.
 *
 * If subtractMean is specified and true, the mean of the input data will be subtracted before the FFT is applied.
 *
 * \note The code is based on fft.cpp but too different from it to make joining both in one main file feasible.
 *
 * \author Berkels
 */

#include <convolution.h>

int main ( int argc, char* argv [] ) {
  try {
    if ( ( argc < 2 ) || ( argc > 4 ) ) {
      cerr << "Usage: " << argv [0] << " input [p] [subtractMean]\n";
      return EXIT_FAILURE;
    }

    const bool subtractMean = ( argc >= 4 ) ? ( atoi ( argv[3] ) != 0 ): false;

    qc::ScalarArray<double, qc::QC_2D> data ( argv [1] );
    qc::ScalarArray<double, qc::QC_2D> modulus ( data, aol::STRUCT_COPY );

    if ( subtractMean )
      data.addToAll( -data.getMeanValue() );

    qc::computeLogFFTModulus<double> ( data, modulus, ( argc >= 3 ) ? atof ( argv[2] ) : 0 );
    modulus.scaleValuesTo01();
    modulus.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
    modulus.savePNG ( ( string ( argv [1] ) + "_modulus.png" ).c_str()  );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
