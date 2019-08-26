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
 * \brief Adds equally distributed noise to an image.
 *
 * \author Berkels
 */

#include <aol.h>
#include <op.h>
#include <FEOpInterface.h>
#include <scalarArray.h>

typedef double RType;
int main ( int argc, char **argv ) {

  try {
    char inputFileName[1024];
    char outputFileName[1024];
    RType variance = 0.;

    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <variance>" << endl;
      return EXIT_FAILURE;
    }
    if ( argc == 3 ) {
      sprintf ( inputFileName, "%s",  argv[1] );
      variance = atof ( argv[2] );
    }
    cerr << "Reading ScalarArray from " << inputFileName << endl;

    qc::ScalarArray<RType, qc::QC_2D> a ( inputFileName );
    a /= a.getMaxValue();

    cerr << "Mean value before adding noise : " << a.getMeanValue() << endl;

    aol::NoiseOperator<RType> noiseOP ( aol::NoiseOperator<RType>::EQUALLY_DISTRIBUTED, -variance, variance );
    noiseOP.applyAdd ( a, a );

    cerr << "Mean value after adding noise : " << a.getMeanValue() << endl;

    sprintf ( outputFileName, "%s_noise_%.1f.dat.bz2", inputFileName, variance );

    a.setQuietMode ( true );
    a.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
    a.save ( outputFileName, qc::PGM_DOUBLE_BINARY );

#ifdef __MINGW32_VERSION
    system ( "PAUSE" );
#endif
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  return 1;
}
