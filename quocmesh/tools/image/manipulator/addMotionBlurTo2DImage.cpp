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
 * \brief Adds motion blur to a 2D image by convolving it with the motion blur
 *        kernel corresponding to the selected velocity.
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>
#include <convolution.h>

typedef double RType;
#ifdef USE_LIB_FFTW

int main ( int argc, char **argv ) {
  try {
    char inputFileName[1024];
    aol::Vec2<double> v;

    if ( argc != 4 ) {
      cerr << "USAGE: " << argv[0] << " <input_file> <x_veclocity> <y_velocity>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }
    if ( argc == 4 ) {
      sprintf ( inputFileName, "%s",  argv[1] );
      v[0] = atof ( argv[2] );
      v[1] = atof ( argv[3] );
    }

    qc::ScalarArray<double, qc::QC_2D> a ( inputFileName );
    qc::addMotionBlurToArray ( v, a, a );

    char outputFileBaseName[1024];
    sprintf ( outputFileBaseName, "%sMotionBlur%.1f-%.1f", inputFileName, v[0], v[1] );
    string outFileNameString = outputFileBaseName;
    outFileNameString += ".pgm";
    a.save ( outFileNameString.c_str(), qc::PGM_UNSIGNED_CHAR_BINARY );

    outFileNameString = outputFileBaseName;
    outFileNameString += ".dat.bz2";
    a.save ( outFileNameString.c_str(), qc::PGM_DOUBLE_BINARY );
    aol::callSystemPauseIfNecessaryOnPlatform();
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  return 1;
}

#else

int main ( int, char** ) {
  cerr << "Needs to be compiled with -DUSE_LIB_FFTW\n";
  return 0;
}

#endif
