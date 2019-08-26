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
 * \brief Averages multiple 2D ScalarArrays and saves the result with double precision.
 *
 * Usage: averageImages InputFile1 ... InputFileN
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>

int main ( int argc, char **argv ) {

  try {
    string outFileName;

    if ( argc < 2 ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile1> ... <InputFileN>" << endl;
      return EXIT_FAILURE;
    }

    qc::ScalarArray<double, qc::QC_2D>  averageArray ( argv[1] );

    const int numberOfQuocArrays = argc - 1;

    for ( int i = 1; i < numberOfQuocArrays; ++i ) {
      qc::ScalarArray<double, qc::QC_2D>  tempArray ( argv[1+i] );
      averageArray += tempArray;
    }

    averageArray /= numberOfQuocArrays;
    averageArray.save ( "average.dat.bz2", qc::PGM_DOUBLE_BINARY );

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
