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
 * \brief Loads 3D quoc dataset, normalizes it and saves it as 8 bit 3D PGM.
 *
 * \author Berkels
 */

#include <aol.h>
#include <FEOpInterface.h>
#include <scalarArray.h>

typedef float RType;
int main ( int argc, char **argv ) {

  try {
    char inputFileName[1024];
    char outputFileName[1024];

    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file>" << endl;
      return EXIT_FAILURE;
    }
    if ( argc == 2 ) {
      sprintf ( inputFileName, "%s",  argv[1] );
    }
    cerr << "Reading 3D PGM from " << inputFileName << endl;

    qc::ScalarArray<RType, qc::QC_3D> a ( inputFileName );
    a /= a.getMaxValue();
    sprintf ( outputFileName, "%s_8.raw", inputFileName );

    /*
        qc::ScalarArray<RType, qc::QC_3D> b( 257, 257, 257 );
        b.scale( a );
        b.save( outputFileName, qc::PGM_UNSIGNED_CHAR_BINARY);
    */
    a.save ( outputFileName, qc::PGM_UNSIGNED_CHAR_BINARY );

#ifdef __MINGW32_VERSION
    system ( "PAUSE" );
#endif
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  return 1;
}
