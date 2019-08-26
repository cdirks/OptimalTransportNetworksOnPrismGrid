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
 * \brief Converts a STM (simple terrain map) file with unsigned short precision to a 2D ScalarArray.
 *
 * Usage: convertQuocToSTM InputFile
 *
 * \author Berkels
 */

#include <scalarArray.h>

int main ( int argc, char **argv ) {

  try {
    if ( argc < 2 ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile>" << endl;
      return EXIT_FAILURE;
    }

    const string inFileName = argv[1];
    const string outFileName = aol::strprintf ( "%s.dat.bz2", inFileName.c_str() );

    const unsigned long MAGIC_VALUE = 0x04030201;
    union { unsigned int val;  char byte[4]; } test;

    std::ifstream in ( inFileName.c_str(), std::ios::binary );
    aol::checkNextLineOrString ( in, "STM", false );
    int width, height;
    in >> width;
    in >> height;
    in.ignore();
    for ( int i = 0; i < 4; ++i )
      test.byte[i] = in.get();
    in.ignore();

    qc::ScalarArray<unsigned short, qc::QC_2D> stmArray ( width, height );
    stmArray.loadRaw ( in, qc::PGM_UNSIGNED_SHORT_BINARY, width, height );
    in.close();

    if ( test.val != MAGIC_VALUE )
      stmArray.swapByteOrder();

    qc::ScalarArray<unsigned short, qc::QC_2D> quocArray ( stmArray, aol::STRUCT_COPY );
    // STM stores the rows to be in reverse order.
    quocArray.flipFrom ( stmArray, qc::QC_Y );
    quocArray.save ( outFileName.c_str(), qc::PGM_DOUBLE_BINARY );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
