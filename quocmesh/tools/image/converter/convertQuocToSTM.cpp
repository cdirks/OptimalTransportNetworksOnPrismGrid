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
 * \brief Converts a 2D ScalarArray to the STM (simple terrain map) format with unsigned short precision.
 *
 * Usage: convertQuocToSTM InputFile
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>

int main ( int argc, char **argv ) {

  try {
    if ( argc < 2 ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile>" << endl;
      return EXIT_FAILURE;
    }

    const string inFileName = argv[1];
    const string outFileName = aol::strprintf ( "%s.stm", inFileName.c_str() );

    qc::ScalarArray<float, qc::QC_2D> quocArrayOrig ( inFileName.c_str() );
    qc::ScalarArray<float, qc::QC_2D> quocArray ( quocArrayOrig, aol::STRUCT_COPY );
    // STM expects the rows to be in reverse order.
    quocArray.flipFrom ( quocArrayOrig, qc::QC_Y );

    ofstream out ( outFileName.c_str(), ios::binary );

    const unsigned long MAGIC_VALUE = 0x04030201;
    union { unsigned int val; unsigned char byte[4]; } test;
    test.val = MAGIC_VALUE ;

    out << aol::strprintf ( "STM %d %d %c%c%c%c\t", quocArray.getNumX(), quocArray.getNumY(), test.byte[0], test.byte[1], test.byte[2], test.byte[3]);
    quocArray.saveRaw ( out, qc::PGM_UNSIGNED_SHORT_BINARY, quocArray.getMinValue(), quocArray.getMaxValue() );

    out.close();

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
