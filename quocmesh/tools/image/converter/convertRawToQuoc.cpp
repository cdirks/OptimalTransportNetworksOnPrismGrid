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
 * \brief Converts 2D or 3D raw data to a ScalarArray.
 *
 * Usage: convertRawToQuoc inputfile numX numY numZ InType OutType
 *
 * A negative numZ value will try to load 3D raw input data from 2D raw slice files.
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>

int main ( int argc, char **argv ) {

  try {
    if ( ( argc < 6 ) || ( argc > 9 ) ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile> <numX> <numY> <numZ> <InType> [OutType] [GuessHeaderSize] [SwapByteOrder]" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    string inputFileName = argv[1];
    const string outputFileName = ( ( aol::fileNameEndsWith ( inputFileName.c_str(), ".raw" ) ? inputFileName.substr ( 0, inputFileName.size() - 4 ) : inputFileName ) + ".dat.bz2" );
    const int numX = atoi ( argv[2] );
    const int numY = atoi ( argv[3] );
    const int numZ = atoi ( argv[4] );
    const qc::SaveType inType = static_cast <qc::SaveType> ( atoi ( argv[5] ) );
    const qc::SaveType outType = ( argc > 6 ) ? static_cast <qc::SaveType> ( atoi ( argv[6] ) ) : inType;
    const bool guessHeaderSize = ( argc > 7 ) ? !!( atoi ( argv[7] ) ) : false;
    const bool swapByteOrder = ( argc > 8 ) ? !!( atoi ( argv[8] ) ) : false;

    if ( numZ >= 0 ) {
      cerr << "Reading raw data from " << inputFileName;

      ifstream in ( inputFileName.c_str(), ios::binary );

      if ( guessHeaderSize ) {
        const int expectedDataSize = numX * numY * aol::Max ( numZ, 1 ) * qc::getSizeOfSaveType ( inType );
        if ( expectedDataSize < aol::getSizeOfFile ( inputFileName ) )
          in.seekg( -expectedDataSize, ios::end );
      }

      // Note: double should have enough precision to store all formats supported by loadRaw.
      if ( numZ > 0 ) {
        qc::ScalarArray<double, qc::QC_3D> dataArray ( numX, numY, numZ );
        dataArray.loadRaw ( in, inType, numX, numY, numZ );
        if ( swapByteOrder )
          dataArray.swapByteOrder();
        cerr << " done.\nSaving quoc array";
        dataArray.save ( outputFileName.c_str(), outType );
      }
      else {
        qc::ScalarArray<double, qc::QC_2D> dataArray ( numX, numY );
        dataArray.loadRaw ( in, inType, numX, numY );
        if ( swapByteOrder )
          dataArray.swapByteOrder();
        cerr << " done.\nSaving quoc array";
        dataArray.save ( outputFileName.c_str(), outType );
      }
      cerr << " done." << endl;
    }
    else {
      const int numSlices = -numZ;
      qc::ScalarArray<double, qc::QC_3D> dataArray ( numX, numY, numSlices );
      qc::ScalarArray<double, qc::QC_2D> sliceArray ( numX, numY );
      for ( int i = 0; i < numSlices; ++i ) {
        string sliceFileName = aol::strprintf ( inputFileName.c_str(), i );
        ifstream in ( sliceFileName.c_str(), ios::binary );
        cerr << "Reading raw data from " << sliceFileName << endl;
        sliceArray.loadRaw ( in, inType, numX, numY );
        dataArray.putSlice ( qc::QC_Z, i, sliceArray );
      }
      if ( swapByteOrder )
        dataArray.swapByteOrder();
      dataArray.save ( outputFileName.c_str(), outType );
    }
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
