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
 * \brief Converts a 2D CIMG float or double data set (supports scalar and vector valued input data) to the quoc array format.
 *
 * Usage: convertCIMGToQuoc inputfile
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>

template <typename DataType>
void readRawAndSave ( const int NumX, const int NumY, const int NumZ, const int ImageDim, const string &BaseFileName, ifstream &In ) {
  qc::ScalarArray<DataType, qc::QC_2D> a ( NumX, NumY );
  for ( int slice = 0; slice < NumZ; ++slice ) {
    for ( int component = 0; component < ImageDim; ++component ) {
      aol::readBinaryData<DataType> ( In, a.getData(), NumX*NumY );
      string outputFileName = aol::strprintf ( "%s_%02d_%02d.bz2", BaseFileName.c_str(), slice, component );
      a.save ( outputFileName.c_str(), qc::SaveTypeTrait<DataType>::BinarySaveType );
    }
  }
}

enum RawFormat {
  RAW_FLOAT,
  RAW_DOUBLE
};

int main ( int argc, char **argv ) {

  try {
    string inputFileName;

    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile>" << endl;
      return EXIT_FAILURE;
    }
    if ( argc == 2 ) {
      inputFileName = argv[1];
    }
    cerr << "Reading CIMG from " << inputFileName << endl;

    ifstream in ( inputFileName.c_str(), ios::binary );
    aol::checkNextLineOrString ( in, "1", false );

    string temp;
    in >> temp;
    RawFormat format = RAW_FLOAT;
    if ( temp.compare ( "float" ) == 0 ) {
      if ( in.peek() == 32 ) {
        aol::checkNextLineOrString ( in, "little_endian", false );
      }
    } else {
      if ( temp.compare ( "double" ) == 0 ) {
        aol::checkNextLineOrString ( in, "little_endian", false );
        format = RAW_DOUBLE;
      } else
        throw aol::Exception ( "Unknown CIMG format", __FILE__, __LINE__ );
    }

    int numX, numY, numZ, imageDim;
    in >> numX >> numY >> numZ >> imageDim;
    // Eat the next char, it closes the header.
    in.get();

    cerr << numX << " " << numY << " "  << numZ << " "  << imageDim << endl;

    switch ( format ) {
      case RAW_FLOAT:
        readRawAndSave<float> ( numX, numY, numZ, imageDim, inputFileName, in );
        break;

      case RAW_DOUBLE:
        readRawAndSave<double> ( numX, numY, numZ, imageDim, inputFileName, in );
        break;

      default:
        throw aol::UnimplementedCodeException ( "Unsopported format", __FILE__, __LINE__ );
        break;
    };

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
