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
 * \brief Reads a 3D rectilinear VTK data file and saves it in the quoc array format.
 *
 * \author Berkels
 */

#include <aol.h>
#include <multiArray.h>

void readCoordinates ( ifstream &In, aol::Vector<float> &Coords, const int Num, const char *VTKCoordName ) {
  aol::checkNextLineOrString ( In, VTKCoordName, false );
  int localNum;
  In >> localNum;
  if ( Num != localNum )
    throw aol::Exception ( aol::strprintf ( "%s mismatch: expected %d, got %d\n", VTKCoordName, Num, localNum ).c_str(), __FILE__, __LINE__ );

  aol::checkNextLineOrString ( In, " float" );

  Coords.reallocate( Num );
  aol::readBinaryData<float> ( In, Coords.getData(), Num );
}

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

    string outputFileNameBase = aol::getFileName ( inputFileName.c_str() );
    // If the file name ends with ".vtk", remove that suffix.
    if ( aol::fileNameEndsWith ( outputFileNameBase.c_str(), ".vtk" ) )
      outputFileNameBase = outputFileNameBase.substr ( 0, outputFileNameBase.length() - 4 );

    cerr << "Reading VTK rectilinear data from " << inputFileName << endl;

    ifstream in ( inputFileName.c_str(), ios::binary );

    aol::checkNextLineOrString ( in, "# vtk DataFile Version 2.0" );
    aol::checkNextLineOrString ( in, "nav" );
    aol::checkNextLineOrString ( in, "BINARY" );
    aol::checkNextLineOrString ( in, "DATASET RECTILINEAR_GRID" );
    aol::checkNextLineOrString ( in, "DIMENSIONS", false );

    int numX, numY, numZ;
    in >> numX >> numY >> numZ;

    cerr << "Input size is " << numX << "x" << numY << "x"  << numZ << endl;

    aol::Vector<float> a;
    readCoordinates ( in, a, numX, "X_COORDINATE" );
    readCoordinates ( in, a, numY, "Y_COORDINATE" );
    readCoordinates ( in, a, numZ, "Z_COORDINATE" );

    aol::checkNextLineOrString ( in, "POINT_DATA", false );
    int numPoints;
    in >> numPoints;
    // Read the next char and make sure that it is a new line.
    if ( in.get() != '\n' )
      throw aol::Exception ( "Parsing error: Didn't get expected '\\n'.\n", __FILE__, __LINE__ );

    if ( numPoints != numX*numY*numZ )
      throw aol::Exception ( aol::strprintf ( "POINT_DATA mismatch: expected %d, got %d\n", numX*numY*numZ, numPoints ).c_str(), __FILE__, __LINE__ );

    string temp;
    getline ( in, temp );

    if ( temp.compare ( "SCALARS l float" ) == 0 ) {
      aol::checkNextLineOrString ( in, "lookup_table default" );

      qc::ScalarArray<float, qc::QC_3D> dataArray ( numX, numY, numZ );
      aol::readBinaryData<float> ( in, dataArray.getData(), numPoints );
      if ( in.fail() )
        throw aol::Exception ( "Failed to read from the input stream.\n", __FILE__, __LINE__ );
      // The byte order in the VTK file seems to be different from ours, so swap it after reading the binary values.
      dataArray.swapByteOrder();
      dataArray.save ( aol::strprintf( "%s.bz2", outputFileNameBase.c_str() ).c_str(), qc::SaveTypeTrait<float>::BinarySaveType );
    }
    else if ( temp.compare ( "VECTORS vectors float" ) == 0 ) {
      a.reallocate ( 3 * numPoints );
      aol::readBinaryData<float> ( in, a.getData(), 3*numPoints );
      if ( in.fail() )
        throw aol::Exception ( "Failed to read from the input stream.\n", __FILE__, __LINE__ );
      a.swapByteOrder();
      qc::MultiArray<float, qc::QC_3D> dataMultiArray ( qc::GridSize<qc::QC_3D> ( numX, numY, numZ ) );
      for ( int i = 0; i < 3 * numPoints; ++i )
        dataMultiArray[i%3][i/3] = a[i];
      dataMultiArray.save ( aol::strprintf( "%s_%%d.bz2", outputFileNameBase.c_str() ).c_str(), qc::SaveTypeTrait<float>::BinarySaveType );
    }
    else {
      string errorMessage = aol::strprintf ( "Unknown header \"%s\".\n", temp.c_str() );
      throw aol::Exception ( errorMessage, __FILE__, __LINE__ );
    }
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
