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
 * \brief Converts a Brainlab VMR file to a 3D ScalarArray.
 *
 * Usage: convertVMRToQuoc InputFile
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

    if ( !aol::fileNameEndsWith ( inFileName.c_str(), ".vmr" ) )
      throw aol::FileFormatException ( "Filenames of VMR files have to end with \".vmr\"!", __FILE__, __LINE__ );

    const string outFileName = aol::strprintf ( "%s.pgm", inFileName.substr ( 0, inFileName.size() - 4 ).c_str() );

    std::ifstream in ( inFileName.c_str(), std::ios::binary );

    const int version = aol::readBinaryData<short, int> ( in );
    cerr << "VMR format " << version << endl;

    aol::Vec3<int> size;
    for ( int i = 0; i < 3; ++i )
      size[i] = aol::readBinaryData<short, int> ( in );

    cerr << "size = " << size << endl;

    qc::ScalarArray<unsigned short, qc::QC_3D> quocArray ( size );
    quocArray.loadRaw ( in, qc::PGM_UNSIGNED_CHAR_BINARY, size[0], size[1], size[2] );
    quocArray.save ( outFileName.c_str(), qc::PGM_UNSIGNED_CHAR_BINARY );
    in.close();
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
