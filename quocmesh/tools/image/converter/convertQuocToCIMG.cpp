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
 * \brief Converts (multiple) 2D ScalarArrays to the CIMG format with float precision.
 *
 * Usage: convertQuocToCIMG InputFile1 ... InputFileN
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

    const int numberOfQuocArrays = argc - 1;
    qc::ArrayHeader header;
    {
      // Make sure that in is destroyed by this additional scope, so that the file argv[1]
      // can be opened again later.
      aol::Bzipifstream in ( argv[1] );
      qc::ReadArrayHeader ( in, header );
    }
    aol::MultiVector<float> quocVectors ( numberOfQuocArrays, header.numX * header.numY * header.numZ );

    for ( int i = 0; i < numberOfQuocArrays; ++i ) {
      qc::ScalarArray<float, qc::QC_2D>  tempArray ( quocVectors[i], header.numX, header.numY, aol::FLAT_COPY );
      tempArray.load ( argv[1+i] );
    }

    outFileName = aol::strprintf ( "%s.cimg", argv[1] );

    ofstream out ( outFileName.c_str(), ios::binary );
    out << "1 float\n";
    out << header.numX << " " << header.numY << " " << header.numZ << " " << numberOfQuocArrays << endl;


    for ( int i = 0; i < numberOfQuocArrays; ++i )
      out.write ( reinterpret_cast<const char*> ( quocVectors[i].getData() ), header.numX * header.numY * header.numZ * sizeof ( float ) );

    out.close();

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
