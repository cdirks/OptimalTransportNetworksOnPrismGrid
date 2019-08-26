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
 * \brief Reads a 3D quoc array and saves it as 2D slices.
 *
 * \author Berkels
 */

#include <aol.h>
#include <FEOpInterface.h>
#include <scalarArray.h>

typedef float RType;
int main ( int argc, char **argv ) {

  try {
    if ( ( argc < 2 ) || ( argc > 5 ) ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> [dir] [normalize] [OutType]" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];

    cerr << "Reading 3D PGM from " << inputFileName << endl;

    const qc::Comp direction = ( argc >= 3 ) ? static_cast<qc::Comp> ( atoi ( argv[2] ) ) : qc::QC_X;
    const bool normalize = ( argc >= 4 ) ? static_cast<bool> ( atoi ( argv[3] ) ) : true;
    const qc::SaveType outType = ( argc >= 5 ) ? static_cast <qc::SaveType> ( atoi ( argv[4] ) ) : qc::PNG_2D;

    qc::ScalarArray<RType, qc::QC_3D> a ( inputFileName.c_str() );

    if ( normalize )
      a.scaleValuesTo01 ();

    string outbasename;
    if ( aol::fileNameEndsWith ( inputFileName.c_str(), ".dat.bz2" ) )
      outbasename = aol::getBaseFileName ( inputFileName.substr ( 0, inputFileName.size() - 4 ) );
    else
      outbasename = aol::getBaseFileName ( inputFileName );

    const string outputFileName = aol::strprintf ( "%s_%%03d%s", outbasename.c_str(), qc::getDefaulSuffixOfSaveType ( outType ) );
    a.saveSlices ( outputFileName.c_str(), direction, outType, NULL, aol::CLIP_THEN_SCALE, 0, 1 );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
