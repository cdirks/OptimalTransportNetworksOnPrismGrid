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
 * \brief Converts a 3D ScalarArray to the meta image data format (compatible with ParaView) with unsigned short precision.
 *
 * Usage: convertQuoc3dToMHD InputFile
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
    qc::ScalarArray<double, qc::QC_3D> quocArray ( inFileName.c_str() );
    quocArray.saveMetaImageFile ( inFileName.c_str(), qc::PGM_UNSIGNED_SHORT_BINARY );

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
