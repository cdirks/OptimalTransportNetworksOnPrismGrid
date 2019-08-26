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

#include <multiArray.h>
/*
 * author: Effland
 *
 */

int main ( int argc, char **argv ) {

  try {
    if ( argc != 6 ) {
      cerr << "USAGE: " << argv[0] << " <input_file> <cropStartX> <cropStartY> <cropSizeX> <cropSizeY>" << endl;
      return ( EXIT_FAILURE);
    }

    const string inputFileName = argv[1];
    const string outputFileName = aol::getBaseFileName ( inputFileName ) + "_cropped.png";
    const aol::Vec2 < int > cropStart ( atoi ( argv[2] ), atoi ( argv[3] ) );
    const aol::Vec2 < int > cropSize ( atoi ( argv[4] ), atoi ( argv[5] ) );
    qc::MultiArray < double, 2, 3 > a ( inputFileName.c_str () );
    qc::MultiArray < double, 2, 3 > out ( cropSize.x (), cropSize.y () );
    a.copyBlockTo ( cropStart, out );
    out.savePNG ( outputFileName.c_str () );
  }
  catch ( aol::Exception &el ) {
    el.dump ();
  }
  aol::callSystemPauseIfNecessaryOnPlatform ();
  return ( EXIT_SUCCESS);
}
