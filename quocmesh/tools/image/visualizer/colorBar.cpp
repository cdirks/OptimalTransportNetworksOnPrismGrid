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
 * \brief Visualizes any aol::colorTrans, e.g. aol::HSV_BLUE_TO_RED, as PNG.
 *
 * \author Berkels
 */

#include "colorWheel.h"

typedef double RType;

int main ( int argc, char **argv ) {
  try {

    if ( ( argc != 2 ) && ( argc != 3 ) ){
      cerr << "USAGE: " << argv[0] << " numX [ColorTrans]" << endl;
      return EXIT_FAILURE;
    }

    const int numX = atoi ( argv[1] );
    const aol::colorTrans colorTrans = ( argc > 2 ) ? static_cast<aol::colorTrans> ( atoi ( argv[2] ) ) : aol::HSV_BLUE_TO_RED;
    const aol::RGBColorMap<RType> hsvMap ( 0., 1., colorTrans );
    qc::MultiArray<RType, 2, 3> tempMArray ( numX, 1 );
    tempMArray.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );

    for ( int i = 0; i < numX; ++i ) {
      aol::Vec3< RType > color;
      hsvMap.scalarToColor ( static_cast<RType> ( i ) / (numX-1), color );
      for ( int j = 0; j < 3; ++j )
        tempMArray[j][i] = color[j];
    }
    tempMArray.savePNG ( "colorBar.png" );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
