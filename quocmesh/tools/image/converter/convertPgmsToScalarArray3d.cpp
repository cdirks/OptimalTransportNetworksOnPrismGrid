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
 * \brief converts pgm files to ScalarArray
 *
 * Programm laedt eine Liste von pgm-Bildern ein
 * und speichert entweder alles oder einen Ausschnitt
 * in einem ScalarArray, dessen Groesse vorgegeben
 * werden kann!
 *
 * \author Nemitz
*/

#include <iostream>
#include <fstream>
#include <string.h>
#include <scalarArray.h>
#include <quoc.h>
#include <aol.h>
#include <parameterParser.h>

#define RealType double

using namespace aol::color;
using namespace std;


// auxiliary function: copies a (possible cutout of the) image into
// scalararray, using possible offsets
void copyImgToScalarArray ( qc::ScalarArray<RealType, qc::QC_3D>* data, qc::ScalarArray<RealType, qc::QC_2D> &image, int z,
                            bool doThreshold, RealType threshold, bool saveOnly01,
                            int SourceOffsetX, int SourceOffsetY,
                            int SourceSizeX,   int SourceSizeY,
                            int DestOffsetX,   int DestOffsetY,  int DestOffsetZ ) {
  cerr << "DestOffsetX: " << DestOffsetX;
  // now copy the rectangular cutout from the picture into the scalarArray
  RealType value = 1;            // this value has to be stored
  for ( int i = 0; i < SourceSizeX; i++ )
    for ( int j = 0; j < SourceSizeY; j++ ) {
      if ( saveOnly01 ) value = 1.;
      else value = image.get ( SourceOffsetX + i, SourceOffsetY + j );

      if ( doThreshold ) {
        if ( image.get ( SourceOffsetX + i, SourceOffsetY + j ) < threshold )
          data->set ( DestOffsetX + i, DestOffsetY + j, DestOffsetZ + z, 0. );
        else  data->set ( DestOffsetX + i, DestOffsetY + j, DestOffsetZ + z, value );
        //cerr<<value<<", "; }
      } else data->set ( DestOffsetX + i, DestOffsetY + j, DestOffsetZ + z, value );
    }
}



int main ( int argc, char** argv ) {
  try {

    if ( argc <= 2 ) {
      cerr << "Converts a lot of single pgm-pics to one scalarArray3d.\n";
      cerr << "usage: " << argv [0] << " parameterfile picture(s)" << endl;
      return 23;
    }

    aol::ParameterParser parser ( argv[1] );

    // Groesse des ScalarArrays einlesen
    int DestSizeX = parser.getInt ( "DestSizeX" );
    int DestSizeY = parser.getInt ( "DestSizeY" );
    int DestSizeZ = parser.getInt ( "DestSizeZ" );

//     if (argc-2 > DestSizeZ) {
//       cerr << Red << "\nScalarArray is too small for such a lot of pictures\n" << black;
//       return 23;
//     }

    // define the scalarArray
    cerr << blue << "Declare memory for the scalarArray3d ... ";
    qc::ScalarArray<RealType, qc::QC_3D>* data = new qc::ScalarArray<RealType, qc::QC_3D> ( DestSizeX, DestSizeY, DestSizeZ );
    cerr << "clear the memory...";
    data->setAll ( 0. );

    // get the necessary parameters for copying the images into
    // the scalarArray
    int SourceSizeX   = parser.getInt ( "SourceSizeX" );
    int SourceSizeY   = parser.getInt ( "SourceSizeY" );

    int SourceOffsetX = parser.getInt ( "SourceOffsetX" );
    int SourceOffsetY = parser.getInt ( "SourceOffsetY" );
    int SourceOffsetZ = parser.getInt ( "SourceOffsetZ" );

    int DestOffsetX   = parser.getInt ( "DestOffsetX" );
    int DestOffsetY   = parser.getInt ( "DestOffsetY" );
    int DestOffsetZ   = parser.getInt ( "DestOffsetZ" );

    bool doThreshold   = !!parser.getInt ( "doThreshold" );
    RealType threshold  = parser.getDouble ( "threshold" );
    bool saveOnly01    = !!parser.getInt ( "saveOnly01" );


    if ( argc - 2 - SourceOffsetZ > DestSizeZ ) {
      cerr << red << "\nScalarArray is too small for such a lot of pictures\n";
      cerr << "=> not all pictures are used!" << black;
      argc = 2 + SourceOffsetZ + DestSizeZ;
    }

    // load the first image and get hereby the size of them
    cerr << green << " ready. \nLoading image number " << red << 1 + SourceOffsetZ << endl << black;
    qc::ScalarArray<RealType, qc::QC_2D> image ( argv[2+SourceOffsetZ] );
    copyImgToScalarArray ( data, image, 0, doThreshold, threshold, saveOnly01,
                           SourceOffsetX, SourceOffsetY, SourceSizeX, SourceSizeY,
                           DestOffsetX, DestOffsetY, DestOffsetZ );

    // now load all images
    for ( int count = 3 + SourceOffsetZ; count < argc; count++ ) {
      cerr << green << "Loading image number " << red << count - 1 << endl << black;
      image.load ( argv[count] );
      copyImgToScalarArray ( data, image, count - 2 - SourceOffsetZ, doThreshold, threshold, saveOnly01,
                             SourceOffsetX, SourceOffsetY, SourceSizeX, SourceSizeY,
                             DestOffsetX, DestOffsetY, DestOffsetZ );
    }

    // divide the scalarArray by its bigges value
    RealType maxVal = data->getMaxValue();
    cerr << blue << "Scaling to [0,1] (maxValue = " << maxVal << ")...";
    ( *data ) /= maxVal;

    RealType offset = parser.getDouble ( "offset" );
    if ( offset ) {
      cerr << "adding the offset " << offset << "...";
      ( *data ).addToAll ( offset );
    }


    // save the scalarArray
    char saveName[ 1024 ];
    parser.getString ( "saveName", saveName );
    cerr << blue << "\nSaving the ScalarArray<QC_3D> '" << green << saveName << blue << "'...";
    data->save ( saveName, qc::SaveTypeTrait<RealType>::BinarySaveType );

    cerr << "\nReady, thank you for using this program and please pay the shareware-fee!\n" << black;

  }// try
  catch ( aol::Exception e ) {
    e.dump();
    cerr << black;
    return 42;
  }
}

