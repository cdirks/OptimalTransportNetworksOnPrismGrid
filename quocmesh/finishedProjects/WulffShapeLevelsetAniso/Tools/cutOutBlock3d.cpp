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

#include <configurators.h>
/* ********************************************************
 * Datei: cutOutBlock3d.cpp
 * Autor: Oliver Nemitz
 * this program saves a cutout from a quocmesh 3d data set
 * in another scalarArray3d. The size and position (the latter
 * as well in the source as in the destination array) can be
 * determined.
 * ******************************************************** */

#include <quoc.h>
#include <aol.h>
#include <scalarArray.h>
#include <parameterParser.h>


#define RealType double

using namespace std;
using namespace aol::color;



void errorExit( const char *errorMsg )
{
  cerr << endl << red << errorMsg << reset << endl;
  exit(42);
}


int main(int argc, char** argv)
{
  try {

    if (argc <= 1) {
      cerr << "Cuts a block out of a scalarArray3d and saves it as a scalarArray3d.\n";
      cerr << "usage: " << argv [0] << " parameterfile" << endl;
      return 23;
    }

    aol::ParameterParser parser( argv[1] );

    // load the source image
    char imgName[ 1024 ];
    parser.getString( "loadName", imgName );
    cerr<<green<<"\nLoading file "<<red<<imgName<<green<<"...\n\n";
    qc::ScalarArray<double, qc::QC_3D> sourceImg( imgName );

    // Grösse des ScalarArrays einlesen
    int destLevel = parser.getInt("destLevel");
    qc::GridDefinition  grid (destLevel, qc::QC_3D);

    // define the scalarArray
    cerr<<blue<<"Declare memory for the scalarArray3d ... ";
    qc::ScalarArray<RealType, qc::QC_3D> destImg( grid );
    cerr<<"clear the memory...";
    destImg.setAll( 0. );

    // get the necessary parameters for copying the images into
    // the scalarArray
    int SourceSizeX   = parser.getInt("SourceSizeX");
    int SourceSizeY   = parser.getInt("SourceSizeY");
    int SourceSizeZ   = parser.getInt("SourceSizeZ");

    int SourceOffsetX = parser.getInt("SourceOffsetX");
    int SourceOffsetY = parser.getInt("SourceOffsetY");
    int SourceOffsetZ = parser.getInt("SourceOffsetZ");

    int DestOffsetX   = parser.getInt("DestOffsetX");
    int DestOffsetY   = parser.getInt("DestOffsetY");
    int DestOffsetZ   = parser.getInt("DestOffsetZ");

    int doThreshold   = parser.getInt("doThreshold");
    double thresholdLower  = sourceImg.getMinValue();
    double thresholdUpper  = sourceImg.getMaxValue();
    if (doThreshold) {
      thresholdLower  = parser.getDouble("thresholdLower");
      thresholdUpper  = parser.getDouble("thresholdUpper");
    }
    int saveOnly01    = parser.getInt("saveOnly01");
    int scaleTo0_1    = parser.getInt("scaleTo0_1");
    double offset     = parser.getDouble("offset");

    // check the size and offsets before copying
    if (  SourceSizeX+DestOffsetX > destImg.getNumX() ||
          SourceSizeY+DestOffsetY > destImg.getNumY() ||
          SourceSizeZ+DestOffsetZ > destImg.getNumZ()  )
      errorExit("DestImg is too small for the cutout!");



    // ---------------------------------------------------------------
    // now the main function: copy the block
    // ---------------------------------------------------------------

    cerr << blue << "\nCopying block:\n";
    RealType value = 0.;

    for (int i=0; i<SourceSizeX; i++) {
      cerr << ".";
      for (int j=0; j<SourceSizeY; j++) {
        for (int k=0; k<SourceSizeZ; k++) {

          value = sourceImg.get( i+SourceOffsetX, j+SourceOffsetY, k+SourceOffsetZ );

          if (doThreshold) { if ( value < thresholdLower || value > thresholdUpper ) value = 0.; }
          if ( saveOnly01 && value != 0 ) value = 1.;

          destImg.set( i+DestOffsetX, j+DestOffsetY, k+DestOffsetZ, value );

        }
      }
    }

    // -----------------------------------------------------------------


    // if desired, scale to [0,1]
    if (scaleTo0_1) {
      cerr << "scaling to [0,1]...";
      destImg -= destImg.getMinValue();
      destImg /= destImg.getMaxValue();
    }

    // adding the offset
    if (offset != 0) {
      cerr << "adding the offset "<<red<<offset<<blue<<"...";
      destImg += offset;
    }


    // save the scalarArray
    parser.getString( "saveName", imgName );
    cerr<<green<<"\nSaving the ScalarArray<QC_3D> '"<<red<<imgName<<green<<"'..."<<endl;
    destImg.save( imgName );

    cerr<<blue<<"\nReady, thank you for using this program and please pay the shareware-fee!\n"<<reset;

  }// try
  catch(aol::Exception e) {
    e.dump();
    cerr<<reset;
    return 42;
  }
}

