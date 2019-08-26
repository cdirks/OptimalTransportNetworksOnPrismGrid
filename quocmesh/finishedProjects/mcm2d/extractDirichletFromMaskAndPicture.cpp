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

// extractDirichletFromMaskAndPicture.cpp
// as the name says, this programs extracts Dirichlet boundary values
// from the original picture and a 0-1-mask!
// -----------------------------------------------------------------------------

#include <quoc.h>
#include <qmException.h>
#include <aol.h>

#include <parameterParser.h>
#include <scalarArray.h>

#define REAL double
using namespace aol::color;
using namespace std;

double uniformrandom (double alpha = 0, double beta = 1) { return (static_cast<double> (rand ()) / (RAND_MAX-1)) * (beta - alpha) + alpha; }


int main( int argc, char **argv )
{
  try {
    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>\n";
      return EXIT_FAILURE;
    }

    aol::ParameterParser parser( argv[1] );


    // ------------- LOAD IMAGE aND MASK AND GET SAVENAME----------------
    char filename[ 1024 ];
    parser.getString( "image", filename );
    cerr<<green<<"\nLoading image '"<<blue<<filename<<green<<"'...\n\n";
    qc::ScalarArray<double, qc::QC_2D> image( filename );

    parser.getString( "mask", filename );
    cerr<<green<<"\nLoading mask '"<<blue<<filename<<green<<"'...\n\n";
    qc::ScalarArray<double, qc::QC_2D> mask( filename );

//     img.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
    cerr<<"\ndone."<<reset;



    // -------------------- GRID AND BACKGROUND INFORMATION ---------------------------------
    int d = parser.getInt("level");
    int N = parser.getInt("gridSize");
    qc::GridDefinition grid( d, qc::QC_2D );
    qc::ScalarArray<double, qc::QC_2D> startImg( grid );
    qc::ScalarArray<double, qc::QC_2D> boundaryMask( grid );

    // set background of the image
    if ( parser.getInt("backgroundType") == 0 )   // constant value
      startImg.setAll( parser.getDouble("backgroundValue") );
    if ( parser.getInt("backgroundType") == 1 )   // white noise
    {
      for (int i=0; i<startImg.size(); i++)
        startImg[i] = uniformrandom( 0., 255. );
    }

    boundaryMask.setAll( 0. );



    // ------------------- THE LOOP OVER THE DIRICHLET-DATA -------------------------
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
        if ( mask.get(i,j) != 0 ) {
          startImg.set( i,j, image.get(i,j) );
          boundaryMask.set( i,j, 255. );
        }


    // --------------- SAVE IT ---------------------------
    parser.getString( "startImgSaveName", filename );
    cerr << green << "\n>>>>>> saving to file '" <<blue<< filename << green << "' <<<<<<\n\n";
    startImg.save( filename, qc::PGM_UNSIGNED_CHAR_BINARY );

    parser.getString( "boundaryMaskSaveName", filename );
    cerr << green << "\n>>>>>> saving to file '" <<blue<< filename << green << "' <<<<<<\n\n";
    boundaryMask.save( filename, qc::PGM_UNSIGNED_CHAR_BINARY );
    cerr << "done!\n" << reset;


  } catch ( aol::Exception &ex ) {

    ex.dump();
  }
  return EXIT_SUCCESS;
}

