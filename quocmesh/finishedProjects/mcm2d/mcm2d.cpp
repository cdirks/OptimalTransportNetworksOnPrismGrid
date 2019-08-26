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

// -----------------------------------------------------------------------------
// mcm2d.cpp
// easy 2d mean curvature motion for 2^n+1 x 2^n+1 - pictures
// uses Marcs mcmSmooth-op
// -----------------------------------------------------------------------------

#include <quoc.h>
#include <qmException.h>
#include <aol.h>

#include <parameterParser.h>
#include <scalarArray.h>
#include <mcm.h>
#include <levelSetDrawer.h>
#include <configurators.h>

#include <signedDistanceOp.h>

#define REAL double
typedef qc::QuocConfiguratorTraitMultiLin<REAL, qc::QC_2D, aol::GaussQuadrature<REAL,qc::QC_2D,7> > ConfigType;
using namespace aol::color;


int main( int argc, char **argv )
{
  try {
    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>\n";
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    aol::ParameterParser parser( argv[1] );


    // ------------- LOAD IMAGE AND GET SAVENAME----------------
    char imgname[ 1024 ];
    parser.getString( "image", imgname );
    cerr<<green<<"\nLoading file "<<imgname<<"...\n\n";
    qc::ScalarArray<double, qc::QC_2D> img( imgname );
//     img.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
    cerr<<"\ndone."<<reset;

    char saveName[ 1024 ];
    parser.getString( "saveName", saveName );    // see below


    // -------------------- GRID AND TIMESTEP INFORMATION ---------------------------------
    int N = img.getNumX ();
    int d = qc::logBaseTwo (N);
    qc::GridDefinition grid( d, qc::QC_2D );
    qc::ScalarArray<double, qc::QC_2D> tmp( grid );
    qc::ScalarArray<double, qc::QC_2D> img_old( grid );

    // HACK: generate the signed distance function
    img /= img.getMaxValue();
    img.addToAll ( - 0.1 );
    qc::SignedDistanceOp<ConfigType> distOp( grid );
    distOp.apply( img, img );
    img.save("KurveSD.bz2",qc::PGM_DOUBLE_BINARY);

    REAL tau = 0.5 * aol::Sqr( parser.getDouble("tau") * grid.H() );
    REAL epsilon = parser.getDouble( "epsilon" ) *  grid.H();


    cerr << red << "declaring the op..." << reset;

    // -------------------- THE MCMSmoothOP -----------------------------------------------------
    qc::MCMSmoothOp< ConfigType > mcm(grid, tau, epsilon);

    // -------------------- MARKS A SPECIAL LEVELSET --------------------------------------------
    qc::LevelSetDrawer<ConfigType> drawer( grid );


    // ------------------- THE MAIN LOOP --------------------------------------------------
    img_old = img;

    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {

      // apply the operator
      cerr << blue << "Step: "<<iter << reset << ", applying op...\n";
      mcm.apply( img_old, img );
      img_old = img;
      cerr << "ready\n";

      // save every "timeOffset"th timestep
      if ( (iter % parser.getInt("timeOffset")) == 0 ) {
          char filename[ 1024 ];
          sprintf( filename, "%s_%05d.pgm", saveName, iter );     // was .bz2
          cerr << green << "\n>>>>>> saving to file " << filename << " <<<<<<\n\n";

          for (int i=0; i<img.size(); i++)
            if (img[i]>0) tmp[i] = 255.; else tmp[i] = 0.;

//           img.save( filename, 9 );
          tmp.save( filename, qc::PGM_UNSIGNED_CHAR_BINARY );
//           drawer.draw( img, tmp, 0.4 );
//           tmp.save( filename , 9 );   // ,9 for double valued pics
          cerr << "done!\n" << reset;
      }

    }


  } catch ( aol::Exception &ex ) {

    ex.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return EXIT_SUCCESS;
}











