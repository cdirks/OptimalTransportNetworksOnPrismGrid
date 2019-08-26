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
// Reconstruct.cpp
// should reconstruct pictures originated by a few parts with non-zero
// Gauss-curvature. Therefore a projective CG solver is used.
// -----------------------------------------------------------------------------

#include <quoc.h>
#include <qmException.h>
#include <aol.h>

#include <parameterParser.h>
#include <scalarArray.h>
#include <mcm.h>
#include <levelSetDrawer.h>
#include <timestepSaver.h>
#include <configurators.h>

#define RealType double
typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,7> > ConfigType;
using namespace aol::color;


int main( int argc, char **argv )
{
  try {
    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>\n";
      return EXIT_FAILURE;
    }

    aol::ParameterParser parser( argv[1] );

    char saveName[ 1024 ];
    parser.getString( "saveName", saveName );    // see below
    aol::TimestepSaver<RealType> tsSaver( parser.getInt("timeOffset"), saveName );

    // ------------- LOAD IMAGE AND GET SAVENAME----------------
    char imgname[ 1024 ];
    parser.getString( "image", imgname );
    cerr<<"\nLoading file "<<green<<imgname<<reset<<"...\n\n";
    qc::ScalarArray<double, qc::QC_2D> img( imgname );

    // Notice: It is important to scale the values to [0,1], because
    // if the difference of the values is too high, there might not exist
    // minimal surfaces (which are the solution of mcm in the graph case).
    img /= img.getMaxValue();

//     img.save("pics/RekonstrResults/anfang.pgm");

    parser.getString( "mask", imgname );
    cerr<<"\nLoading mask "<<green<<imgname<<reset<<"...\n\n";
    qc::ScalarArray<double, qc::QC_2D> mask( imgname );

    // the DofMask for the Dirichlet-Values, 0 = inner node, 1 = boundary node
    aol::DofMask boundaryMask;
    boundaryMask.createDofMaskFromVector( mask );
    qc::ScalarArray<double, qc::QC_2D> DirichletVals( img );
    for (int i=0; i<img.size(); i++)
      if (mask[i] == 0) DirichletVals[i]=0;       // take only the masked DirichletVals

//     DirichletVals.save("pics/DirichletVals.pgm");

    cerr<<"\ndone."<<reset;


    // -------------------- GRID AND TIMESTEP INFORMATION ---------------------------------
    int N = img.getNumX ();
    int d = qc::logBaseTwo (N);
    qc::GridDefinition grid( d, qc::QC_2D );
    qc::ScalarArray<double, qc::QC_2D> tmp( grid );
    qc::ScalarArray<double, qc::QC_2D> img_old( grid );

//     RealType tau = 0.5 * aol::Sqr( parser.getDouble("tau") * grid.H() );
    RealType tau = parser.getDouble("tau") * grid.H();
    RealType epsilon;
    if (parser.getInt("graph"))
      epsilon = 1.;
    else
      epsilon = parser.getDouble( "epsilon" ) *  grid.H();

    cerr << red << "declaring the op..." << reset;

    // -------------------- THE MCMSmoothOP -----------------------------------------------------

    qc::MCMDirichletSmoothOp< ConfigType > mcm(grid, boundaryMask, DirichletVals, tau, epsilon);
    // -------------------- MARKS A SPECIAL LEVELSET --------------------------------------------
    qc::LevelSetDrawer<ConfigType> drawer( grid );


    cerr<<cyan<<"\nepsilon: "<<epsilon<<endl<<reset;

    // ------------------- THE MAIN LOOP --------------------------------------------------
    img_old = img;

    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {

      aol::StopWatch watch;           // Stoppuhr
      watch.start();

      // apply the operator
      cerr << blue << "Step: "<<iter << reset << ", applying op...\n";

      mcm.apply( img_old, img );
      img_old = img;
      watch.stop();
      cerr << "ready, time-step took " << red << watch.elapsedCpuTime() << reset << "s.\n";

//       img.save("pics/RekonstrResults/modified.pgm");

      tmp = img;
      tmp *= 255.;
      tsSaver.saveTimestepPGM( iter, tmp );

    }


  } catch ( aol::Exception &ex ) {

    ex.dump();
  }
  cerr<<"vor beenden...";
  return EXIT_SUCCESS;
}











