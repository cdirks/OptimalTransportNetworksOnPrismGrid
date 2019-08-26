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
// ---------------------------------------------------------------------------
// mcmVolPreserv.cpp
// calculates the mean curvature flow with an additional volume preserving
// term, based on the algorithm of Sethian (or Osher?)
// ---------------------------------------------------------------------------

#include <quoc.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <preconditioner.h>
#include <sparseMatrices.h>

#include <mcm.h>
#include <fastUniformGridMatrix.h>

#include <qmException.h>
#include <aol.h>
#include <parameterParser.h>
#include "grapeInterface3d.h"
#include <integrateLevelSet.h>
#include <timestepSaver.h>


typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;

using namespace aol::color;


int main( int argc, char **argv ) {

// read the parameters with an parameter-file
  if ( argc != 2 ) {
    string s = "USAGE: ";
    s += argv[0];
    s += " <parameterfile>";
    throw aol::Exception( s.c_str(), __FILE__, __LINE__  );
  }

  try {
    aol::ParameterParser parser( argv[1] );

    // ------------------- load the image ------------------------
    char imgname[ 1024 ];
    parser.getString( "image", imgname );
    char saveName[ 1024 ];
    parser.getString( "saveName", saveName );    // see below
    aol::TimestepSaver<double> tsSaver(0, parser.getInt("timeOffset"), saveName) ;

    cerr<<"\nLoading file "<<imgname<<"...\n\n";
    qc::ScalarArray<double, qc::QC_3D> img( imgname );
    cerr<<"done.";

    int N = img.getNumX ();
    int d = qc::logBaseTwo (N);

    qc::ScalarArray<double, qc::QC_3D> rhs( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> curvTemp( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> curvature( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> vector1( N,N,N );
    vector1.setAll( 1. );

    qc::GridDefinition grid( d, qc::QC_3D );
    double eps = 5. * grid.H();


    // -------------- Operator - definitions -------------------------------

    qc::MCMLumpedMassOp< ConfigType > mcmMass( grid, eps, aol::DO_NOT_INVERT );
    qc::MCMStiffOp< ConfigType > mcmStiff( grid, aol::ONTHEFLY, eps );
    qc::MCMStiffOp< ConfigType > mcmStiffAss( grid, aol::ASSEMBLED, eps );

    // inverse mass op for calculating the curvature (true = inverse)
    aol::LumpedMassOp< ConfigType > Minv( grid, aol::INVERT );
    aol::LumpedMassOp< ConfigType > M( grid, aol::DO_NOT_INVERT );
    // for integrating a function over one level set
//     qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> curvatureIntegrator( grid, img, curvature );
//     qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> massIntegrator( grid, img, vector1 );
    qc::LevelSetVolume<double> imgLevelSet( grid, img );

    qc::FastUniformGridMatrix<double, qc::QC_3D> mat( grid );

    aol::StopWatch watch;           // Stoppuhr
    watch.start();

    mcmStiff.setImageReference( img );       // Operatoren konfigurieren
    mcmStiffAss.setImageReference( img );

    // set the last time-step as data to calc |grad u|
    mcmMass.setImageReference( img );

    double tau = 0.5 * aol::Sqr( parser.getDouble("tau") * grid.H() );

    // filling the image-array with senseful values
    img /= img.getMaxValue();

    aol::SSORPreconditioner<aol::Vector<double>, qc::FastUniformGridMatrix<double,qc::QC_3D> > precond( mat );
    aol::PCGInverse<aol::Vector<double> > solver( mat, precond, 1e-16, 1000 );




    // ---------------- TIMESTEP-LOOP -----------------------------------

    double integralCurvature = 0;    // for calculating the vol-preserving term
    double integralLevelSet  = 0;
    double kAverage = 0;             // average curvature

    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {
      // Integrators must rebuild their cfe-structure for each time-step
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> curvatureIntegrator( grid, img, curvature );
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> massIntegrator( grid, img, vector1 );


      mcmMass.reset();
      mcmMass.setImageReference( img );    // own

      mat.setZero();
      cerr << blue << "step "<<iter<< reset << ", assembling and applying matrix.";
      mcmStiff.assembleAddMatrix( mat );

      // apply on img for the later calculation of the curvature
      mcmStiff.apply( img, curvTemp );

      cerr << ".";

      mat *= tau;
      cerr << ".";

      mcmMass.assembleAddMatrix( mat );      // this adds to the already assembled matrix
      cerr << ".done, ";

      cerr << "making rhs...";
      mcmMass.apply( img, rhs );

      // now the volume preserving term has to be added to the rhs
      // calculate the area of the levelset

      cerr << "calculating area of the leveset...";
      integralLevelSet  = massIntegrator.integrate();

      if (integralLevelSet != 0) {
        cerr << "done. calculating curvature...";
        // first calc the curvature by applying the inverse massOp
        Minv.apply( curvTemp, curvature);

        curvature *= -1.;
        integralCurvature = curvatureIntegrator.integrate();

        kAverage = integralCurvature / integralLevelSet;

        cerr << "\nIntegral const: "<<integralLevelSet<<", curvature: "<<integralCurvature;
        cerr << "\naverage Curvature: "<<kAverage<<endl;

        M.apply( vector1, curvTemp );
        curvTemp *= (kAverage * tau * (-1));

        // ready, add this to the rhs
        rhs += curvTemp;
      }
      else cerr << "\nATTENTION: area of level set is zero\n";

      cerr << "done.\n";

      // for testing purpose: get the volume of the levelset
//       cerr << "Finally find out the volume enclose by the levelset...";
//       cerr << "Vol = " << imgLevelSet.getVolume();
//       cerr << " !\n";

      solver.apply( rhs, img );


      // ------------------ save every "timeOffset"th timestep ------------------------
      if ( (iter % parser.getInt("timeOffset")) == 0 ) {
        cerr << "\n>>>>>> saving to file ";
        char filename[ 1024 ];
        sprintf( filename, "%s_%03d.bz2", saveName, iter );
        cerr << "\n>>>>>> saving to file " << filename << " <<<<<<\n\n";

        img.save( filename, qc::PGM_DOUBLE_BINARY );
      }
//       cerr<<"speichere...";
//       tsSaver.saveTimestep( iter, img);

    }

    // ----------------- END TIME-STEP-LOOP ---------------------------------


    watch.stop();
    cerr << "elapsed = " << watch.elapsedCpuTime() << "s.\n";

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}
