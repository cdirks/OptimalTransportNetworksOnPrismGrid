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
#include <quoc.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <preconditioner.h>

#include <mcm.h>
#include <fastUniformGridMatrix.h>

#include <anisotropies.h>
#include <Willmore.h>

#include <qmException.h>
#include <aol.h>
#include <parameterParser.h>
#include "grapeInterface3d.h"
#include <timestepSaver.h>
#include <integrateLevelSet.h>
#include <narrow.h>
#include <narrowBandGrid.h>
#include <narrowBandBoundaryIntegration.h>
#include <narrowBandConfigurators.h>

typedef double RealType;

typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;

typedef nb::NarrowBandGrid<qc::GridDefinition, qc::QC_3D> NarrowBandGridType;
typedef nb::NarrowBandConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3>, NarrowBandGridType > NBConfigType;

void testNAN( const aol::Vector<RealType> &Arg, const char *vecName )
{
  cerr<<aol::color::blue<<"\ntesting vector "<<aol::color::green<<vecName<<aol::color::blue<<"..."<<aol::color::reset;
  for (int i=0; i<Arg.size(); i++) {
    if ( aol::isNaN( Arg[i] ) ) {
      cerr << aol::color::red << "\n******* Found NAN!! ******* \n";
      return;
    }
    if ( !isfinite( Arg[i] ) ) {
      cerr << aol::color::red << "\n******* Found INFINITY!! ******* \n";
      return;
    }
  }
}

// now the main program
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



    // ---------- load image into scalar array -------------

    char saveName[ 1024 ];
    parser.getString( "saveName", saveName );

    aol::TimestepSaver<double> tsSaver( parser.getInt("timeOffset"), saveName );
    tsSaver.setStepDigits( 4 );

    char imgname[ 1024 ];
    parser.getString( "image", imgname );
    cerr<<"\nLoading file "<<imgname<<"...\n\n";
    qc::ScalarArray<double, qc::QC_3D> img( imgname );
    cerr<<"done.";

    int N = img.getNumX ();
    int d = qc::logBaseTwo (N);

    qc::GridDefinition grid( d, qc::QC_3D );

    qc::ScalarArray<double, qc::QC_3D> rhs( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> rhs2( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> tmp( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> img_start( img );

    double eps = 1. * grid.H();    // this was 5.



    // -------------- Operator - definitions -------------------------------

    qc::MCMLumpedMassOp< ConfigType > mcmMass( grid, eps, aol::DO_NOT_INVERT );

    // the anisotropy
    qc::RotatedEllipsoidAnisotropy<double> ellipsoid(  10.,1.,1.,   eps);

    // -------------------- loading vector field
    char loadName[ 1024 ];
    parser.getString( "loadRotationName", loadName );
    cerr<<"\nRestoring Vector-field ... \n";
    char filename[ 1024 ];

    sprintf( filename, "%sX.bz2", loadName);
    qc::ScalarArray<double, qc::QC_3D> rotationX( filename );
    sprintf( filename, "%sY.bz2", loadName);
    qc::ScalarArray<double, qc::QC_3D> rotationY( filename );
    sprintf( filename, "%sZ.bz2", loadName);
    qc::ScalarArray<double, qc::QC_3D> rotationZ( filename );

    if (rotationX.getNumX() != N || rotationY.getNumX() != N || rotationZ.getNumX() != N)
    {
      cerr<<"\n\nError: The vector-field doesn't have the same size as the image !!!!\n\n)";
      exit(43);
    }

    cerr<<"Restoring finished! \n\n";


    // further data fields for the volume preserving part
    qc::ScalarArray<double, qc::QC_3D> curvTemp( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> curvature( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> vector1( N,N,N );
    vector1.setAll( 1. );


    // --------------------------- loading finished ----------------------------------------

    qc::LocalRotatedAnisotropyIntegrationOp< ConfigType, qc::RotatedEllipsoidAnisotropy<double> >
      mcmAnisoInt( grid, ellipsoid, &rotationX, &rotationY, &rotationZ );

    qc::LocalRotatedAnisotropyStiffOp< ConfigType, qc::RotatedEllipsoidAnisotropy<double> >
       mcmAnisoStiff( grid, aol::ASSEMBLED, img, ellipsoid, eps,
                      &rotationX, &rotationY, &rotationZ );


    // -------- volume-preserving definitions -------------------
    // inverse mass op for calculating the curvature (true = inverse)
    aol::LumpedMassOp< ConfigType > Minv( grid, aol::INVERT );
    aol::LumpedMassOp< ConfigType > M( grid, aol::DO_NOT_INVERT );


    qc::FastUniformGridMatrix<double, qc::QC_3D> mat( grid );

    aol::StopWatch watch;           // Stoppuhr
    watch.start();

    // set the last time-step as data to calc |grad u|
    mcmAnisoStiff.setImageReference( img );
    mcmMass.setImageReference( img );

    double tau = 0.5 * aol::Sqr( parser.getDouble("tau") * grid.H() );
    double lambda = parser.getDouble("lambda");

    aol::SSORPreconditioner<aol::Vector<double>, qc::FastUniformGridMatrix<double,qc::QC_3D> > precond( mat );
    aol::PCGInverse<aol::Vector<double> > solver( mat, precond, 1e-18, 1000 );



    // ------------------ stuff for the boundary integral -------------------------------------
    NarrowBandGridType narrowGrid( grid );
    NBConfigType nbConfig( narrowGrid );

    // add all elements to iterator
    qc::GridDefinition::OldFullElementIterator elit;
    for ( elit = grid.begin(); elit != grid.end(); ++elit ) narrowGrid.insert( *elit );

    // determine the boundary of the narrow band
    narrowGrid.extractEdges();

    nb::IntegrateAnisotropyOverBoundary< double, NBConfigType, NarrowBandGridType, qc::RotatedEllipsoidAnisotropy<double>, qc::QC_3D >
        nbAnisoBorderInt( nbConfig, narrowGrid, eps, ellipsoid );

    // ------------------ end boundary stuff ---------------------------------------------------




    // ------------------ the loop over the timesteps ----------------------------
    cerr<<endl<<endl;

    double integralCurvature = 0;    // for calculating the vol-preserving term
    double integralLevelSet  = 0;
    double kAverage = 0;             // average curvature

    aol::StopWatch watchTimeStep;           // Stoppuhr


    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {

      watchTimeStep.start();
      // Integrators must rebuild their cfe-structure for each time-step, isosurface is 0-isosurface
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> curvatureIntegrator( grid, img, curvature );
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> massIntegrator( grid, img, vector1 );

      // ---------------- Reskalierung ---------------------------
      if ( parser.getInt( "rescale") ) {
        cerr<<"\nResizing the values to [0,1] ... ";

        img.addToAll ( - img.getMinValue() );
        img /= img.getMaxValue();
        cerr<<"ready!\n";
      }
      // --------------------------------------------------------


      // the integral over gamma_z * phi
      cerr << aol::color::blue  << "\nstep " << iter;
      cerr << aol::color::black << ": Applying rotatedIntegrationOp...";
      mcmAnisoInt.apply( img, rhs );
      rhs *= -tau;

      // now: M + Tau*lambda*L
      mcmMass.reset();

      mat.setZero();

      cerr << "assembling and applying (M + Tau*lambda*L).";
      mcmAnisoStiff.assembleAddMatrix( mat );
      cerr << ".";

      mat *= tau*lambda;
      //mat.scale( lambda );

      mcmMass.assembleAddMatrix( mat );      // this adds to the already assembled matrix
      cerr << ".done, finishing rhs ...";

      //mat.applyAdd( img, rhs );
      cerr << "done.\nSolving ...";

      //cerr<<"\n\nRHS nach mat: \n"<<rhs<<"\n\n\n";



      // ------------------------------------------------------------------------------------
      // now the volume preserving term has to be added to the rhs
      // calculate the area of the levelset
      // ------------------------------------------------------------------------------------
      cerr << "calculating area of the leveset...";
      integralLevelSet  = massIntegrator.integrate();

      if (integralLevelSet != 0) {
        cerr << "done. calculating anisotropy-curvature...";
        // first calc the curvature by applying the inverse massOp
        mcmAnisoInt.apply( img, curvTemp );
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

        // boundary-integral-stuff
        cerr << "integrating over boundary of the narrowband ... ";
        tmp.setZero();
        testNAN( tmp, "tmp before integration...");
        nbAnisoBorderInt.apply( img, tmp );
        testNAN( tmp, "tmp after integration...");

        cerr << ".";
        tmp *= -tau;
        rhs += tmp;
      }

      else cerr << "\nATTENTION: area of level set is zero\n";

      cerr << "done.\n";
      // -------------------------------------------------------------------------------------




      solver.applyAdd( rhs, img );
      watchTimeStep.stop();
      cerr << ".done."<< aol::color::blue  <<" Time-Step took ";
      cerr << watchTimeStep.elapsedCpuTime() <<"s. \n"<< aol::color::black;


      // ------------------ save every k'th th timestep ------------------------
      tsSaver.saveTimestepBZ2( iter, img, grid );

    }    // loop over time-steps


    watch.stop();
    cerr << "elapsed = " << watch.elapsedCpuTime() << "s.\n";


  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}
