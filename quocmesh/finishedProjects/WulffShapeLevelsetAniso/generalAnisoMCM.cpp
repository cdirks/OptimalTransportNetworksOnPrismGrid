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
#include <integrateLevelSet.h>

#include <timestepSaver.h>

#include <narrowBandConfigurators.h>
#include <narrow.h>
#include <narrowBandBoundaryIntegration.h>


typedef nb::NarrowBandGrid<qc::GridDefinition, qc::QC_3D> NarrowBandGridType;

typedef nb::NarrowBandConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3>, NarrowBandGridType > ConfigType;

// typedef qcLpAnisotropy3d<double> AnisoType;
// typedef qcRotatedLpAnisotropy3d<double> AnisoType;
typedef qc::RotatedEllipsoidAnisotropy<double> AnisoType;
// typedef cone3d<double> AnisoType;


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
    qc::ScalarArray<double, qc::QC_3D> img_start( img );

    double eps = 5. * grid.H();




    // -------------------- loading vector field and pre-calculated radii ------------------

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


    cerr<<"\nrestoring radii...";
    parser.getString( "RadiiName", loadName );
    qc::ScalarArray<double, qc::QC_3D> Radii( loadName );
    double radiusFactor = parser.getInt( "radiusFactor" );
    int useRadii = parser.getInt( "useRadii" );

    cerr<<"Restoring finished! \n\n";



    // ----------- further data fields for the volume preserving part ----------------------

    qc::ScalarArray<double, qc::QC_3D> curvTemp( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> tmp( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> curvature( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> localCurvature( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> curvatureIntegral( N-1,N-1,N-1 );
    qc::ScalarArray<double, qc::QC_3D> massIntegral( N-1,N-1,N-1 );
    qc::ScalarArray<double, qc::QC_3D> vector1( N,N,N );
    vector1.setAll( 1. );


    // --------------------------- loading finished ----------------------------------------


    // ----------------------- initialising the narrow band --------------------------------

    NarrowBandGridType narrowGrid( grid );
    ConfigType config( narrowGrid );

    int useNarrowBand = parser.getInt( "useNarrowBand" );

    qc::GridDefinition::OldFullElementIterator elit;
    for ( elit = grid.begin(); elit != grid.end(); ++elit )
    {
      if (useNarrowBand)   // then only insert elements where radius > 0
      {
        if ( Radii.get( (*elit)[0], (*elit)[1], (*elit)[2] ) > 0 ) narrowGrid.insert( *elit );
      } else {             // insert every element
        narrowGrid.insert( *elit );
      }
    }



    // --------------------- define the operators ------------------------------------------

    qc::MCMLumpedMassOp< ConfigType > mcmMass( narrowGrid, eps, aol::DO_NOT_INVERT );

    // --------------------- define the anisotropy   ---------------------------------------

    // ************ ellipsoid ************
    AnisoType aniso(  2.5,1.,1.,   eps);

    // ************ Test Lp-Norm as Wulffshape ***************
//     AnisoType aniso(  0.9090909091, 1.1,  eps);
//     aniso.setRotateDirection( 1., 1., 1. );
//     aniso.printB();

    // ************ cone3d as Wulffshape ***************
//     AnisoType aniso(  15.,5.,   0.01);


    // integrate over gamma_z( grad u )
//     qc::LocalRotatedAnisotropyIntegrationOp< ConfigType, AnisoType >
//       mcmAnisoInt( narrowGrid, aniso, &rotationX, &rotationY, &rotationZ );

    // border-integral of the narrowband
    nb::IntegrateAnisotropyOverBoundary< double, ConfigType, NarrowBandGridType,
                                         AnisoType, qc::QC_3D >
      nbAnisoBorderInt( config, narrowGrid, eps, aniso );

//     qc::LocalRotatedAnisotropyStiffOp< ConfigType, AnisoType >
//        mcmAnisoStiff( narrowGrid, aol::ASSEMBLED, img, aniso, eps,
//                       &rotationX, &rotationY, &rotationZ );




// ########################### HACK!! ######################################

qc::AnisotropyIntegrationOp< ConfigType, AnisoType >
      mcmAnisoInt( narrowGrid, aniso );

    qc::AnisotropyStiffOp< ConfigType, AnisoType >
       mcmAnisoStiff( narrowGrid, aol::ASSEMBLED, aniso, eps );

//############################################################################




    // inverse mass op for calculating the curvature (true = inverse)
    aol::LumpedMassOp< ConfigType > Minv( narrowGrid, aol::INVERT );
    aol::LumpedMassOp< ConfigType > M( narrowGrid, aol::DO_NOT_INVERT );
    nb::SubGridSparseMatrix<double,ConfigType> mat( config );
//     aol::SparseMatrix<double> mat( grid );

    aol::StopWatch watch;           // Stopwatch for the whole time
    watch.start();

    // set the last time-step as data to calc |grad u|
    mcmAnisoStiff.setImageReference( img );
    mcmMass.setImageReference( img );

    double tau = 0.5 * aol::Sqr( parser.getDouble("tau") * grid.H() );
    double lambda = parser.getDouble("lambda");

//     aol::SSORPreconditioner<aol::Vector<double>, aol::SparseMatrix<double> > precond( mat );
    aol::SSORPreconditioner<aol::Vector<double>, nb::SubGridSparseMatrix<double,ConfigType> > precond( mat );
    aol::PCGInverse<aol::Vector<double> > solver( mat, precond, 1e-16, 1000 );



    // --------------------------- defining several Ball-iterators ----------------------------------
    // The elements belonging to this iterators are pre-calculated and needn't be calculated
    // again for every element. But for that the calculation is not as exact as before, because
    // the radii have to be discretized.

    double finestBallStep = 0.001;

    qc::Ball *epsBall[200];
    for (int i=0; i<200; i++)
      epsBall[i] = new qc::Ball( grid, (i+1)*finestBallStep );



    // ---------------------------------------------------------------------------
    // ------------------ the loop over the timesteps ----------------------------
    // ---------------------------------------------------------------------------

    cerr<<endl<<endl;

    double ballMass = 0;                             // mass of the local ball
    double locCurvVal = 0;                           // local curvature in each point

    aol::Vec3<double> n,rotDirection;                 // the normal and the direction of the vector-field

    aol::StopWatch watchTimeStep;                    // Stoppuhr


    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {

      watchTimeStep.start();
      // Integrators must rebuild their cfe-structure for each time-step
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> curvatureIntegrator( grid, img, curvature );
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> massIntegrator( grid, img, vector1 );


      // the integral over gamma_z * phi
      cerr << aol::color::blue << "step " << iter;
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
      mcmMass.assembleAddMatrix( mat );      // this adds to the already assembled matrix



      // ------------------------------------------------------------------------------------
      // now the volume preserving term has to be added to the rhs
      // calculate the area of the levelset
      // ------------------------------------------------------------------------------------

      cerr << "done, calculating anisotropy-curvature...";
      // first calc the curvature by applying the inverse massOp
      mcmAnisoInt.apply( img, curvTemp );
      Minv.apply( curvTemp, curvature);

      curvature *= -1.;


      cerr << "done, integrating over curvature...";
      curvatureIntegral.setZero();
      curvatureIntegrator.integrateSaveEachValue( curvatureIntegral );

      cerr << "done, integrating over mass ...";
      massIntegral.setZero();
      massIntegrator.integrateSaveEachValue( massIntegral );
      cerr << "done.\n";



      aol::StopWatch watchBall;        // this whole stuff is only for time-measure
      watchBall.start();
      double prozent;


      // ----------------------------------------------------------------------------
      // now calculate the local curvature and the mass of the ball
      // for each point. Therefore use the epicenter "coords" of each element
      aol::Vec3<double> coords;

      for (int i=0; i<N; i++) {

        for (int j=0; j<N; j++) {

          prozent = static_cast<double> ( ((i*N)+j) * 100. / (N*N) );
          cerr.precision(4);
          cerr << aol::color::red <<"\rBerechne lokale Massen..."<<prozent<<" %              ";

          for (int k=0; k<N; k++) {

            coords[0] = i * grid.H();
            coords[1] = j * grid.H();
            coords[2] = k * grid.H();

            // we only need to calculate the average curvature for those pixels, that
            // have a radius > 0
            double radius = 1.;
            if (useRadii) radius = Radii.get( i,j,k );

            if (radius > 0)
            {
              // ------------------------------------------------------------------------------------
              // faster method (hopefully) for integrating over a the levelset in a ball
              qc::Ball::iteratorBall ballIt;
              qc::Ball::iteratorBall ballItEnd;
              qc::Ball::iteratorBall ballItBegin;

              int index=-1;
              double ballRad = radius*radiusFactor*grid.H();

              if ( ballRad >= 0.2 )
              {
                qc::Ball *epsBallNeu = new qc::Ball( grid, ballRad );
                ballItBegin = epsBallNeu->begin(coords);
                cerr<<"Radius zu gross!!\n";
              }
              else
              {
                index = static_cast<int>(ballRad / finestBallStep);
                ballItBegin = epsBall[index]->begin(coords);
                ballItEnd = epsBall[index]->end(coords);
              }

              ballMass = 0;
              // add the parts of the pre-calculated integrals of each element
              for ( ballIt = ballItBegin; ballIt != ballItEnd; ballIt++ ) {
                ballMass += massIntegral.get( (*ballIt)[0], (*ballIt)[1],(*ballIt)[2] );
              }

              // only calculate the integral over the curvature if there is part of the
              // iso-surface in the ball.
              if (ballMass != 0)
              {
                locCurvVal = 0;

                for ( ballIt = ballItBegin; ballIt != ballItEnd; ballIt++ ) {
                  locCurvVal += curvatureIntegral.get( (*ballIt)[0], (*ballIt)[1],(*ballIt)[2] );
                }


//                 locCurvVal = curvatureIntegrator.integrateBall(ballRad, coords);
//                 cerr<<aol::color::green;
//                 cerr<<"("<<ballRad<<","<<index<<","<<fastLocCurvVal<<","<< locCurvVal<<")\n";
//                 cerr<<aol::color::red;


                // -----------------------------------------------------------------------------------
                // the quotient is the local average curvature = vol.pres.-correction term
                // this has to be weighted. We want to have vol.pres. if the normal is
                // orthogonal to the structure, if it's parallel, then we want to have no
                // volume preservation
                // 1. get the normal
                qc::Element El( i,j,k, d );
                img.cellGradient(El, n);
                double norm = n.norm();
                if (norm != 0) n /= norm;

                // 2. get the direction of the vector-field
                rotDirection[0] = rotationX.get( i,j,k );
                rotDirection[1] = rotationY.get( i,j,k );
                rotDirection[2] = rotationZ.get( i,j,k );
                norm = rotDirection.norm();
                if (norm != 0) rotDirection /= rotDirection.norm();

                // now weight the quotient of local average curvature and ball mass
                double correctionTerm = (1./* - abs( n*rotDirection )*/) * (locCurvVal / ballMass);

                localCurvature.set( i,j,k, correctionTerm );
              }
              else  localCurvature.set( i,j,k, 0. );
            }
            else  localCurvature.set( i,j,k, 0. );
          }

        }
      }

      watchBall.stop();
      cerr << "\rBerechnung lokaler Krümmungen abgeschlossen " << aol::color::black;
      cerr << "(" << watchBall.elapsedCpuTime() << "s.)\n";

      M.apply( localCurvature, curvTemp );
      curvTemp *= tau * (-1);

      // ready, add this to the rhs
      rhs += curvTemp;

      cerr << "integrating over boundary of the narrowband ... ";
      tmp.setZero();
      nbAnisoBorderInt.apply( img, tmp );
      cerr << ".";
      tmp *= -tau;
      rhs += tmp;


      // ------------ before solving, test the diagonal-entries of the matrix ---------------
//       cerr << "done, checking diagonal ... ";
//       for (int i=0; i<mat.getNumRows(); i++)
//         if ( mat.get(i,i) == 0 ) mat.set(i,i,1.);


      cerr << "done. Solving\n";

      // -------------------------------------------------------------------------------------

      solver.applyAdd( rhs, img );
      watchTimeStep.stop();
      cerr<< ".done." << aol::color::blue <<  "Time-Step took "<<watchTimeStep.elapsedCpuTime() <<"s. \n";
      cerr << aol::color::black;


      // ------------------ save every timeOffset'th timestep ------------------------
      tsSaver.saveTimestepBZ2( iter, img, grid );

    }    // loop over time-steps


    watch.stop();
    cerr << "elapsed = " << watch.elapsedCpuTime() << "s.\n";


  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}
