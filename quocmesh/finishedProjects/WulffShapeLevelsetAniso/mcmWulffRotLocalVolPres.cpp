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
// #include "grapeInterface3d.h"
#include <integrateLevelSet.h>

// for saving the area and the volume in a file
#include "fileOps.h"

typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;

// typedef qcRotatedLpAnisotropy3d<double> AnisoType;
typedef qc::RotatedEllipsoidAnisotropy<double> AnisoType;

double norm(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  double res = 0.;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
        res += Data.get(i,j,k);

  return res;
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

    double eps = 1. * grid.H();    // factor was 5.



    // -------------- Operator - definitions -------------------------------

    qc::MCMLumpedMassOp< ConfigType > mcmMass( grid, eps, aol::DO_NOT_INVERT );

    // the anisotropy
//     qc::RotatedEllipsoidAnisotropy<double> ellipsoid(  5.,1.,1.,   eps);
    //     AnisoType aniso(  1.1, 0.9090909091,  eps);  // Lp-norm
    AnisoType aniso(  5.,1.,1.,   eps );            // ellipsoid

    // -------------------- loading vector field and pre-calculated radii
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
    qc::ScalarArray<double, qc::QC_3D> localCurvature( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> vector1( N,N,N );
    vector1.setAll( 1. );


    // --------------------------- loading finished ----------------------------------------

    // ------- operator to integrate over the anisotropy ----------------------------
    qc::LocalRotatedAnisotropyIntegrationOp< ConfigType, AnisoType >
      mcmAnisoInt( grid, aniso, &rotationX, &rotationY, &rotationZ );

    // -------------------- with anisotropy weighted stiff-op ---------------------------
    qc::LocalRotatedAnisotropyStiffOp< ConfigType, AnisoType >
       mcmAnisoStiff( grid, aol::ASSEMBLED, img, aniso, eps,
                      &rotationX, &rotationY, &rotationZ );


    // -------- volume-preserving definitions (see below) -------------------
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
    aol::PCGInverse<aol::Vector<double> > solver( mat, precond, 1e-16, 1000 );



    // ---------------------------------------------------------------------------
    // ------------------ the loop over the timesteps ----------------------------
    // ---------------------------------------------------------------------------

    cerr<<endl<<endl;

    double epsilon = parser.getDouble("epsilon");    // radius of the ball in which we integrate
    //epsilon *= grid.H();
    double ballMass = 0;                             // mass of the local ball
    double locCurvVal = 0;                           // local curvature in each point

    aol::Vec3<double> n,rotDirection;                 // the normal and the direction of the vector-field

    aol::StopWatch watchTimeStep;                    // Stoppuhr


    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {

      watchTimeStep.start();
      // Integrators must rebuild their cfe-structure for each time-step
      // these ops are for integrating over the whole 0-levelset and for getting the included volume
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> curvatureIntegrator( grid, img, curvature );
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> massIntegrator( grid, img, vector1 );
      qc::LevelSetVolume<double> imgLevelSet( grid, img );



      // ------------------------------------------------------------------
      // here is the part for examining the local volume preservation:
      // the area and the included volume of the level set are computed
      // ------------------------------------------------------------------

      double LSArea, LSVolume;
      cerr << aol::color::green << "calculating area and volume of the leveset...";
      LSArea   = massIntegrator.integrate();
      LSVolume = imgLevelSet.getVolume();
      cerr << "ready!\nArea: " << LSArea << ", Volume: " << LSVolume << aol::color::black << endl;

      char fileName[ 1024 ];
      parser.getString( "saveAreaFileName", fileName );
      appendAreaAndVolume( fileName, iter, LSArea );
      parser.getString( "saveVolumeFileName", fileName );
      appendAreaAndVolume( fileName, iter, LSVolume );

      // -------------------------------------------------------------------
      // end of examining part
      // -------------------------------------------------------------------





      // the integral over gamma_z * phi
      cerr <<aol::color::blue<< "step " << iter <<aol::color::black;
      cerr << ": Applying rotatedIntegrationOp...";
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


//          for (int i=0; i<32; i++)
//         for (int j=0; j<32; j++)
//           cerr<<mat.get(i,j);
//
      // ------------------------------------------------------------------------------------
      // now the volume preserving term has to be added to the rhs
      // calculate the area of the levelset
      // ------------------------------------------------------------------------------------

      cerr << "done. calculating anisotropy-curvature...";
      // first calc the curvature by applying the inverse massOp
      mcmAnisoInt.apply( img, curvTemp );
      Minv.apply( curvTemp, curvature);

      curvature *= -1.;
      curvature.save( "Data/Tests/Curvature.bz2", qc::PGM_DOUBLE_BINARY );

      aol::Vec3<double> coords;



      // now calculate the local curvature and the mass of the ball
      // for each point. Therefore use the epicenter of each element

      aol::StopWatch watchBall;        // this whole stuff is only for time-measure
      watchBall.start();
      double prozent;


      for (int i=0; i<N; i++) {

        for (int j=0; j<N; j++) {

          prozent = static_cast<double> ( ((i*N)+j) * 100. / (N*N) );
          cerr.precision(4);
          cerr<<"\rBerechne lokale Massen..."<<prozent<<" %              ";

          for (int k=0; k<N; k++) {

            coords[0] = i * grid.H();
            coords[1] = j * grid.H();
            coords[2] = k * grid.H();


            // -----------------------------------------------------------------------
            // now test, whether there are isosurfaces of the value 0 or not,
            // if there are not, it doesn't have to be integrated there.
            // This is tested by checking if all values have the same sign or not.
            qc::Ball epsBall(grid, epsilon);
            qc::Ball::iteratorBall ballIt;
            qc::Ball::iteratorBall ballItEnd( epsBall.end(coords) );
            ballIt = epsBall.begin(coords);

            bool allElementsSameSign = true;
            bool positive = true;


            double val = img.get( (*ballIt).x(), (*ballIt).y(), (*ballIt).z() );
            if (val < 0) positive = false;


            while (allElementsSameSign && (ballIt !=ballItEnd)) {
              val = img.get( (*ballIt).x(), (*ballIt).y(), (*ballIt).z() );

              if ((*ballIt).x()<0 || (*ballIt).x()>N || (*ballIt).x()<0 || (*ballIt).x()>N ||
                  (*ballIt).x()<0 || (*ballIt).x()>N )
              { cerr<<"Indizes ausserhalb des Bereichs: "<<(*ballIt).x()<<", "<<(*ballIt).y()<<", ";
                cerr<<(*ballIt).z()<<endl; }


              if (positive) {
                if (val < 0) allElementsSameSign = false;
              } else {
                if (val > 0) allElementsSameSign = false;
              }
              ballIt ++;
            }
            // ------ test ende ----------------------------------------------------


            if (!allElementsSameSign) ballMass = massIntegrator.integrateBall(epsilon, coords);
            else ballMass = 0;

            if (ballMass != 0) {
              locCurvVal = curvatureIntegrator.integrateBall(epsilon, coords);
              // the quotient is the local average curvature = vol.pres.-correction term
              // this has to be weighted. We want to have vol.pres. if the normal is
              // orthogonal to the structure, if it's parallel, then we want to have no
              // volume preservation
              // 1. get the normal
              qc::Element El( i,j,k, d );
              img.cellGradient(El, n);
              if ( n.norm() > 0.001 ) n /= n.norm();

              // 2. get the direction of the vector-field
              rotDirection[0] = rotationX.get( i,j,k );
              rotDirection[1] = rotationY.get( i,j,k );
              rotDirection[2] = rotationZ.get( i,j,k );
              if (rotDirection.norm() > 0.001) rotDirection /= rotDirection.norm();

              // now weight the quotient of local average curvature and ball mass
              double correctionTerm = (1. - abs( n*rotDirection ) ) * (locCurvVal / ballMass);

              localCurvature.set( i,j,k, correctionTerm );
            }  else  localCurvature.set( i,j,k, 0. );
          }

        }
      }

      watchBall.stop();
      cerr << "\rBerechnung lokaler Massen abgeschlossen (" << watchBall.elapsedCpuTime() << "s.)\n";

      M.apply( localCurvature, curvTemp );
      curvTemp *= tau * (-1);

      // ready, add this to the rhs
      rhs += curvTemp;

//       cerr<<rhs;


      cerr << "done. Solving\n";

      // -------------------------------------------------------------------------------------

      solver.applyAdd( rhs, img );
      watchTimeStep.stop();
      cerr<< ".done."<<aol::color::blue<<" Time-Step took "<<watchTimeStep.elapsedCpuTime();
      cerr <<"s. \n"<<aol::color::black;


      // ------------------ save every timeOffset'th timestep ------------------------
      if ( (iter % parser.getInt("timeOffset")) == 0 ) {
          cerr <<aol::color::green << "\n>>>>>> saving to file ";
          char filename[ 1024 ];
          sprintf( filename, "%s_%03d.bz2", saveName, iter );
          cerr << "\n>>>>>> saving to file " << filename << " <<<<<<\n\n"<<aol::color::black;

          img.save( filename, qc::PGM_DOUBLE_BINARY );
      }
    }    // loop over time-steps


    watch.stop();
    cerr << "elapsed = " << watch.elapsedCpuTime() << "s.\n";


  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}
