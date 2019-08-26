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


typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;



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

    double eps = 5. * grid.H();

    // -------------- Operator - definitions -------------------------------

    qc::MCMLumpedMassOp< ConfigType > mcmMass( grid, eps, aol::DO_NOT_INVERT );


    // ******** the operator for the rhs including the anisotropy ***********
    //qcEllipsoidAnisotropy<double> ellipsoid(  1.,1.,1.,   eps);


    /// ATTENTION: TODO: x-factor HAS TO BE 2.
    // _----------------------------------------_/*;)§/(§M"Ä!"("BN
    qc::RotatedEllipsoidAnisotropy<double> ellipsoid(  100.,1.,1.,   eps);
    // Ellipsoid in Diagonale drehen
//     aol::Vec3<double> v(1.,1.,1.);
//     ellipsoid.setRotateDirection(v);

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

//     rotationX.setAll( 1. );
//     rotationY.setAll( 1. );
//     rotationZ.setAll( 1. );

    // -+#-+#-+#-+#-+#-+#-+#-+#-+#-+#-+#-+#-+#-+#-+-#+-#-+#-+-#-+#--+#--+#-+-#-+#-
    // HACK!! TODO: REMOVE IT AFTER TESTING!!!
    for (int i=0; i<N; i++)
      for (int j=0; j<N/2; j++)
        for (int k=0; k<N; k++)
        {
          rotationX.set(i,j,k, 0.);
          rotationY.set(i,j,k, 0.);
          rotationZ.set(i,j,k, 1.);
        }

    for (int i=0; i<N; i++)
      for (int j=N/2; j<N; j++)
        for (int k=0; k<N; k++)
        {
          rotationX.set(i,j,k, 0.);
          rotationY.set(i,j,k, 1.);
          rotationZ.set(i,j,k, 0.);
        }


//
    // --------------------------- loading finished ----------------------------------------
    //loadVectorField(loadName, rotationX, rotationY, rotationZ, N);

    qc::LocalRotatedAnisotropyIntegrationOp< ConfigType, qc::RotatedEllipsoidAnisotropy<double> >
      mcmAnisoInt( grid, ellipsoid, &rotationX, &rotationY, &rotationZ );

    qc::LocalRotatedAnisotropyStiffOp< ConfigType, qc::RotatedEllipsoidAnisotropy<double> >
       mcmAnisoStiff( grid, aol::ASSEMBLED, img, ellipsoid, eps,
                      &rotationX, &rotationY, &rotationZ );


    qc::FastUniformGridMatrix<double, qc::QC_3D> mat( grid );

    aol::StopWatch watch;           // Stoppuhr
    watch.start();

    // set the last time-step as data to calc |grad u|
    mcmAnisoStiff.setImageReference( img );
    mcmMass.setImageReference( img );

    double tau = 0.5 * aol::Sqr( parser.getDouble("tau") * grid.H() );
    double lambda = parser.getDouble("lambda");

    // filling the image-array with senseful values
    //img /= img.getMaxValue();

    aol::SSORPreconditioner<aol::Vector<double>, qc::FastUniformGridMatrix<double,qc::QC_3D> > precond( mat );
    aol::PCGInverse<aol::Vector<double> > solver( mat, precond, 1e-16, 1000 );

    // ------------------ the loop over the timesteps ----------------------------
    cerr<<endl<<endl;

    aol::StopWatch watchTimeStep;           // Stoppuhr

    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {

      watchTimeStep.start();


      // ---------------- Reskalierung ---------------------------
      if ( parser.getInt( "rescale") ) {
        cerr<<"\nResizing the values to [0,1] ... ";

        img.addToAll ( - img.getMinValue() );
        img /= img.getMaxValue();
        cerr<<"ready!\n";
      }
      // --------------------------------------------------------


      // the integral over gamma_z * phi
      cerr <<aol::color::blue << "step " << iter;
      cerr <<aol::color::black<< ": Applying rotatedIntegrationOp...";
      mcmAnisoInt.apply( img, rhs );


      // TESTING: ============0000000000000000000000===================000000000000000
      cerr<<aol::color::red;
      cerr<<"\nRHS-Norm: "<<norm(rhs,N)<<endl<<aol::color::black;
      // ===================00000000000000000000000======================0000000000000

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



      // TESTING: ============0000000000000000000000===================000000000000000
      qc::ScalarArray<double, qc::QC_3D> vec1( N,N,N );
      vec1.setAll(1.);
      qc::ScalarArray<double, qc::QC_3D> testvec( N,N,N );
      cerr<<aol::color::red;
      mat.apply( vec1, testvec );
      cerr<<"\nMass-Multiplied-Norm: "<<norm(testvec,N)<<endl<<aol::color::black;
      // ===================00000000000000000000000======================0000000000000



      // mat.applyAdd( img, rhs );
      cerr << "done.\nSolving ...";

      //cerr<<"\n\nRHS nach mat: \n"<<rhs<<"\n\n\n";

      solver.applyAdd( rhs, img );
      watchTimeStep.stop();
      cerr<< ".done. Time-Step took "<<watchTimeStep.elapsedCpuTime() <<"s. \n";

      // ------------------ save every 10 th timestep ------------------------
      if ( (iter % parser.getInt("timeOffset")) == 0 ) {
          cerr << "\n>>>>>> saving to file ";
          char filename[ 1024 ];
          sprintf( filename, "%s_%03d.bz2", saveName, iter );
          cerr << "\n>>>>>> saving to file " << filename << " <<<<<<\n\n";

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
