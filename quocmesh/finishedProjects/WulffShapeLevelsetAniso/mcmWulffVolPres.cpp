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
#include <preconditioner.h>
#include <solver.h>

#include <mcm.h>
#include <fastUniformGridMatrix.h>

#include <anisotropies.h>
#include <Willmore.h>

#include <qmException.h>
#include <aol.h>
#include <parameterParser.h>
#include <integrateLevelSet.h>
#include <timestepSaver.h>
#include <narrow.h>
#include <narrowBandGrid.h>
#include <narrowBandBoundaryIntegration.h>
#include <narrowBandConfigurators.h>
#include <signedDistanceOp.h>


using namespace aol::color;

typedef nb::NarrowBandGrid<qc::GridDefinition, qc::QC_3D> NarrowBandGridType;
typedef nb::NarrowBandConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3>, NarrowBandGridType > NBConfigType;

typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;

// typedef qcRotatedLpAnisotropy3d<double> AnisoType;
// typedef qcLpAnisotropy3d<double> AnisoType;
// typedef qc::RotatedEllipsoidAnisotropy<double> AnisoType;
// typedef qcL1Anisotropy3d<double> AnisoType;
// typedef qc::Pedestal3d<double> AnisoType;
// typedef qc::Pentagon3d<double> AnisoType;
// typedef qc::SinePolygonRotated3d<double> AnisoType;
// typedef qc::L1Norm3d<double> AnisoType;
// typedef qc::Rotated3dAnisotropy<double, qc::Hexagon2d<double> > AnisoType;    // rotierte LInf-Norm = Kegel
// typedef qc::EllipsoidAnisotropy<double> AnisoType;
// typedef qcIsotropy3d<double> AnisoType;
// typedef qc::DoubleCone3d<double> AnisoType;
// typedef qc::LInfNorm3d<double> AnisoType;
// typedef qc::RotatedHexagon3d<double> AnisoType;
typedef qc::RotatedHexagon3d<double> AnisoType2;
typedef qc::DoubleCone3d<double> AnisoType1;
typedef qc::LinearBlended3dAnisotropy<ConfigType, AnisoType1, AnisoType2> AnisoType;


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

    aol::TimestepSaver<double> tsSaver( parser.getInt("timeOffset"), saveName );
    tsSaver.setStepDigits( 5 );

    cerr<<"done.";

    int N = img.getNumX ();
    int d = qc::logBaseTwo (N);

    qc::GridDefinition grid( d, qc::QC_3D );

    qc::ScalarArray<double, qc::QC_3D> rhs( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> rhs2( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> img_start( img );

    qc::ScalarArray<double, qc::QC_3D> curvTemp( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> tmp( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> curvature( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> vector1( N,N,N );
    vector1.setAll( 1. );

    double eps = 5. * grid.H();

    // -------------- Operator - definitions -------------------------------

    qc::MCMLumpedMassOp< ConfigType > mcmMass( grid, eps, aol::DO_NOT_INVERT );


    // ******** the operator for the rhs including the anisotropy ***********
//     qcEllipsoidAnisotropy<double> ellipsoid(  2.,1.,1.,   eps);
//     qc::RotatedEllipsoidAnisotropy<double> ellipsoid(  1.,1.,1.,   eps);
//     // Ellipsoid in Diagonale drehen
//     aol::Vec3<double> v(1.0,1.0,1.0);
//     ellipsoid.setRotateDirection(v);
//
//     qc::AnisotropyIntegrationOp< ConfigType, qc::RotatedEllipsoidAnisotropy<double> >
//       mcmAnisoInt( grid, ellipsoid );
//
//     qc::AnisotropyStiffOp< ConfigType, qc::RotatedEllipsoidAnisotropy<double> >
//        mcmAnisoStiff( grid, aol::ASSEMBLED, img, ellipsoid, eps );
//
//
//

    // Test Lp-Norm as Wulffshape
//     AnisoType aniso(  1.1, 0.9090909091,  eps);


    // ****** Ellipsoid ******
//     AnisoType aniso(10.,1.,1.,   eps);
//     AnisoType aniso( 1., 1., 1., 5.*grid.H() );
    // ******  L1-Norm ******
//     AnisoType aniso( grid.H() );
    // ****** Pedestal ******
//     AnisoType aniso( 1.047, 0.1, 0.01 );
    // ****** Pentagon ******
//     AnisoType aniso( 0.02, 0.01 );
    // ****** Hexagon *******
//     AnisoType aniso( 1., 1., 3., M_PI/2., 0.02, 0.01 );
        // ******** diamond-like Quadrangle (Kegel) ********
//     AnisoType aniso( 1., 3., 2., M_PI/4., 0.02, 0.01 );
    // ******** diamond-like Quadrangle (Kegel, aus konvexem gamma) ********
//     qc::LInfNorm2d<double> LInfAniso( 0.001 );
//     qc::Hexagon2d<double> HexAniso( 0.001, 0.7853981632500000000440473768570370793896 );
//     qc::Rotated3dAnisotropy<double, qc::Hexagon2d<double> > aniso( HexAniso, grid.H() );
//     qc::Rotated3dAnisotropy<double, qc::LInfNorm2d<double> > aniso( LInfAniso, grid.H() );
    // ******** Double cone ******
//     AnisoType aniso( grid.H() );
    // ******** LInf-Norm 3d ******
//     AnisoType aniso( grid.H() );
    //     aniso.setRotateDirection( 1., 0., 0. );
        // ******** Rotated hexagon ******
//     AnisoType aniso( grid.H(), 0.7853981633974483096282022398515465511082 );    // pi/4
        // ******** double cone ******
    AnisoType1 aniso1( 0.01 );
    // ****** RotatedHexagon3d ******
    AnisoType2 aniso2( 0.01, 0.7853981632500000000440473768570370793896 );
    double xMin = 0.015;
    double xMax = 0.035;
    // ****** LinearBlendedAnisotropy2dGraph ******
    AnisoType aniso( aniso1, aniso2, grid, xMin, xMax, grid.H() );

    // Use a blending between anisotropies (0=aniso1 (DoubleCone), 1=aniso2 (RotatedHexagon3d))
    aniso.setBlendingValue( 0.5 );

//     aniso1.setGraphFlag();
//     aniso2.setGraphFlag();

    // ------------------ definitions of operators -------------------

    qc::AnisotropyIntegrationOp< ConfigType, AnisoType >
      mcmAnisoInt( grid, aniso );

    qc::AnisotropyStiffOp< ConfigType, AnisoType >
       mcmAnisoStiff( grid, aol::ASSEMBLED, aniso, eps );

    qc::FastUniformGridMatrix<double, qc::QC_3D> mat( grid );




    // ------------------ stuff for the boundary integral -------------------------------------
    NarrowBandGridType narrowGrid( grid );
    NBConfigType nbConfig( narrowGrid );

    // add all elements to iterator
    qc::GridDefinition::OldFullElementIterator elit;
    for ( elit = grid.begin(); elit != grid.end(); ++elit ) narrowGrid.insert( *elit );

    // determine the boundary of the narrow band
    narrowGrid.extractEdges();

    nb::IntegrateAnisotropyOverBoundary< double, NBConfigType, NarrowBandGridType, AnisoType, qc::QC_3D >
        nbAnisoBorderInt( nbConfig, narrowGrid, eps, aniso );

    // ------------------ end boundary stuff ---------------------------------------------------





    // -------- volume-preserving definitions -------------------

//     qc::AnisotropyIntegrationOp< ConfigType, qc::RotatedEllipsoidAnisotropy<double> >
//       mcmAnisoInt( grid, ellipsoid );
//     qc::AnisoGradUStiffOp<ConfigType,
//         qc::RotatedEllipsoidAnisotropy<double> > gammaZ(grid, ellipsoid, img, aol::ONTHEFLY);
    // inverse mass op for calculating the curvature (true = inverse)
    aol::LumpedMassOp< ConfigType > Minv( grid, aol::INVERT );
    aol::LumpedMassOp< ConfigType > M( grid, aol::DO_NOT_INVERT );
    // for integrating a function over one level set
//     qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> curvatureIntegrator( grid, img, curvature );
//     qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> massIntegrator( grid, img, vector1 );

    //  -------------------------------------------------------------



    aol::StopWatch watch;           // Stoppuhr
    watch.start();

    // set the last time-step as data to calc |grad u|
    mcmAnisoStiff.setImageReference( img );
    mcmMass.setImageReference( img );

    double tau = parser.getDouble("tau") * aol::Sqr( grid.H() );
    double lambda = parser.getDouble("lambda");

    aol::SSORPreconditioner<aol::Vector<double>, qc::FastUniformGridMatrix<double,qc::QC_3D> > precond( mat );
    aol::PCGInverse<aol::Vector<double> > solver( mat, precond, 1e-16, 1000, aol::STOPPING_ABSOLUTE, cerr  );

    // ------------------ the loop over the timesteps ----------------------------
    cerr<<endl<<endl;

    double integralCurvature = 0;    // for calculating the vol-preserving term
    double integralLevelSet  = 0;
    double kAverage = 0;             // average curvature

    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {
      // Integrators must rebuild their cfe-structure for each time-step, surface is 0-isosurface
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> curvatureIntegrator( grid, img, curvature );
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> massIntegrator( grid, img, vector1 );

      // the integral over gamma_z * phi
      cerr <<aol::color::blue<< "step " << iter <<aol::color::red<< ": Applying integrationOp...";
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
      cerr << "done.\nSolving ..."<<aol::color::black;

      //cerr<<"\n\nRHS nach mat: \n"<<rhs<<"\n\n\n";


      // ------------------------------------------------------------------------------------
      // now the volume preserving term has to be added to the rhs
      // calculate the area of the levelset
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


        cerr << "integrating over boundary of the narrowband ... ";
        tmp.setZero();
        nbAnisoBorderInt.apply( img, tmp );
        cerr << ".";
        tmp *= -tau;
        rhs += tmp;

      }
      else cerr << "\nATTENTION: area of level set is zero\n";

      cerr << "done.\n";
      // -------------------------------------------------------------------------------------



      solver.applyAdd( rhs, img );
      cerr<< ".done\n";

      // ------------------ save every 10 th timestep ------------------------
      tsSaver.saveTimestepBZ2( iter, img, grid );

      // recompute the signed distance function every k th timestep
      if ( parser.getInt("doRedist") && iter % parser.getInt("redistOffset") == 0 ) {
        cerr << red << "Redistancing...";
        qc::SignedDistanceOp3D<double>( grid ).apply( img, img );
        cerr << "ready.\n" << reset;
      }

    }


    watch.stop();
    cerr << "elapsed = " << watch.elapsedCpuTime() << "s.\n";


  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}
