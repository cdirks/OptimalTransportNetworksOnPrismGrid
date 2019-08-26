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
// anisotropic mean curvature motion in 2D


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


#define RealType double

typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,7> > ConfigType;

// typedef qcLpAnisotropy<double> AnisoType;
// typedef qc::Pedestal3d<double> AnisoType;
// typedef qc::Pentagon2d<double> AnisoType;
// typedef qc::Pentagon3d<double> AnisoType;
// typedef qc::SinePolygonRotated3d<double> AnisoType;
// typedef qc::L1Norm3d<double> AnisoType;
// typedef qc::Rotated3dAnisotropy<double, qc::L1Norm2d<double> > AnisoType;
// typedef qc::Rotated3dAnisotropy<double, qc::LInfNorm2d<double> > AnisoType;    // rotierte LInf-Norm = Kegel
// typedef qcEllipsoidAnisotropy<double> AnisoType;
// typedef qc::RotatedHexagon3d<double> AnisoType1;
typedef qc::RotatedHexagon3d<double> AnisoType;
// typedef qc::DoubleCone3d<double> AnisoType;
// typedef qc::DoubleCone3d<double> AnisoType2;
// typedef qc::LinearBlendedAnisotropy3d<ConfigType, double, AnisoType1, AnisoType2> AnisoType;

using namespace aol::color;

// now the main program
int main( int argc, char **argv ) {
  try {
    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>\n";
      return EXIT_FAILURE;
    }

    aol::ParameterParser parser( argv[1] );


    // ------------- LOAD IMAGE AND GET SAVENAME----------------
    char imgname[ 1024 ];
    parser.getString( "image", imgname );
    cerr<<green<<"\nLoading file "<<imgname<<"...\n\n";
    qc::ScalarArray<double, qc::QC_2D> img( imgname );
    if (img.getMaxValue() > 1.) img /= img.getMaxValue();
    aol::Vector<RealType> rhs( img );
    cerr<<"\ndone."<<reset;

    char saveName[ 1024 ];
    parser.getString( "saveName", saveName );    // see below
    aol::TimestepSaver<double> tsSaver( parser.getInt("timeOffset"), saveName );


    // **HACK ******
    qc::ScalarArray<double, qc::QC_2D> cleanImg( "Data2D/CircleNick/circleNick8Invert.bz2" );
    if (cleanImg.getMaxValue() > 1.) cleanImg /= cleanImg.getMaxValue();


    // -------------------- GRID AND TIMESTEP INFORMATION ---------------------------------
    int N = img.getNumX ();
    int d = qc::logBaseTwo (N);
    qc::GridDefinition grid( d, qc::QC_2D );

    RealType tau = parser.getDouble("tau") * aol::Sqr( grid.H() );
    double lambda = parser.getDouble("lambda");
    const RealType epsilon = parser.getInt( "graph" ) ? 1.0 : parser.getDouble( "epsilon" ) *  grid.H();



    // ------------------- The anisotropy and its operators ------------------------------
    // Lp-Norm with p=8
//     AnisoType aniso( 4,4, epsilon );
    // 2d-pedestal

//     double deltaAbs = 0.;
//     double deltaRad = 0.;
//     cout << "DeltaAbs, DeltaRad: ";
//     cin >> deltaAbs >> deltaRad;

/*    double deltaAbs = 0.1;
    double deltaRad = 0.1;
    double alpha = 1.047; */             // = pi/3 = 60°C

    // ******** Pedestal3d ********
//     AnisoType aniso( alpha, deltaAbs, deltaRad );
    // ******** Pentagon3d ********
//     AnisoType aniso( deltaAbs, deltaRad );
    // ******** Hexagon ********
//     AnisoType aniso( deltaAbs, deltaRad );
    // ******** Hexagon ********
//     AnisoType aniso( 1., 1., 3., M_PI/2., 0.02, 0.01 );
    // ******** diamond-like Quadrangle (Kegel, aus nicht konvexem gamma) ********
//     AnisoType aniso( 1., 3., 2., M_PI/4., 0.02, 0.01 );
    // ******** diamond-like Quadrangle (Kegel, aus konvexem gamma) ********
//     qc::LInfNorm2d<double> LInfAniso( epsilon );
//     qc::Rotated3dAnisotropy<double, qc::LInfNorm2d<double> > aniso( LInfAniso, epsilon );
    // ******** L1-Norm ********
//     AnisoType aniso( grid.H() );
    // ******** Rotated L1Norm ********
//     qc::L1Norm2d<double> L1Norm( grid.H() );
//     AnisoType aniso( L1Norm, 0.0001 );
    // ******** Ellipsoid ********
//     AnisoType aniso( 3., 1., 1., 5.*grid.H() );
    // ******** Rotated hexagon ******
//     AnisoType1 aniso1( grid.H(), 0.78 );
//     aniso1.setGraphFlag();
    AnisoType aniso( 4.*grid.H(), 0.78 );
    aniso.setGraphFlag();
    // ******** Double cone ******
//     AnisoType aniso( 4.*grid.H() );
//     aniso.setGraphFlag();
//     AnisoType2 aniso2( grid.H() );
//     aniso2.setGraphFlag();
    // ******** blending between two anisotropies ******
//     AnisoType aniso( aniso1, aniso2, grid, 0.015, 0.085, grid.H() );
//     aniso.setImageReference( img );

    aniso.setGraphFlag();


    // --------------------- Declaring all necessary ops and the matrix -------------------------
    cerr << red << "declaring the op..." << reset;
    qc::FastUniformGridMatrix<RealType, ConfigType::Dim > mat( grid );

    qc::MCMMassOp<ConfigType> mcmMassOp( grid, aol::ONTHEFLY, epsilon );
    mcmMassOp.setImageReference( img );

    aol::SSORPreconditioner<aol::Vector<RealType>,qc::FastUniformGridMatrix<RealType,ConfigType::Dim> > precond( mat );
    aol::PCGInverse<aol::Vector<RealType> > solver( mat, precond, 1e-18, 1000, aol::STOPPING_ABSOLUTE, cerr );

    // anisotropic operators
    qc::AnisotropyIntegrationOp< ConfigType, AnisoType >  mcmAnisoInt( grid, aniso );
//     qc::AnisotropyStiffOpPositionDepending< ConfigType, AnisoType >  mcmAnisoStiff( grid, aol::ASSEMBLED, img, aniso, epsilon );
    qc::AnisotropyStiffOp< ConfigType, AnisoType >  mcmAnisoStiff( grid, aol::ASSEMBLED, aniso, epsilon );
    mcmAnisoStiff.setImageReference( img );

    // ------------------- THE MAIN LOOP --------------------------------------------------

    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {

      // apply the operator
      cerr << blue << "Step: "<<iter << reset << ", building rhs...";

      mcmAnisoInt.apply( img, rhs );
      cerr << "Norm nach mcmAnisoInt: " << rhs.norm();
      rhs *= -tau;

      cerr << ", assembling matrices...";
      mat.setZero();

      mcmAnisoStiff.assembleAddMatrix( mat );
      mat *= tau*lambda;
      mcmMassOp.assembleAddMatrix( mat );

//       aniso.saveImage();

      cerr << ", solving...\n";
      solver.applyAdd( rhs, img );

      cerr << "ready\n";

      // ------------------ save every k th timestep ------------------------
      tsSaver.saveTimestepBZ2( iter, img, grid );

    }

    cerr << "\nDone.";
//     cerr <<blue<<" DeltaAbs: "<<deltaAbs<< ", DeltaRad: "<<deltaRad<<".\n"<<reset;


  } catch ( aol::Exception &ex ) {

    ex.dump();
  }
  return EXIT_SUCCESS;

}
