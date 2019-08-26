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
#include <mcm.h>
#include <anisotropies.h>
#include <Willmore.h>
#include <timestepSaver.h>
#include <levelSet.h>
#include <parameterParser.h>
#include <scalarArray.h>
#include <narrow.h>
#include <anisotropyVisualization.h>

typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;

using namespace aol::color;

// typedef qcRotatedLpAnisotropy3d<double> AnisoType;
// typedef qcLpAnisotropy3d<double> AnisoType;
// typedef qc::RotatedEllipsoidAnisotropy<double> AnisoType;
// typedef qcL1Anisotropy3d<double> AnisoType;
// typedef qc::Pedestal3d<double> AnisoType;
// typedef qcEllipsoidAnisotropy<double> AnisoType;
// typedef qc::SinePolygonRotated3d<double> AnisoType;
typedef qc::Rotated3dAnisotropy<double, qc::LInfNorm2d<double> > AnisoType;
// typedef qc::L1Norm3d<double> AnisoType;
// typedef qc::DoubleCone3d<double> AnisoType;
// typedef qc::RotatedHexagon3d<double> AnisoType;

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
    cerr<<green<<"\nLoading file "<<imgname<<"...\n\n"<<reset;
    qc::ScalarArray<double, qc::QC_3D> img( imgname );

    int N = img.getNumX ();
    int d = qc::logBaseTwo (N);

    qc::GridDefinition grid( d, qc::QC_3D );
    aol::TimestepSaver<double> tsSaver( parser.getInt("timeOffset"), saveName );

    // -------------- define the anisotropy --------------------------------
    // ******** Pedestal ********
//     AnisoType aniso( 1.047, 0.2, 0.1 );
    // ******** Hexagon ********
//     AnisoType aniso( 1., 1., 3., M_PI/2., 0.02, 0.01 );
    // ******** diamond-like Quadrangle (Kegel) ********
//     AnisoType aniso( 1., 3., 2., M_PI/4., 0.02, 0.01 );
    // ******** L1Norm3d ********
//     AnisoType aniso( 0. /*grid.H()*/ );
    // ******** Rotated L1Norm/LInfNorm ********
    qc::LInfNorm2d<double> aniso2D( grid.H() );
    // ******** Rotated Ellipse ********
//     qcEllipseAnisotropy<double> aniso2D( 10., 1., grid.H() );
    // ******** Rotate a 2d-anisotropy ********
    AnisoType aniso( aniso2D, 0.0001 );
    // ******** Double cone ******
//     AnisoType aniso( 0. );
    // ******** Rotated hexagon ******
//     AnisoType aniso( grid.H(), 0.78 );


    // -------------- Operator - definitions -------------------------------

    qc::GrowWulffshapeEO3D<double, AnisoType> EnquistOsherSchemeOp( grid, aniso );
    EnquistOsherSchemeOp.setData( &img );

      // ------------ get the timestep - size -----------------------------
    double tau = parser.getDouble("tau") * grid.H();


    // ------------------ the loop over the timesteps ----------------------------
    cerr<<endl<<endl;

    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {

      // the integral over gamma_z * phi
      cerr << blue << "step " << iter << endl;


      EnquistOsherSchemeOp.timeStepEO( tau );


      // ------------------ save every 10 th timestep ------------------------
      tsSaver.saveTimestepBZ2( iter, img, grid );

    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}
