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
// -----------------------------------------------------------------------------
// smoothVectorField.cpp
// function smoothes a loaded vector-field, saves it and shows it in GRAPE
// -----------------------------------------------------------------------------

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


typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;


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

    // load image into scalar array
    char loadName[ 1024 ];
    parser.getString( "loadName", loadName );
    cerr<<"Restoring image...";
    qc::ScalarArray<double, qc::QC_3D> img_old( loadName );
    qc::ScalarArray<double, qc::QC_3D> img( loadName );
    qc::ScalarArray<double, qc::QC_3D> rhs( loadName );

    cerr<<"nach deklaration";

    int N = img.getNumX ();
    int d = qc::logBaseTwo (N);
    qc::GridDefinition grid( d, qc::QC_3D );

    cerr<<" nach griddef.";

    // Operators
    aol::StiffOp< ConfigType > L( grid );
    aol::LumpedMassOp< ConfigType > M( grid, aol::DO_NOT_INVERT );  // false: do not inverse matrix

    cerr<<" nach ops.";

    double tau = 0.5 * aol::Sqr( parser.getDouble("tau") * grid.H() );
    cerr<<" nach tau";

    //qc::FastUniformGridMatrix<double, qc::QC_3D> mat( grid );

    aol::SparseMatrix<double> mat( grid );
    cerr<<" nach matrix...";
//     aol::SSORPreconditioner<aol::Vector<double>, qc::FastUniformGridMatrix<double,qc::QC_3D> > precond( mat );
    aol::SSORPreconditioner<aol::Vector<double>, aol::SparseMatrix<double> > precond( mat );
    cerr<<" nach precond. ";
    aol::PCGInverse<aol::Vector<double> > solver( mat, precond, 1e-16, 1000 );
    cerr<<" nach pcginverse";


    // ---------------------------------------------------------------------


    aol::StopWatch watch;           // Stoppuhr
    watch.start();

    // ---------------- Reskalierung ---------------------------
    if ( parser.getInt( "rescale01") ) {
      cerr<<"\nResizing the values to [0,1] ... ";

      img.addToAll ( - img.getMinValue() );
      img /= img.getMaxValue();
      cerr<<"ready!\n";
    }

    if ( parser.getInt( "rescaleMinus1_1") ) {
      cerr<<"\nResizing the values to [-1,1] ... ";

      double a = img.getMinValue();
      double b = img.getMaxValue();

      img.addToAll ( -a );
      img /= (b-a);
      img *= 2;        // width of [-1,1]
      img.addToAll ( -1. );

      cerr<<"ready!\n";
    }


    // Schleife mit der eigentlichen Glättung
    for ( int iter=0; iter<parser.getInt("timesteps"); ++iter ) {

      cerr << "assembling and applying (M + Tau*lambda*L).";
      mat.setZero();
      L.assembleAddMatrix( mat );  //
      cerr << ".";
      mat *= tau;
      M.assembleAddMatrix( mat );      // this adds to the already assembled matrix
      cerr << ".done, finishing rhs ...";

      cerr << "img : " << img.norm() << endl;
      M.apply( img, rhs );
      cerr << "rhs : " << rhs.norm() << endl;
      cerr << "done.\nSolving ...";

      solver.apply( rhs, img );
      cerr << "done!\n";

    }

    // Noch Konstante addieren?
    if ( parser.getInt( "addConstant") ) {
      cerr<<"\nAdding a constant the image... ";

      img.addToAll ( parser.getDouble( "Constant" ) );
      cerr<<"ready!\n";
    }


    watch.stop();
    cerr << "elapsed = " << watch.elapsedCpuTime() << "s.\n";


    // ---------------------------------------------------------------------

    // now save the smoothed vector field
    cerr<<"Saving the smoothed image ... \n";
    char saveName[ 1024 ];
    parser.getString( "saveName", saveName );

    img.save( saveName, qc::PGM_DOUBLE_BINARY );
    cerr<<"... sucessfully finished! \n\n";


    // finally view the result in GRAPE

//     GENMESH3D* mesh = quocmesh_convert_to_gmesh3d ( &img_old, "original picture" );
//     addScalarData( mesh, &img, "geglättete Daten" );
//     // and then start GRAPE, thats it!
//     initStartGrape(mesh, "smoothed VF");

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}

