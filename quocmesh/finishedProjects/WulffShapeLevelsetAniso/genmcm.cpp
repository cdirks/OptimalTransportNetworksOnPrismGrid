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
//#include <array3D.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <preconditioner.h>
//#include <derivativeArray3D.h>
#include <sparseMatrices.h>

#include <mcm.h>
#include <fastUniformGridMatrix.h>


#include <qmException.h>
#include <aol.h>
#include <parameterParser.h>
#include "grapeInterface3d.h"



int main( int argc, char **argv ) {


  //testMatrix3D();

  if ( argc != 2 ) {
    cerr << "USAGE: " << argv[0] << " tau" << endl;
    return EXIT_FAILURE;
  }

  try {
    qc::GridDefinition grid( 5, qc::QC_3D );
    int N = grid.getWidth();

    double eps = 5. * grid.H();

    qc::MCMStiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D,
      aol::GaussQuadrature<double,qc::QC_3D,3> > > mcmStiff( grid, aol::ONTHEFLY, eps );
    qc::MCMStiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D,
      aol::GaussQuadrature<double,qc::QC_3D,3> > > mcmStiffAss( grid, aol::ASSEMBLED, eps );
    qc::MCMLumpedMassOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D,
      aol::GaussQuadrature<double,qc::QC_3D,3> > > mcmMass( grid, eps, aol::DO_NOT_INVERT );

    qc::FastUniformGridMatrix<double, qc::QC_3D> mat( grid );
    //aol::SparseMatrix<double> mat( grid );
    //qc::UniformGridSparseMatrix<double> mat( grid );

    qc::ScalarArray<double, qc::QC_3D> img( "Data/Ellipsoid/Ellipsoid6SD.bz2"); //"Data/Cube/CubeOrig5.bz2"  );
    qc::ScalarArray<double, qc::QC_3D> speed( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> rhs( N,N,N );

    // fill img with data
//     img.clear();
//     generateNoisyDiagLine(img, N);
//     generateBlock(img, N);
    cerr<<"Generiere Daten..."<<endl;
    qc::ScalarArray<double, qc::QC_3D> img_start( "Data/Ellipsoid/Ellipsoid6SD.bz2"); //"Data/Cube/CubeOrig5.bz2" );
//     generateNoisyDiagLine(img_start, N);
//     generateBlock(img_start, N);


    aol::StopWatch watch;           // Stoppuhr
    watch.start();

    mcmStiff.setImageReference( img );       // Operatoren konfigurieren
    mcmStiffAss.setImageReference( img );

    // set the last time-step as data to calc |grad u|
    mcmMass.setImageReference( img );

    double tau = 0.5 * aol::Sqr( atof( argv[1] ) * grid.H() );

    // filling the image-array with senseful values
    //img.load( "input.pgm" );
    img /= img.getMaxValue();
    //img.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. ); // only in 2D


    aol::SSORPreconditioner<aol::Vector<double>, qc::FastUniformGridMatrix<double,qc::QC_3D> > precond( mat );
    //aol::SSORPreconditioner<aol::Vector<double>, qc::UniformGridSparseMatrix<double> > precond( mat );
//     aol::SSORPreconditioner<aol::Vector<double>, aol::SparseMatrix<double> > precond( mat );
    aol::PCGInverse<aol::Vector<double> > solver( mat, precond, 1e-16, 1000 );




    // now for debugging tests
    /*mcmMass.assembleAddMatrix( mat );
    qc::ScalarArray<double, qc::QC_3D> test( N,N,N );
    test.setAll(1.);
    mat.apply(test, rhs);
    cerr << "Summe: "<<rhs.sum()<<endl;
    */



    //char filename[1024];
    for ( int i=0; i<7; i++ ) {
      mcmMass.reset();
      mcmMass.setImageReference( img );    // own

      mat.setZero();
      cerr << "step "<<i<<", assembling and applying matrix.";
      mcmStiff.assembleAddMatrix( mat );
      cerr << ".";

      mat *= tau;
      cerr << ".";

      mcmMass.assembleAddMatrix( mat );      // this adds to the already assembled matrix
      cerr << ".done, ";

      cerr << "making rhs...";
      mcmMass.apply( img, rhs );
      cerr << "done.\n";

      solver.apply( rhs, img );

    }
    watch.stop();
    cerr << "elapsed = " << watch.elapsedCpuTime() << "s.\n";

    //cerr<<"Ausgangsdaten: "<<endl<<img<<endl<<endl;

    img /= img.getMaxValue();

#ifdef USE_EXTERNAL_GRAPE
    // now the GRAPE-stuff to view the first and the last picture:
    GENMESH3D* mesh = quocmesh_convert_to_gmesh3d (&img, "berechnet");
    addScalarData(mesh, &img_start, "Ausgangsbild");
    initStartGrape(mesh, "berechnet");
#else
    cerr << "cannot display in GRAPE without grape external" << endl;
#endif



  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}
