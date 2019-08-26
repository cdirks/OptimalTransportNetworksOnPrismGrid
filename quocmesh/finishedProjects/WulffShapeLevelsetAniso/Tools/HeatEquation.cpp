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
/** *************************************************************************
 * small program to smooth an image with the linSmoothOp, by Oliver Nemitz
 * (last date: 10.10.2006)
 * ************************************************************************* */

#include <linearSmoothOp.h>
#include <quoc.h>
#include <FEOpInterface.h>
// #include <solver.h>
// #include <sparseSolver.h>

// #include <fastUniformGridMatrix.h>
#include <qmException.h>
#include <aol.h>
// #include <parameterParser.h>

#include <timestepSaver.h>
// #include <anisoStiffOps.h>
// #include <mcm.h>

typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;

using namespace aol::color;


// now the main program
int main( int argc, char **argv ) {

  // read the parameters with an parameter-file
  if ( argc < 4 ) {
    cerr << "USAGE: " << argv[0] << " image outputTrunc steps [saveOffset=5 filterwidth=1 (*grid.H())]\n";
    return EXIT_FAILURE;
  }

  try {
    qc::ScalarArray<double, qc::QC_2D> image( argv[1] );
    int anzSteps = atoi( argv[3] );
    int timeOffset = 5;
    if (argc > 4 ) timeOffset = atoi( argv[4] );
    double fw = 1.;
    if (argc > 5 ) fw = atof( argv[5] );

    int N = image.getNumX ();
    int d = qc::logBaseTwo (N);
    qc::GridDefinition grid( d, qc::QC_2D );

    // ---------------- the timestep-saver -----------------------------
//     char saveName[1000];
//     sprintf( saveName, "%sHeat", argv[1] );
    aol::TimestepSaver<double> tsSaver( timeOffset, argv[2] );
    tsSaver.setStepDigits( 4 );

    // ------------------ Operators ------------------------------------
    qc::LinearSmoothOp<double> smoothOp;
    smoothOp.setCurrentGrid( grid );
    smoothOp.setSigma( fw*grid.H() );

    cerr << endl << "Now applying " << anzSteps << " of Heat Equation with filterwidth " << fw << "*grid.H(), every ";
    cerr << timeOffset << " image will be saved! "<< endl;
    for (int i=0; i<anzSteps; i++) {
      cerr << "\rSmoothing the image, step "<<red<<i+1<<reset<<" of "<<anzSteps<<"...";
      smoothOp.apply( image, image );
      tsSaver.saveTimestepPGM( i, image, grid );
    }

    cerr<<"ready!\n";

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}

