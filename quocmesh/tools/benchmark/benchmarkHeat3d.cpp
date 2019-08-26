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

/** \file
 *  \brief benchmark solving the heat equation
 *
 *  Computes given number of time steps of heat conduction on given
 *  level (defaults: 7, 50) and measures runtime.  For default
 *  arguments, times are also printed relative to reference
 *  computation (optimized mode on single processor of brahma) With
 *  parameters bench file <filename>, uses default parameters and
 *  appends relative result to benchmark table in given file
 *
 *  Usage: benchmarkHeat3d <level> <numtimesteps> or benchmarkHeat3D bench file <filename>
 *
 *  \author Schwen
 */

#include <aol.h>
#include <quoc.h>
#include <solver.h>
#include <scalarArray.h>
#include <FEOpInterface.h>
#include <configurators.h>

int main ( int argn, char** argv ) {
  int level = 0;
  int numtimesteps = 0;

  bool bench = false;
  string resultfilename = "/dev/null";

  if ( aol::checkForBenchmarkArguments ( argn, argv, resultfilename ) )
    bench = true;

  if ( argn == 1 ) bench = true;

  if ( argn == 3 ) {

    level = atoi ( argv[1] );
    numtimesteps = atoi ( argv[2] );

  } else {

    level = 7;
    numtimesteps = 50;

  };

  try {

    cerr << "This test computes " << numtimesteps << " timesteps of 3D heat conduction on grid level " << level << ", using aol:: assembled operators and cg solver" << endl;

    aol::StopWatch timer;
    timer.start();

    qc::GridDefinition grid ( level, qc::QC_3D );

    aol::MassOp < qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > >  M_op ( grid, aol::ASSEMBLED );
    aol::StiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > >  L_op ( grid, aol::ASSEMBLED );

    aol::LinCombOp< aol::Vector<double>, aol::Vector<double> > HeatOp;

    const double tau = 0.5 * grid.H() ;
    HeatOp.appendReference ( M_op, 1.0 );
    HeatOp.appendReference ( L_op, tau*tau );
    aol::CGInverse< aol::Vector <double> > HeatSolver ( HeatOp, 1e-12, 50 );
    HeatSolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM  );

    qc::ScalarArray<double, qc::QC_3D> img ( grid ), dwarf ( grid );

    // set some initial values
    for ( int i = 0; i < grid.getWidth(); ++i ) {
      for ( int j = 0; j < grid.getWidth(); ++j ) {
        for ( int k = 0; k < grid.getWidth(); ++k ) {
          img.set ( i, j, k,   cos ( static_cast<double> ( i + j + k ) ) );    // pseudo-random data that should be the same on all platforms.
        }
      }
    }

    for ( int t = 0; t < numtimesteps; ++t ) {
      cerr << "Time step " << t << " of " << numtimesteps << endl;
      M_op.apply ( img, dwarf );
      HeatSolver.apply ( dwarf, img );
    }

    timer.stop();

    timer.printReport ( cerr );

    if ( bench ) {

      double nupsi = 100 * 983.14 / timer.elapsedCpuTime(), wupsi = 100 * 983.14 / timer.elapsedWallClockTime();

      // log benchmark results
      aol::logBenchmarkResult ( "benchmarkHeat3d", nupsi, wupsi, resultfilename );
      // Number of Uniform grid heat conduction timesteps Per Second In appropriate units: 100 nupsi = brahma, single job.
      // in terms of cputime and wall clock time
    }

    aol::callSystemPauseIfNecessaryOnPlatform();

  } catch ( aol::Exception &ex ) {

    ex.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();

  };

  return ( EXIT_SUCCESS );

}


