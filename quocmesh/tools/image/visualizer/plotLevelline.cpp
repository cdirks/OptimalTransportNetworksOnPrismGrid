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

/**
 * \file
 * \brief Extracts a level line of a 2D level set function and plots it with gnuplot.
 *
 * \author Berkels
 */

#include <aol.h>
#include <op.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <multiArray.h>
#include <levelSetDrawer.h>
#include <timestepSaver.h>
#include <configurators.h>
#include <quocTimestepSaver.h>

typedef double RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType;

int main ( int argc, char **argv ) {

  try {
    char inputLevelsetFileName[1024];
    char outputFileName[1024];

    if ( argc < 3 ) {
      cerr << "USAGE: " << argv[0] << " <input_file> <out_file> [isovalue] [save-gnuplot-dat] [hideBorder] [second_input_file]" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }
    sprintf ( inputLevelsetFileName, "%s",  argv[1] );
    sprintf ( outputFileName, "%s",  argv[2] );
    const RType isovalue = ( argc > 3 ) ? atof(argv[3]) : 0.;
    const bool saveGnuplotDat = ( argc > 4 ) ? !!( atoi(argv[4]) ) : 0;
    const bool hideBorder = ( argc > 5 ) ? !!( atoi(argv[5]) ) : 0;
    const string secondInputLevelsetFileName = ( argc > 6 ) ? argv[6] : "";

    qc::ScalarArray<RType, qc::QC_2D> levelsetFunction ( inputLevelsetFileName );
    const ConfType::InitType grid ( qc::GridSize2d::createFrom ( levelsetFunction ) );

    const bool secondFunc = ( secondInputLevelsetFileName.size() > 0 );
    qc::ScalarArray<RType, qc::QC_2D> secondLevelsetFunction;
    if ( secondFunc )
      secondLevelsetFunction.load ( secondInputLevelsetFileName.c_str() );

    qc::QuocTimestepSaver<ConfType> saver ( 1, outputFileName );
    saver.saveIsolineWithGnuplot ( -1, -1, isovalue, aol::GNUPLOT_EPS, levelsetFunction, grid, hideBorder, secondFunc ? &secondLevelsetFunction : NULL );

    if ( saveGnuplotDat ) {
      qc::LevelSetDrawer<ConfType> drawer ( grid );
      string outputDatFileName = outputFileName;
      outputDatFileName += ".dat";
      drawer.drawToGnuplot ( levelsetFunction, outputDatFileName.c_str(), isovalue, true );
    }
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
