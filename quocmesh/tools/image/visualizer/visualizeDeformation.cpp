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
 * \brief Visualizes a 2D deformation by showing a grid deformed with the deformation using gnuplot.
 *
 * \author Berkels
 */

#include <aol.h>
#include <op.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <gnuplotter.h>

typedef double RType;
int main ( int argc, char **argv ) {

  try {
    if ( ( argc < 3 ) || ( argc > 6 ) ){
      cerr << "USAGE: " << argv[0] << "  <input_file def_x> <input_file def_y> [scale] [lineDensity] [noBorder]" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileNameDefX = argv[1];
    const string inputFileNameDefY = argv[2];

    string outbasename = "plot";
    if ( aol::fileNameEndsWith ( inputFileNameDefX.c_str(), "_0.bz2" ) && aol::fileNameEndsWith (inputFileNameDefY.c_str(), "_1.bz2" ) )
      outbasename = aol::getBaseFileName ( inputFileNameDefX.substr ( 0, inputFileNameDefX.size() - 6 ) );
    else if ( aol::fileNameEndsWith ( inputFileNameDefX.c_str(), "_0.dat.bz2" ) && aol::fileNameEndsWith (inputFileNameDefY.c_str(), "_1.dat.bz2" ) )
      outbasename = aol::getBaseFileName ( inputFileNameDefX.substr ( 0, inputFileNameDefX.size() - 10 ) );

    qc::ScalarArray<RType, qc::QC_2D> defX ( inputFileNameDefX.c_str() );
    qc::ScalarArray<RType, qc::QC_2D> defY ( inputFileNameDefY.c_str() );

    if ( argc > 3 ) {
      const RType scale = atof ( argv[3] );
      defX *= scale;
      defY *= scale;
    }

    const int lineDensity = ( argc > 4 ) ? atoi ( argv[4] ) : 64;
    const bool noBorder = ( argc > 5 ) ? atoi ( argv[5] ) : false;

    qc::WriteDeformedGrid<RType>( "grid.dat", defX, defY, lineDensity );

    ofstream gnuplotdat ( "gnuplot.dat" );
    gnuplotdat << "set terminal postscript eps color\n";
    gnuplotdat << "set output \"" << outbasename << ".eps\"\n";
    gnuplotdat << "unset xtics\n";
    gnuplotdat << "unset ytics\n";
    if ( noBorder )
      gnuplotdat << "unset border\n";
    if ( defX.getNumX() == defX.getNumY() )
      gnuplotdat << "set size square 1, 1.4285715\n";
    else
      gnuplotdat << "set size ratio " << static_cast<RType> ( defX.getNumY() ) / static_cast<RType> ( defX.getNumX() ) << endl;
    gnuplotdat << "plot \"grid.dat\" w l notitle lt -1\n";
    gnuplotdat.close();

    aol::runGnuplot ( "gnuplot.dat" );
    system ( aol::strprintf ( "epstopdf %s.eps", outbasename.c_str() ).c_str() );
    remove ( "gnuplot.dat" );
    remove ( "grid.dat" );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
