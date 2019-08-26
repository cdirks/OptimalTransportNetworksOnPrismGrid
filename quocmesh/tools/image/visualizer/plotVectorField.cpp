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
 * \brief Plots a 2D vector field as arrows using gnuplot.
 *
 * \author Berkels
 */

#include <aol.h>
#include <op.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <gnuplotter.h>

/**
 * \author Berkels
 */
template <typename RealType>
void plotVectorFieldWithGnuplot ( qc::ScalarArray<RealType, qc::QC_2D> &DefX,
                                  qc::ScalarArray<RealType, qc::QC_2D> &DefY,
                                  const RealType Spacing,
                                  const bool SubtractAverage,
                                  const bool PlotAverage,
                                  const char *OutBaseName,
                                  const bool MarkRegion = false,
                                  const aol::Vec2<int> *RegionStart = NULL,
                                  const aol::Vec2<int> *RegionEnd = NULL ) {
  const aol::Vec2<RealType> meanDef ( DefX.getMeanValue(), DefY.getMeanValue() );

  if ( SubtractAverage ) {
    DefX.addToAll ( -meanDef[0] );
    DefY.addToAll ( -meanDef[1] );
  }

  ofstream vecoff ( "vec.dat" );
  qc::WriteVectorFieldAsGnuplotFile<RealType> ( vecoff, DefX, DefY, Spacing );
  vecoff.close();
  if ( PlotAverage ) {
    ofstream meanOff ( "mean.dat" );
    meanOff << 0.5 << " " << 0.5 << " " << meanDef[0] << " " << -meanDef[1] << endl;
    meanOff.close();
  }
  ofstream gnuplotdat ( "gnuplot.dat" );
  gnuplotdat << "set terminal postscript eps color\n";
  gnuplotdat << "set output \"" << OutBaseName << ".eps\"\n";
  gnuplotdat << "set xrange [0:1]\n";
  gnuplotdat << "set yrange [0:1]\n";
  if ( MarkRegion ) {
    const RealType h = 1 / static_cast<RealType> ( DefX.getNumXYZ() - 1 );
    gnuplotdat << "set object 1 rect from " << RegionStart->get(0) * h << "," << 1 - RegionStart->get(1) * h << " to " << RegionEnd->get(0) * h << "," << 1 - RegionEnd->get(1) * h << " fs empty border rgb \"green\"\n";
  }
  gnuplotdat << "unset xtics\n";
  gnuplotdat << "unset ytics\n";
  if ( DefX.getNumX() == DefX.getNumY() )
    gnuplotdat << "set size square 1, 1.4285715\n";
  else
    gnuplotdat << "set size ratio " << static_cast<RealType> ( DefX.getNumY() ) / static_cast<RealType> ( DefX.getNumX() ) << endl;
  gnuplotdat << "plot \"vec.dat\" w vec notitle lt -1";
  if ( PlotAverage )
    gnuplotdat << ", \"mean.dat\" w vec notitle head filled lt 1 lw 3";
  gnuplotdat << endl;
  gnuplotdat.close();

  aol::runGnuplot ( "gnuplot.dat" );
  if ( system ( aol::strprintf ( "epstopdf %s.eps", OutBaseName ).c_str() ) != EXIT_SUCCESS )
    cerr << "Warning: Calling epstopdf returned an error.\n";
  remove ( "gnuplot.dat" );
  remove ( "vec.dat" );
  if ( PlotAverage )
    remove ( "mean.dat" );
}

typedef double RType;
int main ( int argc, char **argv ) {

  try {
    char inputFileNameDefX[1024];
    char inputFileNameDefY[1024];

    if ( ( argc < 3 ) || ( argc > 8 ) ){
      cerr << "USAGE: " << argv[0] << "  <input_file def_x> <input_file def_y> [scale] [spacing] [subtractAverage] [plotAverage] [outBaseName]" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    if ( argc >= 3 ) {
      sprintf ( inputFileNameDefX, "%s",  argv[1] );
      sprintf ( inputFileNameDefY, "%s",  argv[2] );
    }
    const RType scale = ( argc >= 4 ) ? atof ( argv[3] ) : 1;
    const RType spacing = ( argc >= 5 ) ? atof ( argv[4] ) : 0.1;
    const bool subtractAverage = ( argc >= 6 ) ? ( atoi ( argv[5] ) != 0 ) : false;
    const bool plotAverage = ( argc >= 7 ) ? ( atoi ( argv[6] ) != 0 ) : false;
    const string outBaseName = ( argc >= 8 ) ? ( argv[7] ) : "plot";
    
    qc::ScalarArray<RType, qc::QC_2D> defX ( inputFileNameDefX );
    qc::ScalarArray<RType, qc::QC_2D> defY ( inputFileNameDefY );

    if ( scale != 1 ) {
      defX *= scale;
      defY *= scale;
    }

    plotVectorFieldWithGnuplot<RType> ( defX, defY, spacing, subtractAverage, plotAverage, outBaseName.c_str() );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
