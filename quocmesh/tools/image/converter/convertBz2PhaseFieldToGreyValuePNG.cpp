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
 *  \brief With this tool one can generate a grey value phase field png by a given Bz2 file. Moreover, if one uses this tool with 5 arguments one can generate a eps file with the function graph of two one dimensional functions saved as Bz2.
 *  \author Franken
 */

#include <aol.h>
#include <scalarArray.h>
#include <gnuplotter.h>

typedef float RType;
int main ( int argc, char **argv ) {

  try {
    char inputFileName[1024];
    char outputFileName[1024];
    char oneDInputFileName1[1024];
    char oneDInputFileName2[1024];
    char oneDOutputFileName[1024];

    cerr<<endl<<argc<<endl;

    if ( argc != 3 && argc != 6) {
      cerr << "USAGE: " << argv[0] << "  <input_file>  <output_file>" << endl;
      cerr << "or: " << argv[0] << "  <2dInput_file> <first1dInput_file> <second1dInput_file> <2dOutput_file> <1dOutput_file>" << endl;
      return EXIT_FAILURE;
    }
    if ( argc == 3 ) {
      sprintf ( inputFileName, "%s",  argv[1] );
      sprintf ( outputFileName, "%s",  argv[2] );
    }
    if ( argc == 6 ) {
      sprintf ( inputFileName, "%s",  argv[1] );
      sprintf ( oneDInputFileName1, "%s",  argv[2] );
      sprintf ( oneDInputFileName2, "%s",  argv[3] );
      sprintf ( outputFileName, "%s",  argv[4] );
      sprintf ( oneDOutputFileName, "%s",  argv[5] );
    }
    cerr << "Reading 2D bz2 from " << inputFileName << endl;

    qc::ScalarArray<RType, qc::QC_2D> a ( inputFileName );

    if ( a.getMaxValue() > 1 || a.getMinValue() < -1 ) {
      cerr<<"\nThis is no double well phase field function!\n";
    }

    a.addToAll ( 1 );
    a *= 120;

    a.savePNG ( outputFileName );

    if (argc == 6 ) {

      cerr << "Reading 1D bz2 from " << oneDInputFileName1 << endl;
      cerr << "and " << oneDInputFileName2 << endl;

      qc::ScalarArray<RType, qc::QC_1D> vec ( oneDInputFileName1 );
      qc::ScalarArray<RType, qc::QC_1D> vec2 ( oneDInputFileName2 );

//     aol::Vector<RType> vec ( 100 );
//     aol::Vector<RType> vec2 ( 100 );
//     for ( int i = 0; i < 100; ++i ) {
//       vec[i] = i;
//       vec2[i] = 100-i;
//     }

      aol::PlotDataFileHandler<RType> plotDataFileHandler;
      plotDataFileHandler.generateFunctionPlot ( vec, 0, 1 );
      plotDataFileHandler.generateFunctionPlot ( vec2, 0, 1 );

      aol::Plotter<RType> plotter;

      plotter.addPlotCommand( plotDataFileHandler.getDataFileNames()[0].c_str(), plotDataFileHandler.getDataFileStyles()[0] );
      plotter.addPlotCommand( plotDataFileHandler.getDataFileNames()[1].c_str(), plotDataFileHandler.getDataFileStyles()[1], 7 );


      plotter.set_outfile_base_name ( oneDOutputFileName );
      plotter.setSpecial ( "unset xtics\nunset ytics\nset border 3\nset size ratio 0.25\n" );
      plotter.genGrayEPS();

    }

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
