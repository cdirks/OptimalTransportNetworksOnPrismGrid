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
 *  \brief cross-section through a 2d image
 *
 *  Reads in a 2d image and produces a png and eps plot of a cross-section through the image.
 *
 *  \author Wirth
 */

#include <aol.h>
#include <op.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <gnuplotter.h>

typedef double RType;
int main ( int argc, char **argv ) {

  try {
    char inputFileName[1024];
    char baseOutputFileName[1024];
    char asciiFileName[1024];

    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << " <input_filename> <base_output_filename>" << endl;
      return EXIT_FAILURE;
    }

    sprintf ( inputFileName, "%s",  argv[1] );
    sprintf ( baseOutputFileName, "%s",  argv[2] );

    cerr << "Reading ScalarArray from " << inputFileName << endl;
    qc::ScalarArray<RType, qc::QC_2D> image ( inputFileName );

    aol::PlotDataFileHandler<RType> plotHandler;
    plotHandler.generateCenterLinePlot ( image );

    aol::Plotter<RType> plotter;
    plotter.addPlotCommandsFromHandler ( plotHandler );
    plotter.set_outfile_base_name ( baseOutputFileName );
    plotter.genPNG();
    plotter.genEPS();

#ifdef __MINGW32_VERSION
    system ( "PAUSE" );
#endif
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  return 1;
}
