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

#include <ChanVese.h>
#include <gnuplotter.h>

typedef float RType;

int main( int /*argc*/, char **/*argv*/ ) {
  try {
    typedef aol::ZeroOneIntervalPenaltyFunction<RType> PenaltyFunctionType;
    PenaltyFunctionType penaltyFunction ( 0.1 );
    aol::PlotDataFileHandler<RType> plotDataFileHandler;
    plotDataFileHandler.generateFunctionAndDerivativePlot ( penaltyFunction, -2, 4, 1000 );

    aol::Plotter<RType> plotter;
    plotter.addPlotCommandsFromHandler ( plotDataFileHandler );
    plotter.set_outfile_base_name ( "test" );
    plotter.genPlot ( aol::GNUPLOT_PNG );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
