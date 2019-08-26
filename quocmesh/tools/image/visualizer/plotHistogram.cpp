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
 * \brief Plots the histogram of an 8-bit 2D ScalarArray and optionally includes the
 *        cumulative distribution in the plot (scaled to height of the histogram).
 *
 * \author Berkels
 */

#include <scalarArray.h>
#include <gnuplotter.h>

typedef double RType;
int main ( int argc, char **argv ) {
  try {
    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <plotCDF>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const bool plotCDF = ( atoi ( argv[2] ) != 0 );

    qc::ScalarArray<unsigned char, qc::QC_2D> a ( inputFileName );
    aol::Vector<int> histo;
    a.createHistogram( histo );
    aol::plotHistogram<RType> ( histo, aol::getBaseFileName ( inputFileName ) + "_histo", plotCDF );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
