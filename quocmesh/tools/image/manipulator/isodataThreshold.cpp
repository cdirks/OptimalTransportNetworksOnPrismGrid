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
 * \brief Thresholds an image and determines the threshold automatically with the isodata algorithm.
 *
 * \author Berkels
 */

#include <scalarArray.h>
#include <gnuplotter.h>

int computeMean ( const aol::Vector<int> &Histo, int StartIndex, int EndIndex ) {
  int firstMoment = 0;
  for ( int i = StartIndex; i <= EndIndex; ++i )
    firstMoment += Histo[i];
  return firstMoment;
}

int computeFirstMoment ( const aol::Vector<int> &Histo, int StartIndex, int EndIndex ) {
  int firstMoment = 0;
  for ( int i = StartIndex; i <= EndIndex; ++i )
    firstMoment += i * Histo[i];
  return firstMoment;
}

typedef double RType;
int main ( int argc, char **argv ) {
  try {
    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <numBins>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const int numBins = atoi ( argv[2] );
    qc::ScalarArray<double, qc::QC_2D> a ( inputFileName );

    aol::Vector<int> histo;
    a.createHistogramOfValues ( histo, numBins );

    RType threshold = 0.5 * numBins;
    int iter = 0;
    do {
      ++iter;
      RType oldThreshold = threshold;

      int thresholdInt = aol::Rint ( threshold );
      threshold = static_cast<RType> ( computeFirstMoment ( histo, 0, thresholdInt - 1 ) ) / computeMean ( histo, 0, thresholdInt - 1 );
      threshold += static_cast<RType> ( computeFirstMoment ( histo, thresholdInt, numBins - 1 ) ) / computeMean ( histo, thresholdInt - 1, numBins - 1 );
      threshold *= 0.5;

      if ( ( aol::Abs ( threshold - oldThreshold ) < 0.5 ) || ( iter > 1000 ) )
        break;
    } while ( true );

    // The threshold was determined on the histogram and thus needs to be rescaled to the original range of the input image.
    const aol::Vec2<RType> minMax = a.getMinMaxValue();
    const RType scaledThreshold = ( threshold / numBins ) * ( minMax[1] - minMax[0] ) + minMax[0];

    cerr << "Determined " << scaledThreshold << " as threshold in " << iter << " iterations\n";

    a.threshold ( scaledThreshold, 0, 255 );
    a.savePNG ( ( inputFileName + "_thres.png" ).c_str() );

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
