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
 * \brief Equalizes the histogram of an 8-bit 2D ScalarArray.
 *
 * \author Berkels
 */

#include <scalarArray.h>

typedef double RType;
int main ( int argc, char **argv ) {
  try {
    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    qc::ScalarArray<unsigned char, qc::QC_2D> a ( inputFileName );
    const int numPixels = a.size();

    aol::Vector<int> histo;
    a.createHistogram( histo );
    const int numVals = histo.size();

    aol::Vector<int> transform ( numVals );
    int cumulativeDist = 0;
    for ( int i = 0; i < transform.size(); ++ i ) {
      cumulativeDist += histo[i];
      transform[i] = aol::Rint ( static_cast<double> ( cumulativeDist * ( numVals - 1 ) ) / numPixels );
    }

    for ( int i = 0; i < numPixels; ++ i )
      a[i] = transform[a[i]];

    a.savePNG( ( inputFileName + "_eq.png" ).c_str() );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
