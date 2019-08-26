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
 * \brief Applies Laplace sharpening with different strengths to a 2D ScalarArray.
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
    const string outputFileNameBase = aol::getBaseFileName ( inputFileName );
    qc::ScalarArray<double, qc::QC_2D> a ( inputFileName );
    qc::ScalarArray<double, qc::QC_2D> aSmooth ( a, aol::STRUCT_COPY_UNINIT );
    qc::ScalarArray<double, qc::QC_2D> laplaceA ( a, aol::STRUCT_COPY_UNINIT );
    qc::ScalarArray<double, qc::QC_2D> tmp ( a, aol::STRUCT_COPY_UNINIT );
    const RType sigma = 1.5;
    const int kernelSize = 5;
    qc::GaussKernel2d<RType> kernel ( kernelSize, sigma );
    qc::GaussDiffKernel2d<RType> kernelXX ( kernelSize, sigma, qc::DIFF_XX );
    qc::GaussDiffKernel2d<RType> kernelYY ( kernelSize, sigma, qc::DIFF_YY );
    a.applyLinearFilterTo ( kernel, aSmooth );
    a.applyLinearFilterTo ( kernelXX, laplaceA );
    a.applyLinearFilterTo ( kernelYY, tmp );
    laplaceA += tmp;

    laplaceA.setOverflowHandlingToCurrentValueRange();
    laplaceA.savePNG ( aol::strprintf( "%s_laplace.png", outputFileNameBase.c_str() ).c_str() );

    for ( int i = 1; i < 10; ++i ) {
      tmp.setSum ( aSmooth, laplaceA, -i );
      tmp.savePNG ( aol::strprintf( "%s_laplSharp%d.png", outputFileNameBase.c_str(), i ).c_str() );
    }

    aol::StopWatch watch;
    watch.start();

    cerr << "duration = " << watch.elapsedCpuTimeString() << endl;
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
