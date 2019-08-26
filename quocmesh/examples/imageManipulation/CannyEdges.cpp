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
 * \brief Detects the edges in an image using the Canny edge detector with different parameters.
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
    qc::ScalarArray<double, qc::QC_2D> adX ( a, aol::STRUCT_COPY_UNINIT );
    qc::ScalarArray<double, qc::QC_2D> adY ( a, aol::STRUCT_COPY_UNINIT );
    const RType sigma = 0.5;
    const int kernelSize = 5;
    const RType threshold = 1;
    for ( int i = 1; i < 10; ++i ) {
      const RType curSigma = i * sigma;
      qc::GaussDiffKernel2d<RType> kernelX ( kernelSize, curSigma, qc::DIFF_X );
      qc::GaussDiffKernel2d<RType> kernelY ( kernelSize, curSigma, qc::DIFF_Y );
      a.applyLinearFilterTo ( kernelX, adX );
      a.applyLinearFilterTo ( kernelY, adY );

      qc::ScalarArray<double, qc::QC_2D> gradMag ( a, aol::STRUCT_COPY_UNINIT );
      qc::ScalarArray<double, qc::QC_2D> edgeImage ( a, aol::STRUCT_COPY_UNINIT );
      for ( int y = 0; y < a.getNumY(); ++y )
        for ( int x = 0; x < a.getNumY(); ++x )
          gradMag.set ( x, y, sqrt ( aol::Sqr( adX.get(x,y) ) + aol::Sqr( adY.get(x,y) ) ) );

      // Ignore edges on the boundary.
      for ( int y = 1; y < a.getNumY()-1; ++y ) {
        for ( int x = 1; x < a.getNumY()-1; ++x ) {
          const RType gradAngle = aol::RadiansToDegrees ( atan2 ( adY.get(x,y), adX.get(x,y) ) );
          int gradAngleRounded = 45 * ( aol::Rint ( gradAngle / 45 ) );
          if ( gradAngleRounded < 0 )
            gradAngleRounded += 180;
          aol::Vec2<int> offset;
          switch ( gradAngleRounded ) {
            case 0:
            case 180:
              offset.set ( 1, 0 );
              break;
            case 45:
              offset.set ( 1, 1 );
              break;
            case 90:
              offset.set ( 0, 1 );
              break;
            case 135:
              offset.set ( -1, 1 );
              break;
            default:
              throw aol::Exception( "Invalid angle, this should not happen.", __FILE__, __LINE__);
              break;
          }
          const bool localMaximum = ( gradMag.get ( x, y ) > aol::Max ( gradMag.get ( x + offset[0], y  + offset[1] ) , gradMag.get ( x - offset[0], y  - offset[1] ) ) );

          edgeImage.set ( x, y, ( localMaximum && ( gradMag.get ( x, y ) > threshold ) )  ? 255 : 0 );
        }
      }

      edgeImage.savePNG ( aol::strprintf( "%s_CannyEdges%.1f.png", outputFileNameBase.c_str(), curSigma ).c_str() );
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
