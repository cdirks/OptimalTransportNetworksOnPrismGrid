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
 * \brief Applies a high and low pass filter to a 2D ScalarArray.
 *        Also saves the corresponding kernels as PNGs.
 *
 * \author Berkels
 */

#include <scalarArray.h>
#include <configurators.h>
#include <generator.h>
#include <convolution.h>
#include "frequencyFilters.h"

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main ( int argc, char **argv ) {
  try {
    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <factor>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const int undersampleFactor = atoi ( argv[2] );
    const string outputFileNameBase = aol::strprintf( "%s_%d", aol::getBaseFileName ( inputFileName ).c_str(), undersampleFactor );
    const qc::ScalarArray<RType, qc::QC_2D> a ( inputFileName );

    ConfType::InitType grid ( qc::GridSize<ConfType::Dim>::createFrom ( a ) );
    qc::MultiArray<RType, 2, 2> image ( grid ), imageTransformed ( grid ), imageMultiplied ( grid ), lowPass ( grid ), highPass ( grid );
    qc::ScalarArray<RType, qc::QC_2D> tmp ( a, aol::STRUCT_COPY_UNINIT );
    qc::ScalarArray<RType, qc::QC_2D> tmp2 ( a, aol::STRUCT_COPY_UNINIT );
    image[0] = a;

    qc::DataGenerator<ConfType> generator ( grid );
    const RType centerX = 0.5 * grid.H() * ( grid.getNumX() - 1 );
    const RType centerY = 0.5 * grid.H() * ( grid.getNumY() - 1 );
    generator.generateCircle( tmp, 1, 0.5 / undersampleFactor, centerX, centerY );
    generator.generateCircle( tmp2, 1, 1. / undersampleFactor, centerX, centerY );
    tmp2 -= tmp;

    tmp.setOverflowHandlingToCurrentValueRange();
    tmp.savePNG ( aol::strprintf( "%s_lowPKernelFT.png", outputFileNameBase.c_str() ).c_str() );
    tmp2.setOverflowHandlingToCurrentValueRange();
    tmp2.savePNG ( aol::strprintf( "%s_highPKernelFT.png", outputFileNameBase.c_str() ).c_str() );

    lowPass[0].shiftByOffset ( grid.getNumX() / 2, grid.getNumY() / 2, tmp );
    highPass[0].shiftByOffset ( grid.getNumX() / 2, grid.getNumY() / 2, tmp2 );

    qc::FourierTransform ( image, imageTransformed, qc::FTForward );

    complexMultiply ( imageTransformed, lowPass, imageMultiplied );
    qc::FourierTransform ( imageMultiplied, image, qc::FTBackward );
    // The fourier transform is not unnormalized.
    image[0] /= grid.getNumberOfNodes();

    image[0].setOverflowHandlingToCurrentValueRange();
    image[0].savePNG ( aol::strprintf( "%s_low.png", outputFileNameBase.c_str() ).c_str() );

    // Undersampling test
    saveUndersampledImage ( a, undersampleFactor, outputFileNameBase.c_str() );
    saveUndersampledImage ( image[0], undersampleFactor, aol::strprintf( "%s_low", outputFileNameBase.c_str() ).c_str() );

    image[0] -= a;
    image[0].setOverflowHandlingToCurrentValueRange();
    image[0].savePNG ( aol::strprintf( "%s_highPart.png", outputFileNameBase.c_str() ).c_str() );

    complexMultiply ( imageTransformed, highPass, imageMultiplied );
    qc::FourierTransform ( imageMultiplied, image, qc::FTBackward );
    image[0].setOverflowHandlingToCurrentValueRange();
    image[0].savePNG ( aol::strprintf( "%s_high.png", outputFileNameBase.c_str() ).c_str() );

    qc::FourierTransform ( lowPass, imageMultiplied, qc::FTBackward );
    saveComplexKernel ( imageMultiplied, aol::strprintf( "%s_lowPKernel.png", outputFileNameBase.c_str() ).c_str() );
    qc::FourierTransform ( highPass, imageMultiplied, qc::FTBackward );
    saveComplexKernel ( imageMultiplied, aol::strprintf( "%s_highPKernel.png", outputFileNameBase.c_str() ).c_str() );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
