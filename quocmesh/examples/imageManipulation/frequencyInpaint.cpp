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
 * \brief A very simple algorithm to inpaint a bandwidth limited image using the frequency domain.
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
      cerr << "USAGE: " << argv[0] << "  <input_file> <lowPassRadius>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const RType lowPassRadius = atof ( argv[2] );
    const string outputFileNameBase = aol::strprintf( "%s_%.3f", aol::getBaseFileName ( inputFileName ).c_str(), lowPassRadius );
    const qc::ScalarArray<RType, qc::QC_2D> a ( inputFileName );

    qc::BitArray<qc::QC_2D> mask ( a.getSize() );
    for ( int i = 3 * a.getNumX() / 8; i < 5 * a.getNumX() / 8; ++i )
      for ( int j = 3 * a.getNumY() / 8; j < 5 * a.getNumY() / 8; ++j )
        mask.set ( i, j, true );

    ConfType::InitType grid ( qc::GridSize<ConfType::Dim>::createFrom ( a ) );
    qc::MultiArray<RType, 2, 2> image ( grid ), imageTransformed ( grid ), imageMultiplied ( grid ), lowPass ( grid );
    qc::ScalarArray<RType, qc::QC_2D> tmp ( a, aol::STRUCT_COPY_UNINIT );

    for ( int i = 0; i < a.size(); ++i )
      if ( mask[i] == false )
        image[0][i] = a[i];

    image[0].savePNG ( aol::strprintf( "%s_inpaintFreqInput.png", outputFileNameBase.c_str() ).c_str() );

    qc::DataGenerator<ConfType> generator ( grid );
    const RType centerX = 0.5 * grid.H() * ( grid.getNumX() - 1 );
    const RType centerY = 0.5 * grid.H() * ( grid.getNumY() - 1 );
    generator.generateCircle( tmp, 1, lowPassRadius, centerX, centerY );

    tmp.setOverflowHandlingToCurrentValueRange();
    tmp.savePNG ( aol::strprintf( "%s_lowPKernelFT.png", outputFileNameBase.c_str() ).c_str() );

    lowPass[0].shiftByOffset ( grid.getNumX() / 2, grid.getNumY() / 2, tmp );

    for ( int iter = 0; iter < 5000; ++iter ) {
      qc::FourierTransform ( image, imageTransformed, qc::FTForward );
      complexMultiply ( imageTransformed, lowPass, imageMultiplied );
      qc::FourierTransform ( imageMultiplied, image, qc::FTBackward );
      // The fourier transform is not unnormalized.
      image[0] /= grid.getNumberOfNodes();

      for ( int i = 0; i < a.size(); ++i )
        if ( mask[i] == false )
          image[0][i] = a[i];

      image[0].clamp(0, 255);
      image[1].setZero();

      if ( ( ( iter + 1 ) % 100 ) == 0 )
        image[0].savePNG ( aol::strprintf( "%s_inpaintFreq_%04d.png", outputFileNameBase.c_str(), iter + 1 ).c_str() );
    }
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
