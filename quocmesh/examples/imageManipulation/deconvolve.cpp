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
 * \brief Convolves an image with a Gaussian and then deconvolves the result using the Fourier convolution theorem.
 *
 * \author Berkels
 */

#include <scalarArray.h>
#include <configurators.h>
#include <generator.h>
#include <convolution.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

template <typename RealType>
void complexDivide ( const RealType ARe, const RealType AIm, const RealType BRe, const RealType BIm, RealType &ResRe, RealType &ResIm, const RealType Epsilon ) {
  const RealType a = ARe;
  const RealType b = AIm;
  const RealType c = BRe;
  const RealType d = BIm;
  const RealType scale = 1 / ( aol::Sqr ( c ) + aol::Sqr( d ) + Epsilon );
  ResRe = scale * ( a * c + b * d );
  ResIm = scale * ( b * c - a * d );
}

template <typename RealType>
void complexDivide ( const qc::MultiArray<RealType, 2, 2> &ImageA, const qc::MultiArray<RealType, 2, 2> &ImageB, qc::MultiArray<RealType, 2, 2> &Result, const RealType Epsilon ) {
  const int numPixels = ImageA.getEqualComponentSize();
  for ( int i = 0; i < numPixels; i++ )
    complexDivide ( ImageA[0][i], ImageA[1][i], ImageB[0][i], ImageB[1][i], Result[0][i], Result[1][i], Epsilon );
}

int main ( int argc, char **argv ) {
  try {
    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const string outputFileNameBase = aol::getBaseFileName ( inputFileName );
    const qc::ScalarArray<RType, qc::QC_2D> a ( inputFileName );

    const ConfType::InitType grid ( qc::GridSize<ConfType::Dim>::createFrom ( a ) );
    qc::MultiArray<RType, 2, 2> image ( grid ), imageTransformed ( grid ), imageMultiplied ( grid ), gaussKernel ( grid ), gaussKernelTransformed ( grid );
    qc::ScalarArray<RType, qc::QC_2D> tmp ( a, aol::STRUCT_COPY_UNINIT );

    qc::generateGaussKernel ( 0.005, gaussKernel[0] );
    gaussKernel[0].setOverflowHandlingToCurrentValueRange();
    gaussKernel[0].savePNG ( aol::strprintf( "%s_gaussKernel.png", outputFileNameBase.c_str() ).c_str() );
    qc::Convolution<qc::QC_2D, RType> conv ( qc::GridSize2d ( grid ).getSizeAsVecDim() );
    conv.convolve ( a, gaussKernel[0], tmp );
    tmp.setOverflowHandlingToCurrentValueRange();
    tmp.savePNG ( aol::strprintf( "%s_blurred.png", outputFileNameBase.c_str() ).c_str() );

    for ( int k = 1; k < 10; ++k ) {
      image[0] = tmp;
      image[1].setZero();
      qc::FourierTransform ( image, imageTransformed, qc::FTForward );
      qc::FourierTransform ( gaussKernel, gaussKernelTransformed, qc::FTForward );
      complexDivide ( imageTransformed, gaussKernelTransformed, imageMultiplied, std::pow ( 0.1, k ) );
      qc::FourierTransform ( imageMultiplied, image, qc::FTBackward );
      image[0].setOverflowHandlingToCurrentValueRange();
      image[0].savePNG ( aol::strprintf( "%s_deconv%d.png", outputFileNameBase.c_str(), k ).c_str() );

      image[0] = tmp;
      image[0].scaleValuesToAB( 0, 255 );
      for ( int i = 0; i < grid.getNumberOfNodes(); ++i )
        image[0][i] = static_cast<unsigned char> ( image[0][i] );

      if ( k == 1 ) {
        image[0].setOverflowHandlingToCurrentValueRange();
        image[0].savePNG ( aol::strprintf( "%s_quant.png", outputFileNameBase.c_str() ).c_str() );
      }
      image[1].setZero();
      qc::FourierTransform ( image, imageTransformed, qc::FTForward );
      complexDivide ( imageTransformed, gaussKernelTransformed, imageMultiplied, std::pow ( 0.1, k ) );
      qc::FourierTransform ( imageMultiplied, image, qc::FTBackward );
      image[0] /= grid.getNumberOfNodes();
      image[0].setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 255 );
      image[0].savePNG ( aol::strprintf( "%s_quantDeconv%d.png", outputFileNameBase.c_str(), k ).c_str() );
    }

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
