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
 * \brief Example for the convolution of a 2d-image with a Gaussian kernel.
 *
 *
 *
 * \author Nemitz
 */

#include <quoc.h>
#include <qmException.h>
#include <aol.h>

#include <kernel2d.h>
#include <scalarArray.h>

int main ( int, char** ) {
  try {

    // declaring Gaussian kernels (makeKernel is called in the constructor)
    qc::GaussKernel2d<double> KernelSize3 ( 3, 1. );
    qc::GaussKernel2d<double> KernelSize9 ( 9, 1. );
    qc::GaussDiffKernel2d<double> DiffKernel ( 3, 1., qc::DIFF_X );

    // load the image
    qc::GridDefinition grid ( 7, qc::QC_2D );
    qc::ScalarArray<double, qc::QC_2D> SourceImg ( grid );
    SourceImg.load ( "../testdata/image_129.pgm.bz2" );
    qc::ScalarArray<double, qc::QC_2D> DestImg ( SourceImg, aol::DEEP_COPY );

    // now convolve parts of the image with the single kernels

    // 1. Gauss-Kernel with size 3
    for ( int i = 20; i < 85; i++ )
      for ( int j = 30; j < 95; j++ )
        DestImg.set ( i, j, SourceImg.getConvolveValue ( i, j, KernelSize3 ) );

    DestImg.save ( "blurred_3.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );

    // 2. Gauss-Kernel with size 9
    for ( int i = 20; i < 85; i++ )
      for ( int j = 30; j < 95; j++ )
        DestImg.set ( i, j, SourceImg.getConvolveValue ( i, j, KernelSize9 ) );

    DestImg.save ( "blurred_9.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );

  // 3. GaussDiffKernel
    for ( int i = 20; i < 85; i++ )
      for ( int j = 30; j < 95; j++ )
        DestImg.set ( i, j, SourceImg.getConvolveValue ( i, j, DiffKernel ) );

    DestImg.save ( "blurred_diff.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );


  } catch ( aol::Exception &ex ) {

    ex.dump();
    return ( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );
}











