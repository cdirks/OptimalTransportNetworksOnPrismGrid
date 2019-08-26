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
 * \brief Example for downsampling and upsampling an image using RestrictOp and ProlongOp.
 *
 *
 *
 * \author Nemitz
 */

#include <quoc.h>
#include <qmException.h>
#include <aol.h>

#include <restriction.h>
#include <prolongation.h>
#include <scalarArray.h>

int main( int, char** ) {
  try {
    // load the fine image
    qc::GridDefinition grid( 7, qc::QC_2D );
    qc::ScalarArray<double, qc::QC_2D> SourceImg( grid );
    SourceImg.load( "../testdata/image_129.pgm" );

    // declare the coarse image
    qc::GridDefinition gridCoarse( 6, qc::QC_2D );
    qc::ScalarArray<double, qc::QC_2D> CoarseImg( gridCoarse );

    // downsampling
    qc::RestrictOp<double, qc::STD_QUOC_RESTRICT> Restrict( gridCoarse, grid );
    Restrict.apply( SourceImg, CoarseImg );

    // save the modified picture
    CoarseImg.save( "image_downsampled.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );

    // declare the fine image
    qc::GridDefinition gridFine( 8, qc::QC_2D );
    qc::ScalarArray<double, qc::QC_2D> FineImg( gridFine );

    // upsampling
    qc::ProlongOp<double> Prolong( grid, gridFine );
    Prolong.apply( SourceImg, FineImg );

    // save the modified picture
    FineImg.save( "image_upsampled.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
