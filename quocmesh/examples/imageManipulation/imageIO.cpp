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

/** \file
 *  \brief various methods for loading and saving 2D and 3D image data
 *
 *  \author Schwen
 */

#include<scalarArray.h>

int main ( int, char** ) {
  try {

    // Loading Data
    // ------------

    // We can load data in the constructor of ScalarArrays ...
    qc::ScalarArray< unsigned char, qc::QC_2D > image2D ( "../testdata/image_129.pgm" );

    // ... or by first creating an appropriate ScalarArray and then calling load:
    qc::GridDefinition grid3D ( 4 /* grid level */, qc::QC_3D );
    qc::ScalarArray<double, qc::QC_3D> image3D ( grid3D );
    image3D.load ( "../testdata/volume_17.dat.bz2" );


    // Saving Data
    // -----------

    // Saving works via save methods:
    image2D.save ( "filename.pgm", qc::PGM_UNSIGNED_CHAR_ASCII, "This is a comment." ); // for human-readable pgm graphics, there are other save types available.
#ifdef USE_LIB_PNG
    image2D.savePNG ( "filename.png" ); // save in png format (surprise!)
#endif

    image3D.save ( "filename.dat.bz2", qc::PGM_DOUBLE_BINARY ); // automatically compresses 3D output in binary double format


    // Treating Raw Data
    // -----------------

    // Suppose we're given raw data (without quoc-specific header) by someone else:
    aol::ipfstream instr ( "../testdata/volume_17.raw" );
    image3D.loadRaw ( instr, qc::PGM_DOUBLE_BINARY, 17, 17, 17 );

    // then maybe we need to fix the endianness if this data was created on a different platform
    // not necessary here, so we'll swap twice
    image3D.swapByteOrder(); image3D.swapByteOrder();

    // If we want to save raw data:
    image3D.saveRaw ( "filename.raw", qc::PGM_DOUBLE_BINARY, 0, 0 );


    // Interaction Between Data of Different Size / Dimension
    // ------------------------------------------------------

    // pad image3D into other array (center to center)
    qc::ScalarArray<double, qc::QC_3D> padImage3D ( 15, 17, 20 );
    padImage3D.padFrom ( image3D );

    // extract 2D slice out of 3D dataset
    qc::ScalarArray<double, qc::QC_2D> slice2D ( image3D.getNumX(), image3D.getNumY() );
    image3D.getSlice ( qc::QC_X, 7, slice2D );

    // import 2D slice into 3D dataset
    image3D.putSlice ( qc::QC_Y, 4, slice2D );

    // save slices of 3D dataset in multiple images
    image3D.addToAll ( -image3D.getMinValue() ); image3D *= 255 / image3D.getMaxValue(); // change contrast so that values range from 0 to 255
    image3D.saveSlices( "z_slices_%03d.pgm", qc::QC_Z, qc::PGM_UNSIGNED_CHAR_ASCII );

    // load slices into 3D dataset from multiple images
    image3D.setZero();
    image3D.loadSlices( "z_slices_%03d.pgm", qc::QC_Z, 5, 10 ); // loading only subset results in slices being centered

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
