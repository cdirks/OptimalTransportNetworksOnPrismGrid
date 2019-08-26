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
 *  \brief  comparison of PPM images
 *
 *  Compares two ppm images, computes number of different pixels and
 *  sum of the absolute of the difference for each color channel and
 *  saves the difference image.
 *
 *  Usage: comparePPMImages inimage1 inimage2 outdiffimage
 *
 *  \author Schwen
 */

#include <iostream>
#include "scalarArray.h"

int main ( int argc, char **argv ) {
  if ( argc == 4 ) {
    FILE * indat1_stream;
    FILE * indat2_stream;
    FILE * outdat_stream;

    if ( ! ( indat1_stream = fopen ( argv[1], "r" ) ) ) cerr << "Cannot open file " << argv[1] << " for reading.";
    if ( ! ( indat2_stream = fopen ( argv[2], "r" ) ) ) cerr << "Cannot open file " << argv[2] << " for reading.";
    if ( ! ( outdat_stream = fopen ( argv[3], "w" ) ) ) cerr << "Cannot open file " << argv[3] << " for writing.";

    char tmp[2][256];
    fscanf ( indat1_stream, "%s", tmp[0] );
    fscanf ( indat2_stream, "%s", tmp[1] );

    if ( tmp[0][0] != 'P' || tmp[1][0] != 'P' ) {
      cerr << "Unknown file header." << endl;
      return ( -2 );
    }
    if ( tmp[0][1] != '3' || tmp[1][1] != '3' ) {
      cerr << "Unknown magic number." << endl;
      return ( -3 );
    }

    int wi, he, tmpn;
    fscanf ( indat1_stream, "%d %d\n%d\n", &wi, &he, &tmpn );
    int wi2, he2, tmpn2;
    fscanf ( indat2_stream, "%d %d\n%d\n", &wi2, &he2, &tmpn2 );

    if ( ( wi != wi2 ) || ( he != he2 ) ) {
      cerr << "Incompatible image sizes" << endl;
      return ( -4 );
    }
    if ( ( tmpn != 255 ) || ( tmpn2 != 255 ) ) {
      cerr << "Max intensity is different from 255." << endl;
    }


    fprintf ( outdat_stream, "P3\n%d %d\n255\n", wi, he );

    qc::ScalarArray<int, qc::QC_2D> image1r ( wi, he ), image1g ( wi, he ), image1b ( wi, he ), image2r ( wi, he ), image2g ( wi, he ), image2b ( wi, he );

    int tmp1r, tmp1g, tmp1b, tmp2r, tmp2g, tmp2b;
    for ( int jy = 0; jy < he; jy++ ) {
      for ( int ix = 0; ix < wi; ix++ ) {
        fscanf ( indat1_stream, "%d %d %d", &tmp1r, &tmp1g, &tmp1b );
        fscanf ( indat2_stream, "%d %d %d", &tmp2r, &tmp2g, &tmp2b );
        image1r.set ( ix, jy, tmp1r );
        image1g.set ( ix, jy, tmp1g );
        image1b.set ( ix, jy, tmp1b );
        image2r.set ( ix, jy, tmp2r );
        image2g.set ( ix, jy, tmp2g );
        image2b.set ( ix, jy, tmp2b );
      }
    }

    int different_pixelsr = 0, different_pixelsg = 0, different_pixelsb = 0,  diffsumr = 0, diffsumg = 0, diffsumb = 0;

    for ( int i = 0; i < image1r.size(); i++ ) {
      if ( image1r[i] != image2r[i] ) {
        different_pixelsr++;
        diffsumr += abs ( image1r[i] - image2r[i] );
      }
      fprintf ( outdat_stream, "%d ", abs ( image1r[i] - image2r[i] ) );
      if ( image1g[i] != image2g[i] ) {
        different_pixelsg++;
        diffsumg += abs ( image1g[i] - image2g[i] );
      }
      fprintf ( outdat_stream, "%d ", abs ( image1g[i] - image2g[i] ) );
      if ( image1b[i] != image2b[i] ) {
        different_pixelsb++;
        diffsumb += abs ( image1b[i] - image2b[i] );
      }
      fprintf ( outdat_stream, "%d ", abs ( image1b[i] - image2b[i] ) );

      fprintf ( outdat_stream, "\n" );
    }

    fclose ( indat1_stream );
    fclose ( indat2_stream );
    fclose ( outdat_stream );

    cout << "In total, " << different_pixelsr << " red pixels were different, that is "
    << ( 100.0 * different_pixelsr ) / image1r.size() << " percent." << endl;
    cout << "The total difference is " << diffsumr << ", the average difference is "
    << ( 1.0 * diffsumr ) / different_pixelsr << " among these, "
    << ( 1.0 * diffsumr ) / image1r.size() << " for all red pixels." << endl;

    cout << "In total, " << different_pixelsg << " green pixels were different, that is "
    << ( 100.0 * different_pixelsg ) / image1g.size() << " percent." << endl;
    cout << "The total difference is " << diffsumg << ", the average difference is "
    << ( 1.0 * diffsumg ) / different_pixelsg << " among these, "
    << ( 1.0 * diffsumg ) / image1g.size() << " for all green pixels." << endl;

    cout << "In total, " << different_pixelsb << " blue pixels were different, that is "
    << ( 100.0 * different_pixelsb ) / image1b.size() << " percent." << endl;
    cout << "The total difference is " << diffsumb << ", the average difference is "
    << ( 1.0 * diffsumb ) / different_pixelsb << " among these, "
    << ( 1.0 * diffsumb ) / image1b.size() << " for all blue pixels." << endl;

    return ( 0 );

  } else {
    cerr << "usage: " << endl << "comparePPMImages image_1.ppm image_2.ppm output_difference_image.ppm" << endl;
    return ( 1 );
  }

}
