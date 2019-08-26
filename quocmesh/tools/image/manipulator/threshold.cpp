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
 *  \brief Sets all values smaller than $a$ and bigger than $b$ to 0.
 *
 *  Sets all values smaller than $a$ and bigger than $b$ to 0
 *  (works both in 2d and 3d).
 *  Usage: threshold inputImg outputImg a b [only01] [scale]
 *  Option [only01] sets the image to 0 for values smaller/bigger than a/b, and 1 else. \\
 *  Option [scale] divides the input image by its max value. \\
 *  {\bf Attention:} The output-image will only be saved with double values if the suffix is '.bz2'!  \\
 *
 *  \author Oliver Nemitz
 */

#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>
#include <scalarArray.h>
#include <string.h>

#include <qmException.h>
#include <fstream>
#include <auxiliary.h>

typedef double RealType;

using namespace aol::color;

// get the first character of the header
char getFirstHeaderChar ( char* fileName ) {
  aol::ipfstream file ( fileName );
  qc::ArrayHeader header;
  qc::ReadArrayHeader ( file, header );
  file.close();

  return header.magic[0];
}

int main ( int argc, char **argv ) {
  try {
    if ( argc < 5 ) {
      cerr << "Sets all values smaller than a and bigger than b to 0.\n";
      cerr << "Option [only01] sets the image to 0 for values smaller/bigger than a/b, and 1 else.\n";
      cerr << "Option [scale] divides the input image by its max value.\n";
      cerr << "Attention: If the suffix of the output is '.bz2', it will be saved with double values, otherwise not!\n" << endl;
      cerr << aol::color::red << "usage: " << argv[0] << " input output a b [only01] [scale]\n\n" << aol::color::black;
      return EXIT_FAILURE;
    }

    double a = atof ( argv[3] );
    double b = atof ( argv[4] );

    // should the image consist only of 0 and 1
    bool only01 = false;
    bool scale = false;

    for ( int i = 5; i < argc; i++ ) {
      if ( !strcmp ( argv[i], "only01" ) && !only01 ) {
        only01 = true;
        cerr << red << "Saving only values 0 and 1...\n";
      }
      if ( !strcmp ( argv[i], "scale" ) && !scale ) {
        scale = true;
        cerr << red << "Dividing input data by max value...\n";
      }
    }

    char c = getFirstHeaderChar ( argv[1] );

    // ------------------------- --------------------------------------------
    // ------------------------- 3D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( c == 'Q' ) {
      qc::ScalarArray<RealType, qc::QC_3D> inputImg ( argv[1] );
      qc::ScalarArray<RealType, qc::QC_3D> outputImg ( inputImg );
      cerr << blue << "3D: Thresholding every value lower than " << red << a << blue << " and bigger than " << red << b << blue << "...";

      for ( int i = 0; i < inputImg.size(); i++ ) {
        if ( inputImg[i] < a || inputImg[i] > b ) outputImg[i] = 0.;
        else if ( only01 ) outputImg[i] = 1.;
      }

      cerr << "done! Saving...\n" << aol::color::black;
      outputImg.save ( argv[2], qc::SaveTypeTrait<RealType>::BinarySaveType );
    }

    // ------------------------- --------------------------------------------
    // ------------------------- 2D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( c == 'P' ) {
      qc::ScalarArray<RealType, qc::QC_2D> inputImg ( argv[1] );
      if ( scale )
        inputImg /= inputImg.getMaxValue();
      qc::ScalarArray<RealType, qc::QC_2D> outputImg ( inputImg );
      cerr << blue << "2D: Thresholding every value lower than " << red << a << blue << " and bigger than " << red << b << blue << "...";

      for ( int i = 0; i < inputImg.size(); i++ ) {
        if ( inputImg[i] < a || inputImg[i] > b ) outputImg[i] = 0.;
        else if ( only01 ) outputImg[i] = 1.;
      }

      cerr << "done! Saving...\n" << aol::color::black;
      qc::recognizeEndingAndSave2d ( outputImg, argv[2], only01 || scale );
    }

    cerr << aol::color::blue << "done! Thanx for using threshold, THE standard utility for thresholding QuocMesh-data!\n";
    cerr << aol::color::black;
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
