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
 *  \brief Scales the input-image to $[0,1]$ (works in 2d and in 3d).
 *
 *  Scales the input-image to $[0,1]$ (works in 2d and in 3d).
 *  Usage: scaleImg01 inputImg outputImg
 *  {\bf Attention:} The output-image will only be saved with double values if the suffix is '.bz2'!
 *
 *  \author Nemitz
 */

#include <quoc.h>
#include <scalarArray.h>
#include <scalarArray.h>

#include <qmException.h>
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

/**
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void scaleImage ( const char *InputFilename, const char *OutputFilename ) {
  qc::ScalarArray<RealType, Dim> img ( InputFilename );

  const RealType min = img.getMinValue();
  const RealType max = img.getMaxValue();
  img.addToAll ( - min );
  cerr << blue << "Minimum value in this ";
  if ( Dim == qc::QC_2D )
    cerr << "2";
  else
    cerr << "3";
  cerr << "d-volume was " << red << min << blue << ", max value was " << red << max << blue << ".\n" << reset;
  img /= ( max - min );

  img.save ( OutputFilename, qc::PGM_DOUBLE_BINARY );
}

int main ( int argc, char **argv ) {
  try {
    if ( argc < 3 ) {
      cerr << "Scales image to [0,1] (either in 2d or in 3d).\n";
      cerr << aol::color::red << "usage: " << argv[0] << " input output\n" << aol::color::black;
      return EXIT_FAILURE;
    }

    const char c = getFirstHeaderChar ( argv[1] );


    // ------------------------- --------------------------------------------
    // ------------------------- 3D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( c == 'Q' )
      scaleImage<double, qc::QC_3D> ( argv[1], argv[2] );

    // ------------------------- --------------------------------------------
    // ------------------------- 2D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( c == 'P' )
      scaleImage<double, qc::QC_2D> ( argv[1], argv[2] );

    cerr << aol::color::blue << "done!\n";
    cerr << aol::color::black;
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
