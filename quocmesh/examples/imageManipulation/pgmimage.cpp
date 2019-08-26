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
 * \brief This example show how to read and write PGM (grey value) pictures and how to access its grey values at given pixel coordinates.
 *
 *
 *
 * \author Lenz
 */


#include <scalarArray.h>

int main ( int, char** ) {
  try {
    const int n = 129;

    qc::ScalarArray<double, qc::QC_2D> picture ( n, n );

    picture.load ( "../testdata/image_129.pgm" );

    cerr << picture.get ( 23, 42 ) << endl;

    for ( int i = 50; i < 100 ; i++ ) {
      for ( int j = 50; j < 100 ; j++ ) {
        picture.set ( i, j, 250 );
      }
    }

    cerr << picture.get ( 70, 70 ) << endl;

    picture.save ( "erg.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );

  } catch ( aol::Exception& ex ) {
    ex.dump();
  }
}
