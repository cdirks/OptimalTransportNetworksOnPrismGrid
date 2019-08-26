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

#include <qmException.h>
#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>


void drawLocalEllipsoid(qc::ScalarArray<double, qc::QC_3D> &Data, int N, int x, int y, int z, double a, double b, double c)
{
  int A = static_cast<int>(a);
  int B = static_cast<int>(b);
  int C = static_cast<int>(c);

  for (int i=x-A; i<=x+A; i++)
    for (int j=y-B; j<=y+B; j++)
      for (int k=z-C; k<=z+C; k++)
        if (i >= 0 && i < N &&  j >= 0 && j < N && k >= 0 && k < N )
          if ( aol::Sqr( (x-i)/a ) + aol::Sqr( (y-j)/b ) + aol::Sqr( (z-k)/c ) <= 1. )
            Data.set(i, j, k, 1.);
}

int main( int, char** ) {

  try {
    // generate different ellipsoids and save them
    qc::GridDefinition grid( 5, qc::QC_3D );
    qc::ScalarArray<double, qc::QC_3D> ellipsoid( grid );
    ellipsoid.setZero();
    drawLocalEllipsoid( ellipsoid, 33, 17, 17, 17, 5, 8, 12 );
    ellipsoid.save( "./image3D_1.bz2", qc::PGM_DOUBLE_BINARY );
    ellipsoid.setZero();
    drawLocalEllipsoid( ellipsoid, 33, 17, 17, 17, 8, 5, 12 );
    ellipsoid.save( "./image3D_2.bz2", qc::PGM_DOUBLE_BINARY );
    ellipsoid.setZero();
    drawLocalEllipsoid( ellipsoid, 33, 17, 17, 17, 5, 8, 15 );
    ellipsoid.save( "./image3D_3.bz2", qc::PGM_DOUBLE_BINARY );
  }
  catch ( aol::Exception &el ) {
    // print out error message
    el.dump();
  }

  // pause before closing the program (to let the user first read the output)
  aol::callSystemPauseIfNecessaryOnPlatform();
  return EXIT_SUCCESS;
}
