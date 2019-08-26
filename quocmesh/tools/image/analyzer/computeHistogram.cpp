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
 *
 *  \brief Loads a 3D dataset and computes a histogram of values in it.
 *
 *  Usage: \code ./computeHistogram image min_val max_val steps \endcode
 *
 *  The histogram divides the interval from \code min_val \endcode to \code max_val \endcode into \code steps \endcode subintervals
 *  and counts how many values fall into the respective subintervals. Values below \code min_val \endcode or above \code max_val \endcode are not counted
 *  (so that the sum of the entries in the histogram may be less than the number of voxels).
 *
 *  \author Schwen
 */

#include <scalarArray.h>
#include <qmException.h>

int main ( int argc, char **argv ) {
  try {

    if ( argc != 5 ) {
      cerr << "Usage: " << argv[0] << " image min_val max_val steps" << endl;
      return ( EXIT_FAILURE );
    } else {

      char filename[1024];
      strncpy ( filename, argv[1], 1024 );

      const double
        min_val = atof ( argv[2] ),
        max_val = atof ( argv[3] );

      const int
      steps = atoi ( argv[4] );

      qc::ScalarArray<double, qc::QC_3D> img ( filename );

      aol::Vector<int> histo ( steps );
      histo.setZero();

      for ( int i = 0; i < img.size(); ++i ) {
        const double rounded_down = floor ( ( img[i] - min_val ) * steps / ( max_val - min_val ) );
        const int r_d = static_cast<int> ( rounded_down );
        if ( r_d >= 0 && r_d < steps ) {
          ++ ( histo[ r_d ] );
        }
      }

      cerr << "Min value = " << img.getMinValue() << ", max value = " << img.getMaxValue() << ", sum = " << img.sum() << endl;

      for ( int i = 0; i < steps; ++i ) {
        cout << min_val + i * ( max_val - min_val ) / steps << " " << histo[i] << endl;
      }

    }

    return ( EXIT_SUCCESS  ) ;

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
}
