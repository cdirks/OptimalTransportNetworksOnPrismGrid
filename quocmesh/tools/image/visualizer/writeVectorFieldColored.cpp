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
 * \brief Writes a vector field as colored image: The color encodes the vector direction, the brightness the vector length.
 *
 * Usage: writeVectorFieldColored <inputFileDefX> <inputFileDefY> [outputFileBasename] [subtractAverage]
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>
#include <colorWheel.h>

template <typename RealType, qc::Dimension Dim>
RealType possiblySubtractAverageAndGetMaxOfPointWiseNorm ( qc::ScalarArray<RealType, Dim> &DefX, qc::ScalarArray<RealType, Dim> &DefY, const bool SubtractAverage ) {
  if ( SubtractAverage ) {
    DefX.addToAll ( -DefX.getMeanValue() );
    DefY.addToAll ( -DefY.getMeanValue() );
  }
  qc::MultiArray<RealType, Dim, 2> defArray ( qc::GridSize<Dim> ( DefX ), qc::MultiArray<RealType, Dim, 2>::CREATE_NO_ARRAYS );
  defArray.appendReference ( DefX );
  defArray.appendReference ( DefY );
  return defArray.getMaxOfPointWiseNorm();
}

typedef double RType;
int main ( int argc, char **argv ) {

  try {
    if ( argc < 3 ) {
      cerr << "USAGE: " << argv[0] << "  <inputFileDefX> <inputFileDefY> [outputFileBasename] [subtractAverage]" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const qc::Dimension dim = qc::getDimensionFromArrayFile ( argv[1] );
    const bool subtractAverage = ( argc >= 5 ) ? !!( atoi ( argv[4] ) ) : false;

    if ( dim == qc::QC_3D ) {
      qc::ScalarArray<RType, qc::QC_3D> defX ( argv[1] );
      qc::ScalarArray<RType, qc::QC_3D> defY ( argv[2] );
      const RType maxNorm = possiblySubtractAverageAndGetMaxOfPointWiseNorm<RType, qc::QC_3D> ( defX, defY, subtractAverage );

      const int numX = defX.getNumX();
      const int numY = defX.getNumY();
      const int numSlices = defX.getNumZ();

      qc::ScalarArray<RType, qc::QC_2D> defXSlice ( numX, numY );
      qc::ScalarArray<RType, qc::QC_2D> defYSlice ( numX, numY );

      for ( int sliceNum = 0; sliceNum < numSlices; ++sliceNum ) {
        defX.getSlice ( qc::QC_Z, sliceNum, defXSlice );
        defY.getSlice ( qc::QC_Z, sliceNum, defYSlice );

        string filename = aol::strprintf ( "%s_%d.png", ( argc > 3 ) ? argv[3] : "output", sliceNum );
        qc::writeColorField ( defXSlice, defYSlice, filename.c_str(), maxNorm );
      }
    }
    else if ( dim == qc::QC_2D ) {
      qc::ScalarArray<RType, qc::QC_2D> defX ( argv[1] );
      qc::ScalarArray<RType, qc::QC_2D> defY ( argv[2] );
      const RType maxNorm = possiblySubtractAverageAndGetMaxOfPointWiseNorm<RType, qc::QC_2D> ( defX, defY, subtractAverage );
      const RType maxNormPX = maxNorm * ( defX.getSize().getMaxValue() - 1 );
      qc::writeColorField ( defX, defY,  aol::strprintf ( "%s_%fpx.png", ( argc > 3 ) ? argv[3] : "output", maxNormPX ) );
      cerr << "maxNorm = " << maxNorm << " ( " << maxNormPX << " pixels)\n";
    }

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
