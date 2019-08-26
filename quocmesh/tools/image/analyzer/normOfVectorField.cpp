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
 * \brief Calculates the pointwise norm of a 2D vector field and saves the result as 2D array.
 *
 * Usage: normOfVectorField <inputFileDefX> <inputFileDefY> [outputFile]
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>
#include <colorWheel.h>

typedef double RType;
int main ( int argc, char **argv ) {

  try {
    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <inputFileDefX> <inputFileDefY>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    qc::ScalarArray<RType, qc::QC_2D> defX ( argv[1] );
    qc::ScalarArray<RType, qc::QC_2D> defY ( argv[2] );

    qc::MultiArray<RType, 2, 2> defArray ( qc::GridSize<qc::QC_2D> ( defX ), qc::MultiArray<RType, 2, 2>::CREATE_NO_ARRAYS );
    defArray.appendReference ( defX );
    defArray.appendReference ( defY );

    qc::ScalarArray<RType, qc::QC_2D> norm ( defX, aol::STRUCT_COPY );
    defArray.getPointWiseNorm( norm );
    norm.save ( aol::strprintf ( "%s-norm.dat.bz2", argv[1] ).c_str(), qc::PGM_DOUBLE_BINARY );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
