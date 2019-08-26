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
 * \brief Clamps a given 2D or 3D ScalarArray into [min,max] and then scales the values to [0,1].
 *
 * \author Berkels
 */

#include <aol.h>
#include <FEOpInterface.h>
#include <scalarArray.h>

template <typename RealType, qc::Dimension Dim>
void clampAndSave ( const char* Filename, const RealType Min, const RealType Max ) {
  qc::ScalarArray<RealType, Dim> a ( Filename );
  a.clamp ( Min, Max );
  a.addToAll ( -Min );
  a /= ( Max - Min );
  const string outputFileName = aol::strprintf ( "%s_clamped_%f_%f.raw", Filename, Min, Max );
  a.save( outputFileName.c_str(), qc::PGM_FLOAT_BINARY );
}

typedef float RType;
int main ( int argc, char **argv ) {

  try {
    char inputFileName[1024];

    if ( argc != 4 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <min> <max>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }
    sprintf ( inputFileName, "%s",  argv[1] );
    const RType min = atof ( argv[2] );
    const RType max = atof ( argv[3] );

    cerr << "Reading from " << inputFileName << endl;

    const qc::Dimension dim = qc::getDimensionFromArrayFile ( inputFileName );
    if ( dim == qc::QC_3D )
      clampAndSave<RType, qc::QC_3D> ( inputFileName, min, max );
    else if ( dim == qc::QC_2D )
      clampAndSave<RType, qc::QC_2D> ( inputFileName, min, max );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
