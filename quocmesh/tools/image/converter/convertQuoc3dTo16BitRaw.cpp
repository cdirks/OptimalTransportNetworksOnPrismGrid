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
 * \brief Tool to load 3D quoc data and save them as unsigned short raw
 *
 * \author Schwen
 */

#include <scalarArray.h>

typedef double RType;

int main ( int argc, char **argv ) {

  try {
    char inputFileName[1024];
    char outputFileName[1024];

    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <output_file>" << endl;
      return ( EXIT_FAILURE );
    }

    sprintf ( inputFileName, "%s",  argv[1] );
    sprintf ( outputFileName, "%s",  argv[2] );

    qc::ScalarArray<RType, qc::QC_3D> a ( inputFileName );
    a.saveRaw ( outputFileName, qc::PGM_UNSIGNED_SHORT_BINARY, 0, a.getMaxAbsValue() );

  } catch ( aol::Exception &el ) {
    el.dump();
  }
  return ( EXIT_SUCCESS );
}
