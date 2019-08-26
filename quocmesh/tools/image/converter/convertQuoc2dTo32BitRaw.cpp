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
 * \brief Tool to load 2D quoc data and save them as unsigned int raw. Assumes that the values are already in the unsigned int range.
 *
 * \author Berkels
 */

#include <scalarArray.h>

int main ( int argc, char **argv ) {

  try {
    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <output_file>" << endl;
      return ( EXIT_FAILURE );
    }

    const string inputFileName = argv[1];
    const string outputFileName = argv[2];

    qc::ScalarArray<uint32_t, qc::QC_2D> a ( inputFileName.c_str() );
    aol::Bzipofstream out ( outputFileName.c_str() );
    out.write ( reinterpret_cast<char*> ( a.getData() ), a.size() * sizeof ( uint32_t ) );
  } catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return ( EXIT_SUCCESS );
}
