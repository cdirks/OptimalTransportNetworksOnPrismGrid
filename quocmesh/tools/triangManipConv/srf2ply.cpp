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
 *  \brief Converts a BrainVoyager SRF file to PLY
 *
 *  Usage: srf2ply sourceFile.srf
 *
 *  \author Berkels
 */
#include <triangMesh.h>

typedef double RType;

int main ( int argc, char **argv ) {
  try {
    if ( argc != 2 ) {
      cerr << "usage: " << argv[0] << " sourceFile.srf" << endl;
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const string outputFileName = aol::getBaseFileName ( inputFileName ) + ".ply";
    aol::TriangMesh<RType> triangMesh;
    triangMesh.loadFromSRF ( inputFileName );
    triangMesh.saveAsPLY ( outputFileName.c_str() );
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
