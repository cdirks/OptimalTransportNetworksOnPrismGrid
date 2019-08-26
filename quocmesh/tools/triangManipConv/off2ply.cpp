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
 *  \brief converts off into ply
 *
 *  Converts an off-file into a ply-file (both file formats to represent
 *  polygonal, not neccessarily triangular, surfaces in 3D).
 *
 *  Usage: off2ply sourceFile.off destinationFile.ply
 *
 *  \author Stefan W. von Deylen
 */
#include <aol.h>

typedef double RealType;

#include "off2ply.h"

int main ( int argc, char **argv ) {
  try {
    if ( argc < 3 ) {
      cerr << "Reads in a .off file and exports it as (standard) .ply file." << endl
           << aol::color::red
           << "usage: " << argv[0] << " source-file.off destination-file.ply" << endl
           << aol::color::reset;
      return EXIT_FAILURE;
    }

  return offToPly ( argv[1], argv[2], aol::STANFORD_PLY );

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
