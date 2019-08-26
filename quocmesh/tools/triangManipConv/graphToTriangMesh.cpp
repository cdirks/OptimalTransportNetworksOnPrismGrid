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
 *  \brief Converts the graph that is given by a 2d-image to a triangMesh and saves this as a ply-file.
 *
 *  Converts the graph that is given by a 2d-image to a triangMesh and saves this as a ply-file.
 *  Usage: graphToTriangMesh SourceArray2d Dest.ply [scaling factor]
 *
 *  \author Oliver Nemitz
 */

#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>
#include <triangMesh.h>
#include <graphToTriangMesh.h>
#include <configurators.h>
#include <progressBar.h>

typedef double RealType;

using namespace aol::color;

int main ( int argc, char **argv ) {
  try {
    if ( argc < 3 ) {
      cerr << "Saves a 2d graph as a ply-file.\n";
      cerr << red << "usage: " << argv[0] << " SourceArray2D Dest.ply [scaling factor]\n" << reset;
      return EXIT_FAILURE;
    }

    // load the source-ScalarArray<QC_3D> and get the grid
    cerr << "loading image '" << green << argv[1] << reset << "'..." << endl;
    qc::ScalarArray<RealType, qc::QC_2D> img ( argv[1] );

    // get a possible scaling factor
    double scaling = 1.;
    if ( argc == 4 ) scaling = atof ( argv[3] );

    // Define the destination-surfmesh
    aol::TriangMesh<RealType> dest;

    // convert
    qcsm::GraphToTriangMesh<RealType> converter;
    converter.setScalingFactor ( scaling );
    converter.apply ( img, dest );

    // save the ply-file
    cerr << reset << "\n\nSaving...  " << reset << endl;
    dest.saveAsUDPLY ( argv[2] );

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
