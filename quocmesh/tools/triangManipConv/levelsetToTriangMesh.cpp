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
 *  \brief Converts the $0$-levelset of a ScalarArray<QC\_3D> to a triangMesh and saves this as a ply-file.
 *
 * Converts the $0$-levelset of a ScalarArray<QC\_3D> to a triangMesh and saves
 * this as a ply-file.
 * Usage: levelsetToTriangMesh SourceArray3D Dest.ply [shift]
 *
 *  \author Oliver Nemitz
 */

#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>
#include <triangMesh.h>
#include <levelsetToTriangMesh.h>

typedef double RealType;

using namespace aol::color;

int main ( int argc, char **argv ) {
  try {
    if ( argc < 3 ) {
      cerr << "Extracts the 0-levelset as a ply-file.\n";
      cerr << red << "usage: " << argv[0] << " SourceArray3D Dest.ply [shift] [DataArray3D]\n" << reset;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    // load the source-ScalarArray<QC_3D>
    cerr << "loading image '" << green << argv[1] << reset << "'..." << endl;
    qc::ScalarArray<RealType, qc::QC_3D> image ( argv[1] );

    if ( argc >= 4 ) {
      const RealType shift = atof ( argv [3] );
      if ( shift != 0 ) {
        cerr << "subtracting " << shift << " from the input image\n";
        image.addToAll ( -shift );
      }
    }

    // Define the destination-surfmesh
    aol::TriangMesh<RealType> dest;

    // convert
    cerr << "done. " << blue << "Converting ... " << reset;
    qcsm::LevelsetToTriangMesh<RealType> converter;
    converter.apply ( image, dest );

    // save the ply-file
    cerr << "done. " << reset << endl;

    if ( argc >= 5 ) {
      image.load ( argv[4] );
      converter.addVertexData(dest, image);
    }

    dest.saveAsPLY ( argv[2] );

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
