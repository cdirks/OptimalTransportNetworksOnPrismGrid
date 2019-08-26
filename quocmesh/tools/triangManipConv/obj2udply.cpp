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
 *  \brief converts obj into ply
 *
 *  Converts an obj-file into a ply-file (both file formats to represent meshed surfaces in 3D).
 *
 *  Usage: obj2udply sourceFile.obj destinationFile.ply
 *
 *  \author Wirth
 */

#include <aol.h>
#include <triangMesh.h>

typedef double RealType;

int main ( int argc, char **argv ) {
  try {
    if ( argc < 3 ) {
      cerr << "Reads in a .obj file and exports it as quoc-.ply file.\n";
      cerr << aol::color::red << "usage: " << argv[0] << " source-file destination-file.ply \n" << aol::color::reset;
      return EXIT_FAILURE;
    }

    // read in .obj-file
    cerr << aol::color::blue << "Reading .obj file ... " << aol::color::reset;
    aol::TriangMesh<RealType> triangMesh;
    char line[256];
    ifstream objFile ( argv[1] );
    while ( !objFile.eof() ) {
      // read in a line of the .obj-file
      objFile.getline ( line, 256 );
      // if line describes a vertex...
      if ( !strncmp ( line, "v ", strlen ( "v " ) ) ) {
        aol::Vec3<RealType> vertex;
        // read in the vertex coordinates and add it to the mesh
        strtok ( line, " " );
        for ( int i = 0; i < 3; i++ )
          vertex[i] = strtod ( strtok ( NULL, " " ), NULL );
        triangMesh.pushBackVertex ( vertex );
      }
      // if line describes a triangle...
      if ( !strncmp ( line, "f ", strlen ( "f " ) ) ) {
        aol::Vec3<int> triangle;
        // read in the vertex indices of the triangle and add it to the mesh
        strtok ( line, " " );
        char *entries[256];
        for ( int i = 0; i < 3; i++ )
          entries[i] = strtok ( NULL, " " );
        for ( int i = 0; i < 3; i++ )
          triangle[i] = atoi( strtok( entries[i], "/" ) ) - 1; // note: node indices start from 0 in .ply, but from 1 in .obj
        triangMesh.pushBackTriang( triangle );
      }
    }
    objFile.close();
    cerr << aol::color::blue << "done. " << aol::color::reset << endl;

    // save the ply-file
    cerr << aol::color::blue << "Saving as .ply file ... " << aol::color::reset;
    triangMesh.saveAsUDPLY ( argv[2] );
    cerr << aol::color::blue << "done. " << aol::color::reset << endl;

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
