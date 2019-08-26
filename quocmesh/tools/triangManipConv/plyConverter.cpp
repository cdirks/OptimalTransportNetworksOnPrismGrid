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
 *  \brief Conversion from PLY formats to PLY formats and povray
 *
 *  Usage: plyConverter options in.ply out.ply
 *
 *  \author Schwen, von Deylen
 */

#include <triangMesh.h>

int main ( int argc, char ** argv ) {

  // check if two arguments have been given
  if ( argc != 4 ) {
    cerr << "usage: plyConverter option infile outfile" << endl
         << "   options are: -udply2ply (UD-PLY to standard PLY)" << endl
         << "   options are: -ply2udply (standard PLY to UD-PLY)" << endl
         << "   options are: -ply2pov (standard PLY to povray)" << endl
         << "   options are: -udply2pov (UD-PLY to povray)" << endl
         << "   options are: -ply2vtk (standard PLY to legacy VTK)" << endl
         << "   options are: -udply2vtk (UD-PLY to legacy VTK)" << endl;
    return ( EXIT_FAILURE );
  }
  string mode ( argv[1] );

  // check if input file exists
  ifstream ifs ( argv[2], ios::in | ios::binary );
  if ( !ifs ) {
    cerr << "file " << argv[2] << " does not exist. Program will be terminated." << endl;
    return ( EXIT_FAILURE );
  }
  ifs.close();

  // check if output file can be written
  ofstream ofs ( argv[3], ios::out | ios::binary );
  if ( !ofs ) {
    cerr << "file " << argv[3] << " cannot be opened for writing. Program will be terminated." << endl;
    return ( EXIT_FAILURE );
  }
  ofs.close();

  // load
  aol::TriangMesh<long double> mesh;
  if ( ( mode == "-udply2ply" ) ||  ( mode == "-udply2pov" ) ||  ( mode == "-udply2vtk" ) )
    mesh.loadFromUDPLY ( argv[2] );
  else if ( ( mode == "-ply2udply" ) || ( mode == "-ply2pov" ) || ( mode == "-ply2vtk" ) )
    mesh.loadFromPLY ( argv[2] );

  cout << "loaded aol::TriangleMesh<long double> " << argv[2] << " with "
       << mesh.getNumVertices() << " vertices and "
       << mesh.getNumTriangs()  << " triangles." << endl;

  // save
  if ( ( mode == "-udply2pov" ) || ( mode == "-ply2pov" ) )
    mesh.saveAsPov ( argv[3] );
  else if ( mode == "-udply2ply" )
    mesh.saveAsPLY ( argv[3] );
  else if ( mode == "-ply2udply" )
    mesh.saveAsUDPLY ( argv[3] );
  else if ( ( mode == "-ply2vtk" ) || ( mode == "-udply2vtk" ) )
    mesh.saveAsLegacyVTK ( argv[3] );

  cout << "Sucessfully written output file " << argv[3] << "." << endl;

  return 0;
}
