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

#ifndef __OFF2PLY_H
#define __OFF2PLY_H

#include <triangMesh.h>

int offToPly ( string inFilename, string outFilename, aol::PLY_FORMAT plyFormat ) {
  // check if input file exists
  ifstream ifs ( inFilename.c_str() );
  if ( !ifs ) {
    cerr << "file " << inFilename << " does not exist. Program will be terminated." << endl;
    return EXIT_FAILURE;
  }

  // check if output file can be written
  ofstream ofs ( outFilename.c_str() );
  if ( !ofs ) {
    cerr << "file " << outFilename << " cannot be opened for writing. Program will be terminated." << endl;
    return EXIT_FAILURE;
  }

  // read off header (missing: take care for comment lines)
  cerr << aol::color::blue << "Reading " << inFilename << ", "
          "writing " << outFilename << " ..." << endl << aol::color::reset;
  string line;
  ifs >> line;
  if ( line != "OFF" )
    cerr << "Warning: File should start with \"OFF\", but first line is \"" << line << "\". "
            "Will continue nevertheless." << endl;

  int vertexCount, faceCount, edgeCount;
  ifs >> vertexCount >> faceCount >> edgeCount;

  // get current date/time
  aol::StopWatch timer;
  timer.start();
  timer.stop();
  string dateTime = timer.startedAt();

  // write ply header
  ofs << "ply" << endl
      << "format ascii 1.0" << endl
      << "comment converted from " << inFilename << " "
         "via QuocMesh off2" << (plyFormat == aol::UD_PLY ? "ud" : "" ) << "ply on " << dateTime << endl
      << "element vertex " << vertexCount << endl
      << "property float x" << endl
      << "property float y" << endl
      << "property float z" << endl
      << "element face " << faceCount << endl
      << "property list uchar int vertex_index" << endl
      << "end_header" << endl;

  // copy vertices
  RealType x, y, z;
  for ( int i = 0; i < vertexCount; ++i ) {
    ifs >> x >> y >> z;
    ofs << x << " " << y << " " << z << endl;
  }

  // copy polygons (missing: after each polygon could be RGBA values
  // which are ignored at the moment)
  int n, v;
  for ( int i = 0; i < faceCount; ++i ) {
    line = "";
    while ( line == "" )
      getline ( ifs, line );
    stringstream ls ( line );
    ls >> n;
    if ( plyFormat == aol::STANFORD_PLY )
      ofs << n << " ";
    else if ( n != 3 ) {
      cerr << "Warning: input file contains polygon with more than 3 vertices. "
              "Will skip this polygon." << endl;
      continue;
    }
    for ( int j = 0; j < n; ++j ) {
      ls >> v; ofs << v << " ";
    }
    ofs << endl;
  }

  ifs.close();
  ofs.close();
  cerr << aol::color::blue << "done. " << aol::color::reset << endl;
  return EXIT_SUCCESS;
}

#endif
