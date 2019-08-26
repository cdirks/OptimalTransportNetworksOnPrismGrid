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
 * This tool reads Time-of-Flight (ToF) range data from a *.bin file and
 * generates a 3D mesh. The mesh is scaled to the range [0,1]^3, the mesh
 * center is located at (0.5,0.5,0.5).
 * The tool outputs a scaleFactor and shiftOffset that serve as a basis to
 * transform different range data into the same 3D coordinate system,
 * to do so use tool RangeDataToTriangMesh.
 *
 * author: Bauer
 */

#include <parameterParser.h>
#include <aol.h>
#include <scalarArray.h>
#include <configurators.h>
#include <imageTools.h>
#include <triangMesh.h>

#include <iostream>

#include "TOF.h"


typedef float RType;

int main ( int argc, char **argv ) {
  try {
    // read parameterfile
    char parameterfilename[1024];
    if ( argc == 2 ) {
      sprintf ( parameterfilename, "%s",  argv[1] );
    } else {
      cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
      return EXIT_FAILURE;
    }

    cerr << "Reading parameters from " << parameterfilename << endl;
    aol::ParameterParser parser ( parameterfilename );

    char rangeDataFile[1024];
    char meshDataFile[1024];
    int tofSizeX = 0, tofSizeY = 0;
    int flagSaveRangeDataPNG = 0;
    RType focalLength = 0.;

    if ( parser.hasVariable ( "rangeDataFile" ) ) {
      parser.getString ( "rangeDataFile", rangeDataFile );
    }
    if ( parser.hasVariable ( "meshDataFile" ) ) {
      parser.getString ( "meshDataFile", meshDataFile );
    }
    if ( parser.hasVariable ( "tofSizeX" ) ) {
      tofSizeX = parser.getInt ( "tofSizeX" );
    }
    if ( parser.hasVariable ( "tofSizeY" ) ) {
      tofSizeY = parser.getInt ( "tofSizeY" );
    }
    if ( parser.hasVariable ( "flagSaveRangeDataPNG" ) ) {
      flagSaveRangeDataPNG = parser.getInt ( "flagSaveRangeDataPNG" );
    }
    if ( parser.hasVariable ( "focalLength" ) ) {
      focalLength = parser.getReal<RType> ( "focalLength" );
    }

    // load range data
    std::ifstream RangeDataInputStream;
    RangeDataInputStream.open ( rangeDataFile, std::ios::in | std::ios::binary );
    qc::ScalarArray<RType, qc::QC_2D> rangeData ( tofSizeX, tofSizeY );
    rangeData.loadRaw ( RangeDataInputStream, qc::PGM_FLOAT_BINARY, tofSizeX, tofSizeY );
    RangeDataInputStream.close();

    if ( flagSaveRangeDataPNG == 1 ) {
      rangeData.setOverflowHandling ( aol::CLIP_THEN_SCALE, rangeData.getMinValue(), rangeData.getMaxValue() );
      char rangeDataFilePNG[1024];
      sprintf ( rangeDataFilePNG, "%s.png", rangeDataFile );
      rangeData.savePNG ( rangeDataFilePNG );
    }

    // generate triangle mesh out of range data
    aol::TriangMesh<RType> Mesh;

    // range to world coord transformation
    // for now, perform a generic transformation with scale factor 1.0, shift offset 0.
    RangeToWorldCoordTransformation<RType> transformation ( focalLength, tofSizeX, tofSizeY, 1., 0., 0., 0. );

    for ( int y = 0; y < tofSizeY; y++ ) {
      for ( int x = 0; x < tofSizeX; x++ ) {
        RType r = rangeData.get ( x, y );
        aol::Vec3<RType> g = transformation.evaluateWorldCoord ( x, y, r );
        Mesh.pushBackVertex ( g );
      }
    }

    for ( int y = 0; y < tofSizeY - 1; y++ ) {
      for ( int x = 0; x < tofSizeX - 1; x++ ) {
        aol::Vec3<int> UpperFace ( y * tofSizeX + x + 1, y * tofSizeX + x, ( y + 1 ) * ( tofSizeX ) + x );
        aol::Vec3<int> LowerFace ( ( y + 1 ) * ( tofSizeX ) + x, ( y + 1 ) * ( tofSizeX ) + x + 1, y * tofSizeX + x + 1 );

        Mesh.pushBackTriang ( UpperFace );
        Mesh.pushBackTriang ( LowerFace );
      }
    }

    //std::cout << Mesh.getNumVertices() << std::endl;
    //std::cout << Mesh.getNumTriangs() << std::endl;


    // compute the bounding box
    aol::Vec3<RType> minXYZ, maxXYZ;
    RType maxWidth = Mesh.computeBoundingBox ( minXYZ, maxXYZ );

    std::cout << "maxWidth = " << maxWidth << std::endl;

    // remember information in the z domain to compute the shift offset later
    RType xmin = minXYZ[0];
    RType ymin = minXYZ[1];
    RType zmin = minXYZ[2];
    RType xmaxminusmin = maxXYZ[0] - minXYZ[0];
    RType ymaxminusmin = maxXYZ[1] - minXYZ[1];
    RType zmaxminusmin = maxXYZ[2] - minXYZ[2];

    // shift the center of the mesh to (0,0,0) and scale to the range [0,1]^3
    maxXYZ += minXYZ;
    maxXYZ *= -.5;
    Mesh.shiftByOffset ( maxXYZ );

    RType scaleFactor = 1. / maxWidth;
    minXYZ.setAll ( scaleFactor );
    Mesh.scaleSizeByFactor ( minXYZ );

    // re-shift the mesh center to (0.5,0.5,0.5)
    aol::Vec3<RType> offset;
    offset.setAll ( 0.5 );
    Mesh.shiftByOffset ( offset );

    Mesh.saveAsPLY ( meshDataFile );

    // print scaleFactor and shiftOffset to console
    std::cout << "scaleFactor = " << scaleFactor << std::endl;
    std::cout << "shiftOffsetX = " << 0.5 - ( xmin + xmaxminusmin / 2. ) / maxWidth << std::endl;
    std::cout << "shiftOffsetY = " << 0.5 - ( ymin + ymaxminusmin / 2. ) / maxWidth << std::endl;
    std::cout << "shiftOffsetZ = " << 0.5 - ( zmin + zmaxminusmin / 2. ) / maxWidth << std::endl;
  } catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
}





