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
 *  \brief Reads standard ply files, first scales and then shifts their vertex coordinates
 *         or can apply a full affine transform when used with a parameter file.
 *
 * Usage: affineTransformPLY in.ply out.ply scale\_x scale\_y scale\_z [shift\_x shift\_y shift\_z]
 *
 *  \author Stefan W. von Deylen
 */

#include <parameterParser.h>
#include <triangMesh.h>

using namespace aol;

typedef long double RealType;

int main ( int argc, char ** argv ) {

  if (argc != 2 && argc != 6 && argc != 9) {
    cerr << "usage: " << argv[0] << " <input file> <output file> <scale x y z> [<shift x y z>]" << endl
         << "or     " << argv[0] << " <parameter file>" << endl;
    exit ( 1 );
  }

  try {
    Vec3<RealType> offset;
    Matrix33<RealType> transfMatrix;
    string inFilename;
    string outFilename;

    // *** read arguments ***

    // case 1: scaling (and possibly offset given)
    if (argc >= 6) {
      transfMatrix[0][0] = atof(argv[3]);
      transfMatrix[1][1] = atof(argv[4]);
      transfMatrix[2][2] = atof(argv[5]);

      inFilename = argv[1];
      outFilename = argv[2];
    }

    // case 1a: offset and scaling given
    if (argc == 9) {
      offset[0] = atof(argv[6]);
      offset[1] = atof(argv[7]);
      offset[2] = atof(argv[8]);
    }

    // case 2: parameter file given
    if (argc == 2) {
      ParameterParser parser ( argv[1] );
      inFilename = parser.getString ( "input" );
      outFilename = parser.getString ( "output" );

      Vector<RealType> offsetVector ( 3 );
      parser.getRealVec ( "offset", offsetVector );

      if ( parser.hasVariable ( "transfMatrix" ) ) {
        MultiVector<RealType> transfMultiVec ( 3, 3 );
        parser.getRealMultiVec ( "transfMatrix", transfMultiVec );
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j)
            transfMatrix[i][j] = transfMultiVec[i][j];
          offset[i] = offsetVector[i];
        }
      }
      else if ( parser.hasVariable ( "rotation" ) ) {
        Vector<RealType> rotation ( 3 );
        parser.getRealVec ( "rotation", rotation );
        Matrix33<RealType> rotX, rotY, rotZ;
        rotX.setRotationAboutX ( rotation[0] );
        rotY.setRotationAboutY ( rotation[1] );
        rotZ.setRotationAboutZ ( rotation[2] );
        transfMatrix = rotX;
        transfMatrix *= rotY;
        transfMatrix *= rotZ;
      }
      else {
        transfMatrix[0][0] = 1.;
        transfMatrix[1][1] = 1.;
        transfMatrix[2][2] = 1.;
      }
    }

    // check if input file exists
    ifstream ifs ( inFilename.c_str(), ios::in | ios::binary );
    if (!ifs) {
      cerr << "file " << inFilename << " does not exist. Program will be terminated." << endl;
      exit ( 1 );
    }
    ifs.close();

    // check if output file can be written
    ofstream ofs ( outFilename.c_str(), ios::out | ios::binary );
    if (!ofs) {
      cerr << "file " << outFilename << " cannot be opened for writing. Program will be terminated." << endl;
      exit ( 1 );
    }
    ofs.close();

    // load
    aol::TriangMesh<RealType> mesh;
    mesh.loadBasedOnSuffix ( inFilename );

    clog << "loaded UD PLY " << inFilename << " with "
         << mesh.getNumVertices() << " vertices and "
         << mesh.getNumTriangs()  << " triangles." << endl;

    if ( transfMatrix.isExactlyDiagonal() ) {
      Vec3<RealType> scaling ( transfMatrix[0][0], transfMatrix[1][1], transfMatrix[2][2] );
      mesh.scaleSizeByFactor ( scaling );
    }
    else
      mesh.transformLinear ( transfMatrix );

    mesh.shiftByOffset ( offset );

    // save
    mesh.saveAsPLY ( outFilename );

    clog << "Sucessfully written output file " << outFilename << "." << endl;
  }
  catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return 0;
}
