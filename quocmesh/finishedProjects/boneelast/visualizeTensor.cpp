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

#include <levelsetToTriangMesh.h>
#include <shapeLevelsetGenerator.h>
#include <tpCFELevelsets.h>

#include "VoigtTensor.h"

const double scaleFactor = 1.0; // currently .8 for human, .2 for animals

int main ( int argc, char** argv ) {
  if ( argc < 3 ) {
    cerr << "Usage: " << argv[0] << " tensor_file ply_file [shadow ply_file]" << endl;
    return ( EXIT_FAILURE );
  }

  try {

    const int level = 8;

    qc::GridDefinition grid ( level, qc::QC_3D );
    qc::ScalarArray<double, qc::QC_3D> levelset ( grid );
    qc::ShapeLevelsetGenerator<double>::generateBallLevelset ( levelset, 1.0 / 3.0 );

    aol::TriangMesh<double> tmesh;
    qcsm::LevelsetToTriangMesh<double> converter;
    converter.apply ( levelset, tmesh );

    // shift ball to origin and rescale
    for ( int pt = 0; pt < tmesh.getNumVertices(); ++pt ) {
      aol::Vec3<double> coord = tmesh.getVertex ( pt );
      coord += aol::Vec3<double> ( -0.5, -0.5, -0.5 );
      coord *= 3.0;
      tmesh.setVertex ( pt, coord );
    }

    aol::Vector<double> sigmaValues ( tmesh.getNumVertices() ), traceValues ( tmesh.getNumVertices() );

    VoigtElastOp<double> inTensor( argv[1] );
    aol::Matrix33<double> Tensor[3][3];

    inTensor.get4thOrderTensor ( Tensor );

    for ( int pt = 0; pt < tmesh.getNumVertices(); ++pt ) {
      aol::Vec3<double> coord = tmesh.getVertex ( pt );

      coord /= coord.norm(); // make sure coord really has unit length and not just approximately so.

      aol::Matrix33<double> N, S;
      N.makeTensorProduct ( coord, coord );

      for ( short i = 0; i < 3; ++i ) {
        for ( short j = 0; j < 3; ++j ) {
          for ( short k = 0; k < 3; ++k ) {
            for ( short l = 0; l < 3; ++l ) {
              S[i][j] += Tensor[i][j][k][l] * N[k][l];
            }
          }
        }
      }

      sigmaValues[pt] = N.ddprod ( S );

      traceValues[pt] = S.tr();

    }

    cerr << "Sigma range: " << sigmaValues.getMinValue() << " " << sigmaValues.getMaxValue() << endl;

    cerr << "Trace range: " << traceValues.getMinValue() << " " << traceValues.getMaxValue() << endl;

    tmesh.createVertexData();

    for ( int pt = 0; pt < tmesh.getNumVertices(); ++pt ) {
      aol::Vec3<double> coord = tmesh.getVertex ( pt );

      coord *= sigmaValues[pt];

      tmesh.setVertex ( pt, coord );
      tmesh.getVertexData()[pt] = traceValues[pt];
    }

    for ( int pt = 0; pt < tmesh.getNumVertices(); ++pt ) {
      aol::Vec3<double> coord = tmesh.getVertex ( pt );
      coord *= scaleFactor;
      tmesh.setVertex ( pt, coord );
    }


    tmesh.saveAsUDPLY ( argv[2] );

    if ( argc == 4 ) {
      for ( int pt = 0; pt < tmesh.getNumVertices(); ++pt ) {
        aol::Vec3<double> coord = tmesh.getVertex ( pt );
        coord[2] = -1.5;
        tmesh.setVertex ( pt, coord );
        tmesh.getVertexData()[pt] = 0;
      }

      tmesh.saveAsUDPLY ( argv[3] );
    }

    return ( EXIT_SUCCESS );

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return  ( EXIT_FAILURE );
  }
}
