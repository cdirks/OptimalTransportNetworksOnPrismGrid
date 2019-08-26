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
 * \file
 * \brief Rotates a mesh
 *
 * Rotates an aol::TriangMesh around the center of its bounding box (angles need to be specified in degrees)
 * and saves the result as PLY. All mesh types compatible with aol::TriangMesh::loadBasedOnSuffix can be used
 * as input.
 *
 * Usage: rotateTriangMesh inputMesh yaw pitch roll
 *
 * \author Berkels
 */

#include <triangMesh.h>
#include <configurators.h>
#include <paramReg.h>

typedef double RType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_3D, aol::GaussQuadrature<RType,qc::QC_3D,3> > VolConfType;

int main( int argc, char **argv ) {
  try {

    if ( argc != 5 ) {
      cerr << "USAGE: " << argv[0] << " <inputfile> <yaw> <pitch> <roll>" << endl;
      return EXIT_FAILURE;
    }

    const std::string inputFileName = argv[1];
    const aol::Vec3<RType> angles ( atof ( argv[2] ), atof ( argv[3] ), atof ( argv[4] ) );

    aol::TriangMesh<RType> inputMesh;
    inputMesh.loadBasedOnSuffix ( inputFileName.c_str() );

    typedef qc::ParametricRigidBodyMotion3D<VolConfType> ParametricDefType;
    aol::MultiVector<RType> deformParameters ( ParametricDefType::getDeformParametersSize() );
    for ( int i = 0; i < 3; ++i )
      deformParameters[1][i] = aol::DegreesToRadians ( angles[i] );

    // The grid is only needed for the syntax (since transformParametric only calls ParametricDefType::evaluateDeformationOn01.
    VolConfType::InitType volGrid ( 2, qc::QC_3D );
    ParametricDefType parDef ( volGrid, deformParameters );

    aol::Vec3<RType> center = inputMesh.centerOfBoundingBox();
    cerr << center << endl;
    // ParametricRigidBodyMotion3D uses the center of the unit cube as center for the rotation.
    // Cancel this by adding this center to all points.
    inputMesh.shiftByOffset ( aol::Vec3<RType> ( 0.5, 0.5, 0.5 ) );
    // Use the center of the bounding box as center instead.
    inputMesh.shiftByOffset ( -center );
    inputMesh.transformParametric ( parDef );
    inputMesh.shiftByOffset ( center );
    inputMesh.saveAsPLY ( aol::strprintf( "%s_rotated_%.2f_%.2f_%.2f.ply", aol::getBaseFileName ( inputFileName ).c_str(), angles[0], angles[1], angles[2] ).c_str() );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
