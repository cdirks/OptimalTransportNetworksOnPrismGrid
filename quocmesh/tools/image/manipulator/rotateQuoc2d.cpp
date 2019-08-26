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
 * \brief Rotates a 2D quoc array around the center (angle needs to be specified in radians) and saves the result with double precision.
 *
 * Usage: rotateQuoc2d input angle
 *
 * \author Berkels
 */

#include <configurators.h>
#include <generator.h>
#include <paramReg.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main ( int argc, char **argv ) {
  try {
    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <angle>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform(); 
      return EXIT_FAILURE;
    }
    const string inputFileName = argv[1];
    const RType angle = atof ( argv[2] );
    const string outFileName = aol::strprintf ( "%s_rotated_%.2f.dat.bz2", aol::getBaseFileName ( inputFileName ).c_str(), angle );

    qc::ScalarArray<RType, qc::QC_2D> inputArray ( inputFileName.c_str() );
    qc::ScalarArray<RType, qc::QC_2D> rotatedArray ( inputArray, aol::STRUCT_COPY );

    ConfType::InitType grid ( qc::GridSize2d::createFrom ( inputArray ) );
    qc::DataGenerator<ConfType> generator ( grid );
    qc::MultiArray<RType, ConfType::Dim> deformation ( grid );
    typedef qc::ParametricRigidBodyMotion2D<ConfType> ParametricDefType;
    aol::MultiVector<RType> deformParameters ( ParametricDefType::getDeformParametersSize() );
    deformParameters[1][0] = angle;
    ParametricDefType parDef ( grid, deformParameters );
    generator.generateDeformationFromParametricDeformation<ParametricDefType, false> ( parDef, deformation );
    qc::DeformImage<ConfType> ( inputArray, grid, rotatedArray, deformation, false );
    rotatedArray.save ( outFileName.c_str(), qc::PGM_DOUBLE_BINARY );

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform(); 
  return 1;
}
