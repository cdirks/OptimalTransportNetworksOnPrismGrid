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
 * \brief Rotates a 3D quoc array around the center (angle needs to be specified in radians), shifts the array and saves the result with double precision.
 *
 * Usage: rotateAndTranslateQuoc3d.cpp  input angle1 angle2 angle3 xshift yshift zshift
 *
 * \author Berkels
 */

#include <configurators.h>
#include <generator.h>
#include <paramReg.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_3D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3>, qc::CellCenteredCubicGrid<DimensionChoice> > ConfType;
//typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main ( int argc, char **argv ) {
  try {
    if ( argc != 8 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file> <angle1> <angle2> <angle3> <xshift> <yshift> <zshift>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform(); 
      return EXIT_FAILURE;
    }
    const string inputFileName = argv[1];

    typedef qc::ParametricRigidBodyMotion3D<ConfType> ParametricDefType;
    aol::MultiVector<RType> deformParameters ( ParametricDefType::getDeformParametersSize() );
    for ( int i = 0; i < 3; ++i ) {
      deformParameters[0][i] = atof ( argv[5+i] );
      deformParameters[1][i] = atof ( argv[2+i] );
    }

    const string outFileName = aol::strprintf ( "%s_rottrans.dat.bz2", aol::getBaseFileName ( inputFileName ).c_str() );

    qc::ScalarArray<RType, qc::QC_3D> inputArray ( inputFileName.c_str() );
    qc::ScalarArray<RType, qc::QC_3D> rotatedArray ( inputArray, aol::STRUCT_COPY );

    ConfType::InitType grid ( qc::GridSize<ConfType::Dim>::createFrom ( inputArray ) );
    qc::DataGenerator<ConfType> generator ( grid );
    qc::MultiArray<RType, ConfType::Dim> deformation ( grid );
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
