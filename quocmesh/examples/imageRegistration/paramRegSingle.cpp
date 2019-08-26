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
 *
 * \brief Registers two images using normalized cross-correlation as similarity measure and a parametric
 *        deformation model without multiscale minimization.
 *
 * \author Berkels
 */

#include <paramReg.h>

typedef float RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {
  try {

    if ( argc > 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
      return EXIT_FAILURE;
    }

    aol::StopWatch watch;
    watch.start();

    aol::ParameterParser parser( argc, argv, "paramReg.par" );

    ConfType::ArrayType referenceImage ( parser.getString ( "reference" ).c_str() );
    ConfType::ArrayType templateImage ( parser.getString ( "template" ).c_str() );

    ConfType::InitType grid ( referenceImage.getSize() );
    typedef qc::ParametricRigidBodyMotion2D<ConfType> ParametricDefType;

    aol::MultiVector<RType> unknowns ( ParametricDefType::getDeformParametersSize() );
    ParametricDefType::setIdentityDeformationParameters ( unknowns );
    qc::updateDeformParameters<ConfType, qc::ParametricSSDEnergy<ConfType, ParametricDefType> > ( grid, referenceImage, templateImage, unknowns );

    qc::DataGenerator<ConfType> generator ( grid );
    qc::MultiArray<RType, ConfType::Dim> deformation ( grid );
    ParametricDefType parDef ( grid, unknowns );
    generator.generateDeformationFromParametricDeformation<ParametricDefType, false> ( parDef, deformation );

    ConfType::ArrayType checkImage ( grid );
    qc::DeformImage<ConfType> ( templateImage, grid, checkImage, deformation );
    generator.generateCheckView ( checkImage, referenceImage, checkImage, parser.getInt ( "checkboxWidth" ), true );
    checkImage.savePNG ( "checkView.png" );

    watch.stop();
    cerr << "duration = " << watch.elapsedCpuTimeString() << endl;
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
