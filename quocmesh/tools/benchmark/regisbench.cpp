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
 * \brief Benchmark gradient descent by registering two images using NCC and Dirichlet regularization.
 *
 * Usage: regisbench <level> or regisbench bench file <filename>
 *
 * \author Berkels
 */

#include <registration.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {
  try {
    aol::StopWatch watch;
    watch.start();

    string resultFilename;

    const int bench = aol::checkForBenchmarkArguments ( argc, argv, resultFilename );

    const int gridDepth = ( bench || argc <= 1 ) ? 9 : atoi ( argv[1] );

    ConfType::InitType grid ( gridDepth, ConfType::Dim );

    qc::ScalarArray<RType, ConfType::Dim> referenceImage ( grid );
    qc::ScalarArray<RType, ConfType::Dim> templateImage ( grid );
    qc::DataGenerator<ConfType> generator ( grid );

    // Generate two input images, each one a noisy circle with same size but different center.
    ConfType::VecType center;
    center.setAll ( 0.5 );
    generator.generateSphereLevelset ( referenceImage, center, 0.3 );
    referenceImage.threshold ( 0, 0, 1 );
    center.setAll ( 0.4 );
    generator.generateSphereLevelset ( templateImage, center, 0.3 );
    templateImage.threshold ( 0, 0, 1 );
    aol::NoiseOperator<RType> noiseOp ( aol::NoiseOperator<RType>::GAUSSIAN, 0, 0.1 );
    noiseOp.applyAddSingle ( referenceImage );
    noiseOp.applyAddSingle ( templateImage );

    // Instantiate the registration classes for the registration on a single level.
    // Note that for real world applications the registration should be done using a multi level scheme,
    // but this is just intended to benchmark the gradient descent algorithm.
    qc::NCCRegistrationConfigurator<ConfType> regisConf;
    qc::DirichletRegularizationConfigurator<ConfType> regulConf ( grid );
    qc::StandardRegistration<ConfType, ConfType::ArrayType, qc::NCCRegistrationConfigurator<ConfType>, qc::DirichletRegularizationConfigurator<ConfType>, aol::H1GradientDescent<ConfType, aol::MultiVector<ConfType::RealType>, qc::LinearSmoothOp<ConfType::RealType, qc::MultilevelArrayTrait<ConfType::RealType, ConfType::InitType>::GridTraitType > > > registration ( qc::GridSize<ConfType::Dim> ( grid ), regisConf, regulConf, 5 );
    registration.setTemplImageReference ( templateImage );
    registration.setRefImageReference ( referenceImage );

    // Do the actual registration on a single level.
    qc::MultiArray<RType, ConfType::Dim> deformation ( grid );
    registration.findTransformation ( deformation );

    // Generate and write a stipe view to check the quality of the registration.
    if ( bench == false ) {
      qc::ScalarArray<RType, ConfType::Dim> match ( grid );
      qc::DeformImage<ConfType> ( registration.getTemplImageReference(), grid, match, deformation );
      qc::ScalarArray<RType, ConfType::Dim> checkBoxImage ( grid );
      generator.generateCheckView ( checkBoxImage, registration.getRefImageReference(), match, 16, true );
      checkBoxImage.setOverflowHandlingToCurrentValueRange ( );
      checkBoxImage.savePNG ( "match.png" );
    }

    watch.stop();
    watch.printReport ( cerr );

    // Still needs to be calibrated.
    const RType nupsi = 100 * 131.96 / watch.elapsedCpuTime(), wupsi = 100 * 131.96 / watch.elapsedWallClockTime();
    aol::logBenchmarkResult ( "regisbench", nupsi, wupsi, resultFilename );

    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  return 0;
}
