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
 * \brief Registers two images using mutual information as similarity measure (without multiscale minimization) and
 *        serves as template on how to include the registration in custom code.
 *
 * \author Berkels
 */

#include <mutualInformation.h>

typedef float RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {
  try {
    char parameterfilename[1024];

    if ( argc > 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
      return EXIT_FAILURE;
    }
    if ( argc == 2 ) {
      sprintf( parameterfilename, "%s",  argv[1] );
    }
    if ( argc == 1 ) {
      if( ConfType::Dim == 3)
        sprintf( parameterfilename, "ParameterFiles/MIReg_3D.par" );
      else
        sprintf( parameterfilename, "ParameterFiles/MIReg_2D.par" );
    }
    cerr << "Reading parameters from " << parameterfilename << endl;
    aol::ParameterParser parser( parameterfilename );

    aol::StopWatch watch;
    watch.start();

    qc::ScalarArray<RType, ConfType::Dim> referenceImage ( parser.getString ( "reference" ).c_str() );
    qc::ScalarArray<RType, ConfType::Dim> templateImage ( parser.getString ( "template" ).c_str() );
    referenceImage.scaleValuesTo01();
    templateImage.scaleValuesTo01();

    // Initialize the paramters.
    const int numberOfIntensityLevels = parser.getInt ( "numberOfIntensityLevels" );
    const RType beta = static_cast<RType> ( parser.getDouble ( "beta_hist" ) );
    const RType lambda = static_cast<RType> ( parser.getDouble ( "lambda" ) );
    const int numX = referenceImage.getNumX();
    const int numY = referenceImage.getNumY();

    // Set pointers to the C-arrays for both input images.
    RType *_pReferenceData = referenceImage.getData();
    RType *_pTemplateData = templateImage.getData();

    qc::GridSize<qc::QC_2D> gridSize( numX, numY );

    const ConfType::InitType grid ( gridSize );

    qc::MIRegistrationConfigurator<ConfType> regisConf ( numberOfIntensityLevels, beta );
    qc::DirichletRegularizationConfigurator<ConfType> regulConf ( grid, parser );
    // If the third template argument is omitted, the registration doesn't work on non-dyadic grids but is faster.
    qc::StandardRegistration<ConfType, ConfType::ArrayType, qc::MIRegistrationConfigurator<ConfType>, qc::DirichletRegularizationConfigurator<ConfType>, aol::H1GradientDescent<ConfType, aol::MultiVector<ConfType::RealType> > > miRegistration ( gridSize, regisConf, regulConf, lambda );

    // Initialize the template and reference image in the registration class by supplying C-array pointers.
    miRegistration.setTemplImageReference ( _pTemplateData );
    miRegistration.setRefImageReference ( _pReferenceData );

    qc::MultiArray<RType, ConfType::Dim> deformation ( grid );
    // Do the registration. Note: The "true" as second argument supresses the console output.
    miRegistration.findTransformation ( deformation /*, true */ );

    qc::ScalarArray<RType, ConfType::Dim> match ( grid );
    qc::DeformImage<ConfType> ( miRegistration.getTemplImageReference(), grid, match, deformation );
    // Get a pointer to the C-array that contains the deformed template image, i.e. deformed such that it corresponds to the reference image.
    //RType *_pDeformedTemplateData = match.getData();

    // Generate and write a stipe view to check the quality of the registration.
    qc::DataGenerator<ConfType> generator ( grid );
    qc::ScalarArray<RType, ConfType::Dim> checkBoxImage ( grid );
    generator.generateCheckView ( checkBoxImage, miRegistration.getRefImageReference(), match, parser.getInt ( "checkboxWidth" ), true );
    checkBoxImage.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    checkBoxImage.savePNG ( "match.png" );

    watch.stop();
    cerr << "duration = " << watch.elapsedCpuTimeString() << endl;

    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  return 0;
}
