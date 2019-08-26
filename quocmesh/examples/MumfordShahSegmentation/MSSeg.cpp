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
 * \brief Segments an image into piecewise constant regions with a Mumford Shah model.
 *
 * Segments an image into piecewise constant regions with a Mumford Shah model. A global MS
 * minimizer is obtained by using a quadratic Esedoglu model.
 *
 * \author Berkels
 */

#include <aol.h>
#include <parameterParser.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <ChanVese.h>
#include <generator.h>
#include <chanVeseDescent.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;
typedef aol::IdentityFunction<RType> HeavisideFuncType;

int main( int argc, char **argv ) {

  try {
    aol::ParameterParser parser( argc, argv, "MSSeg.par" );

    aol::StopWatch watch;
    watch.start();

    typedef RType RealType;
    typedef ConfType ConfiguratorType;
    typedef HeavisideFuncType HeavisideFunctionType;

    HeavisideFunctionType heavisideFunction;

    const int numberOfLevelsetFunctions = 1;
    const RealType epsilonMCM = parser.getDouble( "epsilonMCM" );
    const RealType gamma = parser.getDouble( "gamma" );
    const bool usingThresholdModel = true;

    string saveDirectory = parser.getString( "saveDirectory" );
    aol::makeDirectory( saveDirectory.c_str() );
    parser.dumpToFile( "parameter-dump.txt", saveDirectory.c_str() );

    ConfiguratorType::InitType grid ( qc::getSizeFromArrayFile( parser.getString( "input-image" ) ) );
    ConfiguratorType::ArrayType u0 ( grid );
    u0.load( parser.getString( "input-image" ).c_str() );
    u0 /= u0.getMaxValue();
    aol::MultiVector<RealType> u0MVec( 0, 0 );
    u0MVec.appendReference( u0 );

    aol::MultiVector<RealType> meanGrayValues( 1<<numberOfLevelsetFunctions, 1 );
    // Initialize the mean gray values by distributing them equidistant in [0,1].
    for ( int i = 0; i < 1<<numberOfLevelsetFunctions; ++i )
      meanGrayValues[i][0] = i / static_cast<RealType>( (1<<numberOfLevelsetFunctions) - 1 );

    aol::ClassicalChanVeseEnergyMulti<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      tempE( grid, heavisideFunction, u0 );

    aol::HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      Ereg( grid, heavisideFunction, epsilonMCM );
    aol::VariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, numberOfLevelsetFunctions>
      DEreg( grid, epsilonMCM, gamma );

    aol::ChanVeseEnergyLevelsetPart<RealType> Efidelity( tempE, meanGrayValues );

    aol::ClassicalChanVeseEnergyMultiPhiVariation<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      DEfidelity( grid, heavisideFunction, u0, meanGrayValues );

    aol::LinCombOp<aol::MultiVector<RealType>, aol::Scalar<RealType> > E;
    E.appendReference( Efidelity );
    E.appendReference( Ereg, gamma );
    aol::LinCombOp<aol::MultiVector<RealType> > DE;
    DE.appendReference( DEfidelity );
    DE.appendReference( DEreg );

    aol::MultiVector<RealType> levelsetFunctions ( numberOfLevelsetFunctions, grid.getNumberOfNodes() );
    // Since we can't make any reasonable assumption on the segmentation we are looking for,
    // just keep the initial levelsetFunctions at zero everywhere. Instead we initialize the
    // mean gray values (see above) and don't do any post processing in the gradient descent
    // before the itrations starts (and not in the first couple of steps). This seems to work
    // a lot better than trying to initialize the segments as below and to start with post
    // processed gray values.
    //qc::DataGenerator<ConfiguratorType> generator ( grid );
    //generator.generateChanVeseInitialization( levelsetFunctions );

    // This thresholding is only appropriate when converting the zero levelline, so
    // that it can be used in a Esedoglu like model.
    if ( usingThresholdModel )
      levelsetFunctions.threshold ( 0, 0, 1 );

    // ClassicalChanVeseGradientDescent solves for the levelset functions and
    // for the mean gray values (see class descriptions).
    typedef qc::ClassicalChanVeseGradientDescent<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions, 1> GradientDescentType;
    GradientDescentType gradientDescent_solver( grid, E, DE, heavisideFunction, u0MVec, meanGrayValues, parser.getInt( "numSteps") );
    gradientDescent_solver.setConfigurationFlags ( GradientDescentType::USE_NONLINEAR_CG | GradientDescentType::DO_NOT_POSTPROCESS_STARTING_VALUE );
    gradientDescent_solver.setPostProcessingOffset ( 100 );
    gradientDescent_solver.setFlags( GradientDescentType::DO_NOT_SAVE_TIMESTEPS_AS_ISOLINE_PNG );
    if ( usingThresholdModel )
      gradientDescent_solver.configureForQuadraticEsedogluModel();

    gradientDescent_solver.activateSavingAndConfigure( 2, parser.getInt( "saveOffset"), "phi", saveDirectory.c_str() );
    aol::MultiVector<RealType> mtmp( levelsetFunctions, aol::STRUCT_COPY );

    gradientDescent_solver.apply( levelsetFunctions, mtmp );
    levelsetFunctions = mtmp;

    gradientDescent_solver.setSaveName ( "phi_final" );
    for ( int i = 0; i < numberOfLevelsetFunctions; i++ ){
      gradientDescent_solver.saveTimestepBZ2( i, levelsetFunctions[i], grid );
      gradientDescent_solver.savePNG( levelsetFunctions[i], grid, i );
    }
    gradientDescent_solver.savePiecewiseConstantImageSegmented ( meanGrayValues, levelsetFunctions, grid, gradientDescent_solver.getDisplayedIsovalue() );

    watch.stop();
    watch.printReport( cerr );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
