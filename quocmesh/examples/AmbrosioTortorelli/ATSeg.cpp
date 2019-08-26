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
 * \brief Joint edge dectection and denoising using the Ambrosio-Tortorelli approximation of the Mumford-Shah model.
 *
 * \author Berkels
 */


#include <AmbrosioTortorelli.h>
#include <quocTimestepSaver.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {

  try {
    aol::ParameterParser parser( argc, argv, "ATSeg.par" );

    aol::StopWatch watch;
    watch.start();

    typedef RType RealType;
    typedef ConfType ConfiguratorType;

    const RealType alpha = parser.getDouble( "alpha" );
    const RealType beta = parser.getDouble( "beta" );
    const RealType epsilon = parser.getDouble( "epsilon" );
    const RealType numSteps = parser.getDouble( "numSteps" );

    qc::DefaultArraySaver<RType, ConfType::Dim> stepSaver;
    stepSaver.initFromParser ( parser, true );

    ConfiguratorType::ArrayType u0 ( parser.getString( "input-image" ) );
    u0 /= u0.getMaxValue();
    ConfiguratorType::InitType grid ( qc::GridSize<ConfType::Dim>::createFrom ( u0 ) );

    ConfiguratorType::ArrayType u ( grid );
    ConfiguratorType::ArrayType v ( grid );

    qc::ATSegmentation<ConfiguratorType> atSegmentation ( epsilon, alpha, beta, u0 );
    atSegmentation.segment ( u, v, numSteps );

    stepSaver.saveStep ( u, -1, "u" );
    stepSaver.saveStep ( v, -1, "v" );

    watch.stop();
    watch.printReport( cerr );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
