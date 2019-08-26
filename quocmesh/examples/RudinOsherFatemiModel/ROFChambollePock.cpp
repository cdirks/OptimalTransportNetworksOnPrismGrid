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
 * \brief  Denoises an image with the Rudon Osher Fatemi model. The minimization approach used here
 *         is based on a primal / dual algorithm propsed by Antonin Chambolle and Thomas Pock
 *         in {A first-order primal-dual algorithm for convex problems with applications to imaging}.
 *
 * \author Berkels
 */

#include <parameterParser.h>
#include <configurators.h>
#include <timestepSaver.h>
#include <firstOrderTVAlgos.h>

typedef double RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType;

int main( int argc, char **argv ) {

  try {
    aol::ParameterParser parser( argc, argv, "ROFChambollePock.par" );

    aol::StopWatch watch;
    watch.start();

    typedef RType RealType;
    typedef ConfType ConfiguratorType;

    const RealType gamma = parser.getDouble( "gamma" );

    string saveDirectory = parser.getString( "saveDirectory" );
    aol::makeDirectory( saveDirectory.c_str() );
    parser.dumpToFile( "parameter-dump.txt", saveDirectory.c_str() );

    ConfiguratorType::InitType grid ( qc::getSizeFromArrayFile( parser.getString( "input-image" ) ) );
    ConfiguratorType::ArrayType u0 ( grid );
    u0.load( parser.getString( "input-image" ).c_str() );
    u0 /= u0.getMaxValue();

    qc::FirstOrderPrimalDualROFMinimizer<ConfiguratorType> tvAlgo ( grid, gamma, u0, parser.getInt( "numSteps") );
    ConfiguratorType::ArrayType u ( u0 );
    tvAlgo.minimize ( u );
    u.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
    u.savePNG ( aol::strprintf ( "%su.png", saveDirectory.c_str() ).c_str() );

    watch.stop();
    watch.printReport( cerr );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  } 
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
