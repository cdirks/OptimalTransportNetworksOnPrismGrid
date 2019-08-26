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
 * \brief Registers two 3D images using normalized cross-correlation as similarity measure, a parametric
 *        deformation model and a cascadic multiscale minimization algorithm.
 *
 * \author Berkels
 */

#include <paramReg.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_3D;
typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {
  try {

    if ( argc > 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
      return EXIT_FAILURE;
    }

    aol::StopWatch watch;
    watch.start();

    aol::ParameterParser parser( argc, argv, "paramReg.par" );

    typedef qc::ParametricRigidBodyMotion3D<ConfType> ParametricDefType;
    typedef qc::ParametricRegistrationMultilevelDescent<ConfType, ParametricDefType, qc::ParametricNCCEnergy<ConfType, ParametricDefType> > RegisType;
    RegisType mld ( parser, true );

    mld.solve();

    watch.stop();
    cerr << "duration = " << watch.elapsedCpuTimeString() << endl;
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
