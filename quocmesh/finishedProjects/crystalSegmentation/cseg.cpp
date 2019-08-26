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

#include <aol.h>
#include <parameterParser.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include "cseg.h"

typedef double RType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D, 3> > ConfType;
typedef aol::ArcTanHeavisideFunction<RType> HeavisideFuncType;

int main( int argc, char **argv ) {

  try {
    aol::ParameterParser parser( argc, argv, "ParameterFiles/cseg.par" );

    aol::StopWatch watch;
    watch.start();

    CrystalSegmentation<ConfType, HeavisideFuncType> crystalSegmentation( parser );
    crystalSegmentation.computeTimesteps();

    watch.stop();
    watch.printReport( cerr );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  } 
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
