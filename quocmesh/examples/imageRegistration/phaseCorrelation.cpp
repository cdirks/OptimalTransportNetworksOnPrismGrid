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
 * \brief Registers two images using phase correlation to find the optimal integer shift between the two images.
 *
 * \author Berkels
 */

#include <paramReg.h>
#include <convolution.h>

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

    const ConfType::ArrayType referenceImage ( parser.getString ( "reference" ).c_str() );
    const ConfType::ArrayType templateImage ( parser.getString ( "template" ).c_str() );
    const ConfType::InitType grid ( referenceImage.getSize() );

    qc::CoordType shift = qc::PhaseCorrelationRegistration<ConfType>::registerImages ( referenceImage, templateImage );

    ConfType::ArrayType tmp ( templateImage, aol::STRUCT_COPY );

    tmp.shiftByOffsetFrom ( shift, templateImage );
    tmp.savePNG ( "templShifted.png" );

    watch.stop();
    cerr << "duration = " << watch.elapsedCpuTimeString() << endl;
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
