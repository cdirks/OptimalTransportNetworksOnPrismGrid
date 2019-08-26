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
 * \brief Numerically inverts the deformation (given as displacement) and saves the inverse (as displacement) with double precision.
 *
 * \author Berkels
 */

#include <aol.h>
#include <op.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <multiArray.h>
#include <levelSetDrawer.h>
#include <timestepSaver.h>
#include <configurators.h>
#include <deformations.h>

typedef double RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType, qc::QC_2D, 3> > ConfType;

int main ( int argc, char **argv ) {

  try {
    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << " <input_file def_x> <input_file def_y>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileNameDefX = argv[1];
    const string inputFileNameDefY = argv[2];

    const qc::MultiArray<RType, qc::QC_2D> deformation ( inputFileNameDefX, inputFileNameDefY );
    const ConfType::InitType grid ( deformation[0].getSize() );
    qc::MultiArray<RType, qc::QC_2D> inverse ( grid );
   
    qc::approxInverseDeformation<ConfType> ( grid, deformation, inverse );
    inverse.save ( "inverse_%d.dat.bz2", qc::PGM_DOUBLE_BINARY );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
