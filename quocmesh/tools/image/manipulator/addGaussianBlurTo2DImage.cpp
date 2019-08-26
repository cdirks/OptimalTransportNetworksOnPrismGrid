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
 * \brief Adds Gaussian blur to a 2D image by solving the heat equation with multi linear Finite Elements.
 *
 * \author Berkels
 */

#include <configurators.h>
#include <linearSmoothOp.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main ( int argc, char **argv ) {
  try {
    aol::Vec2<double> v;

    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << " <input_file> <sigma>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const RType sigma = atof ( argv[2] );

    qc::ScalarArray<RType, ConfType::Dim> image ( inputFileName.c_str() );
    ConfType::InitType grid ( qc::GridSize<ConfType::Dim>::createFrom ( image ) );

    const string outputBaseFileName = aol::strprintf ( "%s_blurred_%.2f", inputFileName .c_str(), sigma );
    image.savePNG ( "input.png" );
    qc::GeneralLinearSmoothOp<ConfType> smoother ( grid, sigma * grid.H() );
    smoother.applySingle ( image );
    image.save ( ( outputBaseFileName + ".dat.bz2" ).c_str(), qc::PGM_DOUBLE_BINARY );
    image.savePNG ( ( outputBaseFileName + ".png" ).c_str() );;
    image.scaleValuesTo01();

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
