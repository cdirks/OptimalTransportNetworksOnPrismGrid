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
 * \brief Loads the x and y component of a 2D deformation from two files and saves it
 *        into a single file with qc::SaveMultiVector.
 *
 * \warning Only works on dyadic grids.
 *
 * \author Berkels
 */

#include <aol.h>
#include <op.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <configurators.h>
#include <deformations.h>

typedef float RType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType, qc::QC_2D, 3> > ConfType;
int main ( int argc, char **argv ) {

  try {
    char inputFileNameDefX[1024];
    char inputFileNameDefY[1024];

    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <input_file def_x> <input_file def_y>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    if ( argc == 3 ) {
      sprintf ( inputFileNameDefX, "%s",  argv[1] );
      sprintf ( inputFileNameDefY, "%s",  argv[2] );
    }

    qc::GridDefinition grid ( qc::getGridLevelFromArrayFile ( inputFileNameDefX ), ConfType::Dim );

    qc::ScalarArray<RType, qc::QC_2D> defX ( inputFileNameDefX );
    qc::ScalarArray<RType, qc::QC_2D> defY ( inputFileNameDefY );
    aol::MultiVector<RType> def ( 0, 0 );
    def.appendReference ( defX );
    def.appendReference ( defY );

    qc::SaveMultiVector<ConfType> ( grid , def, "output.dat" );

    aol::callSystemPauseIfNecessaryOnPlatform();
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  return 1;
}
