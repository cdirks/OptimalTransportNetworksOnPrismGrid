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
 * \brief Extracts the zero level set of a levelset function and overlays it over a PNG image.
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

typedef double RType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType, qc::QC_2D, 3> > ConfType;

int main ( int argc, char **argv ) {

  try {
    char inputImageFileName[1024];
    char inputLevelsetFileName[1024];
    char outputFileName[1024];

    if ( argc != 4 && argc != 5 ) {
      cerr << "USAGE: " << argv[0] << " <input_png> <input_levelset> <out_file> [color_channel]" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }
    if ( argc >= 4 ) {
      sprintf ( inputImageFileName, "%s",  argv[1] );
      sprintf ( inputLevelsetFileName, "%s",  argv[2] );
      sprintf ( outputFileName, "%s",  argv[3] );
    }
    int colorChannel = 2;
    if ( argc == 5 )
      colorChannel = atoi ( argv[4] );

    qc::ScalarArray<RType, qc::QC_2D> levelsetFunction ( inputLevelsetFileName );
    const int depth = qc::logBaseTwo ( levelsetFunction.getNumX() );

    qc::GridDefinition grid ( depth, qc::QC_2D );
    qc::MultiArray< RType, 2, 3 > a ( grid );
    a.loadPNG ( inputImageFileName );

    qc::LevelSetDrawer<ConfType> drawer ( grid );
    drawer.draw ( levelsetFunction, a, 0., colorChannel );

    a.setQuietMode ( true );
    a.savePNG ( outputFileName );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
