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

#include <configurators.h>
#include <quoc.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <preconditioner.h>

#include <qmException.h>
#include <aol.h>
#include <integrateLevelSet.h>

// for saving the area and the volume in a file
#include "fileOps.h"


#define REAL double

int main( int argc, char **argv)
{
  try
    {
      if ( argc < 4 )
        {
          cerr <<aol::color::red<< "usage: "<<argv[0];
          cerr<<" data AreaOutputFile VolOutputFile\n"<<aol::color::reset;
          return EXIT_FAILURE;
        }

        // handle all args
        for (int count=1; count<argc-2; count++)
        {
          // define the scalar-array
          qc::ScalarArray<REAL, qc::QC_3D> img( argv[count]);
          int N = img.getNumX ();
          int d = qc::logBaseTwo (N);
          qc::GridDefinition grid( d, qc::QC_3D );
          qc::ScalarArray<double, qc::QC_3D> vector1( grid );
          vector1.setAll( 1. );

          // Integrators must rebuild their cfe-structure for each time-step
          // these ops are for integrating over the whole 0-levelset and for getting the included volume
          qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> massIntegrator( grid, img, vector1 );
          qc::LevelSetVolume<double> imgLevelSet( grid, img );

          // calculate the area and the volume
          double LSArea, LSVolume;
          cerr << aol::color::green << "calculating area and volume of the leveset...";
          LSArea   = massIntegrator.integrate();
          LSVolume = imgLevelSet.getVolume();
          cerr << "ready!\nArea: " << LSArea << ", Volume: " << LSVolume << aol::color::black << endl;

          appendAreaAndVolume( argv[argc-2], count-1, LSArea );
          appendAreaAndVolume( argv[argc-1], count-1, LSVolume );
        }

        cerr<<aol::color::blue<<"ready, area-values stored in "<<argv[argc-2];
        cerr<<", volume-values in "<<argv[argc-1]<<"!\n"<<aol::color::reset;

    }
  catch ( aol::Exception &ex )
    {
      ex.dump();
    }
}
