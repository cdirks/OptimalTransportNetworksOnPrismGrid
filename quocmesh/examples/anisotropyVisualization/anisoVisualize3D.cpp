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
 * \brief Example for visualizing the Wulff shape or the Frank diagram of a given anisotropy
 *
 * Example for visualizing the Wulff shape or the Frank diagram of a given
 * anisotropy defined in anisotropies.h in 3d. A 3d image is created which
 * has to be examined by an external viewer.
 *
 * \author Oliver Nemitz
 */


#include <configurators.h>
#include <anisotropies.h>
#include <anisotropyVisualization.h>

// a typedef for the configurator which contains the RealType, the dimension, the kind of quadrature etc.
typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;

using namespace aol::color;


int main( void ) {

  try {

    // define an anistropy of your choice (here an 2d-L1Norm, that is rotated around the z-axis)
    qc::L1Norm2d<double> L1Norm( 0.0000001 );
    typedef qc::Rotated3dAnisotropy<double, qc::L1Norm2d<double> > AnisoType;
    AnisoType aniso( L1Norm, 0.02 );

    // initialize the Visualizer with image size 32x35x43
    qc::AnisotropyVisualizer3d<double, AnisoType > visualizer( aniso, 32, 35, 43 );
    // Resize the image (just because it is possible, remember, this is an example...;-)
    visualizer.setSizeOfOutputImg( 6 );     // level 6 => image size = 65^3

    // now generate the Wulff shape (derivatives of anisotropy have to implemented)
    // default args = stepsize for getting the derivative maximum and the final traverse-stepsize of the sphere
    visualizer.generateWulffShapeByGammaZ3d( 0.01, 0.01 );
    // save the Wulff shape
    visualizer.saveImage( "WulffShape3d.bz2" );

    // generate the Frank diagram
    visualizer.setVerbose( false );         // no output
    visualizer.generateFrankDiagramByLevelLines3d();
    visualizer.saveImage( "FrankDiagram3d.bz2" );

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}

