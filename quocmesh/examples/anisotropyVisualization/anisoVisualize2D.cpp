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
 * anisotropy defined in anisotropies.h in 2d. A 2d image is created which
 * has to be examined by an external viewer.
 *
 * \author Oliver Nemitz
 */


#include <configurators.h>
#include <anisotropies.h>
#include <anisotropyVisualization.h>

// a typedef for the configurator which contains the RealType, the dimension, the kind of quadrature etc.
typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D,3> > ConfigType;

using namespace aol::color;


int main( void ) {

  try {

    // define an anistropy of your choice (here an ellipse with half-axis 4,1 and a regularization parameter delta=0.001)
    qc::EllipseAnisotropy<double> aniso( 4., 1., 0.001 );

    // initialize the Visualizer with output-image size 32x35
    qc::AnisotropyVisualizer2d<double, qc::EllipseAnisotropy<double> > visualizer( aniso, 32, 35 );
    // Resize the image (just because it is possible, remember, this is an example...;-)
    visualizer.setSizeOfOutputImg( 7 );     // level 7 => image size = 129^3

    // now generate the Wulff shape (derivatives of anisotropy have to implemented)
    // default args = stepsize for getting the derivative maximum and the final traverse-stepsize of the sphere
    visualizer.generateWulffShapeByGammaZ2d( 0.01, 0.01 );
    // finally save the Wulff shape
    visualizer.saveImage( "WulffShape2d.bz2", qc::PGM_DOUBLE_BINARY );

    // generate the Frank diagram
    visualizer.setVerbose( false );         // no output
    visualizer.generateFrankDiagramByLevelLines2d();
    visualizer.saveImage( "FrankDiagram2d.bz2", qc::PGM_DOUBLE_BINARY );

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}

