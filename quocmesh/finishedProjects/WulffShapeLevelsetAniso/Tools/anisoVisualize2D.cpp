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

/** *************************************************************************
 * Program for visualizing the Wulff shape or the Frank diagram of a given
 * anisotropy out of anisotropies.h in 2d.
 * \author Nemitz.
 * (last date: 13.10.2006)
 * ************************************************************************* */

#include <configurators.h>
#include <anisotropies.h>
#include <anisotropyVisualization.h>

typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D,3> > ConfigType;

using namespace aol::color;


// now the main program
int main( void ) {

  try {

    // define an anistropy of your choice
//     qcEllipseAnisotropy<double> aniso( 4., 1., 0.001 );
    qc::Hexagon2d<double> aniso( 0.001, 0.7853981632500000000440473768570370793896 );

    // initialize the Visualizer with image size 32x35x43
    qc::AnisotropyVisualizer2d<double, qc::Hexagon2d<double> > visualizer( aniso, 32, 35 );
    // Resize the image (just because it is possible, remember, this is an example...;-)
    visualizer.setSizeOfOutputImg( 8 );     // level 6 => image size = 65^3

    // now generate the Wulff shape (derivatives of anisotropy have to implemented)
    visualizer.generateWulffShapeBySupFormula2d( 3. );
//     visualizer.generateWulffShapeByGammaZ2d( 0.01, 0.01 );
    // save the Wulff shape
    visualizer.saveImage( "WulffShape2d.bz2", 9 );

    // generate the Frank diagram
    visualizer.setVerbose( false );         // no output
    visualizer.generateFrankDiagramByLevelLines2d();
    visualizer.saveImage( "FrankDiagram2d.bz2", 9 );

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}

