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

/** \file
 *  \brief Analyses an image (ScalarArrayNd) and prints out the: dimension, min, max, mean, median, standard deviation.
 *
 *  Usage: analyseImg inputImg
 *  Program works both in 2d and 3d.
 *
 *  \author Nemitz
 */

#include <quoc.h>
#include <scalarArray.h>
#include <scalarArray.h>

#include <qmException.h>
#include <auxiliary.h>


using namespace aol::color;

/**
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void printImageInfo ( const char* Filename ) {
  qc::ScalarArray<RealType, Dim> img ( Filename );
  if ( Dim == qc::QC_3D )
    cerr << "3D-Image with dimension (" << img.getNumX() << "," << img.getNumY() << "," << img.getNumZ() << ")\n";
  else if ( Dim == qc::QC_2D )
    cerr << "2D-Image with dimension (" << img.getNumX() << "," << img.getNumY() << ")\n";
  else
    cerr << "Unimplemented image dimension\n";
  cerr << "Minimal Value: " << img.getMinValue() << endl;
  cerr << "Maximal Value: " << img.getMaxValue() << endl;
  cerr << "Mean Value: " << img.getMeanValue() << endl;
  cerr << "Median Value: " << img.getMedianValue() << endl;
  cerr << "Standard Deviation Value: " << img.getStdDev() << endl;
}

typedef double RType;

int main ( int argc, char **argv ) {
  try {
    if ( argc < 2 ) {
      cerr << "Analyses an image (dimension, min, max - feel free to extend the analysis).\n";
      cerr << aol::color::red << "usage: " << argv[0] << " image\n" << aol::color::black;
      return EXIT_FAILURE;
    }

    const qc::Dimension dim = qc::getDimensionFromArrayFile ( argv[1] );

    cerr << blue;
    // ------------------------- --------------------------------------------
    // ------------------------- 3D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( dim == qc::QC_3D )
      printImageInfo<RType, qc::QC_3D> ( argv[1] );

    // ------------------------- --------------------------------------------
    // ------------------------- 2D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( dim == qc::QC_2D )
      printImageInfo<RType, qc::QC_2D> ( argv[1] );

    cerr << aol::color::reset << "done! I hope this information helps :-)\n";
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
