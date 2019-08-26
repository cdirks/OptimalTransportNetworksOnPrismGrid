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
 *  \brief Computes the signed distance function to a level-set (default is 0) in 2D.
 *
 *    Computes the signed distance function to a level-set (default is 0) in 2D.
 *    Usage: signedDistance2d inputImg outputImg [isovalue]
 *
 *  \author Nemitz.
 */

#include <eikonalNA.h>
#include <quoc.h>
#include <aol.h>
#include <scalarArray.h>
#include <signedDistanceOp.h>
#include <configurators.h>

typedef double RealType;

int main ( int argc, char **argv ) {
  try {
    if ( argc < 3 ) {
      cerr << "Computes the signed distance function to a level-set (default is 0) in 2D.\n";
      cerr << "usage: signedDistance input output [isovalue]\n";
      return EXIT_FAILURE;
    }
    qc::ScalarArray<RealType, qc::QC_2D> inputimg ( argv[1] );
    int N = inputimg.getNumX ();
    int d = qc::logBaseTwo ( N );
    qc::GridDefinition grid ( d, qc::QC_2D );
    qc::ScalarArray<RealType, qc::QC_2D> outputimg ( grid );

    //! \todo Check if this is the appropriate configurator.
    typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfType;
    qc::SignedDistanceOp<ConfType> distOp ( grid );

    if ( argc == 4 ) {
      double isoValue = aol::convert<double> ( argv[3] );
      distOp.setIsoValue ( isoValue );
    } else distOp.setIsoValue ( 0. );

    distOp.apply ( inputimg, outputimg );

    outputimg.save ( argv[2], qc::SaveTypeTrait<RealType>::BinarySaveType );
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
