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
 *  \brief Computes the integral of a 3D-function over the $0$-level-set of another 3D-function.
 *
 *  Computes the integral of a 3D-function over the $0$-level-set of another 3D-function.
 *  Both functions have to be given as a ScalarArray<QC\_3D>.
 *  Usage: integrateFunctionOverLevelSet LevelSet3D Function3D
 *
 *  \author Oliver Nemitz
 */


#include <aol.h>
#include <eikonalNA.h>
#include <quoc.h>
#include <scalarArray.h>
#include <integrateLevelSet.h>

typedef double RealType;

int main ( int argc, char **argv ) {
  try {
    if ( argc < 3 ) {
      cerr << "Computes the integral of a function over the 0 level-set in 3D.\n";
      cerr << "usage: " << argv[0] << " levelSet function" << endl;
      return EXIT_FAILURE;
    }
    qc::ScalarArray<RealType, qc::QC_3D> levelSet ( argv[1] );
    qc::ScalarArray<RealType, qc::QC_3D> function ( argv[2] );

    int N = levelSet.getNumX ();
    int d = qc::logBaseTwo ( N );

    qc::GridDefinition grid ( d, qc::QC_3D );

    qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> integrator ( grid, levelSet, function );
    cerr << "Integral: " << integrator.integrate() << endl;
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
