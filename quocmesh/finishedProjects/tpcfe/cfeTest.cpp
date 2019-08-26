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

// TODO: this is a self test of the tpcfes that should be moved to quocmesh/selfTest/tpcfe at some point when considered sufficiently stable

#include <tpCFEStandardOp.h>
#include <tpCFETester.h>
#include "affine.h"

static const int DEPTH   = 3;
static const int WIDTH   = ((1<<DEPTH)+1);
static const int WIDTHm2 = (WIDTH-2);

static const double DIFF_COEFF_PLUS  = 10.0;
static const double DIFF_COEFF_MINUS = 10.0;
static const double TANG_A           =  0.0;
static const double TANG_B           =  0.0;

static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS; // we do not want to test CFE_LIEHR by default ...
typedef tpcfe::CFEGrid< double, AT > GridType;
typedef tpcfe::CFEConfigurator < GridType, aol::SparseMatrix<double> > ConfiguratorType;
typedef ConfiguratorType::MatrixType MatrixType;
typedef tpcfe::CFEMassOp < ConfiguratorType > MassOpType;
typedef tpcfe::CFEStiffOp < ConfiguratorType > StiffOpType;
typedef tpcfe::CFEMassOpWI < ConfiguratorType > WMassOpType;
typedef tpcfe::CFEStiffOpWI < ConfiguratorType > WStiffOpType;

int main ( void ) {
  try {
    GridType grid ( DEPTH );
    qc::ScalarArray<double, qc::QC_3D> myInterface ( grid ), affine ( grid );
    qc::AArray<double, qc::QC_3D> diffCoeff ( grid );
    tpcfe::AffineFunction<double> affineFct;

    affineFct.setCoeffs ( DIFF_COEFF_PLUS, DIFF_COEFF_MINUS, TANG_A, TANG_B );
    affineFct.setType( tpcfe::AffineFunction<double>::SPHERICAL_INTERFACE ); affineFct.setRadius(1.0/3.0);
    affineFct.setType ( tpcfe::AffineFunction<double>::ALIGNED_PLANAR_INTERFACE_X );
    affineFct.setWidth ( WIDTH );

    std::cerr << affineFct << std::endl;

    affineFct.createInterface ( &myInterface, &affine, &diffCoeff );

#ifdef LOCAL_ASSEMBLY_TESTER
    tpcfe::CFETester::CFELocalAssemblyTesterLIEHR<MassOpType>();
    tpcfe::CFETester::CFELocalAssemblyTesterLIEHR<StiffOpType>();
#endif

    tpcfe::CFETester::CFEStandardOpTester<MassOpType, DEPTH, AT> ( myInterface, diffCoeff, 1e-15 );
    tpcfe::CFETester::CFEStandardOpTester<StiffOpType, DEPTH, AT> ( myInterface, diffCoeff, 1e-15 );

    tpcfe::CFETester::CFEStandardConvergenceTester<MassOpType, StiffOpType, DEPTH, AT> ( myInterface, diffCoeff, 1e-15 );

  } catch ( aol::Exception &ex ) {
    ex.dump();
    exit ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}
