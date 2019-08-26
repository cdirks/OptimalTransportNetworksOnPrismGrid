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
 *  \brief Misuse of CFE for affine FE computations
 *
 *  Misuse of CFEs without complicated domain for affine FE computations in 3D.
 *
 *  \author Schwen
 */

#include <tpCFEGrid.h>
#include <tpCFEStandardOp.h>

const tpcfe::ConstraintType CT = tpcfe::CFE_CD; // complicated domain, except we don't set one ...
typedef tpcfe::CFEGrid< double, CT >                                         GridType;

typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEBandMatrix<GridType> >  ConfiguratorType;
typedef ConfiguratorType::MatrixType                                         MatrixType;

typedef tpcfe::CFEMassOp  < ConfiguratorType >                               MassOpType;
typedef tpcfe::CFEStiffOp  < ConfiguratorType >                              StiffOpType;

typedef tpcfe::CFEStiffOpWI  < ConfiguratorType >                            WeightedStiffOpType;



int main ( int, char** ) {
  try {

    const int depth = 4;

    GridType grid ( depth );

    // do not set any domain, no need to detect or init virtual nodes.

    MassOpType massOp ( grid, aol::ONTHEFLY );
    MatrixType dummyMat( grid );
    massOp.assembleAddMatrix( dummyMat );

    StiffOpType stiffOp ( grid, aol::ONTHEFLY);
    stiffOp.assembleAddMatrix( dummyMat );


    // now consider weighted ops, e.g. for a non-constant diffusion coefficient
    qc::AArray<double, qc::QC_3D> weights ( grid );

    for ( qc::RectangularIterator<qc::QC_3D> pit ( grid ); pit.notAtEnd(); ++pit ) {
      weights.getRef ( *pit ) = static_cast<double> ( (*pit)[2] );
    }

    WeightedStiffOpType weightedStiffOp ( weights, grid, aol::ONTHEFLY );
    weightedStiffOp.assembleAddMatrix ( dummyMat );


  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
