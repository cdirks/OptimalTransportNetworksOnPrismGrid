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

// include all multigrid header files
#include <boomerAMG.h>
#include <multigrid.h>
#include <quocMultigrid.h>
#include <smoother.h>

#include <FEOpInterface.h>
#include <configurators.h>

int main( int, char** ){

  bool ok = true;

  {
    qc::GridDefinition grid( 4, qc::QC_3D );


    aol::MassOp < qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > >  M_op ( grid, aol::ONTHEFLY );
    qc::UniformGridSparseMatrix<double> M_mat( grid );
    M_op.assembleAddMatrix( M_mat );

    qc::ScalarArray<double, qc::QC_3D> orig( grid ), M_orig( grid ), soln( grid );

    for( int i = 0; i < orig.size(); ++i ){
      orig.set( i, 1.3 * i );
    };

    orig /= orig.norm();

    M_mat.apply( orig, M_orig );

    mg::QuocMultigrid< double, qc::UniformGridSparseMatrix<double> > mg_solver ( M_mat, grid, 2, 3, 3, 1.0, mg::MULTIGRID_V_CYCLE, 1e-16, 100 );
    mg_solver.setStopping( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

    mg_solver.apply ( M_orig, soln );

    soln -= orig;

    ok &= ( soln.norm() < 1.0e-6 );
  }

  if( ok ) {
    aol::printSelfTestSuccessMessage ( "--                  MULTIGRID Self Test Successful                            --" );
    return( EXIT_SUCCESS );
  } else {
    aol::printSelfTestFailureMessage ( "!!                  MULTIGRID Self Test FAILED                                !!" );
    return( EXIT_FAILURE );
  }
}
