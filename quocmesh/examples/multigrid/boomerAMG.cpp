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
 *  \brief using boomerAMG
 *
 *   Example for using the boomerAMG solver, the algebraic multigrid
 *   solver of the HYPRE library.
 *
 *  \author Schwen
 */

#ifdef USE_EXTERNAL_HYPRE

#include <boomerAMG.h>
#include "quocMatrices.h"
#include "configurators.h"
#include "preconditioner.h"

int main ( int, char** ) {

  // skip the introduction ......

  const int depth = 4;

  const int steps = 500;
  const double eps = 1.0e-16;

  qc::GridDefinition grid ( depth, qc::QC_3D );

  aol::StiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > > stiffOp ( grid, aol::ONTHEFLY );
  aol::MassOp<  qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > >  massOp ( grid, aol::ONTHEFLY );

  qc::UniformGridSparseMatrix<double> heatMat ( grid );

  const double tau = 0.5 *  grid.H();

  // set up the matrix we need to invert: M + tau L

  stiffOp.assembleAddMatrix ( heatMat );
  heatMat *= tau;
  massOp.assembleAddMatrix ( heatMat );

  qc::ScalarArray<double, qc::QC_3D> img ( grid ), rhs ( grid );


  // we fill our image with some random data to make things interesting ...
  for ( int i = 0; i < img.size(); ++i ) {
    img[i] = static_cast<double> ( rand() ) / static_cast<double> ( RAND_MAX );
  }

  // set up right hand side
  massOp.apply ( img, rhs );

  // First compute solution via "reliable solver"

  qc::ScalarArray<double, qc::QC_3D> dpcg_soln ( grid );

  aol::DiagonalPreconditioner< aol::Vector<double> > diag_prec ( heatMat );
  aol::PCGInverse< aol::Vector<double> > dpcg_solver ( heatMat, diag_prec, eps, steps );
  dpcg_solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

  dpcg_solver.apply ( rhs, dpcg_soln );

  qc::ScalarArray<double, qc::QC_3D> bamg_soln ( grid );

  // ...... here it becomes interesting: set up and use boomerAMG solver.

  mg::BoomerAMGSolver< qc::UniformGridSparseMatrix<double>, double > bas ( heatMat ); // set up

  bas.apply ( rhs, bamg_soln ); // apply

  qc::ScalarArray<double, qc::QC_3D> difference ( dpcg_soln, aol::DEEP_COPY );
  difference -= bamg_soln;

  cerr << "Comparison: " << difference.norm() << " relative to " << dpcg_soln.norm() << endl;
  if ( difference.norm() / dpcg_soln.norm() < 1e-6 )
    return ( 0 );
  else
    return ( 1 );
}

#else

int main ( int, char** ) {
  return ( 0 );
}

#endif
