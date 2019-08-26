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
 *  \brief using QuocMultigrid
 *
 *   Example for using the QuocMultigrid class, a multigrid solver for
 *   problems on quoc grids with (implicit) homogeneous Neumann
 *   boundary conditions
 *
 *  \author Schwen
 */

#include <quocMultigrid.h>
#include <FEOpInterface.h>
#include <configurators.h>
#include <preconditioner.h>
#include <solver.h>

/** This example illustrates the use of mg::QuocMultigrid derived from mg::ExplicitOperatorHierarchyMultigrid.
 *  That class is a multigrid solver that constructs the coarsened operators (matrices) from the finest one by evaluating them as Restriction o fine Operator o Prolongation
 *  Of course, in this simple example, we could directly specify the coarse operators. If this is possible, it will typically be faster than the operator composition (matrix multiplication) above.
 *  \author Schwen
 */


int main ( int, char** ) {
  try {

    const int depth = 3;

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



    qc::ScalarArray<double, qc::QC_3D>
      img ( grid ),
      rhs ( grid ),
      mg_soln ( grid );


    // we fill our image with some random data to make things interesting ...
    for ( int i = 0; i < img.size(); ++i ) {
      img[i] = static_cast<double> ( rand() ) / static_cast<double> ( RAND_MAX );
    }

    // set up right hand side
    massOp.apply ( img, rhs );

    // now we set up the multigrid solver: just plug in the operator on the finest grid:
    // we coarsen up to grid level 2, do 3 pre- and postsmoothing steps, use a relaxation factor of 1.0 (no overrelaxation), V-cycles and the accuracy parameters defined above.
    mg::QuocMultigrid< double, qc::UniformGridSparseMatrix<double> > mg_solver ( heatMat, grid, 1, 3, 3, 1.0, mg::MULTIGRID_V_CYCLE, eps, steps );
    mg_solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    mg_solver.setVerboseMode ( 2 ); // verbosity mode can be increased up to 6, which produces a lot of output

    // apply solver:
    mg_solver.apply ( rhs, mg_soln );

    // done.



    // We now compare the result to a solution obtained by a standard diagonally preconditioned conjugate gradient solver:

    qc::ScalarArray<double, qc::QC_3D> dpcg_soln ( grid );

    aol::DiagonalPreconditioner< aol::Vector<double> > diag_prec ( heatMat );
    aol::PCGInverse< aol::Vector<double> > dpcg_solver ( heatMat, diag_prec, eps, steps );
    dpcg_solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

    dpcg_solver.apply ( rhs, dpcg_soln );


    // Finally, another application: use exactly one multigrid V-cycle as a preconditioner in a conjugate gradient solver:

    qc::ScalarArray<double, qc::QC_3D> mgpcg_soln ( grid );

    // here, the multigrid is not supposed to solve within one iteration, it is merely used for preconditioning

    mg::QuocMultigrid< double, qc::UniformGridSparseMatrix<double> > mg_prec ( heatMat, grid, 1, 3, 3, 1.0, mg::MULTIGRID_V_CYCLE, 0.99, 1 );
    mg_prec.setVerboseMode ( 2 );
    mg_prec.setStopping( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    aol::PCGInverse< aol::Vector<double> > mgpcg_solver ( heatMat, mg_prec, eps, steps );
    mgpcg_solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

    mgpcg_solver.apply ( rhs, mgpcg_soln );



    // finally, compare the results:
    qc::ScalarArray<double, qc::QC_3D> tmp( grid );

    tmp = mg_soln;
    tmp -= dpcg_soln;

    cerr << "DPCG vs MG:    " << dpcg_soln.norm() << " " << mg_soln.norm() << " " << tmp.norm() <<  endl;


    tmp = mgpcg_soln;
    tmp -= dpcg_soln;

    cerr << "DPCG vs MGPCG: " << dpcg_soln.norm() << " " << mgpcg_soln.norm() << " " << tmp.norm() << endl;

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
