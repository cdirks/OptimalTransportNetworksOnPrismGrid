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
 *  \brief using bi-quadratic finite elements
 *
 *  This example demonstrates the use of bi-quadratic finite elements for the Dirichlet problem:
 *  \f[     - \Delta u = 0                              \mbox{ in } [0,1]^2           \f]
 *  \f[              u = \cos(x 2 \pi) \cos(y 2 \pi)    \mbox{ on } \partial [0,1]^2  \f]
 *
 *  A reference solution is computed on a fine grid.
 *  Then bi-linear and bi-quadratic finite elements are used for approximations on successively refined grids.
 *  Finally errors to the reference solution and the experimental order of convergence (EOC) are computed.
 *  Comments in the code shall point out the differences when using higher order finite elements.
 *
 *  \author Geihe
 */

#include <configurators.h>
#include <solver.h>


int main ( int /*argc*/, char** /*argv*/ ) {
  try {
    // here is the most important difference
    typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double, qc::QC_2D, 3> >  ConfiguratorTypeLinear;
    typedef qc::QuocConfiguratorTraitMultiQuad<double, qc::QC_2D, aol::GaussQuadrature<double, qc::QC_2D, 7> > ConfiguratorTypeQuadratic;

    // first compute a reference solution on a maximal fine grid
    cout << endl << "Computing reference solution, please be patient..." << endl;

    // create the finest grid
    const int maxLevel = 8;
    qc::GridDefinition fineGrid ( maxLevel, qc::QC_2D );
    const int fineN = fineGrid.getWidth();
    const double fineH = fineGrid.H();

    qc::ScalarArray<double, qc::QC_2D> fine_solution ( fineGrid ), fine_boundary_values ( fineGrid ), fine_rhs ( fineGrid );

    // fill in nodal values of exact solution
    qc::GridDefinition::OldFullBoundaryNodeIterator iter;
    for ( iter = fineGrid.begin(); iter != fineGrid.end(); ++iter ) {
      fine_boundary_values.set ( *iter, cos ( fineH * ( *iter ) [0] * 2. * aol::NumberTrait<double>::pi ) * cos ( fineH * ( *iter ) [1] * 2. * aol::NumberTrait<double>::pi ) );
    }

    // create stiffness matrix
    aol::StiffOp< ConfiguratorTypeLinear > stiffOp ( fineGrid, aol::ONTHEFLY );
    qc::UniformGridSparseMatrix<double> fine_matrix ( fineGrid );
    stiffOp.assembleAddMatrix ( fine_matrix );

    // specify the Dirichlet boundary (i.e. the whole boundary here)
    qc::BitArray<qc::QC_2D> fine_DirichletMask ( qc::GridSize<qc::QC_2D>::createFrom ( fineGrid ) );
    fine_DirichletMask.setAll ( false );
    for ( int i = 0; i < fineN; ++i ) {
      fine_DirichletMask.set ( i,         0,         true );
      fine_DirichletMask.set ( i,         fineN - 1, true );
      fine_DirichletMask.set ( 0,         i,         true );
      fine_DirichletMask.set ( fineN - 1, i,         true );
    }

    // transform to zero Dirichlet boundary conditions, i.e. move boundary values part to right hand side
    fine_matrix.apply ( fine_boundary_values, fine_rhs );
    fine_rhs *= -1.;
    // as the stiffness matrix will be modified to keep boundary values fixed set them to desired values here
    fine_rhs.assignMaskedFrom ( fine_boundary_values, fine_DirichletMask );

    // solve the system
    fine_DirichletMask.invert();
    aol::RestrOp<double> fine_matrix_masked ( fine_matrix, fine_DirichletMask );
    aol::CGInverse< aol::Vector<double> > fine_solver ( fine_matrix_masked );
    fine_solver.setStopping ( aol::STOPPING_ABSOLUTE );
    fine_solver.setAccuracy ( 1e-28 );
    fine_solver.apply ( fine_rhs, fine_solution );
    fine_DirichletMask.invert();

    // visualize
    fine_solution.setOverflowHandling ( aol::SCALE, 0, 255 );
    fine_solution.save ( "fine_solution.pgm", qc::PGM_UNSIGNED_CHAR_ASCII );



    // linear finite elements
    {
      double loglasterror = 0.;
      double linf_error = 0.;
      double l2_error = 0.;
      cout << endl;
      cout << "Bi-Linear Finite Elements:" << endl;
      cout << "==========================" << endl;
      cout << setw ( 10 ) << left << "Level" << setw ( 20 ) << left << "Linf error" << setw ( 20 ) << left << "L_2 error" << setw ( 20 ) << left << "EOC" << endl;
      cout << "----------------------------------------------------------" << endl;
      for ( int level = 1; level < maxLevel - 2; ++level ) {

        // create the grid
        qc::GridDefinition grid ( level, qc::QC_2D );

        // create stiffness matrix
        aol::StiffOp< ConfiguratorTypeLinear > stiffOp ( grid, aol::ONTHEFLY );
        qc::UniformGridSparseMatrix<double> matrix ( grid );
        stiffOp.assembleAddMatrix ( matrix );

        // create vector for solution etc.
        qc::BitArray<qc::QC_2D> DirichletMask ( qc::GridSize<qc::QC_2D>::createFrom ( grid ) );
        qc::ScalarArray<double, qc::QC_2D> restricted_fine ( grid ), rhs ( grid ), boundary_values ( grid ), solution ( grid );
        int offset = 1 << ( maxLevel - level );
        for ( int x = 0; x < restricted_fine.getNumX(); ++x )  {
          for ( int y = 0; y < restricted_fine.getNumY(); ++y )  {
            boundary_values.set ( x, y, fine_boundary_values.get ( offset*x, offset*y ) );
            restricted_fine.set ( x, y, fine_solution.get ( offset*x, offset*y ) );
            DirichletMask.set   ( x, y, !fine_DirichletMask.get ( offset*x, offset*y ) );
          }
        }

        // transform to zero Dirichlet boundary conditions, i.e. move boundary values part to right hand side
        matrix.apply ( boundary_values, rhs );
        rhs *= -1.;
        // as the stiffness matrix will be modified to keep boundary values fixed set them to desired values here
        rhs.assignMaskedFrom ( boundary_values, DirichletMask, true );

        // solve the system
        aol::RestrOp<double> matrix_masked ( matrix, DirichletMask );
        aol::CGInverse< aol::Vector<double> > solver ( matrix_masked );
        solver.setStopping ( aol::STOPPING_ABSOLUTE );
        solver.setAccuracy ( 1e-28 );
        solver.setQuietMode ( true );
        solver.apply ( rhs, solution );

        // compute difference and linf error
        solution -= restricted_fine;
        linf_error = solution.getMaxAbsValue();
        // compute L_2 error (needs integration)
        aol::MassOp< ConfiguratorTypeLinear > massOp ( grid );
        massOp.apply ( solution, rhs );
        l2_error = solution * rhs;
        cout << setw ( 10 ) << left << level << setw ( 20 ) << left << linf_error << setw ( 20 ) << left << l2_error << setw ( 20 ) << left;
        if ( level > 1 ) cout << loglasterror - log ( l2_error );
        cout << endl;
        loglasterror = log ( l2_error );
      }
    }



    // quadratic finite elements
    {
      double loglasterror = 0.;
      double linf_error = 0.;
      double l2_error = 0.;
      cout << endl;
      cout << "Bi-Quadratic Finite Elements:" << endl;
      cout << "=============================" << endl;
      cout << setw ( 10 ) << left << "Level" << setw ( 20 ) << left << "Linf error" << setw ( 20 ) << left << "L_2 error" << setw ( 20 ) << left << "EOC" << endl;
      cout << "----------------------------------------------------------" << endl;
      for ( int level = 1; level < maxLevel - 2; ++level ) {

        // create the grid
        qc::GridDefinition grid ( level, qc::QC_2D );

        // *****************
        // *** ATTENTION ***  to define quadratic finite elements 9 nodes per element are required!
        // *****************
        int numNodes = grid.getNumberOfNodes();          // this gives the number of real nodes on the grid (i.e. 4 per element)
        numNodes += 0; // avoid warning
        ConfiguratorTypeQuadratic configurator ( grid );
        int numDofs = configurator.getNumGlobalDofs();   // this gives the number of degrees of freedom (i.e. 9 per element)
        int N = sqrt ( numDofs );

        // create stiffness matrix
        aol::StiffOp< ConfiguratorTypeQuadratic > stiffOp ( grid, aol::ONTHEFLY );
        // *****************
        // *** ATTENTION ***  constructor from grids do not work anymore since they assume nodes and degrees of freedom to be the same
        // *** ATTENTION ***  sparsity structure also changes, so specialized matrices do not work anymore
        // *****************
        aol::SparseMatrix<double> matrix ( numDofs, numDofs );
        stiffOp.assembleAddMatrix ( matrix );

        // create vector for solution etc.
        // *****************
        // *** ATTENTION ***  same problem here
        // *****************
        qc::BitArray<qc::QC_2D> DirichletMask ( N, N );
        qc::ScalarArray<double, qc::QC_2D> restricted_fine ( N, N ), rhs ( N, N ), boundary_values ( N, N ), solution ( N, N );

        // *****************
        // *** ATTENTION ***  also here care must be taken
        // *****************
        int offset = 1 << ( maxLevel - level - 1 );
        for ( int x = 0; x < restricted_fine.getNumX(); ++x )  {
          for ( int y = 0; y < restricted_fine.getNumY(); ++y )  {
            boundary_values.set ( x, y, fine_boundary_values.get ( offset*x, offset*y ) );
            restricted_fine.set ( x, y, fine_solution.get ( offset*x, offset*y ) );
            DirichletMask.set   ( x, y, !fine_DirichletMask.get ( offset*x, offset*y ) );
          }
        }

        // transform to zero Dirichlet boundary conditions, i.e. move boundary values part to right hand side
        matrix.apply ( boundary_values, rhs );
        rhs *= -1.;
        // as the stiffness matrix will be modified to keep boundary values fixed set them to desired values here
        rhs.assignMaskedFrom ( boundary_values, DirichletMask, true );

        // solve the system
        aol::RestrOp<double> matrix_masked ( matrix, DirichletMask );
        aol::CGInverse< aol::Vector<double> > solver ( matrix_masked );
        solver.setStopping ( aol::STOPPING_ABSOLUTE );
        solver.setAccuracy ( 1e-26 );
        solver.setQuietMode ( true );
        solver.apply ( rhs, solution );

        // compute difference and linf error
        solution -= restricted_fine;
        linf_error = solution.getMaxAbsValue();
        // compute L_2 error (needs integration)
        aol::MassOp< ConfiguratorTypeQuadratic > massOp ( grid );
        massOp.apply ( solution, rhs );
        l2_error = solution * rhs;
        cout << setw ( 10 ) << left << level << setw ( 20 ) << left << linf_error << setw ( 20 ) << left << l2_error << setw ( 20 ) << left;
        if ( level > 1 ) cout << loglasterror - log ( l2_error );
        cout << endl;
        loglasterror = log ( l2_error );
      }
      cout << endl;
    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}
