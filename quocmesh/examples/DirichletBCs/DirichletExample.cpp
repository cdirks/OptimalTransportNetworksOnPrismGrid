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
 *  \brief using non-zero Dirichlet boundary conditions
 *
 *  This example shows how nonzero Dirichlet boundary conditions can be treated in a FE computation:
 *  \f[     - \Delta u = 0    \mbox{ in } \Omega                   \f]
 *  \f[              u = 0    \mbox{ on } \{ z = 0 \}              \f]
 *  \f[              u = 1+x  \mbox{ on } \{ z = 1 \}              \f]
 *  \f[ \partial_\nu u = 0  \mbox{ elswehere on } \partial \Omega  \f]
 *  where z is the second coordinate in 2D and the third coordinate in 3D.
 *  The system is first transformed to zero Dirichlet boundary conditions (for which we need the full
 *  operator first, but not necessarily assembled), then zero Dirichlet boundary conditions are imposed
 *  by
 *  a) implicitely modifying it, masking out Dirichlet entries
 *  b) assembling the system of equations leaving out rows and columns for Dirichlet entries
 *
 *  \author Schwen
 */

#include <FEOpInterface.h>
#include <configurators.h>
#include <bitVector.h>
#include <solver.h>
#include <quocMatrices.h>
#include <maskedVector.h>
#include <scalarArray.h>


int main ( int /*argc*/, char** /*argv*/ ) {
  try {

    { // first the 2D case:

      const int depth = 5;
      qc::GridDefinition grid ( depth, qc::QC_2D );
      const int N = grid.getWidth();

      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double, qc::QC_2D, 3> > ConfiguratorType;
      aol::StiffOp< ConfiguratorType > stiffOp ( grid, aol::ONTHEFLY );

      qc::ScalarArray<double, qc::QC_2D> boundary_conditions ( grid ), rhs ( grid );

      qc::BitArray<qc::QC_2D> DirichletMask ( qc::GridSize<qc::QC_2D>::createFrom ( grid ) );
      DirichletMask.setAll ( false );

      for ( int i = 0; i < N; ++i ) {
        DirichletMask.set ( i, 0 , true );    // bottom line: zero Dirichlet boundary conditions
        boundary_conditions.set ( i,  0 , aol::NumberTrait<double>::zero );

        DirichletMask.set ( i, N - 1, true );  // top line: nonzero and nonconstant Dirichlet boundary conditions
        boundary_conditions.set ( i, N - 1,  1 + grid.H() * i );   // for boundary condition, need to translate image coordinate to world coordinate
      }

      // set rhs = - L * boundary_conditions: transform system to zero Dirichlet boundary conditions
      // for this purpose, the operator may but need not be assembled
      stiffOp.apply ( boundary_conditions, rhs );
      rhs *= -1.0;

      qc::ScalarArray<double, qc::QC_2D> solution_onthefly ( grid ), solution_assembled ( grid );

      // first we solve the problem by using MaskedVectors and RestrOps
      {
        // all non-Dirichlet nodes are DOF nodes.
        aol::BitVector nonDirichletMask ( DirichletMask.size() );
        for ( int i = 0; i < DirichletMask.size(); ++i ) {
          nonDirichletMask.set ( i, !DirichletMask.get ( i ) );
        }
        aol::DofMask DOFMask ( nonDirichletMask );

        aol::MaskedVector< double > rhs_masked ( rhs.size(), DOFMask ), solution_masked ( solution_onthefly.size(), DOFMask );
        rhs_masked = rhs; // rhs_masked automatically masks out unused entries

        aol::RestrOp<double> stiffOp_masked ( stiffOp, nonDirichletMask );

        aol::CGInverse< aol::Vector<double> > solver_masked ( stiffOp_masked );
        solver_masked.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

        solver_masked.apply ( rhs_masked, solution_masked );

        solution_onthefly = solution_masked;
        solution_onthefly += boundary_conditions;
      }


      // now we solve the problem by assembling the system matrix for DOF entries only
      {
        qc::UniformGridSparseMatrix<double> matrix ( grid );
        qc::ScalarArray<double, qc::QC_2D> rhs_modified ( rhs, aol::DEEP_COPY );

        stiffOp.assembleAddMatrix ( matrix, &DirichletMask, true );
        // what is the last argument "setDirichletNodes = true" for?
        // when you build a system matrix as sum of FE matrices,
        // you will of course want to have diagonal entries "1"
        // for the Dirichlet nodes. Therefore, the "1" there
        // should only be added once, say: In the first summand.
        // So you give "true" for this summand, "false" for every other.



        // This assembly routine does not consider rows and cols corresponding to Dirichlet entries
        // but sets a one on the diagonal. Hence, it is not necessary to first assemble the whole
        // matrix and modify it afterwards.

        // enforce identity equations for zero Dirichlet boundary conditions
        for ( int i = 0; i < DirichletMask.size(); ++i ) {
          if ( DirichletMask[i] == true ) {
            rhs_modified[i] = 0;
          }
        }

        aol::CGInverse< aol::Vector<double> > solver_modified ( matrix );
        solver_modified.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

        solver_modified.apply ( rhs_modified, solution_assembled );

        solution_assembled += boundary_conditions;
      }


      // finally, we compare the two solutions and let the program fail if difference is too big
      {
        qc::ScalarArray<double, qc::QC_2D> difference ( solution_assembled, aol::DEEP_COPY );
        difference -= solution_onthefly;

        if ( difference.norm() > 1.0e-8 )
          return ( EXIT_FAILURE );

        solution_assembled *= 255;
        solution_assembled.save ( "out.pgm", qc::PGM_UNSIGNED_CHAR_ASCII );
      }

    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}
