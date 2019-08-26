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

/**
 * \file
 * \brief how nonzero Dirichlet boundary conditions can be treated in a FE computation
 *
 * This example shows how nonzero Dirichlet boundary conditions can be treated in a FE computation:
 * \f[     - \Delta u = 0    \mbox{ in } \Omega                     \f]
 * \f[              u = 0    \mbox{ on } \{ z = 0 \}                \f]
 * \f[              u = 1+x  \mbox{ on } \{ z = 1 \}                \f]
 * \f[ \partial_\nu u = 0  \mbox{ elswehere on } \partial \Omega    \f]
 * We prescribe the Dirichlet-values by performing a (only virtual) splitting of the coefficient
 * vector $U = EU^{int} + U^{ext}$, where $E$ is an extension of the vector of interior nodal values
 * to all nodes by 0 and $U^{ext}$ has only non-zero values on the boundary nodes. Inserting this into
 * the matrix equation $LU=0$ gives
 * $$ L(EU^{int} + U^{ext}=0 \quad \Rightarrow \quad R L E U^{int} = - R L U^{ext}. $$
 * The operator $R$ restricts to the space of interior nodes again, because this is the space where
 * we want to solve our problem. By applying $L$ on the rhs, the boundary values have an implicit
 * influence to the solution of our problem. \\
 * The matrix $L$ is as usual assembled on the whole domain, no treatment of lines and columns is
 * necessary. For the extension and restriction we use a maskedOp, which calls an applyMasked of the
 * operator. This applyMasked can work on masked regions and thus simulate the extension and restriction.
 *
 * \author Oliver Nemitz
*/

#include <FEOpInterface.h>
#include <configurators.h>
#include <bitVector.h>
#include <solver.h>
#include <quocMatrices.h>
#include <maskedVector.h>
#include <maskedOp.h>
#include <scalarArray.h>

typedef qc::UniformGridSparseMatrix<double>   MatrixType;
typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double, qc::QC_2D, 3> > ConfiguratorType;


int main ( int /*argc*/, char** /*argv*/ ) {
  try {

    {

      const int depth = 5;
      qc::GridDefinition grid ( depth, qc::QC_2D );
      const int N = grid.getWidth();

      // define two vectors, the solution, which will also contain
      // the Dirichlet-values, and the rhs.
      qc::ScalarArray<double, qc::QC_2D> solution( grid ), rhs ( grid );

      // define a mask that defines the boundary nodes and fill the solution vector at these
      // nodes with the desired Dirichlet-values.
      // The mask is 1 for inner nodes (degrees of freedom) and 0 for nodes on the Dirichlet boundary.
      qc::BitArray<qc::QC_2D> DOFMask ( qc::GridSize<qc::QC_2D>::createFrom ( grid ) );
      DOFMask.setAll ( true );

      for ( int i = 0; i < N; ++i ) {
        DOFMask.set ( i, 0 , false );              // bottom line: zero Dirichlet boundary conditions
        solution.set ( i,  0 , aol::NumberTrait<double>::zero );

        DOFMask.set ( i, N - 1, false );           // top line: nonzero and nonconstant Dirichlet boundary conditions
        solution.set ( i, N - 1,  1 + grid.H() * i );    // for boundary condition, need to translate image coordinate to world coordinate
      }

      // Define the maskedOp which will just contain the stiffness-matrix. The default-mode is
      // INCLUDE_INT_WRITE_INT, i.e. it includes only interior nodes in the computation and writes
      // the result only to interior nodes. This is what we need in the cg-solver.
      aol::MaskedOp<MatrixType> maskedStiffMatrix ( grid, DOFMask, aol::INCLUDE_INT_WRITE_INT );

      // The assembly works as usual, the MaskedOp contains a matrix of MatrixType, which is assembled now.
      aol::StiffOp<ConfiguratorType> stiffOp ( grid );
      stiffOp.assembleAddMatrix ( maskedStiffMatrix );

      // Now build the rhs. Here we have to include boundary nodes, but write the result only to interior nodes.
      // Instead of calling applyMasked, we could set the includeWriteMode and just call the usual apply.
      maskedStiffMatrix.applyMasked( solution, rhs, DOFMask, aol::INCLUDE_BD_WRITE_INT );
      rhs *= -1.;

      // define the solver
      aol::CGInverse< aol::Vector<double> > solver ( maskedStiffMatrix );
      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

      // and solve, this should now only affect the interior nodes.
      solver.apply( rhs, solution );

      // and finally save
      solution *= 255;
      solution.save ( "outSplitting.pgm", qc::PGM_UNSIGNED_CHAR_ASCII );

    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}
