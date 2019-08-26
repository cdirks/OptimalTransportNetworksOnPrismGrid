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

#ifndef __BOOMERAMG_H
#define __BOOMERAMG_H

#ifdef USE_EXTERNAL_HYPRE

#include <FEOpInterface.h>
#include <configurators.h>
#include <preconditioner.h>
#include <solver.h>

#include <hypreIncludes.h>

namespace mg {

/** Very crude version of solver using HYPRE's boomerAMG, an algebraic multigrid solver
 *  \author Schwen
 */

template< typename OpType, typename RealType >
class BoomerAMGSolver : public aol::InverseOp< aol::Vector<RealType> > {
private:
  //! Union allowing to access same object as void** and HYPRE_ParCSRMatrix without warning about breaking strict aliasing rules
  union HypreParCSRMatrixUnion {
    HYPRE_ParCSRMatrix*  hyprematp;
    void**               voidpp;
  };

  union HypreParVectorUnion {
    HYPRE_ParVector*    hyprevecp;
    void**              voidpp;
  };

protected:
  typedef aol::Vector<double> VectorType;

  HYPRE_IJMatrix         hMat;
  HYPRE_ParCSRMatrix     pMatObj;
  HypreParCSRMatrixUnion pMatUnion;
  HYPRE_Solver           solver;

  const RealType         _epsilon;
  const int              _maxIter;
  mutable int            _numLastIter;

public:
  BoomerAMGSolver ( const OpType &Op, const RealType eps = 1.0e-16, const int maxIter = 100 ) : _epsilon ( eps ), _maxIter ( maxIter ), _numLastIter ( -1 ) {
    pMatUnion.hyprematp = &pMatObj;

    int rows = Op.getNumRows();

    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, rows-1, 0, rows-1, &hMat); // must be one less
    HYPRE_IJMatrixSetObjectType( hMat, HYPRE_PARCSR );
    HYPRE_IJMatrixInitialize( hMat );

    for ( int i = 0; i < rows; ++i ) {
      std::vector< aol::Row<double>::RowEntry > revec;
      Op.makeRowEntries ( revec, i );
      int nnz = revec.size();
      int* cols = new int[nnz];
      double* values = new double[nnz];

      for ( int j = 0; j < nnz; ++j ) {
        cols[j] = revec[j].col;
        values[j] = revec[j].value;
      }

      HYPRE_IJMatrixSetValues( hMat, 1, &nnz, &i, cols, values);

      delete[] cols;
      delete[] values;
    }

    HYPRE_IJMatrixAssemble(hMat);
    HYPRE_IJMatrixGetObject( hMat, pMatUnion.voidpp );

    /* Create solver */
    HYPRE_BoomerAMGCreate(&solver);

    /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
    HYPRE_BoomerAMGSetCoarsenType(solver, 6); /* Falgout coarsening */
    HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* G-S/Jacobi hybrid relaxation */
    HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeeps on each level */
    HYPRE_BoomerAMGSetMaxLevels(solver, 25);  /* maximum number of levels */
    HYPRE_BoomerAMGSetTol(solver, sqrt ( _epsilon ) );      /* conv. tolerance */
    HYPRE_BoomerAMGSetMaxIter(solver, _maxIter );      /* conv. tolerance */

    HYPRE_BoomerAMGSetup( solver, *(pMatUnion.hyprematp), NULL, NULL ); // it seems like the Setup step does not care about the vector pointers being passed here
  }

  ~BoomerAMGSolver ( ) {
    HYPRE_BoomerAMGDestroy ( solver );
    HYPRE_IJMatrixDestroy ( hMat );
  }

  virtual void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType result ( Dest, aol::DEEP_COPY ); // so that initial guess can be used
    apply ( Arg, result );
    Dest += result;
  }

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    // redesign me! maybe a lot of this can happen in the constructor during setup phase.
    const int n = Arg.size();
    int* argIndices = new int[ n ];
    int* destIndices = new int[ n ];
    double* argValues = new double[ n ];
    double* destValues = new double[ n ];

    for ( int i = 0; i < n; ++i ) {
      argIndices[i] = i;
      destIndices[i] = i;
      argValues[i] = Arg[i];
      destValues[i] = Dest[i];
    }

    HYPRE_IJVector hArg;
    HYPRE_IJVector hDest;
    HYPRE_IJVectorCreate( MPI_COMM_WORLD, 0, n-1, &hArg );
    HYPRE_IJVectorSetObjectType( hArg, HYPRE_PARCSR );
    HYPRE_IJVectorInitialize( hArg );
    HYPRE_IJVectorSetValues( hArg, n, argIndices, argValues );
    HYPRE_IJVectorAssemble( hArg );

    HYPRE_IJVectorCreate( MPI_COMM_WORLD, 0, n-1, &hDest );
    HYPRE_IJVectorSetObjectType( hDest, HYPRE_PARCSR );
    HYPRE_IJVectorInitialize( hDest );
    HYPRE_IJVectorSetValues( hDest, n, destIndices, destValues );
    HYPRE_IJVectorAssemble( hDest );

    HYPRE_ParVector pArg, pDest;
    HypreParVectorUnion pArgU, pDestU;
    pArgU.hyprevecp = &pArg;
    pDestU.hyprevecp = &pDest;

    HYPRE_IJVectorGetObject( hArg, pArgU.voidpp );
    HYPRE_IJVectorGetObject( hDest, pDestU.voidpp );


    /* Now setup and solve! */
    // HYPRE_BoomerAMGSetup( solver, *(pMatUnion.hyprematp), *(pArgU.hyprevecp), *(pDestU.hyprevecp) ); // since this method apparently needs the vectors, we cannot do it in the constructor. Or can we, passing NULL pointers?
    HYPRE_BoomerAMGSolve( solver, *(pMatUnion.hyprematp), *(pArgU.hyprevecp), *(pDestU.hyprevecp) );

    /* Run info - needed logging turned on */
    double final_res_norm;

    HYPRE_BoomerAMGGetNumIterations(solver, &_numLastIter);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

    cerr << endl << "Iterations = " << _numLastIter << endl << "Final Relative Residual Norm = " << final_res_norm << endl;

    /* copy back values */
    HYPRE_IJVectorGetValues ( hDest, n, destIndices, destValues );
    for ( int i = 0; i < n; ++i ) {
      Dest[i] = destValues[i];
    }

    delete[] argIndices;
    delete[] destIndices;
    delete[] argValues;
    delete[] destValues;
    HYPRE_IJVectorDestroy(hArg);
    HYPRE_IJVectorDestroy(hDest);
  }

  int getCount ( ) const {
    return ( _numLastIter );
  }

};

}

#endif

#endif
