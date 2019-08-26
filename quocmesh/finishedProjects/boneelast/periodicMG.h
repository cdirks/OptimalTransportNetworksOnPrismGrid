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

#ifndef __PERIODICMG_H
#define __PERIODICMG_H

#ifdef USE_EXTERNAL_HYPRE

#include<boomerAMG.h>

namespace mg {

  //! this probably only makes sense for double ...
template< class BlockOpType, typename RealType >
class BoomerAMGBlockSolver : public aol::InverseOp< aol::MultiVector<RealType> > {
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
  typedef aol::MultiVector<RealType> MultiVectorType;

  HYPRE_IJMatrix         hMat;
  HYPRE_ParCSRMatrix     pMatObj;
  HypreParCSRMatrixUnion pMatUnion;
  HYPRE_Solver           solver;

  const RealType         _epsilon;
  const int              _maxIter;
  mutable int            _numLastIter;

public:
  BoomerAMGBlockSolver ( const BlockOpType &Op, const RealType eps = 1.0e-16, const int maxIter = 100 ) : _epsilon ( eps ), _maxIter ( maxIter ), _numLastIter ( -1 ) {
    if ( sizeof ( RealType ) != sizeof ( double ) ) // poor man's type comparison ...
      throw aol::UnimplementedCodeException ( "BoomerAMGBlockSolver only implemented for double", __FILE__, __LINE__ );

    pMatUnion.hyprematp = &pMatObj;

    int rows = Op.getTotalNumRows();

    HYPRE_IJMatrixCreate ( MPI_COMM_WORLD, 0, rows - 1, 0, rows - 1, &hMat ); // must be one less
    HYPRE_IJMatrixSetObjectType ( hMat, HYPRE_PARCSR );
    HYPRE_IJMatrixInitialize ( hMat );

    int i = 0;
    for ( int blockrow = 0; blockrow < Op.getNumRows(); ++blockrow ) {
      for ( int row = 0; row < Op.getReference ( blockrow, 0 ).getNumRows(); ++row, ++i ) {
        std::vector< typename aol::Row<RealType>::RowEntry > revec;
        Op.makeUnblockedRowEntries ( revec, blockrow, row );
        int nnz = revec.size();
        int* cols = new int[nnz];
        double* values = new double[nnz];

        for ( int j = 0; j < nnz; ++j ) {
          cols[j] = revec[j].col;
          values[j] = revec[j].value;
        }

        HYPRE_IJMatrixSetValues ( hMat, 1, &nnz, &i, cols, values );

        delete[] cols;
        delete[] values;
      }
    }

    HYPRE_IJMatrixAssemble ( hMat );
    HYPRE_IJMatrixGetObject ( hMat, pMatUnion.voidpp );

    /* Create solver */
    HYPRE_BoomerAMGCreate ( &solver );

    /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_BoomerAMGSetPrintLevel ( solver, 3 );  /* print solve info + parameters */
    HYPRE_BoomerAMGSetCoarsenType ( solver, 6 ); /* Falgout coarsening */
    HYPRE_BoomerAMGSetRelaxType ( solver, 3 );   /* G-S/Jacobi hybrid relaxation */
    HYPRE_BoomerAMGSetNumSweeps ( solver, 1 );   /* Sweeeps on each level */
    HYPRE_BoomerAMGSetMaxLevels ( solver, 25 );  /* maximum number of levels */
    HYPRE_BoomerAMGSetTol ( solver, _epsilon );      /* conv. tolerance */
    HYPRE_BoomerAMGSetMaxIter ( solver, _maxIter );   /* max. number of iterations */

  }

  ~BoomerAMGBlockSolver ( ) {
    HYPRE_BoomerAMGDestroy ( solver );
    HYPRE_IJMatrixDestroy ( hMat );
  }

  virtual void applyAdd ( const MultiVectorType &Arg, MultiVectorType &Dest ) const {
    MultiVectorType result ( Dest, aol::DEEP_COPY ); // so that initial guess can be used
    apply ( Arg, result );
    Dest += result;
  }

  virtual void apply ( const MultiVectorType &Arg, MultiVectorType &Dest ) const {
    // redesign me! maybe a lot of this can happen in the constructor during setup phase.
    const int n = Arg.getTotalSize();
    int* argIndices = new int[ n ];
    int* destIndices = new int[ n ];
    double* argValues = new double[ n ];
    double* destValues = new double[ n ];

    for ( int bi = 0, i = 0; bi < Arg.numComponents(); ++bi ) {
      for ( int bii = 0; bii < Arg[bi].size(); ++bii, ++i ) {
        argIndices[i] = i;
        destIndices[i] = i;
        argValues[i] = Arg[bi][bii];
        destValues[i] = Dest[bi][bii];
      }
    }

    HYPRE_IJVector hArg;
    HYPRE_IJVector hDest;
    HYPRE_IJVectorCreate ( MPI_COMM_WORLD, 0, n - 1, &hArg );
    HYPRE_IJVectorSetObjectType ( hArg, HYPRE_PARCSR );
    HYPRE_IJVectorInitialize ( hArg );
    HYPRE_IJVectorSetValues ( hArg, n, argIndices, argValues );
    HYPRE_IJVectorAssemble ( hArg );

    HYPRE_IJVectorCreate ( MPI_COMM_WORLD, 0, n - 1, &hDest );
    HYPRE_IJVectorSetObjectType ( hDest, HYPRE_PARCSR );
    HYPRE_IJVectorInitialize ( hDest );
    HYPRE_IJVectorSetValues ( hDest, n, destIndices, destValues );
    HYPRE_IJVectorAssemble ( hDest );

    HYPRE_ParVector pArg, pDest;
    HypreParVectorUnion pArgU, pDestU;
    pArgU.hyprevecp = &pArg;
    pDestU.hyprevecp = &pDest;

    HYPRE_IJVectorGetObject ( hArg, pArgU.voidpp );
    HYPRE_IJVectorGetObject ( hDest, pDestU.voidpp );


    /* Now setup and solve! */
    HYPRE_BoomerAMGSetup ( solver, * ( pMatUnion.hyprematp ), * ( pArgU.hyprevecp ), * ( pDestU.hyprevecp ) ); // since this method apparently needs the vectors, we cannot do it in the constructor.
    HYPRE_BoomerAMGSolve ( solver, * ( pMatUnion.hyprematp ), * ( pArgU.hyprevecp ), * ( pDestU.hyprevecp ) );

    /* Run info - needed logging turned on */
    double final_res_norm;

    HYPRE_BoomerAMGGetNumIterations(solver, &_numLastIter);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm ( solver, &final_res_norm );

    cerr << endl << "Iterations = " << _numLastIter << endl << "Final Relative Residual Norm = " << final_res_norm << endl;

    /* copy back values */
    HYPRE_IJVectorGetValues ( hDest, n, destIndices, destValues );
    for ( int bi = 0, i = 0; bi < Arg.numComponents(); ++bi ) {
      for ( int bii = 0; bii < Arg[bi].size(); ++bii, ++i ) {
        Dest[bi][bii] = destValues[i];
      }
    }

    delete[] argIndices;
    delete[] destIndices;
    delete[] argValues;
    delete[] destValues;
    HYPRE_IJVectorDestroy ( hArg );
    HYPRE_IJVectorDestroy ( hDest );
  }

};

}
#endif


#endif
