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

#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>

#include <bitVector.h>
#include <boomerAMG.h>
#include <shapeLevelsetGenerator.h>

#include "levelsets.h"



#ifdef USE_EXTERNAL_HYPRE

#include <boomerAMG.h>

/** Very crude version of solver for CFE_CD problems using HYPRE's boomerAMG that "compresses" the matrix by removing identity rows and columns corresponding to non-DOFs
 *  \author Schwen
 */

template< typename GridType, typename OpType, typename RealType >
class CFECDAMGSolver : public aol::InverseOp< aol::Vector<RealType> > {
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

  std::vector<int>       _AMG2rNodeIndices;

public:
  CFECDAMGSolver ( const GridType &Grid, const OpType &Op, const RealType eps = 1.0e-16, const int maxIter = 100 ) : _epsilon ( eps ), _maxIter ( maxIter ), _numLastIter ( -1 ) {
    int rNodes = Op.getNumRows(), amgNodes = Grid.getNumDOFs();

    pMatUnion.hyprematp = &pMatObj;

    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, amgNodes-1, 0, amgNodes-1, &hMat); // must be one less WHY???
    HYPRE_IJMatrixSetObjectType( hMat, HYPRE_PARCSR );
    HYPRE_IJMatrixInitialize( hMat );

    int rNodeCounter = 0, amgCounter = 0;
    std::map<int,int>      r2AMGNodeIndices;

    for ( rNodeCounter = 0; rNodeCounter < rNodes; ++rNodeCounter ) {
      if ( Grid.isDOFNode ( rNodeCounter ) ) {
        _AMG2rNodeIndices.push_back ( rNodeCounter );
        r2AMGNodeIndices[ rNodeCounter ] = amgCounter;
        ++amgCounter;
      }

    }

    for ( rNodeCounter = 0; rNodeCounter < rNodes; ++rNodeCounter ) {
      if ( Grid.isDOFNode ( rNodeCounter ) ) {
        std::vector< aol::Row<double>::RowEntry > revec;
        Op.makeNonZeroRowEntries ( revec, rNodeCounter );

        int nnz = revec.size();
        int* cols = new int[nnz];
        double* values = new double[nnz];

        for ( int j = 0; j < nnz; ++j ) {
          cols[j] = r2AMGNodeIndices [ revec[j].col ];
          values[j] = revec[j].value;
        }

        HYPRE_IJMatrixSetValues( hMat, 1, &nnz, &r2AMGNodeIndices [ rNodeCounter ], cols, values);

        delete[] cols;
        delete[] values;
      } else {
        std::vector< aol::Row<double>::RowEntry > revec;
        Op.makeRowEntries ( revec, rNodeCounter );
      }
    }

    HYPRE_IJMatrixAssemble(hMat);
    HYPRE_IJMatrixGetObject( hMat, pMatUnion.voidpp );

    /* Create solver */
    HYPRE_BoomerAMGCreate(&solver);

    /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_BoomerAMGSetPrintLevel(solver, 3);            /* print solve info + parameters */
    HYPRE_BoomerAMGSetCoarsenType(solver, 6);           /* Falgout coarsening */
    HYPRE_BoomerAMGSetRelaxType(solver, 3);             /* G-S/Jacobi hybrid relaxation */
    HYPRE_BoomerAMGSetNumSweeps(solver, 1);             /* Sweeeps on each level */
    HYPRE_BoomerAMGSetMaxLevels(solver, 25);            /* maximum number of levels */
    HYPRE_BoomerAMGSetTol(solver, sqrt ( _epsilon ) );  /* solver accuracy */
    HYPRE_BoomerAMGSetMaxIter(solver, _maxIter );       /* maximal number of iterations */

  }

  ~CFECDAMGSolver ( ) {
    HYPRE_BoomerAMGDestroy ( solver );
    HYPRE_IJMatrixDestroy ( hMat );
  }

  virtual void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType result ( Dest, aol::DEEP_COPY ); // so that initial guess can be used
    apply ( Arg, result );
    Dest += result;
  }

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {

    int* argIndices = new int[ _AMG2rNodeIndices.size() ];
    int* destIndices = new int[ _AMG2rNodeIndices.size() ];
    double* argValues = new double[ _AMG2rNodeIndices.size() ];
    double* destValues = new double[ _AMG2rNodeIndices.size() ];

    for ( int i = 0; i < static_cast<int>( _AMG2rNodeIndices.size() ); ++i ) {
      argIndices[i] = i;
      destIndices[i] = i;

      argValues[i] = Arg[ _AMG2rNodeIndices[i] ];
      destValues[i] = Dest[ _AMG2rNodeIndices[i] ];
    }

    HYPRE_IJVector hArg;
    HYPRE_IJVector hDest;
    HYPRE_IJVectorCreate( MPI_COMM_WORLD, 0,  _AMG2rNodeIndices.size() -1, &hArg );
    HYPRE_IJVectorSetObjectType( hArg, HYPRE_PARCSR );
    HYPRE_IJVectorInitialize( hArg );
    HYPRE_IJVectorSetValues( hArg,  _AMG2rNodeIndices.size() , argIndices, argValues );
    HYPRE_IJVectorAssemble( hArg );

    HYPRE_IJVectorCreate( MPI_COMM_WORLD, 0,  _AMG2rNodeIndices.size() -1, &hDest );
    HYPRE_IJVectorSetObjectType( hDest, HYPRE_PARCSR );
    HYPRE_IJVectorInitialize( hDest );
    HYPRE_IJVectorSetValues( hDest,  _AMG2rNodeIndices.size() , destIndices, destValues );
    HYPRE_IJVectorAssemble( hDest );

    HYPRE_ParVector pArg, pDest;
    HypreParVectorUnion pArgU, pDestU;
    pArgU.hyprevecp = &pArg;
    pDestU.hyprevecp = &pDest;

    HYPRE_IJVectorGetObject( hArg, pArgU.voidpp );
    HYPRE_IJVectorGetObject( hDest, pDestU.voidpp );


    /* Now setup and solve! */
    HYPRE_BoomerAMGSetup( solver, *(pMatUnion.hyprematp), *(pArgU.hyprevecp), *(pDestU.hyprevecp) ); // since this method apparently needs the vectors, we cannot do it in the constructor.
    HYPRE_BoomerAMGSolve( solver, *(pMatUnion.hyprematp), *(pArgU.hyprevecp), *(pDestU.hyprevecp) );

    /* Run info - needed logging turned on */
    double final_res_norm;

    HYPRE_BoomerAMGGetNumIterations(solver, &_numLastIter);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

    cerr << endl << "Iterations = " << _numLastIter << endl << "Final Relative Residual Norm = " << final_res_norm << endl;

    /* copy back values */
    HYPRE_IJVectorGetValues ( hDest, _AMG2rNodeIndices.size(), destIndices, destValues );
    for ( unsigned int i = 0; i < _AMG2rNodeIndices.size() ; ++i ) {
      Dest[ _AMG2rNodeIndices[i] ] = destValues[i];
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

#endif

typedef double                                            RealType;
typedef tpcfe::CFEGrid< RealType, tpcfe::CFE_CD >         GridType;
typedef tpcfe::CFEBandMatrix<GridType>                    MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
typedef tpcfe::CFEStiffOp  < ConfiguratorType >           StiffOpType;


const int scenario = 0;

// 0: solid cube, top/bottom Dirichlet BCs
// 1: solid cube with 1/7, 0.8 slot

const char* MATRIX_FILENAME = "schwierig_33.csr";
const char* RHS_FILENAME = "schwierig_33_rhs.vec";
const char* SOLN_FILENAME = "schwierig_33_soln.vec";


void setLevelsetNumber ( qc::ScalarArray<RealType, qc::QC_3D> &levelset, const int sc ) {
  switch ( sc ) {
  case 0:
    levelset.setAll ( -1.0 );
    break;

  case 1:
    qc::ShapeLevelsetGenerator<RealType>::generateSlotLevelset( levelset, 1.0/14.0, 0.8 );
    break;

  default:
    cerr << "wrong sc" << endl;
  }
}


void setDirichletBCNumber ( GridType &grid, qc::ScalarArray<RealType, qc::QC_3D> &bcond, qc::BitArray<qc::QC_3D> &DirichletMask, const int sc ) {
  const int width = grid.getNumX();
  switch ( sc ) {
  case 0: {
    for ( int j = 0; j < width; ++j ) {
      for ( int k = 0; k < width; ++k ) {
        if ( grid.isDomainNode ( qc::CoordType ( 0, j, k ) ) ) {
          bcond.set ( 0, j, k, -1.0 );
          DirichletMask.set( 0, j, k, true );
        }

        if ( grid.isDomainNode ( qc::CoordType ( width-1, j, k ) ) ) {
          bcond.set ( width-1, j, k, +1.0 );
          DirichletMask.set( width-1, j, k, true );
        }
      }
    }
  } break;

  case 1: {
    // set boundary conditions and Dirichlet Mask
    for ( int j = 0; j < width; ++j ) {
      for ( int k = 0; k < width; ++k ) {
        if ( grid.isDomainNode ( qc::CoordType ( 0, j, k ) ) ) {
          bcond.set ( 0, j, k, ( 2 * k > width ? 1.0 : -1.0 ) );
          DirichletMask.set( 0, j, k, true );
        }
      }
    }
  } break;

  default:
    cerr << "wrong sc" << endl;
  }
}

template< typename MatrixType, typename DataType >
void writeCSRBinary ( const MatrixType& mat, const char* filename ) {
  unsigned int n;
  std::vector< unsigned int > iA, jA; // we don't know the size in advance, so we need a data structure with variable size
  std::vector< DataType > A;

  n = mat.getNumRows();

  iA.reserve ( n + 1 ); // exact size unless there are zero rows in the matrix
  jA.reserve ( 15 * n ); // upper bound to avoid repeated reallocation
  A.reserve ( 15 * n );

  iA.push_back ( 0 );
  for ( unsigned int row = 0; row < n; ++row ) {
    std::vector< typename aol::Row<DataType>::RowEntry > rowEntries;
    mat.makeRowEntries ( rowEntries, row );
    unsigned short int nnz = 0;

    for ( unsigned int e = 0; e < rowEntries.size(); ++e ) {
      const DataType value = rowEntries[e].value;
      if ( !( value == aol::NumberTrait<DataType>::zero ) ) {
        jA.push_back ( rowEntries[e].col );
        A.push_back ( value );
        ++nnz;
      }
    }
    iA.push_back ( iA[ iA.size() - 1 ] + nnz );
  }

  std::ofstream os( filename, std::ios::out | std::ios::binary );
  if ( !os ) {
    char errmsg[1024];
    sprintf ( errmsg, "writeCSRBinary: could not open %s for output", filename );
    throw aol::Exception ( errmsg, __FILE__, __LINE__ );
  }

#ifdef VERBOSE
   cerr << "n = " << n << endl;
#endif
  os.write ( reinterpret_cast< char* >( &n ), sizeof( unsigned int ) );

#ifdef VERBOSE
  cerr << "iA: is " << iA.size() << ", should be " << n+1 << endl;
#endif
  for ( unsigned int i = 0; i < iA.size(); ++i ) {
#ifdef VERBOSE
    cerr << iA[i] << endl;
#endif
    os.write ( reinterpret_cast< char* >( &iA[i] ), sizeof( unsigned int ) );
  }

#ifdef VERBOSE
  cerr << "jA: is " << jA.size() << ", should be " << iA[n] << endl;
#endif
  for ( unsigned int i = 0; i < jA.size(); ++i ) {
#ifdef VERBOSE
    cerr << jA[i] << endl;
#endif
    os.write ( reinterpret_cast< char* >( &jA[i] ), sizeof( unsigned int ) );
  }

#ifdef VERBOSE
  cerr << "A: is " << A.size() << ", should be " << iA[n] << endl;
#endif
  for ( unsigned int i = 0; i < A.size(); ++i ) {
#ifdef VERBOSE
    cerr << A[i] << endl;
#endif
    os.write ( reinterpret_cast< char* >( &A[i] ), sizeof( DataType ) );
  }
}


//! untested!
template< typename MatrixType, typename DataType >
void readCSRBinary ( MatrixType& mat, const char* filename ){
  mat.setZero();

  std::ifstream is( filename, std::ios::in | std::ios::binary );

  if ( !is ) {
    char errmsg[1024];
    sprintf ( errmsg, "readCSRBinary: could not open %s for input", filename );
    throw aol::Exception ( errmsg, __FILE__, __LINE__ );
  }

  unsigned int n = 0;
  is.read ( reinterpret_cast< char* >( &n ), sizeof ( unsigned int ) );

  unsigned int *iA = new unsigned int [ n+1 ];
  is.read ( reinterpret_cast< char* > ( iA ), (n+1) * sizeof ( unsigned int ) );

  unsigned int *jA = new unsigned int [ iA[n] ];
  is.read ( reinterpret_cast< char* > ( jA ), iA[n] * sizeof ( unsigned int ) );

  DataType* A = new DataType[ iA[n] ];
  is.read ( reinterpret_cast< char* > ( A ), iA[n] * sizeof ( DataType ) );

  for ( unsigned int i = 0; i < n; ++i ) {
    for ( unsigned int j = iA[i]; j < iA[i+1]; ++j ) {
      mat.set ( i, jA[j], A[j] );
    }
  }

  delete[] iA;
  delete[] jA;
  delete[] A;
}

template< typename DataType >
void writeVectorBinary ( const aol::Vector<DataType> &vec, const char* filename ) {

  std::ofstream os ( filename,  std::ios::out | std::ios::binary );

  if ( !os ) {
    char errmsg[1024];
    sprintf ( errmsg, "writeCSRBinary: could not open %s for output", filename );
    throw aol::Exception ( errmsg, __FILE__, __LINE__ );
  }

  const unsigned int n = vec.size();
  os.write ( reinterpret_cast< const char* >( &n ), sizeof( unsigned int ) );

  os.write ( reinterpret_cast< char* >( vec.getData() ), n * sizeof ( DataType) );

}


int main ( int, char** ) {

  const int depth = 7;

  const RealType eps = 1.0e-16;

  try {

    qc::ScalarArray<RealType, qc::QC_2D>::quietMode = 1;

    GridType grid ( depth );

    qc::ScalarArray<RealType, qc::QC_3D>
      levelset ( grid ),
      bcond ( grid ),
      cgsoln ( grid ),
      L_bcond ( grid ),
      rhs ( grid );

    qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
    DirichletMask.setAll ( false );

    setLevelsetNumber ( levelset, scenario );

    // Initialize the grid
    grid.setDomainFrom ( levelset );
    grid.detectAndInitVirtualNodes();


    setDirichletBCNumber ( grid, bcond, DirichletMask, scenario );


    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    // computation
    StiffOpType stiffOp ( grid, aol::ASSEMBLED );
    stiffOp.assembleMatrix();

    restrictNonDomainEntries ( grid, stiffOp.getMatrixRef(), 1.0 );

    grid.restrictToDomain ( bcond );
    stiffOp.apply ( bcond, L_bcond );

    rhs -= L_bcond;

    restrictDirichletEntries ( grid, stiffOp.getMatrixRef(), 1.0 );
    grid.restrictDirichletNodes ( rhs );

    // writeCSRBinary<MatrixType, RealType> ( stiffOp.getMatrixRef(), MATRIX_FILENAME );

    // writeVectorBinary<RealType> ( rhs, RHS_FILENAME );

    aol::StopWatch timer;

#if 0
    {
      timer.start();
      aol::CGInverse< aol::Vector<RealType> > cgsolver ( stiffOp.getMatrixRef(), eps, 10000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      cgsolver.apply ( rhs, cgsoln );
      timer.stop();
      cerr << aol::color::pink << "CG solver took " << timer.elapsedWallClockTime() << " seconds" << aol::color::reset << endl;
      // writeVectorBinary<RealType> ( cgsoln, SOLN_FILENAME );
    }
#endif

#if 1
    {
      aol::Vector<RealType> pcgsoln ( cgsoln, aol::STRUCT_COPY );
      timer.start();
      aol::DiagonalPreconditioner< aol::Vector<RealType> > precond ( stiffOp.getMatrixRef() );
      aol::PCGInverse< aol::Vector<RealType> > pcgsolver ( stiffOp.getMatrixRef(), precond, eps, 10000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      pcgsolver.apply ( rhs, pcgsoln );
      timer.stop();
      cerr << aol::color::pink << "DiagonalPCG solver took " << timer.elapsedWallClockTime() << " seconds" << aol::color::reset << endl;
      //       aol::Vector<RealType> difference ( cgsoln, aol::DEEP_COPY );
      //       difference -= pcgsoln;
      //       cerr << ", norm of difference to CG solution: " << difference.norm() << endl;
      cgsoln = pcgsoln;
    }
#endif

#if 0
    // SSOR works, but it is slow
    {
      aol::Vector<RealType> pcgsoln ( cgsoln, aol::STRUCT_COPY );
      timer.start();
      aol::SSORPreconditioner<aol::Vector<RealType>, MatrixType> precond ( stiffOp.getMatrixRef() );
      aol::PCGInverse< aol::Vector<RealType> > pcgsolver ( stiffOp.getMatrixRef(), precond, eps, 10000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      pcgsolver.apply ( rhs, pcgsoln );
      timer.stop();
      cerr << aol::color::pink << "SSOR(1.2)PCG solver took " << timer.elapsedWallClockTime() << " seconds" << aol::color::reset;
      aol::Vector<RealType> difference ( cgsoln, aol::DEEP_COPY );
      difference -= pcgsoln;
      cerr << ", norm of difference to CG solution: " << difference.norm() << endl;
    }
#endif

#if 0
    // ILU 0 does not work.
    {
      aol::Vector<RealType> pcgsoln ( cgsoln, aol::STRUCT_COPY );
      timer.start();
      aol::ILU0Preconditioner<RealType, MatrixType> precond ( stiffOp.getMatrixRef() );
      aol::PCGInverse< aol::Vector<RealType> > pcgsolver ( stiffOp.getMatrixRef(), precond, eps, 10000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      pcgsolver.apply ( rhs, pcgsoln );
      timer.stop();
      cerr << aol::color::pink << "ILU0PCG solver took " << timer.elapsedWallClockTime() << " seconds" << aol::color::reset;
      aol::Vector<RealType> difference ( cgsoln, aol::DEEP_COPY );
      difference -= pcgsoln;
      cerr << ", norm of difference to CG solution: " << difference.norm() << endl;
    }
#endif

#if 0
    {
       const int min_depth = 2;
       const int NUM_SMOOTH = 2;
       const int VMODE = 1;
       const int mgcycles = 500;

      aol::Vector<RealType> mgsoln ( cgsoln, aol::STRUCT_COPY );
      timer.start();
      tpcfe::CFEMultigrid< StiffOpType, aol::ASSEMBLED > mgsolver ( grid, stiffOp.getMatrixRef(), min_depth, NUM_SMOOTH, NUM_SMOOTH, aol::GAUSS_SEIDEL_FORWARD, 1.0, mg::MULTIGRID_V_CYCLE, eps, mgcycles );
      mgsolver.setVerboseMode ( VMODE );
      mgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      mgsolver.setCoarseSolverSteps( 50000 );
      timer.stop();
      cerr << aol::color::pink << " MG solver setup took " << timer.elapsedWallClockTime() << " seconds" << aol::color::reset << endl;

      timer.start();
      mgsolver.apply ( rhs, mgsoln );
      timer.stop();
      cerr << aol::color::pink << " MG solver solve took " << timer.elapsedWallClockTime() << " seconds" << aol::color::reset;
      aol::Vector<RealType> difference ( cgsoln, aol::DEEP_COPY );
      difference -= mgsoln;
      cerr << ", norm of difference to CG solution: " << difference.norm() << endl;

    }
#endif

#if 1
#ifdef USE_EXTERNAL_HYPRE
    {
      aol::Vector<RealType> amgsoln ( cgsoln, aol::STRUCT_COPY );
      timer.start();
      mg::BoomerAMGSolver< MatrixType, RealType > amgsolver ( stiffOp.getMatrixRef() );
      timer.stop();
      cerr << aol::color::pink << "AMG solver setup took " << timer.elapsedWallClockTime() << " seconds" << aol::color::reset << endl;

      timer.start();
      amgsolver.apply ( rhs, amgsoln );
      timer.stop();
      cerr << aol::color::pink << "AMG solver solve took " << timer.elapsedWallClockTime() << " seconds" << aol::color::reset;
      aol::Vector<RealType> difference ( cgsoln, aol::DEEP_COPY );
      difference -= amgsoln;
      cerr << ", norm of difference to CG solution: " << difference.norm() << endl;
    }
#endif
#endif

#if 1
#ifdef USE_EXTERNAL_HYPRE
    {
      aol::Vector<RealType> amgsoln ( cgsoln, aol::STRUCT_COPY );
      timer.start();
      CFECDAMGSolver< GridType, MatrixType, RealType > amgsolver ( grid, stiffOp.getMatrixRef() );
      timer.stop();
      cerr << aol::color::pink << "CFECDAMG solver setup took " << timer.elapsedWallClockTime() << " seconds" << aol::color::reset << endl;

      timer.start();
      amgsolver.apply ( rhs, amgsoln );
      timer.stop();
      cerr << aol::color::pink << "CFECDAMG solver solve took " << timer.elapsedWallClockTime() << " seconds" << aol::color::reset;
      aol::Vector<RealType> difference ( cgsoln, aol::DEEP_COPY );
      difference -= amgsoln;
      cerr << ", norm of difference to CG solution: " << difference.norm() << endl;
    }
#endif
#endif




  } catch ( aol::Exception &exc ) {
    exc.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}
