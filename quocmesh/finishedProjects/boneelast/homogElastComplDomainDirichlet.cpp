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

#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMultigrid.h>
#include <tpCFEUtils.h>

#include <multiArray.h>

#include "computeForce.h"


typedef double RealType;

typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD > GridType;
typedef tpcfe::CFEHybridMatrix< GridType > MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType > ConfiguratorType;

typedef tpcfe::CFEElastOp< ConfiguratorType >  ElastOpType;

static const aol::GaussSeidelSweepingMode gsmode = aol::GAUSS_SEIDEL_ZEBRA2; // for ||ization

int main ( int, char** argv ) {
  try {
    aol::StopWatch::_suppressOpenMPWarning = true;

    // solver parameters
    const int      solversteps = 10000;
    const RealType epsilon = 1.0e-16;

    // parameters for the compression experiment
    const RealType
      E             = 13.0, // Young's modulus
      nu            = 0.32, // Poisson's ratio
      lambda        = tpcfe::computeLambda ( E, nu ), // Lamé-Navier parameters
      mu            = tpcfe::computeMu ( E, nu );


    // levelset description of interface geometry
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( argv[1] );

    aol::Matrix33<RealType> sigmas[5][3][3], epsilons[5][3][3];

    // different loading cases
    for ( int fix_dir = 0; fix_dir < 3; ++fix_dir ) {
      for ( int shift_dir = 0; shift_dir < 3; ++shift_dir ) {

        qc::GridSize<qc::QC_3D> gSize ( levelset );
        GridType grid ( gSize );
        grid.setDomainFrom ( levelset );
        grid.detectAndInitVirtualNodes();

        // affine Dirichlet boundary values
        qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
        DirichletMask.setAll ( false );

        qc::MultiArray< RealType, qc::QC_3D > DirichletBCs ( grid ), rhs ( grid ), soln ( grid );

        for ( qc::RectangularBoundaryIterator<qc::QC_3D> bbit(grid); bbit.notAtEnd(); ++bbit ) {
          DirichletMask.set ( *bbit, true );
          DirichletBCs[ shift_dir ].set ( *bbit, grid.H() * (*bbit)[ fix_dir ] );
        }

        // import Dirichlet data to grid
        grid.setDirichletMask ( DirichletMask );
        grid.setDOFMaskFromDirichletAndDomainNodeMask();

        // FE elasticity matrix
        tpcfe::CFEElastOp< ConfiguratorType > elastOp ( grid, lambda, mu );

        // transformation to zero Dirichlet BC
        soln -= DirichletBCs;

        elastOp.restrictNonDomainEntries();

        grid.restrictToDomain ( soln );

        elastOp.apply ( soln, rhs );

        elastOp.restrictDirichletEntries();
        grid.restrictDirichletNodes ( rhs );

        // matrix memory requirement (once)
        if ( fix_dir == 0 && shift_dir == 0 ) elastOp.getBlockMatrixRef().getReference(0,0).printStatistics();

        aol::GaussSeidelPreconditioner< aol::MultiVector<RealType>, ElastOpType::BlockMatrixType > prec ( elastOp.getBlockMatrixRef() );
        aol::PCGInverse< aol::MultiVector<RealType> > solver ( elastOp.getBlockMatrixRef(), prec, epsilon, solversteps );

        solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

        aol::StopWatch timer_solve;
        timer_solve.start();
        solver.apply( rhs, soln );
        timer_solve.stop();
        cerr << "solver took " << timer_solve.elapsedWallClockTime() << ", " << timer_solve.elapsedCpuTime() << endl;
        cerr << "memusage " << ( aol::memusage() >> 20 ) << " MiB" << endl;

        soln += DirichletBCs;

        for ( short bndL = 0; bndL < 5; ++bndL ) {
          cerr << "Considering boundary layer of " << bndL << " / 16" << endl;

          const aol::Vec3<int>
            lowerBnd (      bndL     *   grid.getNumX() / 16      ,      bndL     *   grid.getNumY() / 16      ,      bndL     *   grid.getNumZ() / 16       ),
            upperBnd ( ( 16 - bndL ) * ( grid.getNumX() - 1 ) / 16, ( 16 - bndL ) * ( grid.getNumY() - 1 ) / 16, ( 16 - bndL ) * ( grid.getNumZ() - 1 ) / 16 );

          cerr << lowerBnd << " " << upperBnd << endl;

          tpcfe::getSigmaEpsilonViaPartialTetTraversal ( grid, soln, 1.0, lambda, mu, lowerBnd, upperBnd, sigmas[bndL][ fix_dir ][ shift_dir ], epsilons[bndL][ fix_dir ][ shift_dir ] );

          cerr << "Evaluated on elements: " << lowerBnd << " --- " << upperBnd << endl
               << "Sigma = " << endl << sigmas[bndL][ fix_dir ][ shift_dir ] << endl
               << "Epsilon = " << endl << epsilons[bndL][ fix_dir ][ shift_dir ] << endl;
        }
      }

    }


    for ( short bndL = 0; bndL < 5; ++bndL ) {
      cerr << "Considering boundary layer of " << bndL << " / 16" << endl;
      tpcfe::convertAverageAndDumpTensors ( sigmas[bndL], epsilons[bndL] );
    }

  } catch(aol::Exception &ex) {
    ex.dump();
    return( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}

