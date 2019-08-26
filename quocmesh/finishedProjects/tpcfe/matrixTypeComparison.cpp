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

#include <shapeLevelsetGenerator.h>
#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>


void printSuccess ( const bool success ) {
  if ( success ) {
    cerr << aol::color::green << "OK." << aol::color::reset << endl;
  } else {
    cerr << aol::color::red << "<----- failed!" << endl << aol::color::reset;
  }
}

template < typename GridType, typename RealType, typename MatrixTypeSparse, typename MatrixTypeBand, typename MatrixTypeHybr, typename OpTypeSparse, typename OpTypeBand, typename OpTypeHybr >
bool scalarMatrixComparison ( const GridType &grid, MatrixTypeSparse &sparseMat, MatrixTypeBand &bandMat, MatrixTypeHybr &hybrMat, const RealType compareTolerance ) {
  bool opOK = true;

  // compare assembled matrices

  tpcfe::restrictNonDOFEntries ( grid, sparseMat, 1.0 );
  tpcfe::restrictNonDOFEntries ( grid, bandMat, 1.0 );
  tpcfe::restrictNonDOFEntries ( grid, hybrMat, 1.0 );

  //   cerr << aol::rowwiseOpLinfDifference< MatrixTypeSparse, MatrixTypeBand, RealType >( sparseMat, bandMat, grid.getNumberOfNodes() ) << " "
  //        << aol::rowwiseOpLinfDifference<   MatrixTypeBand, MatrixTypeHybr, RealType >( bandMat, hybrMat, grid.getNumberOfNodes() )  << " < " <<  compareTolerance << endl;

  opOK &= aol::rowwiseOpLinfDifference< MatrixTypeSparse, MatrixTypeBand, RealType > ( sparseMat, bandMat, grid.getNumberOfNodes() ) < compareTolerance;
  opOK &= aol::rowwiseOpLinfDifference<   MatrixTypeBand, MatrixTypeHybr, RealType > ( bandMat, hybrMat, grid.getNumberOfNodes() ) < compareTolerance;

  return ( opOK );
}


template < typename GridType, typename RealType, typename MatrixTypeSparse, typename MatrixTypeBand, typename MatrixTypeHybr, typename OpTypeSparse, typename OpTypeBand, typename OpTypeHybr >
bool scalarSolveComparison ( const GridType &grid, MatrixTypeSparse &sparseMat, MatrixTypeBand &bandMat, MatrixTypeHybr &hybrMat, const RealType compareTolerance, const RealType solutionTolerance ) {
  bool opOK = true;

  // compare multigrid solvers

  qc::ScalarArray<RealType, qc::QC_3D> orig ( grid ), rhs ( grid ), soln0 ( grid ), soln1 ( grid ), soln2 ( grid ), soln3 ( grid );
  for ( int i = 0; i < orig.size(); ++i )
    orig[i] = cos ( 1.0 * i );

  grid.restrictToDofs ( orig );

  bandMat.apply ( orig, rhs );

  grid.restrictToDofs ( rhs ); // probably unnecessary, but does no harm

  tpcfe::CFEMultigrid< OpTypeSparse, aol::ONTHEFLY, mg::ONTHEFLY_COARSENING > mgSolver0 ( grid, sparseMat, 2, 3, 3, aol::GAUSS_SEIDEL_FORWARD, 1.0, mg::MULTIGRID_V_CYCLE, 1e-16, 1000 );
  mgSolver0.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
  mgSolver0.apply ( rhs, soln0 );

  tpcfe::CFEMultigrid< OpTypeBand, aol::ASSEMBLED, mg::ONTHEFLY_COARSENING > mgSolver1 ( grid, bandMat, 3 /*for TPOS 2 not possible*/, 3, 3, aol::GAUSS_SEIDEL_FORWARD, 1.0, mg::MULTIGRID_V_CYCLE, 1e-16, 1000 );
  mgSolver1.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
  mgSolver1.apply ( rhs, soln1 );

  tpcfe::CFEMultigrid< OpTypeHybr, aol::ONTHEFLY, mg::MATRIXMULT_COARSENING > mgSolver2 ( grid, hybrMat, 2, 3, 3, aol::GAUSS_SEIDEL_FORWARD, 1.0, mg::MULTIGRID_V_CYCLE, 1e-16, 1000 );
  mgSolver2.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
  mgSolver2.apply ( rhs, soln2 );

  tpcfe::CFEMultigrid< OpTypeHybr, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgSolver3 ( grid, hybrMat, 2, 3, 3, aol::GAUSS_SEIDEL_FORWARD, 1.0, mg::MULTIGRID_V_CYCLE, 1e-16, 1000 );
  mgSolver3.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
  mgSolver3.apply ( rhs, soln3 );

  soln1 -= soln0;
  soln2 -= soln0;
  soln3 -= soln0;

  soln1.addToAll ( - soln1.sum() / grid.getNumDOFs() ); grid.restrictToDofs ( soln1 ); // subtract possible constant offset: does not hurt in MassOp case, helps in StiffOpCase
  soln2.addToAll ( - soln2.sum() / grid.getNumDOFs() ); grid.restrictToDofs ( soln2 );
  soln3.addToAll ( - soln3.sum() / grid.getNumDOFs() ); grid.restrictToDofs ( soln3 );

  //   cerr << soln1.norm() / soln0.norm() << " " << soln2.norm() / soln0.norm() << " " << soln3.norm() / soln0.norm() << " <  " <<  compareTolerance << ", " << solutionTolerance << endl;

  opOK &= ( ( soln1.norm() / soln0.norm() ) < solutionTolerance ); // different explicit level in MG solver -> can only expect solver accuracy
  opOK &= ( ( soln2.norm() / soln0.norm() ) < compareTolerance );  // should be "same" solver
  opOK &= ( ( soln3.norm() / soln0.norm() ) < compareTolerance );

  const GridType& smallGrid = mgSolver0.getGridRef ( grid.getGridDepth() - 1 );

  //   cerr << aol::rowwiseOpLinfDifference< MatrixTypeSparse, MatrixTypeBand, RealType >( mgSolver0.getOperatorRef( grid.getGridDepth() - 1 ), mgSolver1.getOperatorRef( grid.getGridDepth() - 1 ), smallGrid.getNumberOfNodes() ) << ", "
  //        << aol::rowwiseOpLinfDifference<   MatrixTypeBand, MatrixTypeHybr, RealType >( mgSolver1.getOperatorRef( grid.getGridDepth() - 1 ), mgSolver2.getOperatorRef( grid.getGridDepth() - 1 ), smallGrid.getNumberOfNodes() ) << ", "
  //        << aol::rowwiseOpLinfDifference<   MatrixTypeHybr, MatrixTypeHybr, RealType >( mgSolver2.getOperatorRef( grid.getGridDepth() - 1 ), mgSolver3.getOperatorRef( grid.getGridDepth() - 1 ), smallGrid.getNumberOfNodes() ) << " <  "
  //        << compareTolerance << endl;

  opOK &= aol::rowwiseOpLinfDifference< MatrixTypeSparse, MatrixTypeBand, RealType > ( mgSolver0.getOperatorRef ( grid.getGridDepth() - 1 ), mgSolver1.getOperatorRef ( grid.getGridDepth() - 1 ), smallGrid.getNumberOfNodes() ) < compareTolerance;
  opOK &= aol::rowwiseOpLinfDifference<   MatrixTypeBand, MatrixTypeHybr, RealType > ( mgSolver1.getOperatorRef ( grid.getGridDepth() - 1 ), mgSolver2.getOperatorRef ( grid.getGridDepth() - 1 ), smallGrid.getNumberOfNodes() ) < compareTolerance;
  opOK &= aol::rowwiseOpLinfDifference<   MatrixTypeHybr, MatrixTypeHybr, RealType > ( mgSolver2.getOperatorRef ( grid.getGridDepth() - 1 ), mgSolver3.getOperatorRef ( grid.getGridDepth() - 1 ), smallGrid.getNumberOfNodes() ) < compareTolerance;

  return ( opOK );
}


template < typename GridType, typename RealType, typename OpTypeSparse, typename OpTypeBand, typename OpTypeHybr >
bool elastSolveComparison ( const GridType &grid, OpTypeSparse &sparseOp, OpTypeBand &bandOp, OpTypeHybr &hybrOp, const RealType compareTolerance, const RealType elastSolutionTolerance ) {
  bool opOK = true;

  qc::MultiArray<RealType, 3> orig ( grid ), rhs ( grid ), soln0 ( grid ), soln1 ( grid ), soln2 ( grid );
  for ( int j = 0; j < 3; ++j ) {
    for ( int i = 0; i < orig[0].size(); ++i ) {
      orig[j][i] = cos ( 1.0 * i + j );
    }
  }

  grid.restrictToDofs ( orig );

  bandOp.getBlockMatrixRef().apply ( orig, rhs );

  // use the four different ones in the scalar case, not necessarily also in the vector-valued case
  tpcfe::CFEBlockMultigrid< OpTypeSparse, aol::ASSEMBLED, mg::ONTHEFLY_COARSENING > mgSolver0 ( grid, sparseOp.getBlockMatrixRef(), 2, 3, 3, aol::GAUSS_SEIDEL_FORWARD, 1.0, mg::MULTIGRID_V_CYCLE, 1e-16, 1000 );
  mgSolver0.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

  tpcfe::CFEBlockMultigrid< OpTypeBand, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgSolver1 ( grid, bandOp.getBlockMatrixRef(), 3 /*see above*/, 3, 3, aol::GAUSS_SEIDEL_FORWARD, 1.0, mg::MULTIGRID_V_CYCLE, 1e-16, 1000 );
  mgSolver1.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

  tpcfe::CFEBlockMultigrid< OpTypeHybr, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgSolver2 ( grid, hybrOp.getBlockMatrixRef(), 2, 3, 3, aol::GAUSS_SEIDEL_FORWARD, 1.0, mg::MULTIGRID_V_CYCLE, 1e-16, 1000 );
  mgSolver2.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

  mgSolver0.apply ( rhs, soln0 );
  mgSolver1.apply ( rhs, soln1 );
  mgSolver2.apply ( rhs, soln2 );

  soln1 -= soln0;
  soln2 -= soln0;

  //   cerr << soln1.norm() / soln0.norm() << ", " << soln2.norm() / soln0.norm() << " <  " << elastSolutionTolerance << ", " << compareTolerance << endl;

  opOK &= ( ( soln1.norm() / soln0.norm() ) < elastSolutionTolerance );
  opOK &= ( ( soln2.norm() / soln0.norm() ) < compareTolerance );


  for ( short i = 0; i < 3; ++i ) {
    for ( short j = 0; j < 3; ++j ) {
      const int depth = grid.getGridDepth();
      const GridType& smallGrid = mgSolver1.getGridRef ( depth - 1 );
      typedef aol::SparseMatrix < RealType>       MatrixTypeSparse;
      typedef tpcfe::CFEBandMatrix< GridType >    MatrixTypeBand;
      typedef tpcfe::CFEHybridMatrix< GridType >  MatrixTypeHybr;

      //       cerr << aol::rowwiseOpLinfDifference< MatrixTypeSparse, MatrixTypeBand, RealType >( mgSolver0.getOperatorRef( depth - 1 ).getReference(i, j),
      //                                                                                           mgSolver1.getOperatorRef( depth - 1 ).getReference(i, j),
      //                                                                                           smallGrid.getNumberOfNodes() )
      //            << ", "
      //            << aol::rowwiseOpLinfDifference< MatrixTypeBand,   MatrixTypeHybr, RealType >( mgSolver1.getOperatorRef( depth - 1 ).getReference(i, j),
      //                                                                                           mgSolver2.getOperatorRef( depth - 1 ).getReference(i, j),
      //                                                                                           smallGrid.getNumberOfNodes() )
      //            << " < " << compareTolerance  << endl;

      opOK &= aol::rowwiseOpLinfDifference< MatrixTypeSparse, MatrixTypeBand, RealType > ( mgSolver0.getOperatorRef ( depth - 1 ).getReference ( i, j ),
                                                                                           mgSolver1.getOperatorRef ( depth - 1 ).getReference ( i, j ),
                                                                                           smallGrid.getNumberOfNodes() )                              < compareTolerance;
      opOK &= aol::rowwiseOpLinfDifference< MatrixTypeBand,   MatrixTypeHybr, RealType > ( mgSolver1.getOperatorRef ( depth - 1 ).getReference ( i, j ),
                                                                                           mgSolver2.getOperatorRef ( depth - 1 ).getReference ( i, j ),
                                                                                           smallGrid.getNumberOfNodes() )                              < compareTolerance;
    }
  }

  return ( opOK );
}


template < typename GridType, typename RealType, typename OpTypeSparse, typename OpTypeBand, typename OpTypeHybr >
bool scalarOpComparison ( const GridType &grid, const RealType compareTolerance, const RealType solutionTolerance ) {
  typedef aol::SparseMatrix < RealType>       MatrixTypeSparse;
  typedef tpcfe::CFEBandMatrix< GridType >    MatrixTypeBand;
  typedef tpcfe::CFEHybridMatrix< GridType >  MatrixTypeHybr;

  MatrixTypeSparse sparseMat ( grid );
  {
    OpTypeSparse opSparse ( grid, aol::ONTHEFLY );
    opSparse._quietMode = true;
    opSparse.assembleAddMatrix ( sparseMat );
  }

  MatrixTypeBand bandMat ( grid );
  {
    OpTypeBand opBand ( grid, aol::ONTHEFLY );
    opBand._quietMode = true;
    opBand.assembleAddMatrix ( bandMat );
  }

  MatrixTypeHybr hybrMat ( grid );
  {
    OpTypeHybr opHybr ( grid, aol::ONTHEFLY );
    opHybr._quietMode = true;
    opHybr.assembleAddMatrix ( hybrMat );
  }

  return ( scalarMatrixComparison< GridType, RealType, MatrixTypeSparse, MatrixTypeBand, MatrixTypeHybr, OpTypeSparse, OpTypeBand, OpTypeHybr > ( grid, sparseMat, bandMat, hybrMat, compareTolerance ) &&
           scalarSolveComparison < GridType, RealType, MatrixTypeSparse, MatrixTypeBand, MatrixTypeHybr, OpTypeSparse, OpTypeBand, OpTypeHybr > ( grid, sparseMat, bandMat, hybrMat, compareTolerance, solutionTolerance ) );
}


template < typename GridType, typename RealType, typename OpTypeSparse, typename OpTypeBand, typename OpTypeHybr >
bool scalarWeightedOpComparison ( const qc::AArray<RealType, qc::QC_3D> &coeff, const GridType &grid, const RealType compareTolerance, const RealType solutionTolerance ) {
  typedef aol::SparseMatrix < RealType>       MatrixTypeSparse;
  typedef tpcfe::CFEBandMatrix< GridType >    MatrixTypeBand;
  typedef tpcfe::CFEHybridMatrix< GridType >  MatrixTypeHybr;

  MatrixTypeSparse sparseMat ( grid );
  {
    OpTypeSparse opSparse ( coeff, grid, aol::ONTHEFLY );
    opSparse._quietMode = true;
    opSparse.assembleAddMatrix ( sparseMat );
  }

  MatrixTypeBand bandMat ( grid );
  {
    OpTypeBand opBand ( coeff, grid, aol::ONTHEFLY );
    opBand._quietMode = true;
    opBand.assembleAddMatrix ( bandMat );
  }

  MatrixTypeHybr hybrMat ( grid );
  {
    OpTypeHybr opHybr ( coeff, grid, aol::ONTHEFLY );
    opHybr._quietMode = true;
    opHybr.assembleAddMatrix ( hybrMat );
  }

  return ( scalarMatrixComparison< GridType, RealType, MatrixTypeSparse, MatrixTypeBand, MatrixTypeHybr, OpTypeSparse, OpTypeBand, OpTypeHybr > ( grid, sparseMat, bandMat, hybrMat, compareTolerance ) &&
           scalarSolveComparison < GridType, RealType, MatrixTypeSparse, MatrixTypeBand, MatrixTypeHybr, OpTypeSparse, OpTypeBand, OpTypeHybr > ( grid, sparseMat, bandMat, hybrMat, compareTolerance, solutionTolerance ) );
}


template < typename GridType, typename RealType, typename OpTypeSparse, typename OpTypeBand, typename OpTypeHybr >
bool elastOpComparison ( const GridType &grid, const RealType compareTolerance, const RealType elastSolutionTolerance ) {

  bool elastOK = true;
  const RealType lambda = tpcfe::computeLambda ( 1.0, 0.33 );
  const RealType mu     = tpcfe::computeMu    ( 1.0, 0.33 );

  OpTypeSparse sparseElastOp ( grid, lambda, mu );
  OpTypeBand bandElastOp ( grid, lambda, mu );
  OpTypeHybr hybrElastOp ( grid, lambda, mu );

  typedef aol::SparseMatrix < RealType>       MatrixTypeSparse;
  typedef tpcfe::CFEBandMatrix< GridType >    MatrixTypeBand;
  typedef tpcfe::CFEHybridMatrix< GridType >  MatrixTypeHybr;

  sparseElastOp.restrictNonDOFEntries();
  bandElastOp.restrictNonDOFEntries();

  hybrElastOp.restrictNonDOFEntries();

  for ( int i = 0; i < 3; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      elastOK &= scalarMatrixComparison<GridType, RealType, MatrixTypeSparse, MatrixTypeBand, MatrixTypeHybr, OpTypeSparse, OpTypeBand, OpTypeHybr > ( grid, sparseElastOp.getBlockMatrixRef().getReference ( i, j ),
                                                                                                                                                       bandElastOp.getBlockMatrixRef().getReference ( i, j ),
                                                                                                                                                       hybrElastOp.getBlockMatrixRef().getReference ( i, j ),
                                                                                                                                                       compareTolerance );
    }
  }

  elastOK &= elastSolveComparison< GridType, RealType, OpTypeSparse, OpTypeBand, OpTypeHybr > ( grid, sparseElastOp, bandElastOp, hybrElastOp, compareTolerance, elastSolutionTolerance );

  return ( elastOK );
}


template < typename GridType, typename RealType, typename OpTypeSparse, typename OpTypeBand, typename OpTypeHybr >
bool weightedElastOpComparison ( const qc::AArray<RealType, qc::QC_3D> &coeff, const GridType &grid, const RealType compareTolerance, const RealType elastSolutionTolerance ) {

  bool elastOK = true;
  const RealType lambda = tpcfe::computeLambda ( 1.0, 0.33 );
  const RealType mu     = tpcfe::computeMu    ( 1.0, 0.33 );

  OpTypeSparse sparseElastOp ( coeff, grid, lambda, mu );
  OpTypeBand bandElastOp ( coeff, grid, lambda, mu );
  OpTypeHybr hybrElastOp ( coeff, grid, lambda, mu );

  typedef aol::SparseMatrix < RealType>       MatrixTypeSparse;
  typedef tpcfe::CFEBandMatrix< GridType >    MatrixTypeBand;
  typedef tpcfe::CFEHybridMatrix< GridType >  MatrixTypeHybr;

  sparseElastOp.restrictNonDOFEntries();
  bandElastOp.restrictNonDOFEntries();

  hybrElastOp.restrictNonDOFEntries();

  for ( int i = 0; i < 3; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      elastOK &= scalarMatrixComparison<GridType, RealType, MatrixTypeSparse, MatrixTypeBand, MatrixTypeHybr, OpTypeSparse, OpTypeBand, OpTypeHybr > ( grid, sparseElastOp.getBlockMatrixRef().getReference ( i, j ),
                                                                                                                                                       bandElastOp.getBlockMatrixRef().getReference ( i, j ),
                                                                                                                                                       hybrElastOp.getBlockMatrixRef().getReference ( i, j ),
                                                                                                                                                       compareTolerance );
    }
  }

  elastOK &= elastSolveComparison< GridType, RealType, OpTypeSparse, OpTypeBand, OpTypeHybr > ( grid, sparseElastOp, bandElastOp, hybrElastOp, compareTolerance, elastSolutionTolerance );

  return ( elastOK );
}


template < typename RealType >
bool complicatedDomainTest ( const int depth, const bool useDirichlet, const RealType compareTolerance, const RealType solutionTolerance, const RealType elastSolutionTolerance ) {
  bool ret = true;

  // Check whether mass and stiff op work with both types of matrix without Dirichlet BCs
  const tpcfe::ConstraintType CT = tpcfe::CFE_CD;
  typedef tpcfe::CFEGrid< RealType, CT >      GridType;
  typedef aol::SparseMatrix < RealType>       MatrixTypeSparse;
  typedef tpcfe::CFEBandMatrix< GridType >    MatrixTypeBand;
  typedef tpcfe::CFEHybridMatrix< GridType >  MatrixTypeHybr;

  typedef tpcfe::CFEConfigurator < GridType, MatrixTypeSparse >  ConfiguratorTypeSparse;
  typedef tpcfe::CFEConfigurator < GridType, MatrixTypeBand >    ConfiguratorTypeBand;
  typedef tpcfe::CFEConfigurator < GridType, MatrixTypeHybr >    ConfiguratorTypeHybr;

  GridType grid ( depth );

  qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
  qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );

  qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
  if ( useDirichlet ) {
    for ( int i = 0; i < grid.getNumX(); ++i ) {
      for ( int j = 0; j < grid.getNumY(); ++j ) {
        DirichletMask.set ( i, j, 0, true );
        DirichletMask.set ( i, j, grid.getWidth() - 1, true );
      }
    }
  } else {
    DirichletMask.setAll ( false );
  }

  // Initialize the grid
  grid.setDomainFrom ( levelset );
  grid.detectAndInitVirtualNodes();

  grid.setDirichletMask ( DirichletMask );
  grid.setDOFMaskFromDirichletAndDomainNodeMask();

  cerr << "Comparing Mass Matrices ... " << tpcfe::ConstraintTypeNames[CT] << endl;
  typedef tpcfe::CFEMassOp  < ConfiguratorTypeSparse >                                     MassOpTypeSparse;
  typedef tpcfe::CFEMassOp  < ConfiguratorTypeBand >                                       MassOpTypeBand;
  typedef tpcfe::CFEMassOp  < ConfiguratorTypeHybr >                                       MassOpTypeHybr;
  const bool mmOk = scalarOpComparison < GridType, RealType, MassOpTypeSparse, MassOpTypeBand, MassOpTypeHybr > ( grid, compareTolerance, solutionTolerance );
  ret = ( ret && mmOk );
  printSuccess ( mmOk );

  cerr << "Comparing Stiff Matrices ... " << tpcfe::ConstraintTypeNames[CT] << endl;
  typedef tpcfe::CFEStiffOp  < ConfiguratorTypeSparse >                                    StiffOpTypeSparse;
  typedef tpcfe::CFEStiffOp  < ConfiguratorTypeBand >                                      StiffOpTypeBand;
  typedef tpcfe::CFEStiffOp  < ConfiguratorTypeHybr >                                      StiffOpTypeHybr;
  const bool smOk = scalarOpComparison < GridType, RealType, StiffOpTypeSparse, StiffOpTypeBand, StiffOpTypeHybr > ( grid, compareTolerance, solutionTolerance );
  ret = ( ret && smOk );
  printSuccess ( smOk );

  if ( useDirichlet ) { // always do this?
    cerr << "Comparing ElastOp " << tpcfe::ConstraintTypeNames[CT] << endl;
    typedef tpcfe::CFEElastOp < ConfiguratorTypeSparse >                                     ElastOpTypeSparse;
    typedef tpcfe::CFEElastOp < ConfiguratorTypeBand >                                       ElastOpTypeBand;
    typedef tpcfe::CFEElastOp < ConfiguratorTypeHybr >                                       ElastOpTypeHybr;

    const bool elOk = elastOpComparison < GridType, RealType, ElastOpTypeSparse, ElastOpTypeBand, ElastOpTypeHybr > ( grid, compareTolerance, elastSolutionTolerance );
    ret = ( ret && elOk );
    printSuccess ( elOk ) ;
  }

  return ( ret );
}


template < typename RealType, tpcfe::ConstraintType CT >
bool jumpingCoefficientTest ( const int depth, const bool useDirichlet, const RealType compareTolerance, const RealType solutionTolerance ) {
  bool ret = true;

  typedef tpcfe::CFEGrid< RealType, CT >      GridType;
  typedef aol::SparseMatrix < RealType>       MatrixTypeSparse;
  typedef tpcfe::CFEBandMatrix< GridType >    MatrixTypeBand;
  typedef tpcfe::CFEHybridMatrix< GridType >  MatrixTypeHybr;

  typedef tpcfe::CFEConfigurator < GridType, MatrixTypeSparse >  ConfiguratorTypeSparse;
  typedef tpcfe::CFEConfigurator < GridType, MatrixTypeBand >    ConfiguratorTypeBand;
  typedef tpcfe::CFEConfigurator < GridType, MatrixTypeHybr >    ConfiguratorTypeHybr;

  GridType grid ( depth );

  qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
  qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );

  qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
  if ( useDirichlet ) {
    for ( int i = 0; i < grid.getNumX(); ++i ) {
      for ( int j = 0; j < grid.getNumY(); ++j ) {
        DirichletMask.set ( i, j, 0, true );
        DirichletMask.set ( i, j, grid.getWidth() - 1, true );
      }
    }
  } else {
    DirichletMask.setAll ( false );
  }

  qc::AArray<RealType, qc::QC_3D> coeff ( grid );
  for ( int i = 0; i < coeff.size(); ++i )
    coeff[i] = ( ( levelset[i] < 0 )  ?  1.0  :  2.0 );

  // Initialize the grid
  grid.addStructureFrom ( levelset );

  grid.detectAndInitVirtualNodes ( coeff );

  grid.setDirichletMask ( DirichletMask );
  grid.setDOFMaskFromDirichletAndDomainNodeMask();

  cerr << "Comparing Mass Matrices ... " << tpcfe::ConstraintTypeNames[CT] << endl;
  typedef tpcfe::CFEMassOp  < ConfiguratorTypeSparse >                                     MassOpTypeSparse;
  typedef tpcfe::CFEMassOp  < ConfiguratorTypeBand >                                       MassOpTypeBand;
  typedef tpcfe::CFEMassOp  < ConfiguratorTypeHybr >                                       MassOpTypeHybr;
  const bool mmOk = scalarOpComparison < GridType, RealType, MassOpTypeSparse, MassOpTypeBand, MassOpTypeHybr > ( grid, compareTolerance, solutionTolerance );
  ret = ( ret && mmOk );
  printSuccess ( mmOk );

  cerr << "Comparing Stiff Matrices ... " << tpcfe::ConstraintTypeNames[CT] << endl;
  typedef tpcfe::CFEStiffOp  < ConfiguratorTypeSparse >                                    StiffOpTypeSparse;
  typedef tpcfe::CFEStiffOp  < ConfiguratorTypeBand >                                      StiffOpTypeBand;
  typedef tpcfe::CFEStiffOp  < ConfiguratorTypeHybr >                                      StiffOpTypeHybr;
  const bool smOk = scalarOpComparison < GridType, RealType, StiffOpTypeSparse, StiffOpTypeBand, StiffOpTypeHybr > ( grid, compareTolerance, solutionTolerance );
  ret = ( ret && smOk );
  printSuccess ( smOk );

  cerr << "Comparing Weighted Stiff Matrices ... " << tpcfe::ConstraintTypeNames[CT] << endl;
  typedef tpcfe::CFEStiffOpWI  < ConfiguratorTypeSparse >                                  WStiffOpTypeSparse;
  typedef tpcfe::CFEStiffOpWI  < ConfiguratorTypeBand >                                    WStiffOpTypeBand;
  typedef tpcfe::CFEStiffOpWI  < ConfiguratorTypeHybr >                                    WStiffOpTypeHybr;
  const bool wsmOk = scalarWeightedOpComparison < GridType, RealType, WStiffOpTypeSparse, WStiffOpTypeBand, WStiffOpTypeHybr > ( coeff, grid, compareTolerance, solutionTolerance );
  ret = ( ret && wsmOk );
  printSuccess ( wsmOk );

  if ( useDirichlet && CT == tpcfe::CFE_CD ) { // always do this?
    cerr << "Comparing ElastOp " << tpcfe::ConstraintTypeNames[CT] << endl;
    typedef tpcfe::CFEElastOp < ConfiguratorTypeSparse >                                   ElastOpTypeSparse;
    typedef tpcfe::CFEElastOp < ConfiguratorTypeBand >                                     ElastOpTypeBand;
    typedef tpcfe::CFEElastOp < ConfiguratorTypeHybr >                                     ElastOpTypeHybr;

    const bool elOk = elastOpComparison < GridType, RealType, ElastOpTypeSparse, ElastOpTypeBand, ElastOpTypeHybr > ( grid, compareTolerance, solutionTolerance );
    ret = ( ret && elOk );
    printSuccess ( elOk ) ;
  }

  if ( useDirichlet && CT == tpcfe::CFE_TPOS ) { // always do this?
    cerr << "No WeightedElastOpComparison yet for CFE_TPOS" << endl;
  }

  return ( ret );
}


int main ( int, char** ) {
  typedef double RealType;

  const int depth = 4;

  bool allOk = true;

  aol::StopWatch timer;
  timer.start();

  allOk = ( allOk && complicatedDomainTest< RealType > ( depth, false, 1.0e-15, 1.0e-07, -1.0e-0 ) );
  allOk = ( allOk && complicatedDomainTest< RealType > ( depth, true , 3.0e-15, 2.0e-07, 3.0e-06 ) );

  allOk = ( allOk && jumpingCoefficientTest< RealType, tpcfe::CFE_LIEHR > ( depth, false, 2.0e-15, 7.0e-08 ) );
  allOk = ( allOk && jumpingCoefficientTest< RealType, tpcfe::CFE_LIEHR > ( depth, true,  2.0e-15, 1.0e-07 ) );

  allOk = ( allOk && jumpingCoefficientTest< RealType, tpcfe::CFE_TPOS  > ( depth, false, 2.0e-15, 1.0e-07 ) );
  allOk = ( allOk && jumpingCoefficientTest< RealType, tpcfe::CFE_TPOS  > ( depth, true , 2.0e-15, 4.0e-08 ) );

  timer.stop();

  cerr << "Comparison took " << timer.elapsedCpuTime() << " seconds (cputime): " << 20086.8 / timer.elapsedCpuTime() << " nupsi." << endl;

  if ( allOk ) {
    cerr << aol::color::ok
         << " ============== " << endl
         << " ===   OK   === " << endl
         << " ============== " << endl
         << aol::color::reset;
    return ( EXIT_SUCCESS );
  } else {
    cerr << aol::color::error
         << " ============== " << endl
         << " === FAILED === " << endl
         << " ============== " << endl
         << aol::color::reset;
    return ( EXIT_FAILURE );
  }
}
