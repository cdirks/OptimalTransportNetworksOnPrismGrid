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
#include <tpCFEPeriodicBC.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>

typedef double RealType;

int main ( int, char** ) {
  try {
    typedef tpcfe::CFEGrid< RealType, tpcfe::CFE_CD >                           GridType;

    typedef aol::SparseMatrix < RealType>                               MatrixTypeSparse;
    typedef tpcfe::CFEConfigurator < GridType, MatrixTypeSparse > ConfiguratorTypeSparse;
    typedef tpcfe::CFEMassOp  < ConfiguratorTypeSparse >                MassOpTypeSparse;
    typedef tpcfe::CFEStiffOp < ConfiguratorTypeSparse >               StiffOpTypeSparse;
    typedef tpcfe::CFEElastOp < ConfiguratorTypeSparse >               ElastOpTypeSparse;

    typedef tpcfe::CFEPeriodicHybridMatrix< GridType >                    MatrixTypeHybr;
    typedef tpcfe::CFEConfigurator < GridType, MatrixTypeHybr >     ConfiguratorTypeHybr;
    typedef tpcfe::CFEMassOp  < ConfiguratorTypeHybr >                    MassOpTypeHybr;
    typedef tpcfe::CFEStiffOp < ConfiguratorTypeHybr >                   StiffOpTypeHybr;
    typedef tpcfe::CFEElastOp < ConfiguratorTypeHybr >                   ElastOpTypeHybr;

    const int depth = 3;

    GridType grid ( depth );

    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );
    qc::BitArray<qc::QC_3D> dummy_DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

    // Initialize the grid
    grid.setDomainFrom ( levelset );

    grid.detectAndInitVirtualNodes ();

    grid.setDirichletMask( dummy_DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    tpcfe::CFEPeriodicityHandler < GridType > periodicityHandler( grid );

    { // Check whether assembled matrices for mass and stiff op coincide for both types
      bool matOK = true;

      MatrixTypeSparse massMatSparse ( grid ), stiffMatSparse ( grid );
      MatrixTypeHybr massMatHybr ( grid ), stiffMatHybr ( grid );
      {
        MassOpTypeSparse massOpSparse ( grid, aol::ONTHEFLY );
        massOpSparse.assembleAddMatrix ( massMatSparse );
        StiffOpTypeSparse stiffOpSparse ( grid, aol::ONTHEFLY );
        stiffOpSparse.assembleAddMatrix ( stiffMatSparse );

        MassOpTypeHybr massOpHybr ( grid, aol::ONTHEFLY );
        massOpHybr.assembleAddMatrix ( massMatHybr );
        StiffOpTypeHybr stiffOpHybr ( grid, aol::ONTHEFLY );
        stiffOpHybr.assembleAddMatrix ( stiffMatHybr );
      }

      qc::ScalarArray<RealType, qc::QC_3D> uSmooth ( grid ), source ( grid ), rhs ( grid ), soln0 ( grid ), soln1 ( grid );

      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        uSmooth.set ( *bit, 1.0 * grid.H() * (*bit)[2] );
      }
      grid.restrictToDomain( uSmooth ); // here, use grid's restriction

      stiffMatHybr.applyAdd ( uSmooth, rhs );                       // note: unperiodized stiffMat.
      periodicityHandler.collapsePeriodicBC ( rhs );

      rhs *= - aol::NumberTrait<RealType>::one;

      // set source term (periodic)
      qc::CoordType p1 ( grid.getWidth() / 3, grid.getWidth() / 3, grid.getWidth() / 2 ), p2 ( 2 * grid.getWidth() / 3, 2 * grid.getWidth() / 3, grid.getWidth() / 2 );
      source.set ( p1, +1.0 );
      source.set ( p2, -1.0 );

      massMatHybr.applyAdd ( source, rhs );

      periodicityHandler.periodicallyCollapseMatrix( massMatSparse );
      periodicityHandler.periodicallyCollapseMatrix( stiffMatSparse );
      periodicityHandler.restrictNonPresentDOFEntries ( stiffMatSparse, 1.0 );

      periodicityHandler.periodicallyCollapseMatrix( massMatHybr );
      periodicityHandler.periodicallyCollapseMatrix( stiffMatHybr );
      periodicityHandler.restrictNonPresentDOFEntries ( stiffMatHybr, 1.0 );

      matOK &= aol::compareOps< RealType > ( massMatSparse, massMatHybr, grid.getNumberOfNodes(), grid.getNumberOfNodes(), 1.0e-19 );
      matOK &= aol::compareOps< RealType > ( stiffMatSparse, stiffMatHybr, grid.getNumberOfNodes(), grid.getNumberOfNodes(), 1.0e-15 );

      aol::RandomAccessContainer< aol::Vector<RealType> > allOnes ( 1 );
      allOnes[0].reallocate ( grid.getNumberOfNodes() );
      allOnes[0].setAll ( 1.0 );
      periodicityHandler.restrictToPresentDOFs ( allOnes[0] );
      periodicityHandler.restrictPeriodicBC ( allOnes[0] );

      aol::RandomAccessContainer< aol::Vector<RealType> > constrVec ( 1 );
      constrVec[0].reallocate( grid.getNumberOfNodes() );
      massMatHybr.apply ( allOnes[0], constrVec[0] );                 // note: periodized massMat
      const RealType volumeFactor = constrVec[0] * allOnes[0];
      constrVec[0] /= volumeFactor;
      periodicityHandler.restrictToPresentDOFs ( constrVec[0] );  // this should not be necessary
      periodicityHandler.restrictPeriodicBC ( constrVec[0] );     // this should not be necessary

      tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::Vector<RealType> >::checkCorrectionResiduumNeutrality ( stiffMatSparse, allOnes[0] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );
      tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::Vector<RealType> >::checkCorrectionResiduumNeutrality ( stiffMatHybr, allOnes[0] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );

      tpcfe::CFEMultigridProjectAvgConstr< StiffOpTypeSparse, MassOpTypeSparse > solver0 ( grid, stiffMatSparse, periodicityHandler, 1, 2, 2, 1.0, mg::MULTIGRID_V_CYCLE, 1.0e-16, 500 );
      solver0.setVerboseMode ( 1 );
      solver0.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver0.setProjectThreshold ( 1.e-7 );
      solver0.apply ( rhs, soln0 );

      tpcfe::CFEMultigridProjectAvgConstr< StiffOpTypeHybr, MassOpTypeHybr > solver1 ( grid, stiffMatHybr, periodicityHandler, 1, 2, 2, 1.0, mg::MULTIGRID_V_CYCLE, 1.0e-16, 500 );
      solver1.setVerboseMode ( 1 );
      solver1.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver1.setProjectThreshold ( 1.e-7 );
      solver1.apply ( rhs, soln1 );

      matOK &= aol::compareOps< RealType > ( solver0.getOperatorRef(2), solver1.getOperatorRef(2), solver0.getGridRef(2).getNumberOfNodes(), solver1.getGridRef(2).getNumberOfNodes(), 1.0e-15 );
      matOK &= aol::compareOps< RealType > ( solver0.getOperatorRef(1), solver1.getOperatorRef(1), solver0.getGridRef(1).getNumberOfNodes(), solver1.getGridRef(1).getNumberOfNodes(), 1.0e-15 );

      soln1 -= soln0;
      matOK &= ( ( soln1.norm() / soln0.norm() ) < 1.0e-15 );

      if ( matOK )
        cerr << aol::color::green << "Comparison of periodic matrices for scalar problem OK" << aol::color::reset << endl;
    }

    {
      bool matOK = true;

      MatrixTypeSparse massMatSparse ( grid );
      {
        MassOpTypeSparse massOpSparse ( grid, aol::ONTHEFLY );
        massOpSparse.assembleAddMatrix ( massMatSparse );
        periodicityHandler.periodicallyCollapseMatrix ( massMatSparse );
      }
      aol::DiagonalBlockOp< RealType > BlockMassMatSparse ( massMatSparse );

      // ... set up elasticity operator
      const RealType
        E                 = 100,
        nu                = 0.33,
        lambda            = tpcfe::computeLambda ( E, nu ),
        mu                = tpcfe::computeMu ( E, nu );

      ElastOpTypeSparse elastOpSparse ( grid, lambda, mu );
      ElastOpTypeHybr elastOpHybr ( grid, lambda, mu );

      qc::MultiArray< RealType, qc::QC_3D > uSmooth ( grid ), rhs ( grid ), soln0 ( grid ), soln1 ( grid );

      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        uSmooth[ 0 ].set ( *bit, 1.0 * grid.H() * (*bit)[ 0 ] );
      }

      grid.restrictToDomain ( uSmooth );

      elastOpSparse.applyAdd ( uSmooth, rhs );
      periodicityHandler.collapsePeriodicBC ( rhs );

      rhs *= - aol::NumberTrait<RealType>::one;

      periodicityHandler.periodicallyCollapseBlockMatrix ( elastOpSparse.getBlockMatrixRef() );
      periodicityHandler.restrictNonPresentDOFEntries ( elastOpSparse.getBlockMatrixRef() );

      periodicityHandler.periodicallyCollapseBlockMatrix ( elastOpHybr.getBlockMatrixRef() );
      periodicityHandler.restrictNonPresentDOFEntries ( elastOpHybr.getBlockMatrixRef() );

      aol::RandomAccessContainer< aol::MultiVector<RealType> > neutralFunctions(3);
      for ( int i = 0; i < 3; ++i ) {
        neutralFunctions[i].reallocate ( 3, grid.getNumberOfNodes() );
        neutralFunctions[i][i].setAll ( 1.0 );
        periodicityHandler.restrictToPresentDOFs ( neutralFunctions[i] );
        periodicityHandler.restrictPeriodicBC ( neutralFunctions[i] );
      }

      aol::RandomAccessContainer< aol::MultiVector<RealType> > constrVec ( 3 );

      for ( int i = 0; i < 3; ++i ) {
        constrVec[i].reallocate ( 3, grid.getNumberOfNodes() );
        BlockMassMatSparse.apply ( neutralFunctions[i], constrVec[i] );
        const RealType volumeFactor = constrVec[i] * neutralFunctions[i];
        constrVec[i] /= volumeFactor;
        periodicityHandler.restrictToPresentDOFs ( constrVec[i] ); // this should not be necessary
        periodicityHandler.restrictPeriodicBC ( constrVec[i] );    // this should not be necessary
      }


      for ( int i = 0; i < 3; ++i ) {
        tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::MultiVector<RealType> >::checkCorrectionResiduumNeutrality ( elastOpSparse.getBlockMatrixRef(), neutralFunctions[i] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );
        tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::MultiVector<RealType> >::checkCorrectionResiduumNeutrality ( elastOpHybr.getBlockMatrixRef(), neutralFunctions[i] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );
      }

      tpcfe::CFEBlockMultigridProjectAvgConstr< ElastOpTypeSparse, MassOpTypeSparse > solver0 ( grid, elastOpSparse.getBlockMatrixRef(), periodicityHandler, 1, 3, 3, 1.0, mg::MULTIGRID_V_CYCLE, 1.0e-16, 1000 );
      solver0.setVerboseMode ( 1 );
      solver0.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver0.setProjectThreshold ( 1.0e-5 );
      solver0.apply ( rhs, soln0 );

      tpcfe::CFEBlockMultigridProjectAvgConstr< ElastOpTypeHybr, MassOpTypeHybr > solver1 ( grid, elastOpHybr.getBlockMatrixRef(), periodicityHandler, 1, 3, 3, 1.0, mg::MULTIGRID_V_CYCLE, 1.0e-16, 1000 );
      solver1.setVerboseMode ( 1 );
      solver1.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver1.setProjectThreshold ( 1.0e-8 );
      solver1.apply ( rhs, soln1 );

      soln1 -= soln0;
      matOK &= ( ( soln1.norm() / soln0.norm() ) < 1.0e-14 );

      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          matOK &= aol::compareOps< RealType > ( solver0.getOperatorRef(2).getReference(i,j), solver1.getOperatorRef(2).getReference(i,j), solver0.getGridRef(2).getNumberOfNodes(), solver1.getGridRef(2).getNumberOfNodes(), 1.0e-13 );
          matOK &= aol::compareOps< RealType > ( solver0.getOperatorRef(1).getReference(i,j), solver1.getOperatorRef(1).getReference(i,j), solver0.getGridRef(1).getNumberOfNodes(), solver1.getGridRef(1).getNumberOfNodes(), 1.0e-13 );
        }
      }

      if ( matOK )
        cerr << aol::color::green << "Comparison of periodic matrices for elasticity problem OK" << aol::color::reset << endl;
    }

  } catch ( aol::Exception &ex ) {

    ex.dump();
    return ( EXIT_FAILURE );

  }



  try {

    const tpcfe::ConstraintType CT = tpcfe::CFE_TPOS;
    typedef tpcfe::CFEGrid< RealType, CT >                                      GridType;

    typedef aol::SparseMatrix < RealType>                               MatrixTypeSparse;
    typedef tpcfe::CFEConfigurator < GridType, MatrixTypeSparse > ConfiguratorTypeSparse;
    typedef tpcfe::CFEMassOp  < ConfiguratorTypeSparse >                MassOpTypeSparse;
    typedef tpcfe::CFEStiffOpWI < ConfiguratorTypeSparse >            WStiffOpTypeSparse;

    typedef tpcfe::CFEPeriodicHybridMatrix< GridType >                    MatrixTypeHybr;
    typedef tpcfe::CFEConfigurator < GridType, MatrixTypeHybr >     ConfiguratorTypeHybr;
    typedef tpcfe::CFEMassOp  < ConfiguratorTypeHybr >                    MassOpTypeHybr;
    typedef tpcfe::CFEStiffOpWI < ConfiguratorTypeHybr >                WStiffOpTypeHybr;

    const int depth = 5;

    GridType grid ( depth );

    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );

    // Initialize the grid
    grid.addStructureFrom ( levelset );
    qc::AArray<RealType, qc::QC_3D> coeff ( grid );
    tpcfe::setCoeffForLevelset ( coeff, levelset, 1.0, 0.01 );
    grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-15, 1.0e-15 );

    tpcfe::CFEPeriodicityHandler < GridType > periodicityHandler( grid );

    { // Check whether assembled matrices for mass and stiff op coincide for both types
      bool matOK = true;

      MatrixTypeSparse massMatSparse ( grid ), stiffMatSparse ( grid );
      MatrixTypeHybr massMatHybr ( grid ), stiffMatHybr ( grid );
      {
        MassOpTypeSparse massOpSparse ( grid, aol::ONTHEFLY );
        massOpSparse.assembleAddMatrix ( massMatSparse );
        WStiffOpTypeSparse stiffOpSparse ( coeff, grid, aol::ONTHEFLY );
        stiffOpSparse.assembleAddMatrix ( stiffMatSparse );

        MassOpTypeHybr massOpHybr ( grid, aol::ONTHEFLY );
        massOpHybr.assembleAddMatrix ( massMatHybr );
        WStiffOpTypeHybr stiffOpHybr ( coeff, grid, aol::ONTHEFLY );
        stiffOpHybr.assembleAddMatrix ( stiffMatHybr );
      }

      qc::ScalarArray<RealType, qc::QC_3D> uSmooth ( grid ), rhs ( grid ), soln0 ( grid ), soln1 ( grid );

      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        uSmooth.set ( *bit, 1.0 * grid.H() * (*bit)[2] );
      }

      stiffMatHybr.applyAdd ( uSmooth, rhs );                       // note: unperiodized stiffMat.
      periodicityHandler.collapsePeriodicBC ( rhs );

      rhs *= - aol::NumberTrait<RealType>::one;

      periodicityHandler.periodicallyCollapseMatrix( massMatSparse );
      periodicityHandler.periodicallyCollapseMatrix( stiffMatSparse );
      periodicityHandler.restrictNonPresentDOFEntries ( stiffMatSparse, 1.0 );

      periodicityHandler.periodicallyCollapseMatrix( massMatHybr );
      periodicityHandler.periodicallyCollapseMatrix( stiffMatHybr );
      periodicityHandler.restrictNonPresentDOFEntries ( stiffMatHybr, 1.0 );

      matOK &= aol::rowwiseOpLinfDifference< MatrixTypeSparse, MatrixTypeHybr, RealType > ( massMatSparse, massMatHybr, grid.getNumberOfNodes() ) < 1.0e-16;
      matOK &= aol::rowwiseOpLinfDifference< MatrixTypeSparse, MatrixTypeHybr, RealType > ( stiffMatSparse, stiffMatHybr, grid.getNumberOfNodes() ) < 1.0e-16;

      aol::RandomAccessContainer< aol::Vector<RealType> > allOnes ( 1 );
      allOnes[0].reallocate ( grid.getNumberOfNodes() );
      allOnes[0].setAll ( 1.0 );
      periodicityHandler.restrictToPresentDOFs ( allOnes[0] );
      periodicityHandler.restrictPeriodicBC ( allOnes[0] );

      aol::RandomAccessContainer< aol::Vector<RealType> > constrVec ( 1 );
      constrVec[0].reallocate( grid.getNumberOfNodes() );
      massMatHybr.apply ( allOnes[0], constrVec[0] );                 // note: periodized massMat
      const RealType volumeFactor = constrVec[0] * allOnes[0];
      constrVec[0] /= volumeFactor;

      tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::Vector<RealType> >::checkCorrectionResiduumNeutrality ( stiffMatSparse, allOnes[0] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );
      tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::Vector<RealType> >::checkCorrectionResiduumNeutrality ( stiffMatHybr, allOnes[0] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );

      aol::DiagonalPreconditioner< aol::Vector<RealType> > precSparse ( stiffMatSparse );
      aol::PCGInverseProjectEqConstr< aol::Vector<RealType> > solverSparse ( stiffMatSparse, precSparse, constrVec, allOnes, 1e-16, 1000 );
      solverSparse.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solverSparse.setProjectThreshold ( 1.e-9 );
      solverSparse.apply ( rhs, soln0 );

      aol::DiagonalPreconditioner< aol::Vector<RealType> > precHybr ( stiffMatHybr );
      aol::PCGInverseProjectEqConstr< aol::Vector<RealType> > solverHybr ( stiffMatHybr, precHybr, constrVec, allOnes, 1e-16, 1000 );
      solverHybr.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solverHybr.setProjectThreshold ( 1.e-9 );
      solverHybr.apply ( rhs, soln1 );

      soln1 -= soln0;
      matOK &= ( ( soln1.norm() / soln0.norm() ) < 1.0e-15 );

      if ( matOK )
        cerr << aol::color::green << "Comparison of periodic matrices for scalar problem OK" << aol::color::reset << endl;
    }

  } catch ( aol::Exception &ex ) {

    ex.dump();
    return ( EXIT_FAILURE );

  }


  return ( EXIT_SUCCESS );

}

