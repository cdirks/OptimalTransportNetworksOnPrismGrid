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

#include <tpCFEPeriodicBC.h>

#include <configurators.h>
#include <FEOpInterface.h>
#include <quocMatrices.h>
#include <anisoStiffOps.h>
#include <shapeLevelsetGenerator.h>

#include <tpCFEStandardOp.h>
#include <tpCFEElastOp.h>

#include <tpCFELevelsets.h>

typedef double RealType;


#define DO_REPR_TEST_CUBE
#define DO_REPR_TEST_BALL
#define DO_CFE_SCALAR
#define DO_CFE_VECTOR

#define TROUBLESOME_CASE

int main ( int, char** ) {
  try {

#if defined DO_REPR_TEST_CUBE || defined DO_REPR_TEST_BALL

    typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD >                                  GridType;
    typedef aol::SparseMatrix< RealType >                                               MatrixType;
    typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;

    const int refDepth = 5;

    cerr << aol::color::blue << "Testing periodic CFE restriction and prolongation operators" << aol::color::reset << endl;

#ifdef DO_REPR_TEST_CUBE
    {
      GridType grid ( refDepth  );
      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
      levelset.setAll ( -1 );

      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes ( );

      grid.setDOFMaskFromDirichletAndDomainNodeMask();

      GridType coarseGrid ( refDepth - 1  );
      tpcfe::coarsenDomainNodeMask ( grid, coarseGrid );
      tpcfe::coarsenDirichletMask ( grid, coarseGrid );
      tpcfe::coarsenDOFMask ( grid, coarseGrid );

      tpcfe::CFEPeriodicityHandler < GridType > periodicityHandler( grid ), coarsePeriodicityHandler ( coarseGrid );

      qc::ScalarArray<RealType, qc::QC_3D> dummy ( grid ), coarseDummy ( coarseGrid ), dummyR( grid ), coarseDummyR ( coarseGrid );

      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        const RealType
          x = grid.H() * (*bit)[0],          y = grid.H() * (*bit)[1],          z = grid.H() * (*bit)[2],
          sx = sin ( aol::NumberTrait<RealType>::pi * x ),
          sy = sin ( aol::NumberTrait<RealType>::pi * y ),
          sz = sin ( aol::NumberTrait<RealType>::pi * z );

        dummyR.set ( *bit, 1 + sx * sy * sz );
      }

      for ( qc::RectangularIterator<qc::QC_3D> bit ( coarseGrid ); bit.notAtEnd(); ++bit ) {
        const RealType
          x = coarseGrid.H() * (*bit)[0],          y = coarseGrid.H() * (*bit)[1],          z = coarseGrid.H() * (*bit)[2],
          sx = sin ( aol::NumberTrait<RealType>::pi * x ),
          sy = sin ( aol::NumberTrait<RealType>::pi * y ),
          sz = sin ( aol::NumberTrait<RealType>::pi * z );

        coarseDummyR.set ( *bit, 1 + sx * sy * sz );
      }

      periodicityHandler.restrictToPresentDOFs ( dummyR );
      coarsePeriodicityHandler.restrictToPresentDOFs ( coarseDummyR );

      {
        tpcfe::CFEPeriodicProlongOp<GridType> prOp ( coarseGrid, grid, coarsePeriodicityHandler, periodicityHandler );
        prOp.apply ( coarseDummyR, dummy);

        dummy -= dummyR;

        RealType difference = dummy.getMaxAbsValue();
        tpcfe::smallOrDie ( difference, 2e-2, "PeriodicProlongation difference?", __FILE__, __LINE__ );
      }

      coarseDummyR *= 8;

      {
        tpcfe::CFEPeriodicRestrictOp<GridType> reOp ( coarseGrid, grid, coarsePeriodicityHandler, periodicityHandler, aol::ASSEMBLED );
        reOp.apply ( dummyR, coarseDummy );

        coarseDummy -= coarseDummyR;

        RealType difference = coarseDummy.getMaxAbsValue();
        tpcfe::smallOrDie ( difference, 5e-1, "PeriodicRestriction difference?", __FILE__, __LINE__ );
      }
    }
#endif

#ifdef DO_REPR_TEST_BALL
    {
      GridType grid ( refDepth  );
      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
      qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );

      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes ( );

      grid.setDOFMaskFromDirichletAndDomainNodeMask();

      GridType coarseGrid ( refDepth - 1  );
      tpcfe::coarsenDomainNodeMask ( grid, coarseGrid );
      tpcfe::coarsenDirichletMask ( grid, coarseGrid );
      tpcfe::coarsenDOFMask ( grid, coarseGrid );

      tpcfe::CFEPeriodicityHandler < GridType > periodicityHandler( grid ), coarsePeriodicityHandler ( coarseGrid );

      qc::ScalarArray<RealType, qc::QC_3D> dummy ( grid ), dummyR( grid ), coarseDummyR ( coarseGrid );

      dummyR.setAll ( 1.0 );
      periodicityHandler.restrictToPresentDOFs ( dummyR );

      coarseDummyR.setAll ( 1.0 );
      coarsePeriodicityHandler.restrictToPresentDOFs ( coarseDummyR );

      {
        tpcfe::CFEPeriodicProlongOp<GridType> prOp ( coarseGrid, grid, coarsePeriodicityHandler, periodicityHandler );
        prOp.apply ( coarseDummyR, dummy);

        dummy -= dummyR;

        RealType difference = dummy.getMaxAbsValue();
        tpcfe::smallOrDie ( difference, 1e-9, "PeriodicProlongation (ball) difference?", __FILE__, __LINE__ );
      }

      // restriction cannot be tested in a simple way.

    }
#endif

#endif

#ifdef DO_CFE_SCALAR
    cerr << aol::color::blue << "Running scalar composite FE test for periodic BC and projecting multigrid solver" << aol::color::reset << endl;
    // CFE scalar test problem on ball geometry with macroscopic part and source term
    {
      typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD >                                  GridType;
      typedef tpcfe::CFEPeriodicHybridMatrix< GridType >                           MatrixType;
      typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
      typedef tpcfe::CFEMassOp< ConfiguratorType >                                        MassOpType;
      typedef tpcfe::CFEStiffOp< ConfiguratorType >                                       StiffOpType;

      // set up grid
#ifndef TROUBLESOME_CASE
      GridType grid ( 4 );
#else
      GridType grid ( 6 );
#endif

      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

#ifndef TROUBLESOME_CASE
      qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );
#else
      qc::ShapeLevelsetGenerator<RealType>::generatePeriodicAnisoRandom3DRodsLevelset ( levelset, 1,
                                                         0.4, 0.3, 0.2,
                                                         0.0, 0.0, 0.0,
                                                         0 );
#endif

      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes ( );

      grid.setDOFMaskFromDirichletAndDomainNodeMask();
      tpcfe::CFEPeriodicityHandler < GridType > periodicityHandler( grid );

      // set MassOp and StiffOp
      MatrixType massMat ( grid ), stiffMat ( grid );
      {
        MassOpType massOp ( grid, aol::ONTHEFLY );
        massOp.assembleAddMatrix ( massMat );
        StiffOpType stiffOp ( grid, aol::ONTHEFLY );
        stiffOp.assembleAddMatrix ( stiffMat );
      }

      qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), pcgsoln ( grid ), mgsoln ( grid ), uSmooth ( grid ), source ( grid );

      // set macroscopic part
      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        uSmooth.set ( *bit, 1.0 * grid.H() * (*bit)[2] );
      }
      grid.restrictToDomain( uSmooth ); // here, use grid's restriction

      stiffMat.applyAdd ( uSmooth, rhs );                       // note: unperiodized stiffMat.
      periodicityHandler.collapsePeriodicBC ( rhs );

      rhs *= - aol::NumberTrait<RealType>::one;

      // set source term (periodic)
      qc::CoordType p1 ( grid.getWidth() / 3, grid.getWidth() / 3, grid.getWidth() / 2 ), p2 ( 2 * grid.getWidth() / 3, 2 * grid.getWidth() / 3, grid.getWidth() / 2 );
      source.set ( p1, +1.0 );
      source.set ( p2, -1.0 );

      massMat.applyAdd ( source, rhs );



      // peroidize massMat and stiffMat.
      periodicityHandler.periodicallyCollapseMatrix( massMat );
      periodicityHandler.periodicallyCollapseMatrix( stiffMat );

      periodicityHandler.restrictNonPresentDOFEntries ( stiffMat, 1.0 );


      // set neutral functions
      aol::RandomAccessContainer< aol::Vector<RealType> > allOnes ( 1 );
      allOnes[0].reallocate ( grid.getNumberOfNodes() );
      allOnes[0].setAll ( 1.0 );
      periodicityHandler.restrictToPresentDOFs ( allOnes[0] );
      periodicityHandler.restrictPeriodicBC ( allOnes[0] );

      tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::Vector<RealType> >::checkCorrectionResiduumNeutrality ( stiffMat, allOnes[0] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );


      // set average constraint
      aol::RandomAccessContainer< aol::Vector<RealType> > constrVec ( 1 );
      constrVec[0].reallocate( grid.getNumberOfNodes() );
      massMat.apply ( allOnes[0], constrVec[0] );                 // note: periodized massMat
      const RealType volumeFactor = constrVec[0] * allOnes[0];
      constrVec[0] /= volumeFactor;
      periodicityHandler.restrictToPresentDOFs ( constrVec[0] );  // this should not be necessary
      periodicityHandler.restrictPeriodicBC ( constrVec[0] );     // this should not be necessary


      // PCG solver with projection
      aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( stiffMat );
      aol::PCGInverseProjectEqConstr< aol::Vector<RealType> > solver ( stiffMat, prec, constrVec, allOnes, 1e-16, 1000 );

      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver.setProjectThreshold ( 1.e-9 );
      solver.apply ( rhs, pcgsoln );

      tpcfe::smallOrDie ( constrVec[0] * pcgsoln / pcgsoln.size(), 1e-10, "Constraint satisfied?", __FILE__, __LINE__ );

      aol::Vector<RealType> dummy ( rhs, aol::STRUCT_COPY );
      stiffMat.apply ( pcgsoln, dummy );
      dummy -= rhs;
      tpcfe::smallOrDie ( dummy.norm() / dummy.size(), 1e-10, "PCG solution solves system?", __FILE__, __LINE__ );

      // periodic extension and addition of macroscopic part
      periodicityHandler.extendPeriodicBC ( pcgsoln );

      pcgsoln += uSmooth;


      tpcfe::CFEMultigridProjectAvgConstr< StiffOpType, MassOpType > mgsolver ( grid, stiffMat, periodicityHandler, 1, 2, 2, 1.0, mg::MULTIGRID_V_CYCLE, 1.0e-24, 1000 );
      mgsolver.setVerboseMode ( 2 );
      mgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      mgsolver.setProjectThreshold ( 1.e-7 );
      mgsolver.apply ( rhs, mgsoln );

      tpcfe::smallOrDie ( constrVec[0] * mgsoln / mgsoln.size(), 1e-10, "Constraint satisfied?", __FILE__, __LINE__ );

      dummy.setZero();
      stiffMat.apply ( mgsoln, dummy );
      dummy -= rhs;
      tpcfe::smallOrDie ( dummy.norm() / dummy.size(), 1e-10, "MG solution solves system?", __FILE__, __LINE__ );

      // periodic extension and addition of macroscopic part
      periodicityHandler.extendPeriodicBC ( mgsoln );

      mgsoln += uSmooth;

      dummy = mgsoln;
      dummy -= pcgsoln;

      tpcfe::smallOrDie ( dummy.norm() / dummy.getTotalSize(), 1e-10, "PCG and MG solutions sufficiently similar?", __FILE__, __LINE__ );
    }
#endif


#ifdef DO_CFE_VECTOR
    cerr << aol::color::blue << "Running CFE elasticity test for homogenization" << aol::color::reset << endl;
    // vector-valued case (elasticity) of periodic boundary conditions (CFE) without source term
    {
      typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD >       GridType;
      typedef tpcfe::CFEPeriodicHybridMatrix< GridType >       MatrixType;
      typedef tpcfe::CFEConfigurator < GridType, MatrixType >  ConfiguratorType;
      typedef tpcfe::CFEMassOp< ConfiguratorType >             MassOpType;
      typedef tpcfe::CFEElastOp< ConfiguratorType >            ElastOpType;

#ifndef TROUBLESOME_CASE
      GridType grid ( 4 );
#else
      GridType grid ( 5 );
#endif

      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

#ifndef TROUBLESOME_CASE
      qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );
#else
      qc::ShapeLevelsetGenerator<RealType>::generatePeriodicAnisoRandom3DRodsLevelset ( levelset, 1,
                                                         0.4, 0.3, 0.2,
                                                         0.0, 0.0, 0.0,
                                                         0 );
#endif

      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes ( );

      grid.setDOFMaskFromDirichletAndDomainNodeMask();

      tpcfe::CFEPeriodicityHandler < GridType > periodicityHandler( grid );

      // set up block mass matrix ...
      MatrixType massMat ( grid );
      {
        MassOpType massOp ( grid, aol::ONTHEFLY );
        massOp.assembleAddMatrix ( massMat );
        periodicityHandler.periodicallyCollapseMatrix ( massMat );
      }
      aol::DiagonalBlockOp< RealType > BlockMassMat ( massMat );

      // ... set up elasticity operator
      const RealType
        E                 = 100,
        nu                = 0.33,
        lambda            = tpcfe::computeLambda ( E, nu ),
        mu                = tpcfe::computeMu ( E, nu );

      ElastOpType elastOp ( grid, lambda, mu );


      qc::MultiArray< RealType, qc::QC_3D > rhs ( grid ), pcgsoln ( grid ), mgsoln ( grid ), amgsoln ( grid ), uSmooth ( grid );

      // set macroscopic part
      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        uSmooth[ 0 ].set ( *bit, 1.0 * grid.H() * (*bit)[ 0 ] );
      }
      // non-oscillatory (macroscopic) part to right hand side
      grid.restrictToDomain ( uSmooth );

      elastOp.applyAdd ( uSmooth, rhs );
      periodicityHandler.collapsePeriodicBC ( rhs );

      // no source term

      rhs *= - aol::NumberTrait<RealType>::one;


      // BlockMassMat has already been collapsed
      periodicityHandler.periodicallyCollapseBlockMatrix ( elastOp.getBlockMatrixRef() );
      periodicityHandler.restrictNonPresentDOFEntries ( elastOp.getBlockMatrixRef() );

      // set neutral functions
      aol::RandomAccessContainer< aol::MultiVector<RealType> > neutralFunctions(3);
      for ( int i = 0; i < 3; ++i ) {
        neutralFunctions[i].reallocate ( 3, grid.getNumberOfNodes() );
        neutralFunctions[i][i].setAll ( 1.0 );
        periodicityHandler.restrictToPresentDOFs ( neutralFunctions[i] );
        periodicityHandler.restrictPeriodicBC ( neutralFunctions[i] );
      }

      for ( int i = 0; i < 3; ++i )
        tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::MultiVector<RealType> >::checkCorrectionResiduumNeutrality ( elastOp.getBlockMatrixRef(), neutralFunctions[i] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );


      // set average constraint
      aol::RandomAccessContainer< aol::MultiVector<RealType> > constrVec ( 3 );

      for ( int i = 0; i < 3; ++i ) {
        constrVec[i].reallocate ( 3, grid.getNumberOfNodes() );
        BlockMassMat.apply ( neutralFunctions[i], constrVec[i] );
        const RealType volumeFactor = constrVec[i] * neutralFunctions[i];
        constrVec[i] /= volumeFactor;
        periodicityHandler.restrictToPresentDOFs ( constrVec[i] ); // this should not be necessary
        periodicityHandler.restrictPeriodicBC ( constrVec[i] );    // this should not be necessary
      }

      for ( int i = 0; i < 3; ++i )
        tpcfe::smallOrDie ( constrVec[i] * rhs, 1e-10, "Constraint satisfied by uSmooth part?", __FILE__, __LINE__ );


      // CG solver with projection
      aol::BlockDiagonalPreconditioner< RealType, qc::QC_3D > prec ( elastOp.getBlockMatrixRef() );
      aol::PCGInverseProjectEqConstr< aol::MultiVector<RealType> > solver ( elastOp.getBlockMatrixRef(), prec, constrVec, neutralFunctions, 1e-16, 1000 );

      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver.setProjectThreshold ( 1.0e-5 );
      solver.apply ( rhs, pcgsoln );


      // test whether constraint is satisfied and soln solves system.
      for ( int i = 0; i < 3; ++i )
        tpcfe::smallOrDie ( constrVec[i] * pcgsoln / pcgsoln.getTotalSize(), 1e-9, "Constraint satisfied?", __FILE__, __LINE__ );

      qc::MultiArray< RealType, qc::QC_3D > dummy( pcgsoln, aol::STRUCT_COPY );
      elastOp.apply ( pcgsoln, dummy );
      dummy -= rhs;
      tpcfe::smallOrDie ( dummy.norm() / dummy.getTotalSize(), 1e-10, "Solves system?", __FILE__, __LINE__ );

      periodicityHandler.extendPeriodicBC ( pcgsoln );
      pcgsoln += uSmooth;

      // MG solver with projection
      tpcfe::CFEBlockMultigridProjectAvgConstr< ElastOpType, MassOpType > mgsolver ( grid, elastOp.getBlockMatrixRef(), periodicityHandler, 1, 3, 3, 1.0, mg::MULTIGRID_V_CYCLE, 1.0e-24, 1000 );
      mgsolver.setVerboseMode ( 2 );
      mgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      mgsolver.setProjectThreshold ( 1.0e-5 );
      mgsolver.apply ( rhs, mgsoln );

      // test whether constraint is satisfied and soln solves system.
      cerr << "Constraint satisfied? ";
      for ( int i = 0; i < 3; ++i )
        tpcfe::smallOrDie ( constrVec[i] * mgsoln / mgsoln.getTotalSize(), 1e-9, "Constraint satisfied?", __FILE__, __LINE__ );

      elastOp.apply ( mgsoln, dummy );
      dummy -= rhs;
      tpcfe::smallOrDie ( dummy.norm() / dummy.getTotalSize(), 1e-10, "Solves system?", __FILE__, __LINE__ );

      periodicityHandler.extendPeriodicBC ( mgsoln );
      mgsoln += uSmooth;

      dummy = mgsoln;
      dummy -= pcgsoln;

      tpcfe::smallOrDie ( dummy.norm() / dummy.getTotalSize(), 1e-7, "PCG and MG solutions sufficiently similar?", __FILE__, __LINE__ );
    }
#endif

  } catch ( aol::Exception &ex ) {
    ex.dump();    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
