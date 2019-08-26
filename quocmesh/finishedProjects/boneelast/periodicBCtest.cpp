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


#define DO_POISSON
#define DO_MFE_SCALAR
#define DO_CFE_SCALAR
#define DO_CFE_VECTOR


int main ( int, char** ) {
  try {

#ifdef DO_POISSON
    // simplest case: Poisson problem made uniquely solvable by enforcing integral average of unknown.
    cerr << aol::color::invert << "Running Poisson test for Projecting solver" << aol::color::reset << endl;
    {
      typedef aol::SparseMatrix< RealType > MatrixType;

      qc::GridDefinition grid ( 4, qc::QC_3D );

      MatrixType massMat ( grid ), stiffMat ( grid );

      {
        aol::MassOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > > massOp ( grid, aol::ONTHEFLY );
        massOp.assembleAddMatrix ( massMat );
        aol::StiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > > stiffOp ( grid, aol::ONTHEFLY );
        stiffOp.assembleAddMatrix( stiffMat );
      }

      qc::ScalarArray<RealType, qc::QC_3D> dummy ( grid ), rhs ( grid ), soln ( grid );

      // construct right hand side for the system of equations
      for ( qc::RectangularIterator<qc::QC_3D> bit ( dummy ); bit.notAtEnd(); ++bit ) {
        dummy.set( *bit, sin ( aol::NumberTrait<RealType>::pi * (*bit)[0] / (grid.getNumX()-1) ) * sin ( 0.5 * aol::NumberTrait<RealType>::pi * (*bit)[1] / (grid.getNumY() - 1) ));
      }

      aol::RandomAccessContainer< aol::Vector<RealType> > allOnes ( 1 );
      allOnes[0].reallocate ( grid.getNumberOfNodes() );
      allOnes[0].setAll ( 1.0 );

      stiffMat.apply ( dummy, rhs );


      // set average constraint
      aol::RandomAccessContainer< aol::Vector<RealType> > constrVec ( 1 );
      constrVec[0].resize( grid.getNumberOfNodes() );
      massMat.apply ( allOnes[0], constrVec[0] );

      tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::Vector<RealType> >::checkCorrectionResiduumNeutrality ( stiffMat, allOnes[0] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );

      // CG solver with projection
      // aol::CGInverseProjectEqConstr< aol::Vector<RealType> > solver ( stiffMat, constrVec, allOnes );
      aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( stiffMat );
      aol::PCGInverseProjectEqConstr< aol::Vector<RealType> > solver ( stiffMat, prec, constrVec, allOnes );

      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver.setProjectThreshold ( 1.0e-7 );
      solver.apply ( rhs, soln );

      // test whether constraint is satisfied and soln solves system.
      tpcfe::smallOrDie ( constrVec[0] * soln / soln.size(), 1e-10, "Constraint satisfied?", __FILE__, __LINE__ );

      stiffMat.apply ( soln, dummy );
      dummy -= rhs;
      tpcfe::smallOrDie ( dummy.norm() / dummy.size(), 1e-10, "Solves system?", __FILE__, __LINE__ );

    }
#endif


#ifdef DO_MFE_SCALAR
    cerr << aol::color::invert << "Running scalar multilinearFE test for periodic BC and projecting solver" << aol::color::reset << endl;
    // multilinear FE test problem on full cube with source term and macroscopic part
    {
      typedef qc::GridDefinition                                                                                GridType;
      typedef aol::SparseMatrix< RealType >                                                                     MatrixType;
      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > ConfiguratorType;
      typedef aol::MassOp< ConfiguratorType >                                                                   MassOpType;
      typedef aol::StiffOp< ConfiguratorType >                                                                  StiffOpType;

      GridType grid ( 4, qc::QC_3D );

      qc::QuocPeriodicityHandler < GridType, RealType, qc::QC_3D > periodicityHandler( grid );

      // set MassOp and StiffOp
      MatrixType massMat ( grid ), stiffMat ( grid );
      {
        MassOpType massOp ( grid, aol::ONTHEFLY );
        massOp.assembleAddMatrix ( massMat );
        StiffOpType stiffOp ( grid, aol::ONTHEFLY );
        stiffOp.assembleAddMatrix ( stiffMat );
      }


      qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), soln ( grid ), uSmooth ( grid ), source ( grid );

      // set macroscopic part
      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        uSmooth.set ( *bit, 1.0 * grid.H() * (*bit)[2] );
      }

      stiffMat.applyAdd ( uSmooth, rhs );   // note: unperiodized stiffMat.
      periodicityHandler.collapsePeriodicBC ( rhs );

      rhs *= - aol::NumberTrait<RealType>::one;

      // set source term (periodic)
      qc::CoordType p1 ( grid.getWidth() / 3, grid.getWidth() / 3, grid.getWidth() / 3 ), p2 = static_cast<short>( 2 ) * p1;
      source.set ( p1, +1.0 );
      source.set ( p2, -1.0 );

      massMat.applyAdd ( source, rhs );


      // peroidize massMat and stiffMat.
      periodicityHandler.periodicallyCollapseMatrix( massMat );
      periodicityHandler.periodicallyCollapseMatrix( stiffMat );


      // set neutral functions
      aol::RandomAccessContainer< aol::Vector<RealType> > allOnes ( 1 );
      allOnes[0].reallocate ( grid.getNumberOfNodes() );
      allOnes[0].setAll ( 1.0 );
      periodicityHandler.restrictPeriodicBC ( allOnes[0] );


      // set average constraint
      aol::RandomAccessContainer< aol::Vector<RealType> > constrVec ( 1 );
      constrVec[0].resize( grid.getNumberOfNodes() );
      massMat.apply ( allOnes[0], constrVec[0] );  // note: periodized massMat
      const RealType volumeFactor = constrVec[0] * allOnes[0];
      constrVec[0] /= volumeFactor;
      periodicityHandler.restrictPeriodicBC ( constrVec[0] );  // necessary?


      tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::Vector<RealType> >::checkCorrectionResiduumNeutrality ( stiffMat, allOnes[0] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );

      // PCG solver with projection
      aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( stiffMat );
      aol::PCGInverseProjectEqConstr< aol::Vector<RealType> > solver ( stiffMat, prec, constrVec, allOnes, 1e-16, 1000 );

      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver.setProjectThreshold ( 1.0e-6 );
      solver.apply ( rhs, soln );


      // test whether constraint is satisfied and soln solves system.
      tpcfe::smallOrDie ( constrVec[0] * soln / soln.size(), 1e-10, "Constraint satisfied?", __FILE__, __LINE__ );

      aol::Vector<RealType> dummy ( rhs, aol::STRUCT_COPY );
      stiffMat.apply ( soln, dummy ); // misuse of source as dummy vector
      dummy -= rhs;
      tpcfe::smallOrDie ( dummy.norm() / dummy.size(), 1e-10, "Solves system?", __FILE__, __LINE__ );


      // periodic extension and addition of macroscopic part
      periodicityHandler.extendPeriodicBC ( soln );
      soln += uSmooth;

    }
#endif


#ifdef DO_CFE_SCALAR
    cerr << aol::color::invert << "Running scalar composite FE test for periodic BC and projecting solver" << aol::color::reset << endl;
    // CFE scalar test problem on ball geometry with macroscopic part and source term
    {
      typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD >        GridType;
      typedef aol::SparseMatrix< RealType >                     MatrixType;
      typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
      typedef tpcfe::CFEMassOp< ConfiguratorType >              MassOpType;
      typedef tpcfe::CFEStiffOp< ConfiguratorType >             StiffOpType;

      // set up grid
      GridType grid ( 4  );
      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
      qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );

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

      qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), soln ( grid ), uSmooth ( grid ), source ( grid );

      // set macroscopic part
      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        uSmooth.set ( *bit, 1.0 * grid.H() * (*bit)[2] );
      }
      grid.restrictToDomain( uSmooth ); // here, use grid's restriction

      stiffMat.applyAdd ( uSmooth, rhs );                       // note: unperiodized stiffMat.
      periodicityHandler.collapsePeriodicBC ( rhs );

      // set source term (periodic)
      qc::CoordType p1 ( grid.getWidth() / 3, grid.getWidth() / 3, grid.getWidth() / 2 ), p2 ( 2 * grid.getWidth() / 3, 2 * grid.getWidth() / 3, grid.getWidth() / 2 );
      source.set ( p1, +1.0 );
      source.set ( p2, -1.0 );

      massMat.applyAdd ( source, rhs );

      rhs *= - aol::NumberTrait<RealType>::one;


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
      solver.apply ( rhs, soln );


      tpcfe::smallOrDie ( constrVec[0] * soln / soln.size(), 1e-10, "Constraint satisfied?", __FILE__, __LINE__ );

      aol::Vector<RealType> dummy ( rhs, aol::STRUCT_COPY );
      stiffMat.apply ( soln, dummy );
      dummy -= rhs;
      tpcfe::smallOrDie ( dummy.norm() / dummy.size(), 1e-10, "Solves system?", __FILE__, __LINE__ );


      // periodic extension and addition of macroscopic part
      periodicityHandler.extendPeriodicBC ( soln );

      soln += uSmooth;

    }
#endif


#ifdef DO_CFE_VECTOR
    cerr << aol::color::invert << "Running CFE elasticity test for homogenization" << aol::color::reset << endl;
    // vector-valued case (elasticity) of periodic boundary conditions (CFE) without source term
    {
      typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD >       GridType;
      typedef tpcfe::CFEPeriodicHybridMatrix< GridType >       MatrixType;
      typedef tpcfe::CFEConfigurator < GridType, MatrixType >  ConfiguratorType;
      typedef tpcfe::CFEMassOp< ConfiguratorType >             MassOpType;
      typedef tpcfe::CFEElastOp< ConfiguratorType >            ElastOpType;

      GridType grid ( 4 );

      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
      qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );

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


      qc::MultiArray< RealType, qc::QC_3D > rhs ( grid ), soln ( grid ), uSmooth ( grid );

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
      solver.apply ( rhs, soln );


      // test whether constraint is satisfied and soln solves system.
      for ( int i = 0; i < 3; ++i )
        tpcfe::smallOrDie ( constrVec[i] * soln / soln.getTotalSize(), 1e-10, "Constraint satisfied?", __FILE__, __LINE__ );

      qc::MultiArray< RealType, qc::QC_3D > dummy( soln, aol::STRUCT_COPY );
      elastOp.apply ( soln, dummy );
      dummy -= rhs;
      tpcfe::smallOrDie ( dummy.norm() / dummy.getTotalSize(), 1e-10, "Solves system?", __FILE__, __LINE__ );

      periodicityHandler.extendPeriodicBC ( soln );
      soln += uSmooth;

    }
#endif

  } catch ( aol::Exception &ex ) {
    ex.dump();    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
