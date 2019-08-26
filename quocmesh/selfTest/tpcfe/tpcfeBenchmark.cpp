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

#include <aol.h>
#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>
#include <shapeLevelsetGenerator.h>

typedef double RealType;


int main ( int argc, char** argv ) {
  bool bench = false;
  string resultfilename = "/dev/null";
  int jcLevel;
  int cdLevel;

  if ( aol::checkForBenchmarkArguments ( argc, argv, resultfilename ) ) {
    bench = true;

    jcLevel = 5;
    cdLevel = 6;
  }
  else if ( argc == 4 ) {
    bench = true;

    jcLevel = atoi( argv[1] );
    cdLevel = atoi( argv[2] );
    resultfilename = argv[3];
  }
  else {
    jcLevel = 2;
    cdLevel = 2;
  }

  try {

    aol::StopWatch timer;
    timer.start();

    { // jump coeff heat cond test

      typedef double RealType;

      static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;
      typedef tpcfe::CFEGrid < RealType, AT >                   GridType;
      typedef tpcfe::CFEHybridMatrix<GridType>                  MatrixType;

      typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
      typedef tpcfe::CFEMassOp  < ConfiguratorType >            MassOpType;
      typedef tpcfe::CFEStiffOpWI < ConfiguratorType >          WStiffOpType;

      GridType grid ( jcLevel );

      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
      qc::ShapeLevelsetGenerator<RealType>::generateConeLevelset ( levelset, 0.231, 0.421 );

      grid.addStructureFrom ( levelset );

      qc::AArray< RealType, qc::QC_3D > coeffLambda ( grid );
      tpcfe::setCoeffForLevelset ( coeffLambda, levelset, 1.0, 23.0 );
      grid.detectAndInitVirtualNodes ( coeffLambda );

      MatrixType massMat ( grid ), stiffMat ( grid ), sysMat ( grid );
      {
        MassOpType massOp ( grid, aol::ONTHEFLY );
        WStiffOpType stiffOp ( coeffLambda, grid, aol::ONTHEFLY );

        massOp.assembleAddMatrix ( massMat );
        stiffOp.assembleAddMatrix ( stiffMat );
      }

      sysMat += stiffMat;
      sysMat *= grid.H();
      sysMat += massMat;

      aol::SSORPreconditioner < aol::Vector<RealType>, MatrixType > precond ( sysMat );
      aol::PCGInverse < aol::Vector<RealType> > solver ( sysMat, precond, 1.0e-16, 1000 );
      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

      qc::ScalarArray<RealType, qc::QC_3D> temperature ( grid ), rhs ( grid );

      aol::NoiseOperator<RealType> noiseOp;
      noiseOp.applySingle ( temperature );

      massMat.apply ( temperature, rhs );
      solver.apply ( rhs, temperature );
    }

    { // complicated domain elasticity test

      const RealType
        E       = 1.0,
        nu      = 0.33,
        lambda  = tpcfe::computeLambda ( E, nu ),
        mu      = tpcfe::computeMu ( E, nu );


      typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD > GridType;
      typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> > ConfiguratorType;

      static const aol::GaussSeidelSweepingMode gsmode = aol::GAUSS_SEIDEL_ZEBRA2;

      GridType grid ( cdLevel );
      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

      qc::ShapeLevelsetGenerator<RealType>::generateConeLevelset ( levelset, 0.231, 0.421 );

      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes();

      qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
      DirichletMask.setAll ( false );

      qc::MultiArray< RealType, qc::QC_3D > dirichletBCs ( grid ), u( grid ), rhs ( grid );

      tpcfe::setGeneralShearing ( grid, dirichletBCs, DirichletMask, qc::QC_Z, qc::QC_X, 1.0 );

      grid.setDirichletMask ( DirichletMask );
      grid.setDOFMaskFromDirichletAndDomainNodeMask();

      tpcfe::CFEElastOp< ConfiguratorType > elastop ( grid, lambda, mu );

      u -= dirichletBCs;

      elastop.restrictNonDomainEntries();

      grid.restrictToDomain ( u );

      elastop.apply ( u, rhs );

      elastop.restrictDirichletEntries();
      grid.restrictDirichletNodes ( rhs );

      tpcfe::CFEBlockMultigrid< tpcfe::CFEElastOp< ConfiguratorType >, aol::ONTHEFLY, mg::ONTHEFLY_COARSENING > mgsolver ( grid, elastop.getBlockMatrixRef(),
                                                                                                                           2, 3, 3,
                                                                                                                           gsmode, 1.0, mg::MULTIGRID_V_CYCLE,
                                                                                                                           1.0e-16, 1000 );
      mgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      mgsolver.setCoarseSolverSteps( 10000 );

      mgsolver.setVerboseMode( 1 );

      mgsolver.apply( rhs, u );
    }

    timer.stop();

    if ( bench ) {
      const RealType
        defaultRuntime = 117.319;

      const int
        nupsi = static_cast<int>( 100 * defaultRuntime / timer.elapsedCpuTime() + 0.5 ),
        wupsi = static_cast<int>( 100 * defaultRuntime / timer.elapsedWallClockTime() + 0.5 );

      aol::logBenchmarkResult ( "tpcfeBenchmark", nupsi, wupsi, resultfilename );
    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}


