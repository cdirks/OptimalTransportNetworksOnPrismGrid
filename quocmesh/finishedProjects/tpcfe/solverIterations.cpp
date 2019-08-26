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

#include <boomerAMG.h>
#include <preconditioner.h>
#include <solver.h>
#include <shapeLevelsetGenerator.h>

#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>


typedef double RealType;

static const RealType eps = 1.0e-16;


void pcgRun ( const aol::Vector<RealType> &rhs, const aol::Vector<RealType> &soln, const aol::Op< aol::Vector<RealType> > &sysMat, const aol::Op< aol::Vector<RealType> > &prec, const char* const slText, aol::StopWatch &timer ) {
  aol::Vector<RealType> pcgsoln ( rhs, aol::STRUCT_COPY );
  aol::PCGInverse< aol::Vector<RealType> > pcgsolver ( sysMat, prec, eps, 1000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
  // pcgsolver.setQuietMode ( true );
  pcgsolver.apply ( rhs, pcgsoln );
  timer.stop();
  pcgsoln -= soln;
  cerr << slText << aol::detailedFormat ( timer.elapsedCpuTime() ) << " seconds for " << aol::intFormat ( pcgsolver.getCount() ) << " iterations, difference " << aol::detailedFormat ( pcgsoln.norm() / pcgsoln.size() ) <<  endl;
}


void pcgRun ( const aol::MultiVector<RealType> &rhs, const aol::MultiVector<RealType> &soln, const aol::Op< aol::MultiVector<RealType> > &sysMat, const aol::Op< aol::MultiVector<RealType> > &prec, const char* const slText, aol::StopWatch &timer ) {
  aol::MultiVector<RealType> pcgsoln ( rhs, aol::STRUCT_COPY );
  aol::PCGInverse< aol::MultiVector<RealType> > pcgsolver ( sysMat, prec, eps, 10000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
  // pcgsolver.setQuietMode ( true );
  pcgsolver.apply ( rhs, pcgsoln );
  timer.stop();
  pcgsoln -= soln;
  cerr << slText << aol::detailedFormat ( timer.elapsedCpuTime() ) << " seconds for " << aol::intFormat ( pcgsolver.getCount() ) << " iterations, difference " << aol::detailedFormat ( pcgsoln.norm() / pcgsoln.getTotalSize() ) <<  endl;
}


template< typename ConfiguratorType, typename GridType,  typename MatrixType >
void solverTests ( const GridType& grid, const MatrixType &sysMat ) {
  aol::StopWatch timer;
  aol::Vector<RealType> soln ( grid ), rhs ( grid );
  aol::NoiseOperator<RealType> nO (aol::NoiseOperator<RealType>::EQUALLY_DISTRIBUTED ); // nO.randomize();
  nO.applySingle ( soln );
  sysMat.apply ( soln, rhs );

  {
    aol::Vector<RealType> cgsoln ( grid );
    timer.start();
    aol::CGInverse< aol::Vector<RealType> > cgsolver ( sysMat, eps, 10000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    cgsolver.apply ( rhs, cgsoln );
    timer.stop();
    cgsoln -= soln;
    cerr << "CG solver took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << " seconds for " << aol::intFormat ( cgsolver.getCount() ) << " iterations, difference " << aol::detailedFormat ( cgsoln.norm() / cgsoln.size() ) <<  endl;
  }
  {
    timer.start();
    aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( sysMat );
    timer.stop();
    cerr << "Diagonal Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, sysMat, prec, "Diagonal-PCG solver took ", timer );
  }
  {
    timer.start();
    aol::GeometricScalingPreconditioner< aol::Vector<RealType> > prec ( sysMat );
    timer.stop();
    cerr << "GeometricScaling Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, sysMat, prec, "GeomScaling-PCG solver took ", timer );
  }
  {
    timer.start();
    aol::SSORPreconditioner<aol::Vector<RealType>, MatrixType > prec ( sysMat );
    timer.stop();
    cerr << "SSOR Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, sysMat, prec, "SSOR-PCG solver took ", timer );
  }
  {
    timer.start();
    aol::SparseMatrix< RealType > sysMatSparse ( sysMat.getNumRows(), sysMat.getNumCols() );
    sysMatSparse += sysMat;
    aol::ILU0Preconditioner< RealType, aol::SparseMatrix<RealType> > prec ( sysMatSparse );
    timer.stop();
    cerr << "ILU0(*) Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, sysMatSparse, prec, "ILU0-PCG(Sparse) solver took ", timer );
  }
  {
    timer.start();
    tpcfe::SlopeInterface< tpcfe::CFEGrid< RealType, tpcfe::CFE_TPOS > >::_useIncorrectStandardCoarsening = false;
    tpcfe::CFEMultigridPreconditioner< tpcfe::CFEStiffOp < ConfiguratorType >, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgPrecond ( grid, sysMat, 2, 3, 3 );
    timer.stop();
    cerr << "MG Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, sysMat, mgPrecond, "MG3V(2)3-PCG solver took ", timer );
  }
  {
    aol::Vector<RealType> mgsoln ( grid );
    timer.start();
    tpcfe::SlopeInterface< tpcfe::CFEGrid< RealType, tpcfe::CFE_TPOS > >::_useIncorrectStandardCoarsening = false;
    tpcfe::CFEMultigrid< tpcfe::CFEStiffOp < ConfiguratorType >, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgsolver ( grid, sysMat, 2, 3, 3,
                                                                                                                        aol::GAUSS_SEIDEL_FORWARD, 1.0,
                                                                                                                        mg::MULTIGRID_V_CYCLE, eps, 1000,
                                                                                                                        aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    mgsolver.setVerboseMode ( 1 );
    timer.stop();
    cerr << "MG Solver setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    mgsolver.apply ( rhs, mgsoln );
    timer.stop();
    mgsoln -= soln;
    cerr << "CFE-MG solver took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << " seconds for " << aol::intFormat ( mgsolver.getCount() ) << " iterations, difference " << aol::detailedFormat ( mgsoln.norm() / mgsoln.size() ) <<  endl;
  }

  {
    timer.start();
    tpcfe::SlopeInterface< tpcfe::CFEGrid< RealType, tpcfe::CFE_TPOS > >::_useIncorrectStandardCoarsening = true;
    tpcfe::CFEMultigridPreconditioner< tpcfe::CFEStiffOp < ConfiguratorType >, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgPrecond ( grid, sysMat, 2, 3, 3 );
    timer.stop();
    cerr << "MG Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, sysMat, mgPrecond, "MG3V(2)3-PCG solver took ", timer );
  }
  {
    aol::Vector<RealType> mgsoln ( grid );
    timer.start();
    tpcfe::SlopeInterface< tpcfe::CFEGrid< RealType, tpcfe::CFE_TPOS > >::_useIncorrectStandardCoarsening = true;
    tpcfe::CFEMultigrid< tpcfe::CFEStiffOp < ConfiguratorType >, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgsolver ( grid, sysMat, 2, 3, 3,
                                                                                                                        aol::GAUSS_SEIDEL_FORWARD, 1.0,
                                                                                                                        mg::MULTIGRID_V_CYCLE, eps, 1000,
                                                                                                                        aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    mgsolver.setVerboseMode ( 1 );
    timer.stop();
    cerr << "MG Solver setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    mgsolver.apply ( rhs, mgsoln );
    timer.stop();
    mgsoln -= soln;
    cerr << "CFE-MG solver took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << " seconds for " << aol::intFormat ( mgsolver.getCount() ) << " iterations, difference " << aol::detailedFormat ( mgsoln.norm() / mgsoln.size() ) <<  endl;
  }

#ifdef USE_EXTERNAL_HYPRE
  {
    aol::Vector<RealType> amgsoln ( grid );
    timer.start();
    mg::BoomerAMGSolver< MatrixType, RealType > amgsolver ( sysMat, eps, 5000 );
    timer.stop();
    cerr << "AMG Solver setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    amgsolver.apply ( rhs, amgsoln );
    timer.stop();
    amgsoln -= soln;
    cerr << "AMG solver took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << " seconds for " << aol::intFormat ( amgsolver.getCount() ) << " iterations, difference " << aol::detailedFormat ( amgsoln.norm() / amgsoln.size() ) <<  endl;
  }
#endif
}


template< typename ConfiguratorType, typename GridType, typename ElastOpType >
void solverTestsElast ( const GridType& grid, const ElastOpType &elastOp ) {
  aol::StopWatch timer;
  aol::MultiVector<RealType> soln ( grid ), rhs ( grid );
  aol::NoiseOperator<RealType> nO (aol::NoiseOperator<RealType>::EQUALLY_DISTRIBUTED ); // nO.randomize();
  nO.applySingle ( soln );
  grid.restrictToDofs ( soln );
  elastOp.apply ( soln, rhs );

  {
    aol::MultiVector<RealType> cgsoln ( grid );
    timer.start();
    aol::CGInverse< aol::MultiVector<RealType> > cgsolver ( elastOp, eps, 10000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    cgsolver.apply ( rhs, cgsoln );
    timer.stop();
    cgsoln -= soln;
    cerr << "CG solver took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << " seconds for " << aol::intFormat ( cgsolver.getCount() ) << " iterations, difference " << aol::detailedFormat ( cgsoln.norm() / cgsoln.getTotalSize() ) <<  endl;
  }

  {
    timer.start();
    aol::DiagonalPreconditioner< aol::MultiVector<RealType> > prec ( elastOp.getBlockMatrixRef() );
    timer.stop();
    cerr << "Diagonal Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, elastOp.getBlockMatrixRef(), prec, "Diagonal-PCG solver took ", timer );
  }

  {
    timer.start();
    aol::BlockDiagonalPreconditioner< RealType, qc::QC_3D > prec ( elastOp.getBlockMatrixRef() );
    timer.stop();
    cerr << "BlockDiagonal Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, elastOp.getBlockMatrixRef(), prec, "BlockDiagonal-PCG solver took ", timer );
  }

  {
    timer.start();
    aol::GeometricScalingPreconditioner< aol::MultiVector<RealType> > prec ( elastOp.getBlockMatrixRef() );
    timer.stop();
    cerr << "GeometricScaling Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, elastOp.getBlockMatrixRef(), prec, "GeomScaling-PCG solver took ", timer );
  }

  {
    timer.start();
    aol::SSORPreconditioner<aol::MultiVector<RealType>, typename ElastOpType::MatrixType > prec ( elastOp.getBlockMatrixRef() );
    timer.stop();
    cerr << "SSOR Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, elastOp.getBlockMatrixRef(), prec, "SSOR-PCG solver took ", timer );
  }

  {
    timer.start();
    aol::BlockSSORPreconditioner< RealType, typename ElastOpType::BlockMatrixType, qc::QC_3D > prec ( elastOp.getBlockMatrixRef() );
    timer.stop();
    cerr << "BlockSSOR Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, elastOp.getBlockMatrixRef(), prec, "BlockSSOR-PCG solver took ", timer );
  }

  {
    timer.start();
    aol::SparseMatrix< RealType > sysMatSparse ( 3 * grid.getNumberOfNodes(), 3 * grid.getNumberOfNodes() );
    elastOp.getBlockMatrixRef().addUnblockedMatrixTo ( sysMatSparse );
    aol::ILU0Preconditioner< RealType, aol::SparseMatrix<RealType> > prec ( sysMatSparse );
    aol::Vector<RealType> rhsUnblocked ( rhs ), solnUnblocked ( soln );
    timer.stop();
    cerr << "ILU0(*) Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhsUnblocked, solnUnblocked, sysMatSparse, prec, "ILU0-PCG(Sparse) solver took ", timer );
  }


  {
    timer.start();
    tpcfe::SlopeInterface< tpcfe::CFEGrid< RealType, tpcfe::CFE_TPOS > >::_useIncorrectStandardCoarsening = false;
    tpcfe::CFEBlockMultigridPreconditioner< tpcfe::CFEElastOp < ConfiguratorType >, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgPrecond ( grid, elastOp.getBlockMatrixRef(), 2, 3, 3 );
    timer.stop();
    cerr << "MG Preconditioner setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    pcgRun ( rhs, soln, elastOp.getBlockMatrixRef(), mgPrecond, "MG3V(2)3-PCG solver took ", timer );
  }

  {
    aol::MultiVector<RealType> mgsoln ( grid );
    timer.start();
    tpcfe::SlopeInterface< tpcfe::CFEGrid< RealType, tpcfe::CFE_TPOS > >::_useIncorrectStandardCoarsening = false;
    tpcfe::CFEBlockMultigrid< tpcfe::CFEElastOp < ConfiguratorType >, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgsolver ( grid, elastOp.getBlockMatrixRef(), 2, 3, 3,
                                                                                                                             aol::GAUSS_SEIDEL_FORWARD, 1.0,
                                                                                                                             mg::MULTIGRID_V_CYCLE, eps, 1000,
                                                                                                                             aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    mgsolver.setVerboseMode ( 1 );
    timer.stop();
    cerr << "MG Solver setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    mgsolver.apply ( rhs, mgsoln );
    timer.stop();
    mgsoln -= soln;
    cerr << "CFE-BlockMG solver took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << " seconds for " << aol::intFormat ( mgsolver.getCount() ) << " iterations, difference " << aol::detailedFormat ( mgsoln.norm() / mgsoln.getTotalSize() ) <<  endl;
  }

#ifdef USE_EXTERNAL_HYPRE
  {
    aol::SparseMatrix< RealType > sysMatSparse ( 3 * grid.getNumberOfNodes(), 3 * grid.getNumberOfNodes() );
    elastOp.getBlockMatrixRef().addUnblockedMatrixTo ( sysMatSparse );

    aol::Vector<RealType> amgsoln ( 3 * grid.getNumberOfNodes() );
    timer.start();
    mg::BoomerAMGSolver< aol::SparseMatrix<RealType>, RealType > amgsolver ( sysMatSparse, eps, 5000 );
    aol::Vector<RealType> rhsUnblocked ( rhs ), solnUnblocked ( soln );
    timer.stop();
    cerr << "AMG Solver setup took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << endl;
    timer.start();
    amgsolver.apply ( rhsUnblocked, amgsoln );
    timer.stop();
    amgsoln -= solnUnblocked;
    cerr << "AMG solver took " << aol::detailedFormat ( timer.elapsedCpuTime() ) << " seconds for " << aol::intFormat ( amgsolver.getCount() ) << " iterations, difference " << aol::detailedFormat ( amgsoln.norm() / amgsoln.size() ) <<  endl;
  }
#endif
}




int main ( int, char** ) {
  try {

    const int level = 6;

    {
      cerr << "Testing CFE_CD" << endl;
      const tpcfe::ConstraintType CT = tpcfe::CFE_CD;

      typedef tpcfe::CFEGrid< RealType, CT >                                         GridType;
      typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> >  ConfiguratorType;
      typedef ConfiguratorType::MatrixType                                           MatrixType;

      typedef tpcfe::CFEMassOp  < ConfiguratorType >                                 MassOpType;
      typedef tpcfe::CFEStiffOp  < ConfiguratorType >                                StiffOpType;

      cerr << "--- CFE_CD scalar solver tests" << endl;
      {
        GridType grid ( level );
        qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
        // levelset.load ( "datasets/961_T1Po2_cubeA_L6_070mu.dat.bz2" );
        qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );
        // qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeRippleLevelset ( levelset, 2.0/3.0, aol::Vec3<RealType> ( 1./40., 1./60., 1./110. ), aol::Vec3<RealType> ( 4., 7., 19. ) );

        grid.setDomainFrom ( levelset );
        grid.detectAndInitVirtualNodes();
        grid.setDOFMaskFromDirichletAndDomainNodeMask();

        MatrixType sysMat ( grid );
        {
          StiffOpType StiffOp ( grid, aol::ONTHEFLY );
          StiffOp._quietMode = true;
          MassOpType MassOp ( grid, aol::ONTHEFLY );
          MassOp._quietMode = true;
          StiffOp.assembleAddMatrix ( sysMat );
          sysMat *= grid.H();
          MassOp.assembleAddMatrix ( sysMat );
        }
        restrictNonDomainEntries ( grid,  sysMat, 1.0 );

        solverTests< ConfiguratorType, GridType, MatrixType > ( grid, sysMat );
      }

      {
        GridType grid ( level-1 );
        qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
        // levelset.load ( "datasets/961_T1Po2_cubeA_L5_070mu.dat.bz2" );
        qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );

        grid.setDomainFrom ( levelset );
        grid.detectAndInitVirtualNodes();
        grid.setDOFMaskFromDirichletAndDomainNodeMask();

        cerr << "--- CFE_CD elasticity solver tests" << endl;
        qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
        DirichletMask.setAll ( false );

        for ( int i = 0; i < grid.getWidth(); ++i ) {
          for ( int j = 0; j < grid.getWidth(); ++j ) {
            DirichletMask.set ( i, j, 0, true );
            DirichletMask.set ( i, j, grid.getWidth() - 1, true );
          }
        }

        grid.setDirichletMask ( DirichletMask );
        grid.setDOFMaskFromDirichletAndDomainNodeMask();

        typedef tpcfe::CFEElastOp< ConfiguratorType > ElastOpType;
        ElastOpType ElastOp ( grid, tpcfe::computeLambda ( 1.0, 0.33 ), tpcfe::computeMu ( 1.0, 0.33 ) );

        ElastOp.restrictNonDomainEntries();
        ElastOp.restrictDirichletEntries();

        solverTestsElast< ConfiguratorType, GridType, ElastOpType > ( grid, ElastOp );
      }
    }

    {
      cerr << "Testing CFE_TPOS for kappa = 42" << endl;

      const tpcfe::ConstraintType CT = tpcfe::CFE_TPOS;

      typedef tpcfe::CFEGrid< RealType, CT >                                         GridType;
      typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> >  ConfiguratorType;
      typedef ConfiguratorType::MatrixType                                           MatrixType;

      typedef tpcfe::CFEMassOp  < ConfiguratorType >                                 MassOpType;
      typedef tpcfe::CFEStiffOp  < ConfiguratorType >                                StiffOpType;

      GridType grid ( level );
      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
      // levelset.load ( "datasets/961_T1Po2_cubeA_L6_070mu.dat.bz2" );
      qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );

      grid.addStructureFrom ( levelset );
      qc::AArray< RealType, qc::QC_3D > coeff ( grid );
      tpcfe::setCoeffForLevelset ( coeff, levelset, 42.0, 1.0 );
      grid.relaxedDetectAndInitVirtualNodes( coeff, 1.0e-14, 1.0e-14 );
      grid.setDOFMaskFromDirichletAndDomainNodeMask();

      MatrixType sysMat ( grid );
      {
        StiffOpType StiffOp ( grid, aol::ONTHEFLY );
        StiffOp._quietMode = true;
        MassOpType MassOp ( grid, aol::ONTHEFLY );
        MassOp._quietMode = true;
        StiffOp.assembleAddMatrix ( sysMat );
        sysMat *= grid.H();
        MassOp.assembleAddMatrix ( sysMat );
      }

      solverTests< ConfiguratorType, GridType, MatrixType > ( grid, sysMat );
    }

    {
      cerr << "Testing CFE_TPOSELAST for ( 70.0, 0.35 ); ( 3.0, 0.38 )" << endl;
      const tpcfe::ConstraintType CT = tpcfe::CFE_TPOSELAST;

      typedef tpcfe::CFEGrid < RealType, CT, tpcfe::IsotropicElasticityCoefficient<RealType> >  GridType;
      typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> >             ConfiguratorType;
      typedef ConfiguratorType::MatrixType                                                      MatrixType;

      GridType grid ( level - 1 );
      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
      // levelset.load ( "datasets/961_T1Po2_cubeA_L5_070mu.dat.bz2" );
      qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );

      grid.addStructureFrom ( levelset );
      qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );
      const GridType::NodalCoeffType ENuMinus ( 70.0, 0.35 ), ENuPlus ( 3.0, 0.38 );
      tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
      grid.relaxedDetectAndInitVirtualNodes ( coeff, 5.0e-15, 1.0e-15 );

      qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
      DirichletMask.setAll ( false );

      for ( int i = 0; i < grid.getWidth(); ++i ) {
        for ( int j = 0; j < grid.getWidth(); ++j ) {
          DirichletMask.set ( i, j, 0, true );
          DirichletMask.set ( i, j, grid.getWidth() - 1, true );
        }
      }
      grid.setDirichletMask ( DirichletMask );
      grid.setDOFMaskFromDirichletAndDomainNodeMask();

      typedef tpcfe::CFEJCElastOp< ConfiguratorType > ElastOpType;
      ElastOpType ElastOp ( grid, coeff );

      solverTestsElast< ConfiguratorType, GridType, ElastOpType > ( grid, ElastOp );

    } // end CFE_TPOSELAST

  } catch ( aol::Exception &exc ) {
    exc.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
