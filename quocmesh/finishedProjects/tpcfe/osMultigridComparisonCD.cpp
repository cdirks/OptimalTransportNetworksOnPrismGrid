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

#include <scalarArray.h>
#include <indexMapper.h>
#include <parameterParser.h>
#include <boomerAMG.h>
#include <quocMultigrid.h>
#include <shapeLevelsetGenerator.h>

#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>

#include "affine.h"
#include "jumpCoeffMG.h"
#include "tpcfe_utils.h"
#include "osTestfunction.h"


typedef double RealType;
static const tpcfe::ConstraintType AT = tpcfe::CFE_CD;
typedef tpcfe::CFEGrid < RealType, AT > GridType;

void doEggBox ( ) {
  aol::StopWatch timer;

  const int level = 8;
  const int numSmooth = 3;

  const RealType eps = 1.0e-16;

  GridType grid ( level );
  qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

  qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeRippleLevelset ( levelset, 2.0/3.0, aol::Vec3<RealType> ( 1./15., 1./25., 1./40. ), aol::Vec3<RealType> ( 6., 13., 25. ) );

  grid.setDomainFrom ( levelset );
  grid.detectAndInitVirtualNodes();

  {
    typedef tpcfe::CFEHybridMatrix< GridType >                MatrixType;
    typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
    typedef tpcfe::CFEMassOp < ConfiguratorType >             MassOpType;
    typedef tpcfe::CFEStiffOp < ConfiguratorType >            StiffOpType;

    qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), bcond ( grid ), Lbcond ( grid );
    qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

    for ( qc::RectangularBoundaryIterator<qc::QC_3D> bnit ( grid ); bnit.notAtEnd(); ++bnit ) {
      if ( ( *bnit ) [2] == 0 ) {
        DirichletMask.set ( *bnit, true );
        bcond.set ( *bnit, 0.0 );
      }

      if ( ( *bnit ) [2] == grid.getNumZ() - 1 ) {
        DirichletMask.set ( *bnit, true );
        bcond.set ( *bnit, 1.0 );
      }
    }

    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    cerr << "Setting up stiffOp" << endl;
    StiffOpType stiffOp ( grid, aol::ONTHEFLY );

    cerr << "Assembling matrix" << endl;
    MatrixType sysMat ( grid );
    stiffOp.assembleAddMatrix ( sysMat );
    restrictNonDomainEntries ( grid, sysMat, 1.0 );

    grid.restrictToDomain ( bcond );

    cerr << "Setting up RHS" << endl;
    sysMat.apply ( bcond, Lbcond );
    rhs -= Lbcond;

    // modify matrix for boundary conditions
    tpcfe::enforceZeroDirichletBCs ( grid, sysMat, rhs );

    // 0. SSOR-PCG
    {
      qc::ScalarArray<RealType, qc::QC_3D> soln ( grid );
      timer.start();
      aol::SSORPreconditioner<aol::Vector<RealType>, MatrixType > prec ( sysMat );
      aol::PCGInverse< aol::Vector<RealType> > pcgsolver ( sysMat, prec, eps, 5000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      pcgsolver.apply ( rhs, soln );
      timer.stop();
      cerr << "TIME   SSORPCG " << aol::intFormat ( 23 ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( pcgsolver.getCount() ) << " iterations" <<  endl;
    }

    // I. CFE+AMG
#ifdef USE_EXTERNAL_HYPRE
    if ( level < 8 ) {
      qc::ScalarArray<RealType, qc::QC_3D> soln ( grid );
      timer.start();
      mg::BoomerAMGSolver< MatrixType, RealType > AMGSolver ( sysMat, eps, 5000 );
      AMGSolver.apply ( rhs, soln );

      timer.stop();
      cerr << "TIME CFE   AMG " << aol::intFormat ( 42 ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( AMGSolver.getCount() ) << " iterations" <<  endl;
    }
#endif

    // II. CFE+CFEMG
    { const short explicitLevel = 2;
      qc::ScalarArray<RealType, qc::QC_3D> soln ( grid );
      timer.start();
      const RealType relax = 1.0;

      tpcfe::CFEMultigrid< tpcfe::CFEMassOp<ConfiguratorType>, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgsolver ( grid, sysMat, explicitLevel, numSmooth, numSmooth, aol::GAUSS_SEIDEL_FORWARD, relax,
                                                                                                                      mg::MULTIGRID_V_CYCLE, eps, 5000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      mgsolver.setVerboseMode ( 6 );
      mgsolver.apply ( rhs, soln );
      timer.stop();
      cerr << "TIME CFE CFEMG " << aol::intFormat ( explicitLevel ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( mgsolver.getCount() ) << " iterations" <<  endl;
    }

  }

}


void doSlot ( ) {
  aol::StopWatch timer;

  const int level = 7;
  const int numSmooth = 3;

  const RealType eps = 1.0e-16;

  GridType grid ( level );
  qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

  qc::ShapeLevelsetGenerator<RealType>::generateSlotLevelset( levelset, 1.0/14.0, 0.8 );
  levelset.save ( "out/slotLevelset.dat.bz2", qc::PGM_FLOAT_BINARY );
  // levelset.setAll ( -1.0 );

  grid.setDomainFrom ( levelset );
  grid.detectAndInitVirtualNodes();

  {
    typedef tpcfe::CFEHybridMatrix< GridType >                MatrixType;
    typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
    typedef tpcfe::CFEMassOp < ConfiguratorType >             MassOpType;
    typedef tpcfe::CFEStiffOp < ConfiguratorType >            StiffOpType;

    qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), bcond ( grid ), Lbcond ( grid );
    qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

    for ( qc::RectangularBoundaryIterator<qc::QC_3D> bnit ( grid ); bnit.notAtEnd(); ++bnit ) {
      if ( ( *bnit ) [0] == 0 && levelset.get ( *bnit ) < 0 ) {
        DirichletMask.set ( *bnit, true );
        bcond.set ( *bnit, ( (*bnit)[2] < grid.getNumZ() / 2 ? -1 : 1 ) );
      }
    }

    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    cerr << "Setting up stiffOp" << endl;
    StiffOpType stiffOp ( grid, aol::ONTHEFLY );

    cerr << "Assembling matrix" << endl;
    MatrixType sysMat ( grid );
    stiffOp.assembleAddMatrix ( sysMat );
    restrictNonDomainEntries ( grid, sysMat, 1.0 );

    grid.restrictToDomain ( bcond );

    cerr << "Setting up RHS" << endl;
    sysMat.apply ( bcond, Lbcond );
    rhs -= Lbcond;

    // modify matrix for boundary conditions
    tpcfe::enforceZeroDirichletBCs ( grid, sysMat, rhs );

    // II. CFE standard multigrid
    for ( short explicitLevel = level - 1 ; explicitLevel >= 0; --explicitLevel ) {
      qc::ScalarArray<RealType, qc::QC_3D> soln ( grid );
      timer.start();
      const RealType relax = 1.0;

      tpcfe::CFEMultigrid< tpcfe::CFEMassOp<ConfiguratorType>, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgsolver ( grid, sysMat, explicitLevel, numSmooth, numSmooth, aol::GAUSS_SEIDEL_FORWARD, relax,
                                                                                                                      mg::MULTIGRID_V_CYCLE, eps, 5000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      mgsolver.setVerboseMode ( 2 );
      mgsolver.apply ( rhs, soln );
      timer.stop();
      cerr << "TIME CFE CFEMG " << aol::intFormat ( explicitLevel ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( mgsolver.getCount() ) << " iterations" <<  endl;

      soln += bcond;

      soln.save ( "out/slotSolution.dat.bz2", qc::PGM_FLOAT_BINARY );

    }

  }
}


int main ( int, char** ) {
  try {

    doEggBox();
    // doSlot();

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
