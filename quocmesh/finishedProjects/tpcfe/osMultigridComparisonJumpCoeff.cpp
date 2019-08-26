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
#include <shapeLevelsetGenerator.h>

#include <boomerAMG.h>
#include <quocMultigrid.h>

#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>


#include "affine.h"
#include "jumpCoeffMG.h"
#include "tpcfe_utils.h"
#include "osTestfunction.h"
#include "levelsets.h"

typedef double RealType;

static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;
typedef tpcfe::CFEGrid < RealType, AT > GridType;


void doEggBox ( ) {
  //   tpcfe::SlopeInterface< GridType >::_performNAKRelaxationAtCircumferentialEdges = false;  // default: false
  //   tpcfe::SlopeInterface< GridType >::_performNAKRelaxationAtAllRemainingEdges = true;      // default: true

  const char* levFN = "out/eggBox.dat.bz2";
  const char* nakQN = NULL; // "out/nakQN_RH1_%02d.out";
  const char* nakQS = NULL; // "out/nakQS_RH1_%02d.out";


  aol::StopWatch timer;

  const RealType DC_PLUS = 1.0, DC_MINUS = 10.0;

  const int level = 8;
  const int numSmooth = 3;

  // const int numRods = 3; // or holes in the cheese ...

  const RealType eps = 1.0e-16;

  GridType grid ( level );
  qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

  // qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeRippleLevelset ( levelset, 2.0/3.0, aol::Vec3<RealType> ( 1./60., 1./90., 1./120. ), aol::Vec3<RealType> ( 14., 19., 25. ) );
  qc::ShapeLevelsetGenerator<RealType>::generateFourPlaneRippleLevelset ( levelset,
                                                                          aol::Vec3<RealType> ( 1./40, 1./80, 1./160 ), aol::Vec3<RealType> ( 1./40, 1./80, 1./160 ), aol::Vec3<RealType> ( 1./40, 1./60, 1./80 ), aol::Vec3<RealType> ( 1./40, 1./80, 1./160 ),
                                                                          aol::Vec3<RealType> ( 2., 4., 8. ),           aol::Vec3<RealType> ( 2.5, 4.5, 8.5 ),        aol::Vec3<RealType> ( 8., 4., 2. ),          aol::Vec3<RealType> ( 3., 5., 9. ) );

  levelset.save ( levFN, qc::PGM_FLOAT_BINARY );

  grid.addStructureFrom ( levelset );

  qc::AArray<RealType, qc::QC_3D> coeff ( grid );
  for ( int i = 0; i < levelset.size(); ++i )
    coeff[i] = ( levelset[i] < 0 ? DC_MINUS : DC_PLUS );

  grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-15, 1.0e-15 );

  {
    typedef tpcfe::CFEHybridMatrix< GridType >                MatrixType;
    typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
    typedef tpcfe::CFEMassOp < ConfiguratorType >             MassOpType;
    typedef tpcfe::CFEStiffOpWI < ConfiguratorType >          WStiffOpType;

    qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), bcond ( grid ), Lbcond ( grid );
    qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

    for ( qc::RectangularBoundaryIterator<qc::QC_3D> bnit ( grid ); bnit.notAtEnd(); ++bnit ) {
      if ( ( *bnit ) [1] == 0 ) {
        DirichletMask.set ( *bnit, true );
        bcond.set ( *bnit, 0.0 );
      }

      if ( ( *bnit ) [1] == grid.getNumY() - 1 ) {
        DirichletMask.set ( *bnit, true );
        bcond.set ( *bnit, 1.0 );
      }
    }

    grid.setDirichletMask ( DirichletMask );


    cerr << "Setting up stiffOp" << endl;
    WStiffOpType stiffOp ( coeff, grid, aol::ONTHEFLY );

    cerr << "Assembling matrix" << endl;
    MatrixType sysMat ( grid );
    stiffOp.assembleAddMatrix ( sysMat );

    cerr << "Setting up RHS" << endl;
    sysMat.apply ( bcond, Lbcond );
    rhs -= Lbcond;

    // modify matrix for boundary conditions
    tpcfe::enforceZeroDirichletBCs ( grid, sysMat, rhs );

#if 1
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
#endif
#ifdef USE_EXTERNAL_HYPRE
    // I. AMG
    if ( level < 8 ) {
      qc::ScalarArray<RealType, qc::QC_3D> soln ( grid );
      timer.start();
      mg::BoomerAMGSolver< MatrixType, RealType > AMGSolver ( sysMat, eps, 5000 );
      AMGSolver.apply ( rhs, soln );

      // comparing solutions would now require adding the boundary conditions ...

      timer.stop();
      cerr << "TIME CFE   AMG " << aol::intFormat ( 42 ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( AMGSolver.getCount() ) << " iterations" <<  endl;
    }
#endif
#if 1
    // II. CFE+CFEMG
    // for ( short explicitLevel = level; explicitLevel >= 3; --explicitLevel ) {
    { const short explicitLevel = 2;
      qc::ScalarArray<RealType, qc::QC_3D> soln ( grid );
      timer.start();
      const RealType relax = 1.0;

      tpcfe::SlopeInterface< GridType >::_useIncorrectStandardCoarsening = false;
      tpcfe::SlopeInterface< GridType >::_checkNAKQualityAfterCorrection = nakQN;
      tpcfe::CFEMultigrid< tpcfe::CFEMassOp<ConfiguratorType>, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgsolver ( grid, sysMat, explicitLevel, numSmooth, numSmooth, aol::GAUSS_SEIDEL_FORWARD, relax,
                                                                                                                      mg::MULTIGRID_V_CYCLE, eps, 5000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      mgsolver.setVerboseMode ( 6 );
      mgsolver.apply ( rhs, soln );
      timer.stop();
      cerr << "TIME CFE CFEMG " << aol::intFormat ( explicitLevel ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( mgsolver.getCount() ) << " iterations" <<  endl;
    }
#endif
#if 1
    // III. CFE+StdMG
    // for ( short explicitLevel = level; explicitLevel >= 3; --explicitLevel ) {
    { const short explicitLevel = 2;
      qc::ScalarArray<RealType, qc::QC_3D> soln ( grid );
      timer.start();
      const RealType relax = 1.0;

      tpcfe::SlopeInterface< GridType >::_useIncorrectStandardCoarsening = true;
      tpcfe::SlopeInterface< GridType >::_checkNAKQualityAfterCorrection = nakQS;

      tpcfe::CFEMultigrid< tpcfe::CFEMassOp<ConfiguratorType>, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgsolver ( grid, sysMat, explicitLevel, numSmooth, numSmooth, aol::GAUSS_SEIDEL_FORWARD, relax,
                                                                                                                      mg::MULTIGRID_V_CYCLE, eps, 5000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      mgsolver.setVerboseMode ( 6 );
      mgsolver.apply ( rhs, soln );
      timer.stop();
      cerr << "TIME CFE StdMG " << aol::intFormat ( explicitLevel ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( mgsolver.getCount() ) << " iterations" <<  endl;
    }
#endif
#if 0
    // I. CFE-MG as preconditioner in CG
    {
      qc::ScalarArray<RealType, qc::QC_3D> soln ( grid );
      timer.start();

      tpcfe::CFEMultigrid< tpcfe::CFEMassOp<ConfiguratorType>, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > prec ( grid, sysMat, 0, numSmooth, numSmooth, aol::GAUSS_SEIDEL_FORWARD, 1.0 /*relax*/,
                                                                                                                  mg::MULTIGRID_V_CYCLE, 0.999999, 1 /*step*/, aol::STOPPING_ABSOLUTE );
      prec.setVerboseMode ( 0 ); // shut up!
      aol::PCGInverse< aol::Vector<RealType> > PCGSolver ( sysMat, prec, eps, 1000 );
      PCGSolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      PCGSolver.apply ( rhs, soln );

      timer.stop();
      cerr << "TIME CFE-PCG   " << aol::intFormat ( 42 ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( PCGSolver.getCount() ) << " iterations" <<  endl;
    }

#endif
  }

}

void doSlot ( ) {
  // no spectacular results
  aol::StopWatch timer;

  const RealType DC_PLUS = 1.0, DC_MINUS = 42.0;

  const int level = 7;
  const int numSmooth = 3;

  const RealType eps = 1.0e-16;

  GridType grid ( level );
  qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

  qc::ShapeLevelsetGenerator<RealType>::generateSlotLevelset( levelset, 1.0/14.0, 0.8 );

  grid.addStructureFrom ( levelset );

  qc::AArray<RealType, qc::QC_3D> coeff ( grid );
  for ( int i = 0; i < levelset.size(); ++i )
    coeff[i] = ( levelset[i] < 0 ? DC_MINUS : DC_PLUS );

  grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-15, 1.0e-15 );

  typedef tpcfe::CFEHybridMatrix< GridType >                MatrixType;
  typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
  typedef tpcfe::CFEMassOp < ConfiguratorType >             MassOpType;
  typedef tpcfe::CFEStiffOpWI < ConfiguratorType >          WStiffOpType;

  qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), bcond ( grid ), Lbcond ( grid );
  qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

  for ( qc::RectangularBoundaryIterator<qc::QC_3D> bnit ( grid ); bnit.notAtEnd(); ++bnit ) {
    if ( ( *bnit ) [0] == 0 ) {
      if ( levelset.get ( *bnit ) < 0 ) {
        DirichletMask.set ( *bnit, true );
        bcond.set ( *bnit, ( (*bnit)[2] < grid.getNumZ() / 2 ? -1 : 1 ) );
      } else {
        // do nothing
      }
    }
  }

  grid.setDirichletMask ( DirichletMask );


  cerr << "Setting up stiffOp" << endl;
  WStiffOpType stiffOp ( coeff, grid, aol::ONTHEFLY );

  cerr << "Assembling matrix" << endl;
  MatrixType sysMat ( grid );
  stiffOp.assembleAddMatrix ( sysMat );

  cerr << "Setting up RHS" << endl;
  sysMat.apply ( bcond, Lbcond );
  rhs -= Lbcond;

  // modify matrix for boundary conditions
  tpcfe::enforceZeroDirichletBCs ( grid, sysMat, rhs );

  // std CFE MG
  for ( short explicitLevel = level - 1; explicitLevel >= 0; --explicitLevel ) {
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
