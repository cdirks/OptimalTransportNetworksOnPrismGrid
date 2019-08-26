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

#include <indexMapper.h>
#include <parameterParser.h>
#include <scalarArray.h>

#include <boomerAMG.h>
#include <quocMultigrid.h>

#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>


#include "tpcfe_utils.h"

typedef double RealType;

static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;
typedef tpcfe::CFEGrid < RealType, AT >                        CFEGridType;
typedef tpcfe::CFEHybridMatrix<CFEGridType>                    CFEMatrixType;

typedef tpcfe::CFEConfigurator < CFEGridType, CFEMatrixType >  CFEConfiguratorType;
typedef tpcfe::CFEStiffOpWI < CFEConfiguratorType >            CFEWStiffOpType;
typedef tpcfe::CFEMassOpWI < CFEConfiguratorType >             CFEWMassOpType;

int main ( int argc, char** argv ) {
  try {

    if ( argc != 2 ) {
      cerr << "usage: osMultigridComparisonAlPmmADSHC <parameterfile>" << endl;
      return( EXIT_FAILURE );
    }

    aol::ParameterParser params ( argv[1] );
    params.dump();

    // "numerics parameters"
    const int
      level             = params.getInt ( "level" );

    const int
      numSmooth         = params.getInt ( "numSmooth" );


    const RealType tauFast = params.getDouble ( "tauFast" );

    CFEGridType grid ( level );

    char levelsetFilename[1024];
    params.getString ( "levelsetFilename", levelsetFilename );
    qc::ScalarArray<RealType, qc::QC_3D> dataset ( levelsetFilename );

    grid.addStructureFrom ( dataset );

    const RealType
      sigma             = params.getDouble ( "sigma" ),         // size of bounding box; in m
      lambda_Al         = params.getDouble ( "lambda_Al" ),     // corresponding to negative level set values
      lambda_PmmA       = params.getDouble ( "lambda_PmmA" ),
      rhoC_Al           = params.getDouble ( "rhoC_Al" ),
      rhoC_PmmA         = params.getDouble ( "rhoC_PmmA" );


    // set up weightInterface
    qc::AArray<RealType, qc::QC_3D> coeffLambda ( grid );
    for ( int i = 0; i < dataset.size(); ++i )
      coeffLambda[i] = ( dataset[i] < 0 ? lambda_Al : lambda_PmmA );

    grid.relaxedDetectAndInitVirtualNodes ( coeffLambda, 1.0e-15, 5.0e-16 );

    qc::AArray<RealType, qc::QC_3D> coeffRhoC ( grid );
    for ( int i = 0; i < dataset.size(); ++i )
      coeffRhoC[i] = ( dataset[i] < 0 ? rhoC_Al : rhoC_PmmA );


    qc::ScalarArray<RealType, qc::QC_3D> temperatureOld ( grid ), temperatureBC ( grid ), rhs ( grid ), rhsModified ( grid ), dwarf ( grid );
    qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
    DirichletMask.setAll ( false );

    for ( int i = 0; i < grid.getNumX(); ++i ) {
      for ( int j = 0; j < grid.getNumX(); ++j ) {
        temperatureBC.set ( i, j,              0    , params.getDouble ( "DirichletValueBottom" ) );
        temperatureBC.set ( i, j, grid.getNumZ() - 1, params.getDouble ( "DirichletValueTop"    ) );

        DirichletMask.set ( i, j,              0    , true );
        DirichletMask.set ( i, j, grid.getNumZ() - 1, true );
      }
    }

    grid.setDirichletMask ( DirichletMask );


    CFEMatrixType massMat ( grid ), stiffMat ( grid ), sysMat ( grid ), sysMatModified ( grid );
    {
      CFEWMassOpType massOp ( coeffRhoC, grid, aol::ASSEMBLED );
      CFEWStiffOpType stiffOp ( coeffLambda, grid, aol::ASSEMBLED );
      massOp.assembleAddMatrix ( massMat );
      stiffOp.assembleAddMatrix ( stiffMat );
    }

    massMat *= aol::Cub ( sigma ); // correct scaling of units
    stiffMat *= sigma;

    temperatureOld.setAll ( params.getDouble ( "InitialValue" ) );

    sysMat += stiffMat; // M + tau L
    sysMat *= tauFast;
    sysMat += massMat;

    sysMatModified += sysMat;
    tpcfe::enforceZeroDirichletBCs ( grid, sysMatModified );

    massMat.apply ( temperatureOld, rhs );
    sysMat.apply ( temperatureBC, dwarf );
    rhs -= dwarf;

    rhsModified = rhs;
    tpcfe::enforceZeroDirichletBCs ( grid, rhsModified );


    aol::StopWatch timer;
    // I. CFE+AMG
#ifdef USE_EXTERNAL_HYPRE
    {
      qc::ScalarArray<RealType, qc::QC_3D> temperatureNew ( grid );
      timer.start();
      mg::BoomerAMGSolver< CFEMatrixType, RealType > AMGSolver ( sysMatModified, 1.0e-16, 1000 );
      AMGSolver.apply ( rhsModified, temperatureNew );
      timer.stop();
      cerr << "TIME CFE   AMG " << aol::intFormat ( 42 ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( AMGSolver.getCount() ) << " iterations" <<  endl;
    }
#endif

    // II. CFE+CFEMG

    for ( short explicitLevel = level; explicitLevel >= params.getInt ( "explicitLevel" ); --explicitLevel ) {
      qc::ScalarArray<RealType, qc::QC_3D> temperatureNew ( grid );
      timer.start();

      tpcfe::SlopeInterface< CFEGridType >::_useIncorrectStandardCoarsening = false;

      tpcfe::CFEMultigrid< /*misuse: */ CFEWStiffOpType, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > solver ( grid, sysMatModified, explicitLevel, numSmooth, numSmooth,
                                                                                                              aol::GAUSS_SEIDEL_FORWARD, 1.0,
                                                                                                              mg::MULTIGRID_V_CYCLE, 1.0e-16, 1000,
                                                                                                              aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver.setVerboseMode( 6 );
      solver.apply ( rhsModified, temperatureNew );
      timer.stop();
      cerr << "TIME CFE CFEMG " << aol::intFormat ( explicitLevel ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( solver.getCount() ) << " iterations" <<  endl;
    }


    // III. CFE+StdMG

    for ( short explicitLevel = level; explicitLevel >= params.getInt ( "explicitLevel" ); --explicitLevel ) {
      qc::ScalarArray<RealType, qc::QC_3D> temperatureNew ( grid );
      timer.start();

      tpcfe::SlopeInterface< CFEGridType >::_useIncorrectStandardCoarsening = true;

      tpcfe::CFEMultigrid< /*misuse: */ CFEWStiffOpType, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > solver ( grid, sysMatModified, explicitLevel, numSmooth, numSmooth,
                                                                                                              aol::GAUSS_SEIDEL_FORWARD, 1.0,
                                                                                                              mg::MULTIGRID_V_CYCLE, 1.0e-16, 1000,
                                                                                                              aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver.setVerboseMode( 6 );
      solver.apply ( rhsModified, temperatureNew );
      timer.stop();
      cerr << "TIME CFE StdMG " << aol::intFormat ( explicitLevel ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( solver.getCount() ) << " iterations" <<  endl;
    }

    // missing: MFE calculation!

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}

