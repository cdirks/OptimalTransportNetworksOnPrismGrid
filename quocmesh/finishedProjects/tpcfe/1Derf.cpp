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

#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>

#include "tpcfe_utils.h"

// strange observation .....
// reusing makes iterations slower, but always converge. using "discontinuous" boundary condition decomposition makes iterations faster but not converge sometimes (limited number of observations)
#define REUSE_TIMESTEP_FOR_BC 1


typedef double RealType;

static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;
typedef tpcfe::CFEGrid < RealType, AT >                   GridType;
typedef tpcfe::CFEHybridMatrix<GridType>                  MatrixType;

typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
typedef tpcfe::CFEStiffOpWI < ConfiguratorType >          WStiffOpType;
typedef tpcfe::CFEMassOpWI < ConfiguratorType >           WMassOpType;

int main ( int, char** ) {
  try {

    const short level = 7;

    GridType grid ( level );

    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    levelset.setAll ( -1.0 );

    grid.addStructureFrom ( levelset );

    const RealType
      sigma     = 1.0,
      lambda    = 1.0,
      rhoC      = 1.0;

    // set up weightInterface
    qc::AArray<RealType, qc::QC_3D> coeffLambda ( grid );
    coeffLambda.setAll ( lambda );

    grid.relaxedDetectAndInitVirtualNodes ( coeffLambda, 1.0e-15, 5.0e-16 );

    qc::AArray<RealType, qc::QC_3D> coeffRhoC ( grid );
    coeffRhoC.setAll ( rhoC );

    qc::ScalarArray<RealType, qc::QC_3D> temperatureNew ( grid ), temperatureOld ( grid ), temperatureBC ( grid ), rhs ( grid ), rhsModified ( grid ), dwarf ( grid );
    qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
    DirichletMask.setAll ( false );

    for ( int i = 0; i < grid.getNumX(); ++i ) {
      for ( int j = 0; j < grid.getNumX(); ++j ) {
        temperatureBC.set ( i, j,              0    , 1.0 );
        DirichletMask.set ( i, j,              0    , true );
      }
    }

    grid.setDirichletMask ( DirichletMask );


    MatrixType massMat ( grid ), stiffMat ( grid );
    {
      WMassOpType massOp ( coeffRhoC, grid, aol::ASSEMBLED );
      WStiffOpType stiffOp ( coeffLambda, grid, aol::ASSEMBLED );
      massOp.assembleAddMatrix ( massMat );
      stiffOp.assembleAddMatrix ( stiffMat );
    }

    massMat *= aol::Cub ( sigma ); // correct scaling of units
    stiffMat *= sigma;

    temperatureOld.setAll ( 0.0 );

    RealType realTime = 0;

    const RealType tau = aol::Sqr ( grid.H() );

    MatrixType sysMat ( grid ), sysMatModified ( grid );

    sysMat += stiffMat; // M + tau L
    sysMat *= tau;
    sysMat += massMat;

    sysMatModified += sysMat;
    tpcfe::enforceZeroDirichletBCs ( grid, sysMatModified );

    tpcfe::CFEMultigrid< /*misuse: */ WStiffOpType, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > solver ( grid, sysMatModified, /*explicitLevel*/ 1, /*numSmooth*/ 3, /*numSmooth*/ 3,
                                                                                                         aol::GAUSS_SEIDEL_FORWARD, 1.0,
                                                                                                         mg::MULTIGRID_V_CYCLE, 1.0e-16, 1000,
                                                                                                         aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

    for ( int timestep = 0; timestep < 1024; ++timestep ) {
      realTime += tau;
      cerr << "time = " << realTime << " seconds for timestep " << timestep << endl;

      massMat.apply ( temperatureOld, rhs );
      sysMat.apply ( temperatureBC, dwarf );
      rhs -= dwarf;

      rhsModified = rhs;
      tpcfe::enforceZeroDirichletBCs ( grid, rhsModified );

      solver.apply ( rhsModified, temperatureNew );

      temperatureNew += temperatureBC;

      temperatureOld = temperatureNew;

      temperatureBC  = temperatureNew;

      const short mx = grid.getNumX() / 2, my = grid.getNumY() / 2;
      cerr << "UVALUES    ";
      for ( short i = 0; i < grid.getNumZ(); ++i )
        cerr << aol::detailedFormat ( temperatureNew.get ( mx, my, i ) ) << " ";
      cerr << endl;

      cerr << "THEOVALUES ";
      for ( short i = 0; i < grid.getNumZ(); ++i ) {
        const RealType theoVal = 1.0 - aol::Erf ( sigma * grid.H() * i / ( 2 * sqrt ( lambda / rhoC ) * sqrt ( realTime ) ) );
        cerr << aol::detailedFormat ( theoVal ) << " ";
      }
      cerr << endl;

    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}

