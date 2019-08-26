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


// reusing makes iterations slower, but always converge. using "discontinuous" boundary condition decomposition makes iterations faster but not converge sometimes (limited number of observations)
#define REUSE_TIMESTEP_FOR_BC 1

// toggle saving ply files, always save temperature datasets
// #define SAVE_PLY 1

// use PCG or Multigrid?
// #define USE_MULTIGRID_SOLVER 1


// move this function somewhere useful!!
template< typename GridType, typename ArrayType >
void writeTemperatureSlicePLYs ( const GridType &grid, const ArrayType &temperature, const char* fnmask ) {
  typedef typename GridType::RealType RealType;
  for ( short diri = 0; diri < 3; ++diri ) {
    const qc::Comp dir = static_cast<qc::Comp>( diri );
    // for ( short pos = 0; pos < grid.getNumXYZ(); ++ pos ) {
    // short pos = grid.getNumXYZ() / 2; {
    for ( short pos = 5; pos < 15; ++ pos ) {

      tpcfe::CFEInterfaceTriangulationWithScalarDataGenerator< GridType > itg ( grid, temperature );
      itg.determineSliceTriangulation ( dir, pos );

      char filename[1024];
      sprintf ( filename, fnmask, dir, pos );
      itg.saveToPLYFile ( filename );
    }

    for ( short pos = grid.getNumXYZ() - 16; pos < grid.getNumXYZ() - 6; ++ pos ) {

      tpcfe::CFEInterfaceTriangulationWithScalarDataGenerator< GridType > itg ( grid, temperature );
      itg.determineSliceTriangulation ( dir, pos );

      char filename[1024];
      sprintf ( filename, fnmask, dir, pos );
      itg.saveToPLYFile ( filename );
    }

  }
}


typedef double RealType;

static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;
typedef tpcfe::CFEGrid < RealType, AT >                   GridType;
typedef tpcfe::CFEHybridMatrix<GridType>                  MatrixType;

typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
typedef tpcfe::CFEStiffOpWI < ConfiguratorType >          WStiffOpType;
typedef tpcfe::CFEMassOpWI < ConfiguratorType >           WMassOpType;

int main ( int argc, char** argv ) {
  try {

    if ( argc != 2 ) {
      cerr << "usage: AlPmmADSHC <parameterfile>" << endl;
      return( EXIT_FAILURE );
    }

    aol::ParameterParser params ( argv[1] );
    params.dump();

    // "numerics parameters"
    const int
      level             = params.getInt ( "level" ),
      numTimestepsSlow  = params.getInt ( "numTimestepsSlow" ),
      numTimestepsFast  = params.getInt ( "numTimestepsFast" ),
      saveTSEverySlow   = params.getInt ( "saveTSEverySlow" ),
      saveTSEveryFast   = params.getInt ( "saveTSEveryFast" );


    const RealType
      tauSlow           = params.getDouble ( "tauSlow" ),
      tauFast           = params.getDouble ( "tauFast" );

    GridType grid ( level );

    char levelsetFilename[1024], plyFNMask[1024], dataFNMask[1024];
    params.getString ( "levelsetFilename", levelsetFilename );
    params.getString ( "plyFNMask", plyFNMask );
    params.getString ( "dataFNMask", dataFNMask );

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


    qc::ScalarArray<RealType, qc::QC_3D> temperatureNew ( grid ), temperatureOld ( grid ), temperatureBC ( grid ), rhs ( grid ), rhsModified ( grid ), dwarf ( grid );
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


    MatrixType massMat ( grid ), stiffMat ( grid );
    {
      WMassOpType massOp ( coeffRhoC, grid, aol::ASSEMBLED );
      WStiffOpType stiffOp ( coeffLambda, grid, aol::ASSEMBLED );
      massOp.assembleAddMatrix ( massMat );
      stiffOp.assembleAddMatrix ( stiffMat );
    }

    massMat *= aol::Cub ( sigma ); // correct scaling of units
    stiffMat *= sigma;

    temperatureOld.setAll ( params.getDouble ( "InitialValue" ) );


    {
#ifdef SAVE_PLY
      char filenamemask[1024];
      sprintf ( filenamemask, plyFNMask, 0 );
      writeTemperatureSlicePLYs ( grid, temperatureOld, filenamemask );
#endif

      char filename [1024];
      sprintf ( filename, dataFNMask, 0 );
      temperatureOld.save ( filename, qc::SaveTypeTrait<RealType>::BinarySaveType );
    }

    RealType realTime = 0;

    // perform slow time steps, then fast time steps
    for ( short speed = 0; speed < 2; ++speed ) {

      const RealType tau = ( speed == 0 ? tauSlow : tauFast );

      MatrixType sysMat ( grid ), sysMatModified ( grid );

      sysMat += stiffMat; // M + tau L
      sysMat *= tau;
      sysMat += massMat;

      sysMatModified += sysMat;
      tpcfe::enforceZeroDirichletBCs ( grid, sysMatModified );

#ifdef USE_MULTIGRID_SOLVER
      const int
        explicitLevel     = params.getInt ( "explicitLevel" ),
        numSmooth         = params.getInt ( "numSmooth" );

      tpcfe::CFEMultigrid< /*misuse: */ WStiffOpType, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > solver ( grid, sysMatModified, explicitLevel, numSmooth, numSmooth,
                                                                                                           aol::GAUSS_SEIDEL_FORWARD, 1.0,
                                                                                                           mg::MULTIGRID_V_CYCLE, 1.0e-16, 1000,
                                                                                                           aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
#else
      aol::SSORPreconditioner< aol::Vector<RealType>, MatrixType > prec ( sysMatModified );
      // tpcfe::CFEMultigridPreconditioner< WStiffOpType, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgPrecond ( grid, sysMatModified, explicitLevel, numSmooth, numSmooth );

      aol::PCGInverse< aol::Vector<RealType> > solver ( sysMatModified, prec, 1.0e-16, 100000 );
      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
#endif

      const int timestepStart = ( speed == 0 ? 0 : numTimestepsSlow );
      const int timestepStop  = ( speed == 0 ? numTimestepsSlow : numTimestepsSlow + numTimestepsFast );
      const int saveTSEvery = ( speed == 0 ? saveTSEverySlow : saveTSEveryFast );

      for ( int timestep = timestepStart; timestep < timestepStop; ++timestep ) {
        realTime += tau;
        cerr << "time = " << realTime << " seconds for timestep " << timestep << endl;

        massMat.apply ( temperatureOld, rhs );
        sysMat.apply ( temperatureBC, dwarf );
        rhs -= dwarf;

        rhsModified = rhs;
        tpcfe::enforceZeroDirichletBCs ( grid, rhsModified );
        tpcfe::enforceZeroDirichletBCs ( grid, temperatureNew ); // start with "fair" residuum

        solver.apply ( rhsModified, temperatureNew );

        cerr << "Memusage = " << ( aol::memusage() >> 20 ) << " MiB" << endl;

        temperatureNew += temperatureBC;

        temperatureOld = temperatureNew;

#ifdef REUSE_TIMESTEP_FOR_BC
        temperatureBC  = temperatureNew;
#endif

        if ( (timestep+1) % saveTSEvery == 0 ) {
#ifdef SAVE_PLY
          char filenamemask[1024];
          sprintf ( filenamemask, plyFNMask, timestep + 1 );
          writeTemperatureSlicePLYs ( grid, temperatureNew, filenamemask );
#endif

          char filename [1024];
          sprintf ( filename, dataFNMask, timestep + 1 );
          temperatureNew.save ( filename, qc::SaveTypeTrait<RealType>::BinarySaveType );
        }

      }
    }
  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}

