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

// a parameter study how removing random individual rods from a 3D grid of rods affects stiffness of the system

#include <tpCFEElastOp.h>
#include <multiArray.h>
#include <parameterParser.h>
#include <shapeLevelsetGenerator.h>
#include <tpCFEMultigrid.h>
#include <tpCFEUtils.h>
#include <tpCFELevelsets.h>

#include "computeForce.h"

typedef double RealType;
typedef qc::BitArray<qc::QC_3D> bitArray;
typedef qc::MultiArray<RealType, qc::QC_3D> MultiArrayType;

typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD, tpcfe::IsotropicElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> > ConfiguratorType;

static const aol::GaussSeidelSweepingMode gsmode = aol::GAUSS_SEIDEL_FORWARD;
// static const aol::GaussSeidelSweepingMode gsmode = aol::GAUSS_SEIDEL_ZEBRA2; // useful if we parallelize

int main ( int /*argc*/, char **argv ) {
  try {

    aol::StopWatch::_suppressOpenMPWarning = true;

    aol::ParameterParser params ( argv[1] );
    params.dump();

    const int
      solversteps       = params.getInt ( "solversteps" ),
      vmode             = params.getInt ( "vmode" ),
      num_smooth        = params.getInt ( "num_smooth" ),
      coarse_level      = params.getInt ( "coarse_level");

    const RealType
      relax             = 1.0,
      epsilon           = params.getDouble ( "epsilon" );


    // parameters for the compression experiment
    const RealType
      shift             = -1.0,
      E                 = params.getDouble( "E" ),
      nu                = params.getDouble( "nu" ),
      lambda            = tpcfe::computeLambda ( E, nu ),
      mu                = tpcfe::computeMu ( E, nu );


    // parameters for the multigrid solver
    const int
      level             = params.getInt("level"),
      rand_steps        = params.getInt("rand_steps"),
      num_diam          = params.getDimSize( "d_x" ),
      num_perc          = params.getDimSize( "p_x" ),
      num_fix           = params.getDimSize( "fix_dirs" ),
      num_shift         = params.getDimSize( "shift_dirs" ),
      num_rodnums       = params.getDimSize( "num_rodnums" );


    aol::MultiVector<RealType> diameters(3,0), percentages(3,0);
    params.getRealVec ( "d_x", diameters[0] );
    params.getRealVec ( "d_y", diameters[1] );
    params.getRealVec ( "d_z", diameters[2] );

    params.getRealVec ( "p_x", percentages[0] );
    params.getRealVec ( "p_y", percentages[1] );
    params.getRealVec ( "p_z", percentages[2] );

    aol::Vector<int> fix_dirs, shift_dirs, rodnums;
    params.getIntVec ( "fix_dirs", fix_dirs );
    params.getIntVec ( "shift_dirs", shift_dirs );

    params.getIntVec ( "num_rodnums", rodnums );

    for ( int i_diam = 0; i_diam < num_diam; ++i_diam ) {
      for ( int i_perc = 0; i_perc < num_perc; ++i_perc ) {

        for ( int i_rodnum = 0; i_rodnum < num_rodnums; ++i_rodnum ) {
          const int num_rods = rodnums [ i_rodnum ];

          for ( int i_fix = 0; i_fix < num_fix; ++i_fix ) {
            for ( int i_shift = 0; i_shift < num_shift; ++i_shift ) {

              for ( int random_step = 0; random_step < rand_steps; ++random_step ) {

                GridType grid ( level );

                qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

                if ( params.hasVariable ( "inputLevelset" ) ) {

                  char levelsetFilename[1024];
                  params.getString ( "inputLevelset", levelsetFilename );
                  cerr << "Loading levelset from file " << levelsetFilename << endl;
                  levelset.load( levelsetFilename );
                  // assumes that poempels have already been removed

                } else {

                  aol::RandomGenerator rnd ( random_step );
                  const unsigned int random_seed = rnd.rUnsignedInt();

                  cerr << "Using random seed " << random_seed << endl;

                  cerr << "Diameters: "
                       << diameters[0][i_diam] / ( 2.0 * num_rods ) << " "
                       << diameters[1][i_diam] / ( 2.0 * num_rods ) << " "
                       << diameters[2][i_diam] / ( 2.0 * num_rods ) << endl
                       << "Removal percentages: "
                       << percentages[0][i_perc] << " " << percentages[1][i_perc] << " " <<  percentages[2][i_perc] << endl
                       << "Number of rods: " << num_rods << endl
                       << "Fixing " << fix_dirs[ i_fix ] << ", shifting " << shift_dirs [ i_fix ] << endl;

                  qc::ShapeLevelsetGenerator<RealType>::generateAnisoRandomRodsRemoved3DRodsLevelset ( levelset, num_rods,
                                                                                                       diameters[0][i_diam] / ( 2.0 * num_rods ), diameters[1][i_diam] / ( 2.0 * num_rods ), diameters[2][i_diam] / ( 2.0 * num_rods ),
                                                                                                       percentages[0][i_perc], percentages[1][i_perc], percentages[2][i_perc],
                                                                                                       random_seed );

                  tpcfe::AllFacePoempelRemover< RealType > poempelRemover;
                  poempelRemover.applySingle( levelset );

                  levelset.save ( "out/levelset.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

                }

                grid.setDomainFrom ( levelset );
                grid.detectAndInitVirtualNodes();

                bitArray DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
                DirichletMask.setAll ( false );

                MultiArrayType dirichletBCs ( grid ), u( grid ), rhs ( grid );

                tpcfe::setGeneralShearing ( grid, dirichletBCs, DirichletMask, fix_dirs[ i_fix ], shift_dirs[ i_shift ], shift );

                grid.setDirichletMask ( DirichletMask );
                grid.setDOFMaskFromDirichletAndDomainNodeMask();

                tpcfe::CFEElastOp< ConfiguratorType > elastop ( grid, lambda, mu );

                cerr << "Using 3 x " << grid.getNumDOFs() << " = " << 3 * grid.getNumDOFs() << " degrees of freedom, total number of nodes in cube: " << grid.getNumberOfNodes() << ", total volume (volume percentage): " << grid.getTotalInnerVolume() << endl;

                u -= dirichletBCs;

                elastop.restrictNonDomainEntries();

                grid.restrictToDomain ( u );

                elastop.apply ( u, rhs );

                elastop.restrictDirichletEntries();
                grid.restrictDirichletNodes ( rhs );

                aol::StopWatch timer_mg_solve;

                elastop.getBlockMatrixRef().getReference(0,0).printStatistics();

                tpcfe::CFEBlockMultigrid< tpcfe::CFEElastOp< ConfiguratorType >, aol::ONTHEFLY, mg::ONTHEFLY_COARSENING > mgsolver ( grid, elastop.getBlockMatrixRef(),
                                                                                                                                     coarse_level, num_smooth, num_smooth,
                                                                                                                                     gsmode, relax, mg::MULTIGRID_V_CYCLE,
                                                                                                                                     epsilon, solversteps );
                mgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

                mgsolver.setCoarseSolverSteps( 10000 );

                mgsolver.setVerboseMode( vmode );

                timer_mg_solve.start();

                mgsolver.apply( rhs, u );

                timer_mg_solve.stop();

                cerr << "MG solver took " << timer_mg_solve.elapsedWallClockTime() << ", " << timer_mg_solve.elapsedCpuTime() << endl;

                u += dirichletBCs;

                // perform force analysis
                cerr << "Force analysis for fixed direction " << fix_dirs[ i_fix ] << ", shift_direction " << shift_dirs[ i_shift ] << endl;
                aol::Vec3<RealType> normal ( 0.0, 0.0, 0.0 );
                normal[ fix_dirs[ i_fix ] ] = 1.0;
                cerr << ", normal = " << normal << endl;

                aol::Vector<RealType> fst_dummy;

                tpcfe::IsotropicElasticityCoefficient<RealType> ENuMinus ( E, nu ), ENuPlus ( aol::NumberTrait<RealType>::NaN, aol::NumberTrait<RealType>::NaN );
                qc::AArray< tpcfe::IsotropicElasticityCoefficient<RealType>, qc::QC_3D > isoCoeff ( grid );
                tpcfe::setCoeffForLevelset ( isoCoeff, levelset, ENuMinus, ENuPlus );

                computeForce ( grid, u, fix_dirs[ i_fix ], normal, 1e-16, -1, 1.0, isoCoeff, cerr, fst_dummy );

              }
            }
          }
        }
      }
    }

  } catch(aol::Exception &ex) {

    ex.dump();
    return( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );
}
