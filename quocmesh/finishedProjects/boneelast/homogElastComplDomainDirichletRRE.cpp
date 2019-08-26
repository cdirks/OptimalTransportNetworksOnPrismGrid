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

// a parameter study how removing random individual rods from a 3D grid of rods affects the macroscopic orthotropic elasticity tensor

#include <tpCFEElastOp.h>
#include <multiArray.h>
#include <parameterParser.h>
#include <shapeLevelsetGenerator.h>
#include <tpCFEMultigrid.h>
#include <tpCFEUtils.h>

#include <tpCFELevelsets.h>

#include "computeForce.h"


typedef double RealType;

typedef qc::MultiArray<RealType, qc::QC_3D> MultiArrayType;

typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD > GridType;
typedef tpcfe::CFEPeriodicHybridMatrix< GridType > MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType > ConfiguratorType;

typedef tpcfe::CFEElastOp< ConfiguratorType >  ElastOpType;

static const aol::GaussSeidelSweepingMode gsmode = aol::GAUSS_SEIDEL_FORWARD;
// static const aol::GaussSeidelSweepingMode gsmode = aol::GAUSS_SEIDEL_ZEBRA2; // useful if we parallelize


int main ( int argc, char **argv ) {
  if ( argc != 2 ) {
    cerr << "usage: " << argv[0] << " <parameter file>" << endl;
    return ( EXIT_FAILURE );
  }

  try {

    aol::ParameterParser params ( argv[1] );

    const int
      solversteps       = ( params.hasVariable("solversteps") ? params.getInt ( "solversteps" ) : 1000 );

    const RealType
      epsilon           = ( params.hasVariable("epsilon") ? params.getDouble ( "epsilon" ) : 1.0e-16 );

    // parameters for the compression experiment
    const RealType
      E                 = params.getDouble( "E" ),
      nu                = params.getDouble( "nu" ),
      lambda            = tpcfe::computeLambda ( E, nu ),
      mu                = tpcfe::computeMu ( E, nu );


    const int
      level             = params.getInt("level"),
      num_rods          = params.getInt("num_rods"),
      rand_steps        = params.getInt("rand_steps"),
      num_diam          = params.getDimSize( "d_x" ),
      num_perc          = params.getDimSize( "p_x" );


    const int
      vmode             = 1, // verbosity
      num_smooth        = 3, // get from parameter file?
      coarse_level      = params.getInt("coarse_level");

    aol::MultiVector<RealType> diameters(3,0), percentages(3,0);
    params.getRealVec ( "d_x", diameters[0] );
    params.getRealVec ( "d_y", diameters[1] );
    params.getRealVec ( "d_z", diameters[2] );

    params.getRealVec ( "p_x", percentages[0] );
    params.getRealVec ( "p_y", percentages[1] );
    params.getRealVec ( "p_z", percentages[2] );


    for ( int i_diam = 0; i_diam < num_diam; ++i_diam ) {
      for ( int i_perc = 0; i_perc < num_perc; ++i_perc ) {

        for ( int random_step = 0; random_step < rand_steps; ++random_step ) {

          aol::RandomGenerator rnd ( random_step );
          const unsigned int random_seed = rnd.rUnsignedInt();

          cerr << "Using random seed " << random_seed << endl;

          GridType outer_grid ( level );

          qc::ScalarArray<RealType, qc::QC_3D> levelset ( outer_grid );

          cerr << diameters[0][i_diam] / ( 2.0 * num_rods ) << " "
               << diameters[1][i_diam] / ( 2.0 * num_rods ) << " "
               << diameters[2][i_diam] / ( 2.0 * num_rods ) << endl;

          qc::ShapeLevelsetGenerator<RealType>::generatePeriodicAnisoRandom3DRodsLevelset ( levelset, num_rods,
                                                                                            diameters[0][i_diam] / ( 2.0 * num_rods ), diameters[1][i_diam] / ( 2.0 * num_rods ), diameters[2][i_diam] / ( 2.0 * num_rods ),
                                                                                            percentages[0][i_perc], percentages[1][i_perc], percentages[2][i_perc],
                                                                                            random_seed );

          tpcfe::AllFacePoempelRemover<RealType> poempelRemover;
          poempelRemover.applySingle( levelset );

          levelset.save ( "out/levelset.dat.bz2", qc::PGM_DOUBLE_BINARY );

          aol::Matrix33<RealType> sigmas[4][3][3], epsilons[4][3][3];

          for ( int fix_dir = 0; fix_dir < 3; ++fix_dir ) {
            for ( int shift_dir = 0; shift_dir < 3; ++shift_dir ) {

              GridType grid ( level );

              grid.setDomainFrom ( levelset );
              grid.detectAndInitVirtualNodes();
              grid.setDOFMaskFromDirichletAndDomainNodeMask();

              qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
              DirichletMask.setAll ( false );

              qc::MultiArray< RealType, qc::QC_3D > DirichletBCs ( grid ), rhs ( grid ), soln ( grid );

              for ( qc::RectangularBoundaryIterator<qc::QC_3D> bbit ( grid ); bbit.notAtEnd(); ++bbit ) {
                DirichletMask.set ( *bbit, true );
                DirichletBCs[ shift_dir ].set ( *bbit, grid.H() * (*bbit)[ fix_dir ] );
                // other entries of DirichletBCs are zero
              }

              grid.setDirichletMask ( DirichletMask );
              grid.setDOFMaskFromDirichletAndDomainNodeMask();

              tpcfe::CFEElastOp< ConfiguratorType > elastop ( grid, lambda, mu );

              soln -= DirichletBCs;

              elastop.restrictNonDomainEntries();

              grid.restrictToDomain ( soln );

              elastop.apply ( soln, rhs );

              elastop.restrictDirichletEntries();
              grid.restrictDirichletNodes ( rhs );

              aol::StopWatch timer_mg_solve;

              elastop.getBlockMatrixRef().getReference(0,0).printStatistics();

              tpcfe::CFEBlockMultigrid< tpcfe::CFEElastOp< ConfiguratorType >, aol::ONTHEFLY, mg::ONTHEFLY_COARSENING > mgsolver ( grid, elastop.getBlockMatrixRef(),
                                                                                                                                   coarse_level, num_smooth, num_smooth,
                                                                                                                                   gsmode, 1.0, mg::MULTIGRID_V_CYCLE,
                                                                                                                                   epsilon, solversteps );
              mgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
              mgsolver.setCoarseSolverSteps( 10000 );
              mgsolver.setVerboseMode( vmode );
              timer_mg_solve.start();

              mgsolver.apply( rhs, soln );

              timer_mg_solve.stop();

              cerr << "MG solver took " << timer_mg_solve.elapsedWallClockTime() << ", " << timer_mg_solve.elapsedCpuTime() << endl;

              soln += DirichletBCs;

              for ( short bdlayer = 0; bdlayer < 4; ++bdlayer ) {
                const aol::Vec3<int>
                  lowerBnd (       bdlayer   *   grid.getNumX() / 8      ,       bdlayer   *   grid.getNumY() / 8      ,       bdlayer   *   grid.getNumZ() / 8       ),
                  upperBnd ( ( 8 - bdlayer ) * ( grid.getNumX() - 1 ) / 8, ( 8 - bdlayer ) * ( grid.getNumY() - 1 ) / 8, ( 8 - bdlayer ) * ( grid.getNumZ() - 1 ) / 8 );
                tpcfe::getSigmaEpsilonViaPartialTetTraversal ( grid, soln, 1.0, lambda, mu, lowerBnd, upperBnd, sigmas[ bdlayer ][ fix_dir ][ shift_dir ], epsilons[ bdlayer ][ fix_dir ][ shift_dir ] );

                cerr << "Evaluated on elements: " << lowerBnd << " --- " << upperBnd << endl
                     << "Sigma = " << endl << sigmas[ bdlayer ][ fix_dir ][ shift_dir ] << endl
                     << "Epsilon = " << endl << epsilons[ bdlayer ][ fix_dir ][ shift_dir ] << endl;
              }

              char fnmask[1024];
              sprintf ( fnmask, "out/RREdeformationDirichletBC_%d_%d_%%d.dat.bz2", fix_dir, shift_dir );
              soln.save ( fnmask, qc::SaveTypeTrait<RealType>::BinarySaveType );

            }
          }

          for ( short bdlayer = 0; bdlayer < 4; ++bdlayer ) {
            cerr << "Tensor evaluate with boundary layer " << bdlayer << " / 8" << endl;

            tpcfe::convertAverageAndDumpTensors ( sigmas[bdlayer], epsilons[bdlayer] );
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
