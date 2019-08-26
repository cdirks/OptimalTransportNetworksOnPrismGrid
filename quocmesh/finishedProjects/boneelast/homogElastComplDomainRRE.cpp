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
#include <tpCFEPeriodicBC.h>

// #define PRINT_FIT_DATA

#include <tpCFELevelsets.h>

#include "computeForce.h"



typedef double RealType;

typedef qc::BitArray<qc::QC_3D> bitArray;
typedef qc::MultiArray<RealType, qc::QC_3D> MultiArrayType;

typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD > GridType;
typedef tpcfe::CFEPeriodicHybridMatrix< GridType > MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType > ConfiguratorType;

typedef tpcfe::CFEMassOp< ConfiguratorType >   MassOpType;
typedef tpcfe::CFEElastOp< ConfiguratorType >  ElastOpType;

static const aol::GaussSeidelSweepingMode gsmode = aol::GAUSS_SEIDEL_FORWARD;
// static const aol::GaussSeidelSweepingMode gsmode = aol::GAUSS_SEIDEL_ZEBRA2; // useful if we parallelize

int main ( int argc, char **argv ) {
  if ( argc != 2 ) {
    cerr << "usage: rod_remove_elast <parameter file>" << endl;
    return ( EXIT_FAILURE );
  }

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
      epsilon           = params.getDouble ( "epsilon" );

    // parameters for the compression experiment
    const RealType
      E                 = params.getDouble( "E" ),
      nu                = params.getDouble( "nu" ),
      lambda            = tpcfe::computeLambda ( E, nu ),
      mu                = tpcfe::computeMu ( E, nu );


    const int
      level             = params.getInt ( "level"),
      rand_steps        = params.getInt ( "rand_steps"),
      num_rods          = params.getInt ( "num_rods" ),
      num_diam          = params.getDimSize ( "d_x" ),
      num_perc          = params.getDimSize ( "p_x" );

    char fnmask[1024];

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

          GridType outer_grid ( level );
          qc::ScalarArray<RealType, qc::QC_3D> levelset ( outer_grid );

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


            cerr << diameters[0][i_diam] / ( 2.0 * num_rods ) << " "
                 << diameters[1][i_diam] / ( 2.0 * num_rods ) << " "
                 << diameters[2][i_diam] / ( 2.0 * num_rods ) << endl
                 << "Removal percentages: "
                 << percentages[0][i_perc] << " " << percentages[1][i_perc] << " " <<  percentages[2][i_perc] << endl
                 << "Number of rods: " << num_rods << endl;

            qc::ShapeLevelsetGenerator<RealType>::generatePeriodicAnisoRandom3DRodsLevelset ( levelset, num_rods,
                                                                                              diameters[0][i_diam] / ( 2.0 * num_rods ), diameters[1][i_diam] / ( 2.0 * num_rods ), diameters[2][i_diam] / ( 2.0 * num_rods ),
                                                                                              percentages[0][i_perc], percentages[1][i_perc], percentages[2][i_perc],
                                                                                              random_seed );

            tpcfe::AllFacePoempelRemover<RealType> poempelRemover;
            poempelRemover.applySingle( levelset );

            levelset.save ( "out/levelset.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

          }

          aol::Matrix33<RealType> sigmas[3][3], epsilons[3][3];

          for ( int fix_dir = 0; fix_dir < 3; ++fix_dir ) {
            for ( int shift_dir = 0; shift_dir < 3; ++shift_dir ) {

              cerr << "Fixing " << fix_dir << ", shifting " << shift_dir << endl;

              GridType grid ( level );

              grid.setDomainFrom ( levelset );
              grid.detectAndInitVirtualNodes();
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


              ElastOpType elastOp ( grid, lambda, mu );


              qc::MultiArray< RealType, qc::QC_3D > rhs ( grid ), soln ( grid ), uSmooth ( grid );

              // set macroscopic part
              for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
                uSmooth[ shift_dir ].set ( *bit, 1.0 * grid.H() * (*bit)[ fix_dir ] );
              }
              sprintf ( fnmask, "out/smoothDeformation_%d_%d_%%d.dat.bz2", fix_dir, shift_dir );
              uSmooth.save ( fnmask, qc::SaveTypeTrait<RealType>::BinarySaveType );

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
                tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::MultiVector<RealType> >::checkCorrectionResiduumNeutrality ( elastOp.getBlockMatrixRef(), neutralFunctions[i] ), 1e-8, "Correction direction neutral for residuum?", __FILE__, __LINE__ );


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
                cerr << "constraint violation by uSmooth part: " << constrVec[i] * rhs << endl;

              cerr << "Using 3 x " << grid.getNumDOFs() << " = " << 3 * grid.getNumDOFs() << " degrees of freedom, total number of nodes in cube: " << grid.getNumberOfNodes() << ", total volume (volume percentage): " << grid.getTotalInnerVolume() << endl;
              elastOp.getBlockMatrixRef().getReference(0,0).printStatistics();

              {
                aol::StopWatch timer_solve;
                timer_solve.start();

#if 0
                // BDB preconditioned CG solver with projection
                aol::BlockDiagonalPreconditioner< RealType, MatrixType > prec ( elastOp.getBlockMatrixRef() );
                aol::PCGInverseProjectEqConstr< aol::MultiVector<RealType> > solver ( elastOp.getBlockMatrixRef(), prec, constrVec, neutralFunctions, epsilon, solversteps );
#else

                const RealType
                  relax             = 1.0;

                tpcfe::CFEBlockMultigridProjectAvgConstr< ElastOpType, MassOpType > solver ( grid, elastOp.getBlockMatrixRef(), periodicityHandler, coarse_level, num_smooth, num_smooth, relax, mg::MULTIGRID_V_CYCLE, epsilon, solversteps );
                solver.setVerboseMode ( vmode );
#endif
                solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
                solver.apply ( rhs, soln );
                timer_solve.stop();

                cerr << "Solver took " << timer_solve.elapsedWallClockTime() << " seconds wall clock time, " << timer_solve.elapsedCpuTime() << " seconds cpu time." << endl;
              }

              cerr << "Constraint satisfied? ";
              for ( int i = 0; i < 3; ++i )
                cerr << i << ": " << constrVec[i] * soln / soln.getTotalSize() << endl;

              periodicityHandler.extendPeriodicBC ( soln );


              sprintf ( fnmask, "out/oscillatoryDeformation_%d_%d_%%d.dat.bz2", fix_dir, shift_dir );
              soln.save ( fnmask, qc::SaveTypeTrait<RealType>::BinarySaveType );

              {
                qc::ScalarArray<RealType, qc::QC_3D> osciNorms ( grid );
                soln.getPointWiseNorm ( osciNorms );
                sprintf ( fnmask, "out/oscillatoryDeformationNorm_%d_%d.dat.bz2", fix_dir, shift_dir );
                osciNorms.save ( fnmask, qc::SaveTypeTrait<RealType>::BinarySaveType );
              }

              soln += uSmooth;


              tpcfe::getSigmaEpsilonViaFullTetTraversal ( grid, soln, 1.0, lambda, mu, sigmas[ fix_dir ][ shift_dir ], epsilons[ fix_dir ][ shift_dir ] );


              cerr << "Sigma = " << endl << sigmas[ fix_dir ][ shift_dir ] << endl;
              cerr << "Epsilon = " << endl << epsilons[ fix_dir ][ shift_dir ] << endl;

              sprintf ( fnmask, "out/sumDeformation_%d_%d_%%d.dat.bz2", fix_dir, shift_dir );
              soln.save ( fnmask, qc::SaveTypeTrait<RealType>::BinarySaveType );

              if ( 0 ) {
                qc::ScalarArray<RealType, qc::QC_3D> sumNorms ( grid );
                soln.getPointWiseNorm ( sumNorms );
                sprintf ( fnmask, "out/sumDeformationNorm_%d_%d.dat.bz2", fix_dir, shift_dir );
                sumNorms.save ( fnmask, qc::SaveTypeTrait<RealType>::BinarySaveType );
              }

            }
          }

          tpcfe::convertAverageAndDumpSigmaTensor ( sigmas );

        }
      }
    }

  } catch(aol::Exception &ex) {

    ex.dump();
    return( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );
}
