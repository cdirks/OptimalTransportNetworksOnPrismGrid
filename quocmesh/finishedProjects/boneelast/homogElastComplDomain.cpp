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

int main ( int, char ** ) {

  try {

    aol::StopWatch::_suppressOpenMPWarning = true;

    const int        solversteps = 1000;
    const RealType   epsilon     = 1.0e-16;

    // parameters for the compression experiment
    const RealType
      E                 = 8,
      nu                = 0.38,
      lambda            = tpcfe::computeLambda ( E, nu ),
      mu                = tpcfe::computeMu ( E, nu );

    char fnmask[1024];

    const int nCells = 14;

    qc::RectangularGrid<qc::QC_3D> setupGrid( qc::GridSize<qc::QC_3D>( 129, 129, 129 ) );

    const int mLocal = 2 * ( static_cast<int> ( floor ( 2 * nCells / sqrt(3.0) ) ) / 2 ); // "floor" to even number
    const int yoffset = static_cast<int> ( setupGrid.getNumY() * ( 1 - mLocal * sqrt(3.0) / ( 2 * nCells ) ) / 2 );
    cerr << "yoffset = " << yoffset << endl;

    qc::ScalarArray< RealType, qc::QC_3D > setupLevelset ( setupGrid );
    aol::Mat<2,2,RealType> thetaMat;
    thetaMat[0][0] = 0.04;    thetaMat[0][1] = 0.03;
    thetaMat[1][0] = 0.02;    thetaMat[1][1] = 0.01;

    qc::ShapeLevelsetGenerator<RealType>::generateGeomOrthoHoneycombLevelset ( setupLevelset, nCells, thetaMat );


    aol::Matrix33<RealType> sigmas[3][3], epsilons[3][3];

    for ( int fix_dir = 0; fix_dir < 3; ++fix_dir ) {
      for ( int shift_dir = 0; shift_dir < 3; ++shift_dir ) {

        cerr << "Fixing " << fix_dir << ", shifting " << shift_dir << endl;

        // not nice but works
        GridType grid( qc::GridSize<qc::QC_3D>( ( 3 * setupGrid.getNumX() ) / 7 + 1 , ( setupGrid.getNumY() - 2 * yoffset ) / 2, setupGrid.getNumZ() / 2 ) );
        qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

        levelset.padFrom ( setupLevelset );

        levelset.save ( "out/geomOrthoHoneycombLevelset.dat.bz2", qc::PGM_FLOAT_BINARY );

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

          aol::BlockGaussSeidelPreconditioner< aol::MultiVector<RealType>, ElastOpType::BlockMatrixType > prec ( elastOp.getBlockMatrixRef() );
          aol::PCGInverseProjectEqConstr< aol::MultiVector<RealType> > solver ( elastOp.getBlockMatrixRef(), prec, constrVec, neutralFunctions, epsilon, solversteps );

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

      }
    }

    tpcfe::convertAverageAndDumpSigmaTensor ( sigmas );

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
