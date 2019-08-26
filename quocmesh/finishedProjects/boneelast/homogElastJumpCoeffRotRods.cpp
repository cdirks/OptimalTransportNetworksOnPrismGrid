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

#include <multiArray.h>
#include <parameterParser.h>
#include <shapeLevelsetGenerator.h>

#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMultigrid.h>
#include <tpCFEPeriodicBC.h>
#include <tpCFEUtils.h>

#include "computeForce.h"

typedef double RealType;

typedef qc::BitArray<qc::QC_3D> bitArray;
typedef qc::MultiArray<RealType, qc::QC_3D> MultiArrayType;

static const tpcfe::ConstraintType CT = tpcfe::CFE_TPOSELAST;

typedef tpcfe::CFEGrid < RealType, CT, tpcfe::IsotropicElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEPeriodicHybridMatrix< GridType > MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType > ConfiguratorType;

typedef tpcfe::CFEJCEMassOp< ConfiguratorType >   MassOpType;
typedef tpcfe::CFEJCElastOp< ConfiguratorType >  ElastOpType;


int main ( int, char** ) {

  try {

    const RealType
      EMinus  = 10.0,
      nuMinus = 0.1,
      EPlus   = 1.0,
      nuPlus  = 0.3;

    const int level = 6, nRods = 5;

    aol::StopWatch timer;

    aol::Matrix33<RealType> sigmas[3][3], epsilons[3][3];

    for ( short fix_dir = 0; fix_dir < 3; ++fix_dir ) {
      for ( short shift_dir = 0; shift_dir < 3; ++shift_dir ) {

        cerr << "Fixing " << fix_dir << ", shifting " << shift_dir << endl;

        GridType grid ( level );

        qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

        qc::ShapeLevelsetGenerator<RealType>::generateXRotated3DRodsLevelset ( levelset, nRods, aol::Vec3<RealType>( 0.12 / nRods, 0.12 / nRods, 0.12 / nRods ), atan ( 1. / nRods ) );
        levelset.save ( "out/xRotatedRodsInterface.dat.bz2", qc::PGM_DOUBLE_BINARY );

        grid.addStructureFrom ( levelset );

        qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );
        cerr << EMinus << " " << nuMinus << " " << EPlus << " " << nuPlus << endl;
        const GridType::NodalCoeffType ENuMinus ( EMinus, nuMinus ), ENuPlus ( EPlus, nuPlus );
        tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
        grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-13, 1.0e-13 );
        grid.setDOFMaskFromDirichletAndDomainNodeMask();

        tpcfe::CFEPeriodicityHandler < GridType > periodicityHandler( grid );

        timer.start();
        // set up block mass matrix ...
        MassOpType massOp ( grid );
        periodicityHandler.periodicallyCollapseBlockMatrix ( massOp.getBlockMatrixRef() );
        // ... and elastOp
        ElastOpType elastOp ( grid, coeff );
        timer.stop();
        cerr << "Assembling matrices took " << timer.elapsedWallClockTime() << " seconds." << endl;

        qc::MultiArray< RealType, qc::QC_3D > rhs ( grid ), soln ( grid ), uSmooth ( grid );

        // set macroscopic part
        for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
          uSmooth[ shift_dir ].set ( *bit, 1.0 * grid.H() * (*bit)[ fix_dir ] );
        }

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
          massOp.apply ( neutralFunctions[i], constrVec[i] );
          const RealType volumeFactor = constrVec[i] * neutralFunctions[i];
          cerr << "Volume factor = " << aol::detailedFormat ( volumeFactor ) << endl;
          constrVec[i] /= volumeFactor;
          periodicityHandler.restrictToPresentDOFs ( constrVec[i] ); // this should not be necessary
          periodicityHandler.restrictPeriodicBC ( constrVec[i] );    // this should not be necessary
        }

        for ( int i = 0; i < 3; ++i )
          cerr << "constraint violations by uSmooth part: " << constrVec[i][0] * rhs[0] << " " << constrVec[i][1] * rhs[1] << " " << constrVec[i][2] * rhs[2] << endl;

        for ( int i = 0; i < 3; ++i )
          cerr << "constraint violation by uSmooth part: " << constrVec[i] * rhs << endl;


        elastOp.getBlockMatrixRef().getReference(0,0).printStatistics();

        {
          timer.start();

          aol::BlockGaussSeidelPreconditioner< aol::MultiVector<RealType>, ElastOpType::BlockMatrixType > prec ( elastOp.getBlockMatrixRef() );
          aol::PCGInverseProjectEqConstr< aol::MultiVector<RealType> > solver ( elastOp.getBlockMatrixRef(), prec, constrVec, neutralFunctions, 1.0e-16, 10000 );

          solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
          solver.apply ( rhs, soln );

          timer.stop();
          cerr << "Solving took " << timer.elapsedWallClockTime() << " seconds." << endl;
        }

        cerr << "Constraint satisfied? ";
        for ( int i = 0; i < 3; ++i )
          cerr << i << ": " << constrVec[i] * soln / soln.getTotalSize() << endl;

        periodicityHandler.extendPeriodicBC ( soln );

        soln += uSmooth;

        tpcfe::getSigmaEpsilonViaFullTetTraversal ( grid, soln, coeff, 1.0, sigmas[ fix_dir ][ shift_dir ], epsilons[ fix_dir ][ shift_dir ] );

        cerr << "Sigma = " << endl << sigmas[ fix_dir ][ shift_dir ] << endl;
        cerr << "Epsilon = " << endl << epsilons[ fix_dir ][ shift_dir ] << endl;

      }
    }

    tpcfe::convertAverageAndDumpSigmaTensor ( sigmas );

  } catch(aol::Exception &ex) {

    ex.dump();
    return( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );
}
