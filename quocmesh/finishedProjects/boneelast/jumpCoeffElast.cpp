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
#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMultigrid.h>
#include <tpCFEUtils.h>

#include "computeForce.h"

#include "periodicMG.h"


typedef double RealType;
static const tpcfe::ConstraintType CT = tpcfe::CFE_TPOSELAST;
typedef tpcfe::CFEGrid < RealType, CT, tpcfe::IsotropicElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEHybridMatrix<GridType> MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType > ConfiguratorType;

typedef tpcfe::CFEJCElastOp< ConfiguratorType > ElastOpType;


int main ( const int, const char**  ) {

  try {

    aol::StopWatch timer;

#if 0
    const RealType
      EMinus  = 70.0, // Al; in GPa
      nuMinus = 0.35,
      EPlus   = 3.0,  // PMMA; in GPa
      nuPlus  = 0.38;

    const char* sfnmask = "out/AlPmmA2_L6_Def%d.dat.bz2";
    const int level = 6;
    GridType grid ( level, qc::QC_3D );

    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    levelset.load ( "datasets/AlPmmA2_L6.dat.bz2" );

    const int explicitLevel = 4;
#endif

#if 0
    const RealType
      EMinus  = 10.0, // Al; in GPa
      nuMinus = 0.1,
      EPlus   = 1.0,  // PMMA; in GPa
      nuPlus  = 0.3;

    // const int level = 4;
    // GridType grid ( level, qc::QC_3D );
    GridType grid ( qc::GridSize<qc::QC_3D> ( 25, 25, 40 ), qc::QC_3D );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset, 1./3. );

    // levelset.save ( "out/CubeSphere6.dat.bz2", qc::PGM_DOUBLE_BINARY );
    const char* sfnmask = "out/CubeSphere6_%d.dat.bz2";

    const int explicitLevel = 1;
#endif

#if 0
    const RealType
      EMinus  = 70.0, // Al; in GPa
      nuMinus = 0.35,
      EPlus   = 3.0,  // PMMA; in GPa
      nuPlus  = 0.38;

    //   const char* sfnmask = "out/AlPmmABrick060060096_Def%d.dat.bz2";
    const char* sfnmask = "out/AlPmmABrick075075120_Def%d.dat.bz2";

    // GridType grid ( qc::GridSize<qc::QC_3D> ( 60, 60, 96 ), qc::QC_3D );
    GridType grid ( qc::GridSize<qc::QC_3D> ( 75, 75, 120 ), qc::QC_3D );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    //  levelset.load ( "datasets/AlPMMA1Brick060060096.dat.bz2" );
    levelset.load ( "datasets/AlPMMA1Brick075075120.dat.bz2" );
    // grid.setAdjustVirtualNodes( 1.0e-2 );                          ???
#endif

#if 0
    const RealType
      EMinus  = 70.0, // Al; in GPa
      nuMinus = 0.35,
      EPlus   = 3.0,  // PMMA; in GPa
      nuPlus  = 0.38;

    const char* sfnmask = "out/AlPMMA5_60u_L7_Def%d.dat.bz2";

    GridType grid ( 7, qc::QC_3D );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    levelset.load ( "datasets/AlPMMA5_60u_L7.dat.bz2" );
#endif

#if 0
    const RealType
      EMinus  = 13.0, // Bone, GPa
      nuMinus = 0.32,
      EPlus   = 3.0,  // PMMA; in GPa
      nuPlus  = 0.38;

    const char* sfnmask = "out/961_T1Po2_070070100_035mu_Def%d.dat.bz2";

    GridType grid ( qc::GridSize<qc::QC_3D> ( 70, 70, 100 ), qc::QC_3D );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    levelset.load ( "datasets/961_T1Po2_070070100_035mu.dat.bz2" );
#endif

#if 1
    const RealType
      EMinus  = 13.0, // Bone, GPa
      nuMinus = 0.32,
      EPlus   = 3.0,  // PMMA; in GPa
      nuPlus  = 0.38;

    const char* sfnmask = "out/961_T1Po2_143143214_035mu_Def%d.dat.bz2";

    GridType grid ( qc::GridSize<qc::QC_3D> ( 143, 143, 214 )  );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    levelset.load ( "datasets/961_T1Po2_143143214_035mu.dat.bz2" );
#endif

    grid.addStructureFrom ( levelset );

    qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );
    cerr << EMinus << " " << nuMinus << " " << EPlus << " " << nuPlus << endl;
    const GridType::NodalCoeffType ENuMinus ( EMinus, nuMinus ), ENuPlus ( EPlus, nuPlus );
    tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
    grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-12, 1.0e-12 );

    qc::MultiArray< RealType, qc::QC_3D > dirichletBCs ( grid ), rhs ( grid ), soln ( grid ), dwarf ( grid );

    qc::BitArray<qc::QC_3D> DirichletMask ( grid );

    tpcfe::setGeneralShearing ( grid, dirichletBCs, DirichletMask, 2, 2, -0.01 );

    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    // tpcfe::CFEJCEMassOp< ConfiguratorType > massOp ( grid );
    ElastOpType elastOp ( grid, coeff );

    dwarf -= dirichletBCs;

    elastOp.apply ( dwarf, rhs );

    elastOp.restrictDirichletEntries();
    grid.restrictDirichletNodes ( rhs );

    timer.start();

    // tpcfe::CFEBlockMultigridPreconditioner< ElastOpType, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > prec ( grid, elastOp.getBlockMatrixRef(), explicitLevel, 3, 3 );
    // prec.setCoarseSolverSteps ( 0 );
    // aol::DiagonalPreconditioner < aol::MultiVector<RealType> > prec ( elastOp.getBlockMatrixRef() );
    // aol::BlockDiagonalPreconditioner < RealType, qc::QC_3D > prec ( elastOp.getBlockMatrixRef() );
    // aol::GeometricScalingPreconditioner< aol::MultiVector<RealType> > prec ( elastOp.getBlockMatrixRef() );
    // aol::SSORPreconditioner< aol::MultiVector<RealType>, MatrixType > prec ( elastOp.getBlockMatrixRef() );
    aol::BlockGaussSeidelPreconditioner < aol::MultiVector<RealType>, ElastOpType::BlockMatrixType > prec ( elastOp.getBlockMatrixRef() );
    prec.setSweepingModes ( aol::GAUSS_SEIDEL_ZEBRA2_FORWARD, aol::GAUSS_SEIDEL_ZEBRA2_BACKWARD );

    aol::PCGInverse< aol::MultiVector<RealType> > solver ( elastOp.getBlockMatrixRef(), prec, 1.0e-16, 200000 );
    solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );



    solver.apply ( rhs, soln );

    cerr << "Memusage = " << ( aol::memusage() >> 20 ) << " MiB" << endl;

    timer.stop();
    cerr << "Solver (including setup) took " << timer.elapsedWallClockTime() << " seconds" << endl;

    soln += dirichletBCs;

    soln.save ( sfnmask, qc::PGM_DOUBLE_BINARY );

    const RealType aleph = 6.81e-03;

    const unsigned int surfaceDirection = 2;
    const RealType tolerance = 1.0e-8;

    aol::Vec3<RealType> surfaceOuterNormal;
    surfaceOuterNormal [ surfaceDirection ] = 1.0;

    const unsigned int surfacePosition = -1;

    aol::Vector<double> dummy;

    tpcfe::computeForce ( grid, soln, surfaceDirection, surfaceOuterNormal, tolerance, surfacePosition, aleph, coeff, cerr, dummy );

  } catch ( aol::Exception &exc ) {
    exc.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
