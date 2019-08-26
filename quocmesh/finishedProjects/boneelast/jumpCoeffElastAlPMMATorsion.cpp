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

  const int level = 6;
  const char* datasetFN = "datasets/AlPMMA5_120u_L6.dat.bz2";
  const char* sfnmask = "out/AlPMMA5_120u_L6_Def%d.dat.bz2";

  const int solverSteps = 200000;
  const RealType solverEps = 1.0e-16;

  const RealType relaxedDAIVNStart = 1.0e-13, relaxedDAIVNStep = 1.0e-13;

  try {

    aol::StopWatch timer;

    const RealType
      EMinus  = 70.0, // Al; in GPa
      nuMinus = 0.35,
      EPlus   = 3.0,  // PMMA; in GPa
      nuPlus  = 0.38;

    GridType grid ( level );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

    levelset.load ( datasetFN );

    grid.addStructureFrom ( levelset );

    qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );
    cerr << EMinus << " " << nuMinus << " " << EPlus << " " << nuPlus << endl;
    const GridType::NodalCoeffType ENuMinus ( EMinus, nuMinus ), ENuPlus ( EPlus, nuPlus );
    tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
    grid.relaxedDetectAndInitVirtualNodes ( coeff, relaxedDAIVNStart, relaxedDAIVNStep );

    qc::MultiArray< RealType, qc::QC_3D > dirichletBCs ( grid ), rhs ( grid ), soln ( grid ), dwarf ( grid );

    qc::BitArray<qc::QC_3D> DirichletMask ( grid );

    tpcfe::setTorsion ( grid, dirichletBCs, DirichletMask, aol::DegreesToRadians ( 1.0 ) );

    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    ElastOpType elastOp ( grid, coeff );

    dwarf -= dirichletBCs;

    elastOp.apply ( dwarf, rhs );

    elastOp.restrictDirichletEntries();
    grid.restrictDirichletNodes ( rhs );

    timer.start();

    aol::BlockGaussSeidelPreconditioner < aol::MultiVector<RealType>, ElastOpType::BlockMatrixType > prec ( elastOp.getBlockMatrixRef() );
    prec.setSweepingModes ( aol::GAUSS_SEIDEL_ZEBRA2_FORWARD, aol::GAUSS_SEIDEL_ZEBRA2_BACKWARD );

    aol::PCGInverse< aol::MultiVector<RealType> > solver ( elastOp.getBlockMatrixRef(), prec, solverEps, solverSteps );
    solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );


    solver.apply ( rhs, soln );

    cerr << "Memusage = " << ( aol::memusage() >> 20 ) << " MiB" << endl;

    timer.stop();
    cerr << "Solver (including setup) took " << timer.elapsedWallClockTime() << " seconds" << endl;

    soln += dirichletBCs;

    soln.save ( sfnmask, qc::PGM_DOUBLE_BINARY );

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
