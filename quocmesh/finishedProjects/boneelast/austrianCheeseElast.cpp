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

// elasticity on solid cubes with periodic holes (closed pores), like in a swiss cheese :-)

#include <tpCFEElastOp.h>
#include <multiArray.h>
#include <parameterParser.h>
#include <shapeLevelsetGenerator.h>
#include <tpCFEMultigrid.h>
#include <tpCFEUtils.h>
#include <tpCFELevelsets.h>


typedef qc::BitArray<qc::QC_3D> bitArray;
typedef qc::MultiArray<double, 3> MultiArrayType;

typedef tpcfe::CFEGrid < double, tpcfe::CFE_CD > GridType;
typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEBandMatrix<GridType> > ConfiguratorType;

static const aol::GaussSeidelSweepingMode gsmode = aol::GAUSS_SEIDEL_FORWARD;

int main ( int, char** ) {
  try {
    qc::ScalarArray<double, qc::QC_2D>::quietMode = 1;

    const int
      num_holes             = 16,
      level                 = 8,
      coarse_level          = 2,
      vmode                 = 2, // verbosity
      num_smooth            = 3,
      mgcycles              = 10000;

    const double
      radius            = 1.0 / ( 4 * (num_holes-1) ),
      relax             = 1.0,
      mg_epsilon        = 1e-16,
      zshift            = -0.01,
      E                 = 7e10,
      nu                = 0.33,
      lambda            = tpcfe::computeLambda ( E, nu ),
      mu                = tpcfe::computeMu ( E, nu );

    GridType grid ( level );

    qc::ScalarArray<double, qc::QC_3D> levelset ( grid );

    qc::ShapeLevelsetGenerator<double>::generateAustrianCheeseLevelset ( levelset, num_holes, radius );

    char levelset_filename[1024];
    sprintf( levelset_filename, "out/austrian_cheese_8.dat.bz2" );
    levelset.save ( levelset_filename, qc::PGM_DOUBLE_BINARY );

    grid.setDomainFrom ( levelset );
    grid.detectAndInitVirtualNodes();

    tpcfe::CFEElastOp< ConfiguratorType > elastop ( grid, lambda, mu );

    bitArray DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
    DirichletMask.setAll ( false );

    MultiArrayType dirichletBCs ( grid ), u( grid ), rhs ( grid );

    setZShift ( grid, dirichletBCs, DirichletMask, zshift ); // compression

    char dirichletBCs_fnmask[1024];
    sprintf( dirichletBCs_fnmask, "out/austrian_cheese_8_BC.dat.bz2" );
    dirichletBCs.save ( dirichletBCs_fnmask, qc::PGM_DOUBLE_BINARY );

    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    cerr << "Using 3 x " << grid.getNumDOFs() << " = " << 3 * grid.getNumDOFs() << " degrees of freedom, total number of nodes in cube: " << grid.getNumberOfNodes() << endl;

    u -= dirichletBCs;

    elastop.restrictNonDomainEntries();

    grid.restrictToDomain ( u );

    elastop.apply ( u, rhs );

    elastop.restrictDirichletEntries();
    grid.restrictDirichletNodes ( rhs );

    aol::StopWatch timer_mg_solve;

    tpcfe::CFEBlockMultigrid< tpcfe::CFEElastOp< ConfiguratorType >, aol::ONTHEFLY, mg::ONTHEFLY_COARSENING > mgsolver ( grid, elastop.getBlockMatrixRef(),
                                                                                                                         coarse_level, num_smooth, num_smooth,
                                                                                                                         gsmode, relax, mg::MULTIGRID_V_CYCLE,
                                                                                                                         mg_epsilon, mgcycles );
    mgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

    mgsolver.setCoarseSolverSteps( 20000 );

    mgsolver.setVerboseMode( vmode );

    timer_mg_solve.start();

    mgsolver.apply( rhs, u );

    timer_mg_solve.stop();

    cerr << "MG solver took " << timer_mg_solve.elapsedCpuTime()  << endl;

    u += dirichletBCs;

    cerr << "Starting save" << endl;
    char deformations_fnmask[1024];
    sprintf( deformations_fnmask, "out/austrian_cheese_8_def%%d.dat.bz2" );
    u.save ( deformations_fnmask, qc::PGM_DOUBLE_BINARY );
    cerr << "Saved." << endl;

  } catch(aol::Exception &ex) {

    ex.dump();
    return( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );
}
