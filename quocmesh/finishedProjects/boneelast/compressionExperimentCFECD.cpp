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
#include <bitVector.h>
#include <parameterParser.h>

#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMultigrid.h>
#include <tpCFEUtils.h>



typedef qc::BitArray<qc::QC_3D> bitArray;

typedef tpcfe::CFEGrid < double, tpcfe::CFE_CD > GridType;
typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> > ConfiguratorType;

static const aol::GaussSeidelSweepingMode GSMODE = aol::GAUSS_SEIDEL_FORWARD;

int main ( int argc, char **argv ) {
  if ( argc != 2 ) {
    cerr << "usage: compression_experiment_cfecd PARAMETERFILE" << endl;
    return( EXIT_FAILURE );
  }

  try {
    qc::ScalarArray<double, qc::QC_2D>::quietMode = 1;

    aol::ParameterParser parser( argv[1] );

    parser.dump();

    const int
      COARSE_LEVEL          = 2, // CHANGED
      VMODE                 = parser.getInt( "vmode" ),
      NUM_SMOOTH            = parser.getInt( "numsmooth" ),
      depth                 = parser.getInt( "depth" ),
      mgcycles              = parser.getInt( "mgcycles" ),
      number_epsilon_steps  = parser.getInt( "number_epsilon_steps" );

    const double
      RELAX             = parser.getDouble( "relax" ),
      // eps               = parser.getDouble( "eps" ),
      step_epsilon      = parser.getDouble( "step_epsilon" ),
      aleph             = parser.getDouble( "aleph" ),                             // conversion from quoc-1 cube to SI units
      zshift            = aleph * parser.getDouble( "zshift" ),                    // in m
      E                 = parser.getDouble( "Emodulus" ),                          // in Pa = N / m^2
      nu                = parser.getDouble( "nu" ),                                // in 1
      lambda            = tpcfe::computeLambda ( E, nu ),
      mu                = tpcfe::computeMu ( E, nu );

    char deformations_fnmask[1024], intermediate_def_fnm[1024], dirichletBCs_fnm[1024];
    parser.getString ( "deformations_fnmask", deformations_fnmask );
    parser.getString ( "intermediate_def_fnm", intermediate_def_fnm );
    parser.getString ( "dirichletBCs_fnmask", dirichletBCs_fnm );

    cerr << "lambda = " << lambda << ", mu = " << mu << endl;

    GridType grid ( depth );

    qc::ScalarArray<double, qc::QC_3D> levelset ( grid );

    char geom_filename[1024];
    parser.getString( "geometry_filename", geom_filename );

    levelset.load ( geom_filename );

    grid.setDomainFrom ( levelset );
    grid.detectAndInitVirtualNodes();

    bitArray DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
    DirichletMask.setAll ( false );

    qc::MultiArray<double, 3> dirichletBCs ( grid ), u ( grid ), rhs ( grid ), dbpcgsoln( grid ), mgsoln( grid ), mgpcgsoln( grid );

    setZShift ( grid, dirichletBCs, DirichletMask, zshift ); // compression
    //    setXShift ( grid, dirichletBCs, DirichletMask, zshift ); // shearing

    dirichletBCs.save ( dirichletBCs_fnm, qc::PGM_DOUBLE_BINARY );

    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    cerr << "Using 3 x " << grid.getNumDOFs() << " = " << 3 * grid.getNumDOFs() << " degrees of freedom, total number of nodes in cube: " << grid.getNumberOfNodes() << endl;

    cerr << "Setting up TP's elasticity operator (based on mixed derivative ops)" << endl;
    tpcfe::CFEElastOp< ConfiguratorType > elastop ( grid, lambda, mu );

    u -= dirichletBCs;

    elastop.restrictNonDomainEntries();

    grid.restrictToDomain ( u );

    elastop.apply ( u, rhs );

    elastop.restrictDirichletEntries();
    grid.restrictDirichletNodes ( rhs );

    aol::StopWatch timer_mg_solve;

    //     elastop.getBlockMatrixRef().getReference(0,0).printStatistics();

    tpcfe::CFEBlockMultigrid< tpcfe::CFEElastOp< ConfiguratorType >, aol::ONTHEFLY, mg::ONTHEFLY_COARSENING > mgsolver ( grid, elastop.getBlockMatrixRef(),
                                                                                                                         COARSE_LEVEL, NUM_SMOOTH, NUM_SMOOTH,
                                                                                                                         GSMODE, RELAX, mg::MULTIGRID_V_CYCLE,
                                                                                                                         1.0, mgcycles / number_epsilon_steps );

    mgsolver.setCoarseSolverSteps( 50000 ); // CHANGED

    // using SolutionStepSaver, we automatically use relative stopping
    mgsolver.setVerboseMode( VMODE );

    timer_mg_solve.start();

    aol::SolutionStepSaver< tpcfe::CFEBlockMultigrid< tpcfe::CFEElastOp< ConfiguratorType >, aol::ONTHEFLY, mg::ONTHEFLY_COARSENING >, qc::MultiArray<double, 3> > sss( mgsolver, intermediate_def_fnm, step_epsilon, number_epsilon_steps );
    sss.apply( rhs, mgsoln );

    timer_mg_solve.stop();
    cerr << "MG solver took " << timer_mg_solve.elapsedCpuTime()  << endl;

    mgsoln += dirichletBCs;

    cerr << "Starting save" << endl;
    mgsoln.save ( deformations_fnmask, qc::PGM_DOUBLE_BINARY );
    cerr << "Saved." << endl;

  } catch(aol::Exception &ex) {

    ex.dump();
    return( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );

}
