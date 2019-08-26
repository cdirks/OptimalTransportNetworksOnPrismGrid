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

// #undef USE_EXTERNAL_AHMED
// Uncomment this to use full matrices even if AHMED is available

#include "ahmed.h"
#include "elasticOps.h"
#include "parametric.h"
#include "matrixInverse.h"
#include "solver.h"
#include "ctrlCCatcher.h"
#include <fstream>

//! Test circle with isotropic expansion Dirichlet boundary
int main ( int argc, char* argv[] ) {

  aol::StopWatch timer;
  timer.start ();
  bool bench = false;
  string resultfilename = "/dev/null";

  try {

#ifdef USE_EXTERNAL_AHMED
#if AHMED_INTERFACE_VERSION == 0
    cerr << aol::color::blue << "--- Testing bm::LinElasticSolution for expanding circle with Dirichlet boundary (using old AHMED interface for hierarchical matrices) ... " << aol::color::reset << endl;
#else
    cerr << aol::color::blue << "--- Testing bm::LinElasticSolution for expanding circle with Dirichlet boundary (using new AHMED interface for hierarchical matrices) ... " << aol::color::reset << endl;
#endif
#else
    cerr << aol::color::blue << "--- Testing bm::LinElasticSolution for expanding circle with Dirichlet boundary (using fully populated matrices) ... " << aol::color::reset << endl;
#endif

    if ( aol::checkForBenchmarkArguments ( argc, argv, resultfilename ) ) bench = true;
    else if (argc == 2 && string (argv [1]) == string ("bench")) bench = true;

    // Set QuocMesh CTRL-C handler
    signal ( InterruptSignal, aol::ctrlCHandler );

    // Compute coefficients of Green's function
    double lambda = 1, mu = 1, isoOffset = 1e-8;
    aol::ElasticTensor<double> elast ( lambda, mu, aol::ElasticTensor<double>::LAMENAVIER );
    elast.isoOffset ( isoOffset );
    bm::ElasticGreen<double> green ( elast );

    for ( int side = 0; side < 2; ++side ) {

      bool exterior = ( side != 0 );

      if ( exterior ) cerr << "--- Testing exterior problem ... " << endl;
      else cerr << "--- Testing interior problem ... " << endl;

      double errold = 0;
      aol::resetCtrlCState();
#ifdef USE_EXTERNAL_AHMED
      for ( int n = 16; n <= (bench ? 2048 : 128) && !aol::getCtrlCState (); n *= 2 ) {
#else
      for ( int n = 16; n <= (bench ? 1024 : 128) && !aol::getCtrlCState (); n *= 2 ) {
#endif

  // Geometry
  double radius = 1;
  bm::Boundary<bm::ParaParticle<double> > bnd;
  bnd.insert ( bm::ParaParticle<double> ( aol::Vec2<double> ( 0, 0 ), radius, radius, n ) );
  aol::MultiVector<double> pos ( 2, n );
  bnd.getPositions ( pos );
  aol::MultiVector<double> norm ( 2, n );
  bnd.getNormals ( norm );

#ifdef USE_EXTERNAL_AHMED
  bm::AHMEDParams Hpar;
  Hpar.clusterEta = 1.25;
  Hpar.minCluster = 10;
  Hpar.lowRankEps = 1E-12;
  Hpar.precondEps = 1E-4;
  Hpar.maxRank = 20;
  Hpar.solverEps = 1E-12;
  Hpar.maxIter = 256;
  // agglomerate fails when efence is enabled
  Hpar.agglomerate = false;
#if AHMED_INTERFACE_VERSION == 0
  // Segment cluster tree for nodal basis functions
  bm::nodalClusterTree<bm::ParaParticle<double> > tree ( bnd, Hpar.clusterEta, Hpar.minCluster );
  tree.quadruplicate ();

  // Hierarchical matrices
  bm::HMatrixOp<bm::nodalMatrixEntry<bm::LinElasticSolution<bm::ParaParticle<double> > > >
    sop ( 2, green, tree, Hpar.lowRankEps, Hpar.maxRank, bm::RENORMALIZE_DIAGONAL_NOT, bm::RESORT_VECTOR_BY_SEGMENT_TREE, exterior ? 1 : -1 );
  bm::HMatrixOp<bm::nodalMatrixEntry<bm::LinElasticSolutionDiff<bm::ParaParticle<double> > > >
    dop ( 2, green, tree, Hpar.lowRankEps, Hpar.maxRank, bm::RENORMALIZE_DIAGONAL_RIGID_BODY_ARGUMENT, bm::RESORT_VECTOR_BY_SEGMENT_TREE, exterior ? -1 : 1 );
#else
  // Neumann values are unknown everywhere
  aol::BitVector neumannIndicator( n );
  neumannIndicator.setAll( true );
  bm::ElasticHMatrixOp< bm::ParaParticle<double>, bm::nodalCollNode<typename bm::ParaParticle<double>::ConstSegmentType>,
                        bm::LinElasticSolution<bm::ParaParticle<double> >, bm::LinElasticSolutionDiff<bm::ParaParticle<double> > >
    hmOp ( bnd, neumannIndicator, green, Hpar, bm::RENORMALIZE_DIAGONAL_RIGID_BODY_ARGUMENT, !exterior );
#endif
#else
  // Fully populated BEM operators
  bm::ElasticOperator<bm::ParaParticle<double>,bm::LinElasticSolution> sopf ( bnd, green, !exterior );
  bm::ElasticOperator<bm::ParaParticle<double>,bm::LinElasticSolutionDiff> dopf ( bnd, green, !exterior );
#endif

  // Boundary values
  double expansion = 1;
  aol::MultiVector<double> mdisp ( pos ), mtrac ( 2, n );
  mdisp *= expansion;
  aol::Vector<double> disp ( mdisp ), trac ( 2*n ), temp ( 2*n );

  // Solve
#ifdef USE_EXTERNAL_AHMED
#if AHMED_INTERFACE_VERSION == 0
  dop.apply ( disp, temp );
#else
  hmOp.mapBtoH( disp, trac );  // trac used as auxiliary vector
  hmOp.applyRhsOp ( trac, temp );
#endif
#else
  dopf.apply ( disp, temp );
#endif
  temp *= -1; // move to right hand side

#ifdef USE_EXTERNAL_AHMED
#if AHMED_INTERFACE_VERSION == 0
  bm::HLUPreconditioner<double> sopai ( sop.getHMatrix (), Hpar.precondEps, Hpar.maxRank );
  aol::PBiCGStabInverse<aol::Vector<double> > sopi ( sop, sopai, aol::Sqr(Hpar.solverEps), Hpar.maxIter, aol::STOPPING_RELATIVE_TO_RIGHT_HAND_SIDE );
  sopi.setQuietMode ( true );
  sopi.apply ( temp, trac );
#else
  hmOp.BiCGStabSolve( temp, disp );
  hmOp.mapHtoB( disp, trac );  // disp used as auxiliary vector
#endif
#else
  aol::QRInverse<double> sopi ( sopf );
  sopi.apply ( temp, trac );
#endif
  mtrac.copySplitFrom ( trac );

  // Compare
  double err = 0;
  for ( int i = 0; i < n; ++i ) {
    aol::Vec2<double> nu ( norm  [0][i], norm  [1][i] ), t_bem ( mtrac [0][i], mtrac [1][i] ), t_ex = exterior ? - 2 * mu * nu : 2 * ( lambda + mu ) * nu, e = t_bem - t_ex;
    double en = e.norm ();
    if ( err < en ) err = en;
  }
  std::cerr << "L^inf-Error in boundary forces " << aol::detailedFormat ( err ) << " when discretized by " << aol::intFormat ( n ) << " points";
  if ( errold ) std::cerr << ", order of convergence is " << aol::mixedFormat ( ( log ( errold ) - log ( err ) ) / log ( 2. ) );
  std::cerr << std::endl;
  errold = err;
      }

      if ( exterior ? ( errold > 2e-3 ) : ( errold > 4e-4 ) ) throw aol::Exception ( "Error too large", __FILE__, __LINE__ );
    }
  } catch ( aol::Exception& e ) {
    e.dump ();

    aol::printSelfTestFailureMessage ( "!!                   BEM dirichlet Test FAILED                                !!" );

    aol::callSystemPauseIfNecessaryOnPlatform();
    return 23;
  }

  aol::printSelfTestSuccessMessage ( "--                   BEM dirichlet Test Successful                            --" );

  timer.stop ();
  timer.printReport ( cerr );

  if (bench) {
    double cputime = timer.elapsedCpuTime();
    double walltime = timer.elapsedWallClockTime();
#ifdef USE_EXTERNAL_AHMED
    double default_runtime = 26*60+41;
#else
    double default_runtime = 35*60+50;
#endif

    int nupsi = 100 * default_runtime / cputime + 0.5, wupsi = 100 * default_runtime / walltime + 0.5;

    cerr << aol::color::blue << "Congratulations: your computer has achieved " << nupsi<< " nupsi and " << wupsi << " wupsi!" << aol::color::reset << endl;
    // cpu time and wall clock time performance, normalized so that 100 nupsi = brahma, single job

    // log results to file, for make bench
    aol::logBenchmarkResult ( "dirichlet", nupsi, wupsi, resultfilename );
  }

  aol::callSystemPauseIfNecessaryOnPlatform();

  return 0;
}
