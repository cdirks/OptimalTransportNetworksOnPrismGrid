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

#include "elasticOps.h"
#include "parametric.h"
#include "matrixInverse.h"

template <template <class ParticleType> class SingleLayerOp, template <class ParticleType> class DoubleLayerOp>
void test ( double maxDefErr, double maxTracErr, int nmin = 16, int nmax = 128,
      double lambdaI = 6, double muI = 2, double lambdaE = 1, double  muE = 1,
      double epsilonI = 1, double isoOffset = 1e-8,
      double radius = 0.5 ) {

    // Compute coefficients of Green's function
    aol::ElasticTensor<double> elastI ( lambdaI, muI, aol::ElasticTensor<double>::LAMENAVIER );
    elastI.isoOffset ( isoOffset );
    bm::ElasticGreen<double> greenI ( elastI );
    aol::ElasticTensor<double> elastE ( lambdaE, muE, aol::ElasticTensor<double>::LAMENAVIER );
    elastE.isoOffset ( isoOffset );
    bm::ElasticGreen<double> greenE ( elastE );

    // Store errors for comparison
    double defErr = 0, tracErr = 0;

    // Geometry
    for ( int n = nmin; n <= nmax; n *= 2 ) {

      bm::Boundary<bm::ParaParticle<double> > bnd;
      bnd.insert ( bm::ParaParticle<double> ( aol::Vec2<double> ( 0, 0 ), radius, radius, n ) );

      // BEM operators for interior and exterior
      bm::ElasticOperator<bm::ParaParticle<double>,SingleLayerOp> sopI ( bnd, greenI, true );
      bm::ElasticOperator<bm::ParaParticle<double>,DoubleLayerOp> dopI ( bnd, greenI, true );
      aol::Matrix22<double> misfit ( epsilonI, 0, 0, epsilonI );
      bm::ElasticEigenstress<bm::ParaParticle<double> > pre ( bnd, greenI, misfit, bm::OperatorTrait<SingleLayerOp<bm::ParaParticle<double> > >::nodecentered, false );
      aol::Vector<double> rhsI ( 2*n ); sopI.apply ( pre, rhsI );

      bm::ElasticOperator<bm::ParaParticle<double>,SingleLayerOp> sopE ( bnd, greenE, false );
      bm::ElasticOperator<bm::ParaParticle<double>,DoubleLayerOp> dopE ( bnd, greenE, false );
      bm::ElasticEigenstress<bm::ParaParticle<double> > rhsE ( bnd, greenE );

      // Assemble and solve full system
      aol::FullMatrix<double> mat ( 4*n, 4*n );
      aol::Vector<double> sol ( 4*n ), rhs ( 4*n );
      mat.setBlock (   0, 0, sopI ); mat.setBlock (   0, 2*n, dopI ); rhs.setBlock (   0, rhsI );
      mat.setBlock ( 2*n, 0, sopE ); mat.setBlock ( 2*n, 2*n, dopE ); rhs.setBlock ( 2*n, rhsE );

      aol::QRInverse<double> inv ( mat );
      inv.apply ( rhs, sol );

      aol::MultiVector<double> defN ( 2, n );
      aol::MultiVector<double> tracN ( 2, n );
      sol.getBlock (   0, tracN [0] );
      sol.getBlock (   n, tracN [1] );
      sol.getBlock ( 2*n, defN [0] );
      sol.getBlock ( 3*n, defN [1] );

      // Compute analytical solution
      double c = epsilonI * ( lambdaI + muI ) / ( lambdaI + muI + muE );

      aol::MultiVector<double> position ( 2, n );
      aol::MultiVector<double> normal ( 2, n );
      aol::MultiVector<double> defA ( 2, n );
      aol::MultiVector<double> tracA ( 2, n );
      bnd.getPositions ( position, bm::OperatorTrait<SingleLayerOp<bm::ParaParticle<double> > >::nodecentered ? bm::ParaParticle<double>::START : bm::ParaParticle<double>::CENTER );
      bnd.getNormals ( normal, ! bm::OperatorTrait<SingleLayerOp<bm::ParaParticle<double> > >::nodecentered );

      defA = position; defA *= c;
      tracA = normal; tracA *= - 2 * c; // I and E just have different sign

      // Compute error
      defA += defN; tracA += tracN;
      double newDefErr = defA.getMaxAbsValue (), newTracErr = tracA.getMaxAbsValue ();
      cerr << "L^inf error in displacement " << aol::mixedFormat ( newDefErr  ) << " when discretized by " << aol::mixedFormat ( n ) << " points";
      if ( defErr ) { cerr << ", order of convergence was " << aol::mixedFormat ( ( log ( defErr ) - log ( newDefErr ) ) / log ( 2. ) ) << endl; }
      else { cerr << endl; };
      cerr << "L^inf error in traction     " << aol::mixedFormat ( newTracErr ) << " when discretized by " << aol::mixedFormat ( n ) << " points";
      if ( tracErr ) { cerr << ", order of convergence was " << aol::mixedFormat ( ( log ( tracErr ) - log ( newTracErr ) ) / log ( 2. ) ) << endl; }
      else { cerr << endl; };

      defErr = newDefErr; tracErr = newTracErr;
    }

    if ( defErr > maxDefErr || tracErr > maxTracErr ) throw aol::Exception ( "Error too large", __FILE__, __LINE__ );
}

//! Test Eshelby's solution for elastic inclusion, circular case
int main ( int /*argc*/, char* /*argv*/[] ) {
  try {

    cerr << "--- Testing bm::ElasticSolution against Eshelby's solution for circular inclusion ... " << endl;
    test<bm::ElasticSolution,bm::ElasticSolutionDiff> ( 6e-4, 3e-2, 16, 64 );

    // Why worse??
    cerr << "--- Testing bm::LinElasticSolution against Eshelby's solution for circular inclusion ... " << endl;
    test<bm::LinElasticSolution,bm::LinElasticSolutionDiff> ( 3e-4, 5e-3, 16, 64 );

  } catch ( aol::Exception& e ) {
    e.dump ();

    aol::printSelfTestFailureMessage ( "!!                     BEM eshelby Test FAILED                                !!" );
    aol::callSystemPauseIfNecessaryOnPlatform();
    return 23;
  }

  aol::printSelfTestSuccessMessage ( "--                     BEM eshelby Test Successful                            --" );
  aol::callSystemPauseIfNecessaryOnPlatform();

  return 0;
}
