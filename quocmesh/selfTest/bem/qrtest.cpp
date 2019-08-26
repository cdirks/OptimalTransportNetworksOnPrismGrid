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

#include "matrix.h"

#include "elasticOps.h"
#include "boundary.h"
#include "matrixInverse.h"
#include "parametric.h"
#include "blocksolver.h"

int main ( int, char** ) {
  try {

    cerr << "--- Testing FullMatrix<double> ..." << endl;

    int n = 8, dim = n * n;
    aol::FullMatrix<double> m ( dim, dim );
    for ( int i = 0; i < dim; ++i ) {
      m.set ( i, i, -4 );
      if ( ( i % n ) != 0 ) m.set ( i, i - 1, 1 );
      if ( ( i % n ) != n - 1 ) m.set ( i, i + 1, 1 );
      if ( i >= n ) m.set ( i, i - n, 1 );
      if ( i + n < dim ) m.set ( i, i + n, 1 );
    }
    aol::LUInverse<double> mi1 ( m );
    aol::QRInverse<double> mi2 ( m );

    aol::Vector<double> b ( dim );
    aol::Vector<double> x1 ( dim );
    aol::Vector<double> x2 ( dim );
    for ( int i = 0; i < dim; ++i )
      b [i] = sin ( static_cast <double> ( i ) );

    mi1.apply ( b, x1 );
    mi2.apply ( b, x2 );

    x1 -= x2;
    cerr << "Difference of LU and QR decomposition is " << x1.norm () << endl;

    aol::FullMatrix<double> q = mi2.getQ (), r = mi2.getR ();
    aol::FullMatrix<double> p ( q.getNumRows (), q.getNumCols () );
    p.makeProduct ( q, r );

    double e = 0;
    for ( int i = 0; i < dim; ++i )
      for ( int j = 0; j < dim; ++j ) {
        double e1 = m.get ( i, j ) - p.get ( i, j );
        e += e1 * e1;
      }

    cerr << "Difference of M and QR is " << sqrt ( e ) << endl;

    aol::QRGivensTridiag<aol::FullMatrix<double> > qrit ( m );
    aol::FullMatrix<double> qq (dim, dim), dd (dim, dim);
    qrit.getQ ( qq );
    qrit.getTridiagFull ( dd );
    aol::FullMatrix<double> qqt = qq;
    qqt.transpose ();

    aol::FullMatrix<double> p1 ( q.getNumRows (), q.getNumCols () ), p2 ( q.getNumRows (), q.getNumCols () );
    p2.makeProduct ( qq, dd ); p1.makeProduct ( p2, qqt );

    double ee = 0;
    for ( int i = 0; i < dim; ++i )
      for ( int j = 0; j < dim; ++j ) {
        double e1 = m.get ( i, j ) - p1.get ( i, j );
        ee += e1 * e1;
      }

    cerr << "Difference of M and QDQ^T is " << sqrt ( ee ) << endl;

    aol::Vector<double> bb ( 3 );
    {
      aol::FullMatrix<double> m ( 3, 2 );
      m.ref ( 0, 0 ) = 1; m.ref ( 0, 1 ) = 0;
      m.ref ( 1, 0 ) = 0; m.ref ( 1, 1 ) = 1;
      m.ref ( 2, 0 ) = 1; m.ref ( 2, 1 ) = 1;

      aol::Vector<double> x ( 2 ), b ( 3 );
      b [0] = 1;
      b [1] = 2;
      b [2] = 3;

      aol::QRInverse<double> mi ( m );
      mi.apply ( b, x );
      m.apply ( x, bb );
      bb -= b;

      cerr << "Residual of Least-Square Problem is " << bb.norm () << endl;
    }

    cerr << "--- Testing BlockSolver<double,int> ..." << endl;

    int ni = 20, no_4 = 6, no = 4 * no_4; n = ni + no;

    bm::Boundary<bm::ParaParticle<double> > boundary, outer;
    bm::ParaParticle<double> blob ( aol::Vec2<double> ( 0, 0 ), 0.5, -0.5, ni );
    bm::ParaParticle<double> box ( aol::Vec2<double> ( 0, 0 ), aol::Vec2<double> ( 1, 1 ), no / 4 );
    boundary.push_back ( blob );
    outer.push_back ( blob );
    outer.push_back ( box );

    aol::ElasticTensor<double> elast ( 0.72, 0.2, aol::ElasticTensor<double>::YOUNGPOISSON );
    elast.isoOffset ( 1E-6 );
    bm::ElasticGreen<double> green ( elast );
    aol::Matrix22<double> bareps ( 1, 0, 0, -1 );

    bm::ElasticOperator<bm::ParaParticle<double>, bm::ElasticSolution > singopM ( outer, green, true );
    bm::ElasticOperator<bm::ParaParticle<double>, bm::ElasticSolutionDiff > doubopM ( outer, green, true );
    for ( int i = 0; i < n; ++i ) for ( int j = 0; j < n; ++j ) {
        doubopM.set ( i, j + n, 0 ); doubopM.set ( i + n, j, 0 );
      }
    aol::Vector<double> rhsM ( 2 * n );
    rhsM.setZero ();

    bm::ElasticOperator<bm::ParaParticle<double>, bm::ElasticSolution > singopP ( boundary, green );
    bm::ElasticOperator<bm::ParaParticle<double>, bm::ElasticSolutionDiff > doubopP ( boundary, green );
    for ( int i = 0; i < ni; ++i ) for ( int j = 0; j < ni; ++j ) {
        doubopP.set ( i, j + ni, 0 ); doubopP.set ( i + ni, j, 0 );
      }
    bm::ElasticEigenstress<bm::ParaParticle<double> > eigen ( boundary, green, bareps, false );
    aol::Vector<double> rhsP ( 2 * ni );
    singopP.apply ( eigen, rhsP );

    bm::DisplacementNormalization<bm::ParaParticle<double> > norm ( boundary, true );
    aol::Vector<double> rhsN ( 2 ); rhsN.setZero ();

    aol::Vec2<double> shift ( 1, 0 );

    aol::Vector<double> def ( 2 * ni ), nst ( 2 * no );

    aol::BlockSolver<double> system;
    system.appendBlockRef ( singopM, doubopM, rhsM );
    system.appendBlockRef ( singopP, doubopP, rhsP );
    system.appendBlockRef ( norm, rhsN );

    system.appendEqualBlock ( 0, 0, ni, 1, 0 );
    system.appendEqualBlock ( 0, n, ni + n, 1, ni );
    system.appendEqualBlock ( 0, 2 * n, 2 * n + ni, 1, 2 * ni );
    system.appendEqualBlock ( 0, 3 * n, 3 * n + ni, 1, 3 * ni );
    system.appendEqualBlock ( 1, 2 * ni, 4 * ni, 2, 0 );

    system.appendPeriodicSquare ( 0, ni, no / 4, true );
    system.appendPeriodicSquare ( 0, n + ni, no / 4, true );
    system.appendPeriodicShiftedSquareX ( 0, 2 * n + ni, no / 4, shift [0], false, true );
    system.appendPeriodicShiftedSquareY ( 0, 3 * n + ni, no / 4, shift [1], false, true );

    aol::Vector<double> sol_slow ( 4 * n + 4 * ni + 2 * ni ), sol ( 4 * n + 4 * ni + 2 * ni );

    aol::StopWatch watch;
    watch.start ();
    system.apply ( system.getRHS (), sol_slow );
    watch.stop ();
    double slow = watch.elapsedCpuTime ();

    watch.start ();
    system.apply_fast ( system.getRHS_fast (), sol );
    watch.stop ();
    double fast = watch.elapsedCpuTime ();

    sol -= sol_slow;

    cerr << "Error of block solver is " << sol.norm () << endl;
    cerr << "Fast block solver took " << fast << " seconds, while slow one took " << slow << " seconds." << endl;

    if ( x1.norm () + sqrt ( e ) + sqrt ( ee ) + bb.norm () < 1E-4 ) {
      aol::printSelfTestSuccessMessage ( "--                      BEM qrtest Test Successful                            --" );
    } else {
      aol::printSelfTestFailureMessage ( "!!                      BEM qrtest Test FAILED                                !!" );
      throw aol::Exception ( "Test failed", __FILE__, __LINE__ );
    }
  } catch ( aol::Exception& ) {
    aol::callSystemPauseIfNecessaryOnPlatform();
    return 23;
  }

  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
