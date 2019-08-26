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

#ifndef __TPCFETESTER_H
#define __TPCFETESTER_H

#include <scalarArray.h>

#include <tpCFEElastOp.h>
#include <tpCFEGrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEVirtualNode.h>

#include <solver.h>

namespace tpcfe {

/** A class implementing several methods which test the cfe code for correctness.
 *  When a piece of cfe-code has been rewritten, the basic functionality should be
 *  easily testable with this class.
 */

class CFETester {
public:


  template <tpcfe::ConstraintType CT>
  static bool CFEElastOpTester() {
    typedef tpcfe::CFEConfigurator < tpcfe::CFEGrid<double, CT>, aol::SparseMatrix<double> > ConfiguratorType;

    // We test on single hexahedral elements
    qc::ScalarArray<double, qc::QC_3D>  interfc ( 2, 2, 2 ), u ( 2, 2, 2 ), v ( 2, 2, 2 );
    qc::AArray<double, qc::QC_3D> dc ( 2, 2, 2 );
    dc.setAll ( 1 );

    bool result = true;
    int  i[8];

    // Run through all possible signatures
    for ( i[0] = -1; i[0] < 2; i[0] += 2 ) {
      interfc.set ( 0, i[0] );
      for ( i[1] = -1; i[1] < 2; i[1] += 2 ) {
        interfc.set ( 1, i[1] );
        for ( i[2] = -1; i[2] < 2; i[2] += 2 ) {
          interfc.set ( 2, i[2] );
          for ( i[3] = -1; i[3] < 2; i[3] += 2 ) {
            interfc.set ( 3, i[3] );
            for ( i[4] = -1; i[4] < 2; i[4] += 2 ) {
              interfc.set ( 4, i[4] );
              for ( i[5] = -1; i[5] < 2; i[5] += 2 ) {
                interfc.set ( 5, i[5] );
                for ( i[6] = -1; i[6] < 2; i[6] += 2 ) {
                  interfc.set ( 6, i[6] );
                  for ( i[7] = -1; i[7] < 2; i[7] += 2 ) {
                    interfc.set ( 7, i[7] );

                    cerr << "Checking config: ";
                    for ( int l = 0; l < 8; ++l ) cerr << ( ( i[l] < 0 ) ? '1' : '0' );
                    cerr << endl;

                    typename ConfiguratorType::GridType grid ( 0 );
                    grid.addStructureFrom ( interfc );
                    grid.detectVirtualNodes();
                    cerr << "Initializing vnodes" << flush;
                    grid.initVirtualNodes ( dc );
                    cerr << "done." << flush;

                    for ( int ii = 0; ii < 3; ++ii ) {
                      for ( int jj = 0; jj < 3; ++jj ) {
                        tpcfe::CFEMixedDerivativeOp<ConfiguratorType> op ( ii, jj, grid, aol::ONTHEFLY );
                        u.setAll ( 1 );
                        op.apply ( u, v );

                        if ( u*v > 1e-14 ) result = false;
                      }
                    }

                  }
                }
              }
            }
          }
        }
      }
    }

    return result;
  }

  /** Test CFEStandardOps by comparing the matrix which results from a non-interfaced grid
   *  with the one resulting from a cfe computation. If the entries differ more than
   *  the wanted accuracy the test fails. This test works only on interfaces that have diffusion
   *  coefficient 1 on either side of the interface.
   */
  template <typename OP_TYPE, int DEPTH, tpcfe::ConstraintType CT> static bool CFEStandardOpTester ( qc::ScalarArray<double, qc::QC_3D> &myInterface,
                                                                                                     qc::AArray<double, qc::QC_3D> &coeff,
      double epsilon ) {
    cerr << "Testing whether matrices for non-interfaced and interfaced geometry coincide" << endl;

    tpcfe::CFEGrid < double, CT > grid ( DEPTH );

    // Compute non-interfaced matrices
    OP_TYPE *popStd = new OP_TYPE ( grid, aol::ASSEMBLED );
    OP_TYPE &opStd = *popStd;
    opStd._quietMode = true;
    opStd.assembleMatrix();

    // Now add the structure to the grid
    grid.addStructureFrom ( myInterface );
    grid.detectAndInitVirtualNodes ( coeff );

    // Compute cfe matrices
    OP_TYPE *popCfe = new OP_TYPE ( grid, aol::ASSEMBLED );
    OP_TYPE &opCfe = *popCfe;
    opCfe._quietMode = true;
    opCfe.assembleMatrix();

    // Make comparison
    bool result = ( opCfe.isApproxEqual ( opStd, epsilon ) );

    // If not approx equal dump out the matrices for debugging purposes
    if ( !result ) {
      cerr << "Dumping std matrix...";
      ofstream stdOut ( "stdMatrix.dat" );
      opStd.getMatrixRef().print ( stdOut );
      cerr << "done. " << endl;

      cerr << "Dumping cfe matrix...";
      ofstream cfeOut ( "cfeMatrix.dat" );
      opCfe.getMatrixRef().print ( cfeOut );
      cerr << "done. " << endl;
    }

    cerr << "Test " << ( result ? "passed" : "failed" ) << endl;

    delete popStd;
    delete popCfe;

    return result;
  }

  template <typename OP_TYPE>
  static bool CFELocalAssemblyTesterLIEHR() {
    qc::ScalarArray<double, qc::QC_3D> structure ( 2, 2, 2 );
    qc::AArray<double, qc::QC_3D> diffCoeff ( 2, 2, 2 );


    OP_TYPE *zeroOp = NULL;

    double x[8] = { -1};
    int    counter = 0;

    diffCoeff.setAll ( 1 );

    const double MINUS_VALUE = -0.2;
    const double PLUS_VALUE  =  1.0;

    for ( x[0] = MINUS_VALUE; x[0] <= PLUS_VALUE; x[0] += PLUS_VALUE - MINUS_VALUE ) {
      for ( x[1] = MINUS_VALUE; x[1] <= PLUS_VALUE; x[1] += PLUS_VALUE - MINUS_VALUE ) {
        for ( x[2] = MINUS_VALUE; x[2] <= PLUS_VALUE; x[2] += PLUS_VALUE - MINUS_VALUE ) {
          for ( x[3] = MINUS_VALUE; x[3] <= PLUS_VALUE; x[3] += PLUS_VALUE - MINUS_VALUE ) {
            for ( x[4] = MINUS_VALUE; x[4] <= PLUS_VALUE; x[4] += PLUS_VALUE - MINUS_VALUE ) {
              for ( x[5] = MINUS_VALUE; x[5] <= PLUS_VALUE; x[5] += PLUS_VALUE - MINUS_VALUE ) {
                for ( x[6] = MINUS_VALUE; x[6] <= PLUS_VALUE; x[6] += PLUS_VALUE - MINUS_VALUE ) {
                  for ( x[7] = MINUS_VALUE; x[7] <= PLUS_VALUE; x[7] += PLUS_VALUE - MINUS_VALUE ) {

                    for ( int i = 0; i < 8; i++ ) structure.set ( i, x[i]*1.0 );

                    CFEGrid<double, CFE_LIEHR> grid ( 0 );
                    grid.setVerboseMode ( false );

                    grid.addStructureFrom ( structure );
                    grid.detectAndInitVirtualNodes ( diffCoeff );

                    if ( counter == 0 ) {
                      zeroOp = new OP_TYPE ( grid );
                      zeroOp->assembleMatrix();
                    } else {
                      OP_TYPE op ( grid );
                      op.assembleMatrix();

                      if ( !zeroOp->getMatrixRef().isApproxEqual ( op.getMatrixRef(), 1e-15 ) ) {
                        cerr << "THE FOLLOWING MATRICES DO NOT COINCIDE" << endl;
                        cerr << "ZERO OP MATRIX" << endl;
                        zeroOp->getMatrixRef().print ( cerr );

                        cerr << "OTHER MATRIX" << endl;
                        op.getMatrixRef().print ( cerr );

                        cerr << "PRESS RETURN";
                        getchar();
                        return false;
                      }
                    }

                    ++counter;
                  }
                }
              }
            }
          }
        }
      }
    }

    cerr << "Test passed: All matrices coincide" << endl;
    return true;
  }

  /** Test CFEStandardOps by computing the FE solution of the simple homogeneous Dirichlet problem
   *  \f$  -\Delta u = f \f$, where \f$ f(x,y,z) = 2(1-y)y(1-z)z + 2(1-x)x(1-z)z + 2(1-x)x(1-y)y  \f$,
   *  and returning the \f$ L^2 \f$ error to the true solution \f$ u(x,y,z) = x(1-x)y(1-y)z(1-z) \f$.
   *  This method can be used to check the convergence of the FE solution with decreasing grid size.
   *  To do so, you should have a loop over increasing DEPTH and check whether the return value
   *  of this method decreases with second order, i.e. increasing DEPTH by 1 would decrease the
   *  error by a factor of 4.
   */
  template <typename MASS_OP_TYPE, typename STIFF_OP_TYPE, int DEPTH, tpcfe::ConstraintType CT>
  static double CFEStandardConvergenceTester ( qc::ScalarArray<double, qc::QC_3D> &myInterface,
                                               qc::AArray<double, qc::QC_3D> &coeff,
                                               double /*epsilon*/ ) {

    typedef tpcfe::CFEConfigurator < tpcfe::CFEGrid<double, CT>, aol::SparseMatrix<double> > ConfiguratorType;
    typedef typename ConfiguratorType::MatrixType MatrixType;

    const int width = ( ( 1 << DEPTH ) + 1 );
    const int size = ( width * width * width );

    cerr << "Testing convergence with standard Dirichlet problem" << endl;

    typename ConfiguratorType::GridType  grid ( DEPTH );
    qc::ScalarArray<double, qc::QC_3D> rhs ( grid ), f ( grid );

    // Now add the structure to the grid
    grid.addStructureFrom ( myInterface );
    grid.detectAndInitVirtualNodes ( coeff );

    // Compute non-interfaced matrices
    MASS_OP_TYPE  massOp ( grid, aol::ONTHEFLY );
    STIFF_OP_TYPE stiffOp ( grid, aol::ASSEMBLED );

    stiffOp.assembleMatrix();
    MatrixType &matrix = stiffOp.getMatrixRef();

    // Adjust matrix to BC and set right hand side
    const int tenPercent = size / 10;
    for ( int x = 0, k = 0; x < width; x++ )
      for ( int y = 0; y < width; y++ )
        for ( int z = 0; z < width; z++, k++ ) {
          const double X = x * 1.0 / ( width - 1.0 );
          const double Y = y * 1.0 / ( width - 1.0 );
          const double Z = z * 1.0 / ( width - 1.0 );
          if ( x == 0 || x == width - 1 ||
               y == 0 || y == width - 1 ||
               z == 0 || z == width - 1 ) {
            matrix.setRowColToZero ( k );
            matrix.set ( k, k, 1.0 );
          }
          f.set ( x, y, z, 2* ( 1.0 - Y ) *Y* ( 1.0 - Z ) *Z +
                  2* ( 1.0 - X ) *X* ( 1.0 - Z ) *Z +
                  2* ( 1.0 - X ) *X* ( 1.0 - Y ) *Y );
          if ( k % tenPercent == 0 ) cerr << ( k / tenPercent ) *10 << "%...";

        }

    // Integrate right hand side
    massOp.apply ( f, rhs );

    aol::CGInverse<aol::Vector<double> > cg ( matrix, 1e-20, 400 );
    cg.setStopping ( aol::STOPPING_ABSOLUTE ); // no idea whether this makes sense; prevent warning
    cg.setQuietMode ( false );
    f.setZero();
    cg.apply ( rhs, f );

    // Compute L^2 error
    double sum = 0.0;
    for ( int x = 0, k = 0; x < width; x++ )
      for ( int y = 0; y < width; y++ )
        for ( int z = 0; z < width; z++, k++ ) {

          const double X = x * 1.0 / ( width - 1.0 );
          const double Y = y * 1.0 / ( width - 1.0 );
          const double Z = z * 1.0 / ( width - 1.0 );

          double value = X * ( 1.0 - X ) * Y * ( 1.0 - Y ) * Z * ( 1.0 - Z );

          sum += aol::Sqr ( value - f.get ( k ) );
        }

    cerr << "Error is " << sum * aol::Cub ( grid.H() ) << endl;

    return sum;
  }




};

}

#endif
