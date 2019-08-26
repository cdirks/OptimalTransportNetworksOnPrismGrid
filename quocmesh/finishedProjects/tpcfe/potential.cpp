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

// Solve Lu = 0 with Dirichlet BC
// test program by TP

#define COMPUTE_SOLUTION
// #define USE_TPOS_MULTIGRID

#include <tpCFEStandardOp.h>
#include <preconditioner.h>
#include <fstream>
#include <math.h>
#include "affine.h"

static const int DEPTH   = 5;
static const int WIDTH   = ((1<<DEPTH)+1);
static const int WIDTHm2 = (WIDTH-2);
static const int TOTALSIZE    = (WIDTH*WIDTH*WIDTH);

#ifndef VERBOSE
#define VERBOSE
#include <solver.h>
#undef VERBOSE
#else
#include <solver.h>
#endif


static const double DIFF_COEFF_PLUS  = 10.0;
static const double DIFF_COEFF_MINUS =  1.0;
static const double TANG_A           =  0.0;
static const double TANG_B           =  0.0;

void outputData ( qc::ScalarArray<double, qc::QC_3D> &v, const char *name, const bool scale = false ) {
  char mask[1024];
  sprintf ( mask, "%s-%%03d.pgm", name );
  std::cerr << name << " Min/Max is " << v.getMinValue() << " " << v.getMaxValue() << std::endl;
  v.setQuietMode ( true );
  if ( scale ) {
    v.saveSlices ( mask, qc::QC_Z, qc::PGM_UNSIGNED_CHAR_BINARY, NULL, aol::SCALE );
  } else {
    v.saveSlices ( mask, qc::QC_Z, qc::PGM_UNSIGNED_CHAR_BINARY );
  }

#if 0
  ofstream o ( name );

  qc::ScalarArray<float, qc::QC_3D> tempArray ( v.getNumX(), v.getNumY(), v.getNumZ() );
  for ( int i = 0; i < v.size(); i++ ) {
    tempArray[i] = ( float ) v[i];
  }
  tempArray.saveRaw ( o );
#endif

}


typedef tpcfe::CFEGrid< double, tpcfe::CFE_TPOS >                       GridType;
typedef tpcfe::CFEConfigurator < GridType, aol::SparseMatrix<double> >  ConfiguratorType;
typedef ConfiguratorType::MatrixType                                    MatrixType;
typedef tpcfe::CFEMassOp < ConfiguratorType >                           MassOpType;
typedef tpcfe::CFEStiffOp < ConfiguratorType >                          StiffOpType;
typedef tpcfe::CFEMassOpWI < ConfiguratorType >                         WMassOpType;
typedef tpcfe::CFEStiffOpWI < ConfiguratorType >                        WStiffOpType;

int main ( int, char** ) {
  try {
    GridType grid ( DEPTH );

    qc::ScalarArray<double, qc::QC_3D>
    myInterface ( grid ),
    affine ( grid ),
    bc ( grid ),
    f ( grid );

    qc::ScalarArray<double, qc::QC_3D>
      rhs ( grid ),
      img ( grid );
    qc::AArray<double, qc::QC_3D>
      diffCoeff ( grid );

    WStiffOpType                          stiffOp ( diffCoeff, grid, aol::ASSEMBLED );
    WMassOpType                           massOp ( diffCoeff, grid, aol::ONTHEFLY );

    tpcfe::AffineFunction<double>         affineFct;

    affineFct.setCoeffs ( DIFF_COEFF_PLUS, DIFF_COEFF_MINUS, TANG_A, TANG_B );
    affineFct.setType ( tpcfe::AffineFunction<double>::SPHERICAL_INTERFACE );
    affineFct.setRadius ( 1.0 / 3.0 );
    //      affineFct.setType(AffineFunction<double>::ALIGNED_PLANAR_INTERFACE_X);
    affineFct.setWidth ( WIDTH );

    std::cerr << affineFct << std::endl;

    affineFct.createInterface ( &myInterface, &affine, &diffCoeff );

    /*    outputData ( diffCoeff,    "out/diffCoeff.dat" ); // this does not work any more */
    outputData ( affine,       "out/affine.dat" );
    outputData ( myInterface,  "out/interface.dat" );

    // Initialize the grid
    grid.addStructureFrom ( myInterface );

    grid.detectAndInitVirtualNodes ( diffCoeff );

    // Set the BC
    for ( int x = 0; x < WIDTH; x++ ) {
      for ( int y = 0; y < WIDTH; y++ ) {
        for ( int z = 0; z < WIDTH; z++ ) {
          if ( x == 0 || x == WIDTH - 1 ) bc.set ( x, y, z, x + 1.0 ); // only at two opposite faces??
        }
      }
    }

#ifdef COMPUTE_SOLUTION
    // Compute right hand side
    stiffOp.apply ( bc, rhs );
    rhs *= -1;

    outputData ( rhs, "out/rhs.dat", true );
    outputData ( bc,  "out/bc.dat",  true );

    // Adjust matrix to boundary conditions
    std::cerr << "Adjusting to boundary conditions ... ";
    MatrixType &matrix = stiffOp.getMatrixRef();

    const int tenPercent = TOTALSIZE / 10;
    for ( int i = 0; i < TOTALSIZE; i++ ) {
      if ( bc.get ( i ) != 0 ) {
        matrix.setRowColToZero ( i );
        matrix.set ( i, i, 1.0 );
        rhs.set ( i, 0.0 );
      }
      if ( i % tenPercent == 0 ) std::cerr << ( i / tenPercent ) *10 << " % ... ";
    }
    std::cerr << " done." << std::endl;

    // Invert system of equations
    {
      aol::DiagonalMatrix<double> diag ( TOTALSIZE );
      for ( int k = 0; k < TOTALSIZE; k++ ) {
        const double value = matrix.get ( k, k );
        diag.set ( k, k, value );
      }
      diag.invert ( false );

      aol::PCGInverse<aol::Vector<double> > cg ( matrix, diag, 1e-20, 400 );
      cg.setQuietMode ( false );
      // cg.setVerboseMode( 5 );

      cg.apply ( rhs, img );

#ifdef USE_TPOS_MULTIGRID
      // TODO: test multigrid solver vs. pcg solver

      // how come this works in combination with enforcing boundary conditions??
      qc::ScalarArray<double, qc::QC_3D> img2 ( grid ), img3 ( grid );
      my_tpcfe_multigrid_sparse tpos_mg ( grid, matrix );
      tpos_mg.setVerboseMode ( 5 );
      tpos_mg.apply ( rhs, img2 );

      img3  = img ;
      img3 -= img2;

      cerr << "Norm of cg solution = " << img.norm() / img.size() << ", norm of mg solution = " << img2.norm() / img2.size() << ", norm of difference = " << img3.norm() / img3.size() << endl;
#endif

    }

    // Readjust to boundary conditions
    img += bc;
#endif

    // Output data
    outputData ( img, "out/test.dat" );

  } catch ( aol::Exception &e ) {
    e.dump();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
