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

// this is a test program by TP showing how to use the CFE_DOMAIN mode


#include <tpCFEStandardOp.h>
#include <fstream>
#include <math.h>
#include "affine.h"

static const tpcfe::ConstraintType AT = tpcfe::CFE_DOMAIN;

static const int DEPTH = 6;
static const int WIDTH = ((1<<DEPTH)+1);
static const int WIDTHm2  = (WIDTH-2);
static const double TAU = (5.0/(WIDTH*WIDTH));


#ifndef VERBOSE
#define VERBOSE
#include <solver.h>
#undef VERBOSE
#else
#include <solver.h>
#endif

void outputData ( qc::ScalarArray<double, qc::QC_3D> &v, const char *name, const bool scale = false ) {
  double min, max;
  char mask[1024];
  sprintf ( mask, "%s-%%03d.pgm", name );
  std::cerr << name << " Min/Max is " << ( min = v.getMinValue() ) << " " << ( max = v.getMaxValue() ) << std::endl;
  v.setQuietMode ( true );
  if ( scale ) {
    v.saveSlices ( mask, qc::QC_Z, qc::PGM_UNSIGNED_CHAR_BINARY, NULL, aol::SCALE );
  } else {
    v.saveSlices ( mask, qc::QC_Z, qc::PGM_UNSIGNED_CHAR_BINARY );
  }
}

using namespace tpcfe;

typedef tpcfe::CFEConfigurator < tpcfe::CFEGrid < double, AT>, aol::SparseMatrix<double> > ConfiguratorType;
typedef tpcfe::CFEMassOp < ConfiguratorType > MassOpType;
typedef tpcfe::CFEStiffOp < ConfiguratorType > StiffOpType;
typedef qc::ScalarArray<double, qc::QC_3D> IT;
typedef aol::Vector<double>       VT;

int main ( int, char** )  {
  try {
    tpcfe::CFEGrid < double, AT >     grid ( DEPTH );
    IT                                myInterface ( grid );
    qc::AArray<double, qc::QC_3D>     dc ( grid );
    IT                                globalImg ( grid ), globalRhs ( grid );
    AffineFunction<double>            affineFct;


    // First we create a nice domain
    dc.setAll ( 1.0 );
    const int RADIUS = WIDTH / 3;
    const int FREQ = 4;
    for ( int z = 0; z < WIDTH; ++z ) {
      for ( int y = 0; y < WIDTH; ++y ) {
        for ( int x = 0; x < WIDTH; ++x ) {
          double X = x - WIDTH / 2;
          double Y = y - WIDTH / 2;
          double angle = atan ( X / ( Y + 0.001 ) );
          double iRadSqr = aol::Sqr ( RADIUS * ( 1 + 0.25 * sin ( FREQ * angle ) ) );

          double rad = aol::Sqr ( X ) + aol::Sqr ( Y );

          myInterface.set ( x, y, z, rad - iRadSqr );
        }
      }
    }
    outputData ( myInterface, "out/domain", true );

    // Initialize the grid
    grid.setDomainFrom ( myInterface );
    grid.detectAndInitVirtualNodes ( dc );

    const int numDofs = grid.largestDofIndex();
    std::cerr << "Grid contains " << numDofs << " nodes in CFE_DOMAIN mode." << endl;

    // Now create some initial data
    aol::Vector<double> img ( numDofs ), rhs ( numDofs );


    CFEGrid<double, AT>::DomainNodeIterator it ( grid );
    for ( it.begin(); !it.end(); ++it ) {
      const double rad =  aol::Sqr ( it.coord().x() - WIDTH / 2 ) +
                          aol::Sqr ( it.coord().y() - WIDTH / 2 ) +
                          aol::Sqr ( it.coord().z() - WIDTH / 2 );
      if ( rad < WIDTH / 4 ) img.set ( *it, 255 );
    }
    grid.mapVectorBack ( img, globalImg );
    outputData ( globalImg, "out/start" );

    std::cerr << "Initial data has been created" << endl;

    // Now build an operator for the heat equation
    MassOpType                   massOp ( grid, aol::ASSEMBLED );
    StiffOpType                  stiffOp ( grid, aol::ASSEMBLED );
    aol::LinCombOp<VT, VT>       myOp;
    myOp.appendReference ( massOp );
    myOp.appendReference ( stiffOp, TAU );

    massOp.assembleMatrix();
    stiffOp.assembleMatrix();

    // Create a solver with diagonal scaling
    aol::DiagonalMatrix<double> diag ( numDofs );
    for ( int k = 0; k < numDofs; k++ ) {
      const double value = massOp.getMatrixRef().get ( k, k ) + TAU * stiffOp.getMatrixRef().get ( k, k );
      diag.set ( k, k, value );
    }
    diag.invert ( false );

    aol::PCGInverse<VT > cg ( myOp, diag, 1e-20, 400 );
    cg.setQuietMode ( false );

    // Run some steps of the heat equation
    for ( int step = 0; step < 5; step++ ) {
      std::cerr << "Running heat equation step " << step << endl;
      // Compute right hand side
      massOp.apply ( img, rhs );
      img.setZero();

      grid.mapVectorBack ( rhs, globalImg );
      outputData ( globalImg, "out/rhs", true );

      img.setZero();
      cg.apply ( rhs, img );

      // Output data
      grid.mapVectorBack ( img, globalImg );
      outputData ( globalImg, "out/test" );
      std::cerr << "Press return" << flush;
      getchar();
    }
  } catch ( aol::Exception &e ) {
    e.dump();
    return 1;
  }
  return 0;
}
