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

#ifndef VERBOSE
#define VERBOSE
#endif

#include <tpCFEStandardOp.h>
#include <tpCFEElastOp.h>
#include <tpCFETester.h>
#include <fstream>
#include <math.h>
#include <progressBar.h>
#include <preconditioner.h>

#ifdef USE_TPCFE_GRAPE
#include "GrapeGrid.h"
#endif

static const tpcfe::ConstraintType AT = tpcfe::CFE_DOMAIN;

static const int DEPTH = 4;
static const int WIDTH = ((1<<DEPTH)+1);
static const int WIDTHm2  = (WIDTH-2);
static const double TAU = (5.0/(WIDTH*WIDTH));
static const double H = (1.0/(WIDTH-1));
static const double H3 = (H*H*H);

#define fileName "bone.dat"

#include <solver.h>

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

typedef tpcfe::CFEConfigurator < tpcfe::CFEGrid< double, AT >, aol::SparseMatrix<double> > ConfiguratorType;
typedef tpcfe::CFEMassOp < ConfiguratorType > MassOpType;
typedef tpcfe::CFEStiffOp < ConfiguratorType > StiffOpType;
typedef tpcfe::CFEMixedDerivativeOp < ConfiguratorType >          *PMDOp;
typedef aol::DiagonalMatrix<double>  *PDMT;
typedef ConfiguratorType::MatrixType *PMT;
typedef ConfiguratorType::MatrixType MatrixType;

typedef qc::ScalarArray<double, qc::QC_3D> IT;
typedef aol::Vector<double>       VT;

int main ( int, char** )  {
  double lambda = 28.;
  double mu = 42;

  std::cerr << "lambda = " << lambda << ", mu = " << mu << endl;

  try {
#if defined USE_EXTERNAL_GRAPE && defined USE_TPCFE_GRAPE
    tpcfe::GrapeGrid < double, AT >   grid ( DEPTH );
#else
    tpcfe::CFEGrid< double, AT >      grid ( DEPTH );
#endif
    IT                                myInterface ( grid );
    qc::AArray< double, qc::QC_3D>    dc ( grid );
    IT                                globalImg ( grid ), globalRhs ( grid );

    // First we create a rod shaped domain
    dc.setAll ( 1.0 );
#ifdef fileName
    myInterface.load ( fileName );
    myInterface.addToAll ( - 0.9 );
    myInterface *= -1;
#else

#if 1
    const int RADIUS = (WIDTH/4);
    const int RADIUS2 = (RADIUS*RADIUS);

    for ( int z = 0; z < WIDTH; ++z ) {
      for ( int y = 0; y < WIDTH; ++y ) {
        for ( int x = 0; x < WIDTH; ++x ) {
          double rad = aol::Sqr ( x - WIDTH / 2 ) + aol::Sqr ( y - WIDTH / 2 );

          myInterface.set ( x, y, z, rad - RADIUS2 );
          /* if(fabs(x-WIDTH/2) < RADIUS && fabs(y-WIDTH/2) < RADIUS)
          myInterface.set(x, y, z, -1);
          else myInterface.set(x, y, z, 1);*/
        }
      }
    }
#else
    for ( int z = 0; z < WIDTH; ++z ) {
      for ( int y = 0; y < WIDTH; ++y ) {
        for ( int x = 0; x < WIDTH; ++x ) {
          //    myInterface.set(x, y, z, WIDTH/2-y-0.5);
          myInterface.set ( x, y, z, -1 );
        }
      }
    }
#endif
#endif    // of if 1

    outputData ( myInterface, "out/domain", true );

    // Initialize the grid
    grid.setDomainFrom ( myInterface );
    grid.detectAndInitVirtualNodes ( dc );

    const int numDofs = grid.largestDofIndex();
    aol::MultiVector<double>  u ( 3, numDofs ), rhs ( 3, numDofs ), bc ( 3, numDofs );

    std::cerr << "Problem has " << numDofs << " degrees of freedom" << endl;

    // Create operators
    MassOpType             opMM ( grid, aol::ASSEMBLED );
    PMDOp                  op[3][3];
    PMT                    mat[3][3];
    PDMT                   diagonal[3];

    CFEGrid<double, AT>::DomainNodeIterator it ( grid );

    aol::DiagonalBlockOp<double> blockMM ( opMM );
    opMM.assembleMatrix();

    aol::BlockOp<double>         blockop ( 3, 3 );
    aol::BlockOp<double>         blockDiagonal ( 3, 3 );
    int i, j;
    for ( i = 0; i < 3; ++i ) {
      diagonal[i] = new aol::DiagonalMatrix<double> ( numDofs );
      blockDiagonal.setReference ( i, i, *diagonal[i] );
      for ( j = 0; j < 3; ++j ) {
        op[i][j] = new tpcfe::CFEMixedDerivativeOp< ConfiguratorType > ( i, j, grid, aol::ASSEMBLED );

        mat[i][j] = op[i][j]->getConfig().createNewMatrix();
        blockop.setReference ( i, j, *mat[i][j] );
        mat[i][j]->setZero();
      }
    }


    for ( i = 0; i < 3; ++i ) {
      // Assemble diagonal blocks
      for ( j = 0; j < 3; ++j ) {
        op[j][j]->assembleAddMatrix ( *mat[i][i] );
      }

      *mat[i][i] *= ( mu / ( lambda + mu ) );
      op[i][i]->assembleAddMatrix ( *mat[i][i] );
      *mat[i][i] *= ( lambda + mu );

      // Assemble off-diagonal blocks
      for ( j = 0; j < 3; ++j ) {
        if ( i == j ) continue;
        op[i][j]->assembleAddMatrix ( *mat[i][j] );
        *mat[i][j] *= ( lambda / mu );
        op[j][i]->assembleAddMatrix ( *mat[i][j] );
        *mat[i][j] *= mu;
      }
    }

#define FORCE
#ifdef FORCE
    // Create force vector
    u.setZero();
    u[1].setAll ( -0.2 );
    blockMM.apply ( u, rhs ); // Kraft wird in Vektor rhs geschrieben
    // Use mass lumping
// rhs = u;
// rhs *= H3;
#else
    std::cerr << "Creating Dirichlet nodes..." << flush;
    bc.setZero();
    for ( it.begin(); !it.end(); ++it ) {
      if ( it.coord().z() == WIDTH - 1 ) {
        const int index = *it;
        bc[1][index] = 0.1;
      }
    }
    std::cerr << "done." << endl;
    blockop.apply ( bc, rhs );
    rhs *= -1;
#endif

    // Count dirichlet nodes first
    int numDir = 0;
    for ( it.begin(); !it.end(); ++it ) {

#ifdef FORCE
      if ( it.coord().z() == 0 ) {
#else
      if ( it.coord().z() == 0 || it.coord().z() == WIDTH - 1 ) {
#endif
        ++numDir;
      }
    }

    std::cerr << "Problem has " << numDir << " Dirichlet nodes" << endl;
    aol::ProgressBar<> pb ( "Adjusting matrix to boundary values" );
    pb.start ( numDir );
    for ( it.begin(); !it.end(); ++it ) {
#ifdef FORCE
      if ( it.coord().z() == 0 ) {
#else
      if ( it.coord().z() == 0 || it.coord().z() == WIDTH - 1 ) {
#endif
        pb++;
        const int index = *it;
        rhs[0][index] = 0;
        rhs[1][index] = 0;
        rhs[2][index] = 0;

        // Adjust matrices?
        for ( i = 0; i < 3; ++i ) {
          for ( j = 0; j < 3; ++j ) {
            mat[i][j]->setRowColToZero ( index );
            if ( i == j ) mat[i][j]->set ( index, index, 1.0 );
          }
        }
      }
    }
    std::cerr << endl;

#if 0
    // Fill the diagonal
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < numDofs; ++j ) {
        diagonal[i]->set ( j, j, mat[i][i]->get ( j, j ) );
      }
    }
#endif

    aol::BlockOp<double>        bprecond ( 3, 3 );
#if 1
    aol::SSORPreconditioner<aol::Vector<double>, MatrixType > precond00 ( *mat[0][0] );
    aol::SSORPreconditioner<aol::Vector<double>, MatrixType > precond11 ( *mat[1][1] );
    aol::SSORPreconditioner<aol::Vector<double>, MatrixType > precond22 ( *mat[2][2] );
#else
    aol::ILU0Preconditioner<double, MatrixType > precond00 ( *mat[0][0] );
    aol::ILU0Preconditioner<double, MatrixType > precond11 ( *mat[1][1] );
    aol::ILU0Preconditioner<double, MatrixType > precond22 ( *mat[2][2] );
#endif

    bprecond.setReference ( 0, 0, precond00 );
    bprecond.setReference ( 1, 1, precond11 );
    bprecond.setReference ( 2, 2, precond22 );

    // Solve the system
    aol::PCGInverse< aol::MultiVector<double> > inv ( blockop, bprecond, 1e-12, 400 );
    inv.setQuietMode ( false );
    inv.apply ( rhs, u );

#ifndef FORCE
    u += bc;
#endif

    // Now output or show data.
    std::cerr << "Min/max of u " << endl
    << u[0].getMinValue() << "/" << u[0].getMaxValue() << endl
    << u[1].getMinValue() << "/" << u[1].getMaxValue() << endl
    << u[2].getMinValue() << "/" << u[2].getMaxValue() << endl;

    std::cerr << "Value at virtual dirichlet nodes is ";
    for ( int l = 0; l < 3; ++l ) std::cerr << u[l][26] << " ,";
    std::cerr << std::endl;

    for ( int l = 0; l < 3; ++l ) std::cerr << u[l][29] << " ,";
    std::cerr << std::endl;

    for ( int l = 0; l < 3; ++l ) std::cerr << u[l][33] << " ,";
    std::cerr << std::endl;

    // Save output
    for ( int l = 0; l < 3; ++l ) {
      char filename[1024];
      sprintf ( filename, "/tmp/u-x%d.dat", l );
      ofstream out ( filename );
      out << numDofs << endl;
      for ( i = 0; i < numDofs; ++i ) out << u[l][i] << endl;
    }
#if defined USE_EXTERNAL_GRAPE && defined USE_TPCFE_GRAPE
    grid.showDeformedGrid ( u );

// grid.grapeVis(myInterface);
#endif
  }
  catch ( aol::Exception &e ) {
    e.dump();
    return 1;
  }
  return 0;
}

