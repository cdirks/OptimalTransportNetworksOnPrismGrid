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

#ifndef __TPCFE_UTILS_H
#define __TPCFE_UTILS_H


// this is a collection of utility functions
// probably, some of these includes are useless -- TODO: remove those...

// TODO: think about where to move the useful code

#include <tpCFEMultigrid.h>

#include <FEOpInterface.h>
#include <matrix.h>
#include <solver.h>
#include <progressBar.h>

#include <tpCFEStandardOp.h>
#include <tpCFEElastOp.h>
#include "affine.h"

#include <multigrid.h>
#include <smoother.h>

#include <sparseMatrices.h>

namespace tpcfe {

// a couple of utility functions that maybe should be moved somewhere else


template <typename RealType>
void analyzeMatrixStructure ( const aol::Matrix<RealType> &matrix, const int N ) {
  const RealType thres = 0.0;
  int R = matrix.getNumRows();
  aol::Vector<int> offsis ( 2*R );
  aol::Vector<RealType> offsvals ( 2*R );
  offsis.setZero();
  offsvals.setZero();

  aol::ProgressBar<> pb ("Reading Entries");
  pb.start( matrix.getNumRows() );

  for ( int i = 0; i < matrix.getNumRows(); ++i, pb++ ) {
    std::vector<typename aol::Row< RealType >::RowEntry > vec;
    matrix.makeRowEntries ( vec, i );
    for ( typename vector<typename aol::Row< RealType >::RowEntry >::const_iterator it = vec.begin(); it != vec.end(); ++it ) {
      if ( aol::Abs ( it->value ) > thres ) {
        int offs = it->col - i ;
        offsis[offs+R]++;
        offsvals[offs+R] += aol::Abs ( it->value );
      }
    }
  }

  cerr << endl;

  cerr << "Matrix structure: (attention: z component corresponds to N^2)" << endl;
  int nofs = 0;
  for ( int i = 0; i < 2*R; ++i ) {
    if ( offsis[i] != 0 ) {
      cout << aol::intFormat ( i - R ) << ": " << aol::intFormat ( offsis[i] ) << " " << aol::mixedFormat ( offsvals[i] ) << ", " ;
      const int j = i - R;
      int
      p = j / ( N * N ) ,
          q = ( j % ( N * N ) ) / N,
              r = j % N;
      if ( q > 3 ) {
        q -= N;
        p += 1;
      }
      if ( q < -3 ) {
        q += N;
        p -= 1;
      }

      if ( r > 3 ) {
        r -= N;
        q += 1;
      }
      if ( r < -3 ) {
        r += N;
        q -= 1;
      }

      cout << aol::intFormat ( j ) << " = " << aol::intFormat ( r ) << " + " << aol::intFormat ( q ) << " x N + " << aol::intFormat ( p ) << " x N^2 " << endl;
      nofs++;
    }
  }
  cerr << nofs << " nonzero bands." << endl;
}


//! not for the tpcfe testing cases ???
//! enforces zero Dirichlet boundary conditions on the boundary of the whole bounding cube.
template <typename RealType, typename GridType >
void enforceZeroDirichletBCsWholeCube ( const GridType &grid, aol::Matrix<RealType> &matrix, qc::ScalarArray<RealType, qc::QC_3D> &arr ) {
  if ( grid.getDimOfWorld() == 3 ) {
    aol::ProgressBar<> pb ( "Enforcing zero Dirichlet BCs" );
    pb.start ( grid.getNumberOfBoundaryNodes() );

    qc::FastILexMapper< qc::QC_3D > index_mapper ( grid );

    for ( qc::RectangularBoundaryIterator<qc::QC_3D> fbnit( grid ); fbnit.notAtEnd(); ++fbnit ) {
      pb++;
      const int ind = index_mapper.getGlobalIndex ( ( *fbnit ) );

      matrix.setRowColToZero ( ind );                                        // it is crucial that the matrix type used has an efficient implementation of setRowColToZero!
      matrix.set ( ind, ind, 1.0 );
      arr.set ( ind,  0.0 );
    }

    cerr << endl;
  } else {
    cerr << "enforcing boundary conditions is not implemented for dimension different from 3" << endl;
  }
}

template <typename RealType, typename GridType >
void enforceZeroDirichletBCs ( const GridType  &grid, aol::Matrix<RealType> &matrix, qc::ScalarArray<RealType, qc::QC_3D> &arr ) {

  aol::ProgressBar<> pb ( "Enforcing zero Dirichlet BCs" );
  pb.start ( grid.getDirichletMaskRef().numTrue() );

  for ( qc::RectangularIterator<qc::QC_3D> fnit ( grid ); fnit.notAtEnd(); ++fnit ) {
    if ( grid.isDirichletNode ( *fnit ) ) {
      pb++;
      const int ind = grid.getIndexMapperRef().getGlobalIndex ( ( *fnit ) );

      matrix.setRowColToZero ( ind );
      matrix.set ( ind, ind, 1.0 );
      arr.set ( ind,  0.0 );
    }
  }
  cerr << endl;
}

template <typename RealType, typename GridType >
void enforceZeroDirichletBCs ( const GridType  &grid, aol::Matrix<RealType> &matrix ) {

  aol::ProgressBar<> pb ( "Enforcing zero Dirichlet BCs" );
  pb.start ( grid.getDirichletMaskRef().numTrue() );

  for ( qc::RectangularIterator<qc::QC_3D> fnit ( grid ); fnit.notAtEnd(); ++fnit ) {
    if ( grid.isDirichletNode ( *fnit ) ) {
      pb++;
      const int ind = grid.getIndexMapperRef().getGlobalIndex ( ( *fnit ) );

      matrix.setRowColToZero ( ind );
      matrix.set ( ind, ind, 1.0 );
    }
  }
  cerr << endl;
}

template <typename RealType, typename GridType >
void enforceZeroDirichletBCs ( const GridType  &grid, qc::ScalarArray<RealType, qc::QC_3D> &arr ) {

  aol::ProgressBar<> pb ( "Enforcing zero Dirichlet BCs" );
  pb.start ( grid.getDirichletMaskRef().numTrue() );

  for ( qc::RectangularIterator<qc::QC_3D> fnit ( grid ); fnit.notAtEnd(); ++fnit ) {
    if ( grid.isDirichletNode ( *fnit ) ) {
      pb++;
      const int ind = grid.getIndexMapperRef().getGlobalIndex ( ( *fnit ) );
      arr.set ( ind,  0.0 );
    }
  }
  cerr << endl;
}


template <typename RealType>
double L2_norm ( const aol::Vector<RealType> &arg, const qc::GridDefinition& /*grid*/ ) {
  // of vector
  return ( sqrt ( arg * arg / arg.size() ) );
}



template < typename _RealType >
class QuocCompatGridDefinition : public qc::GridDefinition {
public:
  typedef _RealType RealType;
  QuocCompatGridDefinition ( const int depth ) : qc::GridDefinition ( depth, qc::QC_3D ){};
  int getNumXYZ () const {
    return ( getNumX() );
  }
};

#if 0
//! Compatibility class: CubicGrid (parent of CFE grid) providing RealType typedef
template < typename _RealType >
class QuocCompatCubicGrid : public qc::CubicGrid<qc::QC_3D> {
public:
  typedef _RealType RealType;
  QuocCompatCubicGrid ( const int depth ) : qc::CubicGrid<qc::QC_3D> ( depth ){};
};
#endif
} // end namespace

#endif
