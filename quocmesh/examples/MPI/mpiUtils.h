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

#ifndef __MPIUTILS_H
#define __MPIUTILS_H

#include <multiArray.h>

namespace qc {

void cart_map ( int reorder, int *procgrid, int *myloc, int procneigh[3][2], int *grid2proc );

void scatter_data ( const double* bins, int num_bins[3], double* bins_local, int num_bins_local[3], int my_bins_start[3] );

void gather_data ( double* bins, int num_bins[3], const double* bins_local, int my_bins_start[3], int my_bins_end[3], int root );

template <typename RealType, qc::Dimension Dim>
struct doCalculateCentralFDDivergenceMPI_PBC {
  static void apply ( const qc::MultiArray<RealType, Dim> &MArg, qc::ScalarArray<RealType, Dim> &Divergence );
};

template <typename RealType>
struct doCalculateCentralFDDivergenceMPI_PBC<RealType, qc::QC_3D> {
  static void apply ( const qc::MultiArray<RealType, qc::QC_3D> &MArg, qc::ScalarArray<RealType, qc::QC_3D> &Divergence ) {
    const int numX = Divergence.getNumX();
    const int numY = Divergence.getNumY();
    const int numZ = Divergence.getNumZ();
    // Straightforward and readable, but not most efficient implementation:
    Divergence.setZero();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int z = 1; z < numZ - 1; ++z ) {
      for ( int y = 1; y < numY - 1; ++y ) {
        for ( int x = 1; x < numX - 1; ++x ) {
          Divergence.add ( x, y, z,  MArg[0].get ( x + 1, y  , z ) );
          Divergence.add ( x, y, z, -MArg[0].get ( x - 1, y  , z ) );
          Divergence.add ( x, y, z,  MArg[1].get ( x  , y + 1, z ) );
          Divergence.add ( x, y, z, -MArg[1].get ( x  , y - 1, z ) );
          Divergence.add ( x, y, z,  MArg[2].get ( x  , y  , z + 1 ) );
          Divergence.add ( x, y, z, -MArg[2].get ( x  , y  , z - 1 ) );
          Divergence.set ( x, y, z, Divergence.get ( x, y, z ) * 0.5 );
        }
      }
    }
  }
};

template <typename RealType, qc::Dimension Dim>
void calculateCentralFDDivergenceMPI_PBC ( const qc::MultiArray<RealType, Dim> &MArg, qc::ScalarArray<RealType, Dim> &Divergence ) {
  doCalculateCentralFDDivergenceMPI_PBC<RealType, Dim>::apply ( MArg, Divergence );
}
template <typename RealType, qc::Dimension Dim>
struct doCalculateCentralFDGradientMPI_PBC {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &Arg, qc::MultiArray<RealType, qc::QC_2D> &Gradient );
};

//! \todo Rewrite to get rid of all the if statements inside the for loops.
template <typename RealType>
struct doCalculateCentralFDGradientMPI_PBC<RealType, qc::QC_3D> {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_3D> &Arg, qc::MultiArray<RealType, qc::QC_3D> &Gradient ) {
    const int numX = Arg.getNumX();
    const int numY = Arg.getNumY();
    const int numZ = Arg.getNumZ();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Straightforward and readable, but not most efficient implementation:
    for ( int z = 1; z < numZ - 1; ++z ) {
      for ( int y = 1; y < numY - 1; ++y ) {
        for ( int x = 1; x < numX - 1; ++x ) {
          Gradient[0].set ( x, y, z, 0.5 * ( Arg.get ( x + 1, y  , z  ) - Arg.get ( x - 1, y  , z  ) ) );
          Gradient[1].set ( x, y, z, 0.5 * ( Arg.get ( x  , y + 1, z  ) - Arg.get ( x  , y - 1, z  ) ) );
          Gradient[2].set ( x, y, z, 0.5 * ( Arg.get ( x  , y  , z + 1 ) - Arg.get ( x  , y  , z - 1 ) ) );
        }
      }
    }
  }
};

template <typename RealType, qc::Dimension Dim>
void calculateCentralFDGradientMPI_PBC ( const qc::ScalarArray<RealType, Dim> &Arg, qc::MultiArray<RealType, Dim> &Gradient ) {
  doCalculateCentralFDGradientMPI_PBC<RealType, Dim>::apply ( Arg, Gradient );
}

}
#endif // __MPIUTILS_H
