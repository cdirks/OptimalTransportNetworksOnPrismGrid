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

#ifndef __FASTSWEEPING_H
#define __FASTSWEEPING_H

#include <quoc.h>

namespace qc {

class Queue {
public:
  Queue ( int nx, int ny, int nz ) {
    this->size = 0;
    this->capacity = aol::Max ( static_cast <int> ( nx * ny * nz * 0.1 ), 128 );
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;

    list = new int[3 * capacity];
    in_queue = new short[nx * ny * nz];

    memset ( in_queue, 0, sizeof ( short ) * nx * ny * nz );
  }

  ~Queue() {
    delete[] list;
    delete[] in_queue;
  }

  void add_point ( int idx_x, int idx_y, int idx_z ) {
    int idx = ILexCombine3 ( idx_x, idx_y, idx_z, nx, ny );
    //only add point if it is not in queue already
    if ( in_queue[ idx ] == 0 ) {
      if ( size == capacity - 1 ) {
        //full
        reallocate();
      }
      list[3 * size] = idx_x;
      list[3 * size + 1] = idx_y;
      list[3 * size + 2] = idx_z;
      in_queue[idx] = 1;
      ++size;
    }
  }

  void reset() {
    this->size = 0;
    memset ( in_queue, 0, sizeof ( short ) * nx * ny * nz );
  }

  short* in_queue; //1 if idx is already in queue; 0 otherwise
  int* list;
  int size;
  int nx, ny, nz;
  int capacity;

private:
  void reallocate() {
    int* tmp = new int[3 * capacity * 2];
    memcpy ( tmp, list, sizeof ( int ) * 3 * capacity );

    capacity *= 2;
    delete[] list;
    list = tmp;
  }
};


void fast_sweeping_queue_PBC ( double* grid, int nx, int ny, int nz, double h );

template <typename RealType>
void sort3 ( const RealType a, const RealType b, const RealType c, RealType* sortedCoeff );

/**
 * Solves ((x - a)^+)^2 + ((x - b)^+)^2 + ((x - c)^+)^2 = h for x. With (y)^+ := (y>0) ? y : 0
 */
template <typename RealType>
RealType solve ( const RealType a, const RealType b, const RealType c, const RealType h );

template <typename RealType>
void print_3d_array ( RealType* grid, int nx, int ny, int nz );

template <typename RealType, int dirX, int dirY, int dirZ>
void sweep3D ( RealType * __restrict__ grid, int nx, int ny, int nz, RealType h, int PBC ) {

  if ( PBC ) {
    const int nxy = nx * ny;
    const int NX = nx * 1.5;
    const int NY = ny * 1.5;
    const int NZ = nz * 1.5;

    for ( int i = ( dirZ ? 0 : NZ ); ( dirZ ? i < NZ : i >= 0 ); i += ( dirZ ? 1 : -1 ) ) {
      for ( int j = ( dirY ? 0 : NY ); ( dirY ? j < NY : j >= 0 ); j += ( dirY ? 1 : -1 ) ) {
        for ( int k = ( dirX ? 0 : NX ); ( dirX ? k < NX : k >= 0 ); k += ( dirX ? 1 : -1 ) ) {

          int idx = ( i % nz ) * nxy + ( j % ny ) * nx + ( k % nx );
          if ( grid[idx] != 0.0 ) {
            RealType a = aol::Min ( grid[ ( i - 1 + nz ) % nz * nxy + j % ny * nx + k % nx], grid[ ( i + 1 ) % nz * nxy + j % ny * nx + k % nx] );
            RealType b = aol::Min ( grid[ i % nz * nxy + ( j - 1 + ny ) % ny * nx + k % nx], grid[ i % nz * nxy + ( j + 1 ) % ny * nx + k % nx] );
            RealType c = aol::Min ( grid[ i % nz * nxy + j % ny * nx + ( k - 1 + nx ) % nx], grid[ i % nz * nxy + j % ny * nx + ( k + 1 ) % nx] );

            RealType new_value = solve ( a, b, c, h );
            new_value = aol::Min ( grid[idx], new_value );

            grid[idx] = new_value;
          }
        }
      }
    }
  } else {
    const int nxy = nx * ny;

    for ( int i = ( dirZ ? 1 : nz - 2 ); ( dirZ ? i < nz - 1 : i > 0 ); i += ( dirZ ? 1 : -1 ) ) {
      for ( int j = ( dirY ? 1 : ny - 2 ); ( dirY ? j < ny - 1 : j > 0 ); j += ( dirY ? 1 : -1 ) ) {
        for ( int k = ( dirX ? 1 : nx - 2 ); ( dirX ? k < nx - 1 : k > 0 ); k += ( dirX ? 1 : -1 ) ) {

          int idx = i * nxy + j * nx + k;
          if ( grid[idx] != 0.0 ) {
            RealType a = aol::Min ( grid[ ( i - 1 ) * nxy + j * nx + k], grid[ ( i + 1 ) * nxy + j * nx + k] );
            RealType b = aol::Min ( grid[ i * nxy + ( j - 1 ) * nx + k], grid[ i * nxy + ( j + 1 ) * nx + k] );
            RealType c = aol::Min ( grid[ i * nxy + j * nx + k - 1], grid[ i * nxy + j * nx + k + 1] );

            RealType new_value = solve ( a, b, c, h );
            new_value = aol::Min ( grid[idx], new_value );

            grid[idx] = new_value;
          }
        }
      }
    }
  }

}

template <typename RealType>
void fast_sweeping3D ( RealType * __restrict__ grid, int nx, int ny, int nz, RealType h, int PBC );

//this function compares the solution of two different fast sweeping implementations and
//measures the speedup
void unit_fast_sweeping3D ( int nx, int ny, int nz, int num_points );

}

#endif // __FASTSWEEPING_H
