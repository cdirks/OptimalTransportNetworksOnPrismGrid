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

#include <quoc.h>
#include "mpiUtils.h"
#include "mpiIncludes.h"

namespace qc {

void cart_map ( int reorder, int *procgrid, int *myloc, int procneigh[3][2], int *grid2proc ) {
  int periods[3];
  periods[0] = periods[1] = periods[2] = 1;
  MPI_Comm cartesian;

  MPI_Cart_create ( MPI_COMM_WORLD, 3, procgrid, periods, reorder, &cartesian );
  MPI_Cart_get ( cartesian, 3, procgrid, periods, myloc );
  MPI_Cart_shift ( cartesian, 0, 1, &procneigh[0][0], &procneigh[0][1] );
  MPI_Cart_shift ( cartesian, 1, 1, &procneigh[1][0], &procneigh[1][1] );
  MPI_Cart_shift ( cartesian, 2, 1, &procneigh[2][0], &procneigh[2][1] );

  int coords[3];
  int i, j, k;
  for ( i = 0; i < procgrid[0]; i++ )
    for ( j = 0; j < procgrid[1]; j++ )
      for ( k = 0; k < procgrid[2]; k++ ) {
        coords[0] = i; coords[1] = j; coords[2] = k;
        MPI_Cart_rank ( cartesian, coords, &grid2proc[i * procgrid[2] * procgrid[1] + j * procgrid[2] + k] );
      }

  MPI_Comm_free ( &cartesian );
}

void scatter_data ( const double* bins, int num_bins[3], double* bins_local, int num_bins_local[3], int my_bins_start[3] ) {
  for ( int k = 1; k < num_bins_local[2] - 1; ++k )
    for ( int j = 1; j < num_bins_local[1] - 1; ++j )
      for ( int i = 1; i < num_bins_local[0] - 1; ++i ) {
        bins_local[ ILexCombine3 ( i, j, k, num_bins_local[0], num_bins_local[1] ) ] = bins[ ILexCombine3 ( my_bins_start[0] +  i - 1, ( my_bins_start[1] + j - 1 ), my_bins_start[2] + k - 1, num_bins[0], num_bins[1] ) ];
      }
}

void insert_data ( double* bins, int num_bins[3], const double* buffer, int bin_range[6] ) {
  for ( int k = bin_range[4], k_local = 0; k < bin_range[5]; ++k, ++k_local )
    for ( int j = bin_range[2], j_local = 0; j < bin_range[3]; ++j, ++j_local )
      for ( int i = bin_range[0], i_local = 0; i < bin_range[1]; ++i, ++i_local )
        bins[ ILexCombine3 ( i, j, k, num_bins[0], num_bins[1] ) ] = buffer[ ILexCombine3 ( i_local, j_local, k_local, bin_range[1] - bin_range[0], bin_range[3] - bin_range[2] ) ];
}

void gather_data ( double* bins, int num_bins[3], const double* bins_local, int my_bins_start[3], int my_bins_end[3], int root ) {
  int rank, mpi_size;
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &mpi_size );

  int bin_range[6] = {my_bins_start[0], my_bins_end[0], my_bins_start[1], my_bins_end[1], my_bins_start[2], my_bins_end[2] };
  int numX = my_bins_end[0] - my_bins_start[0] ;
  int numY = my_bins_end[1] - my_bins_start[1] ;
  int numZ = my_bins_end[2] - my_bins_start[2] ;
  int num_elements = numX * numY * numZ;

  const int capacity = 32 * num_elements;
  double *buffer = ( double* ) malloc ( sizeof ( double ) * capacity );

  //put elements into a contigious buffer (i.e. remove ghosts)
  int counter = 0;

  for ( int k = 0; k < numZ; ++k ) {
    for ( int j = 0; j < numY; ++j ) {
      for ( int i = 0; i < numX; ++i ) {
        buffer[counter++] = bins_local[ ILexCombine3 ( i + 1, j + 1, k + 1, numX + 2, numY + 2 ) ];
      }
    }
  }

  int received_messages = 1;

  if ( rank == root ) {
    //insert local data first
    insert_data ( bins, num_bins, buffer, bin_range );

    //insert remote data
    while ( received_messages != mpi_size ) {
      MPI_Status status;
      MPI_Recv ( bin_range, 6, MPI_INT, MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, &status );


      num_elements = ( bin_range[1] - bin_range[0] ) * ( bin_range[3] - bin_range[2] ) * ( bin_range[5] - bin_range[4] );

      MPI_Recv ( buffer, num_elements, MPI_DOUBLE, status.MPI_SOURCE, 10, MPI_COMM_WORLD, &status );
      insert_data ( bins, num_bins, buffer, bin_range );
      received_messages++;
    }

  } else {
    MPI_Send ( bin_range, 6, MPI_INT, 0, 100, MPI_COMM_WORLD );
    MPI_Send ( buffer, num_elements, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD );
  }

  free ( buffer );
}

}
