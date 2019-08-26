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
#include <mpiIncludes.h>

#include "ghostExchange.h"


namespace qc {

/**
 * \param buffer must be of size num_bins_local[] * num_bins_local[]
 */
void pack_buffer ( const double *bins_local, const int* num_bins_local, double* buffer, Face face ) {
  int counter = 0;
  int numX = num_bins_local[0];
  int numY = num_bins_local[1];
  int numZ = num_bins_local[2];
  if ( face == GH_FRONT ) {
    for ( int i = 0; i < numY; ++i ) {
      for ( int j = 0; j < numX; ++j ) {
        buffer[counter++] = bins_local[ ILexCombine3 ( j, i, 1, numX, numY ) ];
      }
    }
  } else if ( face == GH_BACK ) {
    for ( int i = 0; i < numY; ++i ) {
      for ( int j = 0; j < numX; ++j ) {
        buffer[counter++] = bins_local[ ILexCombine3 ( j, i, numZ - 2, numX, numY ) ];
      }
    }
  } else if ( face == GH_BUTTOM ) {
    for ( int i = 0; i < numZ; ++i ) {
      for ( int j = 0; j < numX; ++j ) {
        buffer[counter++] = bins_local[ ILexCombine3 ( j, 1, i, numX, numY ) ] ;
      }
    }
  } else if ( face == GH_TOP ) {
    for ( int i = 0; i < numZ; ++i ) {
      for ( int j = 0; j < numX; ++j ) {
        buffer[counter++] = bins_local[ ILexCombine3 ( j, numY - 2, i, numX, numY ) ] ;
      }
    }
  } else if ( face == GH_LEFT ) {
    for ( int i = 0; i < numZ; ++i ) {
      for ( int j = 0; j < numY; ++j ) {
        buffer[counter++] = bins_local[ ILexCombine3 ( 1, j, i, numX, numY ) ] ;
      }
    }
  } else if ( face == GH_RIGHT ) {
    for ( int i = 0; i < numZ; ++i ) {
      for ( int j = 0; j < numY; ++j ) {
        buffer[counter++] = bins_local[ ILexCombine3 ( numX - 2, j, i, numX, numY ) ] ;
      }
    }
  }
}

void unpack_buffer ( double *bins_local, const int* num_bins_local, const double* buffer, Face face ) {
  int counter = 0;
  int numX = num_bins_local[0];
  int numY = num_bins_local[1];
  int numZ = num_bins_local[2];

  if ( face == GH_FRONT ) {
    for ( int i = 0; i < numY; ++i ) {
      for ( int j = 0; j < numX; ++j ) {
        bins_local[ ILexCombine3 ( j, i, 0, numX, numY ) ] = buffer[counter++];
      }
    }
  } else if ( face == GH_BACK ) {
    for ( int i = 0; i < numY; ++i ) {
      for ( int j = 0; j < numX; ++j ) {
        bins_local[ ILexCombine3 ( j, i, numZ - 1, numX, numY ) ] = buffer[counter++];
      }
    }
  } else if ( face == GH_BUTTOM ) {
    for ( int i = 0; i < numZ; ++i ) {
      for ( int j = 0; j < numX; ++j ) {
        bins_local[ ILexCombine3 ( j, 0, i, numX, numY ) ]  = buffer[counter++];
      }
    }
  } else if ( face == GH_TOP ) {
    for ( int i = 0; i < numZ; ++i ) {
      for ( int j = 0; j < numX; ++j ) {
        bins_local[ ILexCombine3 ( j, numY - 1, i, numX, numY ) ]  = buffer[counter++];
      }
    }
  } else if ( face == GH_LEFT ) {
    for ( int i = 0; i < numZ; ++i ) {
      for ( int j = 0; j < numY; ++j ) {
        bins_local[ ILexCombine3 ( 0, j, i, numX, numY ) ]  = buffer[counter++];
      }
    }
  } else if ( face == GH_RIGHT ) {
    for ( int i = 0; i < numZ; ++i ) {
      for ( int j = 0; j < numY; ++j ) {
        bins_local[ ILexCombine3 ( numX - 1, j, i, numX, numY ) ]  = buffer[counter++];
      }
    }
  }
}


void exchange_faces ( double* bins_local, int* num_bins_local, int count, int nbrL, int nbrR, Face faceL, Face faceR ) {
  double* send_buffer_left   = ( double* ) malloc ( sizeof ( double ) * count );
  double* recv_buffer_left   = ( double* ) malloc ( sizeof ( double ) * count );
  double* send_buffer_right  = ( double* ) malloc ( sizeof ( double ) * count );
  double* recv_buffer_right  = ( double* ) malloc ( sizeof ( double ) * count );

  pack_buffer ( bins_local, num_bins_local, send_buffer_left, faceL );
  pack_buffer ( bins_local, num_bins_local, send_buffer_right, faceR );

  //send to left and right neighbor
  MPI_Request req[4];
  MPI_Isend ( send_buffer_left, count, MPI_DOUBLE, nbrL, 1, MPI_COMM_WORLD, &req[0] );
  MPI_Isend ( send_buffer_right, count, MPI_DOUBLE, nbrR, 0, MPI_COMM_WORLD, &req[1] );

  //receive from left and right neighbor
  MPI_Irecv ( recv_buffer_left, count, MPI_DOUBLE, nbrL, 0, MPI_COMM_WORLD, &req[2] );
  MPI_Irecv ( recv_buffer_right, count, MPI_DOUBLE, nbrR, 1, MPI_COMM_WORLD, &req[3] );

  //block till all the communication is done
  MPI_Status status[4];
  MPI_Waitall ( 4, req, status );

  unpack_buffer ( bins_local, num_bins_local, recv_buffer_left, faceL );
  unpack_buffer ( bins_local, num_bins_local, recv_buffer_right, faceR );

  MPI_Barrier ( MPI_COMM_WORLD );

  free ( send_buffer_left );
  free ( send_buffer_right );
  free ( recv_buffer_left );
  free ( recv_buffer_right );
}

void ghost_exchange ( double *bins_local, int num_bins_local[3], const int ( &procneigh ) [3][2] ) {
  /********************************************
   * Exchange with left & right neighbor
   ********************************************/
  exchange_faces ( bins_local, num_bins_local, num_bins_local[2] * num_bins_local[1], procneigh[0][0], procneigh[0][1], GH_LEFT, GH_RIGHT );

  /********************************************
   * Exchange with top & buttom neighbor
   ********************************************/
  exchange_faces ( bins_local, num_bins_local, num_bins_local[2] * num_bins_local[0], procneigh[1][0], procneigh[1][1], GH_BUTTOM, GH_TOP );

  /********************************************
   * Exchange with front & back neighbor
   ********************************************/
  exchange_faces ( bins_local, num_bins_local, num_bins_local[0] * num_bins_local[1], procneigh[2][0], procneigh[2][1], GH_FRONT, GH_BACK );
}

}

