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

#include <smallVec.h>
#include <fastSweeping.h>

namespace qc {

template <typename RealType>
void sort3 ( const RealType a, const RealType b, const RealType c, RealType* sortedCoeff ) {
  if ( a < b ) {
    if ( a < c ) {
      sortedCoeff[0] = a;
      if ( b < c ) {
        sortedCoeff[1] = b;
        sortedCoeff[2] = c;
      } else {
        sortedCoeff[1] = c;
        sortedCoeff[2] = b;
      }
    } else {
      sortedCoeff[0] = c;
      if ( a < b ) {
        sortedCoeff[1] = a;
        sortedCoeff[2] = b;
      } else {
        sortedCoeff[1] = b;
        sortedCoeff[2] = a;
      }
    }
  } else {
    if ( b < c ) {
      sortedCoeff[0] = b;
      if ( a < c ) {
        sortedCoeff[1] = a;
        sortedCoeff[2] = c;
      } else {
        sortedCoeff[1] = c;
        sortedCoeff[2] = a;
      }
    } else {
      sortedCoeff[0] = c;
      if ( a < b ) {
        sortedCoeff[1] = a;
        sortedCoeff[2] = b;
      } else {
        sortedCoeff[1] = b;
        sortedCoeff[2] = a;
      }
    }
  }
}

/**
 * Solves ((x - a)^+)^2 + ((x - b)^+)^2 + ((x - c)^+)^2 = h*h for x. With (y)^+ := (y>0) ? y : 0
 */
template <typename RealType>
RealType solve ( const RealType a, const RealType b, const RealType c, const RealType h ) {

  /*************************
   * Sort coefficients
   *************************/
  RealType sortedCoeff[3];
  sort3 ( a, b, c, sortedCoeff );

  RealType solution;

  //solution for 1D
  solution = sortedCoeff[0] + h;
  if ( solution <= sortedCoeff[1] ) {
    return solution;
  }

  //solution for 2D
  solution = ( sortedCoeff[0] + sortedCoeff[1] + sqrt ( 2 * h * h - ( sortedCoeff[0] - sortedCoeff[1] ) * ( sortedCoeff[0] - sortedCoeff[1] ) ) ) * 0.5;
  if ( solution <= sortedCoeff[2] ) {
    return solution;
  }

  //solution for 3D
  RealType tmp = ( sortedCoeff[0] + sortedCoeff[1] + sortedCoeff[2] );
  RealType tmp2 = tmp * tmp - 3.0 * ( sortedCoeff[0] * sortedCoeff[0] + sortedCoeff[1] * sortedCoeff[1] + sortedCoeff[2] * sortedCoeff[2] - h * h );
  QUOC_ASSERT ( tmp2 >= 0 );
  return 1.0 / 3.0 * ( sqrt ( tmp2 ) + tmp );

}

template <typename RealType>
void print_3d_array ( RealType* grid, int nx, int ny, int nz ) {
  for ( int i = 0; i < nz; ++i ) {
    for ( int j = 0; j < ny; ++j ) {
      for ( int k = 0; k < nx; ++k ) {
        int idx = i * nx * ny + j * nx + k;
        printf ( "%f ", grid[idx] );
      }
      printf ( "\n" );
    }
    printf ( "\n" );
  }
}

template <typename RealType>
void fast_sweeping3D ( RealType * __restrict__ grid, int nx, int ny, int nz, RealType h, int PBC ) {
  sweep3D<RealType, 0, 0, 0> ( grid, nx, ny, nz, h, PBC );
  sweep3D<RealType, 1, 0, 0> ( grid, nx, ny, nz, h, PBC );
  sweep3D<RealType, 0, 1, 0> ( grid, nx, ny, nz, h, PBC );
  sweep3D<RealType, 1, 1, 0> ( grid, nx, ny, nz, h, PBC );
  sweep3D<RealType, 0, 0, 1> ( grid, nx, ny, nz, h, PBC );
  sweep3D<RealType, 1, 0, 1> ( grid, nx, ny, nz, h, PBC );
  sweep3D<RealType, 0, 1, 1> ( grid, nx, ny, nz, h, PBC );
  sweep3D<RealType, 1, 1, 1> ( grid, nx, ny, nz, h, PBC );
}

template void fast_sweeping3D<double> ( double* __restrict__ grid, int nx, int ny, int nz, double h, int PBC );

void add_nbr_to_queue ( int idx_x, int idx_y, int idx_z, double* grid, int nx, int ny, int nz, Queue *queue, double new_value ) {
  //add front & back nbr
  for ( int i = -1; i <= 1; ++i ) {
    if ( idx_z + i == 0 || idx_z + i == nz - 1 ) continue; //don't add ghost points

    int idx = ILexCombine3 ( idx_x, idx_y, idx_z + i, nx, ny );
    if ( grid[ idx ] > new_value )
      queue->add_point ( idx_x, idx_y, idx_z + i );
  }

  //add top & bottom nbr
  for ( int j = -1; j <= 1; ++j ) {
    if ( idx_y + j == 0 || idx_y + j == ny - 1 ) continue;

    int idx = ILexCombine3 ( idx_x, idx_y + j, idx_z, nx, ny );
    if ( grid[ idx ] > new_value )
      queue->add_point ( idx_x, idx_y + j, idx_z );
  }

  //add left& right nbr
  for ( int k = -1; k <= 1; ++k ) {
    if ( idx_x + k == 0 || idx_x + k == nx - 1 ) continue;

    int idx = ILexCombine3 ( idx_x + k, idx_y, idx_z, nx, ny );
    if ( grid[ idx ] > new_value )
      queue->add_point ( idx_x + k, idx_y, idx_z );
  }
}

void add_nbr_to_queue_PBC ( int idx_x, int idx_y, int idx_z, double* grid, int nx, int ny, int nz, Queue *queue, double new_value ) {
  //add front & back nbr
  for ( int i = -1; i <= 1; ++i ) {
    int tmp_z = idx_z + i;
    if ( tmp_z == -1 )
      tmp_z = nz - 1;
    else if ( tmp_z == nz )
      tmp_z = 0;

    int idx = ILexCombine3 ( idx_x, idx_y, tmp_z, nx, ny );
    if ( grid[ idx ] > new_value )
      queue->add_point ( idx_x, idx_y, tmp_z );
  }

  //add top & bottom nbr
  for ( int j = -1; j <= 1; ++j ) {
    int tmp_y = idx_y + j;

    if ( tmp_y == -1 )
      tmp_y = ny - 1;
    else if ( tmp_y == ny )
      tmp_y = 0;

    int idx = ILexCombine3 ( idx_x, tmp_y , idx_z, nx, ny );
    if ( grid[ idx ] > new_value )
      queue->add_point ( idx_x, tmp_y, idx_z );
  }

  //add left& right nbr
  for ( int k = -1; k <= 1; ++k ) {
    int tmp_x = idx_x + k;

    if ( tmp_x == -1 )
      tmp_x = nx - 1;
    if ( tmp_x == nx )
      tmp_x = 0;

    int idx = ILexCombine3 ( tmp_x, idx_y, idx_z, nx, ny );
    if ( grid[ idx ] > new_value )
      queue->add_point ( tmp_x, idx_y, idx_z );
  }
}

void fast_sweeping_queue_PBC ( double* grid, int nx, int ny, int nz, double h ) {
  size_t size = nx * ny * nz * sizeof ( double );

  // grid_old <- grid
  double* grid_old = new double[nx * ny * nz];
  memcpy ( grid_old, grid, size );

  Queue* queue = new Queue ( nx, ny, nz ) ;
  Queue* queue_old = new Queue ( nx, ny, nz ) ;

  int max_iter = max ( max ( nx, ny ), nz ) * 10;
  double epsilon = h * 0.01;

  /**************************************
   * Add interfacial grid points to queue
   **************************************/
  for ( int i = 0; i < nz; ++i ) {
    for ( int j = 0; j < ny; ++j ) {
      for ( int k = 0; k < nx; ++k ) {

        double tmp_value = grid[ ILexCombine3 ( k, j, i, nx, ny ) ];
        if ( tmp_value == 0.0 ) continue;

        //loop over nbrs
        if ( grid[ ILexCombine3 ( ( k + 1 ) % nx, j, i, nx, ny ) ] != tmp_value ) {
          queue->add_point ( k, j, i );
          continue;
        }
        if ( grid[ ILexCombine3 ( ( k - 1 + nx ) % nx, j, i, nx, ny ) ] != tmp_value ) {
          queue->add_point ( k, j, i );
          continue;
        }
        if ( grid[ ILexCombine3 ( k, ( j + 1 ) % ny, i, nx, ny ) ] != tmp_value ) {
          queue->add_point ( k, j, i );
          continue;
        }
        if ( grid[ ILexCombine3 ( k, ( j - 1 + ny ) % ny, i, nx, ny ) ] != tmp_value ) {
          queue->add_point ( k, j, i );
          continue;
        }
        if ( grid[ ILexCombine3 ( k, j, ( i + 1 ) % nz, nx, ny ) ] != tmp_value ) {
          queue->add_point ( k, j, i );
          continue;
        }
        if ( grid[ ILexCombine3 ( k, j, ( i - 1 + nz ) % nz, nx, ny ) ] != tmp_value ) {
          queue->add_point ( k, j, i );
          continue;
        }
      }
    }
  }

  int nxy = nx * ny;
  for ( int iter = 0; iter < max_iter; ++iter ) {

    double max_change = 0.0;

    for ( int ii = 0; ii < queue->size; ++ii ) {
      int idx_x = queue->list[3 * ii];
      int idx_y = queue->list[3 * ii + 1];
      int idx_z = queue->list[3 * ii + 2];
      int idx = ILexCombine3 ( idx_x, idx_y, idx_z, nx, ny );

      double a = min ( grid[ ( idx_z - 1 + nz ) % nz * nxy + idx_y           * nx + idx_x           ], grid[ ( idx_z + 1 ) % nz * nxy + idx_y       * nx + idx_x      ] );
      double b = min ( grid[ idx_z           * nxy + ( idx_y - 1 + ny ) % ny * nx + idx_x           ], grid[ idx_z     * nxy + ( idx_y + 1 ) % ny   * nx + idx_x      ] );
      double c = min ( grid[ idx_z           * nxy + idx_y           * nx + ( idx_x - 1 + nx ) % nx ], grid[ idx_z     * nxy + idx_y       * nx + ( idx_x + 1 ) % nx] );

      double new_value = solve ( a, b, c, h );

      if ( new_value < grid_old[ idx ] )
        grid_old[ idx ] = new_value;

      double change = fabs ( new_value - grid[ idx ] );
      if ( change > max_change )
        max_change = change;

      //add nbr grid-points to queue
      add_nbr_to_queue_PBC ( idx_x, idx_y, idx_z, grid, nx, ny, nz, queue_old, new_value );

    }

    //swap queues
    Queue* tmp = queue;
    queue      = queue_old;
    queue_old  = tmp;
    queue_old->reset();

    //grid <- grid_old
    memcpy ( grid, grid_old, size );

    if ( max_change < epsilon )
      break;
  }

  delete[] grid_old;
  delete queue_old;
  delete queue;
}

template <typename RealType>
void unit_fast_sweeping3D ( int nx, int ny, int nz, int num_points ) {
  int pbc_flag = 1;
  printf ( "grid size: %d %d %d\n", nx, ny, nz );
  double h = 1.0;
  int nxy = nx * ny;

  RealType* grid = new RealType[nx * ny * nz];
  RealType* grid_copy = new RealType[nx * ny * nz];

  for ( int i = 0; i < nz; ++i ) {
    for ( int j = 0; j < ny; ++j ) {
      for ( int k = 0; k < nx; ++k ) {
        int idx = i * nxy + j * nx + k;
        grid[idx] = nx * ny * nz * h;
        grid_copy[idx] = nx * ny * nz * h;
      }
    }
  }

  /**************************************
   * generate random points
   *************************************/
  int seed = time ( NULL );
  srand ( seed );
  printf ( "seed: %d\n", seed );

  int* points = new int[num_points * 3];
  for ( int i = 0; i < num_points; ++i ) {
    int x = rand() % ( nx - 2 ) + 1;
    int y = rand() % ( ny - 2 ) + 1;
    int z = rand() % ( nz - 2 ) + 1;
    points[i * 3] = x;
    points[i * 3 + 1] = y;
    points[i * 3 + 2] = z;
    int idx = z * nxy + y * nx + x;
    grid[idx] = 0.0;
    grid_copy[idx] = 0.0;
  }

  /**************************************
   * compute distances
   *************************************/
//   double start = omp_get_wtime();
  fast_sweeping3D ( grid, nx, ny, nz, h, pbc_flag );
//   printf("FSM took: %fms\n",(omp_get_wtime() - start) * 1000);

  /**************************************
   * compute distances
   *************************************/
//   start = omp_get_wtime();
  fast_sweeping_queue_PBC ( grid_copy, nx, ny, nz, h );
//   printf("FSM took: %fms\n",(omp_get_wtime() - start) * 1000);
//   printf("Speedup: %f\n",(omp_get_wtime() - start) * 1000);

  double max_difference = 0.0;
  aol::Vec<6, int> idx_max_difference;
  double max_under_est_difference = 0.0;
  aol::Vec<6, int> idx_max_under_est_difference;


  /**************************************
   * Compare computed distances to
   * exact solution
   *************************************/
  if ( pbc_flag ) {
    for ( int i = 0; i < nz; ++i ) {
      for ( int j = 0; j < ny; ++j ) {
        for ( int k = 0; k < nx; ++k ) {
          double min_distance = 100000000;
          aol::Vec3<int> idx_min;

          //find nearest point on surface from this grid point
          for ( int l = 0; l < num_points; ++l ) {

            //loop over all periodic images of this point
            for ( int nnx = -1; nnx <= 1; nnx ++ ) {
              for ( int nny = -1; nny <= 1; nny ++ ) {
                for ( int nnz = -1; nnz <= 1; nnz ++ ) {

                  double dx = h * ( points[l * 3] + nnx * nx - k );
                  double dy = h * ( points[l * 3 + 1] + nny * ny - j );
                  double dz = h * ( points[l * 3 + 2] + nnz * nz - i );
                  double distance = sqrt ( dx * dx + dy * dy + dz * dz );
                  if ( distance < min_distance ) {
                    min_distance = distance;
                    idx_min[0] = points[l * 3];
                    idx_min[1] = points[l * 3 + 1];
                    idx_min[2] = points[l * 3 + 2];
                  }
                }
              }
            }
          }

          //See A FAST SWEEPING METHOD FOR EIKONAL EQUATIONS by HONGKAI ZHAO for error analysis
          //Theorem 4.3: uh(x, Γ) ≤ uh(x, Γ) ≤ d(x, Γ) + O(|h log h|).
          double difference = ( grid[i * nxy + j * nx + k] - min_distance );
          if ( difference > max_difference ) {
            max_difference = difference;
            idx_max_difference[0] = idx_min[0];
            idx_max_difference[1] = idx_min[1];
            idx_max_difference[2] = idx_min[2];
            idx_max_difference[3] = k;
            idx_max_difference[4] = j;
            idx_max_difference[5] = i;
          }
          //               printf("Computed distance and real distance differ too much. Real: %f, computed: %f\n",min_distance,grid[i * nxy + j * ny + k]);
          if ( grid[i * nxy + j * nx + k] <  min_distance ) {
            double under_est_difference = min_distance - grid[i * nxy + j * nx + k] ;
            if ( max_under_est_difference < under_est_difference ) {
              max_under_est_difference = under_est_difference;
              idx_max_under_est_difference[0] = idx_min[0];
              idx_max_under_est_difference[1] = idx_min[1];
              idx_max_under_est_difference[2] = idx_min[2];
              idx_max_under_est_difference[3] = k;
              idx_max_under_est_difference[4] = j;
              idx_max_under_est_difference[5] = i;
            }
            //               printf("Error: Computed distane is too small.(%d %d %d) (%d %d %d) %f %f\n",k,j,i,idx_min[0],idx_min[1],idx_min[2],min_distance,grid[i * nxy + j * ny + k]);
          }

        }
      }
    }
  }

  if ( max_difference > 0 )
    printf ( "Max difference between surface (%d %d %d) and grid (%d %d %d) over estimated by %f (computed value:%f)\n",
             idx_max_difference[0],
             idx_max_difference[1],
             idx_max_difference[2],
             idx_max_difference[3],
             idx_max_difference[4],
             idx_max_difference[5],
             max_difference,
             grid[idx_max_difference[5] * nxy + idx_max_difference[4] * nx + idx_max_difference[3]] );

  if ( max_under_est_difference > 0 )
    printf ( "Max difference between surface (%d %d %d) and grid (%d %d %d) under estimated by %f (computed value:%f)\n",
             idx_max_under_est_difference[0],
             idx_max_under_est_difference[1],
             idx_max_under_est_difference[2],
             idx_max_under_est_difference[3],
             idx_max_under_est_difference[4],
             idx_max_under_est_difference[5],
             max_under_est_difference,
             grid[idx_max_under_est_difference[5] * nxy + idx_max_under_est_difference[4] * nx + idx_max_under_est_difference[3]] );

  delete[] grid;
}

}
