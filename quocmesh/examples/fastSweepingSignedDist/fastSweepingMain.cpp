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

#include <fastSweeping.h>
#include <smallVec.h>

#define FREQUENCY 2300000000
#define get_ticks(var) {             \
      unsigned int __a, __d;             \
      asm volatile("rdtsc" : "=a" (__a), "=d" (__d));      \
      var = (static_cast<unsigned long> ( __a )) | ((static_cast<unsigned long> ( __d )) << 32); \
   } while(0)



typedef double RealType;

int main ( int /*argc*/, char** /*argv*/ ) {
  int nx = 15;
  int ny = 15;
  int nz = 15;
  int num_points = 15;

  int pbc_flag = 1;
  printf ( "grid size: %d %d %d\n", nx, ny, nz );
  double h = 1.0;
  printf ( "grid spacing: %f \n", h );
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
  unsigned long start, stop;
  get_ticks ( start );
  qc::fast_sweeping3D ( grid_copy, nx, ny, nz, h, pbc_flag );
  get_ticks ( stop );
  double time3d = ( stop - start );
  printf ( "time: %f\n", time3d );

  /**************************************
   * compute distances
   *************************************/
  get_ticks ( start );
  qc::fast_sweeping_queue_PBC ( grid, nx, ny, nz, h );
  get_ticks ( stop );
  double time_queue = ( stop - start );
  printf ( "time: %f\n", time_queue );
  printf ( "Speedup: %f\n", time3d / time_queue );




  /**************************************
   * Compare computed distances to
   * exact solution
   *************************************/
  double max_difference = 0.0;
  aol::Vec<6, int> idx_max_difference;
  double max_under_est_difference = 0.0;
  aol::Vec<6, int> idx_max_under_est_difference;

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
  delete[] grid_copy;
}
