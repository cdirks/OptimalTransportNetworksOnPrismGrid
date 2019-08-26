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

#include <scalarArray.h>
#include <firstOrderTVAlgos.h>

#include <mpiUtils.h>
#include "mpiMSSeg.h"

void create_pgm_image ( const char* filename, double *data, int size[3] );

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_3D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType, DimensionChoice, 3> > ConfType;

int main ( int argc, char **argv ) {

  if ( argc < 4 ) {
    printf ( "Usage: %s <num-procs-x> <num-procs-y> <num-procs-z>\n", argv[0] );
    exit ( -1 );
  }
  /***************************************
  * Initialize image (create random input)
  *    - Draw a sphere
  *    - Add noise
  ***************************************/

  int num_bins[3] = {150, 150, 150};
  double bin_extent = 2.0 / num_bins[0];
  double noise_level = 1.2;

  double* data = ( double* ) malloc ( sizeof ( double ) * num_bins[0] * num_bins[1] * num_bins[2] );
  for ( int i = 0; i < num_bins[2]; ++i ) {
    double dz = ( i - num_bins[2] * 0.5 ) * bin_extent * 1.2;
    for ( int j = 0; j < num_bins[1]; ++j ) {
      double dy = ( j - num_bins[1] * 0.5 ) * bin_extent * 1.2;
      for ( int k = 0; k < num_bins[0]; ++k ) {
        double dx = ( k - num_bins[0] * 0.5 ) * bin_extent * 1.2;
        double distance_from_center = sqrt ( dx * dx + dy * dy + dz * dz );
        if ( distance_from_center <= 1.0 )
          data[i * num_bins[1]*num_bins[0] + j * num_bins[0] + k] = 1.0 + noise_level * ( ( ( double ) rand() ) / RAND_MAX );
        else
          data[i * num_bins[1]*num_bins[0] + j * num_bins[0] + k] = noise_level * ( ( ( double ) rand() ) / RAND_MAX );
      }
    }
  }
  data[20 * num_bins[1]*num_bins[0] + 20 * num_bins[0] + 20] = 1.0 + noise_level * ( ( ( double ) rand() ) / RAND_MAX );
  data[20 * num_bins[1]*num_bins[0] + 21 * num_bins[0] + 20] = 1.0 + noise_level * ( ( ( double ) rand() ) / RAND_MAX );
  data[20 * num_bins[1]*num_bins[0] + 20 * num_bins[0] + 21] = 1.0 + noise_level * ( ( ( double ) rand() ) / RAND_MAX );
  data[20 * num_bins[1]*num_bins[0] + 21 * num_bins[0] + 21] = 1.0 + noise_level * ( ( ( double ) rand() ) / RAND_MAX );


  /******************************************
   * Initialize MPI
   *    - Scatter the data across the processes
   ******************************************/
  MPI_Init ( &argc, &argv );

  int procgrid[3]; //number of processors in x-, y- and z-direction
  procgrid[0] = atoi ( argv[1] );
  procgrid[1] = atoi ( argv[2] );
  procgrid[2] = atoi ( argv[3] );

  int my_rank;
  MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );

  if ( my_rank == 0 )
    create_pgm_image ( "initial_image.pgm", data, num_bins );

  //arange processes in a cartesian grid
  int* grid2proc = new int[procgrid[0]* procgrid[1]* procgrid[2]];
  int myloc[3];                     // which proc I am in each dim
  int procneigh[3][2];              // my 6 neighboring procs, 0/1 = left/right
  qc::cart_map ( 1,  procgrid, myloc, procneigh, grid2proc );
  delete[] grid2proc;

  /****************************************************
   * Extract local data from global image
   *    - Domain decomposition
   *    - Each process computes the segmentation of a
   *       subdoain. This requires ghost-cells surrounding
   *       the subdomain of each process.
   ***************************************************/

  //start and end position for each process
  int my_bins_start[3];
  int my_bins_end[3];

  for ( int i = 0; i < 3 ; ++i ) {
    int num_bins_local = num_bins[i] / procgrid[i];
    my_bins_start[i] = myloc[i] * num_bins_local;
    if ( procgrid[i] - 1 == myloc[i] ) //the last process in each dimension has to do all the remaining work
      my_bins_end[i] = num_bins[i];
    else
      my_bins_end[i] = min ( ( myloc[i] + 1 ) * num_bins_local, num_bins[i] );
  }

  int num_bins_local[3];
  num_bins_local[0] = my_bins_end[0] - my_bins_start[0] + 2; //+2 for ghost
  num_bins_local[1] = my_bins_end[1] - my_bins_start[1] + 2;
  num_bins_local[2] = my_bins_end[2] - my_bins_start[2] + 2;
  int num_bins_local_total =  num_bins_local[0] * num_bins_local[1] * num_bins_local[2];

  double* data_local = new double[num_bins_local[0] * num_bins_local[1] * num_bins_local[2]];
  memset ( data_local, 0, sizeof ( double ) * num_bins_local[0] * num_bins_local[1] * num_bins_local[2] );

  qc::scatter_data ( data, num_bins, data_local, num_bins_local, my_bins_start );

  /******************************************
   * Segment Volume
   ******************************************/
  typedef RType RealType;
  typedef ConfType ConfiguratorType;

  const RealType gamma = 0.005 / 0.011250;
  const int imageDimension = 1;

  qc::ScalarArray<RealType, ConfiguratorType::Dim> u0 ( num_bins_local[0], num_bins_local[1], num_bins_local[2], data_local, aol::FLAT_COPY );

  qc::MultiArray<RealType, ConfiguratorType::Dim, imageDimension> u0_m ( u0, aol::FLAT_COPY );

  ConfiguratorType::InitType grid ( qc::GridSize<ConfType::Dim> ( u0.getSize() ) );

  MPIPiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, imageDimension> segmentor ( grid, gamma, u0_m, procneigh );
  qc::ScalarArray<RealType, ConfiguratorType::Dim> temp ( grid );

  segmentor.setMaxIterations ( 1000 );
  segmentor.setStopEpsilon ( 0.001 );
  segmentor.setOuterIterations ( 1 );

  aol::MultiVector<RealType> &initialGrayVaues = segmentor.getMeanValuesReference();
  initialGrayVaues[0][0] = 0.5 * noise_level;
  initialGrayVaues[1][0] = 1.0 + 0.5 * noise_level;

  double *old_segmentation = new double[3 * num_bins_local_total];
  double *ptr[3];
  ptr[0] = old_segmentation;
  ptr[1] = & ( old_segmentation[num_bins_local_total] );
  ptr[2] = & ( old_segmentation[2 * num_bins_local_total] );
  qc::MultiArray<RealType, ConfiguratorType::Dim> pdual ( qc::GridSize<ConfType::Dim> ( num_bins_local[0], num_bins_local[1], num_bins_local[2] ), ptr );

  //call segmentation algorithm
  segmentor.segment ( temp, &pdual );

  for ( int k = 0; k < num_bins_local[0]; ++k )
    for ( int j = 0; j < num_bins_local[1]; ++j )
      for ( int i = 0; i < num_bins_local[2]; ++i )
        if ( temp.get ( k, j, i ) >= 0.5 )
          data_local[ qc::ILexCombine3 ( k, j, i, num_bins_local[0], num_bins_local[1] ) ] = 1;
        else
          data_local[ qc::ILexCombine3 ( k, j, i, num_bins_local[0], num_bins_local[1] ) ] = 0;


  /*************************************
   * Gather the segmented subdomains into a final image
   *************************************/
  int root = 0;
  memset ( data, 0, sizeof ( double ) * num_bins[0] * num_bins[1] * num_bins[2] );
  qc::gather_data ( data, num_bins, data_local, my_bins_start, my_bins_end, root );

  if ( my_rank == 0 )
    create_pgm_image ( "segmented_image.pgm", data, num_bins );

  MPI_Finalize();

  return 0;
}

void create_pgm_image ( const char* filename, double *data, int size[3] ) {
  FILE* pgm = fopen ( filename, "w+" );

  int max = 99;
  fprintf ( pgm, "P2\n" );
  fprintf ( pgm, "%d %d\n", size[2], size[1] );
  fprintf ( pgm, "%d\n", max );

  int idx_x = size[0] / 2;

  double max_value = -10000;
  for ( int i = 0; i < size[1]; ++i ) {
    for ( int j = 0; j < size[2]; ++j ) {
      double value = data[qc::ILexCombine3 ( idx_x, i, j, size[0], size[1] )];
      if ( value > max_value )
        max_value = value;
    }
  }

  for ( int i = 0; i < size[1]; ++i ) {
    for ( int j = 0; j < size[2]; ++j ) {
      int value = ( int ) ( data[qc::ILexCombine3 ( idx_x, i, j, size[0], size[1] )] / max_value * max );
      fprintf ( pgm, "%d ", max - value );
    }
    fprintf ( pgm, "\n" );
  }
  fclose ( pgm );
}
