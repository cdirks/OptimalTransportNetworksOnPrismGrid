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

/**
 *  \file
 * \brief generating 3D signed distance functions
 *
 * Example for generating signed distance functions to different interfaces
 * (plane, sphere, point) in 3d by using qc::DistanceSweeper3d.
 * This example requires knowing the exact distances for initialization which is
 * generally unknown, use SignedDistanceSweepingOp3D if only a level set function
 * is known.
 *
 * \author Ole Schwen, Martina Teusner
*/

#include<sweeping.h>
#include<quoc.h>
#include<scalarArray.h>

int main ( int, char** ) {
  try {

    qc::GridDefinition grid ( 4, qc::QC_3D );

    {
      cerr << "Treating planar aligned interface " << endl;
      qc::ScalarArray<double, qc::QC_3D> sdf ( grid );        // to compute the signed distance function
      qc::ScalarArray<double, qc::QC_3D> distances ( grid );  // true distances as reference solution
      qc::BitArray<qc::QC_3D> seed                            // will mark which points are used as starting points
        ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

      sdf.setAll ( aol::NumberTrait<double>::Inf );  // first, we set all values to plus infinity

      const double position = 1.0 / 3.0;
      const double bandwidth = grid.H();

      for ( int i = 0; i < grid.getWidth(); ++i ) {
        const double x = i * grid.H();
        for ( int j = 0; j < grid.getWidth(); ++j ) {
          for ( int k = 0; k < grid.getWidth(); ++k ) {
            if ( aol::Abs ( x - position ) < bandwidth ) {
              sdf.set ( i, j, k, x - position );     // next, we initialize points next to the zero interface for which we wish to determine the signed distance function
              seed.set ( i, j, k, true );            // and declare those points as seed points
            }
            distances.set ( i, j, k, x - position ); // and we also compute the reference solution ...
          }
        }
      }

      // Now, we set up a DistanceSweeper. It needs to know the seed mask and the underlying grid.
      qc::DistanceSweeper3d<double> distanceSweeper ( seed, grid );
      distanceSweeper.setVerboseMode ( true );

      // For the computation, we need initial values at the seed points and plus infinity values at all other points.
      distanceSweeper.computeDistances ( sdf );

      // Finally, compare the sdf computed by sweeping to the reference solution.
      qc::ScalarArray<double, qc::QC_3D> difference ( grid );
      difference += sdf;
      difference -= distances;

      cerr << "approximate L2 norm of difference = " << difference.norm() / difference.size() << ", Linf norm = " << difference.getMaxAbsValue() << endl;
      //       difference.saveRaw ( "difference_plane.raw", qc::PGM_UNSIGNED_CHAR_BINARY, difference.getMaxAbsValue() );
    }

    {
      cerr << "Treating spherical interface " << endl;

      qc::ScalarArray<double, qc::QC_3D> sdf ( grid );        // to compute the signed distance function
      qc::ScalarArray<double, qc::QC_3D> distances ( grid );  // true distances as reference solution
      qc::BitArray<qc::QC_3D> seed                            // will mark which points are used as starting points
        ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

      sdf.setAll ( aol::NumberTrait<double>::Inf );  // first, we set all values to plus infinity

      const double radius = 1.0 / 3.0;
      const double bandwidth = sqrt ( 3. ) * grid.H();

      for ( int i = 0; i < grid.getWidth(); ++i ) {
        for ( int j = 0; j < grid.getWidth(); ++j ) {
          for ( int k = 0; k < grid.getWidth(); ++k ) {
            aol::Vec3<double> pos ( i, j, k );
            pos *= grid.H();
            if ( aol::Abs ( pos.norm() - radius ) < bandwidth ) {
              sdf.set ( i, j, k, pos.norm() - radius );     // next, we initialize points next to the zero interface for which we wish to determine the signed distance function
              seed.set ( i, j, k, true );            // and declare those points as seed points
            }
            distances.set ( i, j, k, pos.norm() - radius ); // and we also compute the reference solution ...
          }
        }
      }

      // Now, we set up a DistanceSweeper. It needs to know the seed mask and the underlying grid.
      qc::DistanceSweeper3d<double> distanceSweeper ( seed, grid );
      distanceSweeper.setVerboseMode ( true );

      // For the computation, we need initial values at the seed points and plus infinity values at all other points.
      distanceSweeper.computeDistances ( sdf );

      // Finally, compare the sdf computed by sweeping to the reference solution.
      qc::ScalarArray<double, qc::QC_3D> difference ( grid );
      difference += sdf;
      difference -= distances;

      cerr << "approximate L2 norm of difference = " << difference.norm() / difference.size() << ", Linf norm = " << difference.getMaxAbsValue() << endl;
      //       difference.saveRaw ( "difference_sphere.raw", qc::PGM_UNSIGNED_CHAR_BINARY, difference.getMaxAbsValue() );
    }

    {
      cerr << "Treating point spherical interface " << endl;

      qc::ScalarArray<double, qc::QC_3D> sdf ( grid );        // to compute the signed distance function
      qc::ScalarArray<double, qc::QC_3D> distances ( grid );  // true distances as reference solution
      qc::BitArray<qc::QC_3D> seed                            // will mark which points are used as starting points
        ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

      sdf.setAll ( aol::NumberTrait<double>::Inf );  // first, we set all values to plus infinity

      const int ctr = grid.getWidth() / 2;

      for ( int i = 0; i < grid.getWidth(); ++i ) {
        for ( int j = 0; j < grid.getWidth(); ++j ) {
          for ( int k = 0; k < grid.getWidth(); ++k ) {
            aol::Vec3<double> pos ( i, j, k ), rad ( i - ctr, j - ctr, k - ctr );
            double radius = rad.norm() * grid.H();
            distances.set ( i, j, k, radius );
          }
        }
      }
      sdf.set ( ctr, ctr, ctr, 0.0 );
      seed.set ( ctr, ctr, ctr, true );

      // Now, we set up a DistanceSweeper. It needs to know the seed mask and the underlying grid.
      qc::DistanceSweeper3d<double> distanceSweeper ( seed, grid );
      distanceSweeper.setVerboseMode ( true );

      // For the computation, we need initial values at the seed points and plus infinity values at all other points.
      distanceSweeper.computeDistances ( sdf );

      // Finally, compare the sdf computed by sweeping to the reference solution.
      qc::ScalarArray<double, qc::QC_3D> difference ( grid );
      difference += sdf;
      difference -= distances;

      cerr << "approximate L2 norm of difference = " << difference.norm() / difference.size() << ", Linf norm = " << difference.getMaxAbsValue() << endl;
      //       difference.saveRaw ( "difference_pure_sphere.raw", qc::PGM_UNSIGNED_CHAR_BINARY, difference.getMaxAbsValue() );
    }


  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
