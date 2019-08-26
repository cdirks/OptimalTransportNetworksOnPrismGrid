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

/** \file
 *  \brief generating 2D signed distance functions
 *
 *   Example for generating signed distance functions to different
 *   interfaces (line, circle, point) in 2d by using
 *   qc::DistanceSweeper2d
 *
 *  \author Franken, Schwen
 */

#include<sweeping.h>
#include<quoc.h>
#include<scalarArray.h>

int main ( int, char** ) {
  try {

    qc::GridDefinition grid ( 6, qc::QC_2D );

    {
      cerr << "Treating a line as interface " << endl;
      qc::ScalarArray<double, qc::QC_2D> sdf ( grid );        // to compute the signed distance function
      qc::ScalarArray<double, qc::QC_2D> distances ( grid );  // true distances as reference solution
      qc::BitArray<qc::QC_2D> seed                            // will mark which points are used as starting points
        ( qc::GridSize<qc::QC_2D>::createFrom ( grid ) );

      sdf.setAll ( aol::NumberTrait<double>::Inf );  // first, we set all values to plus infinity

      const double position = 1.0 / 3.0;
      const double bandwidth = grid.H();

      for ( int i = 0; i < grid.getWidth(); ++i ) {
        const double x = i * grid.H();
        for ( int j = 0; j < grid.getWidth(); ++j ) {
          if ( aol::Abs ( x - position ) < bandwidth ) {
            sdf.set ( i, j, x - position );     // next, we initialize points next to the zero interface for which we wish to determine the signed distance function
            seed.set ( i, j, true );            // and declare those points as seed points
          }
          distances.set ( i, j, x - position ); // and we also compute the reference solution ...
        }
      }

      // Now, we set up a DistanceSweeper. It needs to know the seed mask and the underlying grid.
      qc::DistanceSweeper2d<double> distanceSweeper ( seed, grid );
      distanceSweeper.setVerboseMode ( true );

      // For the computation, we need initial values at the seed points and plus infinity values at all other points.
      distanceSweeper.computeDistances ( sdf );

      // Finally, compare the sdf computed by sweeping to the reference solution.
      qc::ScalarArray<double, qc::QC_2D> difference ( grid );
      difference += sdf;
      difference -= distances;


      cerr << "approximate L2 norm of difference = " << difference.norm() / difference.size() << ", Linf norm = " << difference.getMaxAbsValue() << endl;

#if 0
      double maxDiff = difference.getMaxAbsValue();

      for ( int i = 0; i < difference.size(); ++i ) {
        difference[i] = 255 * aol::Abs ( difference[i] ) / maxDiff;
      }
      difference.save ( "difference_plane.pgm" );
      double maxSdf = sdf.getMaxAbsValue();

      for ( int i = 0; i < sdf.size(); ++i ) {
        sdf[i] = 255 * aol::Abs ( sdf[i] ) / maxSdf;
      }
      sdf.save ( "sdf_plane.pgm" );
#endif
    }

    {
      cerr << "Treating circular interface " << endl;

      qc::ScalarArray<double, qc::QC_2D> sdf ( grid );        // to compute the signed distance function
      qc::ScalarArray<double, qc::QC_2D> distances ( grid );  // true distances as reference solution
      qc::BitArray<qc::QC_2D> seed                            // will mark which points are used as starting points
        ( qc::GridSize<qc::QC_2D>::createFrom ( grid ) );

      sdf.setAll ( aol::NumberTrait<double>::Inf );  // first, we set all values to plus infinity

      const double radius = 1.0 / 3.0;
      const double bandwidth = sqrt ( 2. ) * grid.H();

      for ( int i = 0; i < grid.getWidth(); ++i ) {
        for ( int j = 0; j < grid.getWidth(); ++j ) {
          aol::Vec2<double> pos ( i, j );
          pos *= grid.H();
          if ( aol::Abs ( pos.norm() - radius ) < bandwidth ) {
            sdf.set ( i, j, pos.norm() - radius );     // next, we initialize points next to the zero interface for which we wish to determine the signed distance function
            seed.set ( i, j, true );            // and declare those points as seed points
          }
          distances.set ( i, j, pos.norm() - radius ); // and we also compute the reference solution ...
        }
      }

      // Now, we set up a DistanceSweeper. It needs to know the seed mask and the underlying grid.
      qc::DistanceSweeper2d<double> distanceSweeper ( seed, grid );
      distanceSweeper.setVerboseMode ( true );

      // For the computation, we need initial values at the seed points and plus infinity values at all other points.
      distanceSweeper.computeDistances ( sdf );

      // Finally, compare the sdf computed by sweeping to the reference solution.
      qc::ScalarArray<double, qc::QC_2D> difference ( grid );
      difference += sdf;
      difference -= distances;

      cerr << "approximate L2 norm of difference = " << difference.norm() / difference.size() << ", Linf norm = " << difference.getMaxAbsValue() << endl;

#if 0
      double maxDiff = difference.getMaxAbsValue();
      for ( int i = 0; i < difference.size(); ++i ) {
        difference[i] = 255 * aol::Abs ( difference[i] ) / maxDiff;
      }
      difference.save ( "difference_circle.pgm" );
#endif
    }

    {
      cerr << "Treating a point in the middle as interface " << endl;

      qc::ScalarArray<double, qc::QC_2D> sdf ( grid );        // to compute the signed distance function
      qc::ScalarArray<double, qc::QC_2D> distances ( grid );  // true distances as reference solution
      qc::BitArray<qc::QC_2D> seed                            // will mark which points are used as starting points
        ( qc::GridSize<qc::QC_2D>::createFrom ( grid ) );

      sdf.setAll ( aol::NumberTrait<double>::Inf );  // first, we set all values to plus infinity

      const double center = 0.5;
      const double bandwidth = sqrt ( 2. ) * grid.H();

      for ( int i = 0; i < grid.getWidth(); ++i ) {
        for ( int j = 0; j < grid.getWidth(); ++j ) {
          aol::Vec2<double> pos ( i, j );
          pos *= grid.H();
          pos[0] -= center;
          pos[1] -= center;
          if ( aol::Abs ( pos.norm() ) < bandwidth ) {
            sdf.set ( i, j, aol::Abs ( pos.norm() ) );    // next, we initialize points next to the zero interface for which we wish to determine the signed distance function
            seed.set ( i, j, true );            // and declare those points as seed points
          }
          distances.set ( i, j, aol::Abs ( pos.norm() ) ); // and we also compute the reference solution ...
        }
      }

      // Now, we set up a DistanceSweeper. It needs to know the seed mask and the underlying grid.
      qc::DistanceSweeper2d<double> distanceSweeper ( seed, grid );
      distanceSweeper.setVerboseMode ( true );

      // For the computation, we need initial values at the seed points and plus infinity values at all other points.
      distanceSweeper.computeDistances ( sdf );

      // Finally, compare the sdf computed by sweeping to the reference solution.
      qc::ScalarArray<double, qc::QC_2D> difference ( grid );
      difference += sdf;
      difference -= distances;

      cerr << "approximate L2 norm of difference = " << difference.norm() / difference.size() << ", Linf norm = " << difference.getMaxAbsValue() << endl;


#if 0
      double maxDiff = difference.getMaxAbsValue();
      for ( int i = 0; i < difference.size(); ++i ) {
        difference[i] = 255 * aol::Abs ( difference[i] ) / maxDiff;
      }
      difference.save ( "difference_point_middle.pgm" );
      double maxSdf = sdf.getMaxAbsValue();
      for ( int i = 0; i < sdf.size(); ++i ) {
        sdf[i] = 255 * aol::Abs ( sdf[i] ) / maxSdf;
      }
      sdf.save ( "sdf_point_middle.pgm" );
#endif

    }
  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
