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

#include <configurators.h>
#include <mcm.h>
#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <integrateLevelSet.h>

#define RealType double

using namespace aol::color;

int main( int argc, char **argv)
{
  try
    {
      if ( argc < 4 )
      {
        cerr << red << "usage: "<<argv[0]<<" input output max_radius [thresholdFactor=1] \n" << reset;
        cerr << "(for 129^2-images a fine max_radius is 0.002)\n";
        return EXIT_FAILURE;
      }

      double tFactor = 1.;      // factor for the threshold

      if (argc == 5) tFactor = atof( argv[4] );

      qc::ScalarArray<RealType, qc::QC_3D> inputimg( argv[1]);
      qc::ScalarArray<RealType, qc::QC_3D> outputimg( inputimg );

      int N = inputimg.getNumX ();
      int d = qc::logBaseTwo (N);
      qc::GridDefinition grid( d, qc::QC_3D );

      RealType maxr = atof(argv[3]);        // delete all particles with radius smaller than maxr


      // generate a characteristic function from the image:
      for (int i=0; i<inputimg.size(); i++)
        if (inputimg[i] > 0) inputimg[i]=1.; else inputimg[i]=0.;


      // computing the mass of the 0-level-set
      qc::ScalarArray<double, qc::QC_3D> massIntegral( N-1,N-1,N-1 );
      qc::ScalarArray<double, qc::QC_3D> vector1( N,N,N );
      vector1.setAll( 1. );
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> massIntegrator( grid, inputimg, vector1 );
      massIntegral.setZero();
      massIntegrator.integrateSaveEachValue( massIntegral );

      // integrate in a ball of radius 3 times bigger than the maximum particle radius
      qc::Ball *scanBall = new qc::Ball( grid, 10.*maxr );
      // delete in a ball that has the size of the radius (this is to avoid deleting the boundary of vessels)
      qc::Ball *removeBall = new qc::Ball( grid, maxr );

      // iterators for the balls
      qc::Ball::iteratorBall ballIt;
      qc::Ball::iteratorBall ballItEnd;
      qc::Ball::iteratorBall ballItBegin;

      // just for displaying the progress
      RealType percent = 0;
      aol::Vec3<double> coords;


      double threshold = 0.;                                        // just count the pixels in the removeBall
      coords[0]=0.5; coords[1]=0.5; coords[2]=0.5;
      ballItEnd = removeBall->end(coords);;
      ballItBegin = removeBall->begin(coords);

      for (ballIt = ballItBegin; ballIt != ballItEnd; ballIt++) threshold++;

      threshold *= tFactor;

      // BEGIN HACK
//       outputimg.setAll( 0. );
//
//       coords[0]=0.2; coords[1]=0.2; coords[2]=0.2;
//       ballItEnd = removeBall->end(coords);;
//       ballItBegin = removeBall->begin(coords);
//
//       for (ballIt = ballItBegin; ballIt != ballItEnd; ballIt++)
//         outputimg.set( (*ballIt).x(), (*ballIt).y(), (*ballIt).z(), 3. );
//
//       coords[0]=0.7; coords[1]=0.7; coords[2]=0.7;
//       ballItEnd = scanBall->end(coords);;
//       ballItBegin = scanBall->begin(coords);
//
//       for (ballIt = ballItBegin; ballIt != ballItEnd; ballIt++)
//         outputimg.set( (*ballIt).x(), (*ballIt).y(), (*ballIt).z(), 3. );
//
//       coords[0]=0.5; coords[1]=0.5; coords[2]=0.5;
//       ballItEnd = removeBall->end(coords);;
//       ballItBegin = removeBall->begin(coords);
//
//       for (ballIt = ballItBegin; ballIt != ballItEnd; ballIt++)
//         outputimg.set( (*ballIt).x(), (*ballIt).y(), (*ballIt).z(), 3. );
//
      // END HACK

      // now start the loop over all elements
      for (int x = 0; x<N; x++)
      {
        // display progress
        percent = static_cast<RealType> ( x * 100. / N );
        cerr.precision(4);
        cerr<<"\rCalculating..."<<percent<<" %            ";

        for (int y = 0; y<N; y++)
        {
          for (int z = 0; z<N; z++)
          {
            coords[0] = x * grid.H();
            coords[1] = y * grid.H();
            coords[2] = z * grid.H();

            // compute the mass in the ball
            ballItEnd   = scanBall->end(coords);;
            ballItBegin = scanBall->begin(coords);
            double mass = 0.;
            for (ballIt = ballItBegin; ballIt != ballItEnd; ballIt++)
              mass += inputimg.get( (*ballIt).x(), (*ballIt).y(), (*ballIt).z() );

            // too less mass => little spot => remove it!
            if (mass < threshold)  {
              ballItEnd   = removeBall->end(coords);;
              ballItBegin = removeBall->begin(coords);
              for (ballIt = ballItBegin; ballIt != ballItEnd; ballIt++)
                outputimg.set( (*ballIt).x(), (*ballIt).y(), (*ballIt).z(), -10. );
            }

          }
        }
      }
      outputimg.save(argv[2], qc::PGM_DOUBLE_BINARY);
    }
  catch ( aol::Exception &ex )
    {
      ex.dump();
    }
}
