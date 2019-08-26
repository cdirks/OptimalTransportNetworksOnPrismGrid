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
#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>

#define REAL double

using namespace aol::color;

/** *************************************************************************
// * differenceNorm.cpp
// * computes the L2- and the L-infinity-norm of the difference of two
// * images
// * Args: image1, image2
// * Author: Oliver Nemitz
// * Date: 24.04.2006
// ************************************************************************* */


int main( int argc, char **argv)
{
  try
    {
      if ( argc < 3 )
      {
        cerr <<"Computes the L1,L2 and L-inifinity-norm of the difference of two images."<<endl;
        cerr <<red<< "usage: "<<argv[0]<<" image1 image2\n"<<reset;
        return EXIT_FAILURE;
      }
      qc::ScalarArray<REAL, qc::QC_3D> image1( argv[1] );
      qc::ScalarArray<REAL, qc::QC_3D> image2( argv[2] );

      cerr << "Computing very complex, complicated and difficult things, please wait...(ca. 1 hour)...";

      double L2Norm = 0.;         //  L^2-norm
      double L1Norm = 0.;         //  L^1-norm
      double LInfNorm = 0.;       //  L^infty-norm
      double Difference = 0.;     //  difference between two pixels/voxels

      for (int i=0; i<image1.length(); i++) {
        Difference = aol::Abs( image1[i] - image2[i] );
        L2Norm += aol::Sqr( Difference );
        L1Norm += Difference;
        if ( Difference > LInfNorm ) LInfNorm = Difference;
      }

      // dividing by the length is equivalent to dividing by getNumX^3.
      L1Norm /= static_cast<double>( image1.length() );
      L2Norm = sqrt( L2Norm / static_cast<double>( image1.length() ) );

      cerr << "finished!\nResults: \n L1-Norm: "<< red << L1Norm << reset <<endl;
      cerr << ", L2-Norm: "<< red << L2Norm << endl<<reset;
      cerr << ", L-Infinity-Norm: "<< red << LInfNorm <<endl<< reset << endl;

    }
  catch ( aol::Exception &ex )
    {
      ex.dump();
    }
}
