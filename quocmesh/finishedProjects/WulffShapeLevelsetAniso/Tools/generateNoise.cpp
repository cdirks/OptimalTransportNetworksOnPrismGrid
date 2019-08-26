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

// *************************************************************************
// * generateNoise.cpp
// * generates some kind of Noise (to be chosen) on the image
// * Args: input, output, noiseNumber
// * Author: Oliver Nemitz
// * Date: 21.04.2006
// *************************************************************************

#define numNoise 15
#define radiusNoise  0.9
#define Epsilon 1

// auxiliary methods

double uniformrandom (double alpha = 0, double beta = 1) { return (static_cast<double> (rand ()) / (RAND_MAX-1)) * (beta - alpha) + alpha; }


// ----------------- generate Noise - methods -----------------------------------

void generateNoiseOnWholeData(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  cerr<<"generating Noise...";
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
    {
        // get the actual value of the point
      double val = Data.get(i,j,k);
        // change it in a range
      val += uniformrandom( -radiusNoise, radiusNoise );
      Data.set(i,j,k, val);
    }
  }
  cerr<<"ready!"<<endl;
}

void generateSaltNPepaNoise( qc::ScalarArray<double, qc::QC_3D> &Data )
{
  cerr<<"generating Noise...";

  aol::NoiseOperator<double> noiseOp;
  noiseOp.apply( Data, Data );

  cerr<<"ready!"<<endl;
}

void generateNoiseOnSupportData(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  cerr<<"generating Noise...";
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
    {
        // get the actual value of the point
      double val = Data.get(i,j,k);
        // change it in a range
      if (val != 0) {
        val += uniformrandom( -radiusNoise, radiusNoise );
        Data.set(i,j,k, val);
      }
    }
  }
  cerr<<"ready!"<<endl;
}

void generateNoiseAroundSupportData(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  cerr<<"generating Noise around supported data...";

  // generate spots on x percent of the points
  double prozent = 50.;  // 50          // only produce noise in about prozent %  of the support
  double prozentNeighbourhood = 30.;  // 40   // and only in about prozentNB % in the nbh of each point
  int radius = 1;     // radius of the neighbourhood

  qc::ScalarArray<double, qc::QC_3D> DataOrig(Data, aol::DEEP_COPY);

  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j) {
      for (int k=0; k<N; ++k) {
        // get the actual value of the point
        double val = DataOrig.get(i,j,k);
        // change it in a range
        if (val > 0.001)
        {
          // only produce noise in about prozent %  of the points
          double number = uniformrandom( 0., 100. );
          if (number < prozent) {
            // now change the values in the neighbourhood
            for (int i2 = i-radius; i2 < i+radius; ++i2)
              for (int j2 = j-radius; j2 < j+radius; ++j2)
                for (int k2 = k-radius; k2 < k+radius; ++k2) {
                    // only produce noise in about prozent %  of the points
                  double numberNBH = uniformrandom( 0., 100. );
                  if (numberNBH < prozentNeighbourhood) {
                    if (i2 >= 0 && i2 < N &&  j2 >= 0 && j2 < N && k2 >= 0 && k2 < N ) {
                      double val2 = DataOrig.get(i2, j2, k2);
                      if (val2 == 0) val2 = uniformrandom( val-radiusNoise, val+radiusNoise );
                      else val2 += uniformrandom( -radiusNoise, radiusNoise );
                      Data.set(i2, j2 ,k2 , val2);
                    }
                  }
                }
          }
        }
      }
    }
  }
  cerr<<"ready!"<<endl;
}

// background is NOT 0, but -1
void generateNoiseAroundSupportDataBGNeg(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  double radNoise = 0.3;
  cerr<<"generating Noise...";
  // generate spots on x percent of the points
  double prozent = 30.;  // 50
  double prozentNeighbourhood = 20.;  // 30

  qc::ScalarArray<double, qc::QC_3D> DataOrig(Data, aol::DEEP_COPY);

  int radius = 1;
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
    {
        // get the actual value of the point
      double val = DataOrig.get(i,j,k);
        // change it in a range
      if (val != -1.)
      {
          // only produce noise in about prozent %  of the points
        double number = uniformrandom( 0., 100. );
        if (number < prozent) {
            // now change the values in the neighbourhood
          for (int i2 = i-radius; i2 < i+radius; ++i2)
          {
            for (int j2 = j-radius; j2 <= j+radius; ++j2)
            {
              for (int k2 = k-radius; k2 <= k+radius; ++k2)
              {
                  // only produce noise in about prozent %  of the points
                double numberNBH = uniformrandom( 0., 100. );
                if (numberNBH < prozentNeighbourhood) {
                  if (i2 >= 0 && i2 < N &&  j2 >= 0 && j2 < N && k2 >= 0 && k2 < N ) {
                    double val2 = DataOrig.get(i2, j2, k2);
                    if (val2 == -1.) val2 = uniformrandom( val-radNoise, val+radNoise );
                    else val2 += uniformrandom( -radNoise, radNoise );
                    Data.set(i2, j2 ,k2 , val2);
                  }
                }
              }  // k2
            }  // j2
          }  // i2
        }  // if number < prozent
      }  // if val...
    }  // k
  }
  cerr<<"ready!"<<endl;
}


void generateSpotsOnSupportData(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  // generate spots on x percent of the points
  double prozent = 5.;

  cerr<<"generating Noise...";
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
    {
        // get the actual value of the point
      double val = Data.get(i,j,k);
        // change it in a range
      if (val != -1.) {
          // only produce spots around prozent %  of the points
        double number = uniformrandom( 0., 100. );
        if (number < prozent) {
            // get a random color for the spot
          double spotColor = uniformrandom( 0.01, 0.9 );
            // generate a spot in a neighbourhood of a random radius
          int radius = static_cast<int>( uniformrandom( 1., 3. ) );
          for (int l = i-radius; l < i+radius; ++l)
            for (int m = j-radius; m < j+radius; ++m)
              for (int n = k-radius; n < k+radius; ++n)
                if (l >= 0 && l < N &&  m >= 0 && m < N && n >= 0 && n < N )
                  if (Data.get(l,m,n) != 0)
                    Data.set(l,m,n, spotColor);

        }
//           val += uniformrandom( -radiusNoise, radiusNoise );
//           Data.set(i,j,k, val);
      }
    }
  }
  cerr<<"ready!"<<endl;
}

// -------------------------------------------------------------------------------
// ------------------------ the main program -------------------------------------
// -------------------------------------------------------------------------------

int main( int argc, char **argv)
{
  try
    {
      if ( argc < 4 )
      {
        cerr <<aol::color::red<< "usage: "<<argv[0]<<" input output\n"<<aol::color::black;
        return EXIT_FAILURE;
      }
      qc::ScalarArray<REAL, qc::QC_3D> img( argv[1] );
      int N = img.getNumX();
//       int noiseNumber = atoi( argv[3] );

      generateNoiseAroundSupportData(img, N);

      cerr<<"done! Saving...\n"<<aol::color::black;

      img.save(argv[2]);

      cerr<<aol::color::blue<<"done! Thanx for using generateNoise, THE standard utility for generating your personal noise!\n";
      cerr<<aol::color::black;
    }
  catch ( aol::Exception &ex )
    {
      ex.dump();
    }
}
