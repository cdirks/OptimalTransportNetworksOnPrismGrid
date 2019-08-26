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
#include <qmException.h>
#include <aol.h>
#include <quoc.h>
#include <op.h>
#include <parameterParser.h>
#include "grapeInterface3d.h"
#include <narrow.h>
#include <signedDistanceOp.h>
#include <scalarArray.h>
#include <gridBase.h>
#include <linearSmoothOp.h>

typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;

using namespace std;
using namespace aol::color;


/* **********************************************************************
 * file: generateCurve1.cpp
 * generates a very easy levelset-curve, namely the line from
 * (0,0) to (1,1).
 * ********************************************************************** */
#define numNoise 15
#define radiusNoise  0.03
#define Epsilon 1

double uniformrandom (double alpha = 0, double beta = 1) { return (static_cast<double> (rand ()) / (RAND_MAX-1)) * (beta - alpha) + alpha; }

double abs2(double a)
{   if (a>=0) return a; else return(-a); }

double max(double a, double b)
{
  if (a>b) return a;
  else return(b);
}

double sqr(double a) { return a*a; }

int isqr(int a) { return a*a; }

int cutBorderCoord(double coord, int N)
{
  int icoord=static_cast<int>( coord );
  if (icoord<0)  icoord = 0;
  if (icoord>=N) icoord = N-1;
  return icoord;
}


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
  cerr<<"generating Noise...";
  // generate spots on x percent of the points
  double prozent = 50.;  // 50
  double prozentNeighbourhood = 40.;  // 40

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
        if (val != 0)
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



// ------------ Generate SignedDistance-Function to given data ---------------
void generateSignedDistanceFunction(qc::ScalarArray<double, qc::QC_3D> &Data, int /*N*/, int d)
{
  // first copy the data into a temporary array
  qc::ScalarArray<double, qc::QC_3D> Temp(Data);
//   Temp.save( "test.bz2" );

  qc::GridDefinition grid( d, qc::QC_3D );
//   qc::SignedDistanceOp<double> BuildSD( grid );

  qc::SignedDistanceOp3D<double>( grid ).apply( Temp, Data );


//   BuildSD.setIsoValue( 0. );

//   BuildSD.apply( Temp, Data );

}


// -------------------- Generate data - methods ---------------------------

void drawLocalBall(qc::ScalarArray<double, qc::QC_3D> &Data, int N, int x, int y, int z, double radius)
{
  int rad = static_cast<int>(radius);

  for (int i=x-rad; i<=x+rad; i++)
    for (int j=y-rad; j<=y+rad; j++)
      for (int k=z-rad; k<=z+rad; k++)
      {
        if (i >= 0 && i < N &&  j >= 0 && j < N && k >= 0 && k < N ) {
          if ( isqr(x-i) + isqr(y-j) + isqr(z-k) <= sqr(radius) )
            Data.set(i, j, k, 1.);
        }
      }
}


void drawLocalEllipsoid(qc::ScalarArray<double, qc::QC_3D> &Data, int N, int x, int y, int z, double a, double b, double c)
{
  int A = static_cast<int>(a);
  int B = static_cast<int>(b);
  int C = static_cast<int>(c);

  for (int i=x-A; i<=x+A; i++)
    for (int j=y-B; j<=y+B; j++)
      for (int k=z-C; k<=z+C; k++)
  {
    if (i >= 0 && i < N &&  j >= 0 && j < N && k >= 0 && k < N ) {
      if ( sqr( (x-i)/a ) + sqr( (y-j)/b ) + sqr( (z-k)/c ) <= 1. )
        Data.set(i, j, k, 1.);
    }
  }
}

void generate2Balls(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setAll(-1.);

  // zwei übereinanderliegende Kugeln
  drawLocalBall(Data, N, N/2, N/4, N/2, 4.);
  drawLocalBall(Data, N, N/2, 3*N/4 + 1, N/2, 4.);
}

void generateBall(qc::ScalarArray<double, qc::QC_3D> &Data, int N, int x, int y, int z, double radius)
{
  Data.setAll(-1.);
  drawLocalBall(Data, N, x, y, z, radius);
}

void generateEllipsoidBlock(qc::ScalarArray<double, qc::QC_3D> &Data, int N, double a, double b, double c)
{
  // Mittelpunkt des Ellipsoids
  int x = N/2;
  int y = N/2;
  int z = N/2;

  // Abmessung für Schleife
  int radx = static_cast<int>(a);
  int rady = static_cast<int>(b);
  int radz = static_cast<int>(c);

  for (int i=x-radx; i<=x+radx; i++)
    for (int j=y-rady; j<=y+rady; j++)
      for (int k=z-radz; k<=z+radz; k++)
      {
        if (i >= 0 && i < N &&  j >= 0 && j < N && k >= 0 && k < N ) {
          double xd = static_cast<double>(x-i);
          double yd = static_cast<double>(y-j);
          double zd = static_cast<double>(z-k);

          if ( sqr(xd / a) + sqr(yd / b) + sqr(zd / c) <= 1. )
            Data.set(i, j, k, 1.);
        }
      }
}

// generates an ellipsoid with half-axis a,b,c
void generateEllipsoid(qc::ScalarArray<double, qc::QC_3D> &Data, int N, double a, double b, double c)
{
  // Mittelpunkt des Ellipsoids
  int x = N/2;
  int y = N/2;
  int z = N/2;

  double h = 1 / static_cast<double>( N-1 );

  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
      {
        double xd = static_cast<double>(x-i) * h;
        double yd = static_cast<double>(y-j) * h;
        double zd = static_cast<double>(z-k) * h;

        Data.set(i, j, k, sqrt( sqr(xd / a) + sqr(yd / b) + sqr(zd / c) ) - 0.15 );
      }
}

void generateCapsule( qc::ScalarArray<double, qc::QC_3D> &Data, int N, double radius )
{
  double h = 1. / static_cast<double>( N-1. );

  // define z-offset of the cylinder to the x,y-planes
  int zOffset = static_cast<int>( 3. * radius * N );

  // traverse all nodes
  for ( int i=0; i<N; i++ ) {
    for ( int j=0; j<N; j++ ) {
      for ( int k=0; k<N; k++ ) {
        double x = (i * h) - 0.5;
        double y = (j * h) - 0.5;
        double z = (k * h) - 0.5;

        double value = -10.;

        // draw a cylinder parallel to the z-axis in the z-interval [zOffset,N-zOffset]
        if ( k>=zOffset && k<N-zOffset ) value = sqrt( x*x + y*y );
        // draw a halfsphere in the z-interval [0,zOffset]
        if ( k<zOffset ) value = sqrt( x*x + y*y + 0.5*aol::Sqr( z + 0.5 - 3.*radius ) );
        // draw a halfsphere in the z-interval [N-zOffset,N]
        if ( k>=N-zOffset ) value = sqrt( x*x + y*y + 0.5*aol::Sqr( z - 0.5 + 3.*radius ) );

        Data.set( i,j,k, value-radius );               // 0-isosurface is sphere with radius sphereRadius
      }
    }
  }
}


void generateOval( qc::ScalarArray<double, qc::QC_3D> &Data, int N, double radius )
{
  double h = 1. / static_cast<double>( N-1. );

  // define z-offset of the cylinder to the x,y-planes
  int zOffset = static_cast<int> ( 3. * radius * N );

  // traverse all nodes
  for ( int i=0; i<N; i++ ) {
    for ( int j=0; j<N; j++ ) {
      for ( int k=0; k<N; k++ ) {
        double x = (i * h) - 0.5;
//         double y = (j * h) - 0.5;
        double z = (k * h) - 0.5;

        double value = -10.;

        // draw a cylinder parallel to the z-axis in the z-interval [zOffset,N-zOffset]
        if ( k>=zOffset && k<N-zOffset ) value = fabs( x );
        // draw a halfsphere in the z-interval [0,zOffset]
        if ( k<zOffset ) value = sqrt( x*x + 0.5*aol::Sqr( z + 0.5 - 3.*radius ) );
        // draw a halfsphere in the z-interval [N-zOffset,N]
        if ( k>=N-zOffset ) value = sqrt( x*x + 0.5*aol::Sqr( z - 0.5 + 3.*radius ) );

        Data.set( i,j,k, value-radius );               // 0-isosurface is sphere with radius sphereRadius
      }
    }
  }
}

// generates a stripe in y-direction with ellipses of halfaxis a and b
void generateEllipseStripe( qc::ScalarArray<double, qc::QC_3D> &Data, int N, double a, double b, double radius  )
{
  double h = 1. / static_cast<double>( N-1. );

  Data.setAll( 10. );

  // traverse all nodes
  for ( int i=0; i<N; i++ ) {
    for ( int j=5; j<N-5; j++ ) {
      for ( int k=0; k<N; k++ ) {
        double x = (i * h) - 0.5;
        double z = (k * h) - 0.5;

        // compute ellipse-value
        double value = sqrt( aol::Sqr( x/a ) + aol::Sqr( z/b ) );

        Data.set( i,j,k, value-radius );        // 0-isosurface is an ellipse with positive radius
      }
    }
  }
}


void generateNoisyBall(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setAll(-1.);

  // eine verrauschte Kugel
  drawLocalBall(Data, N, N/2, N/2, N/2, 7.);
  generateNoiseAroundSupportDataBGNeg( Data, N );
}


void generateVesselPiece(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  double radius = 3;
  Data.setAll(-1.);

  // draw the first cigar-like shape
  for (int i=N/12; i<N/3 + 2; ++i)
  {
    cerr<<".";
    drawLocalBall(Data, N, i,N/2, N/2, radius );
  }
  // draw the second cigar-like shape
  for (int i=2*N/3 - 2; i<11*N/12; ++i)
  {
    cerr<<".";
    drawLocalBall(Data, N, i,N/2, N/2, radius );
  }

}


void generateXAxisLine(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  double radius = 2;
  Data.setAll(-1.);

  // draw the first cigar-like shape
  for (int i=8; i<N-8; ++i)
  {
    cerr<<".";
    drawLocalBall(Data, N, i,N/2, N/2, radius );
  }
}

void generateXAxisEllipsoidalLine(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  double a = 2.;
  double b = 2.;
  double c = 4.;

  // draw the first cigar-like shape
  for (int i=16; i<N-16; ++i)
  {
    cerr<<".";
    drawLocalEllipsoid(Data, N, i, N/2, N/2, a,b,c );
  }
}



void generateDistancePlane(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setAll(-1.);

  // draw the first cigar-like shape
  for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
        Data.set( i,j,k, static_cast<double>(k)/static_cast<double>(N) -0.5 );

}



void generateDiagLine(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        if ( (i==k) && (j==k) )
        {
          // line with thickness 2*Epsilon
          for (int l=i-Epsilon; l<=i+Epsilon; l++)
            for (int m=j-Epsilon; m<=j+Epsilon; m++)
              for (int n=k-Epsilon; n<=k+Epsilon; n++)
                if (l >= 0 && l <N &&  m >= 0 && m <N && n >= 0 && n <N )
                  Data.set(l,m,n, 1.);
        }
        else
          Data.set(i,j,k, -1.);
      }
  }
}



void generateStraitJunction(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  int LocEpsilon = 2;
  Data.setAll(-1.);

  int mid = N/2;
  // generate two lines: one non-stop line and one forking line
  // first the non-stop line (horizontal):
  for (int i=0; i<N; ++i) {
    cerr<<".";
    // line with thickness 2*Epsilon
    for (int m = mid-LocEpsilon; m <= mid+LocEpsilon; m++)
      for (int n = mid-LocEpsilon; n <= mid+LocEpsilon; n++)
        if (m >= 0 && m < N && n >= 0 && n < N )
          if ( (n-mid)*(n-mid) + (m-mid)*(m-mid) <= LocEpsilon*LocEpsilon )
            Data.set(i,m,n, 1.);
  }

  LocEpsilon = 1;
  // now the junction, that is a diagonal line from the center to the border
  for (int i=mid; i<N; ++i)
  {
    cerr<<".";
    for (int j=mid; j<N; ++j)
      for (int k=mid+1; k<N; ++k)
      {
        if ( (i==k) && (j==k) )
        {
          // line with thickness 2*Epsilon
          for (int l = i-LocEpsilon; l <= i+LocEpsilon; l++)
            for (int m = j-LocEpsilon; m <= j+LocEpsilon; m++)
              for (int n = k-LocEpsilon; n <= k+LocEpsilon; n++)
                if (l >= 0 && l <N &&  m >= 0 && m <N && n >= 0 && n <N )
                  Data.set(l-3,m,n, 1.);
        }
      }
  }
}

void generateNoisyStraitJunction(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  generateStraitJunction( Data, N );
  generateNoiseAroundSupportDataBGNeg( Data, N );
}

void generateNoisyDiagLine(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  generateDiagLine(Data, N);
  generateNoiseAroundSupportDataBGNeg(Data, N);
}


void generateL1Norm(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        double x = (static_cast<double>(i)) / (static_cast<double>(N-1));
        double y = (static_cast<double>(j)) / (static_cast<double>(N-1));
        double z = (static_cast<double>(k)) / (static_cast<double>(N-1));

        Data.set(i,j,k, abs2(x-0.5)+abs2(y-0.5)+abs2(z-0.5));
      }
  }
}


void generateSphere(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
    {
      double x = (static_cast<double>(i)) / (static_cast<double>(N-1));
      double y = (static_cast<double>(j)) / (static_cast<double>(N-1));
      double z = (static_cast<double>(k)) / (static_cast<double>(N-1));

      Data.set( i,j,k, sqrt( x*x + y*y + z*z ) );
    }
  }
}

void generateSupNorm(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
    {
      double x = (static_cast<double>(i)) / (static_cast<double>(N-1));
      double y = (static_cast<double>(j)) / (static_cast<double>(N-1));
      double z = (static_cast<double>(k)) / (static_cast<double>(N-1));

      Data.set(i,j,k, max( max( abs(x-0.5), abs(y-0.5) ), abs2(z-0.5) ) );
    }
  }
}


void generateBlock(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  double offset = N/4;        // offset from the block to the border
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
          if (i>=offset && i<=N-offset && j>=offset && j<=N-offset &&  k>=offset && k<=N-offset)
            Data.set(i,j,k, -1. );
          else
            Data.set(i,j,k, 1. );
      }
  }
}

void generateBlockBorder(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  double offset = N/4;        // offset from the block to the border
  int rad = 0;
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        if ( (i>=offset-rad && i<=N-offset+rad && j>=offset-rad && j<=N-offset+rad &&    // Deckel + Boden
               ( (k>=offset-rad && k<=offset+rad) || (k>=N-offset-rad && k<=N-offset+rad) ) ) ||
             (i>=offset-rad && i<=N-offset+rad && k>=offset-rad && k<=N-offset+rad &&
             ( (j>=offset-rad && j<=offset+rad) || (j>=N-offset-rad && j<=N-offset+rad) ) ) ||
             (j>=offset-rad && j<=N-offset+rad && k>=offset-rad && k<=N-offset+rad &&
             ( (i>=offset-rad && i<=offset+rad) || (i>=N-offset-rad && i<=N-offset+rad) ) ) )
            Data.set(i,j,k, -1. );
        else Data.set(i,j,k, 1. );
      }
  }
}



void generateCube(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setAll( 0. );
  int L = N/2;
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
          double x = (static_cast<double>(i-L)) / (static_cast<double>(N-1));
          double y = (static_cast<double>(j-L)) / (static_cast<double>(N-1));
          double z = (static_cast<double>(k-L)) / (static_cast<double>(N-1));

          Data.set(i,j,k, max( max(abs(x),abs(y)), max(abs(x),abs(z)) ) - 0.2 );
      }
  }
}

void generateNoisyCube(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  generateCube(Data, N);
  generateNoiseOnSupportData(Data, N);
}


void generateEllipsoidalLevel(qc::ScalarArray<double, qc::QC_3D> &Data, int N, double a, double b, double c)
{
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        double x = (static_cast<double>(i)) / (static_cast<double>(N-1));
        double y = (static_cast<double>(j)) / (static_cast<double>(N-1));
        double z = (static_cast<double>(k)) / (static_cast<double>(N-1));

        Data.set(i,j,k, sqrt(sqr( (x-0.5)/a )+sqr( (y-0.5)/b )+sqr( (z-0.5)/c )) );
      }
  }
}

void generateTubes(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setZero();
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        double y = (static_cast<double>(j)) / (static_cast<double>(N-1));
        double z = (static_cast<double>(k)) / (static_cast<double>(N-1));

        double rad = sqrt( sqr(y-0.5) + sqr(z-0.5) );

        Data.set(i,j,k, 1./(1.+rad) - 0.2);
      }
  }
}

void generateNoisyHorizontalTubes(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  // first the tubes
  cerr<<"generating tubes...";
  Data.setZero();
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        double y = (static_cast<double>(j)) / (static_cast<double>(N-1));
        double z = (static_cast<double>(k)) / (static_cast<double>(N-1));

        double rad = sqrt( sqr(y-0.5) + sqr(z-0.5) );

        Data.set(i,j,k, 1./(1.+rad) - 0.2);
      }
  }

  // second: put some noise around it
//   generateNoiseOnWholeData(Data, N);
}

void generateAnyNoisyTubes(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  // first define the direction in which the tube lies
  aol::Vec3<double> _v(1.,1.,1.);
  _v /= _v.norm();

  // now calc the matrix that represents the basis transformation
  aol::Matrix33<double> _B;
  _B.setRow(0, _v );
    // the other directions must be orthogonal to _v, but then, they can be
    // chosen arbitrary, first choose y orthogonal to _v
    // TODO: arranging by size to avoid taking the smallest values (=> stability)
    aol::Vec3<double> y(0, 0, 0);
    if (_v[1] != 0) {
      if (_v[2] != 0) {
        y[1] = -_v[2];
        y[2] = _v[1];
        y /= y.norm();
      }
      else y[2] = 1;
    }
    else y[1] = 1;

    _B.setRow(1, y );

    // now z = vector-product of v and y, both normed => z is normed too
    aol::Vec3<double> z( _v.crossProduct(y) );
    _B.setRow(2, z );

  cerr<<"x: "<<_v<<", y: "<<y<<", z: "<<z<<endl;
  cerr<<"B: "<<_B<<endl;


  // first the tubes
  cerr<<"generating tubes...";
  Data.setZero();
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        double x1 = (static_cast<double>(i)) / (static_cast<double>(N-1));
        double y1 = (static_cast<double>(j)) / (static_cast<double>(N-1));
        double z1 = (static_cast<double>(k)) / (static_cast<double>(N-1));

        aol::Vec3<double> temp( x1-0.5,y1-0.5,z1-0.5 );
        _B.mult(temp, z);

        double rad = sqrt( sqr(z[1]) + sqr(z[2]) );

        Data.set(i,j,k, 1./(1.+rad) - 0.9 );
//         Data.add(i,j,k, -0.3);
      }
  }

  // second: put some noise around it
//   generateNoiseAroundSupportData(Data, N);
   generateNoiseOnWholeData(Data, N);
//   generateSpotsOnSupportData(Data, N);
}


void deleteBlock(qc::ScalarArray<double, qc::QC_3D> &Data, int x0, int y0, int z0,
                                                  int x1, int y1, int z1, double val)
{
  for (int i=x0; i<x1; i++)
    for (int j=y0; j<y1; j++)
      for (int k=z0; k<z1; k++)
        Data.set(i,j,k, val);
}

void generateAnyLine(qc::ScalarArray<double, qc::QC_3D> &Data, int N, double x0, double y0, double z0)
{
  // first define the direction in which the tube lies
  aol::Vec3<double> v(1.,0.,0.);

  // draw the line using a parameter t

  double xr,yr,zr;
  int xi,yi,zi;

  for (double t=-N; t<N; t+=0.2)
  {
    xr = x0 + t*v[0];
    yr = y0 + t*v[1];
    zr = z0 + t*v[2];

    xi = static_cast<int> (xr);
    yi = static_cast<int> (yr);
    zi = static_cast<int> (zr);

    if (xi>0 && xi<N && yi>0 && yi<N && zi>0 && zi<N )
      Data.set(xi,yi,zi, 1.);
  }
}

void generateAnyThickLine(qc::ScalarArray<double, qc::QC_3D> &Data, int N, int My )
{
//   Data.setZero();

//   int My = N / 2;
  int Mz = N / 2;

  double radius = 1.;

  for (double j=My-radius; j<My+radius; j+=0.2)
    for (double k=Mz-radius; k<Mz+radius; k+=0.2)
    {
      if ( sqr(My-j) + sqr(Mz-k) < radius )
        generateAnyLine( Data, N, 0., j, k);
    }

}



void generateHorizontalLine(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setZero();
  // first: the thick line
  cerr<<"generating line...";
  double y,z,rad,val;
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        y = (static_cast<double>(j)) / (static_cast<double>(N-1));
        z = (static_cast<double>(k)) / (static_cast<double>(N-1));

        rad = sqrt( sqr(y-0.5) + sqr(z-0.5) );
        val = 1.; //uniformrandom (0, 1);
        if (rad  <= 0.06) Data.set(i,j,k, val);
      }
  }
  /*// second: put some noise around it
  cerr<<"\ngenerating noise...";
  for (int i=0; i<N; ++i)
  {
    cerr<<".";
    i1 = i + static_cast<int>(uniformrandom (-radiusNoise, radiusNoise));

    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        double y = ((double) j) / ((double) N-1);
        double z = ((double) k) / ((double) N-1);

        double rad = sqrt( sqr(y-0.5) + sqr(z-0.5) );
        if (rad  <= 0.06) Data.set(i,j,k, 1.);
      }
  }*/
  cerr<<endl;
}


void generateSimpleArk(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setAll(-1.);
  // first: the thick line
  cerr<<"generating simple ark...";

  // the ark is realized as a piece of the torus (see below)
  // epicenter of the basis slice
  int MxKlein = N/4;                             // epicenter of the little circle
  int MzKlein = N/2;
  int Mx = static_cast<int> (1.3*N);             // epicenter of the whole torus
  int xNeu, yNeu;                                // rotated coordinates

  // radius of the torus
  double r = static_cast<double> (N / 15.);      // radius of the little circle (war 30)
  double dist = 0;
  int ir = static_cast<int> (r + 1.);
  double angle = 0;                              // angle for rotating the slice
  double angleStep = 1. / (30 * N);              // should be enough :-)

  // now generate the basis slice, which will be rotated
  for (int i=MxKlein-ir; i<MxKlein+ir; i++)
    for (int j=MzKlein-ir; j<MzKlein+ir; j++)
    {
      dist = sqrt( sqr(MxKlein-i) + sqr(MzKlein-j) );
      if ( dist < r)
      {
        // rotate the point along a piece of the torus (about 60 deg)
        for (angle = 0; angle < 1.8; angle += angleStep)
        {
          xNeu = Mx - static_cast<int> (cos(angle) * (Mx - i));
          yNeu = static_cast<int> (sin(angle) * (Mx - i));
          if (xNeu >= 0 && xNeu < N && yNeu >= 0 && yNeu < N)
            Data.set(xNeu,yNeu,j, 1.);
        }
      }
    }

}


void generateNoisySimpleArk(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  generateSimpleArk(Data, N);
  generateNoiseAroundSupportDataBGNeg(Data, N);
}


void generateArk(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setZero();
  // first: the thick line
  cerr<<"generating ark...";

  // the ark is realized as a piece of the torus (see below)
  // epicenter of the basis slice
  int MxKlein = N/4;       // epicenter of the little circle
  int MzKlein = N/2;
  int Mx = static_cast<int> (1.3*N);        // epicenter of the whole torus
  int xNeu, yNeu;          // rotated coordinates

  // radius of the torus
  double r = static_cast<double> (N / 15.);    // radius of the little circle (war 30)
  double dist = 0;
  int ir = static_cast<int> (r + 1.);
  double angle = 0;        // angle for rotating the slice
  double angleStep = 1. / (30 * N);    // should be enough :-)

  // now generate the basis slice, which will be rotated
  for (int i=MxKlein-ir; i<MxKlein+ir; i++)
    for (int j=MzKlein-ir; j<MzKlein+ir; j++)
    {
      dist = sqrt( sqr(MxKlein-i) + sqr(MzKlein-j) );
      if ( dist < r)
      {
        // rotate the point along a piece of the torus (about 60 deg)
        for (angle = 0; angle < 1.8; angle += angleStep)
        {
          xNeu = Mx - static_cast<int> (cos(angle) * (Mx - i));
          yNeu = static_cast<int> (sin(angle) * (Mx - i));
          if (xNeu >= 0 && xNeu < N && yNeu >= 0 && yNeu < N)
            Data.set(xNeu,yNeu,j, 1. - dist/r);
        }
      }
    }

}

void generateNoisyArk(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  generateArk(Data, N);
//   generateSpotsOnSupportData( Data, N );
  generateNoiseAroundSupportData(Data, N);
}


void generateTorus(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setAll(0.);
  // first: the thick line
  cerr<<"generating Torus...";

  // epicenter of the basis slice
  int MxKlein = N/4;       // epicenter of the little circle
  int MyKlein = N/2;
  int Mx, Mz;              // epicenter of the whole torus
  Mx = Mz = N/2;
  int xNeu, zNeu;          // rotated coordinates

  // radius of the torus
  double r = static_cast<double> (N / 50.);    // radius of the little circle
  double dist = 0;
  int ir = static_cast<int> (r + 1.);
  double angle = 0;        // angle for rotating the slice
  double angleStep = 1. / (30 * N);    // should be enough :-)

  // now generate the basis slice
  for (int i=MxKlein-ir; i<MxKlein+ir; i++)
    for (int j=MyKlein-ir; j<MyKlein+ir; j++)
    {
      dist = sqrt( sqr(MxKlein-i) + sqr(MyKlein-j) );
      if ( dist < r)
      {
        // rotate the point along the whole torus
        for (angle = 0; angle < 2.*3.141592653; angle += angleStep)
        {
          xNeu = Mx + static_cast<int> (cos(angle) * (i - Mx));
          zNeu = Mz + static_cast<int> (sin(angle) * (i - Mx));
          Data.set(xNeu,j,zNeu, 1. - dist /r );
        }
      }
    }

}

void generateNoisyTorus(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  generateTorus(Data, N);
  generateNoiseOnSupportData(Data, N);
}




void generateSimpleTorus(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setAll(0.);
  // first: the thick line
  cerr<<"generating Torus...";

  // epicenter of the basis slice
  int MxKlein = N/4;       // epicenter of the little circle
  int MyKlein = N/2;
  int Mx, Mz;              // epicenter of the whole torus
  Mx = Mz = N/2;
  int xNeu, zNeu;          // rotated coordinates

  // radius of the torus
  double r = static_cast<double> (N / 50.);    // radius of the little circle
  double dist = 0;
  int ir = static_cast<int> (r + 1.);
  double angle = 0;        // angle for rotating the slice
  double angleStep = 1. / (30 * N);    // should be enough :-)

  // now generate the basis slice
  for (int i=MxKlein-ir; i<MxKlein+ir; i++)
    for (int j=MyKlein-ir; j<MyKlein+ir; j++)
    {
      dist = sqrt( sqr(MxKlein-i) + sqr(MyKlein-j) );
      if ( dist < r)
      {
        // rotate the point along the whole torus
        for (angle = 0; angle < 2.*3.141592653; angle += angleStep)
        {
          xNeu = Mx + static_cast<int> (cos(angle) * (i - Mx));
          zNeu = Mz + static_cast<int> (sin(angle) * (i - Mx));
          Data.set(xNeu,j,zNeu, 1.);
        }
      }
    }

}

void generateNoisySimpleTorus(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  generateSimpleTorus(Data, N);
  generateNoiseAroundSupportDataBGNeg(Data, N);
}



void generateCurvedJunction(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  Data.setAll(0.);
  // first: the thick line
  cerr<<"generating curved Junction...";

  // The curved junction consists of two arks, first the big ark:

  // the ark is realized as a piece of the torus (see below)
  // epicenter of the basis slice
  int MxKlein = N/3;                           // epicenter of the little circle
  int MzKlein = N/2;
  int Mx = static_cast<int> (1.2*N);           // epicenter of the whole torus
  int My;
  int xNeu, yNeu;                              // rotated coordinates

  // radius of the torus
  double r = static_cast<double> (N / 25.);    // radius of the little circle (war 30)
  double dist = 0;
  int ir = static_cast<int> (r + 1.);
  double angle = 0;                            // angle for rotating the slice
  double angleStep = 1. / (30 * N);            // should be enough :-)

  // now generate the basis slice, which will be rotated
  for (int i=MxKlein-ir; i<MxKlein+ir; i++)
    for (int j=MzKlein-ir; j<MzKlein+ir; j++)
    {
      dist = sqrt( sqr(MxKlein-i) + sqr(MzKlein-j) );
      if ( dist < r)
      {
        // rotate the point along a piece of the torus (about 60 deg)
        for (angle = 0.; angle < 2.1; angle += angleStep)
        {
          xNeu = Mx - static_cast<int> (cos(angle) * (Mx - i));
          yNeu = static_cast<int> (sin(angle) * (Mx - i));
          if (xNeu >= 0 && xNeu < N && yNeu >= 0 && yNeu < N)
            Data.set(xNeu,yNeu,j, 1. );
        }
      }
    }


  // now the smaller ark
  MxKlein = N/2;                               // epicenter of the little circle
  MzKlein = N/2;
  // "traditional" epicenter: Mx = 0, Mz = N/2
  // Bifurcation1: Mx = 10, My = 2            // these absolute values hold for N=65
  // Bifurcation2: Mx = 0, My = 32
  // Bifurcation2: Mx = 0, My = 41
  Mx = 0  ; //static_cast<int> (-0.05*N);             // epicenter of the whole torus
  My = 41;

  // radius of the torus
  r = static_cast<double> (N / 30.);           // radius of the little circle (war 30)
  ir = static_cast<int> (r + 1.);

  // now generate the basis slice, which will be rotated
  for (int i=MxKlein-ir; i<MxKlein+ir; i++)
    for (int j=MzKlein-ir; j<MzKlein+ir; j++)
    {
      dist = sqrt( sqr(MxKlein-i) + sqr(MzKlein-j) );
      if ( dist < r)
      {
        // rotate the point along a piece of the torus (about 60 deg)
        // for the "traditional" curve it was angle = [0,2.1]
        // Bifurcation1: angle = [1.,3.1]
        // Bifurcation2: angle = [0.,2.1]
        // Bifurcation3: angle = [-0.5,1.1]
        for (angle = -0.5; angle < 1.1; angle += angleStep)
        {
          xNeu = Mx + static_cast<int> (cos(angle) * (Mx + i));
          yNeu = My + static_cast<int> (sin(angle) * (Mx + i));
          if (xNeu >= 0 && xNeu < N && yNeu >= 0 && yNeu < N)
            Data.set(xNeu,yNeu,j, 1. );
        }
      }
    }

    // this kills a little overlap in the 3. bifurcation
//     for (int z=0; z<N; z++) {
//       Data.set( 35,41,z, 0. );
//       Data.set( 35,42,z, 0. );
//       Data.set( 34,41,z, 0. );
//     }
}



void drawCroppedCone(int i, int j, qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  double h  = static_cast<double>( N/5. );    // height of the cropped cone
  double RB = static_cast<double>( N/30.);    // bigger radius of the cone
  double rs = static_cast<double>( N/100.);   // smaller radius

  double grad = (RB - rs) / h;                // gradient of the walls of the cone

  double offset = N/2.+2.-h;

//   if (i==15 && j==16) cerr<<"\nh: "<<h<<", RB: "<<RB<<", rs: "<<rs<<", grad: "<<grad;
//   if (i==15 && j==16) cerr<<", offset: "<<offset<<", N/2-2: "<<N/2.-2.<<endl;

  // now the loop over the cuboid
  for ( double z=offset; z<N/2.+2.; z+=0.5 )
  {
    double radMax = grad * (z-offset) + rs;

    for ( double x=i-radMax; x<i+radMax; x+=0.5 )
      for ( double y=j-radMax; y<j+radMax; y+=0.5 )
        if (sqr(x-i) + sqr(y-j) < radMax*radMax)
        {
          int xi = cutBorderCoord( x, N );
          int yi = cutBorderCoord( y, N );
          int zi = cutBorderCoord( z, N );
          Data.set( xi,yi,zi, 0. );
        }
  }
}

void generateWrinkleCurve(qc::ScalarArray<double, qc::QC_3D> &Data, int N, int d)
{
  cerr<<aol::color::red<<"\ninitializing data...";
  Data.setAll(0.);
  // first: fill the lower half of the cube with 1.
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0;  k<N/2; k++)
        Data.set( i,j,k, 1. );

  // now load the 2d-picture that contains the way of the wrinkle
  cerr<<"loading way of the wrinkle...\n"<<aol::color::green;
  qc::GridDefinition grid( d, qc::QC_2D );
  qc::ScalarArray<double, qc::QC_2D> theWay( grid );
  theWay.load( "Data/Wrinkle/wrinkleWayIntegral.pgm" );

  // now dig the wrinkle where the picture says to do so
  cerr<<aol::color::red<<"\ndigging the wrinkle...";
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      if ( theWay.get( i,j ) != 0 )
        drawCroppedCone( i,j, Data, N );

  cerr<<aol::color::red<<"ready! Now saving...";

}

void generateHelix(qc::ScalarArray<double, qc::QC_3D> &Data, int N)
{
  cerr<<aol::color::red<<"\nGenerating helix...";

  double Drehungen = 4.;

  int x,y,z;
  double t = 0.;
  double n = static_cast<double> (N);
  // original: double r = 0.4 / ( (2. * Drehungen) * M_PI );
  double r = 0.3 / ( (2. * Drehungen) * M_PI );

  for (t=0.; t<(2. * Drehungen * M_PI); t+=0.01)
  {
    x = static_cast<int> (n * (r * t * sin(t) + 0.5));
    y = static_cast<int> (n * (r * t * cos(t) + 0.5));
    z = static_cast<int> (n * (t + 2.*M_PI) / ((3. * Drehungen + 1.) * M_PI));

    if ( z > 5 && z < N-5 ) drawLocalBall(Data, N, x, y, z, (t+3.)/(2.*Drehungen));   // 4. was previous t+3.

    if (x>N-1 || y>N-1 || z>N-1) cerr << "Something has gone wrong!!!";
  }
}

void generateHalfSpace(qc::ScalarArray<double, qc::QC_3D> &Data, int N, int xOffset)
{
  cerr<<aol::color::red<<"\nGenerating halfspace...";

  Data.setAll( 0. );
  for (int x=xOffset; x<N; x++)
    for (int y=0; y<N; y++)
      for (int z=0; z<N; z++)
        Data.set(y,z,x, 1.);
}

void generateDestroyedDoubleCone( qc::ScalarArray<double, qc::QC_3D> &Data, int N )
{
  Data.load( "Data/Wulffshapes/DoubleCone_PLUS5000_Tau002_1100.bz2" );
  double z0 = 0.2;
  double h  = 1. / static_cast<double>( N-1 );

  for ( int i=0; i<Data.getNumX(); i++) {
    for ( int j=0; j<Data.getNumY(); j++ ) {
      for ( int k=0; k<Data.getNumZ(); k++ ) {
        double z = k*h;
        if ( z < z0 || z > 1.-z0 /*|| Data.get( i,j,k ) > 0.*/ ) Data.set( i,j,k, 1. );
      }
    }
  }



}


void smoothData( qc::ScalarArray<double, qc::QC_3D> &Data, double rho, qc::GridDefinition &Grid )
{
  qc::LinearSmoothOp<double> linSmooth;
  linSmooth.setSigma( rho );
  linSmooth.setCurrentGrid (Grid);

  linSmooth.apply( Data, Data );
}


int main(int argc, char **argv)
{
  if ( argc != 3 ) {
    cerr << "USAGE: " << argv[0] << " depth savename" << endl;
    return EXIT_FAILURE;
  }
  try
  {
    // the depth of the grid
    int d = atoi( argv[1] );
    int N = (1<<d)+1;              // 2^N + 1
    qc::GridDefinition Grid( d, qc::QC_3D );

    qc::ScalarArray<double, qc::QC_3D> Data(N,N,N);
    cerr<<"Berechne skalare Daten...";


    // generate the curve and save it
    Data.setZero();

//     Data.setAll(-1.);

//     generateHelix( Data, N);
//     generateEllipsoidBlock(Data, N, N/4, N/8, N/8);
//     generateWrinkleCurve(Data, N, d);
//     generateAnyThickLine( Data, N, 32 );
//     deleteBlock( Data, 31,30,30, 33,34,34, 0. );
//     generateAnyThickLine( Data, N, 34 );
//     generateCurvedJunction( Data, N );
//     generateBall( Data, N, N/2, N/2, N/2, static_cast<double>(N/8.25) );
//     generateTorus( Data, N );
//     generateNoisyBall( Data, N );
//     generate2Balls( Data, N );
//     generateTubes( Data, N );
//     generateDistancePlane( Data, N );
//     generateXAxisLine( Data, N );
//     generateVesselPiece( Data, N );
//     generateNoisyStraitJunction( Data, N );
//     generateNoisySimpleArk(Data, N);
//     generateNoisySimpleTorus( Data, N );
//     generateCube(Data, N);
//     generateEllipsoid( Data, N, 1., 1., 1. );
//     generateAnyNoisyTubes( Data, N);
//     generateNoisyDiagLine( Data, N);
//     generateBlock(Data, N);
//     generateSphere( Data, N );
//     generateHalfSpace( Data, N, N/2 );
//     generateHorizontalLine( Data, N );
//     generateXAxisEllipsoidalLine( Data, N );

//     generateOval( Data, N, 0.05 );
//     generateEllipseStripe( Data, N, 0.9, 0.1, 0.45 );
//     generateL1Norm(Data, N);
    generateSupNorm(Data, N);
//     generateDestroyedDoubleCone( Data, N );
//     deleteBlock( Data, 4, 4, N/2 - 2, N-4, N-4, N/2 + 2, 1. );
    Data.addToAll ( - 0.1846153846153846154 );
//     generateSignedDistanceFunction( Data, N, d );
//     deleteBlock(Data, 5*N/8, 3*N/8, 3*N/8,                // (x_0,y_0,z_0), (x_1,y_1,z_1),
//                 7*N/8, 5*N/8, 5*N/8, 0.5);
//     Data -= 0.35;
//     smoothData( Data, 4.0 * Grid.H(), Grid );

//        Data -= 0.99999999999999999999;

//     generateSaltNPepaNoise( Data );
//     generateNoiseOnWholeData( Data, N );
//     generateNoiseAroundSupportDataBGNeg(Data, N);
//     generateNoiseAroundSupportData(Data, N);
//     Data *= 2.;

//     cerr << "Computing Signed Distance Function...";
//     generateSignedDistanceFunction( Data, N, d );

    // scale to [0,1]
//     Data -= Data.getMinValue();
//     Data /= Data.getMaxValue();
    cerr << "Saving...";
    Data.save( argv[2], qc::PGM_DOUBLE_BINARY );

    // ****************** Test the norms of the gradient **********************
//     qc::FindMinMaxNormOfGradientOp<ConfigType> FindGradientNorm( Grid );
//     aol::Vec2<double> maxGradientNorm;
//     FindGradientNorm.apply( Data, maxGradientNorm );
//     cerr << green << "Maximaler Gradient: "<< maxGradientNorm[1] << endl << reset;




    cerr<<aol::color::reset;
    // and view it in Grape
#ifdef USE_EXTERNAL_GRAPE
    GENMESH3D* gmeshNew = quocmesh_convert_to_gmesh3d(&Data, "generated data");

    // ************** now start GRAPE ************
    initStartGrape(gmeshNew, "generated data");
#else
    cerr << "cannot display in GRAPE without grape external" << endl;
#endif
    // thats it - try and enjoy it :-)
  }
  catch(aol::Exception e)
  {
    e.dump();
    return 42;
  }

  return 0;
}
