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
#include <quoc.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <preconditioner.h>

#include <qmException.h>
#include <aol.h>

// #include <anisotropies.h>

#include <math.h>
#include <polarCoords.h>

#define REAL double

double gamma( aol::Vec3<double> z ) {
  return( 100.*100.*z[0]*z[0] + z[1]*z[1] + z[2]*z[2] );
}

double supGammaz( aol::Vec3<double> z ) {
  aol::Vec3<double> gammaZ( 100.*z[0], z[1], z[2] );
  gammaZ /= gamma(z);
  return gammaZ.infinityNorm();
}

double supGammazz( aol::Vec3<double> z ) {
  double g = gamma(z);
  aol::Matrix33<double> Gammazz( 100.*g - aol::Sqr( 100.*z[0] )/g,  - 100.*z[0]*z[1]/g, - 100.*z[0]*z[2]/g,
                                 -z[0]*z[1]/g, g-z[1]*z[1]/g, -z[2]*z[1]/g,
                                 -z[0]*z[2]/g, -z[1]*z[2]/g, g-z[2]*z[2]/g );
  Gammazz /= g*g;
  return Gammazz.infinityNorm();
}

double absMax( double a, double b ) {
  if ( aol::Abs(a) > aol::Abs(b) ) return aol::Abs(a);
  else return aol::Abs(b);
}

int main( int /*argc*/, char **/*argv*/)
{
  try
    {
      aol::Vec3<double> z;

      double gammaInf = 10000000000.;
      double gammaZSup = 0.;

      double g,gz,gzz;

      double step = 0.001;
      for ( double alpha = - M_PI/2.; alpha < M_PI/2.; alpha += step ) {
        for ( double beta = 0.; beta < 2.*M_PI; beta += step ) {
          z[0] = cos(beta)*cos(alpha);
          z[1] = cos(beta)*sin(alpha);
          z[2] = sin(beta);

          g = gamma( z );
          gz = supGammaz( z );
          gzz = supGammazz( z );

          if ( g < gammaInf ) gammaInf = g;

          if ( absMax( gz, gzz ) > gammaZSup ) gammaZSup = absMax( gz, gzz );
        }
      }

      cerr << "Inf( gamma ): " << gammaInf << endl;
      cerr << "Sup( gamma_z, gamma_zz ): " << gammaZSup << endl;
      cerr << "Faktor * Sup(...): " << 1./(sqrt(5.)-1)*gammaZSup << endl;
      cerr << "=> lambda > " <<  ( 1./(sqrt(5.)-1)*gammaZSup ) / gammaInf << endl;

//       aol::SphericalVec<double> n[8];
//
//       n[0].initByCartesianCoords( 1.,1.,1. );
//       n[1].initByCartesianCoords( 1.,1.,-1. );
//       n[2].initByCartesianCoords( 1.,-1.,1. );
//       n[3].initByCartesianCoords( 1.,-1.,-1. );
//       n[4].initByCartesianCoords( -1.,1.,1. );
//       n[5].initByCartesianCoords( -1.,1.,-1. );
//       n[6].initByCartesianCoords( -1.,-1.,1. );
//       n[7].initByCartesianCoords( -1.,-1.,-1. );
//
//       aol::Vec3<double> cartCoords;
//
//       for (int i=0; i<8; i++)
//       {
//         n[i].print();
//         cartCoords = n[i].getCartesianCoords();
//         cerr<<", kartesische Koordinaten: "<<cartCoords<<endl;
//       }
//
//       aol::SphericalVec<double> v1(1,2,3);
//       aol::SphericalVec<double> v2(1,2,3);
//       aol::SphericalVec<double> v3(v2);
//       cerr<<"\nv3: ";
//       v3.print();
//
//       v2.setZero();
//       v2.print();
//
//       v3 = v2;
//       cerr<<"\nv3 (sollte 0 sein): ";
//       v3.print();
//
//       v3.setRadius(3.); v3.setTheta(4.); v3.setPhi(5.);
//       cerr<<"\nv3 (sollte 3,4,5 sein): ";
//       v3.print();
//       v3.normalize();v3.print();
//
//       v3.initByCartesianCoords( cartCoords );
//       cerr<<"\nv3 cart: ";
//       v3.print();
//
//       cerr<<"\nv2=v1? "<<(v1 == v2)<<endl;
//       cerr<<"\nv2!=v1? "<<(v1 != v2)<<endl;
//
// //       cerr<<"atan2(y,x): "<<endl;
// //       cerr<<"arctan(1/1): "<<atan2(1.,1.)<<", arctan(-1,1): "<<atan2(-1.,1.)<<endl;
// //       cerr<<"arctan(1/-1): "<<atan2(1.,-1.)<<", arctan(-1,-1): "<<atan2(-1.,-1.)<<endl;
// //       cerr<<"atan()"<<endl;
// //       cerr<<"arctan(1): "<<atan(1.)<<", arctan(-1): "<<atan(-1.)<<endl;
//
//       cerr<<"\nasin()"<<endl;
//       cerr<<"arcsin(1): "<<asin(1.)<<", arcsin(-1): "<<asin(-1.)<<endl;

//       cone3d<double> kegel(1., 1., 0.001);
//       cerr << "gamma( "<<n1<<"): "<<kegel.gamma(n1)<<endl;
//       cerr << "gamma( "<<n2<<"): "<<kegel.gamma(n2)<<endl;
//       cerr << "gamma( "<<n3<<"): "<<kegel.gamma(n3)<<endl;
//       cerr << "gamma( "<<n4<<"): "<<kegel.gamma(n4)<<endl;
//       cerr << "gamma( "<<n5<<"): "<<kegel.gamma(n5)<<endl;
//       cerr << "gamma( "<<n6<<"): "<<kegel.gamma(n6)<<endl;
//       cerr << "gamma( "<<n7<<"): "<<kegel.gamma(n7)<<endl;
//       cerr << "gamma( "<<n8<<"): "<<kegel.gamma(n8)<<endl;
//       cerr << "gamma( "<<n9<<"): "<<kegel.gamma(n9)<<endl;
    }
  catch ( aol::Exception &ex )
    {
      ex.dump();
    }
}
