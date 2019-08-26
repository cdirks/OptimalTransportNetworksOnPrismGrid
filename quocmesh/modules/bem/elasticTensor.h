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

#ifndef __ELASTICTENSOR_H
#define __ELASTICTENSOR_H

#include <complexUtils.h>
#include <tensor.h>
#include <bemesh.h>
#include <typedParser.h>

namespace bm {
//! Generic implementation of Green function
//! Base class for isotropic and anistropic case
template <class DataType, class GreenType>
struct GenElasticGreen {
  aol::ElasticTensor<DataType> C;
  GenElasticGreen ( );
  GenElasticGreen ( aol::ElasticTensor<DataType> c );
  GreenType&       asImp ()       { return * static_cast<GreenType*> (this); }
  const GreenType& asImp () const { return * static_cast<const GreenType*> (this); }
};

//! Original anisotrpic implementation by ML
//! @todo ML/BG first ctor called with isoOffset=0 and c not satisfying anisotropy criterion results in an endless loop
template <class DataType>
struct ElasticGreen : GenElasticGreen< DataType, ElasticGreen<DataType> > {
  aol::Vec2<complex<DataType> > p;
  aol::Matrix22<complex<DataType> > A;
  aol::Tensor222<complex<DataType> > L;
  aol::Matrix22<complex<DataType> > N;
  aol::Matrix22<DataType> d;
  ElasticGreen ( aol::ElasticTensor<DataType> c, DataType isoOffset = 0 );
  ElasticGreen ( const aol::TypedParameterParser& parameter, string phasename );
};

//! Additional isotropic implementation by BG
//! @todo BG this should be removed
template <class DataType>
struct IsoElasticGreen : GenElasticGreen< DataType, IsoElasticGreen<DataType> > {
  DataType lambda;
  DataType mu;
  DataType c1;
  DataType c2;
  IsoElasticGreen ( aol::ElasticTensor<DataType> c );
};

// GenElasticGreen
template <class DataType, class GreenType>
GenElasticGreen<DataType, GreenType>::GenElasticGreen ( )
    : C ( 0, 0, 0 ) {
}

template <class DataType, typename GreenType>
GenElasticGreen<DataType, GreenType>::GenElasticGreen ( aol::ElasticTensor<DataType> c )
    : C ( c ) {
}

// ElasticGreen
template <class DataType>
ElasticGreen<DataType>::ElasticGreen ( const aol::TypedParameterParser& parameter, string phasename )
    : GenElasticGreen<DataType, ElasticGreen> ( ) {
  aol::ElasticTensor<DataType> c ( 0, 0, 0 );
  double isoOffset;
  parameter.get ( phasename, c );
  parameter.get ( "isoOffset", isoOffset );
  ElasticGreen<DataType> green ( c, isoOffset );
  ( *this ) = green;
}

template <class DataType>
ElasticGreen<DataType>::ElasticGreen ( aol::ElasticTensor<DataType> c, DataType isoOffset )
    : GenElasticGreen<DataType, ElasticGreen> ( c ) {

  // this class was written for orthotropic materials in reference configuration, check if this is true for current c
  if (  aol::Abs(c[0][0][0][1])+aol::Abs(c[0][0][1][0])+aol::Abs(c[0][1][0][0])+aol::Abs(c[1][0][0][0])
       +aol::Abs(c[1][1][1][0])+aol::Abs(c[1][1][0][1])+aol::Abs(c[1][0][1][1])+aol::Abs(c[0][1][1][1]) > 1e-8 )
    throw aol::UnimplementedCodeException ( "ElasticGreen: fundamental solution only implemented for aligned orthotropic tensors", __FILE__, __LINE__ );


  int i, j, k, r, alpha;

  this->C.isoOffset ( isoOffset );

  DataType x = ( this->C.C [1][1] * this->C.C [0][0] - this->C.C [0][1] * this->C.C [0][1] - 2 * this->C.C [0][1] * this->C.C44 )
               / ( this->C.C44 * this->C.C [1][1] );
  DataType diskr = 0.25 * x * x - this->C.C [0][0] / this->C.C [1][1];
  complex<DataType> cdiskr = diskr;

  if ( diskr >= 0 ) {
    p [0] = sqrt ( -0.5 * x + sqrt ( cdiskr ) );
    if ( imag ( p [0] ) < 0 )
      p [0] = conj ( p [0] );
    p [1] = sqrt ( -0.5 * x - sqrt ( cdiskr ) );
    if ( imag ( p [1] ) < 0 )
      p [1] = conj ( p [1] );
  } else {
    p [0] = sqrt ( -0.5 * x + sqrt ( cdiskr ) );
    if ( imag ( p [0] ) < 0 )
      p [0] = conj ( p [0] );
    p [1] = -conj ( p [0] );
  }

  for ( alpha = 0; alpha < 2; alpha++ ) {
    A [0][alpha] = this->C [0][0][1][0] + ( this->C [0][0][1][1] + this->C [0][1][1][0] ) * p [alpha] + this->C [0][1][1][1] * p [alpha] * p [alpha];
    A [1][alpha] = - ( this->C [0][0][0][0] + ( this->C [0][0][0][1] + this->C [0][1][0][0] ) * p [alpha] + this->C [0][1][0][1] * p [alpha] * p [alpha] );
  }

  N = A.inverse ();

  for ( i = 0; i < 2; i++ )
    for ( j = 0; j < 2; j++ )
      for ( alpha = 0; alpha < 2; alpha++ ) {
        L [i][j][alpha] = 0;
        for ( k = 0; k < 2; k++ )
          L [i][j][alpha] += ( this->C [i][j][k][0] + p [alpha] * this->C [i][j][k][1] ) * A [k][alpha];
      }

  aol::Matrix22<DataType> Z;

  for ( i = 0; i < 2; i++ )
    for ( r = 0; r < 2; r++ )
      for ( Z [i][r] = 0, alpha = 0; alpha < 2; alpha++ )
        Z [i][r] += real ( N [alpha][r] * ( L [i][0][alpha] / p [alpha] * ( log ( 1 + p[alpha] ) - log ( 1 - p[alpha] ) ) -
                                            L [i][1][alpha] * ( log ( 1 - p[alpha] ) - log ( -1 - p[alpha] ) ) ) ) / aol::NumberTrait<long double>::pi;

  d = Z.inverse ();
}

// IsoElasticGreen
template <class DataType>
IsoElasticGreen<DataType>::IsoElasticGreen ( aol::ElasticTensor<DataType> c )
    : GenElasticGreen<DataType, IsoElasticGreen> ( c ) {
  // Check isotropy
  if (( this->C.C11 != this->C.C22 ) || ( this->C.C11 != this->C.C12 + 2*this->C.C44 ))
    throw ( aol::Exception ( "  !!! IsoElasticTensor is not isotropic !!!", __FILE__, __LINE__ ) );
  lambda = this->C.C12;
  mu     = this->C.C44;
  c1     = ( lambda + mu ) / ( 4*aol::NumberTrait<long double>::pi * ( mu * ( lambda + 2*mu) ) );
  c2     = - ( lambda + 3*mu ) / ( lambda + mu );
}

}

#endif
