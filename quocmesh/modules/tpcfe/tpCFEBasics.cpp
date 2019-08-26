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

#include <aol.h>
#include <tpCFEBasics.h>
#include <tpCFEUtils.h>

namespace tpcfe {

template < typename RealType >
RealType IsotropicElasticityCoefficient<RealType>::getLambda ( ) const {
  return tpcfe::computeLambda ( getE(), getNu() );
}

template < typename RealType >
RealType IsotropicElasticityCoefficient<RealType>::getMu ( ) const {
  return tpcfe::computeMu ( getE(), getNu() );
}

template < typename RealType >
void IsotropicElasticityCoefficient<RealType>::getAnisotropicTensor ( RealType Tensor[3][3][3][3] ) const {
  tpcfe::setLameNavierTensorFromLambdaMu ( getLambda(), getMu(), Tensor );
}

template < typename RealType >
void IsotropicElasticityCoefficient<RealType>::print ( ostream &os ) const {
  os << "E = " << _E << ", nu = " << _nu << endl;
}

//! OS is not quite sure whether this type of averaging makes sense
template < typename RealType >
IsotropicElasticityCoefficient<RealType>& IsotropicElasticityCoefficient<RealType>::operator+= ( const IsotropicElasticityCoefficient<RealType> &other ) {
  _E += other._E;
  _nu += other._nu;
  return ( *this );
}

template < typename RealType >
IsotropicElasticityCoefficient<RealType>& IsotropicElasticityCoefficient<RealType>::operator*= ( const RealType value ) {
  _E *= value;
  _nu *= value;
  return ( *this );
}

template < typename RealType >
IsotropicElasticityCoefficient<RealType>& IsotropicElasticityCoefficient<RealType>::operator/= ( const RealType value ) {
  _E /= value;
  _nu /= value;
  return ( *this );
}


template < typename RealType >
const short VoigtElasticityCoefficient<RealType>::indexM[3][3] = { { 0, 5, 4 },
                                                                   { 5, 1, 3 },
                                                                   { 4, 3, 2 } };

template < typename RealType >
void VoigtElasticityCoefficient<RealType>::getAnisotropicTensor ( RealType Tensor[3][3][3][3] ) const {
  for ( short i = 0; i < 3; ++i )
    for ( short j = 0; j < 3; ++j )
      for ( short k = 0; k < 3; ++k )
        for ( short l = 0; l < 3; ++l )
          Tensor[i][j][k][l] = _VoigtTensor[ indexM[i][j] ][ indexM[k][l] ];
}

template < typename RealType >
void VoigtElasticityCoefficient<RealType>::averageFromAnisotropicTensor ( const RealType Tensor[3][3][3][3] ) {
  _VoigtTensor.setZero();
  aol::Mat<6,6,RealType> denom;

  for ( short i = 0; i < 3; ++i ) {
    for ( short j = 0; j < 3; ++j ) {
      for ( short k = 0; k < 3; ++k ) {
        for ( short l = 0; l < 3; ++l ) {
          _VoigtTensor[ indexM[i][j] ][ indexM[k][l] ] += Tensor[i][j][k][l];
          denom[ indexM[i][j] ][  indexM[k][l] ] += 1.0;
        }
      }
    }
  }

  // averageing:
  for ( int i = 0; i < 6; ++i ) {
    for ( int j = 0; j < 6; ++j ) {
      _VoigtTensor[i][j] /= denom[i][j];
    }
  }
}

template < typename RealType >
void VoigtElasticityCoefficient<RealType>::rotateTensorBy ( const aol::Matrix33<RealType> rotMat ) {
  RealType tens[3][3][3][3] = {}, rotTens[3][3][3][3] = {}; // all-zero initialization
  this->getAnisotropicTensor ( tens );

  for ( short i = 0; i < 3; ++i )
    for ( short j = 0; j < 3; ++j )
      for ( short k = 0; k < 3; ++k )
        for ( short l = 0; l < 3; ++l )
          for ( short m = 0; m < 3; ++m )
            for ( short n = 0; n < 3; ++n )
              for ( short o = 0; o < 3; ++o )
                for ( short p = 0; p < 3; ++p )
                  rotTens[m][n][o][p] += rotMat[m][i] * rotMat[n][j] * rotMat[o][k] * rotMat[p][l] * tens[i][j][k][l];


  this->averageFromAnisotropicTensor ( rotTens );
}


template < typename RealType >
void VoigtElasticityCoefficient<RealType>::print ( ostream &os ) const {
  os << "Tensor in Voigt's Notation:" << endl << _VoigtTensor << endl;
}

template < typename RealType >
void VoigtElasticityCoefficient<RealType>::setFromENu ( const RealType E, const RealType nu ) {
  IsotropicElasticityCoefficient<RealType> isoTensor ( E, nu );
  RealType tens[3][3][3][3] = {};
  isoTensor.getAnisotropicTensor ( tens );
  this->averageFromAnisotropicTensor ( tens );
}


template < typename RealType >
VoigtElasticityCoefficient<RealType>& VoigtElasticityCoefficient<RealType>::operator+= ( const VoigtElasticityCoefficient<RealType> &other ) {
  _VoigtTensor += other._VoigtTensor;
  return ( *this );
}

template < typename RealType >
VoigtElasticityCoefficient<RealType>& VoigtElasticityCoefficient<RealType>::operator*= ( const RealType value ) {
  _VoigtTensor *= value;
  return ( *this );
}

template < typename RealType >
VoigtElasticityCoefficient<RealType>& VoigtElasticityCoefficient<RealType>::operator/= ( const RealType value ) {
  _VoigtTensor /= value;
  return ( *this );
}


template class IsotropicElasticityCoefficient<float>;
template class IsotropicElasticityCoefficient<double>;
template class IsotropicElasticityCoefficient<long double>;

template class VoigtElasticityCoefficient<float>;
template class VoigtElasticityCoefficient<double>;
template class VoigtElasticityCoefficient<long double>;

template<> const float        CFEStructureTrilin<      float>::_threshold = 1.0e-7;
template<> const double       CFEStructureTrilin<     double>::_threshold = 1.0e-15;
template<> const long double  CFEStructureTrilin<long double>::_threshold = 1.0e-15;

}

namespace aol {
template<> const tpcfe::IsotropicElasticityCoefficient<float>       aol::ZTrait< tpcfe::IsotropicElasticityCoefficient<float> >::zero       = tpcfe::IsotropicElasticityCoefficient<float>();
template<> const tpcfe::IsotropicElasticityCoefficient<double>      aol::ZTrait< tpcfe::IsotropicElasticityCoefficient<double> >::zero      = tpcfe::IsotropicElasticityCoefficient<double>();
template<> const tpcfe::IsotropicElasticityCoefficient<long double> aol::ZTrait< tpcfe::IsotropicElasticityCoefficient<long double> >::zero = tpcfe::IsotropicElasticityCoefficient<long double>();

template<> const tpcfe::VoigtElasticityCoefficient<float>       aol::ZTrait< tpcfe::VoigtElasticityCoefficient<float> >::zero       = tpcfe::VoigtElasticityCoefficient<float>();
template<> const tpcfe::VoigtElasticityCoefficient<double>      aol::ZTrait< tpcfe::VoigtElasticityCoefficient<double> >::zero      = tpcfe::VoigtElasticityCoefficient<double>();
template<> const tpcfe::VoigtElasticityCoefficient<long double> aol::ZTrait< tpcfe::VoigtElasticityCoefficient<long double> >::zero = tpcfe::VoigtElasticityCoefficient<long double>();
}
