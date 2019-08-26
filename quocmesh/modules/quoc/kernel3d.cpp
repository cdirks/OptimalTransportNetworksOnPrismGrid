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

#include <kernel3d.h>
#include <estimator2d.h>
#include <matrix.h>
#include <quoc.h>


template<typename RealType>
void qc::GaussKernel3d<RealType>::makeKernel() {
  for ( int X = -this->offset; X <= this->offset; X++ ) {
    for ( int Y = -this->offset; Y <= this->offset; Y++ ) {
      for ( int Z = -this->offset; Z <= this->offset; Z++ ) {

        RealType val = exp ( -0.5 * ( X * X + Y * Y + Z * Z ) / ( sigma * sigma ) );

        setValue ( X, Y, Z, val );

      }
    }
  }
  this->normalize();
}


template<typename RealType>
void qc::GaussDiffKernel3d<RealType>::makeKernel() {

  RealType s = 2 * sigma * sigma;

  for ( int X = -this->offset; X <= this->offset; X++ ) {
    for ( int Y = -this->offset; Y <= this->offset; Y++ ) {
      for ( int Z = -this->offset; Z <= this->offset; Z++ ) {
        RealType val = exp ( -0.5f * ( X * X + Y * Y + Z * Z ) / ( sigma * sigma ) );
        switch ( comp ) {
        case DIFF_X:
          val *= static_cast< RealType > ( X ) / s;
          break;
        case DIFF_Y:
          val *= static_cast< RealType > ( Y ) / s;
          break;
        case DIFF_Z:
          val *= static_cast< RealType > ( Z ) / s;
          break;
        case DIFF_XX:
          val = val * aol::Sqr ( ( static_cast< RealType > ( X ) ) / s ) + val / s;
          break;
        case DIFF_YY:
          val = val * aol::Sqr ( ( static_cast< RealType > ( Y ) ) / s ) + val / s;
          break;
        case DIFF_ZZ:
          val = val * aol::Sqr ( ( static_cast< RealType > ( Y ) ) / s ) + val / s;
          break;
        case DIFF_YZ:
          val *= static_cast< RealType > ( Z ) * static_cast< RealType > ( Y ) / aol::Sqr ( s );
          break;
        case DIFF_XZ:
          val *= static_cast< RealType > ( Z ) * static_cast< RealType > ( X ) / aol::Sqr ( s );
          break;
        case DIFF_XY:
          val *= static_cast< RealType > ( Y ) * static_cast< RealType > ( X ) / aol::Sqr ( s );
          break;
        default:
          throw aol::UnimplementedCodeException ( "Unknown qc::DiffVarType", __FILE__, __LINE__ );
        }
        this->setValue ( X, Y, Z, val );
      }
    }
  }
}

template class qc::GaussDiffKernel3d<long double>;
template class qc::GaussDiffKernel3d<double>;
template class qc::GaussDiffKernel3d<float>;
