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

#ifndef __STATS_H
#define __STATS_H

#include <vec.h>

namespace aol {

template <typename RealType>
RealType MSE ( const aol::Vector<RealType> &Truth, const aol::Vector<RealType> &Estimate ) {
  if ( Truth.size ( ) != Estimate.size ( ) )
    throw aol::Exception ( "Vector of ground truths and estimates must have the same lengths!", __FILE__, __LINE__ );
  
  RealType res = 0;
  for ( int k=0; k<Truth.size ( ) ; ++k )
    res += aol::Sqr<RealType> ( Estimate[k] - Truth[k] );
  return res / Truth.size ( );
}

template <typename RealType>
RealType PSNR ( const aol::Vector<RealType> &Truth, const aol::Vector<RealType> &Estimate, const RealType FixedMaxValue = 0 ) {
  const RealType peak = ( FixedMaxValue > 0 ) ? FixedMaxValue : Truth.getMaxValue ( );
  return 10 * log10 ( aol::Sqr<RealType> ( peak ) / MSE ( Truth, Estimate ) );
}

}

#endif // __STATS_H
