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

#include <l2Projector.h>

template <typename RealType>
aol::ExitStatus qc::Basis2DConst<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DConst<RealType>* > ( &B ) )  {
    Prod = 1.0;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DX<RealType>* > ( &B ) )  {
    Prod = 0.5;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DY<RealType>* > ( &B ) )  {
    Prod = 0.5;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXY<RealType>* > ( &B ) ) {
    Prod = 0.25;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXX<RealType>* > ( &B ) ) {
    Prod = 1. / 3.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 3.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template <typename RealType>
aol::ExitStatus qc::Basis2DX<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DX<RealType>* > ( &B ) )  {
    Prod = 1. / 3.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DY<RealType>* > ( &B ) )  {
    Prod = 0.25;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXY<RealType>* > ( &B ) ) {
    Prod = 1. / 6.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXX<RealType>* > ( &B ) ) {
    Prod = 1. / 4.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 6.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template <typename RealType>
aol::ExitStatus qc::Basis2DY<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DY<RealType>* > ( &B ) )  {
    Prod = 1. / 3.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXY<RealType>* > ( &B ) ) {
    Prod = 1. / 6.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXX<RealType>* > ( &B ) ) {
    Prod = 1. / 6.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 4.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template <typename RealType>
aol::ExitStatus qc::Basis2DXY<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DXY<RealType>* > ( &B ) )  {
    Prod = 1. / 9.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXX<RealType>* > ( &B ) ) {
    Prod = 1. / 8.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 8.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template <typename RealType>
aol::ExitStatus qc::Basis2DXX<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DXX<RealType>* > ( &B ) ) {
    Prod = 1. / 5.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 9.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template <typename RealType>
aol::ExitStatus qc::Basis2DYY<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 5.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template class qc::Basis2DConst<float>;
template class qc::Basis2DX<float>;
template class qc::Basis2DY<float>;
template class qc::Basis2DXY<float>;
template class qc::Basis2DXX<float>;
template class qc::Basis2DYY<float>;

template class qc::L2Projector<float>;


template class qc::Basis2DConst<double>;
template class qc::Basis2DX<double>;
template class qc::Basis2DY<double>;
template class qc::Basis2DXY<double>;
template class qc::Basis2DXX<double>;
template class qc::Basis2DYY<double>;

template class qc::L2Projector<double>;
