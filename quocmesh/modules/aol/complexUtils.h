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

#ifndef __COMPLEXUTILS_H
#define __COMPLEXUTILS_H

#include <aol.h>

namespace aol {

//! Mathematical constant i
const std::complex<double> I ( 0, 1 );

//! Allows int * complex<double>
template <class TypeA, class TypeB>
std::complex<TypeB> operator * ( TypeA a, std::complex<TypeB> b ) {
  return std::complex<TypeB> ( static_cast<TypeB> ( a ) ) * b;
}

//! Allows complex<double> * int
template <class TypeA, class TypeB>
std::complex<TypeA> operator * ( std::complex<TypeA> a, TypeB b ) {
  return a * std::complex<TypeA> ( static_cast<TypeA> ( b ) );
}

//! Allows int + complex<double>
template <class TypeA, class TypeB>
std::complex<TypeB> operator + ( TypeA a, std::complex<TypeB> b ) {
  return std::complex<TypeB> ( static_cast<TypeB> ( a ) ) + b;
}

//! Allows complex<double> + int
template <class TypeA, class TypeB>
std::complex<TypeA> operator + ( std::complex<TypeA> a, TypeB b ) {
  return a + std::complex<TypeA> ( static_cast<TypeA> ( b ) );
}

//! Allows int - complex<double>
template <class TypeA, class TypeB>
std::complex<TypeB> operator - ( TypeA a, std::complex<TypeB> b ) {
  return std::complex<TypeB> ( static_cast<TypeB> ( a ) ) - b;
}

//! Allows complex<double> - int
template <class TypeA, class TypeB>
std::complex<TypeA> operator - ( std::complex<TypeA> a, TypeB b ) {
  return a - std::complex<TypeA> ( static_cast<TypeA> ( b ) );
}

//! Allows int / complex<double>
template <class TypeA, class TypeB>
std::complex<TypeB> operator / ( TypeA a, std::complex<TypeB> b ) {
  return std::complex<TypeB> ( static_cast<TypeB> ( a ) ) / b;
}

//! Allows complex<double> / int
template <class TypeA, class TypeB>
std::complex<TypeA> operator / ( std::complex<TypeA> a, TypeB b ) {
  return a / std::complex<TypeA> ( static_cast<TypeA> ( b ) );
}

}

// Implicit use in operation of complex with int
using aol::operator*;
using aol::operator+;
using aol::operator-;
using aol::operator/;

#endif
