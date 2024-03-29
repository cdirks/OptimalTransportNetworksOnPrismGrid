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

#include <smallVec.h>

namespace aol {

template <> const Vec2<float> ZTrait<Vec2<float> >::zero = Vec2<float> ( 0, 0 );
template <> const Vec2<double> ZTrait<Vec2<double> >::zero = Vec2<double> ( 0, 0 );
template <> const Vec2<long double> ZTrait<Vec2<long double> >::zero = Vec2<long double> ( 0, 0 );

template <> const Vec3<float> ZTrait<Vec3<float> >::zero = Vec3<float> ( 0, 0, 0 );
template <> const Vec3<double> ZTrait<Vec3<double> >::zero = Vec3<double> ( 0, 0, 0 );
template <> const Vec3<long double> ZTrait<Vec3<long double> >::zero = Vec3<long double> ( 0, 0, 0 );
template <> const Vec3<short> ZTrait<Vec3<short> >::zero = Vec3<short> ( 0, 0, 0 );
template <> const Vec3<uint64_t> ZTrait<Vec3<uint64_t> >::zero = Vec3<uint64_t> ( 0, 0, 0 );


template <> const Format* Vec<1, short>::_pFormat = &mixedFormat;
template <> const Format* Vec<1, unsigned short>::_pFormat = &mixedFormat;
template <> const Format* Vec<1, int>::_pFormat = &mixedFormat;
template <> const Format* Vec<1, unsigned int>::_pFormat = &mixedFormat;
template <> const Format* Vec<1, float>::_pFormat = &mixedFormat;
template <> const Format* Vec<1, double>::_pFormat = &mixedFormat;
template <> const Format* Vec<1, long double>::_pFormat = &mixedFormat;
template <> const Format* Vec<1, std::complex<float> >::_pFormat = &mixedFormat;
template <> const Format* Vec<1, std::complex<double> >::_pFormat = &mixedFormat;
template <> const Format* Vec<2, short>::_pFormat = &mixedFormat;
template <> const Format* Vec<2, unsigned short>::_pFormat = &mixedFormat;
template <> const Format* Vec<2, int>::_pFormat = &mixedFormat;
template <> const Format* Vec<2, unsigned int>::_pFormat = &mixedFormat;
template <> const Format* Vec<2, float>::_pFormat = &mixedFormat;
template <> const Format* Vec<2, double>::_pFormat = &mixedFormat;
template <> const Format* Vec<2, long double>::_pFormat = &mixedFormat;
template <> const Format* Vec<2, std::complex<float> >::_pFormat = &mixedFormat;
template <> const Format* Vec<2, std::complex<double> >::_pFormat = &mixedFormat;
template <> const Format* Vec<3, short>::_pFormat = &mixedFormat;
template <> const Format* Vec<3, unsigned short>::_pFormat = &mixedFormat;
template <> const Format* Vec<3, int>::_pFormat = &mixedFormat;
template <> const Format* Vec<3, float>::_pFormat = &mixedFormat;
template <> const Format* Vec<3, double>::_pFormat = &mixedFormat;
template <> const Format* Vec<3, long double>::_pFormat = &mixedFormat;
template <> const Format* Vec<3, uint64_t>::_pFormat = &mixedFormat;
template <> const Format* Vec<3, std::complex<float> >::_pFormat = &mixedFormat;
template <> const Format* Vec<3, std::complex<double> >::_pFormat = &mixedFormat;
template <> const Format* Vec<4, short>::_pFormat = &mixedFormat;
template <> const Format* Vec<4, unsigned short>::_pFormat = &mixedFormat;
template <> const Format* Vec<4, int>::_pFormat = &mixedFormat;
template <> const Format* Vec<4, float>::_pFormat = &mixedFormat;
template <> const Format* Vec<4, double>::_pFormat = &mixedFormat;
template <> const Format* Vec<4, long double>::_pFormat = &mixedFormat;
template <> const Format* Vec<4, std::complex<float> >::_pFormat = &mixedFormat;
template <> const Format* Vec<4, std::complex<double> >::_pFormat = &mixedFormat;
template <> const Format* Vec<5, short>::_pFormat = &mixedFormat;
template <> const Format* Vec<5, unsigned short>::_pFormat = &mixedFormat;
template <> const Format* Vec<5, int>::_pFormat = &mixedFormat;
template <> const Format* Vec<5, float>::_pFormat = &mixedFormat;
template <> const Format* Vec<5, double>::_pFormat = &mixedFormat;
template <> const Format* Vec<5, long double>::_pFormat = &mixedFormat;
template <> const Format* Vec<5, std::complex<float> >::_pFormat = &mixedFormat;
template <> const Format* Vec<5, std::complex<double> >::_pFormat = &mixedFormat;
template <> const Format* Vec<6, short>::_pFormat = &mixedFormat;
template <> const Format* Vec<6, unsigned short>::_pFormat = &mixedFormat;
template <> const Format* Vec<6, int>::_pFormat = &mixedFormat;
template <> const Format* Vec<6, float>::_pFormat = &mixedFormat;
template <> const Format* Vec<6, double>::_pFormat = &mixedFormat;
template <> const Format* Vec<6, long double>::_pFormat = &mixedFormat;
template <> const Format* Vec<6, std::complex<float> >::_pFormat = &mixedFormat;
template <> const Format* Vec<6, std::complex<double> >::_pFormat = &mixedFormat;
template <> const Format* Vec<8, short>::_pFormat = &mixedFormat;
template <> const Format* Vec<8, unsigned short>::_pFormat = &mixedFormat;
template <> const Format* Vec<8, int>::_pFormat = &mixedFormat;
template <> const Format* Vec<8, float>::_pFormat = &mixedFormat;
template <> const Format* Vec<8, double>::_pFormat = &mixedFormat;
template <> const Format* Vec<8, long double>::_pFormat = &mixedFormat;
template <> const Format* Vec<8, std::complex<float> >::_pFormat = &mixedFormat;
template <> const Format* Vec<8, std::complex<double> >::_pFormat = &mixedFormat;

template <> const Format* Vec<9, short>::_pFormat = &mixedFormat;
template <> const Format* Vec<9, unsigned short>::_pFormat = &mixedFormat;
template <> const Format* Vec<9, int>::_pFormat = &mixedFormat;
template <> const Format* Vec<9, float>::_pFormat = &mixedFormat;
template <> const Format* Vec<9, double>::_pFormat = &mixedFormat;
template <> const Format* Vec<9, long double>::_pFormat = &mixedFormat;
template <> const Format* Vec<9, std::complex<float> >::_pFormat = &mixedFormat;
template <> const Format* Vec<9, std::complex<double> >::_pFormat = &mixedFormat;

template <> const Format* Vec<10, short>::_pFormat = &mixedFormat;
template <> const Format* Vec<10, unsigned short>::_pFormat = &mixedFormat;
template <> const Format* Vec<10, int>::_pFormat = &mixedFormat;
template <> const Format* Vec<10, float>::_pFormat = &mixedFormat;
template <> const Format* Vec<10, double>::_pFormat = &mixedFormat;
template <> const Format* Vec<10, long double>::_pFormat = &mixedFormat;
template <> const Format* Vec<10, std::complex<float> >::_pFormat = &mixedFormat;
template <> const Format* Vec<10, std::complex<double> >::_pFormat = &mixedFormat;

template <> const Format* Vec<12, short>::_pFormat = &mixedFormat;
template <> const Format* Vec<12, unsigned short>::_pFormat = &mixedFormat;
template <> const Format* Vec<12, int>::_pFormat = &mixedFormat;
template <> const Format* Vec<12, float>::_pFormat = &mixedFormat;
template <> const Format* Vec<12, double>::_pFormat = &mixedFormat;
template <> const Format* Vec<12, long double>::_pFormat = &mixedFormat;
template <> const Format* Vec<12, std::complex<float> >::_pFormat = &mixedFormat;
template <> const Format* Vec<12, std::complex<double> >::_pFormat = &mixedFormat;

template <> bool Vec<1, short>::prettyFormat = true;
template <> bool Vec<1, unsigned short>::prettyFormat = true;
template <> bool Vec<1, int>::prettyFormat = true;
template <> bool Vec<1, unsigned int>::prettyFormat = true;
template <> bool Vec<1, float>::prettyFormat = true;
template <> bool Vec<1, double>::prettyFormat = false;
template <> bool Vec<1, long double>::prettyFormat = false;
template <> bool Vec<1, std::complex<float> >::prettyFormat = true;
template <> bool Vec<1, std::complex<double> >::prettyFormat = false;
template <> bool Vec<2, short>::prettyFormat = true;
template <> bool Vec<2, unsigned short>::prettyFormat = true;
template <> bool Vec<2, int>::prettyFormat = true;
template <> bool Vec<2, unsigned int>::prettyFormat = true;
template <> bool Vec<2, float>::prettyFormat = true;
template <> bool Vec<2, double>::prettyFormat = false;
template <> bool Vec<2, long double>::prettyFormat = false;
template <> bool Vec<2, std::complex<float> >::prettyFormat = true;
template <> bool Vec<2, std::complex<double> >::prettyFormat = false;
template <> bool Vec<3, short>::prettyFormat = true;
template <> bool Vec<3, unsigned short>::prettyFormat = true;
template <> bool Vec<3, int>::prettyFormat = true;
template <> bool Vec<3, float>::prettyFormat = true;
template <> bool Vec<3, double>::prettyFormat = false;
template <> bool Vec<3, long double>::prettyFormat = false;
template <> bool Vec<3, uint64_t>::prettyFormat = false;
template <> bool Vec<3, std::complex<float> >::prettyFormat = true;
template <> bool Vec<3, std::complex<double> >::prettyFormat = false;
template <> bool Vec<4, short>::prettyFormat = true;
template <> bool Vec<4, unsigned short>::prettyFormat = true;
template <> bool Vec<4, int>::prettyFormat = true;
template <> bool Vec<4, float>::prettyFormat = true;
template <> bool Vec<4, double>::prettyFormat = false;
template <> bool Vec<4, long double>::prettyFormat = false;
template <> bool Vec<4, std::complex<float> >::prettyFormat = true;
template <> bool Vec<4, std::complex<double> >::prettyFormat = false;
template <> bool Vec<5, short>::prettyFormat = true;
template <> bool Vec<5, unsigned short>::prettyFormat = true;
template <> bool Vec<5, int>::prettyFormat = true;
template <> bool Vec<5, float>::prettyFormat = true;
template <> bool Vec<5, double>::prettyFormat = false;
template <> bool Vec<5, long double>::prettyFormat = false;
template <> bool Vec<5, std::complex<float> >::prettyFormat = true;
template <> bool Vec<5, std::complex<double> >::prettyFormat = false;
template <> bool Vec<6, short>::prettyFormat = true;
template <> bool Vec<6, unsigned short>::prettyFormat = true;
template <> bool Vec<6, int>::prettyFormat = true;
template <> bool Vec<6, float>::prettyFormat = true;
template <> bool Vec<6, double>::prettyFormat = false;
template <> bool Vec<6, long double>::prettyFormat = false;
template <> bool Vec<6, std::complex<float> >::prettyFormat = true;
template <> bool Vec<6, std::complex<double> >::prettyFormat = false;
template <> bool Vec<8, short>::prettyFormat = true;
template <> bool Vec<8, unsigned short>::prettyFormat = true;
template <> bool Vec<8, int>::prettyFormat = true;
template <> bool Vec<8, float>::prettyFormat = true;
template <> bool Vec<8, double>::prettyFormat = false;
template <> bool Vec<8, long double>::prettyFormat = false;
template <> bool Vec<8, std::complex<float> >::prettyFormat = true;
template <> bool Vec<8, std::complex<double> >::prettyFormat = false;

template <> bool Vec<9, short>::prettyFormat = true;
template <> bool Vec<9, unsigned short>::prettyFormat = true;
template <> bool Vec<9, int>::prettyFormat = true;
template <> bool Vec<9, float>::prettyFormat = true;
template <> bool Vec<9, double>::prettyFormat = false;
template <> bool Vec<9, long double>::prettyFormat = false;
template <> bool Vec<9, std::complex<float> >::prettyFormat = true;
template <> bool Vec<9, std::complex<double> >::prettyFormat = false;

template <> bool Vec<10, short>::prettyFormat = true;
template <> bool Vec<10, unsigned short>::prettyFormat = true;
template <> bool Vec<10, int>::prettyFormat = true;
template <> bool Vec<10, float>::prettyFormat = true;
template <> bool Vec<10, double>::prettyFormat = false;
template <> bool Vec<10, long double>::prettyFormat = false;
template <> bool Vec<10, std::complex<float> >::prettyFormat = true;
template <> bool Vec<10, std::complex<double> >::prettyFormat = false;

template <> bool Vec<12, short>::prettyFormat = true;
template <> bool Vec<12, unsigned short>::prettyFormat = true;
template <> bool Vec<12, int>::prettyFormat = true;
template <> bool Vec<12, float>::prettyFormat = true;
template <> bool Vec<12, double>::prettyFormat = false;
template <> bool Vec<12, long double>::prettyFormat = false;
template <> bool Vec<12, std::complex<float> >::prettyFormat = true;
template <> bool Vec<12, std::complex<double> >::prettyFormat = false;

} // end namespace
