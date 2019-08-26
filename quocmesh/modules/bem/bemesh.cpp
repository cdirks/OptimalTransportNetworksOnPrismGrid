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

#include <bemesh.h>

namespace {

// Ensure that this code is called in all programs linking this object
class Foo {
public:
  Foo () {
    qc::ScalarArray<double, qc::QC_2D>::quietMode = ! aol::debugging::arrayverbose;
  }
};

Foo foo;

}

int aol::ld ( int x ) {
  return static_cast<int> ( log ( static_cast<double> ( x ) + 0.5 ) / M_LN2 );
}

void aol::rotate ( aol::Matrix22<double>& mat, double d ) {
  aol::Matrix22<double> q ( aol::rotationMatrix ( d ) ), qt ( q ); qt.transpose ();
  aol::Matrix22<double> tmp ( q ); tmp *= mat; tmp *= qt;
  mat = tmp;
}

void aol::rotate ( aol::Vec2<double>& vec, double d ) {
  aol::Matrix22<double> q ( aol::rotationMatrix ( d ) );
  vec = q * vec;
}

void aol::rotate ( aol::Vector<double>& v, int p0, int p1, double d ) {
  int n = ( p1 - p0 ) / 2;
  for ( int i = 0; i < n; ++i ) {
    aol::Vec2<double> tmp ( v [p0+i], v [p0+n+i] );
    rotate ( tmp, d );
    v [p0+i] = tmp [0]; v [p0+n+i] = tmp [1];
  }
}

void aol::rotate ( aol::Vector<double>& v, int p, double d ) {
  rotate ( v, 0, p, d );
  rotate ( v, p, v.size (), d );
}
