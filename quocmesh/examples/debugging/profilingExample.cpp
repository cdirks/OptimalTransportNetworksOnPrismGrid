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

/** \file
 *  \brief using gprof for profiling
 *
 *  compile whole library with LFLAGS += -pg, CFLAGS += -pg,
 *  ./profilingExample, gprof ./profilingExample Program with a very
 *  inefficient method to set vector entry which can be detected using
 *  gprof.
 *
 *  \author Schwen
 */

#include <aol.h>
#include <vec.h>
#include <randomGenerator.h>

// this is an extremely inefficient method
void setValue ( aol::Vector<double> &vec, const int index, const double value ) {
  aol::Vector<double> cpVec ( vec.size() );

  for ( int i = 0; i < vec.size(); ++i ) {
    if ( i == index )
      cpVec[i] = value;
    else
      cpVec[i] = vec[i];
  }

  vec = cpVec;

}

double findMaxValue ( aol::Vector<double> &vec ) {
  double maxSoFar = - aol::NumberTrait<double>::Inf;
  for ( int i = 0; i < vec.size(); ++i )
    if ( vec[i] > maxSoFar )
      maxSoFar = vec[i];

  return ( maxSoFar );
}

int main ( int, char**) {
  const int size = 100000;
  aol::Vector<double> Vect ( size );

  { // fill vector with some random data
    aol::RandomGenerator rg;
    for ( int i = 0; i < size; ++i )
      Vect[i] = rg.rReal<double>();
  }

  for ( int i = 0; i < 200; ++i ) {
    setValue ( Vect, i, sin ( 1.0 * i ) );
    cerr << findMaxValue ( Vect ) << endl;
  }

}
