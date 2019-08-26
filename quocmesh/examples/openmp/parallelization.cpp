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
 *  \brief usage of openmp for parallelization
 *
 *   Illustrates how to use openmp for parallelization. This requires
 *   -fopenmp (compiler flag) and -lgomp (linker flag)
 *
 *  \author Schwen
 */

#include <aol.h>

// #define _OPENMP 1
// this define should already set automatically if the compiler uses openmp

int main ( int, char** ) {

  const int size = 1 << 25;

  double *vecA = new double[size], *vecB = new double[size];

  aol::StopWatch timer;

  {
    // first, let's fill the vectors with some data:
    timer.start();
    for ( int k=0; k<20; ++k ) {
      cerr << k << ",";
      for ( int i = 0; i < size; ++i ){
        vecA[i] = cos ( 1.1 * i );
        vecB[i] = cos ( 0.9 * i );
      }
    }
    timer.stop();
    cerr << "Non-parallel filling of vectors took (CPU-time) " << timer.elapsedCpuTime() << endl;
    cerr << "Non-Parallel filling of vectors took (Wallclock-time) " << timer.elapsedWallClockTime() << endl;
  }

  // But why not fill the vectors in parallel?
  timer.start();
  // This is how we split in two threads:
#ifdef _OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef _OPENMP
#pragma omp section
#endif
    {
      for ( int k=0; k<20; ++k ) {
        cerr << k << ",";
        for ( int i = 0; i < size; ++i ){
          vecA[i] = cos ( 1.1 * i );
        }
      }
    }

#ifdef _OPENMP
#pragma omp section
#endif
    {
      for ( int k=0; k<20; ++k ) {
        cerr << k << ",";
        for ( int i = 0; i < size; ++i ){
          vecB[i] = cos ( 0.9 * i );
        }
      }
    }
  }
  timer.stop();

  // the threads should rejoin automatically (when results are used, automatically wait for slower thread?)
  {
    double tmp = vecA[size-1] + vecB[size-1];    tmp += 1.0;  // use  the result
  }

  cerr << "Parallel filling of vectors took (CPU-time) " << timer.elapsedCpuTime() << endl;
  cerr << "Parallel filling of vectors took (Wallclock-time) " << timer.elapsedWallClockTime() << endl;


  // Next, let's copy vectors.
  timer.start();
  for( int i = 0; i < size; i++ ) {
    vecA[i] = vecB[i];
  }
  timer.stop();
  cerr << "Non-parallel copying of vectors took " << timer.elapsedCpuTime() << endl;


  // Write access is independent of other iterations, so we can parallelize as follows:
  timer.start();
#ifdef _OPENMP
#pragma omp parallel for shared(vecA, vecB)
#endif
  for( int i = 0; i < size; i++ ) {
    vecA[i] = vecB[i];
  }

  timer.stop();
  cerr << "Parallel copying of vectors took " << timer.elapsedCpuTime() << endl;



  // Finally, let's do something similar to computing a vector norm, but using a cosine evaulation that should take some time:
  {
    timer.start();
    double norm = 0.0;
    for( int i = 0; i < size; i++ ) {
      norm += cos ( vecA[i] * vecA[i] );
    }
    timer.stop();
    vecB[0] = norm; // use  the result
    cerr << "Non-parallel vector cos-norm took " << timer.elapsedCpuTime() << endl;
  }

  // Here, write access is not independent but the result is independent of the order in which the loop is performed and independent of possible splitting of the loop.
  {
    timer.start();
    double norm = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: norm)
#endif
    for( int i = 0; i < size; i++ ) {
      norm += cos ( vecA[i] * vecA[i] );
    }

    timer.stop();
    vecB[0] = norm; // use  the result
    cerr << "Parallel vector cos-norm took " << timer.elapsedCpuTime() << endl;
  }

  delete[] ( vecA );
  delete[] ( vecB );

  return( EXIT_SUCCESS );
}
