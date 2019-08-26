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

#include <iostream>
#include <iomanip>
#include <xmmintrin.h>
#include <vector>

static const int N = 40000; // vector lenght
static const int M = 10000;  // number of scalar products to be calculated

using namespace std;

int main( int, char ** ) {
  // compiler flag -funroll-loops decreased the speed advantage of SSE from 3.7x to 1.8x
  // this could be because the loop over j is really unnecessary
  // This speed advantage of -funroll-loops is killed by the following line:
  std::vector<float*> arrays1(0);
  __attribute__((aligned(16))) float array1[N];
  __attribute__((aligned(16))) float array2[N];

  for( int i = 0; i < N; i++ ) {
    array1[i] = 0.5;
    array2[i] = 2.0;
  }

  /*
    __declspec(align(16)) float array1[N];
    __declspec(align(16)) float array2[N];
    __declspec(align(16)) float result[N];
  */
  double runTimeSSE = 1.;
  double runTimeNoSSE = 1.;
  {
    clock_t start = clock();
    const int nLoop = N/4;

    float result;
    for( int j = 1; j < M; j++) {
      __attribute__((aligned(16))) float temp[4];

      __m128 m1;
      __m128* m2 = reinterpret_cast< __m128* >( temp );

      __m128* pSrc1 = reinterpret_cast< __m128* >( array1 );
      __m128* pSrc2 = reinterpret_cast< __m128* >( array2 );
      *m2 = _mm_set_ps1(0.0f); // set all four values of m2 to 0.

      for ( int i = 0; i < nLoop; i++ ) {
        // This code is easy to understand:
        m1 = _mm_mul_ps(*pSrc1, *pSrc2);          // *pDest = *pSrc1 * *pSrc2
        *m2 = _mm_add_ps(m1, *m2);

        pSrc1++;
        pSrc2++;
      }
      result = (temp[0] + temp[1] + temp[2] + temp[3]);
    }
    clock_t stop = clock();

    runTimeSSE = static_cast<double>(stop - start)/CLOCKS_PER_SEC;
    cerr << "Runtime SSE: " << setprecision(20) << runTimeSSE << " seconds\n";
    cerr << result/N << endl;
  }
  {
    clock_t start = clock();

    float result;
    for( int j = 1; j < M; j++) {
      result = 0.;
      float* pSrc1 = array1;
      float* pSrc2 = array2;

      for ( int i = 0; i < N; i++ ) {
        // This code is easy to understand:
        result += (*pSrc1)*(*pSrc2);          // *pDest = *pSrc1 * *pSrc2

        pSrc1++;
        pSrc2++;
      }
    }
    clock_t stop = clock();
    runTimeNoSSE = static_cast<double>(stop - start)/CLOCKS_PER_SEC;
    cerr << "Runtime non-SSE: " << setprecision(20) << runTimeNoSSE << " seconds\n";
    cerr << result/N << endl;
  }
  cerr << "SSE speed increase factor: " << (runTimeNoSSE/runTimeSSE) << "x\n";
#ifdef __MINGW32_VERSION
    system("PAUSE");
#endif
  return 1;
}
