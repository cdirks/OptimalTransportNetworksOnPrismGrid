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

#ifndef __FREQUENCYFILTERS_H
#define __FREQUENCYFILTERS_H
 
#include <multiArray.h>

template <typename RealType>
void complexMultiply ( const RealType ARe, const RealType AIm, const RealType BRe, const RealType BIm, RealType &ResRe, RealType &ResIm ) {
  const RealType a = ARe;
  const RealType b = AIm;
  const RealType c = BRe;
  const RealType d = BIm;
  ResRe = a * c - b * d;
  ResIm = b * c + a * d;
}

template <typename RealType>
void complexMultiply ( const qc::MultiArray<RealType, 2, 2> &ImageA, const qc::MultiArray<RealType, 2, 2> &ImageB, qc::MultiArray<RealType, 2, 2> &Result ) {
  const int numPixels = ImageA.getEqualComponentSize();
  for ( int i = 0; i < numPixels; i++ )
    complexMultiply ( ImageA[0][i], ImageA[1][i], ImageB[0][i], ImageB[1][i], Result[0][i], Result[1][i] );
}

template <typename RealType>
void saveComplexKernel ( const qc::MultiArray<RealType, 2, 2> &Kernel, const char * FileName ) {
  qc::ScalarArray<RealType, qc::QC_2D> tmp ( Kernel[0], aol::STRUCT_COPY_UNINIT );
  qc::ScalarArray<RealType, qc::QC_2D> tmp2 ( Kernel[0], aol::STRUCT_COPY_UNINIT );

  Kernel.getPointWiseNorm ( tmp2 );
  tmp.shiftByOffset ( tmp.getNumX() / 2, tmp.getNumY() / 2, tmp2 );
  tmp.scaleValuesTo01();
  for ( int i = 0; i < tmp.size(); ++i )
    tmp[i] = log ( tmp[i] + 1 ) / log ( static_cast<RealType>(2) );
  tmp.setOverflowHandlingToCurrentValueRange();
  tmp.savePNG ( FileName );
}

template <typename RealType>
void saveUndersampledImage ( const qc::ScalarArray<RealType, qc::QC_2D> &Image, const int UndersampleFactor, const char *BaseFileName ) {
  qc::ScalarArray<RealType, qc::QC_2D> undersampledImage ( Image.getNumX() / UndersampleFactor, Image.getNumY() / UndersampleFactor );
  for ( int j = 0; j < undersampledImage.getNumY(); ++j )
    for ( int i = 0; i < undersampledImage.getNumX(); ++i )
      undersampledImage.set ( i, j, Image.get ( UndersampleFactor*i+UndersampleFactor/2, UndersampleFactor*j+UndersampleFactor/2 ) );
  undersampledImage.savePNG ( aol::strprintf( "%s_undersampled.png", BaseFileName ).c_str() );
}

#endif //  __FREQUENCYFILTERS_H
