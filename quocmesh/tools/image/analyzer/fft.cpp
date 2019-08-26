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
 *
 *  \brief Computes the fourier transform of an input image and writes the output to an image, using qc::FourierTransform()
 *
 *  Usage: \code ./fft inputfile outputfile scale \endcode
 *
 *  Input and output image are interpreted as complex functions. The input image is a grey valued image, grey values are interpreted as real numbers between 0 and 1.
 *  In the fourier transform, the absolute value \code scale \endcode is interpreted as maximum saturation, 0 as no saturation, the color hue encodes the complex angle / phase.
 *  The origin is shifted to the center in the output image.
 *
 *  \author Lenz
 */

#include "convolution.h"

aol::Vec3<int> complexToColor ( std::complex<double> val, double scale ) {
  double radius = min ( abs ( val ) / scale, 1.0 );
  double angle = arg ( val ) / aol::NumberTrait<double>::pi; if ( angle < 0 ) angle += 1; // Now in (0;1]
  double hsv [3], rgb [3] = { 0, 0, 0 };
  aol::Vec3<int> res;
  hsv [0] = angle;
  hsv [1] = radius;
  hsv [2] = 1;
  aol::RGBColorMap<double>::hsv2rgb ( hsv, rgb );
  for ( int i = 0; i < 3; ++i ) res [i] = static_cast<int> ( 255.0 * rgb [i] );
  return res;
}

int main ( int argc, char* argv [] ) {
  if ( argc != 4 ) {
    cerr << "Usage: " << argv [0] << " input output scale" << endl;
    cerr << "Computes Fourier transform of input image and writes to output image in PPM format (color encoding phase)." << endl
    << "The transform is scaled so that the absolute value <scale> is full saturation." << endl;
    return -1;
  }

  // Read data
  double scale = aol::convert<double> ( argv [3] );
  qc::ScalarArray<double, qc::QC_2D> data ( argv [1] );
  int nx = data.getNumX (), ny = data.getNumY ();

  // Transform
  qc::MultiArray<double, 2, 2> function ( nx, ny ), transform ( nx, ny );
  function [0] = data;
  qc::FourierTransform ( function, transform );

  // Compute output, shift zero to center
  qc::MultiArray<int, 2, 3> image ( nx, ny );
  for ( int i = 0; i < nx; ++i )
    for ( int j = 0; j < ny; ++j )
      image.set ( aol::Vec2<short> ( ( i + nx / 2 ) % nx, ( j + ny / 2 ) % ny ), complexToColor ( std::complex<double> ( transform [0].get ( i, j ), transform [1].get ( i, j ) ), scale ) );
  ofstream file ( argv [2] );
  image.save ( file );
  return 0;
}
