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

#ifndef __CONVOLUTION_H
#define __CONVOLUTION_H

#ifdef USE_LIB_FFTW
#include <fftw3.h>
#endif
#include <scalarArray.h>
#include <ChanVese.h>
#include <multiArray.h>

namespace qc {

#ifdef USE_LIB_FFTW

template <Dimension = QC_2D, typename RealType = double>
class ConvolutionTrait {};

template <>
class ConvolutionTrait<QC_2D,double> {
public:
  typedef fftw_complex FFTWComplex;
  typedef fftw_plan FFTWPlan;

  static void* fftwMalloc( size_t n ){
    return fftw_malloc( n );
  }
  static void fftwFree( void* p ){
    fftw_free( p );
  }
  static FFTWPlan fftwPlan_dft_r2c( aol::Vec2<int> nxy, double *in, FFTWComplex *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftw_plan_dft_r2c_2d( nxy[1], nxy[0], in, out, flags );
  }
  static FFTWPlan fftwPlan_dft_c2r( aol::Vec2<int> nxy, FFTWComplex *in, double *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftw_plan_dft_c2r_2d( nxy[1], nxy[0], in, out, flags );
  }
  static FFTWPlan fftwPlan_dft( aol::Vec2<int> nxy, FFTWComplex *in, FFTWComplex *out, int sign, unsigned flags ){
    return fftw_plan_dft_2d( nxy[1], nxy[0], in, out, sign, flags );
  }
  static void fftwDestroy_plan( FFTWPlan p ){
    fftw_destroy_plan( p );
  }
  static void fftwExecute( FFTWPlan p ){
    fftw_execute( p );
  }
};

template <>
class ConvolutionTrait<QC_3D,double> {
public:
  typedef fftw_complex FFTWComplex;
  typedef fftw_plan FFTWPlan;

  static void* fftwMalloc( size_t n ){
    return fftw_malloc( n );
  }
  static void fftwFree( void* p ){
    fftw_free( p );
  }
  static FFTWPlan fftwPlan_dft_r2c( aol::Vec3<int> nxyz, double *in, FFTWComplex *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftw_plan_dft_r2c_3d( nxyz[2], nxyz[1], nxyz[0], in, out, flags );
  }
  static FFTWPlan fftwPlan_dft_c2r( aol::Vec3<int> nxyz, FFTWComplex *in, double *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftw_plan_dft_c2r_3d( nxyz[2], nxyz[1], nxyz[0], in, out, flags );
  }
  static void fftwDestroy_plan( FFTWPlan p ){
    fftw_destroy_plan( p );
  }
  static void fftwExecute( FFTWPlan p ){
    fftw_execute( p );
  }
};

template <>
class ConvolutionTrait<QC_2D,float> {
public:
  typedef fftwf_complex FFTWComplex;
  typedef fftwf_plan FFTWPlan;

  static void* fftwMalloc( size_t n ){
    return fftwf_malloc( n );
  }
  static void fftwFree( void* p ){
    fftwf_free( p );
  }
  static FFTWPlan fftwPlan_dft_r2c( aol::Vec2<int> nxy, float *in, FFTWComplex *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftwf_plan_dft_r2c_2d( nxy[1], nxy[0], in, out, flags );
  }
  static FFTWPlan fftwPlan_dft_c2r( aol::Vec2<int> nxy, FFTWComplex *in, float *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftwf_plan_dft_c2r_2d( nxy[1], nxy[0], in, out, flags );
  }
  static FFTWPlan fftwPlan_dft( aol::Vec2<int> nxy, FFTWComplex *in, FFTWComplex *out, int sign, unsigned flags ){
    return fftwf_plan_dft_2d( nxy[1], nxy[0], in, out, sign, flags );
  }
  static void fftwDestroy_plan( FFTWPlan p ){
    fftwf_destroy_plan( p );
  }
  static void fftwExecute( FFTWPlan p ){
    fftwf_execute( p );
  }
};

template <>
class ConvolutionTrait<QC_3D,float> {
public:
  typedef fftwf_complex FFTWComplex;
  typedef fftwf_plan FFTWPlan;

  static void* fftwMalloc( size_t n ){
    return fftwf_malloc( n );
  }
  static void fftwFree( void* p ){
    fftwf_free( p );
  }
  static FFTWPlan fftwPlan_dft_r2c( aol::Vec3<int> nxyz, float *in, FFTWComplex *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftwf_plan_dft_r2c_3d( nxyz[2], nxyz[1], nxyz[0], in, out, flags );
  }
  static FFTWPlan fftwPlan_dft_c2r( aol::Vec3<int> nxyz, FFTWComplex *in, float *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftwf_plan_dft_c2r_3d( nxyz[2], nxyz[1], nxyz[0], in, out, flags );
  }
  static void fftwDestroy_plan( FFTWPlan p ){
    fftwf_destroy_plan( p );
  }
  static void fftwExecute( FFTWPlan p ){
    fftwf_execute( p );
  }
};

template <>
class ConvolutionTrait<QC_2D,long double> {
public:
  typedef fftwl_complex FFTWComplex;
  typedef fftwl_plan FFTWPlan;

  static void* fftwMalloc( size_t n ){
    return fftwl_malloc( n );
  }
  static void fftwFree( void* p ){
    fftwl_free( p );
  }
  static FFTWPlan fftwPlan_dft_r2c( aol::Vec2<int> nxy, long double *in, FFTWComplex *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftwl_plan_dft_r2c_2d( nxy[1], nxy[0], in, out, flags );
  }
  static FFTWPlan fftwPlan_dft_c2r( aol::Vec2<int> nxy, FFTWComplex *in, long double *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftwl_plan_dft_c2r_2d( nxy[1], nxy[0], in, out, flags );
  }
  static FFTWPlan fftwPlan_dft( aol::Vec2<int> nxy, FFTWComplex *in, FFTWComplex *out, int sign, unsigned flags ){
    return fftwl_plan_dft_2d( nxy[1], nxy[0], in, out, sign, flags );
  }
  static void fftwDestroy_plan( FFTWPlan p ){
    fftwl_destroy_plan( p );
  }
  static void fftwExecute( FFTWPlan p ){
    fftwl_execute( p );
  }
};

template <>
class ConvolutionTrait<QC_3D,long double> {
public:
  typedef fftwl_complex FFTWComplex;
  typedef fftwl_plan FFTWPlan;

  static void* fftwMalloc( size_t n ){
    return fftwl_malloc( n );
  }
  static void fftwFree( void* p ){
    fftwl_free( p );
  }
  static FFTWPlan fftwPlan_dft_r2c( aol::Vec3<int> nxyz, long double *in, FFTWComplex *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftwl_plan_dft_r2c_3d( nxyz[2], nxyz[1], nxyz[0], in, out, flags );
  }
  static FFTWPlan fftwPlan_dft_c2r( aol::Vec3<int> nxyz, FFTWComplex *in, long double *out, unsigned flags ){
    // The order of the dimensions in each direction is reversed so to pretend that our ScalarArrays are rotated,
    // but in row-major order (in reality, they are in column-major order).
    return fftwl_plan_dft_c2r_3d( nxyz[2], nxyz[1], nxyz[0], in, out, flags );
  }
  static void fftwDestroy_plan( FFTWPlan p ){
    fftwl_destroy_plan( p );
  }
  static void fftwExecute( FFTWPlan p ){
    fftwl_execute( p );
  }
};

#endif // USE_LIB_FFTW


#ifdef USE_LIB_FFTW

enum FourierTransformDirection { FTForward = FFTW_FORWARD, FTBackward = FFTW_BACKWARD };

template <typename RealType>
void FourierTransform ( const qc::MultiArray<RealType, 2, 2>& function, qc::MultiArray<RealType, 2, 2>& transform, enum FourierTransformDirection direction = FTForward );

#else
enum FourierTransformDirection { FTForward, FTBackward };

template <typename RealType>
void FourierTransform ( const qc::MultiArray<RealType, 2, 2>& /*function*/, qc::MultiArray<RealType, 2, 2>& /*transform*/, enum FourierTransformDirection /*direction*/ = FTForward );

#endif // USE_LIB_FFTW

/**
 * Transforms an image to the Fourier domain, zeroes all complex Fourier coefficients with a norm
 * smaller then the given threshold and then transforms the thresholded coefficients back to the
 * original domain (discarding the imaginary part).
 *
 * \author Berkels
 */
template <typename RealType>
void fourierThreshold ( qc::ScalarArray<RealType, qc::QC_2D> &Image, const RealType AbsoluteThreshold ) {
  const int nx = Image.getNumX (), ny = Image.getNumY ();
  qc::MultiArray<RealType, 2, 2> function ( nx, ny ), transform ( nx, ny );
  function [0] = Image;
  qc::FourierTransform ( function, transform, qc::FTForward );

  for ( qc::RectangularIterator<qc::QC_2D> it ( Image ); it.notAtEnd(); ++it ) {
    if ( transform.get ( *it ).norm() < AbsoluteThreshold )
      transform.set( *it, aol::Vec2<RealType> ( 0, 0 ) );
  }

  qc::FourierTransform ( transform, function, qc::FTBackward );
  cerr << function[0].getMaxValue();
  function[0] /= ( nx * ny );
  Image = function[0];
}

/**
 * Computes log2 ( 1 + FFT modulus ) pointwise, possible dropping the lowest DropPercentage values.
 *
 * \author Berkels
 */
template <typename RealType>
void computeLogFFTModulus ( const qc::ScalarArray<RealType, qc::QC_2D> &Image, qc::ScalarArray<RealType, qc::QC_2D> &Modulus, const RealType DropPercentage = 0, const bool ApplyLog = true ) {
  const int nx = Image.getNumX (), ny = Image.getNumY ();

  // Transform
  qc::MultiArray<RealType, 2, 2> function ( nx, ny ), transform ( nx, ny );
  function [0] = Image;
  qc::FourierTransform ( function, transform );

  // Compute output, shift zero to center
  for ( int i = 0; i < nx; ++i ) {
    for ( int j = 0; j < ny; ++j ) {
      const aol::Vec2<short> pos ( ( i + nx / 2 ) % nx, ( j + ny / 2 ) % ny );
      Modulus.set ( pos , aol::Vec2<RealType> ( transform [0].get ( i, j ), transform [1].get ( i, j ) ).norm() );
    }
  }

  if ( ApplyLog ) {
    for ( int i = 0; i < Modulus.size(); ++i )
      Modulus[i] = log ( Modulus[i] + 1 ) / log ( static_cast<RealType> ( 2 ) );
  }

  if ( DropPercentage > 0 ) {
    aol::Vector<RealType> temp ( Modulus );
    temp.sortValues();
    const int sizeM1 = temp.size()-1;
    Modulus.clamp ( temp [ aol::Clamp ( static_cast<int> ( DropPercentage * sizeM1 ), 0, sizeM1 ) ], temp[ sizeM1 ]);
  }
}

/**
 * Two and three dimensional circular convolution using discrete Fourier transform (DFT).
 * Needs libfftw (for more details on that library, visit the fftw website www.fftw.org).
 * With zero padding linear convolution can be emulated.
 *
 * \author Berkels, Wirth
 */
#ifdef USE_LIB_FFTW
template<Dimension Dim= QC_2D, typename RealType = double>
class Convolution {

private:
  typedef ConvolutionTrait<Dim,RealType> ConvType;
  typedef typename ConvType::FFTWComplex FFTWComplex;
  typedef typename ConvType::FFTWPlan FFTWPlan;

  // product of sizes of rectangular grid in x, y, z direction, and same for the FFT result (here, the last dimension is always divided by 2 plus 1)
  const int _numPixels;
  const int _numComplexPixels;
  // arrays "_f" for the (real) input of the FFT and "_g1,_g2" for the (complex) output of the FFT (one for the image and one for the kernel)
  RealType* _f;
  FFTWComplex* _g1;
  FFTWComplex* _g2;
  // FFT procedures for the Fourier transformation and its inverse (will be specified in the constructor)
  FFTWPlan _forwardPlan;
  FFTWPlan _backwardPlan;

public:
  Convolution ( const typename aol::VecDimTrait<int,Dim>::VecType NumXYZ ) :
    _numPixels ( NumXYZ.prod() ),
    _numComplexPixels ( _numPixels / NumXYZ[0] * ( NumXYZ[0] / 2 + 1 ) ),
    // allocate sufficient memory for the FFT input and output
    _f ( static_cast<RealType*> ( ConvType::fftwMalloc ( sizeof ( RealType ) * _numPixels ) ) ),
    _g1 ( static_cast<FFTWComplex*> ( ConvType::fftwMalloc ( sizeof ( FFTWComplex ) * _numComplexPixels ) ) ),
    _g2 ( static_cast<FFTWComplex*> ( ConvType::fftwMalloc ( sizeof ( FFTWComplex ) * _numComplexPixels ) ) ) {

    // Define the fftw plan "_forwardPlan" to do a discrete Fourier transform (dft) from real to complex (r2c) in 3d
    // on a _numX by _numY by _numZ grid, taking the array "_f" as (real) input and "_g1" as (complex) output.
    // Define the fftw plan "_backwardPlan" to do an inverse discrete Fourier transform (dft) from real to complex (r2c) in 3d
    // on a _numX by _numY by _numZ grid, taking the array "_g1" as (real) input and "_f" as (complex) output.
    // The flag "FFTW_MEASURE" as opposed to "FFTW_ESTIMATE" makes the routine to try some test calculations to sort out the fastest plan.
    cerr << "preparing Fourier transformations... ";
    _forwardPlan  = ConvType::fftwPlan_dft_r2c ( NumXYZ, _f, _g1, FFTW_MEASURE );
    _backwardPlan = ConvType::fftwPlan_dft_c2r ( NumXYZ, _g1, _f, FFTW_MEASURE );
    cerr << "done.\n";
  }

  ~Convolution() {
    ConvType::fftwDestroy_plan ( _forwardPlan );
    ConvType::fftwDestroy_plan ( _backwardPlan );
    ConvType::fftwFree ( _f );
    ConvType::fftwFree ( _g1 );
    ConvType::fftwFree ( _g2 );
  }

private:
  void inputDataFFT ( const qc::ScalarArray<RealType, Dim> &InputA, const qc::ScalarArray<RealType, Dim> &InputB ) const {
    // compute the Fourier transform of InputA (result will be in "_g1")
    InputA.copyToBuffer ( _f );
    ConvType::fftwExecute ( _forwardPlan );

    // copy the result of the FFT on the ImageA in "_g1" to "_g2"
    memcpy ( _g2, _g1, sizeof ( FFTWComplex ) * _numComplexPixels );

    // compute the Fourier transform of InputB (result will be in "_g1")
    InputB.copyToBuffer ( _f );
    ConvType::fftwExecute ( _forwardPlan );
  }

  void outputDataIFFT ( qc::ScalarArray<RealType, Dim> &Dest ) const {
    ConvType::fftwExecute ( _backwardPlan );
    Dest.readFromBuffer ( _f );

    // normalize the result (the FFT and inverse FFT were only inverse to each other up to a constant factor)
    Dest *= 1.0 / _numPixels;
  }

public:
  void convolve ( const qc::ScalarArray<RealType, Dim> &Image, const qc::ScalarArray<RealType, Dim> &Kernel, qc::ScalarArray<RealType, Dim> &Dest ) const {
    inputDataFFT ( Image, Kernel );

    // compute the pointwise product of the Fourier transforms of image and kernel
    for ( int i = 0; i < _numComplexPixels; i++ ) {
      // Fourier transform of kernel has value a+ib at that point
      const double a = _g1[i][0];
      const double b = _g1[i][1];
      // Fourier transform of image has value c+id at that point
      const double c = _g2[i][0];
      const double d = _g2[i][1];
      // pointwise product is ac-bd+i(bc+ad) at that point
      _g1[i][0] = a * c - b * d;
      _g1[i][1] = b * c + a * d;
    }

    // compute the inverse FFT of the pointwise product (result will be in "_f" and then copied to Dest), which is the convolution
    outputDataIFFT ( Dest );
  }

  void phaseCorrelation ( const qc::ScalarArray<RealType, Dim> &ImageA, const qc::ScalarArray<RealType, Dim> &ImageB, qc::ScalarArray<RealType, Dim> &Dest ) const {
    inputDataFFT ( ImageA, ImageB );

    for ( int i = 0; i < _numComplexPixels; i++ ) {
      // Calculate the complex product of referenceFFT (_g2) and the complex conjugate of templateFFT (_g1).
      const aol::Vec2<RealType> tmp ( _g2[i][0]*_g1[i][0] + _g2[i][1]*_g1[i][1], _g2[i][1]*_g1[i][0] - _g2[i][0]*_g1[i][1] );
      const RealType norm = tmp.norm();
      const RealType scale = aol::appeqAbsolute ( norm, aol::ZTrait<RealType>::zero ) ? 0 : ( 1 / norm );
      for ( int j = 0; j < 2; ++j )
        _g1[i][j] = tmp[j] * scale;
    }

    // compute the inverse FFT of the pointwise product (result will be in "_f" and then copied to Dest), which is the convolution
    outputDataIFFT ( Dest );
  }

};

#else
template<Dimension Dim = QC_2D, typename RealType = double>
class Convolution {
public:
  Convolution ( const typename aol::VecDimTrait<int,Dim>::VecType /*NumXYZ*/ ) {
    throw aol::Exception ( "Convolution needs libfftw! Compile with -DUSE_LIB_FFTW", __FILE__, __LINE__ );
  }
  ~Convolution() {}
  void convolve ( const qc::ScalarArray<RealType, Dim> &/*Image*/, const qc::ScalarArray<RealType, Dim> &/*Kernel*/, qc::ScalarArray<RealType, Dim> &/*Dest*/ ) const {}
  void phaseCorrelation ( const qc::ScalarArray<RealType, Dim> &/*ImageA*/, const qc::ScalarArray<RealType, Dim> &/*ImageB*/, qc::ScalarArray<RealType, Dim> &/*Dest*/ ) const {}
};

#endif // USE_LIB_FFTW

/**
 * Two dimensional linear convolution using discrete Fourier transform (DFT), needs libfftw.
 * Uses zero padding to emulate linear convolution from circular convolution.
 *
 * ATTENTION: This appromimates linear convolution, i.e. integrating over whole R^2. You
 * will notice that the functions are extended to R^2 by zero at the boundary.
 *
 * \author Berkels
 */
template<typename RealType = double>
class LinearConvolution {
  const int _numX, _numY;
  qc::Convolution<QC_2D,RealType> _conv;
public:
  LinearConvolution ( const int NumX, const int NumY )
    : _numX ( NumX ),
      _numY ( NumY ),
      _conv ( aol::Vec2<int>( 2*_numX, 2*_numY ) ){
  }
  void convolve ( const qc::ScalarArray<RealType, qc::QC_2D> &Image, const qc::ScalarArray<RealType, qc::QC_2D> &Kernel, qc::ScalarArray<RealType, qc::QC_2D> &Dest ) const {

    qc::ScalarArray<RealType, qc::QC_2D> extendedImage ( 2*_numX, 2*_numY );
    qc::ScalarArray<RealType, qc::QC_2D> extendedKernel ( 2*_numX, 2*_numY );

    const int xOffset = _numX / 2 + 1;
    const int yOffset = _numY / 2 + 1;
    // First extend the image and the kernel with zero padding.
    // Note: The kernel needs to be treated carefully.
    for ( int i = 0; i < _numX; ++i ){
      for ( int j = 0; j < _numY; ++j ){
        extendedImage.set( i + xOffset, j + yOffset, Image.get( i, j ) );
        extendedKernel.set( i + (i/xOffset)*_numX, j + (j/yOffset)*_numY, Kernel.get( i, j ) );
      }
    }

    // Do the convolution of the extended image with the extended kernel.
    _conv.convolve ( extendedImage, extendedKernel, extendedImage );

    // Copy the values from the extendedImage to the smaller destination image.
    for ( int i = 0; i < _numX; ++i ){
      for ( int j = 0; j < _numY; ++j ){
        Dest.set( i, j, extendedImage.get( i + xOffset, j + yOffset) );
      }
    }
  }

  //! Convoles an image of size MxN with a kernel of size 2Mx2N by extending the image with mirroring to size 2Mx2M.
  void convolveMirrored ( const qc::ScalarArray<RealType, qc::QC_2D> &Image, const qc::ScalarArray<RealType, qc::QC_2D> &ExtendedKernel, qc::ScalarArray<RealType, qc::QC_2D> &Dest ) const {
    qc::ScalarArray<RealType, qc::QC_2D> extendedImage ( 2*_numX, 2*_numY );
    for ( int i = 0; i < _numX; ++i ){
      for ( int j = 0; j < _numY; ++j ){
        const RealType value = Image.get( i, j );
        extendedImage.set( i, j, value );
        extendedImage.set( 2 * _numX - 1 - i, j, value );
        extendedImage.set( i, 2 * _numY - 1 - j, value );
        extendedImage.set( 2 * _numX - 1 - i, 2 * _numY - 1 - j, value );
      }
    }

    // Do the convolution of the extended image with the extended kernel.
    _conv.convolve ( extendedImage, ExtendedKernel, extendedImage );

    // Copy the values from the extendedImage to the smaller destination image.
    for ( int i = 0; i < _numX; ++i ){
      for ( int j = 0; j < _numY; ++j ){
        Dest.set( i, j, extendedImage.get( i, j) );
      }
    }
  }
};

template <typename RealType>
void generateMotionBlurKernel ( const aol::Vec2<RealType> &Velocity, qc::ScalarArray<RealType, qc::QC_2D> &Kernel ) {
  Kernel.setZero();
  aol::PolynomialHeavisideFunction<RealType> heavisideFunction ( 1. );
  aol::PolynomialHeavisideFunction<RealType> heavisideFunctionLength ( 1. );
  const RealType vNorm = Velocity.norm();

  const int numX = Kernel.getNumX();
  const int numY = Kernel.getNumY();
  const int xOffset = numX / 2 + 1;
  const int yOffset = numY / 2 + 1;

  for ( int i = 0; i < numX; i++ ) {
    for ( int j = 0; j < numY; j++ ) {
      const int x = i - numX / 2;
      const int y = j - numY / 2;
      RealType temp = heavisideFunctionLength.evaluate ( vNorm - sqrt ( aol::Sqr ( static_cast<RealType> ( x ) ) + aol::Sqr ( static_cast<RealType> ( y ) ) ) );
      Kernel.set ( ( i + xOffset ) % numX, ( j + yOffset ) % numY , heavisideFunction.evaluateDerivative ( Velocity[0]*static_cast<RealType> ( y ) - Velocity[1]*static_cast<RealType> ( x ) ) * temp );
    }
  }
  // Ideally dividing by vNorm should be enough to normalize the kernel
  Kernel /= Kernel.sum();
}

//! generates the kernel exp(1-sigma^p/(sigma^p-x^p)) (optionally normalized by its l_1-norm)
template <typename RealType>
void generateCompactSupportKernel1D ( const RealType Sigma, const RealType P, qc::ScalarArray<RealType, qc::QC_1D> &Kernel, const bool Normalize = true) {
  Kernel.setZero();
  const int numX = Kernel.getNumX();
  const RealType sigmaToP = pow( Sigma, P );
  const RealType h = aol::NumberTrait<RealType>::one /( numX - 1 );

  for ( int i = 0; i < numX; i++ ) {
    const RealType x = (i - numX / 2) *h;
    if ( aol::Abs( x ) < Sigma )
      Kernel.set (  i , exp( 1. - sigmaToP / ( sigmaToP - pow( x, P ) ) ) );
  }
  if ( Normalize )
    Kernel /= Kernel.sum();
}

//! generates the kernel -sigma^p*(p*x^(p-1))/(sigma^p-x^p)^2*exp(1-sigma^p/(sigma^p-x^p)) (optionally normalized by its l_1-norm)
template <typename RealType>
void generateCompactSupportKernel1DGrad ( const RealType Sigma, const RealType P, qc::ScalarArray<RealType, qc::QC_1D> &Kernel, const bool Normalize = true) {
  Kernel.setZero();
  const int numX = Kernel.getNumX();
  const RealType sigmaToP = pow( Sigma, P );
  const RealType h = aol::NumberTrait<RealType>::one /( numX - 1 );

  for ( int i = 0; i < numX; i++ ) {
    const RealType x = (i - numX / 2) *h;
    const RealType temp = aol::NumberTrait<RealType>::one / ( sigmaToP - pow( x, P ) );
    if ( aol::Abs( x ) < Sigma )
      Kernel.set (  i , - sigmaToP * aol::Sqr( temp ) * P * pow( x, P - 1. ) * exp( 1. - sigmaToP * temp ) );
  }
  if ( Normalize ){
    qc::ScalarArray<RealType,qc::QC_1D> kernel( Kernel, aol::STRUCT_COPY );
    generateCompactSupportKernel1D ( Sigma, P, kernel, false );
    Kernel /= kernel.sum();
  }
}

template <typename RealType>
void generateGaussKernel1D ( const RealType Sigma, qc::ScalarArray<RealType, qc::QC_1D> &Kernel, const bool Normalize = true) {
  Kernel.setZero();
  const int numX = Kernel.getNumX();
  const RealType twoSigmaSqr = 2*aol::Sqr( Sigma );
  const RealType h = aol::NumberTrait<RealType>::one /( numX - 1 );

  for ( int i = 0; i < numX; i++ ) {
      const RealType x = (i - numX / 2) *h;
      RealType temp = exp( - ( aol::Sqr ( x ) ) / twoSigmaSqr );
      Kernel.set (  i , temp );
    }
  if ( Normalize )
    Kernel /= Kernel.sum();
}

template <typename RealType>
void generateGaussKernel1DGrad ( const RealType Sigma, qc::ScalarArray<RealType, qc::QC_1D> &Kernel, const bool Normalize = true) {
  Kernel.setZero();
  const int numX = Kernel.getNumX();
  const RealType twoSigmaSqr = 2*aol::Sqr( Sigma );
  const RealType h = aol::NumberTrait<RealType>::one /( numX - 1 );

  for ( int i = 0; i < numX; i++ ) {
      const RealType x = (i - numX / 2) *h;
      RealType temp = exp( - ( aol::Sqr ( x ) ) / twoSigmaSqr );
      Kernel.set (  i , - temp*x / aol::Sqr(Sigma) );
    }
  if ( Normalize ){
    qc::ScalarArray<RealType,qc::QC_1D> kernel( Kernel, aol::STRUCT_COPY );
    generateGaussKernel1D ( Sigma, kernel, false);
    Kernel /= kernel.sum();
  }
}

template <typename RealType>
void generateGaussKernel ( const RealType Sigma, qc::ScalarArray<RealType, qc::QC_2D> &Kernel, const bool Normalize = true, const bool CenterAt00 = true ) {
  Kernel.setZero();

  const int numX = Kernel.getNumX();
  const int numY = Kernel.getNumY();
  const RealType twoSigmaSqr = 2*aol::Sqr( Sigma );
  const RealType h = aol::NumberTrait<RealType>::one /( aol::Max ( numX - 1, numY - 1 ) );
  const int xOffset = CenterAt00 ? ( numX / 2 + numX % 2 ) : 0;
  const int yOffset = CenterAt00 ? ( numY / 2 + numY % 2 ) : 0;

  for ( int i = 0; i < numX; i++ ) {
    for ( int j = 0; j < numY; j++ ) {
      const RealType x = (i - numX / 2) *h;
      const RealType y = (j - numY / 2) *h;
      RealType temp = exp( - ( aol::Sqr ( x ) + aol::Sqr ( y ) ) / twoSigmaSqr );
      Kernel.set ( ( i + xOffset ) % numX, ( j + yOffset ) % numY , temp );
    }
  }
  if ( Normalize )
    Kernel /= Kernel.sum();
}

template <typename RealType>
void generateGaussKernelGrad ( const RealType Sigma, aol::auto_container<qc::QC_2D,qc::ScalarArray<RealType, qc::QC_2D> > &KernelGrad, const bool Normalize = true, const bool CenterAt00 = true ) {
  for ( int i = 0; i < qc::QC_2D; i++ )
    KernelGrad[i].setZero();

  const int numX = KernelGrad[0].getNumX();
  const int numY = KernelGrad[0].getNumY();
  const RealType twoSigmaSqr = 2*aol::Sqr( Sigma );
  const RealType h = aol::NumberTrait<RealType>::one /( aol::Max ( numX - 1, numY - 1 ) );
  const int xOffset = CenterAt00 ? ( numX / 2 + 1 ) : 0;
  const int yOffset = CenterAt00 ? ( numY / 2 + 1 ) : 0;

  for ( int i = 0; i < numX; i++ ) {
    for ( int j = 0; j < numY; j++ ) {
      const RealType x = (i - numX / 2) *h;
      const RealType y = (j - numY / 2) *h;
      RealType temp = exp( - ( aol::Sqr ( x ) + aol::Sqr ( y ) ) / twoSigmaSqr );
      KernelGrad[0].set ( ( i + xOffset ) % numX, ( j + yOffset ) % numY , - temp * x / aol::Sqr( Sigma ) );
      KernelGrad[1].set ( ( i + xOffset ) % numX, ( j + yOffset ) % numY , - temp * y / aol::Sqr( Sigma ) );
    }
  }
  if ( Normalize ) {
    qc::ScalarArray<RealType,qc::QC_2D> kernel( KernelGrad[0], aol::STRUCT_COPY );
    generateGaussKernel ( Sigma, kernel, false, CenterAt00 );
    for ( int i = 0; i < qc::QC_2D; i++ )
      KernelGrad[i] /= kernel.sum();
  }
}

void addMotionBlurToArray ( const aol::Vec2<double> &Velocity, const qc::ScalarArray<double, qc::QC_2D> &Arg, qc::ScalarArray<double, qc::QC_2D> &Dest );

template <typename RealType>
class SigmaTauOpBase {
private:
  RealType _sigma;
  virtual void sigmaWasChanged ( ) = 0;

public:
  SigmaTauOpBase () : _sigma( aol::NumberTrait<RealType>::zero ) {}

  virtual ~SigmaTauOpBase () {}

  void setSigma ( const RealType Sigma ) {
    if ( _sigma != Sigma ){
      _sigma = Sigma;
      sigmaWasChanged();
    }
  }

  RealType getSigma (  ) const {
    return _sigma;
  }

  void setTau ( const RealType Tau ) {
    setSigma ( sqrt ( 2*Tau ) );
  }
};

/**
 * Smoothens the argument by linear convolution with a Gauss-kernel.
 * _sigma specifies the width of the kernel.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
class LinearConvolutionSmoothOp {
};

// qc::LinearConvolution:apply is not thread-safe, so the parallelization in aol::BiOp has to be turned off.
template <typename RealType>
class LinearConvolutionSmoothOp<RealType, qc::QC_2D> : public aol::BiOp< aol::Vector<RealType>, false>, public qc::SigmaTauOpBase<RealType> {
  const int _numX, _numY;
  qc::LinearConvolution<RealType> _conv;
  qc::ScalarArray<RealType, qc::QC_2D> _kernel;

  void sigmaWasChanged ( ) {
    qc::generateGaussKernel<RealType> ( this->getSigma(), _kernel );
  }
public:
  LinearConvolutionSmoothOp ( const int NumX, const int NumY )
    : _numX ( NumX ),
      _numY ( NumY ),
      _conv ( _numX, _numY ),
      _kernel ( _numX, _numY ) {
  }

  LinearConvolutionSmoothOp ( const qc::GridSize<qc::QC_2D> &Size )
    : _numX ( Size.getNumX() ),
      _numY ( Size.getNumY() ),
      _conv ( _numX, _numY ),
      _kernel ( _numX, _numY ) {
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    const qc::ScalarArray<RealType, qc::QC_2D> argArray( Arg, _numX, _numY, aol::FLAT_COPY);
    qc::ScalarArray<RealType, qc::QC_2D> destArray( Dest, _numX, _numY, aol::FLAT_COPY);
    _conv.convolve ( argArray, _kernel, destArray );
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  using aol::BiOp<aol::Vector<RealType>, false>::applyAdd;
  using aol::BiOp<aol::Vector<RealType>, false>::apply;
};

/**
 * Class to apply \f$ (M + \tau L)^{-1} \f$ approximately using convolution with a Gaussian.
 *
 * \warning Seems to smooth significantly more than qc::CGBasedInverseH1Metric, possibly sigma is not adjusted properly.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class ConvolutionBasedInverseH1Metric : public aol::BiOp< aol::Vector<typename ConfiguratorType::RealType> >, public qc::SigmaTauOpBase<typename ConfiguratorType::RealType> {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::LumpedMassOp<ConfiguratorType> _lumpedMassInv;
  qc::LinearConvolution<RealType> _conv;
  qc::ScalarArray<RealType, qc::QC_2D> _kernel;

  void sigmaWasChanged ( ) {
    // Since the domain of the argument in apply gets extended, we have to adjust sigma in the kernel.
    qc::generateGaussKernel<RealType> ( this->getSigma() / 2, _kernel );
  }
public:
  ConvolutionBasedInverseH1Metric( const typename ConfiguratorType::InitType &Initializer )
    : _grid ( Initializer ),
      _lumpedMassInv ( Initializer, aol::INVERT ),
      _conv ( _grid.getNumX(), _grid.getNumY() ),
      _kernel ( 2*_grid.getNumX(), 2*_grid.getNumY() ) {}

  ~ConvolutionBasedInverseH1Metric() {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest, aol::DEEP_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> mInvArg ( Arg, aol::STRUCT_COPY );
    // Convolution with a Gaussian solves the heat equation, which in FE terms also involves
    // multiplying the input with M. To get rid of this inherent multiplication multiply
    // the input with M^{-1} first.
    _lumpedMassInv.apply( Arg, mInvArg );
    qc::ScalarArray<RealType, qc::QC_2D> argArray( mInvArg, _grid, aol::FLAT_COPY);
    qc::ScalarArray<RealType, qc::QC_2D> destArray( Dest, _grid, aol::FLAT_COPY);

    _conv.convolveMirrored ( argArray, _kernel, destArray );
  }

  using aol::BiOp<aol::Vector<RealType> >::applyAdd;
  using aol::BiOp<aol::Vector<RealType> >::apply;
};

/**
 * Gradient of LinearConvolutionSmoothOp
 *
 * \author Wirth, Heeren, Boerdgen
 */
template <typename RealType, qc::Dimension Dim>
class LinearConvolutionSmoothGradientOp {
};

template <typename RealType>
class LinearConvolutionSmoothGradientOp<RealType, qc::QC_2D> : public aol::Op< aol::Vector<RealType>, aol::MultiVector<RealType> > {
  const int _numX, _numY;
  RealType _sigma;
  qc::LinearConvolution<RealType> _conv;
  aol::auto_container<qc::QC_2D,qc::ScalarArray<RealType, qc::QC_2D> > _kernelGrad;
public:
  LinearConvolutionSmoothGradientOp ( const int NumX, const int NumY )
    : _numX ( NumX ),
      _numY ( NumY ),
      _sigma( aol::NumberTrait<RealType>::zero ),
      _conv ( _numX, _numY ) {
    qc::ScalarArray<RealType, qc::QC_2D> aux ( _numX, _numY );
    for ( int i = 0; i < qc::QC_2D; i++ )
      _kernelGrad.set_copy( i, aux );
  }

  void setSigma ( const RealType Sigma ) {
    if ( _sigma != Sigma ){
      _sigma = Sigma;
      qc::generateGaussKernelGrad<RealType> ( Sigma, _kernelGrad );
    }
  }

  void setTau ( const RealType Tau ) {
    setSigma ( sqrt ( 2*Tau ) );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    const qc::ScalarArray<RealType, qc::QC_2D> argArray( Arg, _numX, _numY, aol::FLAT_COPY );
    for ( int i = 0; i < qc::QC_2D; i++ ) {
      qc::ScalarArray<RealType, qc::QC_2D> destArray( Dest[i], _numX, _numY, aol::FLAT_COPY );
      _conv.convolve ( argArray, _kernelGrad[i], destArray );
    }
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> tmp ( Dest, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

};

} // namespace qc

#endif // __CONVOLUTION_H
