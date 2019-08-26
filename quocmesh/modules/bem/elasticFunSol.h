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

#ifndef __ELASTICFUNSOL_H
#define __ELASTICFUNSOL_H

#include <bemesh.h>
#include <bemOps.h>
#include <elasticTensor.h>

namespace bm {

//! ElasticSolution Class
/** Single Layer Boundary Integral Operator for piecewise constant ansatz functions.
 *  Overwrites LocalOperator::evaluateLocally()
 */
template <class _ParticleType>
class ElasticSolution
      : public LocalOperator<ElasticSolution<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  ElasticSolution ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp )
      : LocalOperator<ElasticSolution<ParticleType>, ParticleType > (), _greenCoeffs ( greenCoeffs ), _evalComp ( evalComp ), _intComp ( intComp ) {}
  using LocalOperator<ElasticSolution<ParticleType>, ParticleType >::evaluateLocally;

  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool, bool, bool ) const {
    int alpha, j;
    std::complex<DataType> jsum = 0, c1, c2, c_start, c_end, log_start, log_end, result = 0;
    DataType t0, arg_start, arg_end;
    aol::Vec2<DataType> dir, dirs, dire;

    for ( alpha = 0; alpha < 2; alpha++ ) {
      for ( j = 0, jsum = 0; j < 2; j++ )
        jsum += _greenCoeffs.A [_intComp][alpha] * _greenCoeffs.N [alpha][j] * _greenCoeffs.d [j][_evalComp];

      dirs = integrate.getStart () - evaluate;
      dire = integrate.getEnd () - evaluate;
      dir = integrate.getDirection ();
      dir.normalize ();
      aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );

      c_start = pav * aol::Vec2< complex< DataType> >( dirs );
      c_end = pav * aol::Vec2< complex< DataType> >( dire );

      log_start = log ( c_start );
      log_end = log ( c_end );
      arg_start = imag ( log_start );
      arg_end = imag ( log_end );

      c1 = c_start;
      c2 = pav * aol::Vec2< complex< DataType> >( dir );
      result += ( ( c_end * log_end - c_start * log_start ) / c2 - integrate.getLength () ) * jsum;
      if ( ( arg_start >= 0 ) && ( arg_end < 0 ) ) { // Jump due to sign change of the imag part of the log argument
        t0 = - imag ( c1 ) / imag ( c2 );
        if ( real ( c1 + c2 * t0 ) < 0 )
          result += 2 * aol::I * aol::NumberTrait<long double>::pi * ( c1 + c2 * t0 ) / c2 * jsum;
      }
      if ( ( arg_start < 0 ) && ( arg_end >= 0 ) ) {
        t0 = - imag ( c1 ) / imag ( c2 );
        if ( real ( c1 + c2 * t0 ) < 0 )
          result -= 2 * aol::I * aol::NumberTrait<long double>::pi * ( c1 + c2 * t0 ) / c2 * jsum;
      }
    }

    if ( aol::debugging::finiteness && !aol::isFinite ( real ( result ) ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;

    return real ( result ) / ( 2 * aol::NumberTrait<long double>::pi );
  }
private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp;
};

//! ElasticSolutionDiff Class
/** Double Layer Boundary Integral Operator for piecewise constant ansatz functions.
 *  Overwrites LocalOperator::evaluateLocally()
 */
template <class _ParticleType>
class ElasticSolutionDiff
      : public LocalOperator<ElasticSolutionDiff<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  ElasticSolutionDiff ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp )
      : LocalOperator<ElasticSolutionDiff<ParticleType>, ParticleType > (), _greenCoeffs ( greenCoeffs ), _evalComp ( evalComp ), _intComp ( intComp ) {}
  using LocalOperator<ElasticSolutionDiff<ParticleType>, ParticleType >::evaluateLocally;
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool, bool isinterior, bool ) const {
    int alpha, j, t;
    complex<DataType> jtsum = 0, c1, c2, c_start, c_end, log_start, log_end, result = 0;
    DataType t0, arg_start, arg_end;
    aol::Vec2<DataType> dir, dirs, dire;

    if ( isinterior ) return 0;

    for ( alpha = 0; alpha < 2; alpha++ ) {
      for ( j = 0, jtsum = 0; j < 2; j++ )
        for ( t = 0; t < 2; t++ )
          jtsum += _greenCoeffs.L [_intComp][j][alpha] * _greenCoeffs.N [alpha][t] * integrate.getNormal () [j] * _greenCoeffs.d [t][_evalComp];

      dirs = integrate.getStart () - evaluate;
      dire = integrate.getEnd () - evaluate;
      dir = integrate.getDirection ();
      dir.normalize ();
      aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );

      c_start = pav * aol::Vec2< complex< DataType> >( dirs );
      c_end = pav * aol::Vec2< complex< DataType> >( dire );

      log_start = log ( c_start );
      log_end = log ( c_end );
      arg_start = imag ( log_start );
      arg_end = imag ( log_end );

      c1 = c_start;
      c2 = pav * aol::Vec2< complex< DataType> >( dir );
      result += ( log_end - log_start ) / c2 * jtsum;
      if ( ( arg_start >= 0 ) && ( arg_end < 0 ) ) { // Jump due to sign change of the imag part of the log argument
        t0 = - imag ( c1 ) / imag ( c2 );
        if ( real ( c1 + c2 * t0 ) < 0 )
          result += 2 * aol::I * aol::NumberTrait<long double>::pi / c2 * jtsum;
      }
      if ( ( arg_start < 0 ) && ( arg_end >= 0 ) ) {
        t0 = - imag ( c1 ) / imag ( c2 );
        if ( real ( c1 + c2 * t0 ) < 0 )
          result -= 2 * aol::I * aol::NumberTrait<long double>::pi / c2 * jtsum;
      }
    }

    if ( aol::debugging::finiteness && !aol::isFinite ( real ( result ) ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;

    return real ( result ) / ( 2 * aol::NumberTrait<long double>::pi );
  }
private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp;
};

//! ElasticSolutionGrad Class
/** Gradient of Single Layer Boundary Integral Operator for piecewise constant ansatz functions.
 *  Overwrites LocalOperator::evaluateLocally()
 */
template <class _ParticleType>
class ElasticSolutionGrad
      : public LocalOperator<ElasticSolutionGrad<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  ElasticSolutionGrad ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp, int diffComp )
      : LocalOperator<ElasticSolutionGrad<ParticleType>, ParticleType > (),
      _greenCoeffs ( greenCoeffs ), _evalComp ( evalComp ), _intComp ( intComp ), _diffComp ( diffComp ) {}
  using LocalOperator<ElasticSolutionGrad<ParticleType>, ParticleType >::evaluateLocally;
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool, bool isinterior, bool ) const {
    int alpha, j;
    complex<DataType> jsum = 0, c1, c2, c_start, c_end, log_start, log_end, result = 0;
    DataType t0, arg_start, arg_end;
    aol::Vec2<DataType> dir, dirs, dire;

    for ( alpha = 0; alpha < 2; alpha++ ) {
      for ( j = 0, jsum = 0; j < 2; j++ )
        jsum += _greenCoeffs.A [_intComp][alpha] * _greenCoeffs.N [alpha][j] * _greenCoeffs.d [j][_evalComp];

      dirs = integrate.getStart () - evaluate;
      dire = integrate.getEnd () - evaluate;
      dir = integrate.getDirection ();
      dir.normalize ();
      aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );

      c_start = pav * aol::Vec2< std::complex<DataType> > ( dirs );
      c_end = pav * aol::Vec2< std::complex<DataType> > ( dire );

      log_start = log ( c_start );
      log_end = log ( c_end );
      arg_start = imag ( log_start );
      arg_end = imag ( log_end );

      c1 = c_start;
      c2 = pav * aol::Vec2< std::complex<DataType> > ( dir );

      if ( !isinterior ) {
        result -= ( log_end - log_start ) / c2 * jsum * pav [_diffComp];

        if ( ( arg_start >= 0 ) && ( arg_end < 0 ) ) { // Jump due to sign change of the imag part of the log argument
          t0 = - imag ( c1 ) / imag ( c2 );
          if ( real ( c1 + c2 * t0 ) < 0 )
            result -= 2 * aol::I * aol::NumberTrait<long double>::pi / c2 * jsum * pav [_diffComp];
        }
        if ( ( arg_start < 0 ) && ( arg_end >= 0 ) ) {
          t0 = - imag ( c1 ) / imag ( c2 );
          if ( real ( c1 + c2 * t0 ) < 0 )
            result += 2 * aol::I * aol::NumberTrait<long double>::pi / c2 * jsum * pav [_diffComp];
        }
      }
    }

    if ( aol::debugging::finiteness && !aol::isFinite ( real ( result ) ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;

    return - real ( result ) / ( 2 * aol::NumberTrait<long double>::pi );
  }

private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp, _diffComp;
};

//! ElasticSolutionDiffGrad Class
/** Gradient of Double Layer Boundary Integral Operator for piecewise constant ansatz functions.
 *  Overwrites LocalOperator::evaluateLocally()
 */
template <class _ParticleType>
class ElasticSolutionDiffGrad
      : public LocalOperator<ElasticSolutionDiffGrad<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  ElasticSolutionDiffGrad ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp, int diffComp )
      : LocalOperator<ElasticSolutionDiffGrad<ParticleType>, ParticleType > (),
      _greenCoeffs ( greenCoeffs ), _evalComp ( evalComp ), _intComp ( intComp ), _diffComp ( diffComp ) {}
  using LocalOperator<ElasticSolutionDiffGrad<ParticleType>, ParticleType >::evaluateLocally;
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool /*isstart*/, bool /*isinterior*/, bool /*isend*/ ) const {
    int alpha, j, t;
    complex<DataType> jtsum = 0, c1, c2, c_start, c_end, log_start, log_end, result = 0;
    aol::Vec2<DataType> dir, dirs, dire;

    for ( alpha = 0; alpha < 2; alpha++ ) {
      for ( j = 0, jtsum = 0; j < 2; j++ )
        for ( t = 0; t < 2; t++ )
          jtsum += _greenCoeffs.L [_intComp][j][alpha] * _greenCoeffs.N [alpha][t] * integrate.getNormal () [j] * _greenCoeffs.d [t][_evalComp];

      dirs = integrate.getStart () - evaluate;
      dire = integrate.getEnd () - evaluate;
      aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );

      c_start = pav * aol::Vec2< complex< DataType> > ( dirs );
      c_end = pav * aol::Vec2< complex< DataType> > ( dire );
      c1 = c_start;

      result += jtsum * integrate.getLength () / ( c1 * c_end ) * pav [_diffComp];
    }

    if ( aol::debugging::finiteness && !aol::isFinite ( real ( result ) ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;

    return real ( result ) / ( 2 * aol::NumberTrait<long double>::pi );
  }
private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp, _diffComp;
};


//! OffCenterElasticSolutionGrad Class
/** Gradient of Single Layer Boundary Integral Operator for piecewise constant ansatz functions
 *  with shifted singularity of the kernel.
 *  Overwrites LocalOperator::evaluateLocally()
 */
template <class _ParticleType>
class OffCenterElasticSolutionGrad
      : public LocalOperator<ElasticSolutionGrad<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  OffCenterElasticSolutionGrad ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp, int diffComp,
                                 double offset = 0.5 )
      : LocalOperator<ElasticSolutionGrad<ParticleType>, ParticleType > (),
      _greenCoeffs ( greenCoeffs ), _evalComp ( evalComp ), _intComp ( intComp ), _diffComp ( diffComp ), _offset ( offset ) {

    for ( int alpha = 0; alpha < 2; alpha++ )
      for ( int j = 0; j < 2; j++ )
        ANd [alpha] += _greenCoeffs.A [_intComp][alpha] * _greenCoeffs.N [alpha][j] * _greenCoeffs.d [j][_evalComp];
  }
  using LocalOperator<ElasticSolutionGrad<ParticleType>, ParticleType >::evaluateLocally;
  DataType evaluateLocally ( const ConstSegmentType& integrate, const ConstSegmentType& evaluate ) const {
    aol::Vec2<DataType> point = evaluate.getStart () + _offset * evaluate.getDirection (); // off center collocation
    return evaluateLocally ( integrate, point, false, integrate == evaluate, false );
  }
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool, bool isinterior, bool ) const {
    int alpha;
    complex<DataType> c1, c2, c_start, c_end, log_start, log_end, result = 0;
    DataType t0, arg_start, arg_end;
    aol::Vec2<DataType> dir, dirs, dire;

    dir = integrate.getDirection ();
    dirs = integrate.getStart () - evaluate;
    dire = integrate.getEnd () - evaluate;

    if ( isinterior ) {
      if ( _offset >= 0.5 ) dire -= dir * 2 * ( 1 - _offset );
      if ( _offset < 0.5 ) dirs += dir * 2 * _offset;
    }

    dir.normalize ();

    for ( alpha = 0; alpha < 2; alpha++ ) {

      aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );

      c_start = pav * aol::Vec2<complex<DataType> >( dirs );
      c_end = pav * aol::Vec2<complex<DataType> >( dire );

      log_start = log ( c_start );
      log_end = log ( c_end );
      arg_start = imag ( log_start );
      arg_end = imag ( log_end );

      c1 = c_start;
      c2 = pav * aol::Vec2<complex<DataType> >( dir );

      result -= ( log_end - log_start ) / c2 * ANd [alpha] * pav [_diffComp];

      if ( arg_start * arg_end <= 0 ) {
  if ( arg_end < 0 ) { // Jump due to sign change of the imag part of the log argument
    t0 = - imag ( c1 ) / imag ( c2 );
    if ( real ( c1 + c2 * t0 ) < 0 )
      result -= 2 * aol::I * aol::NumberTrait<long double>::pi / c2 * ANd [alpha] * pav [_diffComp];
  }
  if ( arg_start < 0 ) {
    t0 = - imag ( c1 ) / imag ( c2 );
    if ( real ( c1 + c2 * t0 ) < 0 )
      result += 2 * aol::I * aol::NumberTrait<long double>::pi / c2 * ANd [alpha] * pav [_diffComp];
  }
      }
    }

    if ( aol::debugging::finiteness && !aol::isFinite ( real ( result ) ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;

    return - real ( result ) / ( 2 * aol::NumberTrait<long double>::pi );
  }

private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp, _diffComp;
  double _offset;
  aol::Vec2<complex<DataType> > ANd;
};


//! LinElasticSolution Class
/** Single Layer Boundary Integral Operator for piecewise linear ansatz functions.
 *  Overwrites LinLocalOperator::evaluateLocally()
 */
template <class _ParticleType>
class LinElasticSolution
      : public LinLocalOperator<LinElasticSolution<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  LinElasticSolution ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp )
      : LinLocalOperator<LinElasticSolution<ParticleType>, ParticleType > ()
      , _greenCoeffs ( greenCoeffs )
      , _evalComp ( evalComp )
      , _intComp ( intComp )
  {}

  using LinLocalOperator<LinElasticSolution<ParticleType>, ParticleType >::evaluateLocally;

  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool, bool isend, bool firstHalf ) const {
    int alpha, j;
    std::complex<DataType> alphasum, elementsum = 0, c1, c2, c2h, dummy1, h, result = 0;
    DataType dirscal = firstHalf ? 1 : -1;

    aol::Vec2<DataType> collpt = firstHalf ? integrate.getEnd () : integrate.getStart ();
    aol::Vec2<DataType> dirc   = evaluate - collpt;
    aol::Vec2<DataType> dir    = integrate.getDirection ();
    dir.normalize ();
    h = integrate.getLength ();

    for ( j = 0; j < 2; ++j ) {
      for ( alpha = 0, alphasum = 0; alpha < 2; ++alpha ) {

        elementsum = 0;
        aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );
        c1 = pav * aol::Vec2< complex< DataType> >( dirc );
        c2 = pav * aol::Vec2< complex< DataType> >( dir ) * dirscal;
        c2h = c2 * h;
        dummy1 = 0.5 * c1 * c1 / c2h;

        if ( !isstart && !isend )
          elementsum += -0.5 * c1 / c2 - 0.75 * h + 0.5 * h * log ( c1 + c2 * h ) + ( c1 + dummy1 ) * log ( 1 + c2 * h / c1 ) / c2;
        else if ( ( !firstHalf && isend ) || ( firstHalf && isstart ) )
          elementsum += -0.25*h + 0.5*h*log ( -c2h );
        else if ( ( firstHalf && isend ) || ( !firstHalf && isstart ) )
          elementsum += -0.75*h + 0.5*h*log ( c2h );

        alphasum += elementsum * _greenCoeffs.A [_intComp][alpha] * _greenCoeffs.N [alpha][j] * _greenCoeffs.d [j][_evalComp];
      }
      result += alphasum;
    }
    if ( aol::debugging::finiteness && !aol::isFinite ( real ( result ) ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;
    return real ( result ) / ( aol::NumberTrait<long double>::pi * 2 );
  }
private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp;
};


//! LinJumpElasticSolution Class
/** LinElasticSolution for boundaries with jumps
 *  Overwrites LinLocalOperator::evaluateLocally()
 *
 *  @todo ML,BG: make nicer / review
 */
template <class _ParticleType>
class LinJumpElasticSolution
      : public LinJumpLocalOperator<LinJumpElasticSolution<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  LinJumpElasticSolution ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp )
      : LinJumpLocalOperator<LinJumpElasticSolution<ParticleType>, ParticleType > (),
      _greenCoeffs ( greenCoeffs ), _evalComp ( evalComp ), _intComp ( intComp ) {}
  using LinJumpLocalOperator<LinJumpElasticSolution<ParticleType>, ParticleType >::evaluateLocally;

  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool, bool isend, bool firstHalf ) const {
    int alpha, j;
    std::complex<DataType> alphasum, elementsum = 0, c1, c2, c2h, dummy1, dummy2, h;
    DataType dirscal, result = 0;

    dirscal = firstHalf ? -1 : 1;

    for ( j = 0; j < 2; ++j ) {
      for ( alpha = 0, alphasum = 0; alpha < 2; ++alpha ) {

        elementsum = 0;
        aol::Vec2<DataType> collpt = firstHalf ? integrate.getEnd () : integrate.getStart ();
        aol::Vec2<DataType> dirc = collpt - evaluate;
        aol::Vec2<DataType> dir = integrate.getDirection ();
        dir.normalize ();
        aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );
        h = integrate.getLength ();
        c1 = pav * aol::Vec2< std::complex<DataType> > ( dirc );
        c2 = pav * aol::Vec2< std::complex<DataType> > ( dir ) * dirscal;
        c2h = c2 * h;
        elementsum -= 0.5 * c1 / c2 + 0.75 * h;
        dummy1 = 0.5 * c1 * c1 / c2h;
        dummy2 = 0.5 * c2h;

        if ( !isstart && !isend )
          elementsum += 0.5 * h * log ( c1 + c2 * h ) + ( c1 + dummy1 ) * log ( 1 + c2 * h / c1 ) / c2;
        else if ( ( !firstHalf && isend ) || ( firstHalf && isstart ) )
          elementsum -= ( c1 + dummy1 ) * log ( c1 ) / c2;
        else if ( ( firstHalf && isend ) || ( !firstHalf && isstart ) )
          elementsum += ( c1 + dummy1 + dummy2 ) * log ( c1 + c2 * h ) / c2;

        alphasum += elementsum * _greenCoeffs.A [_intComp][alpha] * _greenCoeffs.N [alpha][j];
      }
      result += real ( alphasum ) * _greenCoeffs.d [j][_evalComp];
    }
    if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;
    return result / ( aol::NumberTrait<long double>::pi * 2 );
  }
private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp;
};


//! LinElasticSolutionDiff Class
/** Double Layer Boundary Integral Operator for piecewise linear ansatz functions.
 *  Overwrites LinLocalOperator::evaluateLocally()
 *  @warning This only works when the diagonal gets normalized by utilizing the rigid body argument
 *  @see     ElasticOperator
 */
template <class _ParticleType>
class LinElasticSolutionDiff
      : public LinLocalOperator<LinElasticSolutionDiff<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  LinElasticSolutionDiff ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp )
      : LinLocalOperator<LinElasticSolutionDiff<ParticleType>, ParticleType > ()
      , _greenCoeffs ( greenCoeffs )
      , _evalComp ( evalComp )
      , _intComp ( intComp )
  {}

  using LinLocalOperator<LinElasticSolutionDiff<ParticleType>, ParticleType >::evaluateLocally;

  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool, bool isend, bool firstHalf ) const {
    int alpha, j, t;
    std::complex<DataType> alphasum, elementsum = 0, c1, c2, c2h, dummy;
    DataType h, result = 0, dirscal = firstHalf ? 1 : -1;

    aol::Vec2<DataType> collpt = firstHalf ? integrate.getEnd () : integrate.getStart ();
    aol::Vec2<DataType> dirc   = evaluate - collpt;
    aol::Vec2<DataType> norm   = integrate.getNormal ();
    aol::Vec2<DataType> dir    = integrate.getDirection ();
    dir.normalize ();
    h = integrate.getLength ();

    for ( j = 0; j < 2; ++j ) {
      for ( t = 0; t < 2; ++t ) {
        for ( alpha = 0, alphasum = 0; alpha < 2; ++alpha ) {

          elementsum = 0;
          aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );
          c1 = pav * aol::Vec2< complex< DataType> >( dirc );
          c2 = pav * aol::Vec2< complex< DataType> >( dir ) * dirscal;
          c2h = c2 * h;
          dummy = c1 / c2h;

          // -1/c2 can be dropped here
          if ( !isstart && !isend )
            elementsum += ( 1 + dummy ) * log ( 1 + 1 / dummy ) / c2;
          // else 0

          alphasum += elementsum * _greenCoeffs.L [_intComp][j][alpha] * _greenCoeffs.N [alpha][t];
        }
        result += real ( alphasum ) * _greenCoeffs.d [t][_evalComp] * norm [j];
      }
    }
    if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;
    return -result / ( aol::NumberTrait<long double>::pi * 2 );
  }
private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp;
};

//! LinElasticSolutionGrad Class
/** Gradient of Single Layer Boundary Integral Operator for piecewise linear ansatz functions.
 *  Overwrites LinLocalOperator::evaluateLocally()
 *  @warning There are terms missing in the else cases!
 */
template <class _ParticleType>
class LinElasticSolutionGrad
      : public LinLocalOperator<LinElasticSolutionGrad<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  LinElasticSolutionGrad ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp, int diffComp )
      : LinLocalOperator<LinElasticSolutionGrad<ParticleType>, ParticleType > ()
      , _greenCoeffs ( greenCoeffs )
      , _evalComp ( evalComp )
      , _intComp ( intComp )
      , _diffComp ( diffComp )
  {}

  using LinLocalOperator<LinElasticSolutionGrad<ParticleType>, ParticleType >::evaluateLocally;

  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool, bool isend, bool firstHalf ) const {
    int alpha, j;
    std::complex<DataType> alphasum, elementsum = 0, c1, c2, c2h, dummy;
    DataType h, result = 0, dirscal = firstHalf ? 1 : -1;

    aol::Vec2<DataType> collpt = firstHalf ? integrate.getEnd () : integrate.getStart ();
    aol::Vec2<DataType> dirc   = evaluate - collpt;
    aol::Vec2<DataType> dir    = integrate.getDirection ();
    dir.normalize ();
    h = integrate.getLength ();

    for ( j = 0; j < 2; ++j ) {
      for ( alpha = 0, alphasum = 0; alpha < 2; ++alpha ) {

        elementsum = 0;
        aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );
        c1 = pav * aol::Vec2<complex<DataType> >( dirc );
        c2 = pav * aol::Vec2<complex<DataType> >( dir ) * dirscal;
        c2h = c2 * h;
        dummy = c1 / c2h;

        if ( !isstart && !isend )
          elementsum += ( ( 1 + dummy ) * log ( 1 + 1 / dummy ) - 1 ) / c2;  // 1/c2 does not vanish here
        else if ( ( !firstHalf && isend ) || ( firstHalf && isstart ) )
          elementsum -= 1 / c2;
        else if ( ( firstHalf && isend ) || ( !firstHalf && isstart ) )
          elementsum += ( log ( c2h ) - 1 ) / c2;

        alphasum += elementsum * _greenCoeffs.A [_intComp][alpha] * _greenCoeffs.N [alpha][j] * pav [_diffComp];
      }
      result += real ( alphasum ) * _greenCoeffs.d [j][_evalComp];
    }

    if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;
    return result / ( aol::NumberTrait<long double>::pi * 2 );
  }
private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp, _diffComp;
};

//! LinElasticSolutionDiffGrad Class
/** Gradient of Double Layer Boundary Integral Operator for piecewise linear ansatz functions.
 *  Overwrites LinLocalOperator::evaluateLocally()
 *  @warning There are terms missing in the else cases!
 */
template <class _ParticleType>
class LinElasticSolutionDiffGrad
      : public LinLocalOperator<LinElasticSolutionDiffGrad<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType                           ParticleType;
  typedef typename ParticleType::DataType         DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  LinElasticSolutionDiffGrad ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp, int diffComp )
      : LinLocalOperator<LinElasticSolutionDiffGrad<ParticleType>, ParticleType > ()
      , _greenCoeffs ( greenCoeffs )
      , _evalComp ( evalComp )
      , _intComp ( intComp )
      , _diffComp ( diffComp )
  {}

  using LinLocalOperator<LinElasticSolutionDiffGrad<ParticleType>, ParticleType >::evaluateLocally;

  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool /*isinterior*/, bool isend, bool firstHalf ) const {
    int alpha, j, t;
    std::complex<DataType> alphasum, elementsum = 0, c1, c2, c2h;
    DataType h, result = 0, dirscal = firstHalf ? 1 : -1;

    double eps = 0.252646e-5;

    aol::Vec2<DataType> collpt = firstHalf ? integrate.getEnd () : integrate.getStart ();
    aol::Vec2<DataType> dirc   = evaluate - collpt;
    aol::Vec2<DataType> norm   = integrate.getNormal ();
    aol::Vec2<DataType> dir    = integrate.getDirection ();
    dir.normalize ();
    h = integrate.getLength ();

    for ( j = 0; j < 2; ++j ) {
      for ( t = 0; t < 2; ++t ) {
        for ( alpha = 0, alphasum = 0; alpha < 2; ++alpha ) {
          elementsum = 0;
          aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );
          c1 = pav * aol::Vec2< std::complex<DataType> > ( dirc );
          c2 = pav * aol::Vec2< std::complex<DataType> > ( dir ) * dirscal;
          c2h = c2 * h;

          //if ( isinterior )                                                          // ???
          //  {elementsum += 1 / ( c1 * c2 );} // part finit
          if ( !isstart && !isend )                                                    // regular case, OK!
            {elementsum -= ( log ( 1 + c2h / c1 ) / c2h - 1 / c1 ) / c2;}
          else if ( ( !firstHalf && isend ) || ( firstHalf && isstart ) )  // singularity where base funtions becomes 0
            {elementsum += (log ( c2h*c2h )-log(eps*c2h*c2) ) / ( c2h * c2 );}                                  // was log ( c1 ), log(h)
          else if ( ( firstHalf && isend ) || ( !firstHalf && isstart ) )  // singularity at center of base function, OK CPV?
            {elementsum -= ( log ( h ) + 1 ) / ( c2h * c2 );}                          // was log ( c2h )

          alphasum += elementsum * _greenCoeffs.L [_intComp][j][alpha] * _greenCoeffs.N [alpha][t] * pav [_diffComp];
        }
        result += real ( alphasum ) * _greenCoeffs.d [t][_evalComp] * norm[j];
      }
    }
    if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;
    return result / ( aol::NumberTrait<long double>::pi * 2 );
  }
private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp, _diffComp;
};


//! OffCenterLinElasticSolutionGrad Class
/** Gradient of Single Layer Boundary Integral Operator for piecewise linear ansatz functions
 *  with shifted singularity of the kernel.
 *  Overwrites LinLocalOperator::evaluateLocally()
 */
template <class _ParticleType>
class OffCenterLinElasticSolutionGrad
      : public LinLocalOperator<OffCenterLinElasticSolutionGrad<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  OffCenterLinElasticSolutionGrad ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp, int diffComp, double offset = 0.5 )
      : LinLocalOperator<OffCenterLinElasticSolutionGrad<ParticleType>, ParticleType > ()
      , _greenCoeffs ( greenCoeffs )
      , _evalComp ( evalComp )
      , _intComp ( intComp )
      , _diffComp ( diffComp )
      , _offset ( offset )
  {}

  using LinLocalOperator<OffCenterLinElasticSolutionGrad<ParticleType>, ParticleType >::evaluateLocally;

  DataType evaluateLocally ( const ConstSegmentType& integrate, const ConstSegmentType& evaluate, bool firstHalf ) const {
    aol::Vec2<DataType> point = evaluate.getStart () + _offset * evaluate.getDirection (); // off center collocation
    return evaluateLocally ( integrate, point, false, integrate == evaluate , false, firstHalf );
  }

  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool, bool isinterior, bool, bool firstHalf ) const {
    int alpha, j;
    std::complex<DataType> alphasum, elementsum = 0, c1, c2, c2h, dummy, h;
    DataType dirscal, result = 0;

    dirscal = firstHalf ? 1 : -1;
    aol::Vec2<DataType> collpt = firstHalf ? integrate.getEnd () : integrate.getStart ();
    aol::Vec2<DataType> dirc   = evaluate - collpt;
    aol::Vec2<DataType> dir    = integrate.getDirection (); dir.normalize ();
    h = integrate.getLength ();

    for ( j = 0; j < 2; ++j ) {
      for ( alpha = 0, alphasum = 0; alpha < 2; ++alpha ) {

        elementsum = 0;
        aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );
        c1  = pav * aol::Vec2<complex<DataType> > ( dirc );
        c2  = pav * aol::Vec2<complex<DataType> > ( dir ) * dirscal;
        c2h = c2 * h;
        dummy = c1 / c2h;

        if ( !isinterior )
          elementsum += ( ( 1.0 + dummy ) * log ( 1.0 + 1.0 / dummy ) - 1.0 ) / c2;
        else
          elementsum += -1.0/c2; // --> try 0 ?

        alphasum += elementsum * _greenCoeffs.A [_intComp][alpha] * _greenCoeffs.N [alpha][j] * pav [_diffComp];
      }
      result += real ( alphasum ) * _greenCoeffs.d [j][_evalComp];
    }

    if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;
    return result / ( aol::NumberTrait<long double>::pi * 2 );
  }
private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp, _diffComp;
  double _offset;
};

//! OffCenterLinElasticSolutionDiffGrad Class
/** Gradient of Double Layer Boundary Integral Operator for piecewise linear ansatz functions
 *  with shifted singularity of the kernel.
 *  Overwrites LinLocalOperator::evaluateLocally()
 *  @warning There are terms missing in the else case!
 */
template <class _ParticleType>
class OffCenterLinElasticSolutionDiffGrad
      : public LinLocalOperator<OffCenterLinElasticSolutionDiffGrad<_ParticleType>, _ParticleType > {
public:
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  OffCenterLinElasticSolutionDiffGrad ( const ElasticGreen<DataType>& greenCoeffs, int evalComp, int intComp, int diffComp, double offset = 0.5 )
      : LinLocalOperator<OffCenterLinElasticSolutionDiffGrad<ParticleType>, ParticleType > ()
      , _greenCoeffs ( greenCoeffs )
      , _evalComp ( evalComp )
      , _intComp ( intComp )
      , _diffComp ( diffComp )
      , _offset ( offset )
  {}

  using LinLocalOperator<OffCenterLinElasticSolutionDiffGrad<ParticleType>, ParticleType >::evaluateLocally;

  DataType evaluateLocally ( const ConstSegmentType& integrate, const ConstSegmentType& evaluate, bool firstHalf ) const {
    aol::Vec2<DataType> point = evaluate.getStart () + _offset * evaluate.getDirection (); // off center collocation
    return evaluateLocally ( integrate, point, false, integrate == evaluate , false, firstHalf );
  }

  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool, bool isinterior, bool, bool firstHalf ) const {
    int alpha, j, t;
    std::complex<DataType> alphasum, elementsum = 0, c1, c2, c2h, dummy, h;
    DataType dirscal, result = 0;

    dirscal = firstHalf ? 1 : -1;
    aol::Vec2<DataType> collpt = firstHalf ? integrate.getEnd () : integrate.getStart ();
    aol::Vec2<DataType> dirc   = evaluate - collpt;
    aol::Vec2<DataType> norm   = integrate.getNormal ();
    aol::Vec2<DataType> dir    = integrate.getDirection (); dir.normalize ();
    h = integrate.getLength ();

    for ( j = 0; j < 2; ++j ) {
      for ( t = 0; t < 2; ++t ) {
        for ( alpha = 0, alphasum = 0; alpha < 2; ++alpha ) {
          elementsum = 0;
          aol::Vec2<complex<DataType> > pav ( 1, _greenCoeffs.p [alpha] );
          c1  = pav * aol::Vec2<complex<DataType> > ( dirc );
          c2  = pav * aol::Vec2<complex<DataType> > ( dir ) * dirscal;
          c2h = c2 * h;
          dummy = c2h / c1;

          if ( !isinterior )
            elementsum += ( 1.0/c1 - log ( 1 + dummy ) / c2h ) / c2;
          else
            elementsum += -2.0/(c2h*c2);

          alphasum += elementsum * _greenCoeffs.L [_intComp][j][alpha] * _greenCoeffs.N [alpha][t] * pav [_diffComp];
        }
        result += real ( alphasum ) * _greenCoeffs.d [t][_evalComp] * norm[j];
      }
    }

    if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
      cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;
    return result / ( aol::NumberTrait<long double>::pi * 2 );
  }
private:
  ElasticGreen<DataType> _greenCoeffs;
  int _evalComp, _intComp, _diffComp;
  double _offset;
};

}

#endif
