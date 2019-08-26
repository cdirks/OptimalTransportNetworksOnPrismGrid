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

#ifndef __LAPLACEFUNSOL_H
#define __LAPLACEFUNSOL_H

#include <bemOps.h>
#include <boundary.h>

namespace bm {

/************************/
//! Single layer potential
template <class _ParticleType>
class SingleLayerPotential
      : public LocalOperator<SingleLayerPotential<_ParticleType>, _ParticleType > {

public:

  //! Types
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  using LocalOperator<SingleLayerPotential<ParticleType>, ParticleType >::evaluate;
  using LocalOperator<SingleLayerPotential<ParticleType>, ParticleType >::evaluateLocally;

  //! Evaluates the kernel
  DataType evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center ) const;

  //! Evaluates operator locally at a point, integrated over a segment
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool isinterior, bool isend ) const;
};

/*********************************/
//! Single layer potential Gradient
template <class _ParticleType>
class SingleLayerPotentialGrad
      : public LocalOperator<SingleLayerPotentialGrad<_ParticleType>, _ParticleType > {

public:
  SingleLayerPotentialGrad ( int direction ) : _direction ( direction ) { }

  //! Types
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  using LocalOperator<SingleLayerPotentialGrad<ParticleType>, ParticleType >::evaluate;
  using LocalOperator<SingleLayerPotentialGrad<ParticleType>, ParticleType >::evaluateLocally;

  //! Evaluates the kernel
  DataType evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center ) const;

  //! Evaluates operator locally at a point, integrated over a segment
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool isinterior, bool isend ) const;

private:
  const int _direction;
};

/************************/
//! Double layer potential
template <class _ParticleType>
class DoubleLayerPotential
      : public LocalOperator<DoubleLayerPotential<_ParticleType>, _ParticleType > {

public:

  //! Types
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  using LocalOperator<DoubleLayerPotential<ParticleType>, ParticleType >::evaluate;
  using LocalOperator<DoubleLayerPotential<ParticleType>, ParticleType >::evaluateLocally;

  //! Evaluates the kernel
  DataType evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center, aol::Vec2<DataType> normal ) const;

  //! Evaluates operator locally at a point, integrated over a segment
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool isinterior, bool isend ) const;
};

/************************/
//! Single layer potential
template <class _ParticleType>
class LinSingleLayerPotential
      : public LinLocalOperator<LinSingleLayerPotential<_ParticleType>, _ParticleType > {

public:

  //! Types
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  using LinLocalOperator<LinSingleLayerPotential<ParticleType>, ParticleType >::evaluate;
  using LinLocalOperator<LinSingleLayerPotential<ParticleType>, ParticleType >::evaluateLocally;

  //! Evaluates the kernel
  DataType evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center ) const;

  //! Evaluates operator locally at a point, integrated over a segment
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool isinterior, bool isend, bool firstHalf ) const;

private:
  SingleLayerPotential<ParticleType> constOp;
};

/************************/
//! Double layer potential
template <class _ParticleType>
class LinDoubleLayerPotential
      : public LinLocalOperator<LinDoubleLayerPotential<_ParticleType>, _ParticleType > {

public:

  //! Types
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  using LinLocalOperator<LinDoubleLayerPotential<ParticleType>, ParticleType >::evaluate;
  using LinLocalOperator<LinDoubleLayerPotential<ParticleType>, ParticleType >::evaluateLocally;

  //! Evaluates the kernel
  DataType evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center, aol::Vec2<DataType> normal ) const;

  //! Evaluates operator locally at a point, integrated over a segment
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool isinterior, bool isend, bool firstHalf ) const;

private:
  DoubleLayerPotential<ParticleType> constOp;
};

  /************************************/
//! Derivative of single layer potential
//! @remark incomplete, not working
//! @todo ML: complete
template <class _ParticleType>
class LinSingleLayerPotentialGrad
      : public LinLocalOperator<LinSingleLayerPotentialGrad<_ParticleType>, _ParticleType > {

public:

  //! Types
  typedef _ParticleType ParticleType;
  typedef typename ParticleType::DataType DataType;
  typedef typename ParticleType::ConstSegmentType ConstSegmentType;

  using LinLocalOperator<LinSingleLayerPotentialGrad<ParticleType>, ParticleType >::evaluate;
  using LinLocalOperator<LinSingleLayerPotentialGrad<ParticleType>, ParticleType >::evaluateLocally;

  //! Constructor, dir = 0 for x-direction, dir = 1 for y-direction
  LinSingleLayerPotentialGrad ( int dir );

  //! Evaluates the kernel
  DataType evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center ) const;

  //! Evaluates operator locally at a point, integrated over a segment
  DataType evaluateLocally ( const ConstSegmentType& integrate, const aol::Vec2<DataType>& evaluate,
                             bool isstart, bool isinterior, bool isend, bool firstHalf ) const;

  //! Evaluates operator locally at a point, integrated over a segment.
  //! Basis function is zero at point a and 1 at point b.
  //! Considers only the case where coll != a and coll != b.
  DataType evaluateLocallyRegular ( const aol::Vec2<DataType>& coll, aol::Vec2<DataType> a, aol::Vec2<DataType> b ) const;

private:

  // Direction of differentiation
  int _dir;
};



// Class SingleLayerPotential<ParticleType>

template <class ParticleType>
typename ParticleType::DataType SingleLayerPotential<ParticleType>::evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center ) const {
  point -= center;
  return - log ( point.norm () ) * ( 0.5 / aol::NumberTrait<long double>::pi );
}

template <class ParticleType>
typename ParticleType::DataType SingleLayerPotential<ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
    const aol::Vec2<DataType>& evaluate,
    bool isstart, bool isinterior, bool isend ) const {
  if ( ( isstart || isend ) && !isinterior )
    throw aol::ParameterException ( "SingleLayerPotential::evaluateLocally, start or end points are also considered interior points", __FILE__, __LINE__ );

  aol::Vec2<DataType> a, b, dir, normal;
  dir = integrate.getDirection ();
  normal = integrate.getNormal ();
  a = evaluate - integrate.getStart ();
  b = a - dir;

  DataType len = integrate.getLength ();
  DataType result = - len;

  DataType alen = 0, atilde, blen = 0, btilde; // Avoid uninitialized warning
  DataType c, cosphi, phi;

  if ( !isstart ) {
    alen = a.norm ();
    atilde = a * dir / len;

    result += atilde * log ( alen );
  }

  if ( !isend ) {
    blen = b.norm ();
    btilde = - b * dir / len;

    result += btilde * log ( blen );
  }

  if ( !isstart && !isinterior && !isend ) {
    c = fabs ( a * normal );
    cosphi = a * b / ( alen * blen );
    if ( cosphi >= 1.0 ) phi = 0.0;
    else phi = acos ( cosphi );

    result += c * phi;
  }

  if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
    cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;

  return - 0.5 * result / aol::NumberTrait<long double>::pi;
}

// Class SingleLayerPotentialGrad<ParticleType>

template <class ParticleType>
typename ParticleType::DataType SingleLayerPotentialGrad<ParticleType>::evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center ) const {
  point -= center;
  return - 0.5 / aol::NumberTrait<long double>::pi * (point [_direction]) / point.normSqr ();
}

template <class ParticleType>
typename ParticleType::DataType SingleLayerPotentialGrad<ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
    const aol::Vec2<DataType>& evaluate,
    bool isstart, bool isinterior, bool isend ) const {

  if ( isstart || isinterior || isend ) return 0; // really! Symmetry...

  aol::Vec2<double> a = evaluate - integrate.getStart ();
  aol::Vec2<double> b = a - integrate.getDirection();

  aol::Vec2<double> dir = integrate.getDirection() / integrate.getLength ();
  aol::Vec2<double> normal = integrate.getNormal();

  double A = - a [_direction], D = dir [_direction];
  double alen = a.norm(), blen = b.norm ();
  double atilde = a * dir;
  double c = abs ( a * normal );

  double len = integrate.getLength ();

  double result = (A+atilde*D)/c * (atan ( (len - atilde) / c ) - atan ( - atilde / c ))  + D * (log(blen) - log(alen));

  if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
    cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;

  return 0.5 * result / aol::NumberTrait<long double>::pi;
}

// Class DoubleLayerPotential<ParticleType>

template <class ParticleType>
typename ParticleType::DataType DoubleLayerPotential<ParticleType>::evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center, aol::Vec2<DataType> normal ) const {
  point -= center;
  return  0.5 / aol::NumberTrait<long double>::pi * (point * normal) / point.normSqr ();
}

template <class ParticleType>
typename ParticleType::DataType DoubleLayerPotential<ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
    const aol::Vec2<DataType>& evaluate,
    bool isstart, bool isinterior, bool isend ) const {
  if ( ( isstart || isend ) && !isinterior )
    throw aol::ParameterException ( "DoubleLayerPotential::evaluateLocally, start or end points are also considered interior points", __FILE__, __LINE__ );

  if ( isstart || isinterior || isend ) return 0;

  aol::Vec2<DataType> a, dir, normal;
  dir = integrate.getDirection ();
  normal = integrate.getNormal ();
  a = evaluate - integrate.getStart ();

  DataType len = integrate.getLength ();
  DataType atilde = a * dir / len;
  DataType c = abs ( a * normal );
  DataType alpha = atan2 ( - atilde, c ), beta = atan2 ( len - atilde, c );

  DataType result = beta - alpha;

  if ( a * normal > 0 ) result *= -1;

  if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
    cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;

  return 0.5 * result / aol::NumberTrait<long double>::pi;
}

// Class LinSingleLayerPotential<ParticleType>

template <class ParticleType>
typename ParticleType::DataType LinSingleLayerPotential<ParticleType>::evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center ) const {
  point -= center;
  return - log ( point.norm () ) * ( 0.5 / aol::NumberTrait<long double>::pi );
}

template <class ParticleType>
typename ParticleType::DataType LinSingleLayerPotential<ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
                       const aol::Vec2<DataType>& evaluate,
                       bool isstart, bool isinterior, bool isend, bool firstHalf ) const {
  if ( ( isstart || isend ) && !isinterior )
    throw aol::ParameterException ( "LinSingleLayerPotential::evaluateLocally, start or end points are also considered interior points", __FILE__, __LINE__ );

  aol::Vec2<DataType> a, b, dir, normal;
  DataType alen = 0, atilde, blen = 0; // Avoid uninitialized warning
  DataType c, cosphi, phi, len = integrate.getLength ();

  aol::Vec2<DataType> collpt = integrate.getStart ();
  dir = integrate.getDirection ();
  normal = integrate.getNormal ();
  a = evaluate - collpt;
  b = a - dir;
  alen = a.norm ();
  atilde = a * dir / len;
  c = fabs ( a * normal );

  DataType result = - len * ( 2 * atilde + len );

  if ( !isstart ) {

    result -= 2 * ( c * c - atilde * atilde ) * log ( alen );
  }

  if ( !isend ) {
    blen = b.norm ();

    result += 2 * ( c * c - atilde * atilde + len * len ) * log ( blen );
  }

  if ( !isstart && !isinterior && !isend ) {
    cosphi = a * b / ( alen * blen );
    if ( cosphi >= 1.0 ) phi = 0.0;
    else phi = acos ( cosphi );

    result += 4 * atilde * c * phi;
  }

  if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
    cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;

  if ( firstHalf ) return - 0.25 * 0.5 * result / ( aol::NumberTrait<long double>::pi * len );
  else return 0.25 * 0.5 * result / ( aol::NumberTrait<long double>::pi * len ) + constOp.evaluateLocally ( integrate, evaluate, isstart, isinterior, isend );
}

// Class LinDoubleLayerPotential<ParticleType>

template <class ParticleType>
typename ParticleType::DataType LinDoubleLayerPotential<ParticleType>::evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center, aol::Vec2<DataType> normal ) const {
  point -= center;
  return  0.5 / aol::NumberTrait<long double>::pi * (point * normal) / point.normSqr ();
}

template <class ParticleType>
typename ParticleType::DataType LinDoubleLayerPotential<ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
                       const aol::Vec2<DataType>& evaluate,
                       bool isstart, bool isinterior, bool isend, bool firstHalf ) const {
  if ( ( isstart || isend ) && !isinterior )
    throw aol::ParameterException ( "LinDoubleLayerPotential::evaluateLocally, start or end points are also considered interior points", __FILE__, __LINE__ );

  if (isstart || isend) return 0;

  aol::Vec2<DataType> a, b, dir, normal;
  dir = integrate.getDirection ();
  normal = integrate.getNormal ();
  a = evaluate - integrate.getStart ();
  b = evaluate - integrate.getEnd ();

  DataType len = integrate.getLength ();
  DataType atilde = a * dir / len;
  DataType c = abs (a * normal);
  DataType alpha = atan2 ( - atilde, c ), beta = atan2 ( len - atilde, c );

  DataType result = a * normal * (log (b.norm()) - log(a.norm()));
  if ( c > 1e-24) result += a * normal * atilde / c * (beta - alpha);

  if ( a * normal > 0 ) result *= -1;

  if ( aol::debugging::finiteness && !aol::isFinite ( result ) )
    cerr << "Integral over " << integrate << " on " << evaluate << " is not finite." << endl;

  if ( firstHalf ) return 0.5 * result / ( aol::NumberTrait<long double>::pi * len );
  else return - 0.5 * result / ( aol::NumberTrait<long double>::pi * len ) - constOp.evaluateLocally ( integrate, evaluate, isstart, isinterior, isend );
}

// Class LinSingleLayerPotentialGrad<ParticleType>

template <class ParticleType>
LinSingleLayerPotentialGrad<ParticleType>::LinSingleLayerPotentialGrad ( int dir )
  : LinLocalOperator<LinSingleLayerPotentialGrad<ParticleType>,ParticleType> (), _dir ( dir ) {}

template <class ParticleType>
typename ParticleType::DataType LinSingleLayerPotentialGrad<ParticleType>::evaluate ( aol::Vec2<DataType> point, aol::Vec2<DataType> center ) const {
  point -= center;
  return /*- point [_dir] /*/ 1 / point.normSqr () * ( 0.5 / aol::NumberTrait<DataType>::pi );
}

template <class ParticleType>
typename ParticleType::DataType LinSingleLayerPotentialGrad<ParticleType>::evaluateLocally ( const ConstSegmentType& integrate,
                       const aol::Vec2<DataType>& evaluate,
                       bool isstart, bool, bool isend, bool firstHalf ) const {

  if (isstart) { return 0; } //!!
  if (isend)   { return 0; } //!!

  if (firstHalf) return evaluateLocallyRegular ( evaluate, integrate.getStart (), integrate.getEnd ());
  else return evaluateLocallyRegular ( evaluate, integrate.getEnd (), integrate.getStart ());
}

template <class ParticleType>
typename ParticleType::DataType LinSingleLayerPotentialGrad<ParticleType>::evaluateLocallyRegular ( const aol::Vec2<DataType>& coll, aol::Vec2<DataType> a, aol::Vec2<DataType> b ) const {

  aol::Vec2<typename ParticleType::DataType> dir = b - a;
  aol::Vec2<typename ParticleType::DataType> normal = dir; normal.normalize (); normal.rotateRight ();
  a -= coll; b -= coll;
  DataType len = dir.norm (), alen = a.norm (), atilde = a * dir / len, blen = b.norm (), btilde = b * dir / len, c = abs (a * normal), alpha = atan2 ( atilde, c ), beta = atan2 ( btilde, c );
  DataType result = log (blen) - log (alen);

  if ( c > 1e-24) result -= atilde / c * (beta - alpha);
  else result += alen / blen - 1;

  result /= len;

  return result * ( 0.5 / aol::NumberTrait<DataType>::pi );
}

} // namespace bm

#endif
