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

#ifndef __OSTESTFUNCTION_H
#define __OSTESTFUNCTION_H

#include <tpCFEGrid.h>
#include <iterators.h>
#include <matrixInverse.h>

namespace tpcfe {

template< typename RealType >
class InterfaceTestFunction {
protected:
  RealType     _coeffPlus, _coeffMinus; // diff coeff in _normal direction on both sides
  int          _width;                  // width of the grid

public:
  RealType     _FDh;

  InterfaceTestFunction ( ) :
    _coeffPlus ( 1.0 ), _coeffMinus ( 1.0 ), _width ( 1 ), _FDh ( 1.0e-5 ) { // 1 should lead to division by zero
  }

  virtual ~InterfaceTestFunction() {}

  virtual void setCoeffs ( const RealType coeffPlus, const RealType coeffMinus ) {
    _coeffPlus  = coeffPlus;
    _coeffMinus = coeffMinus;
  }

  virtual void setWidth ( const int width ) {
    _width = width;
  }


  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &/*p_wc*/ ) const = 0;
  virtual void     gradientInterfaceValueAt ( const aol::Vec3<RealType> &/*p_wc*/, aol::Vec3<RealType> &/*gradient*/ ) const = 0;
  virtual RealType valueAt          ( const aol::Vec3<RealType> &/*p_wc*/ ) const = 0;


  //! default implementation via difference quotient, overload on derived classes if analytic partial derivative known.
  virtual RealType dx_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
#ifdef VERBOSE
    cerr << "Computing partial x derivative via difference quotient" << endl;
#endif
    const aol::Vec3<RealType> hx ( _FDh, 0, 0 );
    return ( ( valueAt ( p_wc + hx ) - valueAt ( p_wc - hx ) ) / ( 2 * _FDh ) );
  }

  //! default implementation via difference quotient, overload on derived classes if analytic partial derivative known.
  virtual RealType dx2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
#ifdef VERBOSE
    cerr << "Computing second partial x derivative via difference quotient" << endl;
#endif
    const aol::Vec3<RealType> hx ( _FDh, 0, 0 );
    return ( ( valueAt ( p_wc + hx ) - 2 * valueAt ( p_wc ) + valueAt ( p_wc - hx ) ) / ( aol::Sqr( _FDh ) ) );
  }


  //! default implementation via difference quotient, overload on derived classes if analytic partial derivative known.
  virtual RealType dy_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
#ifdef VERBOSE
    cerr << "Computing partial y derivative via difference quotient" << endl;
#endif
    const aol::Vec3<RealType> hy ( 0, _FDh, 0 );
    return ( ( valueAt ( p_wc + hy ) - valueAt ( p_wc - hy ) ) / ( 2 * _FDh ) );
  }

  //! default implementation via difference quotient, overload on derived classes if analytic partial derivative known.
  virtual RealType dy2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
#ifdef VERBOSE
    cerr << "Computing second partial y derivative via difference quotient" << endl;
#endif
    const aol::Vec3<RealType> hy ( 0, _FDh, 0 );
    return ( ( valueAt ( p_wc + hy ) - 2 * valueAt ( p_wc ) + valueAt ( p_wc - hy ) ) / ( aol::Sqr ( _FDh ) ) );
  }


  //! default implementation via difference quotient, overload on derived classes if analytic partial derivative known.
  virtual RealType dz_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
#ifdef VERBOSE
    cerr << "Computing partial z derivative via difference quotient" << endl;
#endif
    const aol::Vec3<RealType> hz ( 0, 0, _FDh );
    return ( ( valueAt ( p_wc + hz ) - valueAt ( p_wc - hz ) ) / ( 2 * _FDh ) );
  }

  //! default implementation via difference quotient, overload on derived classes if analytic partial derivative known.
  virtual RealType dz2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
#ifdef VERBOSE
    cerr << "Computing second partial z derivative via difference quotient" << endl;
#endif
    const aol::Vec3<RealType> hz ( 0, 0, _FDh );
    return ( ( valueAt ( p_wc + hz ) - 2 * valueAt ( p_wc ) + valueAt ( p_wc - hz ) ) / ( aol::Sqr ( _FDh ) ) );
  }

  virtual RealType Laplace_valueAt ( const aol::Vec3<RealType> &p_wc )  const {
    return ( dx2_valueAt ( p_wc ) + dy2_valueAt ( p_wc ) + dz2_valueAt ( p_wc ) );
  }

  RealType interfaceValueAt_global ( const aol::Vec3<RealType> &p ) const {
    const aol::Vec3<RealType> p_wc = p / ( _width - 1.0 );
    return ( interfaceValueAt ( p_wc ) );
  }

  // the following functions allow evaluation at global coordinates by converting to world coordinates
  RealType valueAt_global          ( const aol::Vec3<RealType> &p ) const {
    const aol::Vec3<RealType> p_wc = p / ( _width - 1.0 );
    return ( valueAt ( p_wc ) );
  }

  void gradientInterfaceValueAt_global ( const aol::Vec3<RealType> &p, aol::Vec3<RealType> &gradient ) const {
    const aol::Vec3<RealType> p_wc = p / ( _width - 1.0 );
    gradientInterfaceValueAt ( p_wc, gradient );
  }


  RealType dx_valueAt_global       ( const aol::Vec3<RealType> &p ) const {
    const aol::Vec3<RealType> p_wc = p / ( _width - 1.0 );
    return ( dx_valueAt ( p_wc ) );
  }

  RealType dy_valueAt_global       ( const aol::Vec3<RealType> &p ) const {
    const aol::Vec3<RealType> p_wc = p / ( _width - 1.0 );
    return ( dy_valueAt ( p_wc ) );
  }

  RealType dz_valueAt_global       ( const aol::Vec3<RealType> &p ) const {
    const aol::Vec3<RealType> p_wc = p / ( _width - 1.0 );
    return ( dz_valueAt ( p_wc ) );
  }

  aol::Vec3<RealType> gradient_valueAt_global ( const aol::Vec3<RealType> &p ) const {
    aol::Vec3<RealType> grad;
    grad[0] = dx_valueAt_global ( p );
    grad[1] = dy_valueAt_global ( p );
    grad[2] = dz_valueAt_global ( p );
    return ( grad );
  }

  RealType Laplace_valueAt_global  ( const aol::Vec3<RealType> &p ) const {
    const aol::Vec3<RealType> p_wc = p / ( _width - 1.0 );
    return ( Laplace_valueAt ( p_wc ) );
  }

  virtual RealType coefficient_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( interfaceValueAt ( p_wc ) < 0  ?  _coeffMinus  :  _coeffPlus  );
  }

  RealType coefficientValueAtGlobal ( const aol::Vec3<RealType> &p ) const {
    const aol::Vec3<RealType> p_wc = p / ( _width - 1.0 );
    return ( coefficient_valueAt ( p_wc ) );
  }


  void createInterface ( qc::ScalarArray<RealType, qc::QC_3D> &interfc ) const {
    const RealType h = 1.0 / ( _width - 1.0 );
    if ( _width != interfc.getNumXYZ() )
      throw aol::Exception ( "tpcfe::InterfaceTestFunction<...>::createInterface: _width does not match interface", __FILE__, __LINE__ );

    for ( qc::RectangularIterator<qc::QC_3D> posit ( interfc ); posit.notAtEnd(); ++posit ) {
      const aol::Vec3<RealType> p_wc ( h * (*posit)[0], h * (*posit)[1], h * (*posit)[2] );
      interfc.set ( *posit, interfaceValueAt ( p_wc ) );
    }
  }

  void createInterfaceFunctionCoeffs ( qc::ScalarArray<RealType, qc::QC_3D> &interfc,
                                       qc::ScalarArray<RealType, qc::QC_3D> &function,
                                       qc::AArray<RealType, qc::QC_3D> &coeffs          ) const {

    if ( _width != interfc.getNumXYZ() || _width != function.getNumXYZ() || _width != coeffs.getNumXYZ() )
      throw aol::Exception ( "tpcfe::InterfaceTestFunction<...>::createInterface: _width does not match interface", __FILE__, __LINE__ );


    const RealType w1 = ( _width - 1.0 ) ;

    for ( qc::RectangularIterator<qc::QC_3D> rit ( interfc ); rit.notAtEnd(); ++rit ) {
      const aol::Vec3<RealType> p_wc ( (*rit)[0] / w1, (*rit)[1] / w1, (*rit)[2] / w1 );
      const RealType interfaceValue = interfaceValueAt ( p_wc );

      interfc.set ( *rit, interfaceValue );
      coeffs.getRef ( *rit ) = this->coefficientValueAtGlobal ( aol::Vec3<RealType>( *rit ) );
      function.set ( *rit, valueAt ( p_wc ) );
    }
  }


  void derivativeTest ( const RealType tol = 1.0e-2 ) const {
    for ( qc::RectangularIterator<qc::QC_3D> bit ( aol::Vec3<int>( 0, 0, 0 ), aol::Vec3<int> ( _width, _width, _width ) ); bit.notAtEnd(); ++bit ) {
      const RealType h = 1.0 / ( _width - 1);
      aol::Vec3<RealType> p_wc ( (*bit)[0] * h, (*bit)[1] * h, (*bit)[2] * h );

      const RealType
        dxfd  = InterfaceTestFunction<RealType>::dx_valueAt ( p_wc ),   dx = dx_valueAt ( p_wc ),
        dyfd  = InterfaceTestFunction<RealType>::dy_valueAt ( p_wc ),   dy = dy_valueAt ( p_wc ),
        dzfd  = InterfaceTestFunction<RealType>::dz_valueAt ( p_wc ),   dz = dz_valueAt ( p_wc ),

        dx2fd = InterfaceTestFunction<RealType>::dx2_valueAt ( p_wc ), dx2 = dx2_valueAt ( p_wc ),
        dy2fd = InterfaceTestFunction<RealType>::dy2_valueAt ( p_wc ), dy2 = dy2_valueAt ( p_wc ),
        dz2fd = InterfaceTestFunction<RealType>::dz2_valueAt ( p_wc ), dz2 = dz2_valueAt ( p_wc );

      if ( aol::isNaN ( dx ) || aol::isNaN ( dxfd ) || ( fabs ( dxfd - dx ) > tol * aol::Max ( fabs ( dxfd ), fabs ( dx ) ) ) )
        cerr << "x derivative test failed: " << p_wc << " " << aol::detailedFormat ( dx ) << " " << aol::detailedFormat ( dxfd ) << " "<< aol::detailedFormat ( dx / dxfd ) << endl;
      if ( aol::isNaN ( dy ) || aol::isNaN ( dyfd ) || ( fabs ( dyfd - dy ) > tol * aol::Max ( fabs ( dyfd ), fabs ( dy ) ) ) )
        cerr << "y derivative test failed: " << p_wc << " " << aol::detailedFormat ( dy ) << " " << aol::detailedFormat ( dyfd ) << " "<< aol::detailedFormat ( dy / dyfd ) << endl;
      if ( aol::isNaN ( dz ) || aol::isNaN ( dzfd ) || ( fabs ( dzfd - dz ) > tol * aol::Max ( fabs ( dzfd ), fabs ( dz ) ) ) )
        cerr << "z derivative test failed: " << p_wc << " " << aol::detailedFormat ( dz ) << " " << aol::detailedFormat ( dzfd ) << " "<< aol::detailedFormat ( dz / dzfd ) << endl;
      if ( aol::isNaN ( dx2 ) || aol::isNaN ( dx2fd ) || fabs ( dx2fd - dx2 ) > tol * aol::Max ( fabs ( dx2fd ), fabs ( dx2 ) ) )
        cerr << "xx derivative test failed: " << p_wc << " " << aol::detailedFormat ( dx2 ) << " " << aol::detailedFormat ( dx2fd ) << " " << aol::detailedFormat ( dx2 / dx2fd ) << endl;
      if ( aol::isNaN ( dy2 ) || aol::isNaN ( dy2fd ) || fabs ( dy2fd - dy2 ) > tol * aol::Max ( fabs ( dy2fd ), fabs ( dy2 ) ) )
        cerr << "yy derivative test failed: " << p_wc << " " << aol::detailedFormat ( dy2 ) << " " << aol::detailedFormat ( dy2fd ) << " " << aol::detailedFormat ( dy2 / dy2fd ) << endl;
      if ( aol::isNaN ( dz2 ) || aol::isNaN ( dz2fd ) || fabs ( dz2fd - dz2 ) > tol * aol::Max ( fabs ( dz2fd ), fabs ( dz2 ) ) )
        cerr << "zz derivative test failed: " << p_wc << " " << aol::detailedFormat ( dz2 ) << " " << aol::detailedFormat ( dz2fd ) << " " << aol::detailedFormat ( dz2 / dz2fd) << endl;
    }
  }

};


/** x-planar interface, rotated about y axis */
template< typename RealType >
class PlanarInterfaceTestFunction : public InterfaceTestFunction<RealType> {
protected:
  RealType                 _interfacePosition;  // position of interface in world coordinates
  aol::Matrix33<RealType>  _rotationMatrix, _inverseRotationMatrix;
  RealType                 _tang1, _tang2;          // tangential gradients

public:
  PlanarInterfaceTestFunction ( ) : InterfaceTestFunction<RealType>(), _interfacePosition ( -1.0 ), _tang1 ( 0.0 ), _tang2 ( 0.0 ) {
    _rotationMatrix.setIdentity();
    _inverseRotationMatrix.setIdentity();
  }

  virtual ~PlanarInterfaceTestFunction ( ) { }

public:
  void setInterfacePosition ( const RealType pos ) {
    _interfacePosition = pos;
  }

  void setRotationAngle ( const RealType theta ) { // in degrees
    _rotationMatrix.setRotationAboutY (        - theta * aol::NumberTrait<RealType>::pi / 180 );
    _inverseRotationMatrix.setRotationAboutY ( + theta * aol::NumberTrait<RealType>::pi / 180 );

  }

  void setTangCoeffs ( const RealType tang1, const RealType tang2 ) {
    this->_tang1 = tang1;
    this->_tang2 = tang2;
  }


protected:
  void rotate ( const aol::Vec3<RealType> &arg, aol::Vec3<RealType> &dest ) const {
    aol::Vec3<RealType> tmp;
    tmp = arg - aol::Vec3<RealType> ( 0.5, 0.5, 0.5 );
    dest = _rotationMatrix * tmp ;
    dest += aol::Vec3<RealType> ( 0.5, 0.5, 0.5 );
  }

  void backRotate ( const aol::Vec3<RealType> &arg, aol::Vec3<RealType> &dest ) const {
    aol::Vec3<RealType> tmp;
    tmp = arg - aol::Vec3<RealType> ( 0.5, 0.5, 0.5 );
    dest = _inverseRotationMatrix() * tmp ;
    dest += aol::Vec3<RealType> ( 0.5, 0.5, 0.5 );
  }


  void getGradient ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    const aol::Vec3<RealType> unrotatedGradient ( ( interfaceValueAt ( p_wc ) < 0  ?  this->_coeffPlus  :  this->_coeffMinus ), this->_tang1, this->_tang2 );
    gradient = _inverseRotationMatrix * unrotatedGradient;
  }

public:
  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    aol::Vec3<RealType> p_rot;
    rotate ( p_wc, p_rot );
    return ( p_rot[ 0 ] - _interfacePosition );
  }

  virtual void gradientInterfaceValueAt ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    getGradient ( p_wc, gradient );
  }


  virtual RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType dminus = 0;
    const RealType dplus  = ( this->_coeffPlus - this->_coeffMinus ) * _interfacePosition + dminus;

    aol::Vec3<RealType> p_rot;
    rotate ( p_wc, p_rot );

    const RealType common = this->_tang1 * p_rot[1] + this->_tang2 * p_rot[2];
    if ( interfaceValueAt ( p_wc ) < 0 )
      return ( common + this->_coeffPlus  * p_rot[ 0 ] + dminus );  // _coeffPlus is the slope on the minus side!
    else
      return ( common + this->_coeffMinus * p_rot[ 0 ] + dplus  );  // and vice versa
  }

  virtual RealType dx_valueAt       ( const aol::Vec3<RealType> &p_wc ) const {
    aol::Vec3<RealType> gradient;
    getGradient ( p_wc, gradient );
    return ( gradient[0] );
  }

  virtual RealType dy_valueAt       ( const aol::Vec3<RealType> &p_wc ) const {
    aol::Vec3<RealType> gradient;
    getGradient ( p_wc, gradient );
    return ( gradient[1] );
  }

  virtual RealType dz_valueAt       ( const aol::Vec3<RealType> &p_wc ) const {
    aol::Vec3<RealType> gradient;
    getGradient ( p_wc, gradient );
    return ( gradient[2] );
  }

  virtual RealType dx2_valueAt ( const aol::Vec3<RealType> & ) const {
    return ( aol::NumberTrait<RealType>::zero );
  }

  virtual RealType dy2_valueAt ( const aol::Vec3<RealType> & ) const {
    return ( aol::NumberTrait<RealType>::zero );
  }

  virtual RealType dz2_valueAt ( const aol::Vec3<RealType> & ) const {
    return ( aol::NumberTrait<RealType>::zero );
  }

};


template< typename RealType >
class CylindricalInterfaceTestFunction : public InterfaceTestFunction<RealType> {
protected:
  RealType _radius;  // position of interface in world coordinates
  RealType _tangCoeff;
  aol::Vec2<RealType> _centerWC;

public:
  CylindricalInterfaceTestFunction ( ) : InterfaceTestFunction<RealType>(), _radius ( -1.0 ), _tangCoeff ( 0.0 ) {}

  virtual ~CylindricalInterfaceTestFunction() {}

  void setRadius ( const RealType rad ) {
    _radius = rad;
  }

  void setTangCoeff ( const RealType tangCoeff ) {
    _tangCoeff = tangCoeff;
  }

  void setCenterWC ( const aol::Vec2<RealType> &ctr_wc ) {
    _centerWC = ctr_wc;
  }


  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta;
    getValues ( p_wc, radius, beta );
    return ( radius - _radius );
  }

  virtual void gradientInterfaceValueAt ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    gradient = aol::Vec3<RealType> ( p_wc[0] - _centerWC[0], p_wc[1] - _centerWC[1], 0 );
  }

  virtual RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta;
    getValues ( p_wc, radius, beta );

    RealType value;

    if ( radius < _radius ) {
      value =  - exp ( 1 - aol::Sqr ( radius / _radius ) ) + _tangCoeff * p_wc[2];
    } else {
      value = beta * ( radius - _radius ) - 1.0 + _tangCoeff * p_wc[2];
    }

    /*     cout << aol::detailedFormat ( radius ) << " " << aol::detailedFormat ( value ) << endl; */

    return ( value );
  }

  virtual RealType dx_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta;
    getValues ( p_wc, radius, beta );
    const aol::Vec2<RealType> pmc_wc ( p_wc[0] - _centerWC[0], p_wc[1] - _centerWC[1] );

    if ( radius < _radius ) {
      return ( 2.0 * pmc_wc[0] *  exp ( 1 - aol::Sqr ( radius / _radius ) )  /  aol::Sqr ( _radius ) );
    } else {
      return ( beta * pmc_wc[0] / radius );
    }
  }

  virtual RealType dx2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta;
    getValues ( p_wc, radius, beta );
    const aol::Vec2<RealType> pmc_wc ( p_wc[0] - _centerWC[0], p_wc[1] - _centerWC[1] );

    if ( radius < _radius ) {
      return ( - 4.0 * aol::Sqr ( pmc_wc[0] ) *  exp ( 1 - aol::Sqr ( radius / _radius ) )  /  aol::Sqr ( aol::Sqr ( _radius ) )  +  2 * exp ( 1 - aol::Sqr ( radius / _radius ) )  /  aol::Sqr ( _radius ) );
    } else {
      return ( - beta * aol::Sqr ( pmc_wc[0] ) / aol::Cub ( radius ) + beta / radius );
    }

  }


  virtual RealType dy_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta;
    getValues ( p_wc, radius, beta );
    const aol::Vec2<RealType> pmc_wc ( p_wc[0] - _centerWC[0], p_wc[1] - _centerWC[1] );

    if ( radius < _radius ) {
      return ( 2.0 * pmc_wc[1] *  exp ( 1 - aol::Sqr ( radius / _radius ) )  /  aol::Sqr ( _radius ) );
    } else {
      return ( beta * pmc_wc[1] / radius );
    }
  }

  virtual RealType dy2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta;
    getValues ( p_wc, radius, beta );
    const aol::Vec2<RealType> pmc_wc ( p_wc[0] - _centerWC[0], p_wc[1] - _centerWC[1] );

    if ( radius < _radius ) {
      return ( - 4.0 * aol::Sqr ( pmc_wc[1] ) *  exp ( 1 - aol::Sqr ( radius / _radius ) )  /  aol::Sqr ( aol::Sqr ( _radius ) )  +  2 * exp ( 1 - aol::Sqr ( radius / _radius ) )  /  aol::Sqr ( _radius ) );
    } else {
      return ( - beta * aol::Sqr ( pmc_wc[1] ) / aol::Cub ( radius ) + beta / radius );
    }

  }


  virtual RealType dz_valueAt ( const aol::Vec3<RealType> &/*p_wc*/ ) const {
    return ( _tangCoeff );
  }

  virtual RealType dz2_valueAt ( const aol::Vec3<RealType> &/*p_wc*/ ) const {
    return ( aol::NumberTrait<RealType>::zero );
  }



private:
  void getValues ( const aol::Vec3<RealType> p_wc, RealType &radius, RealType &beta ) const {
    radius = sqrt ( aol::Sqr ( p_wc[0] - _centerWC[0] ) + aol::Sqr ( p_wc[1] - _centerWC[1] ) );
    beta = ( 2.0 / _radius ) * ( this->_coeffMinus / this->_coeffPlus );
  }

};



template< typename RealType >
class SphericalInterfaceTestFunction : public InterfaceTestFunction<RealType> {
protected:
  RealType               _radius;    // position of interface in world coordinates
  aol::Vec3<RealType>    _centerWC;  // center of the sphere in world coordinates
  RealType               _tangFactor;

public:
  SphericalInterfaceTestFunction ( ) : InterfaceTestFunction<RealType>(), _radius ( -1.0 ), _tangFactor ( 0.0 ) {}

  virtual ~SphericalInterfaceTestFunction() {}

public:
  void setRadius ( const RealType rad ) {
    _radius = rad;
  }

  void setCenterWC ( const aol::Vec3<RealType> &ctr_wc ) {
    _centerWC = ctr_wc;
  }

  void setTangFactor ( const RealType tangFactor ) {
    _tangFactor = tangFactor;
  }

public:
  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType radius = ( p_wc - _centerWC ).norm();
    return ( radius - _radius );
  }

  virtual void gradientInterfaceValueAt ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    gradient = ( p_wc - _centerWC ); // need not be normalized!
  }

  virtual RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta, r, s;
    getValues ( p_wc, radius, beta, r, s );

    /*     cout << aol::detailedFormat ( radius ) << " " << aol::detailedFormat ( r*s ) << endl; */

    return ( r * s );
  }

  virtual RealType dx_valueAt       ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta, r, s;
    getValues ( p_wc, radius, beta, r, s );

    RealType xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues ( p_wc, xmc, ymc, zmc, rxy, ryz, rzx );

    const RealType
      partial_x_r = ( ( radius < _radius )
                      ?
                      ( ( 2.0 * xmc / aol::Sqr ( _radius ) ) *  exp ( 1 - aol::Sqr ( radius / _radius ) ) )
                      :
                      ( beta * xmc / radius ) ),

      partial_x_s = (_tangFactor != 0.0 ) ? _tangFactor * ymc * zmc * ( -             1.0  / ( rxy              * ryz * rzx              )
                                                + aol::Sqr ( xmc ) / ( aol::Cub ( rxy ) * ryz * rzx              )
                                                + aol::Sqr ( xmc ) / ( rxy              * ryz * aol::Cub ( rzx ) ) ) : 0.0;

    return ( partial_x_r * s + partial_x_s * r );
  }

  virtual RealType dy_valueAt       ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta, r, s;
    getValues ( p_wc, radius, beta, r, s );

    RealType xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues ( p_wc, xmc, ymc, zmc, rxy, ryz, rzx );

    const RealType
      partial_y_r = ( ( radius < _radius )
                      ?
                      ( ( 2.0 * ymc / aol::Sqr ( _radius ) ) *  exp ( 1 - aol::Sqr ( radius / _radius ) ) )
                      :
                      ( beta * ymc / radius ) ),

      partial_y_s = (_tangFactor != 0.0) ? _tangFactor * xmc * zmc * ( -              1.0 / ( rxy              * ryz              * rzx )
                                                + aol::Sqr ( ymc ) / ( aol::Cub ( rxy ) * rzx              * ryz )
                                                + aol::Sqr ( ymc ) / ( rxy              * aol::Cub ( ryz ) * rzx ) ) : 0.0;

    return ( partial_y_r * s + partial_y_s * r );
  }

  virtual RealType dz_valueAt       ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta, r, s;
    getValues ( p_wc, radius, beta, r, s );

    RealType xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues ( p_wc, xmc, ymc, zmc, rxy, ryz, rzx );

    const RealType
      partial_z_r = ( ( radius < _radius )
                      ?
                      ( ( 2.0 * zmc / aol::Sqr ( _radius ) ) *  exp ( 1 - aol::Sqr ( radius / _radius ) ) )
                      :
                      ( beta * zmc / radius ) ),

      partial_z_s = (_tangFactor != 0.0) ? _tangFactor * xmc * ymc * ( -              1.0 / ( rxy * ryz              * rzx              )
                                                + aol::Sqr ( zmc ) / ( rxy * ryz              * aol::Cub ( rzx ) )
                                                + aol::Sqr ( zmc ) / ( rxy * aol::Cub ( ryz ) * rzx              ) ) : 0.0;

    return ( partial_z_r * s + partial_z_s * r );
  }


private:
  inline void getValues ( const aol::Vec3<RealType> &p_wc, RealType &radius, RealType &beta, RealType &r, RealType &s ) const {
    radius = ( p_wc - _centerWC ).norm();
    beta = this->_coeffMinus / this->_coeffPlus * 2.0 / _radius;

    r = ( ( radius < _radius )  ?  ( - exp ( 1 - aol::Sqr ( radius / _radius ) ) )  :  ( beta * ( radius - _radius ) - 1.0 ) );

    RealType x, y, z, rxy, ryz, rzx;
    getSValues ( p_wc, x, y, z, rxy, ryz, rzx );
    const RealType
      sigmaxy = ( rxy > 1e-12 ? y / rxy : 0 ),
      sigmayz = ( ryz > 1e-12 ? z / ryz : 0 ),
      sigmazx = ( rzx > 1e-12 ? x / rzx : 0 );

    s = 1 - _tangFactor * sigmaxy * sigmayz * sigmazx;
  }

  inline void getSValues ( const aol::Vec3<RealType> &p_wc, RealType &x, RealType &y, RealType &z, RealType &rxy, RealType &ryz, RealType &rzx ) const {
    x = p_wc[0] - _centerWC[0];
    y = p_wc[1] - _centerWC[1];
    z = p_wc[2] - _centerWC[2];
    rxy = sqrt ( x*x + y*y );
    ryz = sqrt ( y*y + z*z );
    rzx = sqrt ( z*z + x*x );
  }

};


template< typename RealType >
class SphericalPlateauInterfaceTestFunction : public InterfaceTestFunction<RealType> {
protected:
  const aol::Vec3<RealType> _centerWC;  // center of the sphere in world coordinates
  const RealType            _expOuterRadius, _linOuterRadius, _transOuterRadius;

public:
  SphericalPlateauInterfaceTestFunction ( ) : InterfaceTestFunction<RealType>(), _centerWC( 0.5, 0.5, 0.5 ), _expOuterRadius ( 0.20 ), _linOuterRadius ( 0.30 ), _transOuterRadius ( 0.45 ) {
  }

  virtual ~SphericalPlateauInterfaceTestFunction() {
  }

  void setCenterWC ( const aol::Vec3<RealType> &ctr_wc ) {
    _centerWC = ctr_wc;
  }

  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType radius = ( p_wc - _centerWC ).norm();
    return ( _expOuterRadius - radius ); // points inside ball are "outside" ^= ball has fixed slope, "inside" has flatter slope as kappa gets large => ball has relatively steep slope
  }

  virtual void gradientInterfaceValueAt ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    gradient = ( _centerWC - p_wc );
  }

  virtual RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType radius = ( p_wc - _centerWC ).norm(), kappaInv = 1. / ( this->_coeffMinus / this->_coeffPlus ), beta = kappaInv * 2 / _expOuterRadius;
    RealType value;

    const RealType
      Re = _expOuterRadius, Rl = _linOuterRadius, Rt = _transOuterRadius,
      Ct = - exp(1.0),
      Cl = kappaInv * exp ( 1 - aol::Sqr ( ( Rt - Rl ) / ( Rt - Rl ) ) ) + kappaInv * Ct - beta * ( Rl - Re ),
      Ce = beta * ( Re - Re ) + Cl - ( - exp ( 1 - aol::Sqr ( Re / Re ) ) );

    if ( radius < Re )
      value = - exp ( 1 - aol::Sqr ( radius / Re ) ) + Ce;
    else if ( radius < Rl )
      value = beta * ( radius - Re ) + Cl;
    else if ( radius < Rt )
      value = kappaInv * exp ( 1 - aol::Sqr ( ( Rt - radius ) / ( Rt - Rl ) ) ) + kappaInv * Ct;
    else
      value = 0;


    /*     cout << aol::detailedFormat ( radius ) << " " << aol::detailedFormat ( value ) << endl;  */
    return ( value );
  }


  virtual RealType Laplace_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType radius = ( p_wc - _centerWC ).norm(), kappa = this->_coeffMinus / this->_coeffPlus;

    RealType value;

    const RealType
      Re = _expOuterRadius, Rl = _linOuterRadius, Rt = _transOuterRadius;

    if ( radius < Re ) {
      value = ( 6 / aol::Sqr ( Re ) - ( 4 * aol::Sqr ( radius ) / aol::Sqr ( aol::Sqr ( Re ) ) ) ) * exp ( 1 - aol::Sqr ( radius / Re ) ); // falsch?
    } else if ( radius < Rl ) {
      value = 4.0 / ( radius * kappa * Re );
    } else if ( radius < Rt ) {
      value = ( ( 4 * Rt - 6 * radius ) / ( radius * kappa * aol::Sqr ( Rt - Rl ) ) + ( 4 * aol::Sqr ( Rt - radius ) ) / ( kappa * aol::Sqr ( aol::Sqr ( Rt - Rl ) ) ) ) * exp ( 1 - aol::Sqr ( ( Rt - radius ) / ( Rt - Rl ) ) );
    } else {
      value = 0;
    }
    /*     cout << aol::detailedFormat ( radius ) << " " << aol::detailedFormat ( value ) << endl; */

    return ( value );
  }

};


template< typename RealType >
class CylindricalPolynomialTestFunction : public InterfaceTestFunction < RealType > {
protected:
  RealType _radius;
  aol::Vec2<RealType> _centerWC;

public:
  CylindricalPolynomialTestFunction ( ) : InterfaceTestFunction<RealType>(), _radius ( 1.0/3.0 ), _centerWC ( 0.5, 0.5 ) {
  }

  virtual ~CylindricalPolynomialTestFunction ( ) {
  }

  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( getRadius ( p_wc ) - _radius );
  }

  virtual void gradientInterfaceValueAt ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    gradient[0] = p_wc[0] - _centerWC[0];
    gradient[1] = p_wc[1] - _centerWC[1];
    gradient[2] = 0;
  }

  virtual RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType radius = getRadius ( p_wc );
    const RealType bPlus  = ( this->_coeffPlus - this->_coeffMinus ) * aol::Sqr ( _radius );
    const RealType value  = ( ( radius < _radius ) ? ( this->_coeffPlus * aol::Sqr ( radius ) ) : ( this->_coeffMinus * aol::Sqr ( radius ) + bPlus ) );

    //    cout << "VAL " << aol::detailedFormat ( radius ) << " " << aol::detailedFormat ( value ) << endl;

    return ( value );
  }


  virtual RealType dx_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( ( ( getRadius ( p_wc ) < _radius ) ? ( this->_coeffPlus ) : ( this->_coeffMinus ) ) * 2 * ( p_wc[0] - _centerWC[0] ) );
  }

  virtual RealType dy_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( ( ( getRadius ( p_wc ) < _radius ) ? ( this->_coeffPlus ) : ( this->_coeffMinus ) ) * 2 * ( p_wc[1] - _centerWC[1] ) );
  }

  virtual RealType dz_valueAt ( const aol::Vec3<RealType> & ) const {
    return ( 0 );
  }

  virtual RealType dx2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( 2 * ( ( getRadius ( p_wc ) < _radius ) ? ( this->_coeffPlus ) : ( this->_coeffMinus ) ) );
  }

  virtual RealType dy2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( 2 * ( ( getRadius ( p_wc ) < _radius ) ? ( this->_coeffPlus ) : ( this->_coeffMinus ) ) );
  }

  virtual RealType dz2_valueAt ( const aol::Vec3<RealType> & ) const {
    return ( 0 );
  }

private:
  inline RealType getRadius ( const aol::Vec3<RealType> &p_wc ) const {
    return ( ( aol::Vec2<RealType> ( p_wc[0], p_wc[1] ) - _centerWC ).norm() );
  }

};


template< typename RealType >
class SphericalPolynomialTestFunction : public InterfaceTestFunction < RealType > {
protected:
  RealType _radius, _tangFactor;
  aol::Vec3<RealType> _centerWC;

public:
  SphericalPolynomialTestFunction ( ) : InterfaceTestFunction<RealType>(), _radius ( 1.0/3.0 ), _tangFactor ( 0 ), _centerWC ( 0.5, 0.5, 0.5 ) {
  }

  virtual ~SphericalPolynomialTestFunction ( ) {
  }

  void setTangFactor ( RealType TangFactor ) {
    _tangFactor = TangFactor;
  }

  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( getRadius ( p_wc ) - _radius );
  }

  virtual void gradientInterfaceValueAt ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    gradient = ( p_wc - _centerWC );
  }

  virtual RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    RealType x, y, z, rxy, ryz, rxz, s;
    getSValues ( p_wc, x, y, z, rxy, ryz, rxz, s );

    // cout << "VAL " << aol::detailedFormat ( getRadius ( p_wc ) ) << " " << aol::detailedFormat ( s * getRValue( p_wc ) ) << endl;

    return ( s * getRValue( p_wc ) );
  }


  virtual RealType dx_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType s, xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues( p_wc, xmc, ymc, zmc, rxy, ryz, rzx, s );

    const RealType
      r   = getRValue ( p_wc ),
      dxr = ( ( getRadius ( p_wc ) < _radius ) ? ( this->_coeffPlus ) : ( this->_coeffMinus ) ) * 2 * ( p_wc - _centerWC )[ 0 ],
      dxs = (_tangFactor == 0.0) ? 0.0 : _tangFactor * ymc * zmc * ( -             1.0  / ( rxy              * ryz * rzx              )
                                        + aol::Sqr ( xmc ) / ( aol::Cub ( rxy ) * ryz * rzx              )
                                        + aol::Sqr ( xmc ) / ( rxy              * ryz * aol::Cub ( rzx ) ) );
    return ( dxr * s + r * dxs );
  }

  virtual RealType dy_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType s, xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues( p_wc, xmc, ymc, zmc, rxy, ryz, rzx, s );

    const RealType
      r   = getRValue ( p_wc ),
      dyr = ( ( getRadius ( p_wc ) < _radius ) ? ( this->_coeffPlus ) : ( this->_coeffMinus ) ) * 2 * ( p_wc - _centerWC )[ 1 ],
      dys = (_tangFactor == 0.0) ? 0.0 : _tangFactor * xmc * zmc * ( -              1.0 / ( rxy              * ryz              * rzx )
                                        + aol::Sqr ( ymc ) / ( aol::Cub ( rxy ) * rzx              * ryz )
                                        + aol::Sqr ( ymc ) / ( rxy              * aol::Cub ( ryz ) * rzx ) );
    return ( dyr * s + r * dys );
  }

  virtual RealType dz_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType s, xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues( p_wc, xmc, ymc, zmc, rxy, ryz, rzx, s );

    const RealType
      r   = getRValue ( p_wc ),
      dzr = ( ( getRadius ( p_wc ) < _radius ) ? ( this->_coeffPlus ) : ( this->_coeffMinus ) ) * 2 * ( p_wc - _centerWC )[ 2 ],
      dzs = (_tangFactor == 0.0) ? 0.0 : _tangFactor * xmc * ymc * ( -              1.0 / ( rxy * ryz              * rzx              )
                                        + aol::Sqr ( zmc ) / ( rxy * ryz              * aol::Cub ( rzx ) )
                                        + aol::Sqr ( zmc ) / ( rxy * aol::Cub ( ryz ) * rzx              ) );
    return ( dzr * s + r * dzs );
  }

private:
  inline RealType getRadius ( const aol::Vec3<RealType> &p_wc ) const {
    return ( ( p_wc - _centerWC ).norm() );
  }

  inline RealType getRValue ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType radius = getRadius ( p_wc );
    const RealType bPlus  = ( this->_coeffPlus - this->_coeffMinus ) * aol::Sqr ( _radius );
    return ( ( radius < _radius ) ? ( this->_coeffPlus * aol::Sqr ( radius ) ) : ( this->_coeffMinus * aol::Sqr ( radius ) + bPlus ) );
  }

  inline void getSValues ( const aol::Vec3<RealType> &p_wc, RealType &x, RealType &y, RealType &z, RealType &rxy, RealType &ryz, RealType &rzx, RealType &s ) const {
    x = p_wc[0] - _centerWC[0];
    y = p_wc[1] - _centerWC[1];
    z = p_wc[2] - _centerWC[2];
    rxy = sqrt ( x*x + y*y );
    ryz = sqrt ( y*y + z*z );
    rzx = sqrt ( z*z + x*x );

    const RealType
      sigmaxy = ( rxy > 1e-12 ? y / rxy : 0 ),
      sigmayz = ( ryz > 1e-12 ? z / ryz : 0 ),
      sigmazx = ( rzx > 1e-12 ? x / rzx : 0 );

    s = 1 - _tangFactor * sigmaxy * sigmayz * sigmazx;
  }

};


template< typename RealType >
class SphericalPolynomialDoubleKinkTestFunction : public InterfaceTestFunction < RealType > {
protected:
  RealType _radius0, _radius1, _tangFactor, _secondKinkFactor;
  aol::Vec3<RealType> _centerWC;

public:
  SphericalPolynomialDoubleKinkTestFunction ( ) : InterfaceTestFunction<RealType>(), _radius0 ( 2.0/9.0 ),  _radius1 ( 4.0/9.0 ), _tangFactor ( 0.0 ), _secondKinkFactor ( sqrt ( 2.0 ) ), _centerWC ( 0.5, 0.5, 0.5 ) {
  }

  virtual ~SphericalPolynomialDoubleKinkTestFunction ( ) {
  }

  void setRadii( const RealType radius0, const RealType radius1){
    _radius0 = radius0;
    _radius1 = radius1;
  }

  void setTangFactor ( RealType TangFactor ) {
    _tangFactor = TangFactor;
  }

  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType rad = getRadius ( p_wc );
    return ( aol::Min ( rad - _radius0, _radius1 - rad ) );
  }

  virtual RealType coefficient_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType rad = getRadius ( p_wc );
    if ( rad < _radius0 ) {
      return ( this->_coeffMinus );
    } else if ( rad < _radius1 ) {
      return ( this->_coeffPlus );
    } else{
      return ( 1.0 / _secondKinkFactor * this->_coeffMinus );
    }
  }

  virtual void gradientInterfaceValueAt ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    const RealType rad = getRadius ( p_wc );
    if ( rad < ( 0.5 * ( _radius0 + _radius1 ) ) ) {
      gradient = p_wc - _centerWC;
    } else {
      gradient = _centerWC - p_wc;
    }
  }

  virtual RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    RealType x, y, z, rxy, ryz, rxz, s;
    getSValues ( p_wc, x, y, z, rxy, ryz, rxz, s );
    // cout << "VAL " << aol::detailedFormat ( getRadius ( p_wc ) ) << " " << aol::detailedFormat ( getRValue( p_wc ) ) << " " << aol::detailedFormat ( s ) << " " << aol::detailedFormat ( s * getRValue( p_wc ) ) << endl;
    return ( s * getRValue( p_wc ) );
  }


  virtual RealType dx_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType s, xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues( p_wc, xmc, ymc, zmc, rxy, ryz, rzx, s );

    const RealType
      r   = getRValue ( p_wc ),
      dxr = getRPartialDerivative ( p_wc, qc::QC_X ),
      trm = _tangFactor * ymc * zmc * ( -             1.0  / ( rxy              * ryz * rzx              )
                                        + aol::Sqr ( xmc ) / ( aol::Cub ( rxy ) * ryz * rzx              )
                                        + aol::Sqr ( xmc ) / ( rxy              * ryz * aol::Cub ( rzx ) ) ),
      dxs = ( aol::isNaN ( trm ) ? 0.0 : trm );
    return ( dxr * s + r * dxs );
  }

  virtual RealType dy_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType s, xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues( p_wc, xmc, ymc, zmc, rxy, ryz, rzx, s );

    const RealType
      r   = getRValue ( p_wc ),
      dyr = getRPartialDerivative ( p_wc, qc::QC_Y ),
      trm = _tangFactor * xmc * zmc * ( -              1.0 / ( rxy              * ryz              * rzx )
                                        + aol::Sqr ( ymc ) / ( aol::Cub ( rxy ) * rzx              * ryz )
                                        + aol::Sqr ( ymc ) / ( rxy              * aol::Cub ( ryz ) * rzx ) ),
      dys = ( aol::isNaN ( trm ) ? 0.0 : trm );
    return ( dyr * s + r * dys );
  }

  virtual RealType dz_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType s, xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues( p_wc, xmc, ymc, zmc, rxy, ryz, rzx, s );

    const RealType
      r   = getRValue ( p_wc ),
      dzr = getRPartialDerivative ( p_wc, qc::QC_Z ),
      trm = _tangFactor * xmc * ymc * ( -              1.0 / ( rxy * ryz              * rzx              )
                                        + aol::Sqr ( zmc ) / ( rxy * ryz              * aol::Cub ( rzx ) )
                                        + aol::Sqr ( zmc ) / ( rxy * aol::Cub ( ryz ) * rzx              ) ),
      dzs = ( aol::isNaN ( trm ) ? 0.0 : trm );
    return ( dzr * s + r * dzs );
  }

  virtual RealType dx2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType s, xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues( p_wc, xmc, ymc, zmc, rxy, ryz, rzx, s );
    const RealType
      r    = getRValue ( p_wc ),
      dxr  = getRPartialDerivative ( p_wc, qc::QC_X ),
      dx2r = getRSecondPartialDerivative ( p_wc, qc::QC_X ),
      trm  = _tangFactor * ymc * zmc * ( -             1.0  / ( rxy              * ryz * rzx              )
                                         + aol::Sqr ( xmc ) / ( aol::Cub ( rxy ) * ryz * rzx              )
                                         + aol::Sqr ( xmc ) / ( rxy              * ryz * aol::Cub ( rzx ) ) ),
      dxs = ( aol::isNaN ( trm ) ? 0.0 : trm ),
      term0 = _tangFactor * ymc * zmc,
      term1 = ryz / ( aol::Sqr ( rxy * rzx * ryz ) ) * ( xmc / rxy * rzx + xmc / rzx * rxy ),
      term2 = ( ( 2 * xmc * aol::Cub ( rxy ) * ryz * rzx ) - aol::Sqr ( xmc ) * ryz * ( 3 * xmc * rxy * rzx + aol::Cub ( rxy ) * xmc / rzx ) ) / aol::Sqr ( aol::Cub ( rxy ) * ryz * rzx ),
      term3 = ( ( 2 * xmc * rxy * ryz * aol::Cub ( rzx ) ) - aol::Sqr ( xmc ) * ryz * ( xmc / rxy * aol::Cub ( rzx ) + rxy * 3 * xmc * rzx ) ) / aol::Sqr ( rxy * ryz * aol::Cub ( rzx ) ),
      dx2s = term0 * ( ( aol::isNaN ( term1 ) ? 0.0 : term1 ) + ( aol::isNaN ( term2 ) ? 0.0 : term2 ) + ( aol::isNaN ( term3 ) ? 0.0 : term3 ) ),
      ret = r * dx2s + 2 * dxr * dxs + dx2r * s;

    return ( ret );
  }

  virtual RealType dy2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType s, xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues( p_wc, xmc, ymc, zmc, rxy, ryz, rzx, s );
    const RealType
      r    = getRValue ( p_wc ),
      dyr  = getRPartialDerivative ( p_wc, qc::QC_Y ),
      dy2r = getRSecondPartialDerivative ( p_wc, qc::QC_Y ),
      trm  = _tangFactor * xmc * zmc * ( -              1.0 / ( rxy              * ryz              * rzx )
                                         + aol::Sqr ( ymc ) / ( aol::Cub ( rxy ) * rzx              * ryz )
                                         + aol::Sqr ( ymc ) / ( rxy              * aol::Cub ( ryz ) * rzx ) ),
      dys = ( aol::isNaN ( trm ) ? 0.0 : trm ),
      term0 = _tangFactor * xmc * zmc,
      term1 = rzx / ( aol::Sqr ( rxy * rzx * ryz ) ) * ( ymc / rxy * ryz + ymc / ryz * rxy ),
      term2 = ( ( 2 * ymc * aol::Cub ( rxy ) * ryz * rzx ) - aol::Sqr ( ymc ) * rzx * ( 3 * ymc * rxy * ryz + aol::Cub ( rxy ) * ymc / ryz ) ) / aol::Sqr ( aol::Cub ( rxy ) * ryz * rzx ),
      term3 = ( ( 2 * ymc * rxy * aol::Cub ( ryz ) * rzx ) - aol::Sqr ( ymc ) * rzx * ( ymc / rxy * aol::Cub ( ryz ) + rxy * 3 * ymc * ryz ) ) / aol::Sqr ( rxy * aol::Cub ( ryz ) * rzx ),
      dy2s = term0 * ( ( aol::isNaN ( term1 ) ? 0.0 : term1 ) + ( aol::isNaN ( term2 ) ? 0.0 : term2 ) + ( aol::isNaN ( term3 ) ? 0.0 : term3 ) ),
      ret = r * dy2s + 2 * dyr * dys + dy2r * s;

    return ( ret );
  }

  virtual RealType dz2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    RealType s, xmc, ymc, zmc, rxy, ryz, rzx;
    getSValues( p_wc, xmc, ymc, zmc, rxy, ryz, rzx, s );

    const RealType
      r    = getRValue ( p_wc ),
      dzr  = getRPartialDerivative ( p_wc, qc::QC_Z ),
      dz2r = getRSecondPartialDerivative ( p_wc, qc::QC_Z ),
      trm  = _tangFactor * xmc * ymc * ( -              1.0 / ( rxy * ryz              * rzx              )
                                         + aol::Sqr ( zmc ) / ( rxy * ryz              * aol::Cub ( rzx ) )
                                         + aol::Sqr ( zmc ) / ( rxy * aol::Cub ( ryz ) * rzx              ) ),
      dzs = ( aol::isNaN ( trm ) ? 0.0 : trm ),
      term0 = _tangFactor * xmc * ymc,
      term1 = rxy / ( aol::Sqr ( rxy * rzx * ryz ) ) * ( zmc / ryz * rzx + zmc / rzx * ryz ),
      term2 = ( ( 2 * zmc * rxy * ryz * aol::Cub ( rzx ) ) - aol::Sqr ( zmc ) * rxy * ( 3 * zmc * rzx * ryz + aol::Cub ( rzx ) * zmc / ryz ) ) / aol::Sqr ( rxy * ryz * aol::Cub ( rzx ) ),
      term3 = ( ( 2 * zmc * rxy * aol::Cub ( ryz ) * rzx ) - aol::Sqr ( zmc ) * rxy * ( zmc / rzx * aol::Cub ( ryz ) + rzx * 3 * zmc * ryz ) ) / aol::Sqr ( rxy * aol::Cub ( ryz ) * rzx ),
      dz2s = term0 * ( ( aol::isNaN ( term1 ) ? 0.0 : term1 ) + ( aol::isNaN ( term2 ) ? 0.0 : term2 ) + ( aol::isNaN ( term3 ) ? 0.0 : term3 ) ),
      ret = r * dz2s + 2 * dzr * dzs + dz2r * s;

    return ( ret );
  }

private:
  inline RealType getRadius ( const aol::Vec3<RealType> &p_wc ) const {
    return ( ( p_wc - _centerWC ).norm() );
  }

  inline RealType getRValue ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType rad = getRadius ( p_wc );
    if ( rad < _radius0 ) {
      return ( this->_coeffPlus  * aol::Sqr ( rad ) );
    } else if ( rad < _radius1 ) {
      return ( this->_coeffMinus * aol::Sqr ( rad ) + ( this->_coeffPlus - this->_coeffMinus ) * aol::Sqr ( _radius0 ) );
    } else{
      return ( _secondKinkFactor * this->_coeffPlus  * aol::Sqr ( rad ) + ( this->_coeffMinus - _secondKinkFactor * this->_coeffPlus ) * aol::Sqr ( _radius1 ) + ( this->_coeffPlus - this->_coeffMinus ) * aol::Sqr ( _radius0 ) );
    }
  }

  inline RealType getRPartialDerivative ( const aol::Vec3<RealType> &p_wc, const qc::Comp Dir ) const {
    const RealType rad = getRadius ( p_wc );
    RealType slope = 0;
    if ( rad < _radius0 ) {
      slope = this->_coeffPlus;
    } else if ( rad < _radius1 ) {
      slope = this->_coeffMinus;
    } else {
      slope = _secondKinkFactor * this->_coeffPlus;
    }
    return ( 2 * slope * ( p_wc - _centerWC )[ Dir ] );
  }

  inline RealType getRSecondPartialDerivative ( const aol::Vec3<RealType> &p_wc, const qc::Comp /*Dir*/ ) const {
    const RealType rad = getRadius ( p_wc );
    RealType slope = 0;
    if ( rad < _radius0 ) {
      slope = this->_coeffPlus;
    } else if ( rad < _radius1 ) {
      slope = this->_coeffMinus;
    } else {
      slope = _secondKinkFactor * this->_coeffPlus;
    }
    return ( 2 * slope );
  }



  inline void getSValues ( const aol::Vec3<RealType> &p_wc, RealType &x, RealType &y, RealType &z, RealType &rxy, RealType &ryz, RealType &rzx, RealType &s ) const {
    x = p_wc[0] - _centerWC[0];
    y = p_wc[1] - _centerWC[1];
    z = p_wc[2] - _centerWC[2];
    rxy = sqrt ( x*x + y*y );
    ryz = sqrt ( y*y + z*z );
    rzx = sqrt ( z*z + x*x );

    const RealType
      sigmaxy = ( rxy > 1e-12 ? y / rxy : 0 ),
      sigmayz = ( ryz > 1e-12 ? z / ryz : 0 ),
      sigmazx = ( rzx > 1e-12 ? x / rzx : 0 );

    s = 1 - _tangFactor * sigmaxy * sigmayz * sigmazx;
  }

};


template< typename RealType >
class SphericalPolynomialPlateauTestFunction : public InterfaceTestFunction < RealType > {
protected:
  RealType _kinkRadius, _transOuterRadius;
  aol::Vec3<RealType> _centerWC;
  aol::Vector<RealType> _polyC;

public:
  SphericalPolynomialPlateauTestFunction ( ) : InterfaceTestFunction<RealType>(), _kinkRadius ( 0.3 ), _transOuterRadius ( 0.45 ), _centerWC ( 0.5, 0.5, 0.5 ), _polyC ( 10 ) {
  }

  virtual void setCoeffs ( const RealType coeffPlus, const RealType coeffMinus ) {
    InterfaceTestFunction<RealType>::setCoeffs ( coeffPlus, coeffMinus );

    const RealType
      K = this->_coeffMinus / this->_coeffPlus,
      R = _kinkRadius,
      T = _transOuterRadius;

    aol::FullMatrix<RealType> PC ( 10, 10 );
    aol::Vector<RealType> rhs ( 10 );
    PC.set ( 0, 0,  1.0 );    rhs[0] = 1.0;
    PC.set ( 1, 1,  1.0 );
    PC.set ( 2, 2,  2.0 );
    PC.set ( 3, 2,  2.0 );    PC.set ( 3, 3, 6*R   );    PC.set ( 3, 4, 12*R*R   );
    PC.set ( 4, 7,  2.0 );    PC.set ( 4, 8, 6*R   );    PC.set ( 4, 9, 12*R*R   );
    PC.set ( 5, 5,  1.0 );    PC.set ( 5, 6,  T    );    PC.set ( 5, 7,  T*T     );    PC.set ( 5, 8,  T*T*T    );    PC.set ( 5, 9, T*T*T*T  );
    PC.set ( 6, 6,  1.0 );    PC.set ( 6, 7, 2*T   );    PC.set ( 6, 8, 3*T*T    );    PC.set ( 6, 9, 4*T*T*T   );
    PC.set ( 7, 7,  2.0 );    PC.set ( 7, 8, 6*T   );    PC.set ( 7, 9, 12*T*T   );
    PC.set ( 8, 0,  1.0 );    PC.set ( 8, 1,  R    );    PC.set ( 8, 2,  R*R     );    PC.set ( 8, 3,  R*R*R    );    PC.set ( 8, 4, R*R*R*R  );
    PC.set ( 8, 5, -1.0 );    PC.set ( 8, 6, -R    );    PC.set ( 8, 7, -R*R     );    PC.set ( 8, 8, -R*R*R    );    PC.set ( 8, 9, -R*R*R*R );
    PC.set ( 9, 1,   K  );    PC.set ( 9, 2, 2*K*R );    PC.set ( 9, 3,  3*K*R*R );    PC.set ( 9, 4, 4*K*R*R*R );
    PC.set ( 9, 6, -1.0 );    PC.set ( 9, 7, -2*R  );    PC.set ( 9, 8, -3*R*R   );    PC.set ( 9, 9, -4*R*R*R  );

    aol::LUInverse<RealType> PCi ( PC );
    PCi.apply ( rhs, _polyC );

    cerr << _polyC << endl;

  }

  virtual ~SphericalPolynomialPlateauTestFunction ( ) {
  }

  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( getRadius ( p_wc ) - _kinkRadius ); // Ball = inside
  }

  virtual void gradientInterfaceValueAt ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    gradient = ( p_wc - _centerWC );
  }

  virtual RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType r = getRadius ( p_wc );

    RealType value;
    if ( r < _kinkRadius ) { // inner part
      value = _polyC[0] + _polyC[1] * r + _polyC[2] * r*r + _polyC[3] * r*r*r + _polyC[4] * r*r*r*r;
    } else if ( r < _transOuterRadius ) { // outer part
      value = _polyC[5] + _polyC[6] * r + _polyC[7] * r*r + _polyC[8] * r*r*r + _polyC[9] * r*r*r*r;
    } else { // zero outside
      value = 0;
    }

    // cout << "VAL " << aol::longScientificFormat ( r ) << " " << aol::longScientificFormat ( value ) << endl;

    return ( value );
  }

  virtual RealType dx_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( getPartialDerivative ( p_wc, 0 ) );
  }

  virtual RealType dy_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( getPartialDerivative ( p_wc, 1 ) );
  }

  virtual RealType dz_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( getPartialDerivative ( p_wc, 2 ) );
  }

  virtual RealType Laplace_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType r = getRadius ( p_wc );

    RealType value;
    if ( r == 0 ) {
      value = 0; // continuous extension since _polyC[1] = 0
    } else if ( r < _kinkRadius ) { // inner part
      value = 2 * _polyC[1] / r + 6 * _polyC[2] + 12 * _polyC[3] * r + 20 * _polyC[4] * r * r;
    } else if ( r < _transOuterRadius ) { // outer part
      value = 2 * _polyC[6] / r + 6 * _polyC[7] + 12 * _polyC[8] * r + 20 * _polyC[9] * r * r;
    } else { // zero outside
      value = 0;
    }

    return ( value );

  }

private:
  inline RealType getRadius ( const aol::Vec3<RealType> &p_wc ) const {
    return ( ( p_wc - _centerWC ).norm() );
  }

  RealType getPartialDerivative ( const aol::Vec3<RealType> &p_wc, const short dir ) const {
    const RealType r = getRadius ( p_wc );

    RealType value = r * dir;
    if ( r == 0 ) {
      value = 0; // continuous extension since _polyC[1] = 0
    } else if ( r < _kinkRadius ) { // inner part
      value = ( _polyC[1] / r + 2 * _polyC[2] + 3 * _polyC[3] * r + 4 * _polyC[4] * r * r ) * ( p_wc - _centerWC )[ dir ];
    } else if ( r < _transOuterRadius ) { // outer part
      value = ( _polyC[6] / r + 2 * _polyC[7] + 3 * _polyC[8] * r + 4 * _polyC[9] * r * r ) * ( p_wc - _centerWC )[ dir ];
    } else { // zero outside
      value = 0;
    }

    return ( value );
  }

};


template< typename RealType >
class PlanarPolynomialTestFunction : public InterfaceTestFunction < RealType > {
protected:
  RealType _thres;
  aol::Matrix33<RealType>  _rotationMatrix, _inverseRotationMatrix;

public:
  PlanarPolynomialTestFunction ( ) : InterfaceTestFunction<RealType>(), _thres ( 1./3. ) {
    _rotationMatrix.setRotationAboutY (          3 * aol::NumberTrait<RealType>::pi / 180 );
    _inverseRotationMatrix.setRotationAboutY ( - 3 * aol::NumberTrait<RealType>::pi / 180 );
  }

  virtual ~PlanarPolynomialTestFunction ( ) {
  }

public:
  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    aol::Vec3<RealType> p_rot = _rotationMatrix * p_wc;
    return ( p_rot[0] - _thres );
  }

  virtual void gradientInterfaceValueAt ( const aol::Vec3<RealType> &/*p_wc*/, aol::Vec3<RealType> &gradient ) const {
    gradient = _inverseRotationMatrix * aol::Vec3<RealType> ( 1.0, 0.0, 0.0 ); // is this correct?
  }

  virtual RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    aol::Vec3<RealType> p_rot = _rotationMatrix * p_wc;
    const RealType bPlus  = ( this->_coeffMinus - this->_coeffPlus ) * aol::Sqr ( _thres );
    const RealType value  = ( ( p_rot[0] < _thres ) ? ( - this->_coeffPlus * aol::Sqr ( p_rot[0] ) ) : ( - this->_coeffMinus * aol::Sqr ( p_rot[0] ) + bPlus ) );
    return ( value );
  }
};


template< typename RealType, typename ITFType >
class PeriodicInterfaceTestFunction : public InterfaceTestFunction<RealType> {
protected:
  RealType FREQ;
  ITFType *_pCellFct;

public:
  explicit PeriodicInterfaceTestFunction ( const RealType Freq ) : InterfaceTestFunction<RealType>(), FREQ ( Freq ), _pCellFct ( NULL ) {
    cerr << "Periodicity " << FREQ << endl;
  }

  void setCellFct ( ITFType &CellFct ) {
    _pCellFct = &CellFct;
  }

  void setFreq ( const RealType freq ) {
    FREQ = freq;
    cerr << "Periodicity " << FREQ << endl;
  }

  RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( _pCellFct->interfaceValueAt ( pPos ( p_wc ) ) );
  }

  void gradientInterfaceValueAt ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    _pCellFct->gradientInterfaceValueAt ( pPos ( p_wc ), gradient );
    gradient *= FREQ;
  }

  RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    return ( _pCellFct->valueAt ( pPos ( p_wc ) ) );
  }

  RealType dx_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( FREQ * _pCellFct->dx_valueAt ( pPos ( p_wc ) ) );
  }

  RealType dx2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( FREQ * FREQ * _pCellFct->dx2_valueAt ( pPos ( p_wc ) ) );
  }


  RealType dy_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( FREQ * _pCellFct->dy_valueAt ( pPos ( p_wc ) ) );
  }

  RealType dy2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( FREQ * FREQ * _pCellFct->dy2_valueAt ( pPos ( p_wc ) ) );
  }


  RealType dz_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( FREQ * _pCellFct->dz_valueAt ( pPos ( p_wc ) ) );
  }

  RealType dz2_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( FREQ * FREQ * _pCellFct->dz2_valueAt ( pPos ( p_wc ) ) );
  }


  RealType Laplace_valueAt ( const aol::Vec3<RealType> &p_wc ) const {
    return ( FREQ * FREQ * _pCellFct->Laplace_valueAt ( pPos ( p_wc ) ) );
  }


  void setCoeffs ( const RealType coeffPlus, const RealType coeffMinus ) {
    InterfaceTestFunction<RealType>::setCoeffs ( coeffPlus, coeffMinus );
    _pCellFct->setCoeffs ( coeffPlus, coeffMinus );
  }

  void setWidth ( const int width ) {
    InterfaceTestFunction<RealType>::setWidth ( width );
    _pCellFct->setWidth ( width );
  }

protected:
  inline aol::Vec3<RealType> pPos ( const aol::Vec3<RealType> &pos ) const {
    aol::Vec3<RealType> tmp ( pos );
    // trying to avoid rounding problems
    for ( short i = 0; i < 3; ++i )
      if ( tmp[i] == 1.0 )
        tmp[i] -= 1.0e-15;

    tmp *= FREQ;

    return ( aol::Vec3<RealType> ( tmp[0] - floor ( tmp[0] ), tmp[1] - floor ( tmp[1] ), tmp[2] - floor ( tmp[2] ) ) );
  }

};

template< typename RealType >
class SwissCheeseInterfaceTestFunction : public PeriodicInterfaceTestFunction< RealType, SphericalPlateauInterfaceTestFunction<RealType> > {
protected:
  SphericalPlateauInterfaceTestFunction<RealType> _sITF;

public:
  explicit SwissCheeseInterfaceTestFunction ( const RealType Freq ) : PeriodicInterfaceTestFunction< RealType, SphericalPlateauInterfaceTestFunction<RealType> >( Freq ), _sITF() {
    this->setCellFct ( _sITF );
  }

  virtual ~SwissCheeseInterfaceTestFunction() {}

  void setRadius ( const RealType rad ) {
    _sITF.setRadius ( rad );
  }

  void setCenterWC ( const aol::Vec3<RealType> &ctr_wc ) {
    _sITF.setCenterWC ( ctr_wc );
  }

};

template< typename RealType >
  class SwissCheesePolynomialInterfaceTestFunction : public PeriodicInterfaceTestFunction< RealType, SphericalPolynomialPlateauTestFunction<RealType> > {
protected:
  SphericalPolynomialPlateauTestFunction<RealType> _sITF;

public:
  explicit SwissCheesePolynomialInterfaceTestFunction ( const RealType Freq ) : PeriodicInterfaceTestFunction< RealType, SphericalPolynomialPlateauTestFunction<RealType> >( Freq ), _sITF() {
    this->setCellFct ( _sITF );
  }

  virtual ~SwissCheesePolynomialInterfaceTestFunction() {}
};


template< typename RealType >
InterfaceTestFunction<RealType>* selectITFByType ( const int TYPE ) {

  const RealType planarPosition = 3.0 / 7.0;

  switch ( TYPE ) {

    // planar interface

  case 0: { // aligned planar interface, tangential gradients = 0

    tpcfe::PlanarInterfaceTestFunction<RealType>* testFct = new tpcfe::PlanarInterfaceTestFunction<RealType>;
    testFct->setInterfacePosition ( planarPosition );
    return ( testFct );
  }
  break;

  case 1: { // aligned planar interface, tangential gradients = 2, 5

    tpcfe::PlanarInterfaceTestFunction<RealType>* testFct = new tpcfe::PlanarInterfaceTestFunction<RealType>;
    testFct->setInterfacePosition ( planarPosition );
    testFct->setTangCoeffs ( 2.0, 5.0 );
    return ( testFct );
  }
  break;

  case 2: { // non-aligned planar interface, tangential gradients = 0

    tpcfe::PlanarInterfaceTestFunction<RealType> *testFct = new tpcfe::PlanarInterfaceTestFunction<RealType>;
    testFct->setInterfacePosition ( planarPosition );
    testFct->setRotationAngle ( 23.0 );
    return ( testFct );
  }
  break;

  case 3: { // non-aligned planar interface, tangential gradients = 2, 5

    tpcfe::PlanarInterfaceTestFunction<RealType> *testFct = new tpcfe::PlanarInterfaceTestFunction<RealType>;
    testFct->setInterfacePosition ( planarPosition );
    testFct->setRotationAngle ( 23.0 );
    testFct->setTangCoeffs ( 2.0, 5.0 );
    return ( testFct );
  }
  break;


  // cylindrical interface

  case 4: { // cylindrical interface centered outside the bounding cube, tangential gradient 0

    CylindricalInterfaceTestFunction<RealType>* testFct = new CylindricalInterfaceTestFunction<RealType>;
    testFct->setRadius ( 3.6 * sqrt ( 2.0 ) );
    testFct->setCenterWC ( aol::Vec2<RealType> ( -3.0, -3.0 ) );
    return ( testFct );
  }
  break;

  case 5: { // cylindrical interface centered outside the bounding cube, tangential gradient 1

    CylindricalInterfaceTestFunction<RealType>* testFct = new CylindricalInterfaceTestFunction<RealType>;
    testFct->setRadius ( 3.6 * sqrt ( 2.0 ) );
    testFct->setCenterWC ( aol::Vec2<RealType> ( -3.0, -3.0 ) );
    testFct->setTangCoeff ( 0.1 );
    return ( testFct );
  }
  break;

  case 6: { // cylindrical interface centered inside the bounding cube, tangential gradient 0

    CylindricalInterfaceTestFunction<RealType>* testFct = new CylindricalInterfaceTestFunction<RealType>;
    testFct->setRadius ( 1.0 / 3.0 );
    testFct->setCenterWC ( aol::Vec2<RealType> ( 0.5, 0.5 ) );
    return ( testFct );
  }
  break;

  case 7: { // cylindrical interface centered inside the bounding cube, tangential gradient 1

    CylindricalInterfaceTestFunction<RealType>* testFct = new CylindricalInterfaceTestFunction<RealType>;
    testFct->setRadius ( 1.0 / 3.0 );
    testFct->setCenterWC ( aol::Vec2<RealType> ( 0.5, 0.5 ) );
    testFct->setTangCoeff ( 0.1 );
    return ( testFct );
  }
  break;


  // spherical interface

  case 8: { // spherical interface centered outside the bounding cube, no tangential gradient

    SphericalInterfaceTestFunction<RealType>* testFct = new SphericalInterfaceTestFunction<RealType>;
    testFct->setRadius ( 3.6 * sqrt ( 3.0 ) );
    testFct->setCenterWC ( aol::Vec3<RealType> ( -3, -3, -3 ) );
    return ( testFct );
  }
  break;

  case 9: { // spherical interface centered outside the bounding cube, tangential gradient present

    SphericalInterfaceTestFunction<RealType>* testFct = new SphericalInterfaceTestFunction<RealType>;
    testFct->setRadius ( 3.6 * sqrt ( 3.0 ) );
    testFct->setCenterWC ( aol::Vec3<RealType> ( -3, -3, -3 ) );
    testFct->setTangFactor ( 0.1 );
    return ( testFct );
  }
  break;

  case 10: { // spherical interface centered inside the bounding cube, no tangential gradient

    SphericalInterfaceTestFunction<RealType>* testFct = new SphericalInterfaceTestFunction<RealType>;
    testFct->setRadius ( 1.0 / 3.0 );
    testFct->setCenterWC ( aol::Vec3<RealType> ( 0.5, 0.5, 0.5 ) );
    return ( testFct );
  }
  break;

  case 11: { // spherical interface centered inside the bounding cube, tangential gradient present

    SphericalInterfaceTestFunction<RealType>* testFct = new SphericalInterfaceTestFunction<RealType>;
    testFct->setRadius ( 1.0 / 3.0 );
    testFct->setCenterWC ( aol::Vec3<RealType> ( 0.5, 0.5, 0.5 ) );
    testFct->setTangFactor ( 0.1 );
    return ( testFct );
  }
  break;


  case 20: { // spherical interface centered inside bounding cube with smooth transition to zero inside the boundary; periodic cell for swiss cheese
    SphericalPlateauInterfaceTestFunction<RealType>* testFct = new SphericalPlateauInterfaceTestFunction<RealType>;
    return ( testFct );
  }
  break;

  case 21: { // periodic version of SphericalPlateauITF
    SwissCheeseInterfaceTestFunction<RealType>* testFct = new SwissCheeseInterfaceTestFunction<RealType> ( 3 );
    return ( testFct );
  }
  break;


  case 30: { // piecewise polynomial 2nd order
    PlanarPolynomialTestFunction<RealType>* testFct = new PlanarPolynomialTestFunction<RealType>;
    return ( testFct );
  }
  break;

  case 100:
  case 16: { // cylindrically polynomial
    CylindricalPolynomialTestFunction<RealType> *testFct = new CylindricalPolynomialTestFunction<RealType>;
    return ( testFct );
  }
  break;

  case 32: { // spherically-symmetric piecewise polynomial 2nd order
    SphericalPolynomialTestFunction<RealType>* testFct = new SphericalPolynomialTestFunction<RealType>;
    return ( testFct );
  }
  break;

  case 101:
  case 33: { // spherically-symmetric piecewise polynomial modulated
    SphericalPolynomialTestFunction<RealType>* testFct = new SphericalPolynomialTestFunction<RealType>;
    testFct->setTangFactor ( 0.1 );
    return ( testFct );
  }
  break;

  case 102: { // periodic spherically-symmetric polynomial-linear-exponential-plateau
    SwissCheesePolynomialInterfaceTestFunction<RealType>* testFct = new SwissCheesePolynomialInterfaceTestFunction<RealType> ( 1 );
    return ( testFct );
  }

  case 103: { // periodic spherically-symmetric polynomial-linear-exponential-plateau
    SwissCheesePolynomialInterfaceTestFunction<RealType>* testFct = new SwissCheesePolynomialInterfaceTestFunction<RealType> ( 3 );
    return ( testFct );
  }

  case 200: { // periodic spherically-symmetric polynomial-linear-exponential-plateau
    SphericalPolynomialDoubleKinkTestFunction<RealType>* testFct = new SphericalPolynomialDoubleKinkTestFunction<RealType>;
    testFct->setTangFactor ( 0.1 );
    return ( testFct );
  }

  default:
    cerr << "Error: invalid type" << endl;
    return ( NULL );
  }
}


template< typename RealType >
class CFEStructureAnalytic : public CFEStructure<RealType> {
private:
  InterfaceTestFunction< RealType > *_pITF;
  qc::OTFILexMapper<qc::QC_3D> _iMap;
  static const int _numSteps = 15;
  static const RealType _threshold;

public:
  CFEStructureAnalytic ( ) : _pITF ( NULL ) { }

  void setITFPointer ( InterfaceTestFunction<RealType> *pITF, const int width ) {
    _pITF = pITF;
    _iMap.resize ( width, width, width );
  }

  RealType getValue ( const qc::CoordType& pos ) const {
    aol::Vec3<RealType> rpos ( pos );
    return ( _pITF->interfaceValueAt_global( rpos ) );
  }

  RealType getCutRelationBetween ( const qc::CoordType& pos0, const qc::CoordType &pos1 ) const {
    // via secant search

    RealType b0 = aol::NumberTrait<RealType>::zero, b1 = aol::NumberTrait<RealType>::one;
    const aol::Vec3<RealType> psr0 ( pos0 ), psr1 ( pos1 );
    RealType v0 = _pITF->interfaceValueAt_global ( psr0 ), v1 = _pITF->interfaceValueAt_global ( psr1 );

    if ( differentSign ( v0, v1 ) ) {
      int i = 0;
      do {
        const RealType b2 = ( b0 * v1 - b1 * v0 ) / ( v1 - v0 );
        b0 = b1;
        v0 = v1;
        b1 = b2;
        {
          const aol::Vec3<RealType> p1 = psr0 + b1 * ( psr1 - psr0 );
          v1 = _pITF->interfaceValueAt_global ( p1 );
        }
        ++i;
      } while ( ( i < _numSteps ) && ( fabs ( v1 ) > _threshold ) );

      return ( b1 );
    } else { // no zero on edge
      return ( aol::NumberTrait<RealType>::NaN );
    }
  }

  aol::Vec3<RealType> getInterfaceNormal ( const aol::Vec3<RealType> &pos ) const {
    aol::Vec3<RealType> n;
    _pITF->gradientInterfaceValueAt_global ( pos, n );
    n.normalize();
    return ( n );
  }

};


template<> const float CFEStructureAnalytic<float>::_threshold = 1.0e-7;
template<> const double CFEStructureAnalytic<double>::_threshold = 1.0e-15;
template<> const long double CFEStructureAnalytic<long double>::_threshold = 1.0e-15;



template< typename RealType >
class SphericalInterfaceBVPTestFunction : public InterfaceTestFunction<RealType> {
protected:
  RealType               _radius;    // position of interface in world coordinates
  aol::Vec3<RealType>    _centerWC;  // center of the sphere in world coordinates
  RealType               _tangFactor;

public:
  SphericalInterfaceBVPTestFunction ( ) : InterfaceTestFunction<RealType>(), _radius ( 1.0/3.0 ), _centerWC( 0.5, 0.5, 0.5 ), _tangFactor ( 0.0 ) {
  }

  virtual ~SphericalInterfaceBVPTestFunction() {
  }

public:

  void setTangFactor ( const RealType tangFactor ) {
    _tangFactor = tangFactor;
  }

public:
  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &p_wc ) const {
    const RealType radius = ( p_wc - _centerWC ).norm();
    return ( radius - _radius );
  }

  virtual void gradientInterfaceValueAt ( const aol::Vec3<RealType> &p_wc, aol::Vec3<RealType> &gradient ) const {
    gradient = ( p_wc - _centerWC ); // need not be normalized!
  }

  virtual RealType valueAt          ( const aol::Vec3<RealType> &p_wc ) const {
    RealType radius, beta, r, s, c;
    getValues ( p_wc, radius, beta, r, s, c );

    return ( c * r * s );
  }

private:
  inline void getValues ( const aol::Vec3<RealType> &p_wc, RealType &radius, RealType &beta, RealType &r, RealType &s, RealType &c ) const {
    radius = ( p_wc - _centerWC ).norm();
    beta = this->_coeffMinus / this->_coeffPlus * 2.0 / _radius;

    r = ( ( radius < _radius )  ?  ( - exp ( 1 - aol::Sqr ( radius / _radius ) ) )  :  ( beta * ( radius - _radius ) - 1.0 ) );

    RealType x, y, z, rxy, ryz, rzx;
    getSValues ( p_wc, x, y, z, rxy, ryz, rzx );
    const RealType
      sigmaxy = ( rxy > 1e-12 ? y / rxy : 0 ),
      sigmayz = ( ryz > 1e-12 ? z / ryz : 0 ),
      sigmazx = ( rzx > 1e-12 ? x / rzx : 0 );

    s = 1 - _tangFactor * sigmaxy * sigmayz * sigmazx;

    c = cutOffValue ( radius );
  }

  inline void getSValues ( const aol::Vec3<RealType> &p_wc, RealType &x, RealType &y, RealType &z, RealType &rxy, RealType &ryz, RealType &rzx ) const {
    x = p_wc[0] - _centerWC[0];
    y = p_wc[1] - _centerWC[1];
    z = p_wc[2] - _centerWC[2];
    rxy = sqrt ( x*x + y*y );
    ryz = sqrt ( y*y + z*z );
    rzx = sqrt ( z*z + x*x );
  }

  inline RealType transFct ( const RealType x ) const {
    return ( exp ( -1.0 / ( 1 - x * x ) ) );
  }

  inline RealType cutOffValue ( const RealType radius ) const {
    const RealType cutOffRadius = 0.49;
    if ( radius <= 0 )
      return ( 1.0 );
    else if ( radius < cutOffRadius )
      return ( exp(1.0) * transFct ( radius / cutOffRadius ) );
    else
      return ( 0.0 );
  }

};



template< typename RealType >
  class PeriodicZeroBallInterfaceTestFunction : public PeriodicInterfaceTestFunction< RealType, SphericalInterfaceBVPTestFunction<RealType> > {
protected:
  SphericalInterfaceBVPTestFunction<RealType> _sITF;

public:
  explicit PeriodicZeroBallInterfaceTestFunction ( const RealType Freq ) : PeriodicInterfaceTestFunction< RealType, SphericalInterfaceBVPTestFunction<RealType> >( Freq ), _sITF() {
    this->setCellFct ( _sITF );
  }

  virtual ~PeriodicZeroBallInterfaceTestFunction() {}

};



template< typename RealType >
class FrickelInterfaceTestFunction : public InterfaceTestFunction<RealType> {
protected:
  const qc::ScalarArray<RealType, qc::QC_3D> &_levelset;

  virtual RealType interfaceValueAt ( const aol::Vec3<RealType> &/*p_wc*/ ) const {
    throw aol::Exception ( "FrickelInterfaceTestFunction::interfaceValueAt must not be called", __FILE__, __LINE__ );
  }
  virtual void     gradientInterfaceValueAt ( const aol::Vec3<RealType> &/*p_wc*/, aol::Vec3<RealType> &/*gradient*/ ) const {
    throw aol::Exception ( "FrickelInterfaceTestFunction::gradientInterfaceValueAt must not be called", __FILE__, __LINE__ );
  }
  virtual RealType valueAt          ( const aol::Vec3<RealType> &/*p_wc*/ ) const{
    throw aol::Exception ( "FrickelInterfaceTestFunction::valueAt must not be called", __FILE__, __LINE__ );
  }

public:
  FrickelInterfaceTestFunction ( const RealType CoeffPlus, const RealType CoeffMinus, const qc::ScalarArray<RealType, qc::QC_3D> &Levelset ) : InterfaceTestFunction<RealType> ( ), _levelset ( Levelset ) {
    this->setCoeffs ( CoeffPlus, CoeffMinus );
  }

  RealType coefficientValueAtGlobal ( const aol::Vec3<RealType> &p ) const {
    return ( _levelset.interpolate ( p ) < 0 ? this->_coeffMinus : this->_coeffPlus );
  }
};


}

namespace aol {

template< typename ConfiguratorType >
class MyMFEWeightedStiffOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, MyMFEWeightedStiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;

public:
  MyMFEWeightedStiffOp ( const typename ConfiguratorType::InitType &Initializer, const tpcfe::InterfaceTestFunction<RealType> &Itf, aol::OperatorType OpType = aol::ONTHEFLY )
    : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, MyMFEWeightedStiffOp<ConfiguratorType> > ( Initializer, OpType ),
      _itf ( Itf ) {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int /*QuadPoint*/, const typename ConfiguratorType::VecType &RefCoord ) const {
    typename ConfiguratorType::VecType imageCOQuadPoint ( El ); // explicit conversion from Vec3<integer> to Vec3<floating>
    imageCOQuadPoint += RefCoord;
    return ( _itf.coefficientValueAtGlobal ( imageCOQuadPoint ) );
  }

protected:
  const tpcfe::InterfaceTestFunction<RealType> &_itf;
};

}

#endif
