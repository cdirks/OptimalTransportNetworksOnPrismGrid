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

#ifndef __AFFINE_H
#define __AFFINE_H

// in this file, we provide test functions and test interfaces for the
// tpCFE tests and a function to compute error norms between computed
// and analytical solutions


#include <matrix.h>
#include <quoc.h>
#include <scalarArray.h>
#include <tpCFEGrid.h>
#include <iostream>
#include <progressBar.h>

// this name is misleading - the spherical function is not affine, not even spherical-affine!

namespace tpcfe {

template <typename RealType>
class AffineFunction {
public:
  enum Type { ALIGNED_PLANAR_INTERFACE_X,
              ALIGNED_PLANAR_INTERFACE_Y,
              ALIGNED_PLANAR_INTERFACE_Z,
              NON_ALIGNED_PLANAR_INTERFACE,
              SPHERICAL_INTERFACE,
              NO_INTERFACE,                    // domain without interface
              SSOS_INTERFACE,
              NO_INTERFACE_SET              }; // interface not set yet

protected:
  int                     _width; // of image
  Type                    _type;
  qc::Comp                _axis;
  RealType                _radius;
  RealType                _normalDCPlus, _normalDCMinus; // diff coeff in _normal direction on either side
  RealType                _tangDC1, _tangDC2;            // diff coeff in _tangential direction
  RealType                _interfacePos;
  RealType                _theta;
  RealType                _factor;
  aol::Matrix33<RealType> A;
  bool                    rotation_matrix_set;

public:
  AffineFunction() :
      _width (              0 ),
      _type (               NO_INTERFACE_SET ),
      _axis (               qc::QC_Z ),
      _radius (             1.0 / 3.0 ),
      _normalDCPlus (       0.0 ),
      _normalDCMinus (      0.0 ),
      _tangDC1 (            0.0 ),
      _tangDC2 (            0.0 ),
      _interfacePos (       1.0 / 3.0 ),
      _theta (              M_PI / 8.0 ),
      _factor (             1.0 ),
      rotation_matrix_set ( false ) {
    A.setZero();
  }

public:
  void setCoeffs ( RealType normalDCPlus, RealType normalDCMinus, RealType tangDC1, RealType tangDC2 ) {
    _normalDCPlus = normalDCPlus;
    _normalDCMinus = normalDCMinus;
    _tangDC1 = tangDC1;
    _tangDC2 = tangDC2;
  }

  void setType ( Type type ) { _type = type; }
  void setWidth ( int width ) { _width = width; }
  void setInterfacePos ( RealType pos ) { _interfacePos = pos; }
  void setAngle ( RealType angle )  { _theta = angle; }
  void setRadius ( RealType radius ) { _radius = radius; } // radius is meant world coordinate units
  void setAxis ( qc::Comp axis ) { _axis = axis; }
  void setFactor ( RealType factor ) { _factor = factor; }

  void setTypeNumber ( int _typeNo ) {
    switch ( _typeNo ) {
    case 0:
      cout << "% Testing planar interface x" << endl;
      setType ( ALIGNED_PLANAR_INTERFACE_X );
      break;
    case 1:
      cout << "% Testing planar interface y" << endl;
      setType ( ALIGNED_PLANAR_INTERFACE_Y );
      break;
    case 2:
      cout << "% Testing planar interface z" << endl;
      setType ( ALIGNED_PLANAR_INTERFACE_Z );
      break;
    case 3:
      cout << "% Testing non-aligned planar interface" << endl;
      setType ( NON_ALIGNED_PLANAR_INTERFACE );
      break;
    case 4:
      cout << "% Testing spherical interface r=1/8" << endl;
      setType ( SPHERICAL_INTERFACE );
      setRadius ( 1.0 / 8.0 ); // this means we need to deal with interface lying on regular grid points
      break;
    case 5:
      cout << "% Testing spherical interface r=1/4" << endl;
      setType ( SPHERICAL_INTERFACE );
      setRadius ( 1.0 / 4.0 ); // this means we need to deal with interface lying on regular grid points
      break;
    case 6:
      cout << "% Testing spherical interface r=1/3" << endl;
      setType ( SPHERICAL_INTERFACE );
      setRadius ( 1.0 / 3.0 );
      break;
    case 7:
      cout << "% Testing no interface" << endl;
      setType ( NO_INTERFACE );
      break;
    case 8:
      cout << "% Testing SSOS trilinear function" << endl;
      setType ( SSOS_INTERFACE );
    default:
      cerr << "Please choose correct interface type" << endl;
      break;
    }
  }

protected:

  void setRotation_matrix() {
    switch ( _axis ) {
    case qc::QC_X:
      A.setRow ( 0, 1, 0, 0 );
      A.setRow ( 1, 0, cos ( _theta ), -sin ( _theta ) );
      A.setRow ( 2, 0, sin ( _theta ),  cos ( _theta ) );
      break;
    case qc::QC_Y:
      A.setRow ( 0,  cos ( _theta ), 0, -sin ( _theta ) );
      A.setRow ( 1,  0, 1, 0 );
      A.setRow ( 2,  sin ( _theta ),  0, cos ( _theta ) );
      break;
    case qc::QC_Z:
      A.setRow ( 0,  cos ( _theta ), -sin ( _theta ), 0 );
      A.setRow ( 1,  sin ( _theta ),  cos ( _theta ), 0 );
      A.setRow ( 2,  0, 0, 1 );
      break;
    default:
      std::cerr << "Rotation axis is not set" << std::endl;
    }
  }

  void B_wc ( aol::Vec3<RealType> p, aol::Vec3<RealType> &Bp ) {
    // expects world coordinates
    if ( !rotation_matrix_set )
      setRotation_matrix();

    p[0] -= 0.5;
    p[1] -= 0.5;
    p[2] -= 0.5;

    Bp = A * p;

#if 0
    // Check how good the back rotation is
    aol::Matrix33<RealType> B ( A );
    B.transpose();

    aol::Vec3<RealType> r = B * Bp;
    cerr.precision ( 15 );
    std::cerr << r << endl << p << endl;
    r -= p;
    std::cerr << r.norm()  << endl;
    getchar();
#endif

    Bp[0] += 0.5;
    Bp[1] += 0.5;
    Bp[2] += 0.5;

  }

  // aligned affine function
  // =======================
  RealType alignedAffineFct ( const aol::Vec3<RealType> &p, const int comp ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    aol::Vec3<RealType> a1, a2;

    a1[comp] = _normalDCPlus; a2[comp] = _normalDCMinus;
    a1[ ( comp+1 ) % 3 ] = a2[ ( comp+1 ) % 3 ] = _tangDC1;
    a1[ ( comp+2 ) % 3 ] = a2[ ( comp+2 ) % 3 ] = _tangDC2;
    const RealType b1 = 0;
    const RealType b2 = ( _normalDCPlus - _normalDCMinus ) * _interfacePos + b1;

    if ( p_wc[comp] < _interfacePos ) {
      return ( a1*p_wc + b1 );
    } else {
      return ( a2*p_wc + b2 );
    }
  }

  // ??? why _tangDC1 and not zero??

  RealType dx_i_alignedAffineFct ( const aol::Vec3<RealType> &p, const int comp, const int i ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    if ( i == comp ) { // derivative normal to the kink
      if ( p_wc[comp] < _interfacePos ) {
        return _normalDCPlus;
      } else {
        return _normalDCMinus;
      }
    } else if ( i == ( comp + 1 ) % 3 ) {
      return _tangDC1;
    } else if ( i  == ( comp + 2 ) % 3 ) {
      return _tangDC2;
    } else {
      return aol::NumberTrait<RealType>::NaN ;
    }
  }

  RealType Laplace_alignedAffineFct ( const aol::Vec3<RealType>& , const int ) {
    return ( 0.0 );
  }

  // non-aligned affine function
  // ===========================
  RealType nonAlignedAffineFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    aol::Vec3<RealType> Bp_wc;
    B_wc ( p_wc, Bp_wc );  // Compute transformed coordinate in world coordinates

    aol::Vec3<RealType> Bp = w1 * Bp_wc;

    return alignedAffineFct ( Bp, 0 );
  }

  RealType partial_x_nonAlignedAffineFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    if ( !rotation_matrix_set )
      setRotation_matrix();

    aol::Vec3<RealType> Bp_wc;
    B_wc ( p_wc, Bp_wc );

    if ( Bp_wc[0] < _interfacePos ) {
      return ( ( _normalDCPlus  * A[0][0] + _tangDC1 * A[1][0] + _tangDC2 * A[2][0] ) );
    } else {
      return ( ( _normalDCMinus * A[0][0] + _tangDC1 * A[1][0] + _tangDC2 * A[2][0] ) );
    }
  }

  RealType partial_y_nonAlignedAffineFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    if ( !rotation_matrix_set )
      setRotation_matrix();

    aol::Vec3<RealType> Bp_wc;
    B_wc ( p_wc, Bp_wc );

    if ( Bp_wc[0] < _interfacePos ) {
      return ( ( _normalDCPlus  * A[0][1] + _tangDC1 * A[1][1] + _tangDC2 * A[2][1] ) );
    } else {
      return ( ( _normalDCMinus * A[0][1] + _tangDC1 * A[1][1] + _tangDC2 * A[2][1] ) );
    }
  }

  RealType partial_z_nonAlignedAffineFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    if ( !rotation_matrix_set )
      setRotation_matrix();

    aol::Vec3<RealType> Bp_wc;
    B_wc ( p_wc, Bp_wc );

    if ( Bp_wc[0] < _interfacePos ) {
      return ( ( _normalDCPlus  * A[0][2] + _tangDC1 * A[1][2] + _tangDC2 * A[2][2] ) );
    } else {
      return ( ( _normalDCMinus * A[0][2] + _tangDC1 * A[1][2] + _tangDC2 * A[2][2] ) );
    }
  }

  RealType Laplace_nonAlignedAffineFct ( const aol::Vec3<RealType> & ) {
    return ( 0.0 );
  }

  // spherical function
  RealType sphericalFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    const aol::Vec3<RealType> center ( 0.5, 0.5, 0.5 );
    const RealType radius = ( p_wc - center ).norm();
    const RealType alpha = 1.0;
    const RealType beta = _normalDCMinus / _normalDCPlus * 2.0 * alpha / _radius;

    if ( radius < _radius ) {
      return ( alpha * exp ( 1 - aol::Sqr ( radius / _radius ) ) );
    } else {
      return ( - beta * ( radius - _radius ) + alpha );
    }
  }

  RealType partial_x_sphericalFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    const aol::Vec3<RealType> center ( 0.5, 0.5, 0.5 );
    aol::Vec3<RealType> pmc_wc = p_wc - center;
    const RealType radius = pmc_wc.norm();
    const RealType alpha = 1.0;
    const RealType beta = _normalDCMinus / _normalDCPlus * 2.0 * alpha / _radius;

    if ( radius < _radius ) {
      return ( -2.0 * ( alpha / aol::Sqr ( _radius ) ) *  exp ( 1 - aol::Sqr ( radius / _radius ) ) * pmc_wc[0] );
    } else {
      return ( - ( beta / radius ) * pmc_wc[0] );
    }

  }

  RealType partial_y_sphericalFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    const aol::Vec3<RealType> center ( 0.5, 0.5, 0.5 );
    aol::Vec3<RealType> pmc_wc = p_wc - center;
    const RealType radius = pmc_wc.norm();
    const RealType alpha = 1.0;
    const RealType beta = _normalDCMinus / _normalDCPlus * 2.0 * alpha / _radius;

    if ( radius < _radius ) {
      return ( -2.0 * ( alpha / aol::Sqr ( _radius ) ) *  exp ( 1 - aol::Sqr ( radius / _radius ) ) * pmc_wc[1] );
    } else {
      return ( - ( beta / radius ) * pmc_wc[1] );
    }

  }

  RealType partial_z_sphericalFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    const aol::Vec3<RealType> center ( 0.5, 0.5, 0.5 );
    aol::Vec3<RealType> pmc_wc = p_wc - center;
    const RealType radius = pmc_wc.norm();
    const RealType alpha = 1.0;
    const RealType beta = _normalDCMinus / _normalDCPlus * 2.0 * alpha / _radius;

    if ( radius < _radius ) {
      return ( -2.0 * ( alpha / aol::Sqr ( _radius ) ) *  exp ( 1 - aol::Sqr ( radius / _radius ) ) * pmc_wc[2] );
    } else {
      return ( - ( beta / radius ) * pmc_wc[2] );
    }

  }

  RealType Laplace_sphericalFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    const aol::Vec3<RealType> center ( 0.5, 0.5, 0.5 );
    aol::Vec3<RealType> pmc_wc = p_wc - center;
    const RealType radius = pmc_wc.norm();
    const RealType alpha = 1.0;
    const RealType beta = _normalDCMinus / _normalDCPlus * 2.0 * alpha / _radius;

    if ( radius < _radius ) {
      return ( 2.0 * ( alpha / aol::Sqr ( _radius ) )  * ( 2.0 * aol::Sqr ( radius ) / aol::Sqr ( _radius ) - 3.0 ) * exp ( 1.0 - aol::Sqr ( radius / _radius ) ) );
    } else {
      return ( - 2.0 * ( beta / radius ) );
    }

  }


  // for the testFct, we have four different cases:
  // it may have zero- or nonzero-Neumann-BC (N0, N1)
  // the function may have zero- or nonzero-Dirichlet-BC (D0, D1)
#define N1D1   bla
  RealType testFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

#if   defined N0D0
    return (  ( 1.0 - cos ( 2 * M_PI * p_wc[0] ) ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[1] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[2] ) ) );
#elif defined N0D1
    return ( 2.0 +  ( 1.0 - cos ( 2 * M_PI * p_wc[0] ) ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[1] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[2] ) ) );
#elif defined N1D0
    return ( 2.0 + sin ( 2 * M_PI * p_wc[0] ) * sin ( 2 * M_PI * p_wc[1] ) * sin ( 2 * M_PI * p_wc[2] ) );
#elif defined N1D1
    return ( aol::Cub ( p_wc[0] ) + aol::Sqr ( p_wc[1] ) + ( p_wc[2] ) );
#else
    return ( 0.0 );
#endif

  }

  RealType partial_x_testFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

#if   defined N0D0
    return ( 2 * M_PI * sin ( 2 * M_PI * p_wc[0] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[1] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[2] ) ) );
#elif defined N0D1
    return ( 2 * M_PI * sin ( 2 * M_PI * p_wc[0] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[1] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[2] ) ) );
#elif defined N1D0
    return ( 2 * M_PI * cos ( 2 * M_PI * p_wc[0] ) * sin ( 2 * M_PI * p_wc[1] ) * sin ( 2 * M_PI * p_wc[2] ) );
#elif defined N1D1
    return ( 3.0 * aol::Sqr ( p_wc[0] ) );
#endif

  }

  RealType partial_y_testFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

#if   defined N0D0
    return ( 2 * M_PI * sin ( 2 * M_PI * p_wc[1] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[0] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[2] ) ) );
#elif defined N0D1
    return ( 2 * M_PI * sin ( 2 * M_PI * p_wc[1] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[0] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[2] ) ) );
#elif defined N1D0
    return ( 2 * M_PI * cos ( 2 * M_PI * p_wc[1] ) * sin ( 2 * M_PI * p_wc[0] ) * sin ( 2 * M_PI * p_wc[2] ) );
#elif defined N1D1
    return ( 2.0 * p_wc[1] );
#else
    return ( 0.0 );
#endif

  }

  RealType partial_z_testFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

#if   defined N0D0
    return ( 2 * M_PI * sin ( 2 * M_PI * p_wc[2] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[0] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[1] ) ) );
#elif defined N0D1
    return ( 2 * M_PI * sin ( 2 * M_PI * p_wc[2] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[0] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[1] ) ) );
#elif defined N1D0
    return ( 2 * M_PI * cos ( 2 * M_PI * p_wc[2] ) * sin ( 2 * M_PI * p_wc[0] ) * sin ( 2 * M_PI * p_wc[1] ) );
#elif defined N1D1
    return ( 1.0 );
#else
    return ( 0.0 );
#endif

  }

  RealType Laplace_testFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

#if   defined N0D0
    return ( 4 * M_PI * M_PI * ( cos ( 2 * M_PI * p_wc[0] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[1] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[2] ) ) +
                                 cos ( 2 * M_PI * p_wc[1] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[0] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[2] ) ) +
                                 cos ( 2 * M_PI * p_wc[2] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[0] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[1] ) ) ) );
#elif defined N0D1
    return ( 4 * M_PI * M_PI * ( cos ( 2 * M_PI * p_wc[0] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[1] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[2] ) ) +
                                 cos ( 2 * M_PI * p_wc[1] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[0] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[2] ) ) +
                                 cos ( 2 * M_PI * p_wc[2] ) *  ( 1.0 - cos ( 2 * M_PI * p_wc[0] ) ) * ( 1.0 - cos ( 2 * M_PI * p_wc[1] ) ) ) );
#elif defined N1D0
    return ( -12 * M_PI * M_PI *  sin ( 2 * M_PI * p_wc[0] ) * sin ( 2 * M_PI * p_wc[1] ) * sin ( 2 * M_PI * p_wc[2] ) );
#elif defined N1D1
    return ( 6.0 * p_wc[0] + 2.0 );
#else
    return ( 0.0 );
#endif

  }

  // SSOS example that will hopefully show bad convergence with TPOS computation of weights
  RealType SSOSFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    if ( p_wc[0] < _interfacePos ) {
      return ( _normalDCPlus  * ( p_wc[0] - _interfacePos ) * ( p_wc[1] + 1 ) * ( p_wc[2] + 1 ) );
    } else {
      return ( _normalDCMinus * ( p_wc[0] - _interfacePos ) * ( p_wc[1] + 1 ) * ( p_wc[2] + 1 ) );
    }
  }

  RealType partial_x_SSOSFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    if ( p_wc[0] < _interfacePos ) {
      return ( _normalDCPlus  * (           1.0           ) * ( p_wc[1] + 1 ) * ( p_wc[2] + 1 ) );
    } else {
      return ( _normalDCMinus * (           1.0           ) * ( p_wc[1] + 1 ) * ( p_wc[2] + 1 ) );
    }
  }

  RealType partial_y_SSOSFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    if ( p_wc[0] < _interfacePos ) {
      return ( _normalDCPlus  * ( p_wc[0] - _interfacePos ) * (     1.0     ) * ( p_wc[2] + 1 ) );
    } else {
      return ( _normalDCMinus * ( p_wc[0] - _interfacePos ) * (     1.0     ) * ( p_wc[2] + 1 ) );
    }
  }

  RealType partial_z_SSOSFct ( const aol::Vec3<RealType> &p ) {
    const RealType w1 = ( _width - 1.0 ) ;
    const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

    if ( p_wc[0] < _interfacePos ) {
      return ( _normalDCPlus  * ( p_wc[0] - _interfacePos ) * ( p_wc[1] + 1 ) * (    1.0     ) );
    } else {
      return ( _normalDCMinus * ( p_wc[0] - _interfacePos ) * ( p_wc[1] + 1 ) * (    1.0     ) );
    }
  }

  RealType Laplace_SSOSFct ( const aol::Vec3<RealType> & ) {
    return ( 0.0 );
  }

public:
  /** Check how close a given point is to the interface
   */
  RealType distanceToInterface ( const aol::Vec3<RealType> &p ) {
    switch ( _type ) {
    case NON_ALIGNED_PLANAR_INTERFACE: {
      const RealType w1 = ( _width - 1.0 ) ;
      const aol::Vec3<RealType> p_wc ( p[0] / w1, p[1] / w1, p[2] / w1 );

      if ( !rotation_matrix_set ) setRotation_matrix();

      aol::Vec3<RealType> r;
      B_wc ( p_wc, r );

      r[0] -= _interfacePos;
      cerr.precision ( 20 );
      std::cerr << "Distance to interface " << r[0] << endl;
      return r[0];
    }
    default:
      std::cerr << "isPointOnInterface not implemented for this type" << std::endl;
    }
    return aol::NumberTrait<RealType>::NaN;
  }

  /** Return the value of the function at coordinate p
   */
  RealType valueAt ( const aol::Vec3<RealType> &p ) {
    switch ( _type ) {
    case ALIGNED_PLANAR_INTERFACE_X:
      return _factor*alignedAffineFct ( p, 0 );  break;
    case ALIGNED_PLANAR_INTERFACE_Y:
      return _factor*alignedAffineFct ( p, 1 );  break;
    case ALIGNED_PLANAR_INTERFACE_Z:
      return _factor*alignedAffineFct ( p, 2 );  break;
    case NON_ALIGNED_PLANAR_INTERFACE:
      return _factor*nonAlignedAffineFct ( p ); break;
    case SPHERICAL_INTERFACE:
      return _factor*sphericalFct ( p );        break;
    case NO_INTERFACE:
      return _factor*testFct ( p );             break;
    case SSOS_INTERFACE:
      return _factor*SSOSFct ( p );             break;
    default:
      return ( aol::NumberTrait<RealType>::NaN );
    }
    return ( aol::NumberTrait<RealType>::NaN );
  }

  RealType partial_x_valueAt ( const aol::Vec3<RealType> &p ) {
    switch ( _type ) {
    case ALIGNED_PLANAR_INTERFACE_X:
      return _factor*dx_i_alignedAffineFct ( p, 0, 0 );   break;
    case ALIGNED_PLANAR_INTERFACE_Y:
      return _factor*dx_i_alignedAffineFct ( p, 1, 0 );   break;
    case ALIGNED_PLANAR_INTERFACE_Z:
      return _factor*dx_i_alignedAffineFct ( p, 2, 0 );   break;
    case NON_ALIGNED_PLANAR_INTERFACE:
      return _factor*partial_x_nonAlignedAffineFct ( p ); break;
    case SPHERICAL_INTERFACE:
      return _factor*partial_x_sphericalFct ( p );        break;
    case NO_INTERFACE:
      return _factor*partial_x_testFct ( p );             break;
    case SSOS_INTERFACE:
      return _factor*partial_x_SSOSFct ( p );             break;
    default:
      return ( aol::NumberTrait<RealType>::NaN );
    }
    return ( aol::NumberTrait<RealType>::NaN );
  }

  RealType partial_y_valueAt ( const aol::Vec3<RealType> &p ) {
    switch ( _type ) {
    case ALIGNED_PLANAR_INTERFACE_X:
      return _factor*dx_i_alignedAffineFct ( p, 0, 1 );   break;
    case ALIGNED_PLANAR_INTERFACE_Y:
      return _factor*dx_i_alignedAffineFct ( p, 1, 1 );   break;
    case ALIGNED_PLANAR_INTERFACE_Z:
      return _factor*dx_i_alignedAffineFct ( p, 2, 1 );   break;
    case NON_ALIGNED_PLANAR_INTERFACE:
      return _factor*partial_y_nonAlignedAffineFct ( p ); break;
    case SPHERICAL_INTERFACE:
      return _factor*partial_y_sphericalFct ( p );        break;
    case NO_INTERFACE:
      return _factor*partial_y_testFct ( p );             break;
    case SSOS_INTERFACE:
      return _factor*partial_y_SSOSFct ( p );             break;
    default:
      return ( aol::NumberTrait<RealType>::NaN );
    }
    return ( aol::NumberTrait<RealType>::NaN );
  }

  RealType partial_z_valueAt ( const aol::Vec3<RealType> &p ) {
    switch ( _type ) {
    case ALIGNED_PLANAR_INTERFACE_X:
      return _factor*dx_i_alignedAffineFct ( p, 0, 2 );   break;
    case ALIGNED_PLANAR_INTERFACE_Y:
      return _factor*dx_i_alignedAffineFct ( p, 1, 2 );   break;
    case ALIGNED_PLANAR_INTERFACE_Z:
      return _factor*dx_i_alignedAffineFct ( p, 2, 2 );   break;
    case NON_ALIGNED_PLANAR_INTERFACE:
      return _factor*partial_z_nonAlignedAffineFct ( p ); break;
    case SPHERICAL_INTERFACE:
      return _factor*partial_z_sphericalFct ( p );        break;
    case NO_INTERFACE:
      return _factor*partial_z_testFct ( p );             break;
    case SSOS_INTERFACE:
      return _factor*partial_z_SSOSFct ( p );             break;
    default:
      return ( aol::NumberTrait<RealType>::NaN );
    }
    return ( aol::NumberTrait<RealType>::NaN );
  }

  // evaluation at position given in world coordinates
  inline RealType Laplace_valueAt ( const aol::Vec3<RealType> &p ) {
    switch ( _type ) {
    case ALIGNED_PLANAR_INTERFACE_X:
      return _factor*Laplace_alignedAffineFct ( p, 0 );  break;
    case ALIGNED_PLANAR_INTERFACE_Y:
      return _factor*Laplace_alignedAffineFct ( p, 1 );  break;
    case ALIGNED_PLANAR_INTERFACE_Z:
      return _factor*Laplace_alignedAffineFct ( p, 2 );  break;
    case NON_ALIGNED_PLANAR_INTERFACE:
      return _factor*Laplace_nonAlignedAffineFct ( p ); break;
    case SPHERICAL_INTERFACE:
      return _factor*Laplace_sphericalFct ( p );        break;
    case NO_INTERFACE:
      return _factor*Laplace_testFct ( p );             break;
    case SSOS_INTERFACE:
      return _factor*Laplace_SSOSFct ( p );             break;
    default:
      return ( aol::NumberTrait<RealType>::NaN );
    }
    return ( aol::NumberTrait<RealType>::NaN );
  }

  // TODO: using pointers here does not really make sense ...
  void createInterface ( qc::ScalarArray<RealType, qc::QC_3D> *myInterface,
                         qc::ScalarArray<RealType, qc::QC_3D> *affine = NULL,
                         qc::AArray<RealType, qc::QC_3D> *diffCoeff = NULL ) {

    if ( myInterface ) myInterface->setZero();

    const RealType w1 = ( _width - 1.0 ) ;

    aol::ProgressBar<> pb ( "Creating affine fct" );
    pb.start ( aol::Cub ( _width ) );

    for ( int z = 0; z < _width; z++ ) {
      for ( int y = 0; y < _width; y++ ) {
        for ( int x = 0; x < _width; x++ ) {
          aol::Vec3<RealType> p ( x, y, z ), p_wc ( x / w1, y / w1, z / w1 );

          RealType interface_value = 0.0;
          pb++;

          switch ( _type ) {
          case SPHERICAL_INTERFACE: {
            const aol::Vec3<RealType> center_wc ( 0.5, 0.5, 0.5 );
            const RealType radius = ( p_wc - center_wc ).norm();

            interface_value = _factor * ( radius - _radius );
            break;
          }
          case NON_ALIGNED_PLANAR_INTERFACE: {
            if ( !rotation_matrix_set )
              setRotation_matrix();

            aol::Vec3<RealType> Bp_wc;
            B_wc ( p_wc, Bp_wc );

            interface_value = _factor * ( Bp_wc[0] - _interfacePos );
            break;
          }
          case ALIGNED_PLANAR_INTERFACE_X: {
            interface_value = _factor * (  p_wc[0] - _interfacePos );
            break;
          }
          case ALIGNED_PLANAR_INTERFACE_Y: {
            interface_value = _factor * (  p_wc[1] - _interfacePos );
            break;
          }
          case ALIGNED_PLANAR_INTERFACE_Z: {
            interface_value = _factor * (  p_wc[2] - _interfacePos );
            break;
          }
          case NO_INTERFACE: {
            interface_value = _factor * ( 1.0 ); // minus ?
            break;
          }
          case SSOS_INTERFACE: {
            interface_value = _factor * (  p_wc[0] - _interfacePos );
            break;
          }
          default:
            cerr << "ERROR: no interface set" << endl;
          }

          if ( myInterface ) {
            myInterface->set ( x, y, z, _factor*interface_value ); // now twice multiplied by _factor -- why at all?
          }

          if ( diffCoeff ) {
            if ( interface_value < 0 ) {
              diffCoeff->getRef ( x, y, z ) = _normalDCMinus;
            } else {
              diffCoeff->getRef ( x, y, z ) = _normalDCPlus;
            }
          }

          if ( affine ) {
            // valueAt already multiplies with _factor
            const RealType affine_value = valueAt ( p );
            affine->set ( x, y, z, affine_value );
          }

        }
      }
    }
    std::cerr << std::endl;
  }

  void dump ( ostream &out ) {
    switch ( _type ) {
    case ALIGNED_PLANAR_INTERFACE_X:
      out << "ALIGNED_PLANAR_INTERFACE_X at pos " << _interfacePos;
      break;
    case ALIGNED_PLANAR_INTERFACE_Y:
      out << "ALIGNED_PLANAR_INTERFACE_Y at pos " << _interfacePos;
      break;
    case ALIGNED_PLANAR_INTERFACE_Z:
      out << "ALIGNED_PLANAR_INTERFACE_Z at pos " << _interfacePos;
      break;
    case NON_ALIGNED_PLANAR_INTERFACE:
      out << "NON_ALIGNED_PLANAR_INTERFACE at pos " << _interfacePos << ", angle " << _theta << ", rotated around ";
      switch ( _axis ) {
      case qc::QC_Z: out << " Z-axis"; break;
      case qc::QC_Y: out << " Y-axis"; break;
      case qc::QC_X: out << " X-axis"; break;
      default:
        out << " not specified axis"; break;
      }
      break;
    case SPHERICAL_INTERFACE:
      out << "SPHERICAL_INTERFACE with radius " << _radius;
      break;
    case NO_INTERFACE:
      out << "NO_INTERFACE";
      break;
    case SSOS_INTERFACE:
      out << "SSOS_INTERFACE";
      break;
    case NO_INTERFACE_SET:
      out << "NO INTERFACE SET";
      break;
    default:
      throw aol::Exception ( "unhandled case", __FILE__, __LINE__ );
    }
    out << " normalDCPlus = " << _normalDCPlus;
    out << " normalDCMinus = " << _normalDCMinus;
    out << std::endl;
  }

};

template <typename RealType>
std::ostream &operator<< ( ostream &out, AffineFunction<RealType> &a ) {
  a.dump ( out );
  return out;
}


#define LEAVE_OUT_VERIFICATION_OF_APPROXIMATION_FUNCTION
// #define SUPPRESS_POINTWISE_PROBLEM_WARNING


// Linf difference is currently evaluated at the points used whereas we use midpoint integration for the L2 and H1 differences
template< typename GridType >
#ifndef SUPPRESS_POINTWISE_PROBLEM_WARNING
void compute_L2_H1_difference_to_utrue ( qc::ScalarArray<typename GridType::RealType, qc::QC_3D> &arg, const GridType &grid, AffineFunction<typename GridType::RealType> &affineFct, typename GridType::RealType &L2_difference, typename GridType::RealType &H1_difference, typename GridType::RealType &Linf_difference, const int DB = -1e6, const int DC = -1e6, const typename GridType::RealType L2_PRINT_THRES = 1.0e-6, const typename GridType::RealType H1_PRINT_THRES = 1.0e-6 ) {
#else
void compute_L2_H1_difference_to_utrue ( qc::ScalarArray<typename GridType::RealType, qc::QC_3D> &arg, const GridType &grid, AffineFunction<typename GridType::RealType> &affineFct, typename GridType::RealType &L2_difference, typename GridType::RealType &H1_difference, typename GridType::RealType &Linf_difference, const int DB = -1e6, const int DC = -1e6 ) {
#endif
  typedef typename GridType::RealType RealType;

  RealType
    _L2_difference = 0.0,
    _Linf_difference = 0.0,
    _H1x_difference = 0.0,
    _H1y_difference = 0.0,
    _H1z_difference = 0.0;

  const int width = grid.getWidth();
  const RealType h = grid.H();

  // Run through all elements

  aol::ProgressBar<> pb ( "Computing error" );
  pb.start ( grid.getNumberOfElements() );

  const qc::GridSize<qc::QC_3D> gridSize ( grid );
  for ( typename GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
    pb++;

    tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

    const int cfetype = el.cfeType();
    el.computeAssembleData ( grid );

    //! \todo fixme
    if ( true ) {
      // Sum up the error for each tetrahedron
      for ( CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {

        const CFETetra<RealType> *t = *tit;
        RealType pointValue[4];
        aol::Vec3<RealType> points[4];

        // Get interpolated values
        for ( int i = 0; i < 4; i++ ) {
          pointValue[i] = 0.0;
          int gIdx0, gIdx1;

          gIdx0 = el._globalIndex[ ( *t ) ( i,0 ) ];
          if ( ( *t ) ( i, 1 ) == 11 ) {
            gIdx1 = -1;
          } else {
            gIdx1 = el._globalIndex[ ( *t ) ( i,1 ) ];
          }

          if ( gIdx1 == -1 ) { // non-interpolated node
            pointValue[i] = arg.get ( gIdx0 );
          } else {           // interpolated node
            const typename GridType::VirtualNodeType &vn = grid.getVirtualNodeRef ( gIdx0, gIdx1 );
#ifdef VERBOSE
            cout << "number of constraints = " << vn._numOfConstraints << endl;
#endif
            for ( unsigned char j = 0; j < vn._numOfConstraints; ++j ) {
#ifdef VERBOSE
              cout << "constraint " << vn._constraint->get ( j ) << ", weight = " << vn._weight->get ( j ) << endl;
#endif
              pointValue[i] += vn._weight->get ( j ) * arg.get ( vn._constraint->get ( j ) ) ;  //  not sure whether this is correct ...
              // TP: This seems to be perfect. It is exactly the way vn.extrapolate does it.
            }
          }
        }

        // Find value at barycenter for midpoint integration - for pcw affine fct, arithmetic mean is exact.
        const RealType interValue = ( pointValue[0] + pointValue[1] + pointValue[2] + pointValue[3] ) / 4.0 ;

        // Get real values
        aol::Vec3<RealType> coord, p; // coord is used for barycenter in global (image) coordinates.
        coord.setZero(); p.setZero();
        for ( int i = 0; i < 4; i++ ) {
          t -> computeGlobalCoordinate ( p, el, i );

          points[i] = p;

          const RealType difference = fabs ( pointValue[i] - affineFct.valueAt ( p ) );
          if ( difference > _Linf_difference ) {
            _Linf_difference = difference;
          }

          coord += p;
        }
        coord /= 4.0;

        const RealType realValue = affineFct.valueAt ( coord );

        // this will be used to verify whether the function interpolate_data_at_position_wc works correctly.
#ifndef LEAVE_OUT_VERIFICATION_OF_APPROXIMATION_FUNCTION
        const aol::Vec3<RealType> wc_coord = coord * grid.H();
        const RealType inter_value_2 = grid.interpolateDataAtPositionWC ( arg, wc_coord );
        cerr << interValue << " - " << inter_value_2 << " = " << interValue - inter_value_2 << endl;
#endif


#ifndef SUPPRESS_POINTWISE_PROBLEM_WARNING
        if ( aol::Sqr ( pointValue[0] - affineFct.valueAt ( points[0] ) ) +
             aol::Sqr ( pointValue[1] - affineFct.valueAt ( points[1] ) ) +
             aol::Sqr ( pointValue[2] - affineFct.valueAt ( points[2] ) ) +
             aol::Sqr ( pointValue[3] - affineFct.valueAt ( points[3] ) )   > L2_PRINT_THRES ) {
          cerr.precision ( 6 );
          cerr << cfetype << ": "
          << pointValue[0] << " "  << affineFct.valueAt ( points[0] ) << "    #    "
          << pointValue[1] << " "  << affineFct.valueAt ( points[1] ) << "    #    "
          << pointValue[2] << " "  << affineFct.valueAt ( points[2] ) << "    #    "
          << pointValue[3] << " "  << affineFct.valueAt ( points[3] ) << " baryc:  "
          << interValue << " - " << realValue << " = " << interValue - realValue << endl;
        }
#endif

        // DB is the number of slices near the boundary we want to ignore in the computation of L2 and H1 norms. DB != 0 is cheated ...
        // DC is the minimum linf difference from the center for points to be considered in the computation. (this is used for leaving out the "peak" of the spherical function.
        if ( ( p[0] >= DB ) && ( p[0] <= width - DB ) && ( p[1] >= DB ) && ( p[1] <= width - DB ) && ( p[2] >= DB ) && ( p[2] <= width - DB ) &&
             ! ( p[0] >= width / 2 - DC && p[0] <= width / 2 + DC && p[1] >= width / 2 - DC && p[1] <= width / 2 + DC && p[2] >= width / 2 - DC && p[2] <= width / 2 + DC ) ) {
          // add local L2 error
          _L2_difference += t->getVolume() * aol::Cub ( grid.H() )  *  aol::Sqr ( realValue - interValue ); // volume is computed for h = 1, so factor h^3 okay.

          // compute and add local H1 error

          // for this purpose, we compute directional derivatives along some edges of the tetrahedron to compute the gradient of our piecewise affine-linear function (at the barycenter)

          RealType len[3];
          aol::Vec3<RealType> dir[3], con[3];
          for ( int i = 0; i < 3; i++ ) {
            con[i]  =  h * ( points[i] - points[3] );
            len[i]  =  con[i].norm();
            dir[i]  =  ( 1.0 / len[i] ) * con[i] ;
          }

          aol::Vec3<RealType>
            directional_derivatives ( ( pointValue[0] - pointValue[3] ) / len[0] ,
                                      ( pointValue[1] - pointValue[3] ) / len[1] ,
                                      ( pointValue[2] - pointValue[3] ) / len[2]   ),
            gradient;


          aol::Matrix33<RealType> A, B;
          A.setRow ( 0, dir[0][0], dir[0][1], dir[0][2] );
          A.setRow ( 1, dir[1][0], dir[1][1], dir[1][2] );
          A.setRow ( 2, dir[2][0], dir[2][1], dir[2][2] );

          try {
            B = A.inverse();                                 // this should work: for non-degenerate tetrahedra, dir[i] are linearly independent, thus A is invertible
          } catch ( aol::Exception &e ) {
            e.dump();
            continue;
          }

          gradient = B * directional_derivatives;

#ifndef SUPPRESS_POINTWISE_PROBLEM_WARNING
          if ( aol::Sqr ( gradient[0] - affineFct.partial_x_valueAt ( coord ) ) +
               aol::Sqr ( gradient[1] - affineFct.partial_y_valueAt ( coord ) ) +
               aol::Sqr ( gradient[2] - affineFct.partial_z_valueAt ( coord ) )    > H1_PRINT_THRES    ) {
            cerr.precision ( 6 );
            cerr << cfetype << " " << _L2_difference <<  " ## "
            << gradient[0] << " / " << affineFct.partial_x_valueAt ( coord ) << " =  " << gradient[0] / affineFct.partial_x_valueAt ( coord ) << ",  "
            << gradient[1] << " / " << affineFct.partial_y_valueAt ( coord ) << " =  " << gradient[1] / affineFct.partial_y_valueAt ( coord ) << ",  "
            << gradient[2] << " / " << affineFct.partial_z_valueAt ( coord ) << " =  " << gradient[2] / affineFct.partial_z_valueAt ( coord ) << endl;
          }
#endif

          _H1x_difference += t->getVolume() * aol::Cub ( grid.H() ) * aol::Sqr ( gradient[0] - affineFct.partial_x_valueAt ( coord ) );
          _H1y_difference += t->getVolume() * aol::Cub ( grid.H() ) * aol::Sqr ( gradient[1] - affineFct.partial_y_valueAt ( coord ) );
          _H1z_difference += t->getVolume() * aol::Cub ( grid.H() ) * aol::Sqr ( gradient[2] - affineFct.partial_z_valueAt ( coord ) );

        } else {
#ifdef VERBOSE
          cerr << "Ignoring point (" << p[0] << ", " << p[1] << ", " << p[2] << ") in comparison." << endl;
#endif
        }

      }
    }
  }

  std::cerr << std::endl;

  L2_difference = sqrt ( _L2_difference );
  H1_difference = sqrt ( _L2_difference + _H1x_difference + _H1y_difference + _H1z_difference );
  Linf_difference = _Linf_difference;
}

} // end namespace tpcfe
#endif
