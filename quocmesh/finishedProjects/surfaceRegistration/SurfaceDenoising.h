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

#ifndef __SURFACEREGISTRATION_H
#define __SURFACEREGISTRATION_H

#include "signedDistanceOp.h"
#include "TOF.h"
#include "deformations.h"
#include <op.h>


template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class SDFBase {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  /// The 2D grid (2D range data, parameter domain)
  const typename ConfiguratorType2D::InitType & _grid2D;
  /// The 3D grid (3D world coord data, SDF of reference surface)
  const typename ConfiguratorType3D::InitType & _grid3D;

  /// The discrete function of the SDF of the reference surface
  const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFuncs;

  /// The converter from range data to 3D world coordinates
  const RangeToWorldCoordTransformation<RealType> & _rangeToWorldCoord;

  /// Count # quad points that are transformed to the outside of the
  /// parameter domain
  //mutable int _cntIn;
  //mutable int _cntOut;

  /// Constructor
  SDFBase (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const typename ConfiguratorType3D::InitType &Grid3D,
    const aol::Vector<RealType> &SDF,
    const RangeToWorldCoordTransformation<RealType> & RangeToWorldCoord )
    : _grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdfDiscFuncs ( Grid3D, SDF ),
      _rangeToWorldCoord ( RangeToWorldCoord ) {}

  /// Compute the deformed world coord for a quad point
  /// x[r](\xi)
  bool evaluateWorldCoordDeformed (
    const aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> >
    &DefDiscFuncs,
    const aol::DiscreteFunctionDefault<ConfiguratorType2D> &RangeDiscFuncs,
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::DomVecType &RefCoord,
    typename ConfiguratorType3D::ElementType &WorldCoordDeformedEl,
    typename ConfiguratorType3D::VecType &WorldCoordDeformedLocalCoord ) const {

    // Coordinates of the point in the parameter domain [0,1]^2
    RealType u = El[0] + RefCoord[0];
    RealType v = El[1] + RefCoord[1];

    // Get range data at quad point by bilinear interpolation of the range data
    RealType r = RangeDiscFuncs.evaluateAtQuadPoint ( El, QuadPoint );

    // Compute 3D world coord for this range value
    typename ConfiguratorType3D::VecType WorldCoord =
      _rangeToWorldCoord.evaluateWorldCoord ( u, v, r );

    typename ConfiguratorType3D::ElementType WorldCoordEl;
    typename ConfiguratorType3D::VecType WorldCoordLocalCoord;

    typename ConfiguratorType3D::VecType Deformation;

    for ( int i = 0; i < 3; ++i ) {
      // Get the 3D deformation vector
      Deformation[i] = DefDiscFuncs[i].evaluateAtQuadPoint ( El, QuadPoint );

      // Get the element IDs and local coords of the world coord
      WorldCoordEl[i] = static_cast<short> ( WorldCoord[i] / _grid3D.H() );
      WorldCoordLocalCoord[i] = WorldCoord[i] / _grid3D.H() - WorldCoordEl[i];
    }

    // Compute deformed world coord
    if ( !qc::transformCoord<ConfiguratorType3D> ( _grid3D, WorldCoordEl,
                                                   WorldCoordLocalCoord, Deformation, WorldCoordDeformedEl,
                                                   WorldCoordDeformedLocalCoord ) ) {

      //_cntOut++;
      return false;
    }
    //_cntIn++;
    return true;
  }
};


template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class SDFEvaluator : public SDFBase<ConfiguratorType2D, ConfiguratorType3D> {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  /// The discrete function of the range data
  const aol::DiscreteFunctionDefault<ConfiguratorType2D> _rangeDiscFuncs;

  /// Constructor
  SDFEvaluator (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const typename ConfiguratorType3D::InitType &Grid3D,
    const aol::Vector<RealType> &Sdf,
    const RangeToWorldCoordTransformation<RealType>
    & RangeToWorldCoord,
    const aol::Vector<RealType> &Range )
    : SDFBase<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D, Grid3D, Sdf,
                                                        RangeToWorldCoord ), _rangeDiscFuncs ( Grid2D, Range ) {}

  SDFEvaluator (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const typename ConfiguratorType3D::InitType &Grid3D,
    const aol::Vector<RealType> &Sdf,
    const RangeToWorldCoordTransformation<RealType>
    & RangeToWorldCoord )
    : SDFBase<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D, Grid3D, Sdf,
                                                        RangeToWorldCoord ) {}

  RealType evaluateSqrtDet (
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::DomVecType &RefCoord ) const {

    // gamma
    // Coordinates of the point in the parameter domain [0,1]^2
    // Compute the transformation for this range value
    RealType u = El[0] + RefCoord[0];
    RealType v = El[1] + RefCoord[1];
    aol::Vec3<RealType> gamma = this->_rangeToWorldCoord.evaluateGamma ( u, v );

    // grad gamma
    aol::Vec3<RealType> grad_gamma1;
    aol::Vec3<RealType> grad_gamma2;
    this->_rangeToWorldCoord.evaluateGammaDerivative ( u, v, grad_gamma1, grad_gamma2 );

    // r
    // Get range data at quad point
    RealType r = _rangeDiscFuncs.evaluateAtQuadPoint ( El, QuadPoint );

    // grad r
    // Compute range value gradient in the parameter domain
    aol::Vec2<RealType> grad_r;
    _rangeDiscFuncs.evaluateGradientAtQuadPoint ( El, QuadPoint, grad_r );

    // a1, a2
    aol::Vec3<RealType> a1, a2;
    a1 = grad_gamma1 * r + gamma * grad_r[0];
    a2 = grad_gamma2 * r + gamma * grad_r[1];

    // \sqrt{ \vert det \big( (Dx[r])^T Dx[r] \big) \vert }
    return sqrt ( a1.normSqr() * a2.normSqr() - aol::Sqr ( a1 * a2 ) );
  }
};


template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class SDFEvaluator2 : public SDFBase<ConfiguratorType2D, ConfiguratorType3D> {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  /// Constructor
  SDFEvaluator2 (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const typename ConfiguratorType3D::InitType &Grid3D,
    const aol::Vector<RealType> &Sdf,
    const RangeToWorldCoordTransformation<RealType>
    & RangeToWorldCoord )
    : SDFBase<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D, Grid3D, Sdf,
                                                        RangeToWorldCoord ) {}

  RealType evaluateSqrtDet (
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::DomVecType &RefCoord,
    const aol::DiscreteFunctionDefault<ConfiguratorType2D> RangeDiscFuncs ) const {

    // gamma
    // Coordinates of the point in the parameter domain [0,1]^2
    // Compute the transformation for this range value
    RealType u = El[0] + RefCoord[0];
    RealType v = El[1] + RefCoord[1];
    aol::Vec3<RealType> gamma = this->_rangeToWorldCoord.evaluateGamma ( u, v );

    // grad gamma
    aol::Vec3<RealType> grad_gamma1;
    aol::Vec3<RealType> grad_gamma2;
    this->_rangeToWorldCoord.evaluateGammaDerivative ( u, v, grad_gamma1, grad_gamma2 );

    // r
    // Get range data at quad point
    RealType r = RangeDiscFuncs.evaluateAtQuadPoint ( El, QuadPoint );

    // grad r
    // Compute range value gradient in the parameter domain
    aol::Vec2<RealType> grad_r;
    RangeDiscFuncs.evaluateGradientAtQuadPoint ( El, QuadPoint, grad_r );

    // a1, a2
    aol::Vec3<RealType> a1, a2;
    a1 = grad_gamma1 * r + gamma * grad_r[0];
    a2 = grad_gamma2 * r + gamma * grad_r[1];

    // \sqrt{ \vert det \big( (Dx[r])^T Dx[r] \big) \vert }
    return sqrt ( a1.normSqr() * a2.normSqr() - aol::Sqr ( a1 * a2 ) );
  }
};


//

template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class MatchingEnergyOfPhi : public aol::FENonlinIntegrationVectorInterface < ConfiguratorType2D,
  MatchingEnergyOfPhi<ConfiguratorType2D, ConfiguratorType3D>, 3 > ,
  public SDFEvaluator<ConfiguratorType2D, ConfiguratorType3D> {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  const RealType _alpha;

  MatchingEnergyOfPhi (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const typename ConfiguratorType3D::InitType &Grid3D,
    const aol::Vector<RealType> &Sdf,
    const RangeToWorldCoordTransformation<RealType> & RangeToWorldCoord,
    const aol::Vector<RealType> &Range,
    const RealType Alpha )
    : aol::FENonlinIntegrationVectorInterface < ConfiguratorType2D,
      MatchingEnergyOfPhi<ConfiguratorType2D, ConfiguratorType3D>, 3 > ( Grid2D ),
      SDFEvaluator<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D, Grid3D, Sdf,
                                                             RangeToWorldCoord, Range ), _alpha ( Alpha ) {}

  // d(x[r]+\hat{\phi})
  RealType evaluateIntegrand (
    const aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> >
    &DefDiscFuncs,
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::DomVecType &RefCoord ) const {

    typename ConfiguratorType3D::ElementType WorldCoordDeformedEl;
    typename ConfiguratorType3D::VecType WorldCoordDeformedLocalCoord;

    // compute deformed world coord
    if ( this->evaluateWorldCoordDeformed ( DefDiscFuncs, this->_rangeDiscFuncs, El, QuadPoint,
                                            RefCoord, WorldCoordDeformedEl, WorldCoordDeformedLocalCoord ) ) {
      // sqr(d(x[r]+\hat{\phi}))
      RealType dSqr =  aol::Sqr ( this->_sdfDiscFuncs.evaluate ( WorldCoordDeformedEl,
                                                                 WorldCoordDeformedLocalCoord ) );
      // \sqrt{ \vert det \big( (Dx[r])^T Dx[r] \big) \vert }

      // with area term
      //RealType sqrtDet = evaluateSqrtDet( El, QuadPoint, RefCoord );
      //return 0.5 * ( dSqr + _alpha ) * sqrtDet;
      //return 0.5 * ( dSqr ) * sqrtDet;

      // without area term
      return 0.5 * ( dSqr );
    } else {
      return 0;
    }
  }
};


template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class VariationOfMatchingEnergyWRTPhi : public aol::FENonlinVectorOpInterface< ConfiguratorType2D, 3, 3, VariationOfMatchingEnergyWRTPhi<ConfiguratorType2D, ConfiguratorType3D> > ,
  public SDFEvaluator<ConfiguratorType2D, ConfiguratorType3D> {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  VariationOfMatchingEnergyWRTPhi (  const typename ConfiguratorType2D::InitType &Grid2D,
                                     const typename ConfiguratorType3D::InitType &Grid3D,
                                     const aol::Vector<RealType> &Sdf,
                                     const RangeToWorldCoordTransformation<RealType> & RangeToWorldCoord,
                                     const aol::Vector<RealType> &Range )
    : aol::FENonlinVectorOpInterface< ConfiguratorType2D, 3, 3, VariationOfMatchingEnergyWRTPhi<ConfiguratorType2D, ConfiguratorType3D> > ( Grid2D ),
      SDFEvaluator<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D,
                                                             Grid3D, Sdf, RangeToWorldCoord, Range ) {}

  void getNonlinearity (
    aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> >
    &DefDiscFuncs,
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::VecType &RefCoord,
    aol::Vec<3 , typename ConfiguratorType2D::RealType> &NL ) const {

    // deformed WorldCoord
    typename ConfiguratorType3D::ElementType WorldCoordDeformedEl;
    typename ConfiguratorType3D::VecType WorldCoordDeformedLocalCoord;

    if ( this->evaluateWorldCoordDeformed ( DefDiscFuncs, this->_rangeDiscFuncs, El,
                                            QuadPoint, RefCoord, WorldCoordDeformedEl, WorldCoordDeformedLocalCoord ) ) {
      // d(x[r]+\hat{\phi})
      RealType d = this->_sdfDiscFuncs.evaluate ( WorldCoordDeformedEl,
                                                  WorldCoordDeformedLocalCoord );
      // NL = \nabla \big( d(x[r]+\hat{\phi}) \big)
      this->_sdfDiscFuncs.evaluateGradient ( WorldCoordDeformedEl,
                                             WorldCoordDeformedLocalCoord, NL );
      // \sqrt{ \vert det \big( (Dx[r])^T Dx[r] \big) \vert }

      // with area term
      /*      RealType sqrtDet = evaluateSqrtDet( El, QuadPoint, RefCoord );
            NL *= d * sqrtDet;    */

      // without area term
      NL *= d;
    } else {
      NL.setZero();
    }
  }
};


template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class MatchingEnergyOfRange : public aol::Op< aol::Vector<typename ConfiguratorType2D::RealType>, aol::Scalar<typename ConfiguratorType2D::RealType> >
    , public SDFBase<ConfiguratorType2D, ConfiguratorType3D> {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  const aol::MultiVector<RealType> _deformation;
  const RealType _alpha;

  MatchingEnergyOfRange (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const typename ConfiguratorType3D::InitType &Grid3D,
    const aol::Vector<RealType> &Sdf,
    const RangeToWorldCoordTransformation<RealType>
    & RangeToWorldCoord,
    const aol::MultiVector<RealType> &Deformation,
    const RealType Alpha )
    : SDFBase<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D, Grid3D, Sdf,
                                                        RangeToWorldCoord ), _deformation ( Deformation ), _alpha ( Alpha ) {}

  void applyAdd ( const aol::Vector<RealType> &Range,
                  aol::Scalar<RealType> &Energy ) const {

    MatchingEnergyOfPhi <ConfiguratorType2D, ConfiguratorType3D> E ( this->_grid2D,
                                                                     this->_grid3D, this->_sdfDiscFuncs.getDofs(),
                                                                     this->_rangeToWorldCoord, Range, _alpha );

    E.applyAdd ( _deformation, Energy );
  }
};


template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class FirstPartOfVariationOfMatchingEnergyWRTRange : public aol::FENonlinOpInterface< ConfiguratorType2D, FirstPartOfVariationOfMatchingEnergyWRTRange<ConfiguratorType2D, ConfiguratorType3D> >
    , public SDFEvaluator2<ConfiguratorType2D, ConfiguratorType3D> {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> >
  _defDiscFunc;

  FirstPartOfVariationOfMatchingEnergyWRTRange (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const typename ConfiguratorType3D::InitType &Grid3D,
    const aol::Vector<RealType> &Sdf,
    const RangeToWorldCoordTransformation<RealType>
    & RangeToWorldCoord,
    const aol::MultiVector<RealType> &Deformation )
    : aol::FENonlinOpInterface< ConfiguratorType2D, FirstPartOfVariationOfMatchingEnergyWRTRange<ConfiguratorType2D, ConfiguratorType3D> > ( Grid2D ),
      SDFEvaluator2<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D,
                                                              Grid3D, Sdf, RangeToWorldCoord ) {

    for ( int c = 0; c < 3; c++ ) {
      _defDiscFunc.set_copy ( c,
                              aol::DiscreteFunctionDefault<ConfiguratorType2D> ( this->_grid2D,
                                                                                 Deformation[c] ) );
    }
  }

  void getNonlinearity (
    const aol::DiscreteFunctionDefault<ConfiguratorType2D> &RangeDiscFuncs,
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::VecType &RefCoord,
    typename ConfiguratorType2D::RealType &NL ) const {

    // deformed world coord
    typename ConfiguratorType3D::ElementType WorldCoordDeformedEl;
    typename ConfiguratorType3D::VecType WorldCoordDeformedLocalCoord;

    if ( this->evaluateWorldCoordDeformed ( _defDiscFunc, RangeDiscFuncs, El, QuadPoint,
                                            RefCoord, WorldCoordDeformedEl, WorldCoordDeformedLocalCoord ) ) {
      // d(x[r]+\hat{\phi})
      RealType d =  this->_sdfDiscFuncs.evaluate ( WorldCoordDeformedEl,
                                                   WorldCoordDeformedLocalCoord );

      // \nabla \big( d(x[r]+\hat{\phi}) \big)
      aol::Vec3<RealType> gradD;
      this->_sdfDiscFuncs.evaluateGradient ( WorldCoordDeformedEl,
                                             WorldCoordDeformedLocalCoord, gradD );

      // \gamma
      RealType u = El[0] + RefCoord[0];
      RealType v = El[1] + RefCoord[1];
      aol::Vec3<RealType> gamma = this->_rangeToWorldCoord.evaluateGamma ( u, v );

      // for now: // \sqrt{ \vert det \big( (Dx[r])^T Dx[r] \big) \vert } = 1
      RealType e = gradD * gamma;
      e *= d;
      NL = e;

      // new: \sqrt{ \vert det \big( (Dx[r])^T Dx[r] \big) \vert }
      //RealType sqrtDet = evaluateSqrtDet( El, QuadPoint, RefCoord, RangeDiscFuncs );
      //NL *= sqrtDet;
    } else {
      NL = 0.;
    }
  }
};


template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class SecondPartOfVariationOfMatchingEnergyWRTRange : public aol::FENonlinOpInterface< ConfiguratorType2D, SecondPartOfVariationOfMatchingEnergyWRTRange<ConfiguratorType2D, ConfiguratorType3D> >
    , public SDFEvaluator2<ConfiguratorType2D, ConfiguratorType3D> {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  const RealType _alpha;

  aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> >
  _defDiscFunc;

  SecondPartOfVariationOfMatchingEnergyWRTRange ( const typename ConfiguratorType2D::InitType &Grid2D,
                                                  const typename ConfiguratorType3D::InitType &Grid3D,
                                                  const aol::Vector<RealType> &Sdf,
                                                  const RangeToWorldCoordTransformation<RealType>
                                                  & RangeToWorldCoord,
                                                  const aol::MultiVector<RealType> &Deformation,
                                                  const RealType Alpha )
    : aol::FENonlinOpInterface< ConfiguratorType2D, SecondPartOfVariationOfMatchingEnergyWRTRange<ConfiguratorType2D, ConfiguratorType3D> > ( Grid2D ),
      SDFEvaluator2<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D,
                                                              Grid3D, Sdf, RangeToWorldCoord ), _alpha ( Alpha ) {

    for ( int c = 0; c < 3; c++ ) {
      _defDiscFunc.set_copy ( c,
                              aol::DiscreteFunctionDefault<ConfiguratorType2D> ( this->_grid2D,
                                                                                 Deformation[c] ) );
    }
  }

  void getNonlinearity (
    const aol::DiscreteFunctionDefault<ConfiguratorType2D> &RangeDiscFuncs,
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::VecType &RefCoord,
    typename ConfiguratorType2D::RealType &NL ) const {

    // deformed world coord
    typename ConfiguratorType3D::ElementType WorldCoordDeformedEl;
    typename ConfiguratorType3D::VecType WorldCoordDeformedLocalCoord;

    if ( evaluateWorldCoordDeformed ( _defDiscFunc, RangeDiscFuncs, El, QuadPoint,
                                      RefCoord, WorldCoordDeformedEl, WorldCoordDeformedLocalCoord ) ) {
      // gamma
      RealType u = El[0] + RefCoord[0];
      RealType v = El[1] + RefCoord[1];
      aol::Vec3<RealType> gamma = this->_rangeToWorldCoord.evaluateGamma ( u, v );

      // grad gamma
      aol::Vec3<RealType> grad_gamma1;
      aol::Vec3<RealType> grad_gamma2;
      this->_rangeToWorldCoord.evaluateGammaDerivative ( u, v, grad_gamma1, grad_gamma2 );

      // r
      RealType r = RangeDiscFuncs.evaluateAtQuadPoint ( El, QuadPoint );

      // grad r
      aol::Vec2<RealType> grad_r;
      RangeDiscFuncs.evaluateGradientAtQuadPoint ( El, QuadPoint, grad_r );

      // a1, a2
      aol::Vec3<RealType> a1, a2;
      a1 = grad_gamma1 * r + gamma * grad_r[0];
      a2 = grad_gamma2 * r + gamma * grad_r[1];

      // h1, h2
      aol::Vec3<RealType> h1, h2;
      h1 = a2.normSqr() * a1 - ( a1 * a2 ) * a2;
      h2 = a1.normSqr() * a2 - ( a1 * a2 ) * a1;

      // vartheta
      RealType vartheta = 2 * ( grad_gamma1 * h1 + grad_gamma2 * h2 );

      RealType sqrtDet = evaluateSqrtDet ( El, QuadPoint, RefCoord, RangeDiscFuncs );

      RealType dSqr =  aol::Sqr ( this->_sdfDiscFuncs.evaluate ( WorldCoordDeformedEl,
                                                                 WorldCoordDeformedLocalCoord ) );

      NL = 0.5 * ( dSqr + _alpha ) / ( 2 * sqrtDet ) * vartheta;
    } else {
      NL = 0.;
    }
  }
};


template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class ThirdPartOfVariationOfMatchingEnergyWRTRange : public aol::FENonlinDiffOpInterface< ConfiguratorType2D, ThirdPartOfVariationOfMatchingEnergyWRTRange<ConfiguratorType2D, ConfiguratorType3D> >
    , public SDFEvaluator2<ConfiguratorType2D, ConfiguratorType3D> {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  const RealType _alpha;

  aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> >
  _defDiscFunc;

  ThirdPartOfVariationOfMatchingEnergyWRTRange ( const typename ConfiguratorType2D::InitType &Grid2D,
                                                 const typename ConfiguratorType3D::InitType &Grid3D,
                                                 const aol::Vector<RealType> &Sdf,
                                                 const RangeToWorldCoordTransformation<RealType>
                                                 & RangeToWorldCoord,
                                                 const aol::MultiVector<RealType> &Deformation,
                                                 const RealType Alpha )
    : aol::FENonlinDiffOpInterface< ConfiguratorType2D, ThirdPartOfVariationOfMatchingEnergyWRTRange<ConfiguratorType2D, ConfiguratorType3D> > ( Grid2D ),
      SDFEvaluator2<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D,
                                                              Grid3D, Sdf, RangeToWorldCoord ), _alpha ( Alpha ) {

    for ( int c = 0; c < 3; c++ ) {
      _defDiscFunc.set_copy ( c,
                              aol::DiscreteFunctionDefault<ConfiguratorType2D> ( this->_grid2D,
                                                                                 Deformation[c] ) );
    }
  }

  void getNonlinearity (
    const aol::DiscreteFunctionDefault<ConfiguratorType2D> &RangeDiscFuncs,
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::VecType &RefCoord,
    aol::Vec<ConfiguratorType2D::Dim, typename ConfiguratorType2D::RealType> &NL ) const {

    // deformed world coord
    typename ConfiguratorType3D::ElementType WorldCoordDeformedEl;
    typename ConfiguratorType3D::VecType WorldCoordDeformedLocalCoord;

    if ( evaluateWorldCoordDeformed ( _defDiscFunc, RangeDiscFuncs, El, QuadPoint,
                                      RefCoord, WorldCoordDeformedEl, WorldCoordDeformedLocalCoord ) ) {
      // gamma
      RealType u = El[0] + RefCoord[0];
      RealType v = El[1] + RefCoord[1];
      aol::Vec3<RealType> gamma = this->_rangeToWorldCoord.evaluateGamma ( u, v );

      // grad gamma
      aol::Vec3<RealType> grad_gamma1;
      aol::Vec3<RealType> grad_gamma2;
      this->_rangeToWorldCoord.evaluateGammaDerivative ( u, v, grad_gamma1, grad_gamma2 );

      // r
      RealType r = RangeDiscFuncs.evaluateAtQuadPoint ( El, QuadPoint );

      // grad r
      aol::Vec2<RealType> grad_r;
      RangeDiscFuncs.evaluateGradientAtQuadPoint ( El, QuadPoint, grad_r );

      // a1, a2
      aol::Vec3<RealType> a1, a2;
      a1 = grad_gamma1 * r + gamma * grad_r[0];
      a2 = grad_gamma2 * r + gamma * grad_r[1];

      // h1, h2
      aol::Vec3<RealType> h1, h2;
      h1 = a2.normSqr() * a1 - ( a1 * a2 ) * a2;
      h2 = a1.normSqr() * a2 - ( a1 * a2 ) * a1;

      // grad vartheta
      RealType grad_vartheta1 = 2 * ( gamma * h1 );
      RealType grad_vartheta2 = 2 * ( gamma * h2 );

      RealType sqrtDet = evaluateSqrtDet ( El, QuadPoint, RefCoord, RangeDiscFuncs );

      RealType dSqr =  aol::Sqr ( this->_sdfDiscFuncs.evaluate ( WorldCoordDeformedEl,
                                                                 WorldCoordDeformedLocalCoord ) );

      NL[0] = 0.5 * ( dSqr + _alpha ) / ( 2 * sqrtDet ) * grad_vartheta1;
      NL[1] = 0.5 * ( dSqr + _alpha ) / ( 2 * sqrtDet ) * grad_vartheta2;
    } else {
      NL.setZero();
    }
  }
};


template <typename VecType>
class DirichletEnergyOp : public aol::Op<VecType, aol::Scalar<typename VecType::DataType> > {
protected:
  const aol::Op<VecType, VecType> &_vectorValuedOp;
public:
  // Constructor
  // VectorValuedOp: L
  DirichletEnergyOp ( const aol::Op<VecType, VecType> &VectorValuedOp )
    : _vectorValuedOp ( VectorValuedOp ) {}

  // Arg: R
  virtual void applyAdd ( const VecType &Arg,
                          aol::Scalar<typename VecType::DataType> &Dest ) const {

    VecType temp ( Arg, aol::STRUCT_COPY );
    // LR
    _vectorValuedOp.apply ( Arg, temp );
    // 0.5 LR \cdot R
    Dest[0] += 0.5 * ( temp * Arg );
  }
};


template <typename ConfiguratorType2D>
class L1DifferenceEnergy : public aol::FENonlinIntegrationScalarInterface<ConfiguratorType2D, L1DifferenceEnergy<ConfiguratorType2D> > {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  /// The discrete functions of the measured range data
  const aol::DiscreteFunctionDefault<ConfiguratorType2D> _rangeDiscFunc;
  /// The converter from range data to 3D world coordinates
  const RangeToWorldCoordTransformation<RealType> & _rangeToWorldCoord;

  const RealType _delta;

  /// Constructor
  L1DifferenceEnergy ( const typename ConfiguratorType2D::InitType &Grid2D,
                       const aol::Vector<RealType> &Range,
                       const RangeToWorldCoordTransformation<RealType>
                       & RangeToWorldCoord,
                       const RealType Delta )
    : aol::FENonlinIntegrationScalarInterface< ConfiguratorType2D, L1DifferenceEnergy<ConfiguratorType2D> > ( Grid2D ),
      _rangeDiscFunc ( Grid2D, Range ), _rangeToWorldCoord ( RangeToWorldCoord ), _delta ( Delta ) {}

  RealType evaluateIntegrand (
    const aol::DiscreteFunctionDefault<ConfiguratorType2D> &RangeDiscFuncs,
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::DomVecType &/*RefCoord*/ ) const {

    // the denoised range value
    RealType r = RangeDiscFuncs.evaluateAtQuadPoint ( El, QuadPoint );
    // the measured range value
    RealType rMeasured = _rangeDiscFunc.evaluateAtQuadPoint ( El, QuadPoint );

    // \vert x \vert _{\delta}
    return std::sqrt ( aol::Sqr ( r - rMeasured ) + aol::Sqr ( _delta ) );
  }
};


template <typename ConfiguratorType2D>
class VariationOfL1DifferenceEnergy : public aol::FENonlinOpInterface< ConfiguratorType2D, VariationOfL1DifferenceEnergy<ConfiguratorType2D> > {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  /// The discrete functions of the measured range data
  const aol::DiscreteFunctionDefault<ConfiguratorType2D> _rangeDiscFunc;
  /// The converter from range data to 3D world coordinates
  const RangeToWorldCoordTransformation<RealType> & _rangeToWorldCoord;

  const RealType _delta;

  VariationOfL1DifferenceEnergy (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const aol::Vector<RealType> &Range,
    const RangeToWorldCoordTransformation<RealType> & RangeToWorldCoord,
    const RealType Delta )
    : aol::FENonlinOpInterface< ConfiguratorType2D, VariationOfL1DifferenceEnergy<ConfiguratorType2D> > ( Grid2D ),
      _rangeDiscFunc ( Grid2D, Range ), _rangeToWorldCoord ( RangeToWorldCoord ), _delta ( Delta ) {}

  void getNonlinearity (
    const aol::DiscreteFunctionDefault<ConfiguratorType2D> &RangeDiscFuncs,
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::VecType &/*RefCoord*/,
    typename ConfiguratorType2D::RealType &NL ) const {

    // the denoised range value
    RealType r = RangeDiscFuncs.evaluateAtQuadPoint ( El, QuadPoint );
    // the measured range value
    RealType rMeasured = _rangeDiscFunc.evaluateAtQuadPoint ( El, QuadPoint );

    // \vert x \vert _{\delta}
    NL = ( r - rMeasured ) /
         ( std::sqrt ( aol::Sqr ( r - rMeasured ) + aol::Sqr ( _delta ) ) );
  }
};


template <typename ConfiguratorType2D>
class L2DifferenceEnergy : public aol::FENonlinIntegrationScalarInterface<ConfiguratorType2D, L2DifferenceEnergy<ConfiguratorType2D> > {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  /// The discrete functions of the measured range data
  const aol::DiscreteFunctionDefault<ConfiguratorType2D> _rangeDiscFunc;
  /// The converter from range data to 3D world coordinates
  const RangeToWorldCoordTransformation<RealType> & _rangeToWorldCoord;

  const RealType _delta;

  /// Constructor
  L2DifferenceEnergy ( const typename ConfiguratorType2D::InitType &Grid2D,
                       const aol::Vector<RealType> &Range,
                       const RangeToWorldCoordTransformation<RealType>
                       & RangeToWorldCoord,
                       const RealType Delta )
    : aol::FENonlinIntegrationScalarInterface< ConfiguratorType2D, L2DifferenceEnergy<ConfiguratorType2D> > ( Grid2D ),
      _rangeDiscFunc ( Grid2D, Range ), _rangeToWorldCoord ( RangeToWorldCoord ), _delta ( Delta ) {}

  RealType evaluateIntegrand (
    const aol::DiscreteFunctionDefault<ConfiguratorType2D> &RangeDiscFuncs,
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::DomVecType &/*RefCoord*/ ) const {

    // the denoised range value
    RealType r = RangeDiscFuncs.evaluateAtQuadPoint ( El, QuadPoint );
    // the measured range value
    RealType rMeasured = _rangeDiscFunc.evaluateAtQuadPoint ( El, QuadPoint );

    // \vert x \vert _{\delta}
    return aol::Sqr ( r - rMeasured );
  }
};



template <typename ConfiguratorType2D>
class VariationOfL2DifferenceEnergy : public aol::FENonlinOpInterface< ConfiguratorType2D, VariationOfL2DifferenceEnergy<ConfiguratorType2D> > {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  /// The discrete functions of the measured range data
  const aol::DiscreteFunctionDefault<ConfiguratorType2D> _rangeDiscFunc;
  /// The converter from range data to 3D world coordinates
  const RangeToWorldCoordTransformation<RealType> & _rangeToWorldCoord;

  const RealType _delta;

  VariationOfL2DifferenceEnergy (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const aol::Vector<RealType> &Range,
    const RangeToWorldCoordTransformation<RealType> & RangeToWorldCoord,
    const RealType Delta )
    : aol::FENonlinOpInterface< ConfiguratorType2D, VariationOfL2DifferenceEnergy<ConfiguratorType2D> > ( Grid2D ),
      _rangeDiscFunc ( Grid2D, Range ), _rangeToWorldCoord ( RangeToWorldCoord ), _delta ( Delta ) {}

  void getNonlinearity (
    const aol::DiscreteFunctionDefault<ConfiguratorType2D> &RangeDiscFuncs,
    const typename ConfiguratorType2D::ElementType &El,
    int QuadPoint,
    const typename ConfiguratorType2D::VecType &/*RefCoord*/,
    typename ConfiguratorType2D::RealType &NL ) const {

    // the denoised range value
    RealType r = RangeDiscFuncs.evaluateAtQuadPoint ( El, QuadPoint );
    // the measured range value
    RealType rMeasured = _rangeDiscFunc.evaluateAtQuadPoint ( El, QuadPoint );

    // \vert x \vert _{\delta}
    NL = 2 * ( r - rMeasured );
  }
};


#endif //__SURFACEREGISTRATION_H






























