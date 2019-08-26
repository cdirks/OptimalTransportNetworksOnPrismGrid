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


/**
 * \author Bauer
 */

template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class SDFBase {
public:
  typedef typename ConfiguratorType2D::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdf;
  const aol::DiscreteFunctionDefault<ConfiguratorType2D> _range;
  const RangeToWorldCoordTransformation<RealType> & _rangeToWorldCoord;
  const typename ConfiguratorType3D::InitType & _grid3D;
  mutable int _cntIn;
  mutable int _cntOut;

  SDFBase ( const typename ConfiguratorType2D::InitType &Grid2D,
            const typename ConfiguratorType3D::InitType &Grid3D,
            const aol::Vector<RealType> &Sdf,
            const aol::Vector<RealType> &Range,
            const RangeToWorldCoordTransformation<RealType> & RangeToWorldCoord )
    : _sdf ( Grid3D, Sdf ), _range ( Grid2D, Range ), _rangeToWorldCoord ( RangeToWorldCoord ), _grid3D ( Grid3D ) { _cntIn = 0; _cntOut = 0; }

  // g[r](u) +\hat{\phi}(u)
  bool evaluateGPlusPhiHat ( const aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> > &DiscFuncs,
                             const typename ConfiguratorType2D::ElementType &El,
                             int QuadPoint, const typename ConfiguratorType2D::DomVecType &RefCoord,
                             typename ConfiguratorType3D::ElementType &GPlusPhiHatEl,
                             typename ConfiguratorType3D::VecType &GPlusPhiHatLocalCoord ) const {

    RealType r = _range.evaluateAtQuadPoint ( El, QuadPoint );

    // sphere_full
    //RType shiftOffset = 2.65;
    // sphere_0.33
    RealType shiftOffset = 2.878;

    // g(u,v)
    RealType u = El[0] + RefCoord[0];
    RealType v = El[1] + RefCoord[1];
    aol::Vec3<RealType> abc = _rangeToWorldCoord.evaluate ( u, v );
    typename ConfiguratorType3D::VecType g = abc * r;
    g[2] -= shiftOffset;
    typename ConfiguratorType3D::VecType phihat;

    typename ConfiguratorType3D::VecType gLocalCoord;
    typename ConfiguratorType3D::ElementType gEl;

    //std::cout << g[0] << " " << g[1] << " " << g[2] << std::endl;

    for ( int i = 0; i < 3; ++i ) {
      phihat[i] = DiscFuncs[i].evaluateAtQuadPoint ( El, QuadPoint );

      //gEl[i] = static_cast<short> ( g[i] );
      //gLocalCoord[i] = g[i] - gEl[i];

      gEl[i] = static_cast<short> ( g[i] / _grid3D.H() );
      gLocalCoord[i] = g[i] / _grid3D.H() - gEl[i];
      //std::cout << gEl[i] << std::endl;
      //std::cout << gLocalCoord[i] << std::endl;

    }

    if ( !qc::transformCoord<ConfiguratorType3D> ( _grid3D, gEl, gLocalCoord, phihat, GPlusPhiHatEl, GPlusPhiHatLocalCoord ) ) {

      _cntOut++;
      return false;
    }
    _cntIn++;
    return true;
  }

  // d(x), x \in \mathbb{R}^3
  RealType evaluateD ( const aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> > &DiscFuncs,
                       const typename ConfiguratorType3D::ElementType &gPlusPhiHatEl,
                       const typename ConfiguratorType3D::DomVecType &gPlusPhiHatLocalCoord ) const {

    return ( _sdf.evaluate ( gPlusPhiHatEl, gPlusPhiHatLocalCoord ) );
  }

  // \nabla d(x), x \in \mathbb{R}^3
  void evaluateGradD ( const aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> > &DiscFuncs,
                       const typename ConfiguratorType3D::ElementType &gPlusPhiHatEl,
                       const typename ConfiguratorType3D::DomVecType &gPlusPhiHatLocalCoord,
                       aol::Vec<3 , typename ConfiguratorType2D::RealType> &NL ) const {

    _sdf.evaluateGradient ( gPlusPhiHatEl, gPlusPhiHatLocalCoord, NL );
  }

  // \vert det Dg[r](u) \vert
  RealType evaluateDetDg (   const aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> > &DiscFuncs,
                             const typename ConfiguratorType2D::ElementType &El,
                             int QuadPoint, const typename ConfiguratorType2D::DomVecType &RefCoord ) const {

    RealType r = _range.evaluateAtQuadPoint ( El, QuadPoint );

    RealType u = El[0] + RefCoord[0];
    RealType v = El[1] + RefCoord[1];

    aol::Vec3<RealType> du;
    aol::Vec3<RealType> dv;

    _rangeToWorldCoord.evaluateDerivative ( u, v, du, dv );

    aol::Vec2<RealType> gradR;
    _range.evaluateGradientAtQuadPoint ( El, QuadPoint, gradR );

    aol::Vec3<RealType> abc = _rangeToWorldCoord.evaluate ( u, v );

    du *= r;
    du.addMultiple ( abc, gradR[0] );

    dv *= r;
    dv.addMultiple ( abc, gradR[1] );

    return sqrt ( du.normSqr() + dv.normSqr() - du * dv );
  }
};




/**
 * \brief  Signed Distance Function (SDF) Energy
 *
 * \author Bauer
 */
template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class SDFEnergy : public aol::FENonlinIntegrationVectorInterface<ConfiguratorType2D, SDFEnergy<ConfiguratorType2D, ConfiguratorType3D>, 3 >, public SDFBase<ConfiguratorType2D, ConfiguratorType3D> {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  SDFEnergy ( const typename ConfiguratorType2D::InitType &Grid2D,
              const typename ConfiguratorType3D::InitType &Grid3D,
              const aol::Vector<RealType> &Sdf,
              const aol::Vector<RealType> &Range,
              const RangeToWorldCoordTransformation<RealType> & RangeToWorldCoord )
    : aol::FENonlinIntegrationVectorInterface< ConfiguratorType2D, SDFEnergy<ConfiguratorType2D, ConfiguratorType3D>, 3 > ( Grid2D ),
      SDFBase<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D, Grid3D, Sdf, Range, RangeToWorldCoord ) {}

  RealType evaluateIntegrand ( const aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> > &DiscFuncs,
                               const typename ConfiguratorType2D::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType2D::DomVecType &RefCoord ) const {

    typename ConfiguratorType3D::ElementType gPlusPhiHatEl;
    typename ConfiguratorType3D::VecType gPlusPhiHatLocalCoord;

    if ( evaluateGPlusPhiHat ( DiscFuncs, El, QuadPoint, RefCoord, gPlusPhiHatEl, gPlusPhiHatLocalCoord ) ) {
      RealType dSqr =  aol::Sqr ( evaluateD ( DiscFuncs, gPlusPhiHatEl, gPlusPhiHatLocalCoord ) );
      RealType DetDg = evaluateDetDg ( DiscFuncs, El, QuadPoint, RefCoord );

      return 0.5 * dSqr * DetDg;
    } else {
      return 0;
    }
  }
};

template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class SDFForce : public aol::FENonlinVectorOpInterface< ConfiguratorType2D, 3, 3, SDFForce<ConfiguratorType2D, ConfiguratorType3D> > , public SDFBase<ConfiguratorType2D, ConfiguratorType3D> {
public:
  typedef typename ConfiguratorType2D::RealType RealType;

  SDFForce (  const typename ConfiguratorType2D::InitType &Grid2D,
              const typename ConfiguratorType3D::InitType &Grid3D,
              const aol::Vector<RealType> &Sdf,
              const aol::Vector<RealType> &Range,
              const RangeToWorldCoordTransformation<RealType> & RangeToWorldCoord )
    : aol::FENonlinVectorOpInterface< ConfiguratorType2D, 3, 3, SDFForce<ConfiguratorType2D, ConfiguratorType3D> > ( Grid2D ),
      SDFBase<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D, Grid3D, Sdf, Range, RangeToWorldCoord ) {}

  void getNonlinearity ( aol::auto_container<3, aol::DiscreteFunctionDefault<ConfiguratorType2D> > &DiscFuncs,
                         const typename ConfiguratorType2D::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType2D::VecType &RefCoord,
                         aol::Vec<3 , typename ConfiguratorType2D::RealType> &NL ) const {

    typename ConfiguratorType3D::ElementType gPlusPhiHatEl;
    typename ConfiguratorType3D::VecType gPlusPhiHatLocalCoord;

    if ( evaluateGPlusPhiHat ( DiscFuncs, El, QuadPoint, RefCoord, gPlusPhiHatEl, gPlusPhiHatLocalCoord ) ) {
      RealType d =  evaluateD ( DiscFuncs, gPlusPhiHatEl, gPlusPhiHatLocalCoord );
      RealType DetDg = evaluateDetDg ( DiscFuncs, El, QuadPoint, RefCoord );
      evaluateGradD ( DiscFuncs, gPlusPhiHatEl, gPlusPhiHatLocalCoord, NL );

      NL *= d * DetDg;
    } else {
      NL.setZero();
    }

  }
};




#endif //__SURFACEREGISTRATION_H
