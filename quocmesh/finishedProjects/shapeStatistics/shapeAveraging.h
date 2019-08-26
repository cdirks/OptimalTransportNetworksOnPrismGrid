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

#ifndef __SHAPEAVERAGING_H
#define __SHAPEAVERAGING_H

#include "testResultSaving.h"
#include "auxiliaryOps.h"
#include "stressOps.h"
#include <aol.h>
#include <hyperelastic.h>
#include <parameterParser.h>
#include <signedDistanceOp.h>
#include <eigenvectors.h>
#include <gradientDescent.h>
#include <Newton.h>
#include <deformations.h>
#include <morphology.h>
#include <simplexConfigurators.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#include <gradientflow.h>

// for resolveSelfPenetration:
#include "../../finishedProjects/geodesicCalculus/geodesicCalculusEnergies.h"

// bool TEST = false;
// bool TEST1 = false;
// int TEST_INDEX = 1000;


/**************************************************************
 * The averaging energy of one single deformation and its variation.
 **************************************************************/

enum MismatchPenaltyType {
  SquaredDifference,
  PhasefieldMismatch
};

template<typename ConfiguratorType>
class StressMatch :
  public aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim*ConfiguratorType::Dim, ConfiguratorType::Dim, StressMatch<ConfiguratorType> > {

public:
  StressMatch( const typename ConfiguratorType::InitType &Grid ) :
    aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim*ConfiguratorType::Dim, ConfiguratorType::Dim, StressMatch<ConfiguratorType> >( Grid ) {}

  void getNonlinearity ( const aol::auto_container<ConfiguratorType::Dim*ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        NL[i][j] = DiscFuncs[i*ConfiguratorType::Dim+j].evaluateAtQuadPoint( El, QuadPoint );
  }
};

template<typename ConfiguratorType>
class PreStressedStressMatch :
  public aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim*ConfiguratorType::Dim, ConfiguratorType::Dim, PreStressedStressMatch<ConfiguratorType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  // the prestressing displacement
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _preStressDisp;
  // the inverse prestress displacement, $\phi^{-1}$-identity, in all Cartesian directions as vector of the dofs and as multilinear interpolation
  aol::MultiVector<RealType> _invDispDOFs;
  const aol::DiscreteVectorFunctionDefault< ConfiguratorType, ConfiguratorType::Dim > _invDisp;
  // array encoding whether $\phi^{-1}$ is defined at the corresponding node
  qc::BitArray<ConfiguratorType::Dim> _invPhiDefined;

public:
  PreStressedStressMatch( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<typename ConfiguratorType::RealType> &PreStressDisp ) :
    aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim*ConfiguratorType::Dim, ConfiguratorType::Dim, PreStressedStressMatch<ConfiguratorType> >( Grid ),
    _preStressDisp( Grid, PreStressDisp ),
    _invDispDOFs( PreStressDisp, aol::STRUCT_COPY ),
    _invDisp( Grid, _invDispDOFs ),
    _invPhiDefined( Grid ) {
    // generate the identity function
    aol::MultiVector<RealType> identity( PreStressDisp, aol::STRUCT_COPY );
    qc::DataGenerator<ConfiguratorType>( Grid ).generateIdentity( identity );
    // compute $\phi^{-1}$-identity
    qc::TransformFunction<RealType, ConfiguratorType::Dim> transformOp( Grid );
    transformOp.setDeformation( PreStressDisp );
    transformOp.transform( identity, _invDispDOFs, _invPhiDefined );
    _invDispDOFs -= identity;
  }

  void getNonlinearity ( const aol::auto_container<ConfiguratorType::Dim*ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                         aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {
    NL.setZero();

    // if the inverse prestress deformation is not defined at this position, return 0
    if ( !_invPhiDefined.elementTrue( El ) )
      return;

    // compute the position before prestress deformation
    typename ConfiguratorType::VecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    if ( !qc::transformCoord<ConfiguratorType> ( this->_initializer, _invDisp, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord ) )
      return;

    // compute stress at the position before prestressing
    typename ConfiguratorType::MatType stress;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        stress[i][j] = DiscFuncs[i*ConfiguratorType::Dim+j].evaluate( transformedEl, transformedLocalCoord );

    // compute the deformation gradient of the prestressing deformation
    typename ConfiguratorType::MatType dphi;
    _preStressDisp.evaluateGradient ( transformedEl, transformedLocalCoord, dphi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dphi[i][i] += 1.;

    NL.makeProductABtransposed ( stress, dphi );
    NL /= aol::Abs( dphi.det() );
  }
};

template <typename ConfiguratorType, typename ElasticEnergyType, MismatchPenaltyType MismatchType>
class EnergyOfD :
  public aol::FENonlinIntegrationVectorInterface< ConfiguratorType, EnergyOfD<ConfiguratorType, ElasticEnergyType, MismatchType>, ConfiguratorType::Dim > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  // the two given shapes/images of which the second will be deformed, both with values {0,1}
  const aol::DiscreteFunctionDefault<ConfiguratorType> _uFixed, _uTemplate;
  // the (hyper-)elastic energy term
  const ElasticEnergyType &_elasticEnergy;
  // the mismatch parameter
  const RealType _gamma_epsilon;
  
  mutable aol::Scalar<RealType> _lastMismatchEnergy, _lastElasticEnergy;

  // if required, the stress for the stress mismatch penalty
  aol::MultiVector<RealType> _stressProd;

private:
  inline RealType evaluateSqrDifference( const RealType A, const RealType B ) const {
    return aol::Sqr( A - B );
  }

  inline RealType evaluatePhasefieldMismatch( const RealType A, const RealType B ) const {
    return aol::Sqr( A * ( 1 - B ) ) + aol::Sqr( ( 1 - A ) * B );
  }

public:
  EnergyOfD( const typename ConfiguratorType::InitType &Grid,
             const aol::Vector<RealType> &UFixed,
             const aol::Vector<RealType> &UTemplate,
             const ElasticEnergyType &ElasticEnergy,
             const RealType Gamma_Epsilon ) :
    aol::FENonlinIntegrationVectorInterface< ConfiguratorType, EnergyOfD<ConfiguratorType, ElasticEnergyType, MismatchType>, ConfiguratorType::Dim >( Grid ),
    _uFixed( Grid, UFixed ),
    _uTemplate( Grid, UTemplate ),
    _elasticEnergy( ElasticEnergy ),
    _gamma_epsilon( Gamma_Epsilon ) {}

  EnergyOfD( const typename ConfiguratorType::InitType &Grid,
             const aol::Vector<RealType> &UFixed,
             const aol::Vector<RealType> &UTemplate,
             const ElasticEnergyType &ElasticEnergy,
             const RealType Gamma_Epsilon,
             const aol::MultiVector<RealType> &Stress ) :
    aol::FENonlinIntegrationVectorInterface< ConfiguratorType, EnergyOfD<ConfiguratorType, ElasticEnergyType, MismatchType>, ConfiguratorType::Dim >( Grid ),
    _uFixed( Grid, UFixed ),
    _uTemplate( Grid, UTemplate ),
    _elasticEnergy( ElasticEnergy ),
    _gamma_epsilon( Gamma_Epsilon ) {
    if ( Stress.numComponents() > 0 ) {
      _stressProd.reallocate( Grid );
      ( StressMatch<ConfiguratorType>( Grid ) ).apply( Stress, _stressProd );
    }
  }

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint,
                              const typename ConfiguratorType::DomVecType &RefCoord) const {
    //! compute \f$\phi(x)\f$
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    qc::transformAndClipCoord<ConfiguratorType> ( *this->_config, DiscFuncs, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord );

    //! compute the value of the deformed shape and the undeformed one
    RealType uTemplate = _uTemplate.evaluate( transformedEl, transformedLocalCoord ),
             uFixed    = _uFixed.evaluateAtQuadPoint( El, QuadPoint );

    // if $x$ was not displaced outside the grid, return result
    switch ( MismatchType ) {
    case SquaredDifference:  return _gamma_epsilon * evaluateSqrDifference( uFixed, uTemplate );
    case PhasefieldMismatch: return _gamma_epsilon * evaluatePhasefieldMismatch( uFixed, uTemplate );
    default:
      throw aol::UnimplementedCodeException ( "Unsupported MismatchType", __FILE__, __LINE__ );
    }
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    // compute the matching term of the energy
    _lastMismatchEnergy.setZero();
    aol::FENonlinIntegrationVectorInterface< ConfiguratorType, EnergyOfD<ConfiguratorType, ElasticEnergyType, MismatchType>, ConfiguratorType::Dim >::applyAdd( Arg, _lastMismatchEnergy );
    Dest += _lastMismatchEnergy;
    
    // compute the hyperelastic energy
    _elasticEnergy.apply( Arg, _lastElasticEnergy );
    Dest += _lastElasticEnergy;

    // if required, compute the stress mismatch penalty
    if ( _stressProd.numComponents() > 0 )
      Dest += _stressProd * Arg;
  }

  RealType getLastMismatchEnergy( bool weighted = true ) const {
    return weighted ? _lastMismatchEnergy[0] : _lastMismatchEnergy[0] / _gamma_epsilon;
  }
  
  RealType getLastElasticEnergy( ) const {
    return _lastElasticEnergy[0];
  }

};

template <typename ConfiguratorType, typename ElasticEnergyType, MismatchPenaltyType MismatchType>
class PreStressedEnergyOfD :
  public aol::FENonlinIntegrationVectorInterface< ConfiguratorType, PreStressedEnergyOfD<ConfiguratorType, ElasticEnergyType, MismatchType>, ConfiguratorType::Dim > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  // the two given shapes/images of which the second will be deformed, both with values {0,1}
  const aol::DiscreteFunctionDefault<ConfiguratorType> _uFixed, _uTemplate;
  // the (hyper-)elastic energy term
  const ElasticEnergyType &_elasticEnergy;
  // the mismatch parameter
  const RealType _gamma_epsilon;
  // the stress for the stress mismatch penalty
  aol::MultiVector<RealType> _stressProd;
  // the prestressing displacement
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _preStressDisp;

private:
  inline RealType evaluateSqrDifference( const RealType A, const RealType B ) const {
    return aol::Sqr( A - B );
  }

  inline RealType evaluatePhasefieldMismatch( const RealType A, const RealType B ) const {
    return aol::Sqr( A * ( 1 - B ) ) + aol::Sqr( ( 1 - A ) * B );
  }

public:
  PreStressedEnergyOfD( const typename ConfiguratorType::InitType &Grid,
                        const aol::Vector<RealType> &UFixed,
                        const aol::Vector<RealType> &UTemplate,
                        const ElasticEnergyType &ElasticEnergy,
                        const RealType Gamma_Epsilon,
                        const aol::MultiVector<RealType> &PreStressDisp,
                        const aol::MultiVector<RealType> &Stress ) :
    aol::FENonlinIntegrationVectorInterface< ConfiguratorType, PreStressedEnergyOfD<ConfiguratorType, ElasticEnergyType, MismatchType>, ConfiguratorType::Dim >( Grid ),
    _uFixed( Grid, UFixed ),
    _uTemplate( Grid, UTemplate ),
    _elasticEnergy( ElasticEnergy ),
    _gamma_epsilon( Gamma_Epsilon ),
    _stressProd( Grid ),
    _preStressDisp( Grid, PreStressDisp ) {
    ( PreStressedStressMatch<ConfiguratorType>( Grid ) ).apply( Stress, _stressProd );
  }

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint,
                              const typename ConfiguratorType::VecType &RefCoord) const {
    //! compute \f$\phi(x)\f$
    typename ConfiguratorType::VecType intermedCoord, transformedLocalCoord, offset;
    typename ConfiguratorType::ElementType intermedEl, transformedEl;
    if( !qc::transformCoord<ConfiguratorType> ( this->_initializer, _preStressDisp, El, QuadPoint, RefCoord, intermedEl, intermedCoord ) )
      // if point $x$ was displaced outside the grid, return 0
      return 0.;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      offset[i] = DiscFuncs[i].evaluate( intermedEl, intermedCoord );
    if( !qc::transformCoord<ConfiguratorType> ( this->_initializer, intermedEl, intermedCoord, offset, transformedEl, transformedLocalCoord ) )
      // if point $x$ was displaced outside the grid, return 0
      return 0.;

    //! compute the value of the deformed shape and the undeformed one
    RealType uTemplate = _uTemplate.evaluate( transformedEl, transformedLocalCoord ),
             uFixed    = _uFixed.evaluateAtQuadPoint( El, QuadPoint );

    // if $x$ was not displaced outside the grid, return result
    switch ( MismatchType ) {
    case SquaredDifference:  return _gamma_epsilon * evaluateSqrDifference( uFixed, uTemplate );
    case PhasefieldMismatch: return _gamma_epsilon * evaluatePhasefieldMismatch( uFixed, uTemplate );
    }
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    // compute the matching term of the energy
    aol::FENonlinIntegrationVectorInterface< ConfiguratorType, EnergyOfD<ConfiguratorType, ElasticEnergyType, MismatchType>, ConfiguratorType::Dim >::applyAdd( Arg, Dest );

    // compute the hyperelastic energy
    _elasticEnergy.applyAdd( Arg, Dest );

    // compute the stress mismatch penalty
    Dest += _stressProd * Arg;
  }
};

template <typename ConfiguratorType, MismatchPenaltyType MismatchType>
class EnergyOfDDefault :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > WeightedEnergyDensityType;
  typedef qc::HyperelasticEnergy<ConfiguratorType,WeightedEnergyDensityType> WeightedHyperelasticEnergyType;

  // the default elastic energy
  const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> _energyDensity;
  const WeightedEnergyDensityType _weightedDensity;
  const WeightedHyperelasticEnergyType _hyperelasticEnergy;

  // the energy term itself
  const EnergyOfD<ConfiguratorType, WeightedHyperelasticEnergyType, MismatchType> _energyOfD;

public:
  EnergyOfDDefault( const typename ConfiguratorType::InitType &Grid,
                    const aol::Vector<RealType> &UFixed,
                    const aol::Vector<RealType> &UTemplate,
                    const RealType LengthEnergyWeight,
                    const RealType VolumeEnergyWeight,
                    const RealType Gamma_Epsilon,
                    const aol::Vector<RealType> &Weight ) :
    _energyDensity( LengthEnergyWeight, 0., VolumeEnergyWeight ),
    _weightedDensity( _energyDensity, Grid, Weight ),
    _hyperelasticEnergy( Grid, _weightedDensity ),
    _energyOfD( Grid, UFixed, UTemplate, _hyperelasticEnergy, Gamma_Epsilon ) {}

  EnergyOfDDefault( const typename ConfiguratorType::InitType &Grid,
                    const aol::Vector<RealType> &UFixed,
                    const aol::Vector<RealType> &UTemplate,
                    const RealType LengthEnergyWeight,
                    const RealType VolumeEnergyWeight,
                    const RealType Gamma_Epsilon,
                    const aol::Vector<RealType> &Weight,
                    const aol::MultiVector<RealType> &Stress ) :
    _energyDensity( LengthEnergyWeight, 0., VolumeEnergyWeight ),
    _weightedDensity( _energyDensity, Grid, Weight ),
    _hyperelasticEnergy( Grid, _weightedDensity ),
    _energyOfD( Grid, UFixed, UTemplate, _hyperelasticEnergy, Gamma_Epsilon, Stress ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    _energyOfD.applyAdd( Arg, Dest );
  }
};

template <typename ConfiguratorType, typename ElasticGradientType, MismatchPenaltyType MismatchType>
class EnergyVariationWRTD :
  public aol::FENonlinVectorOpInterface< ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, EnergyVariationWRTD<ConfiguratorType, ElasticGradientType, MismatchType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  // the two given shapes/images of which the second will be deformed, both with values {0,1}
  const aol::DiscreteFunctionDefault<ConfiguratorType> _uFixed, _uTemplate;
  // the (hyper-)elastic energy term
  const ElasticGradientType &_elasticGradient;
  // the mismatch parameter
  const RealType _gamma_epsilon;

  // if required, the stress for the stress mismatch penalty
  aol::MultiVector<RealType> _stressProd;

private:
  inline RealType evaluateSqrDifferenceDeriv( const RealType A, const RealType B ) const {
    return 2 * ( A - B );
  }

  inline RealType evaluatePhasefieldMismatchDeriv( const RealType A, const RealType B ) const {
    return 2 * ( A * aol::Sqr( 1 - B ) + ( A - 1 ) * aol::Sqr( B ) );
  }

public:
  EnergyVariationWRTD( const typename ConfiguratorType::InitType &Grid,
                       const aol::Vector<RealType> &UFixed,
                       const aol::Vector<RealType> &UTemplate,
                       const ElasticGradientType &ElasticGradient,
                       const RealType Gamma_Epsilon ) :
    aol::FENonlinVectorOpInterface< ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, EnergyVariationWRTD<ConfiguratorType, ElasticGradientType, MismatchType> >( Grid ),
    _uFixed( Grid, UFixed ),
    _uTemplate( Grid, UTemplate ),
    _elasticGradient( ElasticGradient ),
    _gamma_epsilon( Gamma_Epsilon ) {}

  EnergyVariationWRTD( const typename ConfiguratorType::InitType &Grid,
                       const aol::Vector<RealType> &UFixed,
                       const aol::Vector<RealType> &UTemplate,
                       const ElasticGradientType &ElasticGradient,
                       const RealType Gamma_Epsilon,
                       const aol::MultiVector<RealType> &Stress ) :
    aol::FENonlinVectorOpInterface< ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, EnergyVariationWRTD<ConfiguratorType, ElasticGradientType, MismatchType> >( Grid ),
    _uFixed( Grid, UFixed ),
    _uTemplate( Grid, UTemplate ),
    _elasticGradient( ElasticGradient ),
    _gamma_epsilon( Gamma_Epsilon ) {
    if ( Stress.numComponents() > 0 ) {
      _stressProd.reallocate( Grid );
      ( StressMatch<ConfiguratorType>( Grid ) ).apply( Stress, _stressProd );
    }
  }

  void getNonlinearity ( aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                       const typename ConfiguratorType::ElementType &El,
                       int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                       aol::Vec<ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {

    //! compute \f$\phi(x)\f$
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    aol::Vec<ConfiguratorType::Dim,bool> coordinateWithinLimits;
    qc::transformAndClipCoord<ConfiguratorType> ( *this->_config, DiscFuncs, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord, coordinateWithinLimits );

    //! compute the value of the deformed shape and the undeformed one as well as the gradient
    RealType uTemplate = _uTemplate.evaluate( transformedEl, transformedLocalCoord ),
             uFixed    = _uFixed.evaluateAtQuadPoint( El, QuadPoint );
    _uTemplate.evaluateGradient ( transformedEl, transformedLocalCoord, NL );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      if ( coordinateWithinLimits[i] == false )
        NL[i] = 0.;

    //! return the derivative
    switch ( MismatchType ) {
    case SquaredDifference:  NL *= _gamma_epsilon * evaluateSqrDifferenceDeriv( uTemplate, uFixed ); break;
    case PhasefieldMismatch: NL *= _gamma_epsilon * evaluatePhasefieldMismatchDeriv( uTemplate, uFixed ); break;
    default:
      throw aol::UnimplementedCodeException ( "Unsupported MismatchType", __FILE__, __LINE__ );
    }
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // compute the phase field term of the energy variation
    aol::FENonlinVectorOpInterface< ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, EnergyVariationWRTD<ConfiguratorType, ElasticGradientType, MismatchType> >::applyAdd( Arg, Dest );

    // compute the hyperelastic energy variation
    _elasticGradient.applyAdd( Arg, Dest );

    // if required, compute the stress mismatch penalty
    if ( _stressProd.numComponents() > 0 )
      Dest += _stressProd;
  }
};

template <typename ConfiguratorType, typename ElasticGradientType, MismatchPenaltyType MismatchType>
class PreStressedEnergyVariationWRTD :
  public aol::FENonlinVectorOpInterface< ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, PreStressedEnergyVariationWRTD<ConfiguratorType, ElasticGradientType, MismatchType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  // the two given shapes/images of which the second will be deformed, both with values {0,1}
  const aol::DiscreteFunctionDefault<ConfiguratorType> _uFixed, _uTemplate;
  // the (hyper-)elastic energy term
  const ElasticGradientType &_elasticGradient;
  // the mismatch parameter
  const RealType _gamma_epsilon;
  // the stress for the stress mismatch penalty
  aol::MultiVector<RealType> _stressProd;
  // the prestressing displacement
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _preStressDisp;
  // the inverse displacement, $\phi^{-1}$-identity, in all Cartesian directions as vector of the dofs and as multilinear interpolation
  aol::MultiVector<RealType> _invDispDOFs;
  const aol::DiscreteVectorFunctionDefault< ConfiguratorType, ConfiguratorType::Dim > _invDisp;
  // array encoding whether $\phi^{-1}$ is defined at the corresponding node
  qc::BitArray<ConfiguratorType::Dim> _invPhiDefined;

private:
  inline RealType evaluateSqrDifferenceDeriv( const RealType A, const RealType B ) const {
    return 2 * ( A - B );
  }

  inline RealType evaluatePhasefieldMismatchDeriv( const RealType A, const RealType B ) const {
    return 2 * ( A * aol::Sqr( 1 - B ) + ( A - 1 ) * aol::Sqr( B ) );
  }

public:
  PreStressedEnergyVariationWRTD( const typename ConfiguratorType::InitType &Grid,
                                  const aol::Vector<RealType> &UFixed,
                                  const aol::Vector<RealType> &UTemplate,
                                  const ElasticGradientType &ElasticGradient,
                                  const RealType Gamma_Epsilon,
                                  const aol::MultiVector<RealType> &PreStressDisp,
                                  const aol::MultiVector<RealType> &Stress ) :
    aol::FENonlinVectorOpInterface< ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, EnergyVariationWRTD<ConfiguratorType, ElasticGradientType, MismatchType> >( Grid ),
    _uFixed( Grid, UFixed ),
    _uTemplate( Grid, UTemplate ),
    _elasticGradient( ElasticGradient ),
    _gamma_epsilon( Gamma_Epsilon ),
    _stressProd( Grid ),
    _preStressDisp( Grid, PreStressDisp ),
    _invDispDOFs( PreStressDisp, aol::STRUCT_COPY ),
    _invDisp( Grid, _invDispDOFs ),
    _invPhiDefined( Grid ) {
    ( PreStressedStressMatch<ConfiguratorType>( Grid ) ).apply( Stress, _stressProd );
    // generate the identity function
    aol::MultiVector<RealType> identity( PreStressDisp, aol::STRUCT_COPY );
    qc::DataGenerator<ConfiguratorType>( Grid ).generateIdentity( identity );
    // compute $\phi^{-1}$-identity
    qc::TransformFunction<RealType, ConfiguratorType::Dim> transformOp( Grid );
    transformOp.setDeformation( PreStressDisp );
    transformOp.transform( identity, _invDispDOFs, _invPhiDefined );
    _invDispDOFs -= identity;
  }

  void getNonlinearity ( aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                       const typename ConfiguratorType::ElementType &El,
                       int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                       aol::Vec<ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {
    NL.setZero();

    //! if \f$\phi_2^{-1}\f$ is not defined at this position, return 0
    if ( !_invPhiDefined.elementTrue( El ) )
      return;

    //! compute \f$\phi_1(x)\f$ and \f$\phi_2^{-1}(x)\f$
    typename ConfiguratorType::VecType transformedLocalCoord1, transformedLocalCoord2;
    typename ConfiguratorType::ElementType transformedEl1, transformedEl2;
    if( !qc::transformCoord<ConfiguratorType> ( this->_initializer, DiscFuncs, El, QuadPoint, RefCoord, transformedEl1, transformedLocalCoord1 ) )
      // if point $x$ was displaced outside the grid, return 0
      return;
    if( !qc::transformCoord<ConfiguratorType> ( this->_initializer, _invDisp, El, QuadPoint, RefCoord, transformedEl2, transformedLocalCoord2 ) )
      // if point $x$ was displaced outside the grid, return 0
      return;

    //! compute the value of the deformed shape and the undeformed one
    RealType uTemplate = _uTemplate.evaluate( transformedEl1, transformedLocalCoord1 ),
             uFixed    = _uFixed.evaluate( transformedEl2, transformedLocalCoord2 );

    //! if \f$x\f$ was not displaced outside the grid, return result
    _uTemplate.evaluateGradient ( transformedEl1, transformedLocalCoord1, NL );
    typename ConfiguratorType::MatType dphi;
    _invDisp.evaluateGradientAtQuadPoint( El, QuadPoint, dphi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dphi[i][i] += 1.;
    switch ( MismatchType ) {
    case SquaredDifference:  NL *= _gamma_epsilon * evaluateSqrDifferenceDeriv( uTemplate, uFixed ) * aol::Abs( dphi.det() ); break;
    case PhasefieldMismatch: NL *= _gamma_epsilon * evaluatePhasefieldMismatchDeriv( uTemplate, uFixed ) * aol::Abs( dphi.det() ); break;
    }
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // compute the phase field term of the energy variation
    aol::FENonlinVectorOpInterface< ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, PreStressedEnergyVariationWRTD<ConfiguratorType, ElasticGradientType, MismatchType> >::applyAdd( Arg, Dest );

    // compute the hyperelastic energy variation
    _elasticGradient.applyAdd( Arg, Dest );

    // compute the stress mismatch penalty
    Dest += _stressProd;
  }
};

template <typename ConfiguratorType, MismatchPenaltyType MismatchType>
class EnergyVariationWRTDDefault :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > WeightedEnergyDensityType;
  typedef qc::HyperelasticGradient<ConfiguratorType,WeightedEnergyDensityType> WeightedHyperelasticGradientType;

  // the default elastic energy
  const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> _energyDensity;
  const WeightedEnergyDensityType _weightedDensity;
  const WeightedHyperelasticGradientType _hyperelasticGradient;

  // the energy term itself
  const EnergyVariationWRTD<ConfiguratorType, WeightedHyperelasticGradientType, MismatchType> _energyVariationWRTD;

public:
  EnergyVariationWRTDDefault( const typename ConfiguratorType::InitType &Grid,
                              const aol::Vector<RealType> &UFixed,
                              const aol::Vector<RealType> &UTemplate,
                              const RealType LengthEnergyWeight,
                              const RealType VolumeEnergyWeight,
                              const RealType Gamma_Epsilon,
                              const aol::Vector<RealType> &Weight ) :
    _energyDensity( LengthEnergyWeight, 0., VolumeEnergyWeight ),
    _weightedDensity( _energyDensity, Grid, Weight ),
    _hyperelasticGradient( Grid, _weightedDensity ),
    _energyVariationWRTD( Grid, UFixed, UTemplate, _hyperelasticGradient, Gamma_Epsilon ) {}

  EnergyVariationWRTDDefault( const typename ConfiguratorType::InitType &Grid,
                              const aol::Vector<RealType> &UFixed,
                              const aol::Vector<RealType> &UTemplate,
                              const RealType LengthEnergyWeight,
                              const RealType VolumeEnergyWeight,
                              const RealType Gamma_Epsilon,
                              const aol::Vector<RealType> &Weight,
                              const aol::MultiVector<RealType> &Stress ) :
    _energyDensity( LengthEnergyWeight, 0., VolumeEnergyWeight ),
    _weightedDensity( _energyDensity, Grid, Weight ),
    _hyperelasticGradient( Grid, _weightedDensity ),
    _energyVariationWRTD( Grid, UFixed, UTemplate, _hyperelasticGradient, Gamma_Epsilon, Stress ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    _energyVariationWRTD.applyAdd( Arg, Dest );
  }
};

template <typename ConfiguratorType>
class CommonDisplacementEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > WeightedEnergyDensityType;
  typedef qc::HyperelasticDeformEnergy<ConfiguratorType,WeightedEnergyDensityType> WeightedHyperelasticEnergyType;

  // the hyperelastic energy
  const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> _energyDensity;
  std::vector<WeightedEnergyDensityType*> _weightedDensity;
  std::vector<WeightedHyperelasticEnergyType*> _energy;

  // the stress energy
  aol::MultiVector<RealType> _stressProd;

public:
  CommonDisplacementEnergy( const typename ConfiguratorType::InitType &Grid,
                            const aol::RandomAccessContainer<aol::MultiVector<RealType> > &Displacements,
                            const aol::MultiVector<RealType> &Stress,
                            const RealType LengthEnergyWeight,
                            const RealType VolumeEnergyWeight,
                            const aol::MultiVector<RealType> &Weights ) :
    _energyDensity( LengthEnergyWeight, 0., VolumeEnergyWeight ),
    _weightedDensity( Displacements.size() ),
    _energy( Displacements.size() ),
    _stressProd( Grid ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Displacements.size(); i++ ) {
      _weightedDensity[i] = new WeightedEnergyDensityType( _energyDensity, Grid, Weights[i] );
      _energy[i] = new WeightedHyperelasticEnergyType( Grid, *_weightedDensity[i], Displacements[i] );
    }
    ( StressMatch<ConfiguratorType>( Grid ) ).apply( Stress, _stressProd );
  }

  ~CommonDisplacementEnergy() {
    for ( unsigned int i = 0; i < _energy.size(); i++ ) {
      delete _weightedDensity[i];
      delete _energy[i];
    }
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < static_cast<int>( _energy.size() ); i++ ) {
      aol::Scalar<RealType> dest;
      _energy[i]->apply( Arg, dest );
#ifdef _OPENMP
#pragma omp critical
#endif
      { Dest += dest; }
    }
    Dest -= _stressProd * Arg;
  }
};

template <typename ConfiguratorType>
class CommonDisplacementEnergyVariation :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > WeightedEnergyDensityType;
  typedef qc::HyperelasticDeformGradient<ConfiguratorType,WeightedEnergyDensityType> WeightedHyperelasticGradientType;

  // the hyperelastic energy
  const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> _energyDensity;
  std::vector<WeightedEnergyDensityType*> _weightedDensity;
  std::vector<WeightedHyperelasticGradientType*> _energyGradient;

  // the stress energy
  aol::MultiVector<RealType> _stressProd;

public:
  CommonDisplacementEnergyVariation( const typename ConfiguratorType::InitType &Grid,
                                     const aol::RandomAccessContainer<aol::MultiVector<RealType> > &Displacements,
                                     const aol::MultiVector<RealType> &Stress,
                                     const RealType LengthEnergyWeight,
                                     const RealType VolumeEnergyWeight,
                                     const aol::MultiVector<RealType> &Weights ) :
    _energyDensity( LengthEnergyWeight, 0., VolumeEnergyWeight ),
    _weightedDensity( Displacements.size() ),
    _energyGradient( Displacements.size() ),
    _stressProd( Grid ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Displacements.size(); i++ ) {
      _weightedDensity[i] = new WeightedEnergyDensityType( _energyDensity, Grid, Weights[i] );
      _energyGradient[i] = new WeightedHyperelasticGradientType( Grid, *_weightedDensity[i], Displacements[i] );
    }
    ( StressMatch<ConfiguratorType>( Grid ) ).apply( Stress, _stressProd );
  }

  ~CommonDisplacementEnergyVariation() {
    for ( unsigned int i = 0; i < _energyGradient.size(); i++ ) {
      delete _weightedDensity[i];
      delete _energyGradient[i];
    }
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < static_cast<int>( _energyGradient.size() ); i++ ) {
      aol::MultiVector<RealType> dest( Arg, aol::STRUCT_COPY );
      _energyGradient[i]->apply( Arg, dest );
#ifdef _OPENMP
#pragma omp critical
#endif
      { Dest += dest; }
    }
    Dest -= _stressProd;
  }
};

template <typename RealType, typename VectorType, typename QuadraticOpType>
class HilbertSpaceQuadraticEnergyOp:
  public aol::Op<VectorType, aol::Scalar<RealType> > {
protected:
  const QuadraticOpType &_quadraticOp;
  const VectorType &_linearOp;
public:
  HilbertSpaceQuadraticEnergyOp( const QuadraticOpType &QuadraticOp, const VectorType &LinearOp ) :
    _quadraticOp( QuadraticOp ),
    _linearOp( LinearOp ) {}
  void applyAdd( const VectorType &Arg, aol::Scalar<RealType> &Dest ) const {
    VectorType aux( Arg, aol::STRUCT_COPY );
    _quadraticOp.apply( Arg, aux );
    aux += _linearOp;
    Dest += aux * Arg;
  }
};

template <typename RealType, typename VectorType, typename QuadraticOpType>
class HilbertSpaceQuadraticEnergyGradientOp:
  public aol::Op<VectorType> {
protected:
  const QuadraticOpType &_quadraticOp;
  const VectorType &_linearOp;
public:
  HilbertSpaceQuadraticEnergyGradientOp( const QuadraticOpType &QuadraticOp, const VectorType &LinearOp ) :
    _quadraticOp( QuadraticOp ),
    _linearOp( LinearOp ) {}
  void applyAdd( const VectorType &Arg, VectorType &Dest ) const {
    VectorType aux( Arg, aol::STRUCT_COPY );
    _quadraticOp.apply( Arg, aux );
    aux *= 2.;
    Dest += aux;
    Dest += _linearOp;
  }
};


/**************************************************************
 * Classes for averaging of given shapes and related problems.
 **************************************************************/

template <typename ConfiguratorType>
class StressPCAOp {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  // grid over the examined domain
  const typename ConfiguratorType::InitType _grid;
  // number of given shapes
  const int _numberOfShapes;
  aol::MultiVector<RealType> _images;
  // the directory for the results
  char _destDirectory[1024];
  // model parameters
  const RealType _epsFactor;
  // weights of the hyperelastic deformation energy
  const RealType _volumeEnergyWeight, _lengthEnergyWeight;
  const bool _constantElasticParameters;
  // the displacments to the average shape
  aol::RandomAccessContainer<aol::MultiVector<RealType> > _displacements;
  // the average shape (as binary image)
  ArrayType _shape;

public:
  StressPCAOp( aol::ParameterParser &Parser ) :
    _grid( Parser.getInt( "GridDepth" ), ConfiguratorType::Dim ),
    _numberOfShapes( Parser.getInt( "numberOfImages" ) ),
    _images( _numberOfShapes, _grid.getNumberOfNodes() ),
    _epsFactor( Parser.getDouble( "epsFactor" ) ),
    _volumeEnergyWeight( Parser.getDouble( "volumeEnergyWeight" ) ),
    _lengthEnergyWeight( Parser.getDouble( "lengthEnergyWeight" ) ),
    _constantElasticParameters( Parser.getInt( "constantElasticParameters" ) ),
    _displacements( _numberOfShapes ),
    _shape( Parser.getString( "averageImageFile" ) ) {

    Parser.getString( "destDirectory", _destDirectory );

    // load deformations of input shapes to the average shape (or rather the displacements)
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      _displacements[i].reallocate( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        // get filename of j-th component of i-th displacement
        char fileName[1024];
        sprintf( fileName, Parser.getString( "deformationFileNameTemplate" ).c_str(), i, j );
        // load that component
        cerr<<"load " + string( fileName ) + "...\n";
        ( ArrayType( _displacements[i][j], _grid, aol::FLAT_COPY ) ).load( fileName );
      }
    }

    // load the input objects and scale them to the interval [0,1]
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 1; i <= _numberOfShapes; i++ ) {
      // get filename of i-th object
      char fileName[1024];
      sprintf( fileName, Parser.getString( "imageFileNameTemplate" ).c_str(), i );
      // load i-th object
      cerr<<"load " + string( fileName ) + "...\n";
      ArrayType tmpArray( _images[i-1], _grid, aol::FLAT_COPY );
      tmpArray.load( fileName );
      // intensify the contrast
      tmpArray.addToAll( -tmpArray.getMinValue() );
      tmpArray /= tmpArray.getMaxValue();
    }

    // make the shape-image values lie in [0,1]
    _shape.addToAll( -_shape.getMinValue() );
    _shape /= _shape.getMaxValue();
    // if the shape was the first of the input shapes, deform it to the average shape
    /*ArrayType origShape( _shape, aol::DEEP_COPY );
    qc::InvDeformImage<ConfiguratorType,ArrayType>( origShape, _grid, _shape, _d[0] );*/

    // notify the user
    cerr << "Parameters loaded." << endl;
  }

  void execute() {
    cerr << "Compute the signed distance function of the average shape and the weight for local stress averaging..." << endl;
    ArrayType weight( _shape, aol::DEEP_COPY ), shapeSignedDist( _shape, aol::STRUCT_COPY );
    weight.addToAll( -.5 );
    ( typename qc::SignedDistanceOpTrait<ConfiguratorType,ConfiguratorType::Dim>::OpType( _grid ) ).apply( weight, shapeSignedDist );
    weight = shapeSignedDist;
    weight /= _grid.H(); // now each pixel corresponds to an increase/decrease by one
    weight.addToAll( -3.5 * _epsFactor + 2. ); // diminished stresses (due to phasefield) are within 3.5 epsilon band, accuracy of weight is 2 pixel
    weight.clamp( 1.e-8, 1. );

    cerr << "Compute the stresses of the deformed shapes..." << endl;
    aol::RandomAccessContainer<aol::MultiVector<RealType> > stress( _numberOfShapes );
    qc::HyperelasticEnergyDensityDefault<ConfiguratorType> density( _lengthEnergyWeight, 0., _volumeEnergyWeight );
    typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > ElasticEnergyType;
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      ArrayType energyWeight( _images[i], _grid, aol::STRUCT_COPY );
      if ( _constantElasticParameters )
        energyWeight.setAll( 1. );
      else
        energyWeight = _images[i];
      qc::DeformImage<ConfiguratorType>( _shape, _grid, energyWeight, _displacements[i] ); // (assuming that this is the energy density used for averaging)
      ElasticEnergyType weightedDensity( density, _grid, energyWeight );
      CauchyStressOp<ConfiguratorType,ElasticEnergyType> stressOp( _grid, weightedDensity );
      stress[i].reallocate( ConfiguratorType::Dim * ConfiguratorType::Dim, _shape.size() );
      stressOp.apply( _displacements[i], stress[i] );
    }

    cerr << "Perform a PCA of the boundary stresses..." << endl;
    typename ConfiguratorType::ArrayType kernel( _grid );
    ( qc::DataGenerator<ConfiguratorType>( _grid ) ).generatePeriodicSphereCharacFunc( kernel, 7 * _grid.H() * _epsFactor );
    aol::RandomAccessContainer<aol::MultiVector<RealType> > modes;
    aol::Vector<RealType> variances( _numberOfShapes );
    //DiffusedBoundaryL2Product<ConfiguratorType> dotProd( _grid, shapeSignedDist, _shape, kernel );
    //( aol::PCAOp<RealType,aol::MultiVector<RealType> >( dotProd, _numberOfShapes ) ).modesAndVariances( stress, modes, variances );

    //! alternatively to the two last rows:
    // subtract the mean stress from the stresses
    aol::MultiVector<RealType> mean( stress[0], aol::STRUCT_COPY );
    aol::MeanOp<aol::MultiVector<RealType> >().meanAndVariation( stress, mean, stress );
    // compute the (inward) normal stresses (only on the object!)
    StressMultWeightedNormal<ConfiguratorType> normalOp( _grid, shapeSignedDist, _shape );
    aol::RandomAccessContainer<aol::MultiVector<RealType> >normalStress( stress.size() );
    for ( int i = 0; i < stress.size(); i++ ) {
      normalStress[i].reallocate( _grid );
      normalOp.apply( stress[i], normalStress[i] );
    }

    // locally average the normal stresses
    WeightedConvolutionOp<ConfiguratorType> weightedConvOp( weight, kernel );
    for ( int i = 0; i < stress.size(); i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        typename ConfiguratorType::ArrayType aux( normalStress[i][j], _grid, aol::FLAT_COPY );
        weightedConvOp.apply( aux, aux );
      }

    // compute the covariance matrix
    qc::LevelsetVectorMassOp<ConfiguratorType,ConfiguratorType::Dim> levelsetL2Prod( _grid, shapeSignedDist );
    aol::FullMatrix<RealType> covMatrix;
    aol::CovarianceMatrixOp<RealType,aol::MultiVector<RealType> >( levelsetL2Prod ).apply( normalStress, covMatrix );
    // decompose covariance matrix (first row of decompResult contains eigenvalues in decreasing order, further rows contain corresponding eigenvectors)
    aol::MultiVector<RealType> decompResult;
    ( aol::QREigenvectorOp<aol::FullMatrix<RealType> >() ).apply( covMatrix, decompResult );
    // compute modes of variation (note that eigenvalues are ordered from small to large)
    modes.reallocate( _numberOfShapes );
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      modes[i].reallocate( stress[0] );
      for ( int j = 0; j < _numberOfShapes; j++ )
        modes[i].addMultiple( stress[j], decompResult[_numberOfShapes-i][j] );
      variances[i] = decompResult[0][_numberOfShapes-i-1];
      modes[i] /= sqrt( variances[i] );
    }

    cerr << "Saving the results..." << endl;
    for ( int mode = 0; mode < stress.size(); mode++ ) {
      char fileName[1024];
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          sprintf( fileName, "%s/stressMode_%d_%d%d.bz2", _destDirectory, mode + 1, i, j );
          ( typename ConfiguratorType::ArrayType( stress[mode][i*ConfiguratorType::Dim+j], _grid, aol::FLAT_COPY ) ).save( fileName, qc::PGM_DOUBLE_BINARY );
        }
    }
    char fileName[1024];
    sprintf( fileName, "%s/variances.dat", _destDirectory );
    variances.saveAsVector( fileName );
    sprintf( fileName, "%s/variances.txt", _destDirectory );
    variances.saveASCII( fileName );
  }
};


//! Given input shapes O_i and displacements d_i, s.t. the corresponding deformation \phi_i = id + d_i deforms
//! O_i into some (average) shape O, the induced Cauchy stresses  \f$ \sigma_i = \sigma[\phi_i] = W_{,A}(D\phi_i) \text{cof}(D\Phi)^{-1} \f$
//! are computed simultaneously.
//! NOTE The average object O is actually not used in the code!
template <typename ConfiguratorType>
class CauchyStressComputation {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  // grid over the examined domain
  const aol::ParameterParser& _parser;
  const typename ConfiguratorType::InitType _grid;
  // number of given shapes
  const int _numberOfShapes;
  aol::MultiVector<RealType> _images;
  // the directory for the results
  char _destDirectory[1024];
  // model parameters
  const RealType _epsFactor;
  // weights of the hyperelastic deformation energy
  const RealType _volumeEnergyWeight, _lengthEnergyWeight;
  const bool _constantElasticParameters;
  const bool _doNotSubtractMeanStress;
  // the displacments to the average shape
  aol::RandomAccessContainer<aol::MultiVector<RealType> > _displacements;
  // the average shape (as binary image)
  ArrayType _shape;

public:
  CauchyStressComputation( const aol::ParameterParser &Parser ) :
    _parser( Parser ),
    _grid( Parser.getInt( "GridDepth" ), ConfiguratorType::Dim ),
    _numberOfShapes( Parser.getInt( "numberOfImages" ) ),
    _images( _numberOfShapes, _grid.getNumberOfNodes() ),
    _epsFactor( Parser.getDouble( "epsFactor" ) ),
    _volumeEnergyWeight( Parser.getDouble( "volumeEnergyWeight" ) ),
    _lengthEnergyWeight( Parser.getDouble( "lengthEnergyWeight" ) ),
    _constantElasticParameters( Parser.getInt( "constantElasticParameters" ) ),
    _doNotSubtractMeanStress( Parser.hasVariable( "doNotSubtractMeanStress" ) ? Parser.getInt( "doNotSubtractMeanStress" ) : false ),
    _displacements( _numberOfShapes ),
    _shape( Parser.getString( "averageImageFile" ) ) {

    Parser.getString( "destDirectory", _destDirectory );

    // load deformations of input shapes to the average shape (or rather the displacements)
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      _displacements[i].reallocate( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        // get filename of j-th component of i-th displacement
        char fileName[1024];
        sprintf( fileName, Parser.getString( "deformationFileNameTemplate" ).c_str(), i, j );
        cerr<<"load " + string( fileName ) + "...\n";
        // load that component
        ( ArrayType( _displacements[i][j], _grid, aol::FLAT_COPY ) ).load( fileName );
      }
    }

    // load the input objects and scale them to the interval [0,1]
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 1; i <= _numberOfShapes; i++ ) {
      // get filename of i-th object
      char fileName[1024];
      sprintf( fileName, Parser.getString( "imageFileNameTemplate" ).c_str(), i );
      cerr<<"load " + string( fileName ) + "...\n";
      // load i-th object
      ArrayType tmpArray( _images[i-1], _grid, aol::FLAT_COPY );
      tmpArray.load( fileName );
      // intensify the contrast
      tmpArray.addToAll( -tmpArray.getMinValue() );
      tmpArray /= tmpArray.getMaxValue();
    }

    // make the shape-image values lie in [0,1]
    _shape.addToAll( -_shape.getMinValue() );
    _shape /= _shape.getMaxValue();
    // if the shape was the first of the input shapes, deform it to the average shape
    /*ArrayType origShape( _shape, aol::DEEP_COPY );
    qc::InvDeformImage<ConfiguratorType,ArrayType>( origShape, _grid, _shape, _d[0] );*/

    // notify the user
    cerr << "Parameters loaded." << endl;
  }

  void execute() {
    cerr << "Compute the stresses of the deformed shapes..." << endl;
    aol::RandomAccessContainer<aol::MultiVector<RealType> > stress( _numberOfShapes );
    qc::HyperelasticEnergyDensityDefault<ConfiguratorType> density( _lengthEnergyWeight, 0., _volumeEnergyWeight );
    typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > ElasticEnergyType;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      ArrayType energyWeight( _images[i], _grid, aol::STRUCT_COPY );
      if ( _constantElasticParameters )
        energyWeight.setAll( 1. );
      else
        dilateAndClamp<ConfiguratorType>( _images[i], _grid, energyWeight );
      ElasticEnergyType weightedDensity( density, _grid, energyWeight );
      CauchyStressOp<ConfiguratorType,ElasticEnergyType> stressOp( _grid, weightedDensity );
      stress[i].reallocate( ConfiguratorType::Dim * ConfiguratorType::Dim, _shape.size() );
      stressOp.apply( _displacements[i], stress[i] );
    }

    if ( !_doNotSubtractMeanStress ) {
      cerr << "Substract the mean stress (so that stresses due to phasefield contraction do not count)..." << endl;
      aol::MultiVector<RealType> meanStress( stress[0], aol::STRUCT_COPY );
      for ( int i = 0; i < _numberOfShapes; i++ )
        meanStress += stress[i];
      meanStress /= _numberOfShapes;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < _numberOfShapes; i++ )
        stress[i] -= meanStress;
    }

    cerr << "Saving the results..." << endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 0; k < _numberOfShapes; k++ ) {
      char fileName[1024];
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
	  //sprintf( fileName, "%s/stress_%d_%d%d.bz2", _destDirectory, k + 1, i, j );
	  sprintf( fileName, _parser.getString( "stressFileNameTemplate" ).c_str(), k + 1, i, j );
          ( typename ConfiguratorType::ArrayType( stress[k][i*ConfiguratorType::Dim+j], _grid, aol::FLAT_COPY ) ).save( fileName, qc::PGM_DOUBLE_BINARY );
        }
    }
  }
};

template <typename ConfiguratorType>
class PreStressRelaxationOp {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  // finest grid over the examined domain
  const typename ConfiguratorType::InitType _grid;
  // number of input objects and images
  const int _numberOfImages;
  MultiLevelMultiArray<ConfiguratorType> _images;
  // the directory for the results
  char _destDirectory[1024];
  // weights of the hyperelastic deformation energy
  const RealType _lengthEnergyWeight, _volumeEnergyWeight;
  const bool _constantElasticParameters;
  // coarsest and finest level for the multiscale method (only needed, if nonlinear problem is solved, ie. not LINEARIZED_METHOD)
  const int _coarsestLevel, _finestLevel;
  // number of gradient descent steps for the deformation (only needed, if nonlinear problem is solved, ie. not LINEARIZED_METHOD)
  const int _maxSteps;
  // average shape (as binary image so that it can be used as weight)
  qc::MultilevelArray<RealType> _average;
  // prestress displacements
  MultiLevelMultiArray<ConfiguratorType> _displacements;
  // the additional stresses to be imposed
  const int _modeLower, _modeNumber; // if computation shall only be for few imposed stress
  MultiLevelMultiArray<ConfiguratorType> _stresses;
  const RealType _stressAmpl; // (only needed, if nonlinear problem is solved, ie. not LINEARIZED_METHOD)
  // the resulting deformations
  MultiLevelMultiArray<ConfiguratorType> _commonDisplacements;

public:
  PreStressRelaxationOp( const aol::ParameterParser &Parser ) :
    // define the finest grid over the domain
    _grid( Parser.getInt( "GridDepth" ), ConfiguratorType::Dim ),
    // load number of objects and create the image-vector
    _numberOfImages( Parser.getInt( "numberOfImages" ) ),
    _images( _grid, 1, _numberOfImages ),
    // load hyperelastic and phasefield parameters
    _lengthEnergyWeight( Parser.getDouble( "lengthEnergyWeight" ) ),
    _volumeEnergyWeight( Parser.getDouble( "volumeEnergyWeight" ) ),
    _constantElasticParameters( Parser.getInt( "constantElasticParameters" ) ),
    // load coarsest and finest level as well as number of gradient descent steps per level
    _coarsestLevel( Parser.getInt( "coarsestLevel" ) ),
    _finestLevel( Parser.getInt( "GridDepth" ) ),
    _maxSteps( Parser.getInt( "maxSteps" ) ),
    // create the average shape and the displacements
    _average( _grid ),
    _displacements( _grid, _numberOfImages, ConfiguratorType::Dim ),
    // create the additionally imposed stresses
    _modeLower( 1 /*Parser.hasVariable( "StressMode" ) ? Parser.getInt( "StressMode" ) : 1*/ ),
    _modeNumber( Parser.hasVariable( "StressMode" ) ? Parser.getInt( "StressMode" ) : _numberOfImages ),
    _stresses( _grid, _modeNumber, ConfiguratorType::Dim * ConfiguratorType::Dim ),
    _stressAmpl( Parser.getDouble( "stressAmplification" ) ),
    _commonDisplacements( _grid, _modeNumber, ConfiguratorType::Dim ) {

    // read in the directory, where results are to be saved
    Parser.getString( "destDirectory", _destDirectory );

    // load the input objects and scale them to the interval [0,1]
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 1; i <= _numberOfImages; i++ ) {
      // get filename of i-th object
      char fileName[1024];
      sprintf( fileName, Parser.getString( "imageFileNameTemplate" ).c_str(), i );
      // load i-th object
      cerr<<"load " + string( fileName ) + "...\n";
      ArrayType tmpArray( _images.getArray( i - 1 ), aol::FLAT_COPY );
      tmpArray.load( fileName );
      // intensify the contrast
      tmpArray.addToAll( -tmpArray.getMinValue() );
      tmpArray /= tmpArray.getMaxValue();
    }

    // load average shape
    cerr<<"load "<<Parser.getString( "averageImageFile" )<<"..."<<endl;
    ( ArrayType( _average.current(), aol::FLAT_COPY ) ).load( Parser.getString( "averageImageFile" ).c_str() );
    ( _average.current() ).addToAll( -( _average.current() ).getMinValue() );
    _average.current() /= ( _average.current() ).getMaxValue();

    // initialise average shape and displacements (already happened automatically)
    /*_average.current().setZero();
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      aol::MultiVector<RealType> displacement;
      _displacements.appendMultiArrayReferenceTo( i, displacement );
      displacement.setZero();
    }
    aol::MultiVector<RealType> displacement;
    _commonDisplacement.appendMultiArrayReferenceTo( 0, displacement );
    displacement.setZero();*/

   // load the imposed Cauchy stresses
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 0; k < _modeNumber; k++ ) {
      char fileName[1024];
      aol::MultiVector<RealType> cauchyStress;
      _stresses.appendMultiArrayReferenceTo( k, cauchyStress );
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          sprintf( fileName, Parser.getString( "stressFileNameTemplate" ).c_str(), _modeLower + k, i, j );
          cerr<<"load " + string( fileName ) + "...\n";
          ArrayType( cauchyStress[i*ConfiguratorType::Dim+j], _grid, aol::FLAT_COPY ).load( fileName );
        }
      cauchyStress *= _stressAmpl * _numberOfImages;
    }

    // load the prestressing displacements
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numberOfImages; i++ ) {
      // load the displacement
      char fileName[1024];
      aol::MultiVector<RealType> displacement;
      _displacements.appendMultiArrayReferenceTo( i, displacement );
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        sprintf( fileName, Parser.getString( "deformationFileNameTemplate" ).c_str(), _finestLevel, i, j );
	// check if this filename exists
	if ( !fopen( fileName, "rb" ) )
          sprintf( fileName, Parser.getString( "deformationFileNameTemplate" ).c_str(), i, j );	
        cerr<<"load " + string( fileName ) + "...\n";
        ArrayType( displacement[j], _grid, aol::FLAT_COPY ).load( fileName );
      }
    }

    // generate the needed coarser scales of images, phase field, and displacements
    if ( _coarsestLevel < _finestLevel ) {
      _images.levRestrict( _coarsestLevel, _finestLevel );
      _average.levRestrict( _coarsestLevel, _finestLevel );
      _displacements.levRestrict( _coarsestLevel, _finestLevel );
      _stresses.levRestrict( _coarsestLevel, _finestLevel );
      _commonDisplacements.levRestrict( _coarsestLevel, _finestLevel );
    }

/*if ( !_constantElasticParameters ) {
  aol::RandomAccessContainer<aol::MultiVector<RealType> > displacements;
  _displacements.getMultiArrays( displacements );
  for ( int i = 0; i < _numberOfImages; i++ ) {
    ArrayType aux( _grid );
    dilateAndClamp<ConfiguratorType>( _images.getArray( i ), _grid, aux, 3.5, 0., 1. );
    aux *= -1.;
    aux.addToAll( 1. );
    typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,AntiInterpenetrationHyperelasticEnergyDensity<ConfiguratorType> > WeightedEnergyDensityType;
    typedef qc::HyperelasticEnergy<ConfiguratorType,WeightedEnergyDensityType> WeightedHyperelasticEnergyType;
    typedef qc::HyperelasticGradient<ConfiguratorType,WeightedEnergyDensityType> WeightedHyperelasticGradientType;
    AntiInterpenetrationHyperelasticEnergyDensity<ConfiguratorType> energyDensity( 1e-0, 1e2 );
    WeightedEnergyDensityType weightedDensity( energyDensity, _grid, aux );
    WeightedHyperelasticEnergyType energy( _grid, weightedDensity );
    WeightedHyperelasticGradientType grad( _grid, weightedDensity );
    aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > gradientDescent( _grid, energy, grad, 500 );
    gradientDescent.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::USE_NONLINEAR_CG );
    aol::MultiVector<RealType> initialDisp( displacements[i], aol::DEEP_COPY );
    gradientDescent.apply( initialDisp, displacements[i] );
    saveDisplacement<ConfiguratorType>( displacements[i], _grid, _destDirectory, "smoothedDisp", i, false );
  }
  abort();
}
*/
    // notify the user
    cerr << "Parameters loaded." << endl;
  }

  //! Afterwards, "DefEnergyWeights[i]" is one inside the ith object and close to zero outside.
  inline void computeDeformationEnergyWeights( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Images, aol::MultiVector<RealType> &DefEnergyWeights, const RealType PixelOverlap = 3.5 ){
    DefEnergyWeights.reallocate( Images );
    if ( _constantElasticParameters )
      DefEnergyWeights.setAll( 1. );
    else
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < DefEnergyWeights.numComponents(); i++ )
        dilateAndClamp<ConfiguratorType>( Images[i], Grid, DefEnergyWeights[i], PixelOverlap );
  }

  void findMinimizingCommonDisplacement( const typename ConfiguratorType::InitType &Grid,
                                         const aol::MultiVector<RealType> &Images,
                                         const aol::Vector<RealType> &/*Average*/,
                                         const aol::MultiVector<RealType> &Stress,
                                         const aol::RandomAccessContainer<aol::MultiVector<RealType> > &Displacements,
                                         aol::MultiVector<RealType> &CommonDisplacement,
                                         const RealType LengthEnergyWeight,
                                         const RealType VolumeEnergyWeight,
                                         const int NumSteps ) {
    //! define total energy "E" and its variation with respect to the displacement, "DE" (can also be adapted for image averaging)
    aol::MultiVector<RealType> weights( Images );
    computeDeformationEnergyWeights( Grid, Images, weights );

    CommonDisplacementEnergy<ConfiguratorType> E( Grid, Displacements, Stress, LengthEnergyWeight, VolumeEnergyWeight, weights );
    CommonDisplacementEnergyVariation<ConfiguratorType> DE( Grid, Displacements, Stress, LengthEnergyWeight, VolumeEnergyWeight, weights );

    //! perform a few steps of gradient descent
    aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > gradientDescent( Grid, E, DE, NumSteps << 2 * ( _finestLevel - Grid.getGridDepth() ) );
    gradientDescent.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::USE_NONLINEAR_CG );
    aol::MultiVector<RealType> initialDisp( CommonDisplacement, aol::DEEP_COPY );
    gradientDescent.apply( initialDisp, CommonDisplacement );
  }

  void saveResults( const typename ConfiguratorType::InitType &Grid,
                    const aol::Vector<RealType> &Average,
                    const aol::MultiVector<RealType> &/*Images*/,
                    const aol::RandomAccessContainer<aol::MultiVector<RealType> > &/*Displacements*/,
                    const aol::MultiVector<RealType> &CommonDisplacement,
                    const char* DestDirectory,
                    const int Mode,
                    const int Index ) {
    char imageFileExtension[4]; // last character will contain the zero-byte, signifying the end of the string
    sprintf( imageFileExtension, "%s", ( ConfiguratorType::Dim == 2 ) ? "pgm" : "bz2" );
    char fileName[1024], prefix[1024];
    sprintf( prefix, "%s/mode_%d_", DestDirectory, Mode );

    //! save deformed average shape
    ArrayType average( Average, Grid, aol::STRUCT_COPY );
    qc::InvDeformImage<ConfiguratorType,aol::Vector<RealType> >( Average, Grid, average, CommonDisplacement );
    average *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
    sprintf( fileName, "%saverageShape_%d.%s", prefix, Index, imageFileExtension );
    average.save( fileName, ( ConfiguratorType::Dim == 2 ) ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );

    //! save the (single common) displacement
    for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
      sprintf( fileName, "%scommonDisplacement_%d_%d.bz2", prefix, Index, j );
      ArrayType disp( CommonDisplacement[j], Grid, aol::DEEP_COPY );
      disp /= _stressAmpl;
      disp.save( fileName, qc::PGM_DOUBLE_BINARY );
    }
  }

  void execute() {
#define LINEARIZED_METHOD // if system is linearized and linear system is solved
    // if a multilevel method is used, then on each level do... (only needed for alternative D)
    for ( int level = _coarsestLevel; level <= _finestLevel; level++ ) {
      cerr << endl << aol::color::red << "level " << level << " of " << _finestLevel << aol::color::reset << endl;

      //! refine mesh
      typename ConfiguratorType::InitType grid( level, ConfiguratorType::Dim );
      _images.setCurLevel( level );
      _average.setCurLevel( level );
      _displacements.setCurLevel( level );
      _stresses.setCurLevel( level );
      _commonDisplacements.setCurLevel( level );

      //! put all variables into a more convenient form
      aol::MultiVector<RealType> images;
      aol::RandomAccessContainer<aol::MultiVector<RealType> > displacements, stresses, commonDisplacements;
      _displacements.getMultiArrays( displacements );
      _stresses.getMultiArrays( stresses );
      _images.appendMultiArrayReferenceTo( 0, images );
      _commonDisplacements.getMultiArrays( commonDisplacements );

#ifdef LINEARIZED_METHOD
      qc::meanZeroShiftAndRotationProjector<ConfiguratorType> projector( grid, _average.current(), true );
#else
      qc::meanZeroShiftAndRotationProjector<ConfiguratorType> projector( grid, _average.current(), false );
#endif

#ifdef LINEARIZED_METHOD
      //! assemble system matrix
      const aol::MultiVector<RealType> zeros( grid );
      aol::Vector<int> bc( grid );
      typedef typename ConfiguratorType::MatrixType SubMatrixType; // alternative: aol::SparseMatrix<RealType>
      aol::BlockMatrix<SubMatrixType> systemMatrix( grid );
      typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > WeightedEnergyDensityType;
      const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> energyDensity( _lengthEnergyWeight, 0, _volumeEnergyWeight );
      aol::MultiVector<RealType> weights( images );
      computeDeformationEnergyWeights( grid, images, weights, 3.5 );

      aol::BitVector objectMaskAssembly( grid.getNumberOfNodes() );
      qc::BitArray<ConfiguratorType::Dim> objectMaskSave( qc::GridSize<ConfiguratorType::Dim>::createFrom( grid ) );
      if ( !_constantElasticParameters ) {
        aol::Vector<RealType> aux( grid );
        dilateAndClamp<ConfiguratorType>( _average.current(), grid, aux, 8, 0., 1. ); // 5 instead of 8 already hampers the movement!
        objectMaskAssembly.setNonzeroPatternFrom<RealType>( aux );
        objectMaskAssembly.invert();
        dilateAndClamp<ConfiguratorType>( _average.current(), grid, aux, 3.5, 0., 1. );
        for ( int k = 0; k < objectMaskSave.size(); k++ )
          objectMaskSave.set( k, aux[k] > 0. );
      }

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < _numberOfImages; i++ ) {
        cerr<<"assemble stiffness matrix of "<<i<<"th deformed object..."<<endl;
        aol::BlockMatrix<SubMatrixType> localSystemMatrix( grid );
        WeightedEnergyDensityType weightedDensity( energyDensity, grid, weights[i] );
        if ( _constantElasticParameters )
          qc::HyperelasticDeformHessian<ConfiguratorType,WeightedEnergyDensityType,SubMatrixType>( grid, weightedDensity, displacements[i] ).applyAdd( zeros, localSystemMatrix );
        else
          qc::HyperelasticDeformHessian<ConfiguratorType,WeightedEnergyDensityType,SubMatrixType>( grid, weightedDensity, displacements[i] ).applyAdd( zeros, localSystemMatrix, objectMaskAssembly );
        cerr<<"assembly of stiffness matrix of "<<i<<"th deformed object finished"<<endl;
#ifdef _OPENMP
#pragma omp critical // results not reproducible for "critical"? "ordered" would not work since parallel threads are started in wrong order so that parallelization does no longer work
#endif
        systemMatrix += localSystemMatrix;
      }

      /*//! alternative A: set zero displacement bc at the central node (using the splitting method, see examples)
      // careful! this is highly unstable and usually ignored by the solution around this node
      int nodeIndex = grid.getNumberOfNodes() / 2;
      bc[nodeIndex] = 1.;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          systemMatrix.getReference( i, j ).setRowColToZero( nodeIndex );
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        systemMatrix.getReference( i, i ).set( nodeIndex, nodeIndex, 1. );*/

      /*//! alternative B: set zero displacement bc at the boundary (using the splitting method, see examples)
      qc::GridDefinition::OldFullBoundaryNodeIterator it;
      for ( it = grid.begin(); it != grid.end(); ++it ) {
        int nodeIndex = 0;
        for ( int i = ConfiguratorType::Dim - 1; i >= 0 ; i-- )
          nodeIndex = nodeIndex * grid.getNumX() + (*it)[i];
        bc[nodeIndex] = 0.;
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          for ( int j = 0; j < ConfiguratorType::Dim; j++ )
            systemMatrix.getReference( i, j ).setRowColToZero( nodeIndex );
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          systemMatrix.getReference( i, i ).set( nodeIndex, nodeIndex, 1. );
      }*/
#endif // LINEARIZED_METHOD

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < _modeNumber; i++ ) {
#ifdef LINEARIZED_METHOD // if system is linearized and linear system is solved
        //! assemble rhs
        aol::MultiVector<RealType> rhs( grid );
        StressMatch<ConfiguratorType>( grid ).apply( stresses[i], rhs );
        for ( int l = 0; l < bc.size(); l++ )
          if ( bc[l] )
            for ( int k = 0; k < ConfiguratorType::Dim; k++ )
              rhs[k][l] = 0.;

        //! alternative 1: solve system
        aol::MultiVector<RealType> commonDisplacement( grid ), initialGuess( grid );
        const aol::BlockDiagonalPreconditioner<RealType,ConfiguratorType::Dim> preCond( systemMatrix );
        aol::PCGInverse<aol::MultiVector<RealType> >( systemMatrix, preCond, 1e-16, _maxSteps, aol::STOPPING_ABSOLUTE ).apply( rhs, commonDisplacements[i] );

        /*//! alternative 2: minimize the corresponding quadratic energy
        const HilbertSpaceQuadraticEnergyOp<RealType,aol::MultiVector<RealType>,aol::BlockMatrix<SubMatrixType> > e( systemMatrix, rhs );
        const HilbertSpaceQuadraticEnergyGradientOp<RealType,aol::MultiVector<RealType>,aol::BlockMatrix<SubMatrixType> > de( systemMatrix, rhs );
        aol::NewtonInfo<RealType> newtonInfo ( 1E-15, 300, 1E-20, 1000, aol::STOPPING_ABSOLUTE );
        //aol::QuasiNewtonIteration< RealType, aol::MultiVector<RealType>, aol::MultiVector<RealType> >( e, de, newtonInfo ).apply( initialGuess, commonDisplacements[i] );
        aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >( grid, e, de, 150 ).apply( initialGuess, commonDisplacements[i] );
        */

        //! since the displacement outside the shapes is not interesting and should be removed
        if ( !_constantElasticParameters ) {
          aol::MultiVector<RealType> displacement( commonDisplacements[i], aol::DEEP_COPY );
          qc::SmoothlyExtendImage<ConfiguratorType>( displacement, grid, commonDisplacements[i], objectMaskSave );
        }

        //! remove any mean translation or rotation
        projector.apply( commonDisplacements[i], commonDisplacements[i] );

#else // LINEARIZED_METHOD // let the shapes deform nonlinearly and minimize total nonlinear energy

        //! on the coarsest grid level test whether the global minimizing displacement requires a turn by 180° and change sign of stress in that case
        if ( level == _coarsestLevel ) {
          aol::MultiVector<RealType> displacement( grid );
          aol::Scalar<RealType> val1, val2;
          CommonDisplacementEnergy<ConfiguratorType> E( grid, displacements, stresses[i], _lengthEnergyWeight, _volumeEnergyWeight, images );
          E.apply( displacement, val1 );
          qc::DataGenerator<ConfiguratorType>( grid ).generateIdentity( displacement );
          displacement *= -2;
          E.apply( displacement, val2 );
          if ( val1 > val2 )
            for ( int k = _finestLevel; k >= _coarsestLevel; k-- ) {
              _stresses.setCurLevel( k );
              aol::MultiVector<RealType> stress;
              _stresses.appendMultiArrayReferenceTo( i, stress );
              stress *= -1;
          }
        }

        //! approximately solve for the (single common) displacement
        findMinimizingCommonDisplacement( grid, images, _average.current(), stresses[i], displacements, commonDisplacements[i], _lengthEnergyWeight, _volumeEnergyWeight, _maxSteps );

        //! remove any mean translation or rotation
        projector.apply( commonDisplacements[i], commonDisplacements[i] );
#endif // LINEARIZED_METHOD

        //! save data
        cerr << "common displacement " << i << " updated" << endl;
        saveResults( grid, _average[level], images, displacements, commonDisplacements[i], _destDirectory, _modeLower + i, level );
      }

      //! prolongate the results of this level
      _commonDisplacements.levProlongate();
    }
  }
};

template <typename ConfiguratorType>
class ShapeAverageOp {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  // finest grid over the examined domain
  const typename ConfiguratorType::InitType _grid;
  // number of input objects and input images, the objects themselves (as binary images) and their shapes (as phasefields)
  const int _numberOfImages;
  MultiLevelMultiArray<ConfiguratorType> _inputImages, _images, _phasefields;
  // the directory for the results
  char _destDirectory[1024];
  // weights of the hyperelastic deformation energy
  const RealType _lengthEnergyWeight, _volumeEnergyWeight;
  // phasefield parameters for phasefield segmentation and averaging
  RealType _alpha, _beta, _nu, _epsilon, _epsFactor;
  RealType _gamma, _mu;
  const bool _constantElasticParameters;
  const bool _jointSegmentationAveraging;
  const bool _imageSmoothing;
  const bool _imageMatching;
  // coarsest and finest level for the multiscale method
  const int _coarsestLevel, _finestLevel;
  // iteration number of alternate minimisations and gradient descent steps for the deformations per minimisation
  const int _maxIterations, _maxSteps;
  // average shape (as binary image or phase field)
  qc::MultilevelArray<RealType> _average;
  // displacements
  MultiLevelMultiArray<ConfiguratorType> _displacements;
  // if the elastic energies of the shapes are to be weighted differently
  aol::Vector<RealType> _shapeWeights;
  // if a mode of variation is to be computed using a stressmode, then the "_mode"th stress is prescribed
  const int _mode;
  MultiLevelMultiArray<ConfiguratorType> _stresses;
  const RealType _stressAmpl;
  char _deformationFileNameTemplate[1024];
  // if one only wants to calculate the displacements, the average is loaded and fixed
  const bool _updateShapeAverage;

public:
  ShapeAverageOp( const aol::ParameterParser &Parser, bool updateShapeAverage = true ) :
    // define the finest grid over the domain
    _grid( Parser.getInt( "GridDepth" ), ConfiguratorType::Dim ),
    // load number of objects and create the object-vector and shape vector
    _numberOfImages( Parser.getInt( "numberOfImages" ) ),
    _inputImages( _grid, 1, Parser.getInt( "imageSmoothing") ? _numberOfImages : 0 ),
    _images( _grid, 1, _numberOfImages ),
    _phasefields( _grid, 1, _numberOfImages ),
    // load hyperelastic and phasefield parameters
    _lengthEnergyWeight( Parser.getDouble( "lengthEnergyWeight" ) ),
    _volumeEnergyWeight( Parser.getDouble( "volumeEnergyWeight" ) ),
    _alpha( Parser.getDouble( "alpha" ) ),
    _beta( Parser.getDouble( "beta" ) ),
    _nu( Parser.getDouble( "nu" ) ),
    _epsFactor( Parser.getDouble( "epsFactor" ) ),
    _gamma( Parser.getDouble( "gamma" ) ),
    _mu( Parser.getDouble( "mu" ) ),
    _constantElasticParameters( Parser.getInt( "constantElasticParameters" ) ),
    _jointSegmentationAveraging( Parser.getInt( "jointSegmentationAveraging" ) ),
    _imageSmoothing( Parser.getInt( "imageSmoothing" ) ),
    _imageMatching( Parser.hasVariable( "imageMatching" ) ? Parser.getInt( "imageMatching" ) : false ),
    // load coarsest and finest level as well as number of optimisation steps per level and gradient descent steps per optimisation step
    _coarsestLevel( Parser.getInt( "coarsestLevel" ) ),
    _finestLevel( Parser.getInt( "GridDepth" ) ),
    _maxIterations( Parser.getInt( "maxIterations" ) ),
    _maxSteps( Parser.getInt( "maxSteps" ) ),
    // create the average shape and the displacements
    _average( _grid ),
    _displacements( _grid, _numberOfImages, ConfiguratorType::Dim ),
    // create the vector of weights for weighted averaging
    _shapeWeights( _numberOfImages ),
    // create the stress for the mode of variation
    _mode( Parser.getInt( "StressMode" ) ),
    _stresses( _grid, _numberOfImages, _mode ? ConfiguratorType::Dim * ConfiguratorType::Dim : 0 ),
    _stressAmpl( Parser.getDouble( "stressAmplification" ) ),
    _updateShapeAverage( updateShapeAverage ){

    // read in the directory, where results are to be saved
    Parser.getString( "destDirectory", _destDirectory );

    // load the input objects and scale them to the interval [0,1]
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 1; i <= _numberOfImages; i++ ) {
      // get filename of i-th object
      char fileName[1024];
      sprintf( fileName, Parser.getString( "imageFileNameTemplate" ).c_str(), i );
      cerr<<"load " + string( fileName ) + "...\n";
      // load i-th object
      ArrayType tmpArray( _images.getArray( i - 1 ), aol::FLAT_COPY );
      tmpArray.load( fileName );
      // intensify the contrast
      tmpArray.addToAll( -tmpArray.getMinValue() );
      tmpArray /= tmpArray.getMaxValue();
      // if input images and piecewise smooth images are distinguished:
      if ( _imageSmoothing )
        _inputImages.getArray( i - 1 ) = tmpArray;
    }

    // initialise average shape and displacements (already happened automatically)
    /*_average.current().setZero();
    for ( int i = 0; i < _numberOfImages; i++ ) {
      aol::MultiVector<RealType> displacement;
      _displacements.appendMultiArrayReferenceTo( i, displacement );
      displacement.setZero();
    }*/
    
    // if one only wants to calculate the displacements, the average is loaded and fixed
    if( !_updateShapeAverage ){
      cerr<<"load average from " + Parser.getString( "averageFileNameTemplate" ) + "...\n";
      ArrayType tmpArray( _average.current(), aol::FLAT_COPY );
      tmpArray.load( Parser.getString( "averageFileNameTemplate" ).c_str() );
    }

    // read in the weights for the different shapes for a weighted average
    if ( Parser.hasVariable( "shapeWeights" ) )
      Parser.getRealVec( "shapeWeights", _shapeWeights );
    else
      _shapeWeights.setAll( 1. );

    // if initial displacements are proposed, load these
    Parser.getString( "deformationFileNameTemplate", _deformationFileNameTemplate );
    if ( Parser.checkAndGetBool( "loadInitialDeformations" ) ) {      
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < _numberOfImages; i++ ) {
        // load the displacement
        char fileName[1024];
        aol::MultiVector<RealType> displacement;
        _displacements.appendMultiArrayReferenceTo( i, displacement );
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          sprintf( fileName, _deformationFileNameTemplate, _finestLevel, i, j );
          cerr<<"load " + string( fileName ) + "...\n";
          ArrayType( displacement[j], _grid, aol::FLAT_COPY ).load( fileName );
        }
      }
    }

    // if stress modes are imposed, ie. a stress is prescribed, initialise displacements and load the stress
    if ( _mode > 0 ) {
      // load the prescribed Cauchy stress
      char fileName[1024];
      aol::MultiVector<RealType> cauchyStress( ConfiguratorType::Dim * ConfiguratorType::Dim, _grid.getNumberOfNodes() );
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          sprintf( fileName, Parser.getString( "stressFileNameTemplate" ).c_str(), _mode, i, j );
          cerr<<"load " + string( fileName ) + "...\n";
          ArrayType( cauchyStress[i*ConfiguratorType::Dim+j], _grid, aol::FLAT_COPY ).load( fileName );
        }
      if ( Parser.hasVariable( "variancesFileName" ) ) {
        aol::Vector<RealType> variances;
        variances.loadAsVector( Parser.getString( "variancesFileName" ).c_str() );
      }
      cauchyStress *= _stressAmpl / _numberOfImages; // * sqrt( variances[_mode-1] );
      // compute the first Piola Kirchhoff stresses
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < _numberOfImages; i++ ) {
        aol::MultiVector<RealType> firstPiolaKirchhoffStress, displacement;
        _stresses.appendMultiArrayReferenceTo( i, firstPiolaKirchhoffStress );
        _displacements.appendMultiArrayReferenceTo( i, displacement );
        ( Cauchy2FirstPiolaKirchhoff<ConfiguratorType>( _grid, displacement ) ).apply( cauchyStress, firstPiolaKirchhoffStress );
      }
    }

    // generate the needed coarser scales of images, phase field, and displacements
    _inputImages.levRestrict( _coarsestLevel, _finestLevel );
    _images.levRestrict( _coarsestLevel, _finestLevel );
    _phasefields.levRestrict( _coarsestLevel, _finestLevel );
    _average.levRestrict( _coarsestLevel, _finestLevel );
    _displacements.levRestrict( _coarsestLevel, _finestLevel );
    _stresses.levRestrict( _coarsestLevel, _finestLevel );

    // notify the user
    cerr << "Parameters loaded." << endl;
  }

  /**
   * \brief solves \f$Ax=b\f$.
   */
  inline void solveLinearEquation( const typename ConfiguratorType::MatrixType &A, const aol::Vector<RealType> &B, aol::Vector<RealType> &X ) {
    aol::DiagonalPreconditioner< aol::Vector<RealType> > preCond( A );
    aol::PCGInverse<aol::Vector<RealType> > inv( A, preCond, 1e-16, 1000 );
    //aol::CGInverse<aol::Vector<RealType> > inv( A );
    inv.setStopping( aol::STOPPING_ABSOLUTE );
    inv.apply( B, X );
  }

  void findMinimizingAverageImage( const typename ConfiguratorType::InitType &Grid,
                                   const aol::MultiVector<RealType> &Images,
                                   const aol::RandomAccessContainer<aol::MultiVector<RealType> > &Displacements,
                                   aol::Vector<RealType> &Average,
                                   const RealType Gamma,
                                   const RealType Mu,
                                   const RealType /*Epsilon*/ ) {
    //! assemble the system matrix and rhs
    aol::Vector<RealType> rhs( Grid ), aux( Grid ), ones( Grid );
    ones.setAll( 1 );
    typename ConfiguratorType::MatrixType systemMatrix( Grid );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Images.numComponents(); i++ ) {
      // a local system matrix and rhs are only introduced for parallelization
      aol::Vector<RealType> localRhs( Grid );
      typename ConfiguratorType::MatrixType localSystemMatrix( Grid );
      qc::InvDeformedWeightMassOp<ConfiguratorType>( Grid, ones, Displacements[i] ).assembleAddMatrix( localSystemMatrix );
      qc::InvDeformedWeightMassOp<ConfiguratorType>( Grid, Images[i], Displacements[i] ).applyAdd( ones, localRhs );
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        rhs += localRhs;
        systemMatrix += localSystemMatrix;
      }
    }
    // add some small regularization to render the system solvable
    // (since the inverse displacements might not have been defined everywhere)
    aol::MassOp<ConfiguratorType>( Grid ).assembleAddMatrix( systemMatrix, Mu / Gamma );

    //! solve for the average object
    solveLinearEquation( systemMatrix, rhs, Average );
  }

  void findMinimizingAveragePhasefield( const typename ConfiguratorType::InitType &Grid,
                                        const aol::MultiVector<RealType> &Phasefields,
                                        const aol::RandomAccessContainer<aol::MultiVector<RealType> > &Displacements,
                                        aol::Vector<RealType> &Average,
                                        const RealType Gamma,
                                        const RealType Mu,
                                        const RealType Epsilon ) {
    //! assemble the system matrix and rhs
    aol::Vector<RealType> rhs( Grid ), aux( Grid ), ones( Grid );
    ones.setAll( 1 );
    typename ConfiguratorType::MatrixType systemMatrix( Grid );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Phasefields.numComponents(); i++ ) {
      // a local system matrix and rhs are only introduced for parallelization
      typename ConfiguratorType::MatrixType localSystemMatrix( Grid );
      aol::Vector<RealType> localRhs( Grid ), localAux( Phasefields[i] );
      localAux.addToAll( -1. );
      /*! caution! the following mass matrix may not be assembled as M[1-2v]+2M[v^2]! This gives instabilities! Instead: M[(1-v)^2]+M[v^2] */
      /*! caution! in the operators "qc::SquaredInvDeformedWeightMassOp" the nodal interpolant \f$I[\phi^{-1}]\f$ is computed and then during quadrature
       *           concatenated with v. If instead the nodal interpolant \f$I[v\circ\phi^{-1}]\f$ were computed and then used, the algorithm will
       *           produce strong artificial stresses at the phase field! */
      qc::SquaredWeightDeformMassOp<ConfiguratorType>( Grid, Displacements[i], localAux ).assembleAddMatrix( localSystemMatrix );
      qc::SquaredWeightDeformMassOp<ConfiguratorType>( Grid, Displacements[i], Phasefields[i] ).assembleAddMatrix( localSystemMatrix );
      qc::SquaredWeightDeformMassOp<ConfiguratorType>( Grid, Displacements[i], Phasefields[i] ).apply( ones, localRhs );
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        rhs += localRhs;
        systemMatrix += localSystemMatrix;
      }
    }
    aol::MassOp<ConfiguratorType>( Grid ).assembleAddMatrix( systemMatrix, Phasefields.numComponents() * Mu / ( 4. * Gamma ) );
    aol::StiffOp<ConfiguratorType>( Grid ).assembleAddMatrix( systemMatrix, Phasefields.numComponents() * aol::Sqr( Epsilon ) * Mu / Gamma );
    aux.setAll( Phasefields.numComponents() * Mu / ( 4. * Gamma ) );
    aol::MassOp<ConfiguratorType>( Grid ).applyAdd( aux , rhs );

    //! solve for the average shape
    solveLinearEquation( systemMatrix, rhs, Average );
  }

  void findMinimizingDisplacement( const typename ConfiguratorType::InitType &Grid,
                                   const aol::Vector<RealType> &Image,
                                   const aol::Vector<RealType> &Phasefield,
                                   const aol::Vector<RealType> &Average,
                                   const aol::MultiVector<RealType> &Stress,
                                   aol::MultiVector<RealType> &Displacement,
                                   const RealType LengthEnergyWeight,
                                   const RealType VolumeEnergyWeight,
                                   const RealType Gamma,
                                   const RealType Epsilon,
                                   const int NumSteps,
                                   const RealType RigidnessOfMaterialOutside = 1e-3 ) {
    aol::Vector<RealType> weight( Image, aol::STRUCT_COPY );
    if ( _constantElasticParameters )
      weight.setAll( 1. );
    else
      dilateAndClamp<ConfiguratorType>( Image, Grid, weight, 3.5, RigidnessOfMaterialOutside, 1. );
    aol::MultiVector<RealType> initialDisp( Displacement, aol::DEEP_COPY );
    if ( _imageMatching ) {
      //! define total energy "E" and its variation with respect to the displacement, "DE" (can also be adapted for image averaging)
      EnergyOfDDefault<ConfiguratorType,SquaredDifference> E( Grid, Image, Average, LengthEnergyWeight, VolumeEnergyWeight, Gamma / Epsilon, weight, Stress );
      EnergyVariationWRTDDefault<ConfiguratorType,SquaredDifference> DE( Grid, Image, Average, LengthEnergyWeight, VolumeEnergyWeight, Gamma / Epsilon, weight, Stress );

      //! perform a few steps of gradient descent
      aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > gradientDescent( Grid, E, DE, NumSteps << ( _finestLevel - Grid.getGridDepth() ) );
      gradientDescent.apply( initialDisp, Displacement );
    } else {
      //! define total energy "E" and its variation with respect to the displacement, "DE" (can also be adapted for image averaging)
      EnergyOfDDefault<ConfiguratorType,PhasefieldMismatch> E( Grid, Phasefield, Average, LengthEnergyWeight, VolumeEnergyWeight, Gamma / Epsilon, weight, Stress );
      EnergyVariationWRTDDefault<ConfiguratorType,PhasefieldMismatch> DE( Grid, Phasefield, Average, LengthEnergyWeight, VolumeEnergyWeight, Gamma / Epsilon, weight, Stress );

      //! perform a few steps of gradient descent
      aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > gradientDescent( Grid, E, DE, NumSteps << ( _finestLevel - Grid.getGridDepth() ) );
      // with automatic filter width the gradient descent does too few iterations and hence does not properly find rotations
      //aol::GradientDescentWithAutomaticFilterWidth<ConfiguratorType,aol::MultiVector<RealType> > gradientDescent( Grid, E, DE, fac/*NumSteps*/ );
      //gradientDescent.setMaximumNumIterationsPerFilterWidth( 10 );
      //gradientDescent.setStopFilterWidth( 1.e-20 );
      //gradientDescent.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::USE_NONLINEAR_CG );
      gradientDescent.apply( initialDisp, Displacement );
    }
  }

  void findMinimizingPhasefield( const typename ConfiguratorType::InitType &Grid,
                                 const aol::Vector<RealType> &Image,
                                 const aol::Vector<RealType> &AveragePhasefield,
                                 const aol::MultiVector<RealType> &Displacement,
                                 aol::Vector<RealType> &Phasefield,
                                 const RealType Gamma,
                                 const RealType Beta,
                                 const RealType Nu,
                                 const RealType Epsilon ) {
    //! assemble the Ambrosio-Tortorelli part of the system matrix
    typename ConfiguratorType::MatrixType systemMatrix( Grid );
    aol::MassOp<ConfiguratorType>( Grid ).assembleAddMatrix( systemMatrix, Nu / ( 4. * Epsilon ) );
    aol::StiffOp<ConfiguratorType>( Grid ).assembleAddMatrix( systemMatrix, Nu * Epsilon );
    aol::SquaredDiffWeightMassOp<ConfiguratorType>( Grid, Image ).assembleAddMatrix( systemMatrix, Beta / 2. );

    //! assemble the averaging mismatch part of the system matrix
    typename ConfiguratorType::ArrayType average( AveragePhasefield, Grid, aol::DEEP_COPY );
    average.addToAll( -1. );
    qc::SquaredDeformedWeightMassOp<ConfiguratorType>( Grid, average, Displacement ).assembleAddMatrix( systemMatrix, Gamma / Epsilon );
    qc::SquaredDeformedWeightMassOp<ConfiguratorType>( Grid, AveragePhasefield, Displacement ).assembleAddMatrix( systemMatrix, Gamma / Epsilon );

    //! assemble the rhs
    aol::Vector<RealType> rhs( Grid ), aux( Grid );
    aux.setAll( Nu / ( 4. * Epsilon ) );
    aol::MassOp<ConfiguratorType>( Grid ).apply( aux , rhs );
    aux.setAll( Gamma / Epsilon );
    qc::SquaredDeformedWeightMassOp<ConfiguratorType>( Grid, AveragePhasefield, Displacement ).applyAdd( aux, rhs );

    //! solve for the phasefield
    solveLinearEquation( systemMatrix, rhs, Phasefield );
  }

  void saveResults( const typename ConfiguratorType::InitType &Grid,
                    const aol::Vector<RealType> &Average,
                    const aol::MultiVector<RealType> &Images,
                    const aol::MultiVector<RealType> &/*Phasefields*/,
                    const aol::RandomAccessContainer<aol::MultiVector<RealType> > &Displacements,
                    const char* DestDirectory,
                    const int Mode,
                    const int Index1,
                    const int Index2 ) {
    char imageFileExtension[4]; // last character will contain the zero-byte, signifying the end of the string
      sprintf( imageFileExtension, "%s", ( ConfiguratorType::Dim == 2 ) ? "pgm" : "bz2" );
    char fileName[1024], prefix[1024];
    if ( Mode > 0 )
      sprintf( prefix, "%s/mode_%d_", DestDirectory, Mode );
    else
      sprintf( prefix, "%s/", DestDirectory );

    //! save average shape
    ArrayType average( Average, Grid, aol::DEEP_COPY );
    average *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
    sprintf( fileName, "%saverageShape_%d_%d.%s", prefix, Index1, Index2, imageFileExtension );
    average.save( fileName, ( ConfiguratorType::Dim == 2 ) ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );

    //! save the displacements
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Images.numComponents(); i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        char fileName[1024];
        sprintf( fileName, "%sdisplacement_%d_%d_%d.bz2", prefix, Index1, i, j );
        ( ArrayType( Displacements[i][j], Grid, aol::FLAT_COPY ) ).save( fileName, qc::PGM_DOUBLE_BINARY );
      }

    //! save the inversely deformed objects
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Images.numComponents(); i++ ) {
      char fileName[1024];
      ArrayType invDefU( Grid );
      qc::InvDeformImage<ConfiguratorType,aol::Vector<RealType> >( Images[i], Grid, invDefU, Displacements[i] );
      invDefU *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
      sprintf( fileName, "%s/deformedImage_%d_%d_%d.%s", DestDirectory, Index1, Index2, i, imageFileExtension );
      invDefU.save( fileName, ( ConfiguratorType::Dim == 2 ) ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );
    }

    //! save the first inversely deformed object as average image
    if ( Mode == 0 ) {
      char fileName[1024];
      ArrayType invDefU( Grid );
      qc::InvDeformImage<ConfiguratorType,aol::Vector<RealType> >( Images[0], Grid, invDefU, Displacements[0] );
      invDefU *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
      sprintf( fileName, "%s/averageImage.%s", DestDirectory, imageFileExtension );
      invDefU.save( fileName, ( ConfiguratorType::Dim == 2 ) ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );
    }

    /*//! save the inversely deformed phasefields
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Images.numComponents(); i++ ) {
      char fileName[1024];
      ArrayType invDefPhasefield( Grid );
      qc::InvDeformImage<ConfiguratorType,aol::Vector<RealType> >( Phasefields[i], Grid, invDefPhasefield, Displacements[i] );
      invDefPhasefield *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
      sprintf( fileName, "%s/deformedPhasefield_%d_%d_%d.%s", DestDirectory, Index1, Index2, i, imageFileExtension );
      invDefPhasefield.save( fileName, ( ConfiguratorType::Dim == 2 ) ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );
    }*/

    //! save the displacement of the average shape according to a mode of variation
    if ( Mode > 0 ) {
      aol::MultiVector<RealType> linDisp( Grid );
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < Images.numComponents(); i++ ) {
        char fileName[1024];
        aol::MultiVector<RealType> displacement( Grid );
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          sprintf( fileName, _deformationFileNameTemplate, Index1, i, j );
          ArrayType( displacement[j], Grid, aol::FLAT_COPY ).load( fileName );
        }
        displacement -= Displacements[i];
        aol::MultiVector<RealType> auxDisp( Grid );
        qc::InvDeformImage<ConfiguratorType,aol::MultiVector<RealType> >( displacement, Grid, auxDisp, Displacements[i] );
#ifdef _OPENMP
#pragma omp critical
#endif
        { linDisp += auxDisp; }
      }
      linDisp /= Images.numComponents();
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        sprintf( fileName, "%scommonDisplacement_%d_%d.bz2", prefix, Index1, j );
        ( ArrayType( linDisp[j], Grid, aol::FLAT_COPY ).save( fileName, qc::PGM_DOUBLE_BINARY ) );
      }
    }
  }

  void execute() {
    //! on each level of the multilevel method do...
    for ( int level = _coarsestLevel; level <= _finestLevel; level++ ) {
      cerr << endl << aol::color::red << "level " << level << " of " << _finestLevel << aol::color::reset << endl;

      // refine mesh
      typename ConfiguratorType::InitType grid( level, ConfiguratorType::Dim );
      _images.setCurLevel( level );
      _inputImages.setCurLevel( level );
      _phasefields.setCurLevel( level );
      _average.setCurLevel( level );
      _displacements.setCurLevel( level );
      _stresses.setCurLevel( level );
      _epsilon = grid.H() * _epsFactor;

      // put all variables into a more convenient form
      aol::MultiVector<RealType> images, phasefields;
      aol::RandomAccessContainer<aol::MultiVector<RealType> > displacements, stresses;
      _displacements.getMultiArrays( displacements );
      _stresses.getMultiArrays( stresses );
      _images.appendMultiArrayReferenceTo( 0, images );
      _phasefields.appendMultiArrayReferenceTo( 0, phasefields );

      //! solve for the phasefields of the input images (may be left out, if only average image is to be found)
      if ( (!_jointSegmentationAveraging) && (!_imageMatching) ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int i = 0; i < _numberOfImages; i++ ) {
          AmbrosioTortorelliPhasefieldOp<ConfiguratorType> phasefieldOp( grid, _beta, _nu, _epsilon );
          phasefieldOp.apply( images[i], phasefields[i] );
          cerr << "phase field " << i+1 << " of " << _numberOfImages << " updated" << endl;
        }
      }
    
      //! resolve self-penetration
      if( !_updateShapeAverage ){
        cerr << "Resolve self-penetration on level " << level << ": " << endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int i = 1; i < _numberOfImages; i++ )
          resolveSelfPenetration<ConfiguratorType>( displacements[i], grid, 100, 0, true );	
      }

      //! perform optimisation steps on this level
      for ( int iteration = 0; iteration < _maxIterations; iteration++ ) {
        cerr << aol::color::green << "iteration " << iteration + 1 << " of " << _maxIterations << aol::color::reset << endl;

        //! solve for the input phasefields (only for joint segmentation and averaging)
        if ( _jointSegmentationAveraging )
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for ( int i = 0; i < _numberOfImages; i++ ) {
            if ( _imageSmoothing )
              AmbrosioTortorelliImageSmootherOp<ConfiguratorType>( grid, phasefields[i], _alpha, _beta ).apply( _inputImages.getArray( i ), images[i] );
            if ( !_imageMatching )
              findMinimizingPhasefield( grid, images[i], _average.current(), displacements[i], phasefields[i], _gamma, _beta, _nu, _epsilon );
            cerr << "phase field " << i+1 << " of " << _numberOfImages << " updated" << endl;
          }

        //! solve for the average shape (can alternatively be done for the average image)
        if( _updateShapeAverage ){
          if ( _imageMatching )
            findMinimizingAverageImage( grid, images, displacements, _average.current(), _gamma, _mu, _epsilon );
          else
            findMinimizingAveragePhasefield( grid, phasefields, displacements, _average.current(), _gamma, _mu, _epsilon );
          cerr << "average shape updated" << endl;
	}

        //! approximately solve for the displacements
        const RealType rigidnessOutside[10] = { 1, 1, 1, .5, .5, .1, .05, .01, .005, .001 };
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int i = 0; i < _numberOfImages; i++ ) {
          findMinimizingDisplacement( grid, images[i], phasefields[i], _average.current(), stresses[i], displacements[i], _shapeWeights[i] * _lengthEnergyWeight, _shapeWeights[i] * _volumeEnergyWeight, _gamma, _epsilon, _maxSteps, level > _coarsestLevel ? 1e-3 : rigidnessOutside[iteration] );
          cerr << "displacement " << i+1 << " of " << _numberOfImages << " updated" << endl;
        }

        //! save data: if shape average is not updated, we only save the displacements
        if( _updateShapeAverage )
          saveResults( grid, _average.current(), images, phasefields, displacements, _destDirectory, _mode, level, iteration );
	else{
	  //cerr << "Save displacements on level " << level << ": " << endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for ( int i = 0; i < displacements.size(); i++ )
            for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
              char fileName[1024];
              sprintf( fileName, _deformationFileNameTemplate, i, j );
              ( ArrayType( displacements[i][j], grid, aol::FLAT_COPY ) ).save( fileName, qc::PGM_DOUBLE_BINARY );
            }
	}
      }
  
      // check energy contributions
      for ( int i = 0; i < _numberOfImages; i++ )
	evaluateEnergy( grid, images[i], _average.current(), displacements[i], _shapeWeights[i] * _lengthEnergyWeight, _shapeWeights[i] * _volumeEnergyWeight, _gamma, _epsilon, 1e-3, i );

      // prolongate the results of this level
      _displacements.levProlongate();
      if( _updateShapeAverage )
        _average.levProlongate();      
    }
  }

protected:
  void evaluateEnergy( const typename ConfiguratorType::InitType &Grid,
                                   const aol::Vector<RealType> &Image,
                                   const aol::Vector<RealType> &Average,
                                   aol::MultiVector<RealType> &Displacement,
                                   const RealType LengthEnergyWeight,
                                   const RealType VolumeEnergyWeight,
                                   const RealType Gamma,
                                   const RealType Epsilon,
                                   const RealType RigidnessOfMaterialOutside,
				   int number ) const {
    aol::Vector<RealType> weight( Image, aol::STRUCT_COPY );
    if ( _constantElasticParameters )
      weight.setAll( 1. );
    else
      dilateAndClamp<ConfiguratorType>( Image, Grid, weight, 3.5, RigidnessOfMaterialOutside, 1. );
    
    typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > WeightedEnergyDensityType;
    typedef qc::HyperelasticEnergy<ConfiguratorType,WeightedEnergyDensityType> WeightedHyperelasticEnergyType;
    qc::HyperelasticEnergyDensityDefault<ConfiguratorType> energyDensity( LengthEnergyWeight, 0., VolumeEnergyWeight );
    WeightedEnergyDensityType weightedDensity( energyDensity, Grid, weight );
    WeightedHyperelasticEnergyType hyperelasticEnergy( Grid, weightedDensity );

    // the energy term itself
    EnergyOfD<ConfiguratorType, WeightedHyperelasticEnergyType, SquaredDifference> energyOfD( Grid, Image, Average, hyperelasticEnergy, Gamma / Epsilon );
    //EnergyOfDDefault<ConfiguratorType,SquaredDifference> E( Grid, Image, Average, LengthEnergyWeight, VolumeEnergyWeight, Gamma / Epsilon, weight );
    aol::Scalar<RealType> energy;
    energyOfD.apply( Displacement, energy );    
    
    cerr << number << "th displacement has energy  E[phi] = " << energy[0] << " = " << energyOfD.getLastElasticEnergy( ) << " + " << Gamma / Epsilon << " * " << energyOfD.getLastMismatchEnergy( false ) << ", with rel. mismatch = " << Average.size() * energyOfD.getLastMismatchEnergy( false ) / Average.sum()  << endl;    
  }

};

template <typename ConfiguratorType>
class DisplacementVariationOp {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  // grid over the examined domain
  const typename ConfiguratorType::InitType _grid;
  // number of input objects
  const int _numberOfImages;
  aol::MultiVector<RealType> _weights;
  aol::RandomAccessContainer<aol::MultiVector<RealType> > _averageDisp;
  // the directory for the results
  char _destDirectory[1024];
  // weights of the hyperelastic deformation energy
  const RealType _lengthEnergyWeight, _volumeEnergyWeight;
  RealType _epsilon, _epsFactor;
  const bool _constantElasticParameters;
  // average shape (as binary image or phase field)
  ArrayType _average, _averageShape;
  // displacements
  aol::VectorContainer<aol::MultiVector<RealType> > _displacements;
  // visualization flags
  const bool _overlaidVisualization, _visualizationForShape, _maskedDeformationForVisualization;
  const RealType _deformationAmplification;

public:
  DisplacementVariationOp( aol::ParameterParser &Parser ) :
    // define the finest grid over the domain
    _grid( Parser.getInt( "GridDepth" ), ConfiguratorType::Dim ),
    // load number of objects and create the object-vector and shape vector
    _numberOfImages( Parser.getInt( "numberOfImages" ) ),
    _weights( _numberOfImages, _grid.getNumberOfNodes() ),
    _averageDisp( _numberOfImages ),
    // load hyperelastic and phasefield parameters
    _lengthEnergyWeight( Parser.getDouble( "lengthEnergyWeight" ) ),
    _volumeEnergyWeight( Parser.getDouble( "volumeEnergyWeight" ) ),
    _epsFactor( Parser.getDouble( "epsFactor" ) ),
    _constantElasticParameters( Parser.getInt( "constantElasticParameters" ) ),
    // create the average shape and the displacements
    _average( Parser.getString( "averageImageFile" ) ),
    _averageShape( Parser.getString( "averageShapeFile" ) ),
    _displacements( _numberOfImages ),
    _overlaidVisualization( Parser.hasVariable( "overlaidVisualization" ) ? Parser.getInt( "overlaidVisualization" ) == 1 : ConfiguratorType::Dim == qc::QC_2D ),
    _visualizationForShape( Parser.hasVariable( "visualizationForShape" ) ? Parser.getInt( "visualizationForShape" ) == 1 : ConfiguratorType::Dim == qc::QC_2D ),
    _maskedDeformationForVisualization( Parser.hasVariable( "maskedDeformationForVisualization" ) && Parser.getInt( "maskedDeformationForVisualization" ) == 1 ),
    _deformationAmplification( Parser.hasVariable( "deformationAmplification" ) ? Parser.getDouble( "deformationAmplification" ) : 1. ){

    // read in the directory, where results are to be saved
    Parser.getString( "destDirectory", _destDirectory );

    // scale the average image to the interval [0,1]
    _average.addToAll( -_average.getMinValue() );
    _average /= _average.getMaxValue();
    _averageShape.addToAll( -_averageShape.getMinValue() );
    _averageShape /= _averageShape.getMaxValue();

    // load the original images to obtain elastic weights
    if ( _constantElasticParameters )
      _weights.setAll( 1. );
    else
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < _numberOfImages; i++ ) {
        char fileName[1024];
        sprintf( fileName, Parser.getString( "imageFilenameTemplate" ).c_str(), i + 1 );
        cerr<<"load " + string( fileName ) + "...\n";
        ArrayType image( fileName );
        dilateAndClamp<ConfiguratorType>( image, _grid, _weights[i] );
      }

    // load the displacements to the average
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numberOfImages; i++ ) {
      _averageDisp[i].reallocate( _grid );
      if ( Parser.hasVariable( "averageDisplacementFilenameTemplate" ) )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          char fileName[1024];
          sprintf( fileName, Parser.getString( "averageDisplacementFilenameTemplate" ).c_str(), i, j );
          cerr<<"load " + string( fileName ) + "...\n";
          ArrayType( _averageDisp[i][j], _grid, aol::FLAT_COPY ).load( fileName );
        }
    }

    // compute the displacements of the average
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int mode = 0; mode < _numberOfImages; mode++ ) {
      _displacements[mode].reallocate( _grid );
      if ( Parser.hasVariable( "modeDisplacementFilenameTemplate" ) ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int i = 0; i < _numberOfImages; i++ ) {
          // load the displacement of the ith object for the modeth mode
          aol::MultiVector<RealType> displacement( _grid ), auxDisp( _grid );
          char fileName[1024];
          for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
            sprintf( fileName, Parser.getString( "modeDisplacementFilenameTemplate" ).c_str(), mode, i, j );
            cerr<<"load " + string( fileName ) + "...\n";
            ArrayType( displacement[j], _grid, aol::FLAT_COPY ).load( fileName );
          }
          // compute the displacement of the average
          displacement -= _averageDisp[i];
          qc::InvDeformImage<ConfiguratorType,aol::MultiVector<RealType> >( displacement, _grid, auxDisp, _averageDisp[i] );
#ifdef _OPENMP
#pragma omp critical
#endif
          { _displacements[mode] += auxDisp; }
        }
        _displacements[mode] /= _numberOfImages;
      } else {
        char fileName[1024];
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          sprintf( fileName, Parser.getString( "displacementFilenameTemplate" ).c_str(), mode + 1, j );
          cerr<<"load " + string( fileName ) + "...\n";
          ArrayType( _displacements[mode][j], _grid, aol::FLAT_COPY ).load( fileName );
        }
      }
    }

    // notify the user
    cerr << "Parameters loaded." << endl;
  }

  void saveResults( const typename ConfiguratorType::InitType &Grid,
                    const aol::VectorContainer<aol::MultiVector<RealType> > &Displacements,
                    const aol::Vector<RealType> &Variances,
                    const ArrayType &Average,
                    const ArrayType &AverageShape,
                    const char* DestDirectory,
		    RealType factor  ) {
    char imageFileExtension[4]; // last character will contain the zero-byte, signifying the end of the string
      sprintf( imageFileExtension, "%s", ( ConfiguratorType::Dim == 2 ) ? "pgm" : "bz2" );

    //! fix, whether visualization shall be done for the average image or the average shape
    ArrayType average( Average, aol::DEEP_COPY );
    if ( _visualizationForShape )
      average = AverageShape;
    if ( ConfiguratorType::Dim == qc::QC_2D )
      average *= 255;

    //! generate a mask for the deformation (if only interior object shall be deformed)
    qc::BitArray<ConfiguratorType::Dim> mask( qc::GridSize<ConfiguratorType::Dim>::createFrom( Grid ) );
    ArrayType image( Average, aol::DEEP_COPY );
    if ( _maskedDeformationForVisualization ) {
      ArrayType shapeSignedDist( Average, aol::STRUCT_COPY );
      image.addToAll( -.5 );
      ( typename qc::SignedDistanceOpTrait<ConfiguratorType,ConfiguratorType::Dim>::OpType( _grid ) ).apply( image, shapeSignedDist );
      shapeSignedDist.addToAll( -0.5/*3.5*/ * _epsFactor * Grid.H() );
      for ( int i = 0; i < shapeSignedDist.size(); i++ )
        mask.set( i, shapeSignedDist[i] > 0 );
    }
    else
      mask.setAll( true );
    image.setAll( ( ConfiguratorType::Dim == qc::QC_2D ) ? 255 : 0. );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Displacements.size(); i++ ) {
      //! save the displacement
      char fileName[1024];
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        sprintf( fileName, "%sdisplacementMode_%d_%d.bz2", DestDirectory, i, j );
        ( ArrayType( Displacements[i][j], Grid, aol::FLAT_COPY ) ).save( fileName, qc::PGM_DOUBLE_BINARY );
      }

      //! save the average, deformed according to the displacement mode
      ArrayType invDefU( Grid );
      aol::MultiVector<RealType> displacement( Displacements[i], aol::DEEP_COPY );
      displacement *= factor;
      qc::InvDeformImage<ConfiguratorType,aol::Vector<RealType> >( average, Grid, invDefU, displacement );
      sprintf( fileName, "%s/deformedImageMode_%d.%s", DestDirectory, i, imageFileExtension );
      invDefU.save( fileName, ( ConfiguratorType::Dim == qc::QC_2D ) ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );

      //! visualize the displacement modes
      ArrayType overlaidImage( Grid );
      // produce deformations of different strengths
      for ( int j = -3; j < 4; j++ ) {
        aol::MultiVector<RealType> scaledDispl( displacement, aol::DEEP_COPY );
        scaledDispl *= j * .8;
        qc::TransformFunction<RealType,ConfiguratorType::Dim> transformOp( Grid, &mask );
        transformOp.setDeformation( scaledDispl );
        qc::GridSize<ConfiguratorType::Dim> gridSize( Grid );
        qc::BitArray<ConfiguratorType::Dim> aux( gridSize );
        if ( _overlaidVisualization ) {
        // in 2d, produce overlaid phase fields, in 3d save the deformed objects
          transformOp.transform( average, invDefU, aux, image );
          //qc::InvDeformImage<ConfiguratorType,aol::Vector<RealType> >( AverageShape, Grid, invDefU, scaledDispl );
          invDefU *= ( j == 0 ) ? 5 : 1;
          overlaidImage += invDefU;
        } else {
          transformOp.transform( average, invDefU, aux, image );
          sprintf( fileName, "%s/shapeMode_%d_%d.%s", DestDirectory, i, j, imageFileExtension );
          invDefU.save( fileName, ConfiguratorType::Dim == qc::QC_2D ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );
        }
      }
      if ( _overlaidVisualization ) {
        overlaidImage *= ( ConfiguratorType::Dim == qc::QC_2D ? 255 : 1 ) / overlaidImage.getMaxValue();
        sprintf( fileName, "%s/shapeMode_%d.%s", DestDirectory, i, imageFileExtension );
        overlaidImage.save( fileName, ConfiguratorType::Dim == qc::QC_2D ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );
      }
    }

    //! save the variances
    char fileName[1024];
    sprintf( fileName, "%s/variances.dat", DestDirectory );
    Variances.saveAsVector( fileName );
    sprintf( fileName, "%s/variances.txt", DestDirectory );
    Variances.saveASCII( fileName );
    cerr<<"variances "<<Variances<<endl;
  }

  void execute() {
    // compute the weight
    aol::Vector<RealType> weight( _average, aol::DEEP_COPY );
    if ( _constantElasticParameters )
      weight.setAll( 1. );
    else {
      /*// smoothly extend the displacements outside the shape
      qc::BitArray<ConfiguratorType::Dim> valueDefined( _grid );
      for ( int k = 0; k < weight.size(); ++k )
        valueDefined.set( k, weight[k] > 1e-2 );
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < _displacements.size(); i++ ) {
        aol::MultiVector<RealType> displacement( _displacements[i], aol::DEEP_COPY );
        qc::SmoothlyExtendImage<ConfiguratorType>( displacement, _grid, _displacements[i], valueDefined );
      }*/
    }

#if 0 // this should have happened already when computing the displacements
    // make the mean shift and moment of the displacement zero
    qc::meanZeroShiftAndRotationProjector<ConfiguratorType> projector( _grid, weight, true );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _displacements.size(); i++ ) {
      projector.apply( _displacements[i], _displacements[i] );
      // check, whether moment and mean are really zero 2D
      cerr<<"mean shift of displacement "<<i<<": "<<projector.totalShift( _displacements[i] )<<";  "
          <<"mean moment of displacement "<<i<<": "<<projector.totalMoment( _displacements[i] )<<endl;
    }
#endif

    // compute object volume, \int_O 1 dx
    aol::Vector<RealType> ones( _grid.getNumberOfNodes() ), aux( ones, aol::STRUCT_COPY );
    ones.setAll( 1. );
    aol::MassOp<ConfiguratorType>( _grid ).apply( weight, aux );
    RealType volume = aux * ones;

    /*// alternative 1: define L2 product as inner product
    aol::WeightedMassOp<ConfiguratorType> l2Prod( _grid, weight, aol::ASSEMBLED );
    aol::BlockOp<RealType> dotProd( ConfiguratorType::Dim, ConfiguratorType::Dim );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      dotProd.setReference( j, j, l2Prod );*/

    // alternative 2:define elastic inner product
    const aol::MultiVector<RealType> zeros( _grid );
    typedef typename ConfiguratorType::MatrixType SubMatrixType;
    typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > WeightedEnergyDensityType;
    const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> energyDensity( _lengthEnergyWeight, 0, _volumeEnergyWeight );
    aol::BlockMatrix<SubMatrixType> dotProd( _grid );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numberOfImages; i++ ) {
      cerr<<"assemble stiffness matrix of "<<i<<"th deformed object..."<<endl;
      aol::BlockMatrix<SubMatrixType> localMatrix( _grid );
      _weights[i] *= 1. / _numberOfImages;
      WeightedEnergyDensityType weightedDensity( energyDensity, _grid, _weights[i] );
      qc::HyperelasticDeformHessian<ConfiguratorType, WeightedEnergyDensityType,SubMatrixType>( _grid, weightedDensity, _averageDisp[i] ).applyAdd( zeros, localMatrix );
#ifdef _OPENMP
#pragma omp critical
#endif
      dotProd += localMatrix;
    }

    // perform the PCA 
    // caution: covariance operator does not divide by (number of samples - 1) !
    aol::VectorContainer<aol::MultiVector<RealType> > modes;
    aol::Vector<RealType> variances;
    aol::PCAOp<RealType,aol::MultiVector<RealType> > pca( dotProd, _numberOfImages );
    pca.modesAndVariances( _displacements, modes, variances );
    //variances /= volume;


    cerr << "Test Mahalanobis distance with the original displacements:" << endl;
    for( int i = 0; i < _numberOfImages; i++ ){
      aol::Vector<RealType> distances( _numberOfImages - 1 );
      distances.setZero();            
      aol::MultiVector<RealType> temp( _displacements[i], aol::STRUCT_COPY );
      dotProd.apply( _displacements[i], temp );
      for( int j = 0; j < _numberOfImages - 1; j++ ){
	RealType aux = temp * modes[j];
        distances[j] = aux * aux / variances[j];
      }
      // smallest eigenvalue is zero!
      cerr << "sum of " << i << "th mode = " << distances << endl;
      cerr << "Mahalanobis distance = " << sqrt( distances.sum() )<< endl << endl;
    }

    // save the results
    saveResults( _grid, modes, variances, _average, _averageShape, _destDirectory, _deformationAmplification * .01 * sqrt( volume ) );
  }
};

template<typename ConfiguratorType>
class StressTest {

  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType _grid;

public:
  StressTest( const int GridDepth ) :
    _grid( GridDepth, ConfiguratorType::Dim ) {}

  void execute() {
    aol::MultiVector<RealType> displacement( _grid ), stress( ConfiguratorType::Dim*ConfiguratorType::Dim, _grid.getNumberOfNodes() );
    qc::DataGenerator<ConfiguratorType>( _grid ).generateIdentity( displacement );
    displacement *= -.01;
    qc::HyperelasticEnergyDensityDefault<ConfiguratorType> density( 1., 0., 1. );
    FirstPiolaKirchhoffStressOp<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> >( _grid, density ).apply( displacement, stress );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        cerr<<"stress "<<i<<j<<". min: "<<stress[i*ConfiguratorType::Dim+j].getMinValue()<<". max: "<<stress[i*ConfiguratorType::Dim+j].getMaxValue()<<"."<<endl;
    displacement.setZero();

    aol::RandomAccessContainer<aol::MultiVector<RealType> > preStressDisp;
    preStressDisp.pushBack( displacement );
    preStressDisp[0].setZero();
    aol::MultiVector<RealType> weight( 1, _grid.getNumberOfNodes() );
    weight.setAll( 1. );
    CommonDisplacementEnergy<ConfiguratorType> E( _grid, preStressDisp, stress, 1., 1., weight );
    CommonDisplacementEnergyVariation<ConfiguratorType> DE( _grid, preStressDisp, stress, 1., 1., weight );

    //! perform a few steps of gradient descent
    aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > gradientDescent( _grid, E, DE, 4000 );
    gradientDescent.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::USE_NONLINEAR_CG );
    aol::MultiVector<RealType> initialDisp( displacement, aol::STRUCT_COPY );
    /*qc::DataGenerator<ConfiguratorType>( _grid ).generateIdentity( initialDisp );
    initialDisp *= -.01;*/
    gradientDescent.apply( initialDisp, displacement );

/*//if buckling modes are to be examined and rotation shall be prevented, this piece of code has to be included in CommonDisplacementEnergy
typename ConfiguratorType::ArrayType arg1( Arg[0], _grid, aol::FLAT_COPY), arg2( Arg[1], _grid, aol::FLAT_COPY);
for ( int x = 0; x < arg1.getNumX() / 2; x++ )
for ( int y = 0; y < arg1.getNumY() / 2; y++ ) {
arg1.set(arg1.getNumX()-x-1,y,-arg1.get(x,y));
arg1.set(arg1.getNumX()-x-1,arg1.getNumY()-y-1,-arg1.get(x,y));
arg1.set(x,arg1.getNumY()-y-1,arg1.get(x,y));
arg2.set(arg2.getNumX()-x-1,y,arg2.get(x,y));
arg2.set(arg2.getNumX()-x-1,arg2.getNumY()-y-1,-arg2.get(x,y));
arg2.set(x,arg2.getNumY()-y-1,-arg2.get(x,y));
}
*/

    // save the resulting displacement
    displacement[0].addToAll( -displacement[0].getMinValue() );
    displacement[1].addToAll( -displacement[1].getMinValue() );
    saveDisplacement<ConfiguratorType>( displacement, _grid, "testDisp" );

    // save the local energy distribution
    aol::MultiVector<RealType> energy( 1, _grid.getNumberOfNodes() );
    qc::HyperelasticEnergy<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> >( _grid, density ).applyAddIntegrand( displacement, energy );
    save<ConfiguratorType>( energy[0], _grid, "testEnergyTotal" );
    energy.setZero();
    qc::HyperelasticEnergyDensityDefault<ConfiguratorType> densityLength( 1., 0., 0. );
    qc::HyperelasticEnergy<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> >( _grid, densityLength ).applyAddIntegrand( displacement, energy );
    save<ConfiguratorType>( energy[0], _grid, "testEnergyLength" );
    energy.setZero();
    qc::HyperelasticEnergyDensityDefault<ConfiguratorType> densityVolume( 0., 0., 1. );
    qc::HyperelasticEnergy<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> >( _grid, densityVolume ).applyAddIntegrand( displacement, energy );
    save<ConfiguratorType>( energy[0], _grid, "testEnergyVolume" );
  }
};



#endif
