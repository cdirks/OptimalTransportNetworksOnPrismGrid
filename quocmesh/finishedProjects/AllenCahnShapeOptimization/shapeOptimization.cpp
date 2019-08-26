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

#include "../shapeStatistics/testResultSaving.h"
#include <aol.h>
#include <hyperelastic.h>
#include <parameterParser.h>
#include <gradientDescent.h>
#include <Newton.h>
#include <boundaryIntegration.h>
#include <ringBuffer.h>
#include <suiteSparseSolver.h>
#include <trustRegionMethod.h>

// debugging tools
#include <gradientflow.h>
#include <FEOpInterface.h>
#include <gnuplotter.h>
//static bool TEST = false;

// use exactly one of the following four defines
//#define LINEARIZED_ELASTICITY
//#define USE_INTERNAL_ELASTIC_ENERGY
#define USE_EXTERNAL_LOAD_POTENTIAL
//#define USE_ELASTIC_DISSIPATION

//#define DIRECT_SOLVER

template <typename RealType>
class DoubleWellPotential {
public:
  static RealType evaluate( const RealType W ){ return aol::Sqr( 1 - aol::Sqr( W ) ); }
  static RealType evaluateDerivative( const RealType W ){ return 4 * W * ( aol::Sqr( W ) - 1 ); }
};

template <typename ConfiguratorType,typename IntegrandType>
class FENonlinIntegrationOp :
  public aol::FENonlinIntegrationScalarInterface< ConfiguratorType, FENonlinIntegrationOp<ConfiguratorType, IntegrandType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const IntegrandType &_integrand;

public:
  FENonlinIntegrationOp( const typename ConfiguratorType::InitType &Grid, const IntegrandType &Integrand ) :
    aol::FENonlinIntegrationScalarInterface< ConfiguratorType, FENonlinIntegrationOp<ConfiguratorType, IntegrandType> >( Grid ),
    _integrand( Integrand ){}

  RealType evaluateIntegrand( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint,
                              const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    return _integrand.evaluate( DiscFunc.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

template <typename ConfiguratorType,typename IntegrandType>
class FENonlinOp :
  public aol::FENonlinOpInterface< ConfiguratorType, FENonlinOp<ConfiguratorType, IntegrandType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const IntegrandType &_integrand;

public:
  FENonlinOp( const typename ConfiguratorType::InitType &Grid, const IntegrandType &Integrand ) :
    aol::FENonlinOpInterface< ConfiguratorType, FENonlinOp<ConfiguratorType, IntegrandType> >( Grid ),
    _integrand( Integrand ){}

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint,
                        const typename ConfiguratorType::VecType &/*RefCoord*/,
                        typename ConfiguratorType::RealType &NL ) const {
    NL = _integrand.evaluateDerivative( DiscFunc.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class ElasticEnergyWeightDerivative :
  public aol::FENonlinOpInterface< ConfiguratorType, ElasticEnergyWeightDerivative<ConfiguratorType,HyperelasticEnergyDensityType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;
  const HyperelasticEnergyDensityType &_energyDensity;

public:
  ElasticEnergyWeightDerivative( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement, const HyperelasticEnergyDensityType &EnergyDensity ) :
    aol::FENonlinOpInterface< ConfiguratorType, ElasticEnergyWeightDerivative<ConfiguratorType,HyperelasticEnergyDensityType> >( Grid ),
    _displacement( Grid, Displacement ),
    _energyDensity( EnergyDensity ){}

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &/*DiscFunc*/,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint,
                        const typename ConfiguratorType::VecType &/*RefCoord*/,
                        typename ConfiguratorType::RealType &NL ) const {
    // compute the deformation gradient
    typename ConfiguratorType::MatType dphi;
    _displacement.evaluateGradientAtQuadPoint ( El, QuadPoint, dphi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dphi[i][i] += 1.;

    // compute the elastic energy density
    NL = qc::HyperelasticEnergy<ConfiguratorType,HyperelasticEnergyDensityType>::energyDensity( dphi, _energyDensity, El, QuadPoint );
  }
};

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class ElasticLinearizedEnergyWeightDerivative :
  public aol::FENonlinOpInterface< ConfiguratorType, ElasticLinearizedEnergyWeightDerivative<ConfiguratorType,HyperelasticEnergyDensityType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement, _displacementVariation;
  const HyperelasticEnergyDensityType &_energyDensity;

public:
  ElasticLinearizedEnergyWeightDerivative( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement, const aol::MultiVector<RealType> &DisplacementVariation, const HyperelasticEnergyDensityType &EnergyDensity ) :
    aol::FENonlinOpInterface< ConfiguratorType, ElasticLinearizedEnergyWeightDerivative<ConfiguratorType,HyperelasticEnergyDensityType> >( Grid ),
    _displacement( Grid, Displacement ),
    _displacementVariation( Grid, DisplacementVariation ),
    _energyDensity( EnergyDensity ){}

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &/*DiscFunc*/,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint,
                        const typename ConfiguratorType::VecType &/*RefCoord*/,
                        typename ConfiguratorType::RealType &NL ) const {
#ifdef LINEARIZED_ELASTICITY
    // compute the displacement gradient and the one of the displacement variation
    typename ConfiguratorType::MatType dv, du;
    _displacement.evaluateGradientAtQuadPoint ( El, QuadPoint, dv );
    _displacementVariation.evaluateGradientAtQuadPoint ( El, QuadPoint, du );

    // compute the elastic energy density
    NL = 0.;
    typename ConfiguratorType::VecType stress;
    typename ConfiguratorType::MatType tensor, id;
    id.setIdentity();
    for ( int k = 0; k < ConfiguratorType::Dim; k++ )
      for ( int l = 0; l < ConfiguratorType::Dim; l++ ) {
        qc::HyperelasticSubHessian<ConfiguratorType,HyperelasticEnergyDensityType>::elasticityTensor( id, _energyDensity, El, QuadPoint, k, l, tensor );
        stress = tensor * du[l];
        NL += stress * dv[k];
      }
#else
    typename ConfiguratorType::MatType dphi, du, stress;

    // compute the deformation gradient
    _displacement.evaluateGradientAtQuadPoint ( El, QuadPoint, dphi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dphi[i][i] += 1.;

    // compute the deformation variation gradient
    _displacementVariation.evaluateGradientAtQuadPoint ( El, QuadPoint, du );

    // compute the elastic energy density
    qc::HyperelasticGradient<ConfiguratorType,HyperelasticEnergyDensityType>::firstPiolaKirchhoffStress( dphi, _energyDensity, El, QuadPoint, stress );
    NL = stress.dotProduct( du );
#endif
  }
};

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class PhasefieldWeightedHyperelasticEnergyDensity {
private:
  typedef typename ConfiguratorType::RealType RealType;
  const HyperelasticEnergyDensityType &_energy;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _phasefield;
  const RealType _delta;
public:
  PhasefieldWeightedHyperelasticEnergyDensity( const HyperelasticEnergyDensityType &Energy, const typename ConfiguratorType::InitType &Grid, const aol::Vector<RealType> &Phasefield, const RealType Delta ) :
    _energy( Energy ),
    _phasefield( Grid, Phasefield ),
    _delta( Delta ) {}

  inline RealType evaluate ( const RealType I1, const RealType I2, const RealType I3,
                             const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    const RealType w = _phasefield.evaluate( El, RefCoord );
    return ( aol::Sqr( w + 1 ) + _delta ) * _energy.evaluate( I1, I2, I3, El, RefCoord );
  }

  inline RealType evaluateAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                        const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    const RealType w = _phasefield.evaluateAtQuadPoint( El, QuadPoint );
    return ( aol::Sqr( w + 1 ) + _delta ) * _energy.evaluateAtQuadPoint( I1, I2, I3, El, QuadPoint );
  }

  inline aol::Vec3<RealType> evaluateDerivative ( const RealType I1, const RealType I2, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    const RealType w = _phasefield.evaluate( El, RefCoord );
    return ( aol::Sqr( w + 1 ) + _delta ) * _energy.evaluateDerivative( I1, I2, I3, El, RefCoord );
  }

  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    const RealType w = _phasefield.evaluateAtQuadPoint( El, QuadPoint );
    return ( aol::Sqr( w + 1 ) + _delta ) * _energy.evaluateDerivativeAtQuadPoint( I1, I2, I3, El, QuadPoint );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivative ( const RealType I1, const RealType I2, const RealType I3,
                                                            const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    aol::Matrix33<RealType> result( _energy.evaluateSecondDerivative( I1, I2, I3, El, RefCoord ) );
    const RealType w = _phasefield.evaluate( El, RefCoord );
    result *= aol::Sqr( w + 1 ) + _delta;
    return result;
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivativeAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                                                       const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    aol::Matrix33<RealType> result( _energy.evaluateSecondDerivativeAtQuadPoint( I1, I2, I3, El, QuadPoint ) );
    const RealType w = _phasefield.evaluateAtQuadPoint( El, QuadPoint );
    result *= aol::Sqr( w + 1 ) + _delta;
    return result;
  }
};

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class RestrictedDomainHyperelasticEnergyDensity {
private:
  typedef typename ConfiguratorType::RealType RealType;
  const HyperelasticEnergyDensityType &_energy;
  const qc::BitArray<ConfiguratorType::Dim> &_domain;
public:
  RestrictedDomainHyperelasticEnergyDensity( const HyperelasticEnergyDensityType &Energy, const qc::BitArray<ConfiguratorType::Dim> &Domain ) :
    _energy( Energy ),
    _domain( Domain ) {}

  inline RealType evaluate ( const RealType I1, const RealType I2, const RealType I3,
                             const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    if ( _domain.elementTrue( El ) )
      return _energy.evaluate( I1, I2, I3, El, RefCoord );
    else
      return 0.;
  }

  inline RealType evaluateAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                        const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    if ( _domain.elementTrue( El ) )
      return _energy.evaluateAtQuadPoint( I1, I2, I3, El, QuadPoint );
    else
      return 0.;
  }

  inline aol::Vec3<RealType> evaluateDerivative ( const RealType I1, const RealType I2, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    if ( _domain.elementTrue( El ) )
      return _energy.evaluateDerivative( I1, I2, I3, El, RefCoord );
    else
      return aol::Vec3<RealType>( 0., 0., 0. );
  }

  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    if ( _domain.elementTrue( El ) )
      return _energy.evaluateDerivativeAtQuadPoint( I1, I2, I3, El, QuadPoint );
    else
      return aol::Vec3<RealType>( 0., 0., 0. );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivative ( const RealType I1, const RealType I2, const RealType I3,
                                                            const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    if ( _domain.elementTrue( El ) )
      return _energy.evaluateSecondDerivative( I1, I2, I3, El, RefCoord );
    else
      return aol::Matrix33<RealType>( 0., 0., 0., 0., 0., 0., 0., 0., 0. );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivativeAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                                                       const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    if ( _domain.elementTrue( El ) )
      return _energy.evaluateSecondDerivativeAtQuadPoint( I1, I2, I3, El, QuadPoint );
    else
      return aol::Matrix33<RealType>( 0., 0., 0., 0., 0., 0., 0., 0., 0. );
  }
};

template <typename ConfiguratorType,typename qc::BoundaryFaceElementBase<typename ConfiguratorType::RealType>::BoundaryFaceType FaceType = qc::BoundaryFaceElementBase<typename ConfiguratorType::RealType>::Z_UPPER_BOUNDARY>
class BoundaryFaceMassOp :
  public qc::FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,BoundaryFaceMassOp<ConfiguratorType,FaceType> > {

  typedef typename ConfiguratorType::RealType RealType;

public:
  explicit BoundaryFaceMassOp ( const typename ConfiguratorType::InitType &Grid, const aol::OperatorType OpType = aol::ONTHEFLY ) :
    qc::FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,BoundaryFaceMassOp<ConfiguratorType,FaceType> > ( Grid, OpType ) {}

  inline RealType getCoeff ( const qc::BoundaryFaceElement<RealType,ConfiguratorType::Dim> &El, const typename ConfiguratorType::DomVecType &/*Local3DCoord*/ ) const {
    return El.getBoundaryFaceType() == FaceType ? 1. : 0.;
  }
};

#ifdef LINEARIZED_ELASTICITY

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class MechanicalEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  aol::BlockMatrix<typename ConfiguratorType::MatrixType> _stiffnessMatrix;
  const aol::MultiVector<RealType> &_externalEnergy;

public:
  MechanicalEnergy( const typename ConfiguratorType::InitType &Grid, const HyperelasticEnergyDensityType &EnergyDensity, const aol::MultiVector<RealType> &ExternalEnergy ) :
    _stiffnessMatrix( Grid ),
    _externalEnergy( ExternalEnergy ) {
    aol::MultiVector<RealType> zero( Grid );
    qc::HyperelasticHessian<ConfiguratorType,HyperelasticEnergyDensityType,typename ConfiguratorType::MatrixType>( Grid, EnergyDensity ).apply( zero, _stiffnessMatrix );
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector<RealType> aux( Arg, aol::STRUCT_COPY );
    _stiffnessMatrix.apply( Arg, aux );
    Dest += .5 * ( aux * Arg );
    Dest -= _externalEnergy * Arg;
  }

  void apply( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector<RealType> aux( Arg, aol::STRUCT_COPY );
    _stiffnessMatrix.apply( Arg, aux );
    Dest = .5 * ( aux * Arg );
    Dest -= _externalEnergy * Arg;
  }
};

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class MechanicalDerivative :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  aol::BlockMatrix<typename ConfiguratorType::MatrixType> _stiffnessMatrix;
  const aol::MultiVector<RealType> &_externalEnergyDeriv;
  const aol::BitVector &_dirichletNodes;

public:
  MechanicalDerivative( const typename ConfiguratorType::InitType &Grid, const HyperelasticEnergyDensityType &EnergyDensity, const aol::MultiVector<RealType> &ExternalEnergy, const aol::BitVector &DirichletNodes ) :
    _stiffnessMatrix( Grid ),
    _externalEnergyDeriv( ExternalEnergy ),
    _dirichletNodes( DirichletNodes ) {
    aol::MultiVector<RealType> zero( Grid );
    qc::HyperelasticHessian<ConfiguratorType,HyperelasticEnergyDensityType,typename ConfiguratorType::MatrixType>( Grid, EnergyDensity ).apply( zero, _stiffnessMatrix );
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> dest( Dest, aol::STRUCT_COPY );
    apply( Arg, dest );
    Dest += dest;
  }

  void apply( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    _stiffnessMatrix.apply( Arg, Dest );
    Dest -= _externalEnergyDeriv;
    // set values belonging to Dirichlet-nodes to zero
    for ( int j = 0; j < Dest[0].size(); j++ )
      if ( _dirichletNodes[j] )
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          Dest[i][j] = 0.;
  }
};

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType, typename SubMatrixType = typename ConfiguratorType::MatrixType>
class MechanicalHessian :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::BlockMatrix<SubMatrixType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  aol::BlockMatrix<typename ConfiguratorType::MatrixType> _stiffnessMatrix;

public:
  MechanicalHessian( const typename ConfiguratorType::InitType &Grid, const HyperelasticEnergyDensityType &EnergyDensity, const aol::BitVector &DirichletNodes ) :
    _stiffnessMatrix( Grid ) {
    aol::MultiVector<RealType> zero( Grid );
    qc::HyperelasticHessian<ConfiguratorType,HyperelasticEnergyDensityType,typename ConfiguratorType::MatrixType>( Grid, EnergyDensity ).apply( zero, _stiffnessMatrix );
    // set rows and columns belonging to Dirichlet-nodes to the identity
    for ( int nodeIndex = 0; nodeIndex < Grid.getNumberOfNodes(); nodeIndex++ )
      if ( DirichletNodes[nodeIndex] ) {
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          for ( int j = 0; j < ConfiguratorType::Dim; j++ )
            _stiffnessMatrix.getReference( i, j ).setRowColToZero( nodeIndex );
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          _stiffnessMatrix.getReference( i, i ).set( nodeIndex, nodeIndex, 1. );
      }
  }

  void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::BlockMatrix<SubMatrixType> &Dest ) const {
    Dest += _stiffnessMatrix;
  }
};

#else

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class MechanicalEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const qc::HyperelasticEnergy<ConfiguratorType, HyperelasticEnergyDensityType> _internalEnergy;
  const aol::MultiVector<RealType> &_externalEnergy;

public:
  MechanicalEnergy( const typename ConfiguratorType::InitType &Grid, const HyperelasticEnergyDensityType &EnergyDensity, const aol::MultiVector<RealType> &ExternalEnergy ) :
    _internalEnergy( Grid, EnergyDensity ),
    _externalEnergy( ExternalEnergy ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    _internalEnergy.applyAdd( Arg, Dest );
    Dest -= _externalEnergy * Arg;
  }

  void apply( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    _internalEnergy.apply( Arg, Dest );
    Dest -= _externalEnergy * Arg;
  }
};

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class MechanicalDerivative :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const qc::HyperelasticGradient<ConfiguratorType, HyperelasticEnergyDensityType> _internalEnergyDeriv;
  const aol::MultiVector<RealType> &_externalEnergyDeriv;
  const aol::BitVector &_dirichletNodes;

public:
  MechanicalDerivative( const typename ConfiguratorType::InitType &Grid, const HyperelasticEnergyDensityType &EnergyDensity, const aol::MultiVector<RealType> &ExternalEnergy, const aol::BitVector &DirichletNodes ) :
    _internalEnergyDeriv( Grid, EnergyDensity ),
    _externalEnergyDeriv( ExternalEnergy ),
    _dirichletNodes( DirichletNodes ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> dest( Dest, aol::STRUCT_COPY );
    apply( Arg, dest );
    Dest += dest;
  }

  void apply( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    _internalEnergyDeriv.apply( Arg, Dest );
    Dest -= _externalEnergyDeriv;
    // set values belonging to Dirichlet-nodes to zero
    for ( int j = 0; j < Dest[0].size(); j++ )
      if ( _dirichletNodes[j] )
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          Dest[i][j] = 0.;
  }
};

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType, typename SubMatrixType = typename ConfiguratorType::MatrixType>
class MechanicalHessian :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::BlockMatrix<SubMatrixType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const qc::HyperelasticHessian<ConfiguratorType, HyperelasticEnergyDensityType> _internalEnergyHessian;
  const aol::BitVector &_dirichletNodes;

public:
  MechanicalHessian( const typename ConfiguratorType::InitType &Grid, const HyperelasticEnergyDensityType &EnergyDensity, const aol::BitVector &DirichletNodes ) :
    _grid( Grid ),
    _internalEnergyHessian( Grid, EnergyDensity ),
    _dirichletNodes( DirichletNodes ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::BlockMatrix<SubMatrixType> &Dest ) const {
    aol::BlockMatrix<SubMatrixType> dest( _grid );
    apply( Arg, dest );
    Dest += dest;
  }

  void apply( const aol::MultiVector<RealType> &Arg, aol::BlockMatrix<SubMatrixType> &Dest ) const {
    _internalEnergyHessian.apply( Arg, Dest );
    // set rows and columns belonging to Dirichlet-nodes to the identity
    for ( int nodeIndex = 0; nodeIndex < Arg[0].size(); nodeIndex++ )
      if ( _dirichletNodes[nodeIndex] ) {
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          for ( int j = 0; j < ConfiguratorType::Dim; j++ )
            Dest.getReference( i, j ).setRowColToZero( nodeIndex );
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          Dest.getReference( i, i ).set( nodeIndex, nodeIndex, 1. );
      }
  }
};

#endif

template<typename RealType, typename SubMatrixType, typename SolverType, qc::Dimension Dim = qc::QC_2D>
class NewtonMinimizerElasticity :
  public aol::NewtonMinimizationBase< RealType, aol::MultiVector<RealType>, aol::MultiVector<RealType>, aol::BlockMatrix<SubMatrixType> > {

  mutable aol::Op<aol::MultiVector<RealType> > *_pPrecond;

public:
  NewtonMinimizerElasticity( const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &F, const aol::Op<aol::MultiVector<RealType>,aol::MultiVector<RealType> > &DF, const aol::Op<aol::MultiVector<RealType>,aol::BlockMatrix<SubMatrixType> > &D2F, aol::BlockMatrix<SubMatrixType> *PMatD2F, const int MaxIterations = 50, const RealType StopEpsilon = 1.e-6, const bool WriteTimeSteps = false, const char *BaseSaveName = NULL ) :
    aol::NewtonMinimizationBase< RealType, aol::MultiVector<RealType>, aol::MultiVector<RealType>, aol::BlockMatrix<SubMatrixType> >( F, DF, D2F, PMatD2F, MaxIterations, StopEpsilon, WriteTimeSteps, BaseSaveName ),
    _pPrecond( NULL ) {}

  ~NewtonMinimizerElasticity() {
    delete _pPrecond;
  }

  virtual void prepareSolver() const {
#ifdef DIRECT_SOLVER
    if ( this->_pSolver == NULL )
      this->_pSolver = new aol::CholeskyBlockInverseOp<RealType,SubMatrixType>( Dim == qc::QC_2D ? 9 : 27 );
    static_cast<aol::CholeskyBlockInverseOp<RealType,SubMatrixType>*>( this->_pSolver )->setMatrix( *this->_pMatDF );
#else
    delete this->_pSolver;
    delete _pPrecond;
    _pPrecond = new aol::BlockDiagonalPreconditioner<RealType,Dim>( *this->_pMatDF );
    //_pPrecond = new aol::DiagonalPreconditioner< aol::MultiVector<RealType> >( *this->_pMatDF );
    //_pPrecond = new aol::SORBlockPreconditioner<RealType,SubMatrixType>( *this->_pMatDF );
    this->_pSolver = new SolverType ( *this->_pMatDF, *_pPrecond, this->_pInfo->getSolverInfo() );
    this->_pSolver->setStopping( aol::STOPPING_ABSOLUTE ); // the solver-accuracy is automatically set to fit to the Newton-minimizer accuracy
#endif
    //cerr<<endl;
  }
};

template <typename ConfiguratorType>
class ElasticProblemOp :
  public aol::Op< aol::Vector<typename ConfiguratorType::RealType>,aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef PhasefieldWeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > PhasefieldDensityType;
  typedef RestrictedDomainHyperelasticEnergyDensity<ConfiguratorType,PhasefieldDensityType> DensityType;
  typedef typename ConfiguratorType::MatrixType SubMatrixType; // alternative: aol::SparseMatrix<RealType>
#ifdef DIRECT_SOLVER
  typedef NewtonMinimizerElasticity< RealType, SubMatrixType, aol::CholeskyBlockInverseOp<RealType,SubMatrixType>, ConfiguratorType::Dim > NewtonMinimizerType;
#else
  typedef NewtonMinimizerElasticity< RealType, SubMatrixType, aol::PCGInverse<aol::MultiVector<RealType> >, ConfiguratorType::Dim > NewtonMinimizerType;
#endif

  const typename ConfiguratorType::InitType &_grid;
  const RealType _delta;
  const qc::BitArray<ConfiguratorType::Dim> &_domain; // true for all nodes which are inside the computational domain
  aol::BitVector _nonDOF; // true for all nodes which are no degrees of freedom (outside computational domain or Dirichlet boundary)
  const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> _defaultDensity;
  aol::MultiVector<RealType> _surfaceWorkOp;
  mutable aol::auto_container<2,aol::MultiVector<RealType> > _recentValues; // the most resent one at position 0
  const bool _quietMode;

public:
  ElasticProblemOp( const typename ConfiguratorType::InitType &Grid, const qc::BitArray<ConfiguratorType::Dim> &Domain, const aol::BitVector &DirichletNodes, const aol::MultiVector<RealType> &SurfaceLoad, const RealType Lambda, const RealType Mu, const RealType Delta, const bool QuietMode = true ) :
    _grid( Grid ),
    _delta( Delta ),
    _domain( Domain ),
    _nonDOF( Domain ),
    _defaultDensity( Mu, 0, Lambda ),
    _surfaceWorkOp( _grid ),
    _quietMode( QuietMode ) {
    // create linear operator for non-homogeneous Neumann bc on the domain at z=1
    BoundaryFaceMassOp<ConfiguratorType,ConfiguratorType::Dim == qc::QC_3D ? qc::BoundaryFaceElementBase<RealType>::Z_UPPER_BOUNDARY : qc::BoundaryFaceElementBase<RealType>::Y_UPPER_BOUNDARY > massOp( Grid, aol::ASSEMBLED );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      massOp.apply( SurfaceLoad[i], _surfaceWorkOp[i] );
    // initialize the recent value list
    aol::MultiVector<RealType> aux( Grid );
    for ( int i = 0; i < 2; i++ )
      _recentValues.set_copy( i, aux );
    // create BitVector which is true for Dirichlet nodes or nodes outside the computational domain
    _nonDOF.invert();
    _nonDOF |= DirichletNodes;
  }

  const aol::MultiVector<RealType>& getExternalPotentialOp() const {
    return _surfaceWorkOp;
  }

  void getDisplacement( aol::MultiVector<RealType> &Displacement ) const {
    Displacement = _recentValues[0];
  }

  void setDisplacement( const aol::MultiVector<RealType> &Displacement ) const {
    _recentValues[0] = Displacement;
  }

  RealType computeEnergy( const aol::Vector<RealType> &Phasefield, const aol::MultiVector<RealType> &Displacement ) const {
    aol::Scalar<RealType> result;
    PhasefieldDensityType phasefieldDensity( _defaultDensity, _grid, Phasefield, _delta );
    DensityType density( phasefieldDensity, _domain );
    MechanicalEnergy<ConfiguratorType, DensityType>( _grid, density, _surfaceWorkOp ).apply( Displacement, result );
    return result.v;
  }

  void apply( const aol::Vector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // define the hyperelastic energy density
    PhasefieldDensityType phasefieldDensity( _defaultDensity, _grid, Arg, _delta );
    DensityType density( phasefieldDensity, _domain );

    // define the total mechanical energy and its derivatives
    MechanicalEnergy<ConfiguratorType, DensityType> e( _grid, density, _surfaceWorkOp );
    MechanicalDerivative<ConfiguratorType, DensityType> de( _grid, density, _surfaceWorkOp, _nonDOF );
    MechanicalHessian<ConfiguratorType, DensityType> d2e( _grid, density, _nonDOF ); // the Hessian shall have unit rows and columns for all nodes which are not degrees of freedom

    // three alternative optimization routines:
    // 1. trust region method
    aol::TrustRegionMethod<RealType,aol::MultiVector<RealType>,aol::BlockMatrix<SubMatrixType>,aol::BlockDiagonalPreconditioner<RealType,ConfiguratorType::Dim>,SubMatrixType> minimizer(e, de, d2e, new aol::BlockMatrix<SubMatrixType>( _grid ), 5000, 5.e-12 ); // a higher accuracy than |gradient|<5e-7 can sometimes not be reached, since then the energy already varies within the machine epsilon
    minimizer.setQuietMode( _quietMode );
    minimizer.setIterativeMode( false );
    // 2. Newton descent method
    /*NewtonMinimizerType minimizer( e, de, d2e, new aol::BlockMatrix<SubMatrixType>( _grid ), 500, 1.e-7, false );
    minimizer.getNewtonInfo().setQuietMode( _quietMode );
    minimizer.setTimestepController( aol::NewtonInfo<RealType>::ARMIJO );
    minimizer.setSolverQuietMode( true );*/
    // 3. gradient descent method
    /*aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > minimizer( _grid, e, de, 500, 1., 1.e-7 );
    minimizer.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::USE_NONLINEAR_CG | aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::DO_NOT_SMOOTH_DESCENT_DIRECTION );
    minimizer.setQuietMode( _quietMode );*/

    // check which of the recent values is smaller before starting the iteration;
    // this may be important since for most linesearch techniques, the second last energy evaluation usually is the current one...
    // can be switched off by letting the loop go to -1
    int i = 0;
    aol::Scalar<RealType> energy;
    RealType min = aol::NumberTrait<RealType>::Inf;
    for ( int j = 0; j < 2; j++ ) {
      e.apply( _recentValues[j], energy );
      if ( energy.v < min ) {
        i = j;
        min = energy.v;
      }
    }

    minimizer.apply( _recentValues[i], Dest );
    _recentValues[1] = _recentValues[0];
    _recentValues[0] = Dest;
    //cerr<<endl;
  }

  void applyAdd( const aol::Vector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> aux;
    apply( Arg, aux );
    Dest += aux;
  }
};

template <typename ConfiguratorType>
class ShapeOptimizationEnergy :
  public aol::Op< aol::Vector<typename ConfiguratorType::RealType>,aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef PhasefieldWeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > PhasefieldDensityType;
  typedef RestrictedDomainHyperelasticEnergyDensity<ConfiguratorType,PhasefieldDensityType> DensityType;

  const typename ConfiguratorType::InitType &_grid;
  const RealType _lambda, _mu, _alpha, _nu, _eta, _delta, _epsilon;

  const DoubleWellPotential<RealType> _potential;
  const FENonlinIntegrationOp<ConfiguratorType,DoubleWellPotential<RealType> > _doubleWellPotential;
  const aol::StiffOp<ConfiguratorType> _stiffOp;
  const aol::MassOp<ConfiguratorType> _massOp;
  const ElasticProblemOp<ConfiguratorType> &_elasticEnergyMinimizer;
  const qc::FELinBoundaryMassOp<ConfiguratorType> _boundaryMassOp;

  const qc::BitArray<ConfiguratorType::Dim> &_domain;
  const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> _defaultDensity;

public:
  ShapeOptimizationEnergy( const typename ConfiguratorType::InitType &Grid, const qc::BitArray<ConfiguratorType::Dim> &Domain, const RealType Lambda, const RealType Mu, const RealType Alpha, const RealType Nu, const RealType Eta, const RealType Delta, const RealType Epsilon, const ElasticProblemOp<ConfiguratorType> &ElasticEnergyMinimizer ) :
    _grid( Grid ),
    _lambda( Lambda ),
    _mu( Mu ),
    _alpha( Alpha ),
    _nu( Nu ),
    _eta( Eta ),
    _delta( Delta ),
    _epsilon( Epsilon ),
    _potential(),
    _doubleWellPotential( Grid, _potential ),
    _stiffOp( Grid, aol::ASSEMBLED ),
    _massOp( Grid, aol::ASSEMBLED ),
    _elasticEnergyMinimizer( ElasticEnergyMinimizer ),
    _boundaryMassOp( Grid, aol::ASSEMBLED ),
    _domain( Domain ),
    _defaultDensity( Mu, 0, Lambda ) {}

  void computeEnergyComponents( const aol::Vector<RealType> &Arg, aol::Vec<3,RealType> &Energies ) const {
    aol::Scalar<RealType> aux;
    aol::Vector<RealType> auxVec( _grid ), auxVec2( Arg, aol::DEEP_COPY );

    // Allen-Cahn energy
    _doubleWellPotential.apply( Arg, aux );
    Energies[0] = aux.v * _nu / _epsilon;
    _stiffOp.apply( Arg, auxVec );
    Energies[0] += _nu * _epsilon / 2 * ( auxVec * Arg );

    // interface energy at the boundary (to avoid that edges form at the domain boundary)
    auxVec2.addToAll( 1. );
    _boundaryMassOp.apply( auxVec2, auxVec );
    Energies[0] += _nu / 3 * sqrt( 2 ) * ( auxVec * auxVec2 );

    // volume energy
    _massOp.apply( auxVec2, auxVec );
    Energies[1] = _eta / 4 * ( auxVec * auxVec2 );

    // deformation energy
    aol::MultiVector<RealType> displacement( _grid );
    _elasticEnergyMinimizer.apply( Arg, displacement );
    PhasefieldDensityType phasefieldDensity( _defaultDensity, _grid, Arg, _delta );
    DensityType density( phasefieldDensity, _domain );
    aol::Scalar<RealType> dest;
#ifdef LINEARIZED_ELASTICITY
    aol::MultiVector<RealType> zeroAux( displacement, aol::STRUCT_COPY );
    aol::BlockMatrix<typename ConfiguratorType::MatrixType> stiffnessMatrix( _grid );
    qc::HyperelasticHessian<ConfiguratorType,DensityType,typename ConfiguratorType::MatrixType>( _grid, density ).apply( zeroAux, stiffnessMatrix );
    stiffnessMatrix.apply( displacement, zeroAux );
    dest.v = .5 * ( zeroAux * displacement );
#endif
#ifdef USE_INTERNAL_ELASTIC_ENERGY
    qc::HyperelasticEnergy<ConfiguratorType, DensityType>( _grid, density ).apply( displacement, dest );
#endif
#ifdef USE_EXTERNAL_LOAD_POTENTIAL
    dest.v = _elasticEnergyMinimizer.getExternalPotentialOp() * displacement;
#endif
#ifdef USE_ELASTIC_DISSIPATION
    MechanicalEnergy<ConfiguratorType, DensityType>( _grid, density, _elasticEnergyMinimizer.getExternalPotentialOp() ).apply( displacement, dest );
    dest *= -1.;
#endif
    Energies[2] = dest.v * _alpha;
  }

  void applyAdd( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::Vec<3,RealType> energies;
    computeEnergyComponents( Arg, energies );
    for ( int i = 0; i < 3; i++ )
      Dest += energies[i];
  }
};

template <typename ConfiguratorType>
class ShapeOptimizationDerivative :
  public aol::Op< aol::Vector<typename ConfiguratorType::RealType>,aol::Vector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > WeightedDensityType;
  typedef RestrictedDomainHyperelasticEnergyDensity<ConfiguratorType,WeightedDensityType> DensityType;
  typedef PhasefieldWeightedHyperelasticEnergyDensity<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > PhasefieldDensityType;
  typedef RestrictedDomainHyperelasticEnergyDensity<ConfiguratorType,PhasefieldDensityType> DensityType2;

  const typename ConfiguratorType::InitType &_grid;
  const RealType _lambda, _mu, _alpha, _nu, _eta, _delta, _epsilon;

  const DoubleWellPotential<RealType> _potential;
  const FENonlinOp<ConfiguratorType,DoubleWellPotential<RealType> > _doubleWellPotentialDeriv;
  const aol::StiffOp<ConfiguratorType> _stiffOp;
  const aol::MassOp<ConfiguratorType> _massOp;
  const ElasticProblemOp<ConfiguratorType> &_elasticEnergyMinimizer;
  aol::BitVector _nonDOF;
  const qc::FELinBoundaryMassOp<ConfiguratorType> _boundaryMassOp;

  const qc::BitArray<ConfiguratorType::Dim> &_domain;
  const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> _defaultDensity;

  const int _maxSteps;

public:
  ShapeOptimizationDerivative( const typename ConfiguratorType::InitType &Grid, const qc::BitArray<ConfiguratorType::Dim> &Domain, const RealType Lambda, const RealType Mu, const RealType Alpha, const RealType Nu, const RealType Eta, const RealType Delta, const RealType Epsilon, const ElasticProblemOp<ConfiguratorType> &ElasticEnergyMinimizer, const aol::BitVector &DirichletNodes, const int MaxSteps ) :
    _grid( Grid ),
    _lambda( Lambda ),
    _mu( Mu ),
    _alpha( Alpha ),
    _nu( Nu ),
    _eta( Eta ),
    _delta( Delta ),
    _epsilon( Epsilon ),
    _potential(),
    _doubleWellPotentialDeriv( Grid, _potential ),
    _stiffOp( Grid, aol::ASSEMBLED ),
    _massOp( Grid, aol::ASSEMBLED ),
    _elasticEnergyMinimizer( ElasticEnergyMinimizer ),
    _nonDOF( Domain ),
    _boundaryMassOp( Grid, aol::ASSEMBLED ),
    _domain( Domain ),
    _defaultDensity( Mu, 0, Lambda ),
    _maxSteps( MaxSteps ) {
    _nonDOF.invert();
    _nonDOF |= DirichletNodes;
  }

  void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> auxVec( _grid ), auxVec2( Arg, aol::DEEP_COPY );

    // Allen-Cahn derivative
    _doubleWellPotentialDeriv.apply( Arg, Dest );
    Dest *= _nu / _epsilon;
    _stiffOp.apply( Arg, auxVec );
    Dest.addMultiple( auxVec, _nu * _epsilon );

    // interface derivative at the boundary (to avoid that edges form at the domain boundary)
    auxVec2.addToAll( 1. );
    _boundaryMassOp.apply( auxVec2, auxVec );
    Dest.addMultiple( auxVec, _nu * 2 / 3 * sqrt( 2 ) );

    // volume derivative
    _massOp.apply( auxVec2, auxVec );
    Dest.addMultiple( auxVec, _eta / 2 );

    // deformation energy weight derivative
    aol::MultiVector<RealType> displacement( _grid );
    _elasticEnergyMinimizer.apply( Arg, displacement );
    auxVec2 *= 2;
    WeightedDensityType weightedDensity( _defaultDensity, _grid, auxVec2 );
    DensityType density( weightedDensity, _domain );
#ifdef LINEARIZED_ELASTICITY
    ElasticLinearizedEnergyWeightDerivative<ConfiguratorType,DensityType>( _grid, displacement, displacement, density ).apply( Arg, auxVec );
    Dest.addMultiple( auxVec, .5 * _alpha );
#endif
#ifdef USE_INTERNAL_ELASTIC_ENERGY
    ElasticEnergyWeightDerivative<ConfiguratorType,DensityType>( _grid, displacement, density ).apply( Arg, auxVec );
    Dest.addMultiple( auxVec, _alpha );
#endif
#ifdef USE_EXTERNAL_LOAD_POTENTIAL
#endif
#ifdef USE_ELASTIC_DISSIPATION
    ElasticEnergyWeightDerivative<ConfiguratorType,DensityType>( _grid, displacement, density ).apply( Arg, auxVec );
    Dest.addMultiple( auxVec, -_alpha );
#endif

    // compute dual variable
    aol::MultiVector<RealType> p( _grid ), dJdPhi( _grid );
    PhasefieldDensityType phasefieldDensity( _defaultDensity, _grid, Arg, _delta );
    DensityType2 density2( phasefieldDensity, _domain );
    aol::BlockMatrix<typename ConfiguratorType::MatrixType> elasticHessian( _grid );
    MechanicalHessian<ConfiguratorType, DensityType2>( _grid, density2, _nonDOF ).apply( displacement, elasticHessian );
#ifdef LINEARIZED_ELASTICITY
    aol::MultiVector<RealType> zero( displacement, aol::STRUCT_COPY );
    aol::BlockMatrix<typename ConfiguratorType::MatrixType> stiffnessMatrix( _grid );
    qc::HyperelasticHessian<ConfiguratorType, DensityType2, typename ConfiguratorType::MatrixType>( _grid, density2 ).apply( zero, stiffnessMatrix );
    stiffnessMatrix.apply( displacement, dJdPhi );
#endif
#ifdef USE_INTERNAL_ELASTIC_ENERGY
    qc::HyperelasticGradient<ConfiguratorType, DensityType2>( _grid, density2 ).apply( displacement, dJdPhi );
#endif
#ifdef USE_EXTERNAL_LOAD_POTENTIAL
    dJdPhi = _elasticEnergyMinimizer.getExternalPotentialOp();
#endif
#ifdef USE_ELASTIC_DISSIPATION
    qc::HyperelasticGradient<ConfiguratorType, DensityType2>( _grid, density2 ).apply( displacement, dJdPhi );
    dJdPhi *= -1.;
    dJdPhi += _elasticEnergyMinimizer.getExternalPotentialOp();
#endif
    dJdPhi *= -_alpha;
    for ( int j = 0; j < dJdPhi[0].size(); j++ ) // implement the homogeneous Dirichlet conditions for the displacement (also outside computational domain)
      if ( _nonDOF[j] )
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          dJdPhi[i][j] = 0.;
    // 1. alternative: direct solver
    aol::SolverInfo<RealType> info( 1e-16, _maxSteps, aol::STOPPING_ABSOLUTE );
    aol::CholeskyBlockInverseOp<RealType,typename ConfiguratorType::MatrixType> solver( ConfiguratorType::Dim == qc::QC_2D ? 9 : 27 );
    solver.setMatrix( elasticHessian );
    /*// 2. alternative: iterative solver
    const aol::BlockDiagonalPreconditioner<RealType,ConfiguratorType::Dim> preCond( elasticHessian );
    aol::PCGInverse<aol::MultiVector<RealType> > solver( elasticHessian, preCond, 1e-16, _maxSteps, aol::STOPPING_ABSOLUTE ); // _maxSteps~3*matrix dimension
    solver.setQuietMode( true );*/
    solver.apply( dJdPhi, p );

    // compute dual variable part of the derivative
    ElasticLinearizedEnergyWeightDerivative<ConfiguratorType,DensityType>( _grid, displacement, p, density ).applyAdd( Arg, Dest );

    // implement Dirichlet values of the double well at z=1
    aol::Vector<RealType> zeros( ConfiguratorType::Dim == qc::QC_3D ? _grid.getNumX() * _grid.getNumY() : _grid.getNumX() );
    Dest.setBlock( Dest.size() - zeros.size(), zeros );
  }

  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> dest( Dest, aol::STRUCT_COPY );
    apply( Arg, dest );
    Dest += dest;
  }
};

template <typename ConfiguratorType>
class ShapeOptimizationOp {

private:
  typedef typename ConfiguratorType::RealType RealType;

  // numerical parameters
  const int _coarsestLevel, _finestLevel, _maxSteps, _maxIterations, _maxRepetitions;
  // model parameters
  const RealType _lambda, _mu, _delta, _epsFactor;
  aol::Vector<RealType> _alpha, _nu, _eta;
  RealType _epsilon;
  // the surface load and the double well phase field representing the material distribution
  qc::MultiDimMultilevelArray<RealType,qc::Array<RealType>,qc::ProlongOp<RealType>,qc::RestrictOp<RealType,qc::THROW_AWAY_RESTRICT> > _surfaceLoad;
  qc::MultiDimMultilevelArray<RealType> _displacement;
  qc::MultilevelArray<RealType> _dirichletNodes;
  qc::MultilevelArray<RealType> _phasefield;
  // the directory for the results
  char _destDirectory[1024];

public:
  ShapeOptimizationOp( aol::ParameterParser &Parser ) :
    // load numerical parameters
    _coarsestLevel( Parser.getInt( "coarsestLevel" ) ),
    _finestLevel( Parser.getInt( "GridDepth" ) ),
    _maxSteps( Parser.getInt( "maxSteps" ) ),
    _maxIterations( Parser.getInt( "maxIterations" ) ),
    _maxRepetitions( Parser.getInt( "maxRepetitions" ) ),
    // load model parameters
    _lambda( Parser.getDouble( "lambda" ) ),
    _mu( Parser.getDouble( "mu" ) ),
    _delta( Parser.getDouble( "delta" ) ),
    _epsFactor( Parser.getDouble( "epsFactor" ) ),
    _alpha( _finestLevel - _coarsestLevel + 1 ),
    _nu( _finestLevel - _coarsestLevel + 1 ),
    _eta( _finestLevel - _coarsestLevel + 1 ),
    // initialize surface load and phase field
    _surfaceLoad( _finestLevel, ConfiguratorType::Dim, ConfiguratorType::Dim ),
    _displacement( _finestLevel, ConfiguratorType::Dim, ConfiguratorType::Dim ),
    _dirichletNodes( _finestLevel, ConfiguratorType::Dim ),
    _phasefield( _finestLevel, ConfiguratorType::Dim ) {
    // read in the directory, where results are to be saved
    Parser.getString( "destDirectory", _destDirectory );

    // read in numerical parameters
    if ( Parser.getNumDim( "alpha" ) > 0 )
      Parser.getRealVec( "alpha", _alpha );
    else
      _alpha.setAll( Parser.getDouble( "alpha" ) );
    if ( Parser.getNumDim( "nu" ) > 0 )
      Parser.getRealVec( "nu", _nu );
    else
      _nu.setAll( Parser.getDouble( "nu" ) );
    if ( Parser.getNumDim( "eta" ) > 0 )
      Parser.getRealVec( "eta", _eta );
    else
      for ( int i = 0; i < _eta.size(); i++ )
        _eta[i] = ( 1 << i ) * Parser.getDouble( "eta" );

    // define the surface load
    aol::Vector<RealType> loads;
    if ( Parser.hasVariable( "loads" ) )
      Parser.getRealVec( "loads", loads );
    int n = 1 << ( _finestLevel - 4 );
    qc::Array<RealType> xLoad( _surfaceLoad.getArray( 0 ), aol::FLAT_COPY );
    qc::Array<RealType> zLoad( _surfaceLoad.getArray( ConfiguratorType::Dim - 1 ), aol::FLAT_COPY );
    switch ( Parser.getInt( "loadCase" ) ) {
    case 1: // compression
      _surfaceLoad.getArray( ConfiguratorType::Dim - 1 ).setAll( -1.e-1 ); // load is only integrated on the boundary z=1, anyway
      break;
    case 2: // dilation
      _surfaceLoad.getArray( ConfiguratorType::Dim - 1 ).setAll( 1. ); // load is only integrated on the boundary z=1, anyway
      break;
    case 3: // shearing
      _surfaceLoad.getArray( 0 ).setAll( 1. ); // load is only integrated on the boundary z=1, anyway
      break;
    case 4: // "bridge"
      _surfaceLoad.getArray( ConfiguratorType::Dim - 1 ).setAll( loads[0] ); // load is only integrated on the boundary z=1, anyway
      break;
    case 5: // general load (e.g. compression + shearing)
      _surfaceLoad.getArray( ConfiguratorType::Dim - 1 ).setAll( loads[0] ); // load is only integrated on the boundary z=1, anyway
      _surfaceLoad.getArray( 0 ).setAll( loads[1] ); // load is only integrated on the boundary z=1, anyway
      break;
    case 6: // "bridge" with local loads
    case 7: // "chain bridge" with local loads
      for ( int y = ( ( ConfiguratorType::Dim == qc::QC_2D ) ? zLoad.getNumY() - 1 : 0 ); y < zLoad.getNumY(); y++ )
        for ( int i = -n; i < n; i++ ) {
          zLoad.set( 8 * n + i, y, zLoad.getNumZ() - 1, ( loads[0] * ( n - aol::Abs( i ) ) ) / n ); // middle zLoad
          zLoad.set( 4 * n + i, y, zLoad.getNumZ() - 1, ( loads[0] * ( n - aol::Abs( i ) ) ) / n ); // left zLoad
          zLoad.set( 12 * n + i, y, zLoad.getNumZ() - 1, ( loads[0] * ( n - aol::Abs( i ) ) ) / n ); // right zLoad
        }
      break;
    case 8: // cantilever
    case 9: // long cantilever
      for ( int y = ( ( ConfiguratorType::Dim == qc::QC_2D ) ? xLoad.getNumY() - 1 : 0 ); y < xLoad.getNumY(); y++ )
        for ( int i = -n; i < n; i++ )
          xLoad.set( 8 * n + i, y, xLoad.getNumZ() - 1, ( loads[0] * ( n - aol::Abs( i ) ) ) / n ); // middle load
      break;
    case 10: // "3D bridge"
      _surfaceLoad.getArray( ConfiguratorType::Dim - 1 ).setAll( loads[0] ); // load is only integrated on the boundary z=1, anyway
      break;
    default:
      throw aol::Exception ( "Missing default case", __FILE__, __LINE__ );
    };

    // define the Dirichlet data on the coarsest grid, then prolongate
    // so that only new nodes with their closest neighbours being Dirichlet nodes may become Dirichlet nodes
    if ( _finestLevel > _coarsestLevel )
      _dirichletNodes.levRestrict( _coarsestLevel, _finestLevel );
    _dirichletNodes.setCurLevel( _coarsestLevel );
    _dirichletNodes.current().setAll( 1. );
    if ( ( Parser.getInt( "loadCase" ) == 4 ) || ( Parser.getInt( "loadCase" ) == 6 ) ) // bridge: Dirichlet at z=1 and only at the boundary
      for ( int y = ( ( ConfiguratorType::Dim == qc::QC_2D ) ? _dirichletNodes.current().getNumY() - 1 : 0 ); y < _dirichletNodes.current().getNumY(); y++ )
        for ( int x1 = 0, x2 = _dirichletNodes.current().getNumX() - 1; x1 < 1 + ( 1/*#nodes*/ - 1 ) * ( 1 << ( _coarsestLevel - 4 ) ); x1++, x2-- ){
          _dirichletNodes.current().set( x1, y, _dirichletNodes.current().getNumZ() - 1, 0. );
          _dirichletNodes.current().set( x2, y, _dirichletNodes.current().getNumZ() - 1, 0. );
        }
    else if ( Parser.getInt( "loadCase" ) == 7 ) // chain bridge: Dirichlet at x=0 and x=1
      for ( int z = 0; z < _dirichletNodes.current().getNumZ(); z++ )
        for ( int y = 0; y < _dirichletNodes.current().getNumY(); y++ ) {
          _dirichletNodes.current().set( 0, y, z, 0. );
          _dirichletNodes.current().set( _dirichletNodes.current().getNumX() - 1, y, z, 0. );
        }
    else if ( Parser.getInt( "loadCase" ) == 9 ) // long cantilever: Dirichlet at z=0 and .25<x<.75
      for ( int y = 0; y < ( ( ConfiguratorType::Dim == qc::QC_2D ) ? 1 : _dirichletNodes.current().getNumY() ); y++ )
        for ( int x = _dirichletNodes.current().getNumX() / 4; x <= ( _dirichletNodes.current().getNumX() / 4 ) * 3; x++ )
          _dirichletNodes.current().set( x, y, _dirichletNodes.current().getNumZ() - 1, 0. );
    else if ( Parser.getInt( "loadCase" ) == 10 ) // 3D bridge
      for ( int x = 0; x < _dirichletNodes.current().getNumX(); x++ )
        for ( int y1 = 0, y2 = _dirichletNodes.current().getNumY() - 1; y1 < 1 + ( 1/*#nodes*/ - 1 ) * ( 1 << ( _coarsestLevel - 4 ) ); y1++, y2-- ){
          _dirichletNodes.current().set( x, y1, _dirichletNodes.current().getNumZ() - 1, 0. );
          _dirichletNodes.current().set( x, y2, _dirichletNodes.current().getNumZ() - 1, 0. );
        }
    else // Dirichlet at z=0
      for ( int y = 0; y < ( ( ConfiguratorType::Dim == qc::QC_2D ) ? 1 : _dirichletNodes.current().getNumY() ); y++ )
        for ( int x = 0; x < _dirichletNodes.current().getNumX(); x++ )
          _dirichletNodes.current().set( x, y, 0, 0. );
    _dirichletNodes.levProlongate( _coarsestLevel, _finestLevel );
    for ( int level = _coarsestLevel; level <= _finestLevel; level++ ) {
      _dirichletNodes[level] *= - ( 1 << ( level - _coarsestLevel ) );
      _dirichletNodes[level].clamp( -1., 0. );
      _dirichletNodes[level].addToAll( 1. );
    }

    // initialize the phase field and deformation
    if ( Parser.hasVariable( "InitPhasefield" ) ) {
      char filename[1024];
      Parser.getString( "InitPhasefield", filename );
      qc::ScalarArray<RealType,ConfiguratorType::Dim>( _phasefield[_finestLevel], aol::FLAT_COPY ).load( filename );
      if ( ConfiguratorType::Dim == qc::QC_2D ) {
        _phasefield[_finestLevel] *= 2. / 255.;
        _phasefield[_finestLevel].addToAll( -1. );
      }
    }
    if ( Parser.hasVariable( "InitDeformationTemplate" ) ) {
      char filenameTemplate[1024], filename[1024];
      Parser.getString( "InitDeformationTemplate", filenameTemplate );
      for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
        sprintf( filename, filenameTemplate, i );
        qc::ScalarArray<RealType,ConfiguratorType::Dim>( _displacement.getArray( i ), aol::FLAT_COPY ).load( filename );
      }
    }

    // restriction for multiscale approach
    if ( _finestLevel > _coarsestLevel ) {
      _surfaceLoad.levRestrict( _coarsestLevel, _finestLevel );
      _displacement.levRestrict( _coarsestLevel, _finestLevel );
      _phasefield.levRestrict( _coarsestLevel, _finestLevel );
    }

    // initialize the phase field with noise
    _phasefield.setCurLevel( _coarsestLevel );
    if ( !Parser.hasVariable( "InitPhasefield" ) )
      aol::NoiseOperator<RealType> ( aol::NoiseOperator<RealType>::EQUALLY_DISTRIBUTED, -.1, .1 ).apply( _phasefield.current(), _phasefield.current() );
  }

  void execute() {
    for ( int level = _coarsestLevel; level <= _finestLevel; level++ ) {
      cerr << aol::color::green << "Relaxing on level " << level << " of " << _finestLevel << aol::color::reset << endl;
      typename ConfiguratorType::InitType grid( level, ConfiguratorType::Dim );
      _epsilon = _epsFactor * grid.H();
      _surfaceLoad.setCurLevel( level );
      _displacement.setCurLevel( level );
      _phasefield.setCurLevel( level );
      aol::BitVector dirichletNodes( grid.getNumberOfNodes() );
      dirichletNodes.setNonzeroPatternFrom<RealType>( _dirichletNodes[level] );

      for ( int i = 0; i < _maxRepetitions; i++ ) {
        // if needed, reduce computational domain to save computation time
        qc::BitArray<ConfiguratorType::Dim> domain( grid );
        domain.setAll();
        // only compute where the phase field is sufficiently far away from -1
        for ( int x = 0; x < grid.getNumX(); x++ )
          for ( int y = 0; y < grid.getNumY(); y++ )
            for ( int z = 0; z < grid.getNumZ(); z++ ) {
              qc::CoordType coord( x, y, z );
              if ( _phasefield.current().get( coord ) < -.95 )
                domain.set( coord, false );
            }
        for ( int j = 0; j < 3; j++ )
          for ( int d = 0; d < ConfiguratorType::Dim; d++ )
            domain.dilateByOne( d );
        /*// only compute inside a middle stripe
        for ( int x = 0, xE = grid.getNumX() - 1; x < grid.getNumX() / 4; x++, xE-- )
          for ( int y = 0; y < grid.getNumY(); y++ )
            for ( int z = 0; z < grid.getNumZ(); z++ ) {
              qc::CoordType coord( x, y, z );
              domain.set( coord, false );
              coord[0] = xE;
              domain.set( coord, false );
            }*/
        /*// only compute inside the upper half
        if ( ConfiguratorType::Dim == qc::QC_2D )
          for ( int x = 0; x < grid.getNumX(); x++ )
            for ( int y = 0; y < grid.getNumY() / 2; y++ ) {
              qc::CoordType coord( x, y, 0 );
              domain.set( coord, false );
            }
        else
          for ( int x = 0; x < grid.getNumX(); x++ )
            for ( int y = 0; y < grid.getNumY(); y++ )
              for ( int z = 0; z < grid.getNumZ() / 2; z++ ) {
              qc::CoordType coord( x, y, z );
              domain.set( coord, false );
            }*/
        /*// 3D: leave a tube free
        for ( int x = grid.getNumX() * 7 / 16; x < grid.getNumX() * 10 / 16; x++ )
          for ( int y = 0; y < grid.getNumY(); y++ )
            for ( int z = grid.getNumZ() * 12 / 16; z < grid.getNumZ() * 14 / 16; z++ ) {
              qc::CoordType coord( x, y, z );
              domain.set( coord, false );
            }*/

        // set phasefield to 1 at the non-homogeneous Neumann boundary
        aol::Vector<RealType> ones( ConfiguratorType::Dim == qc::QC_2D ? grid.getNumX() : grid.getNumX() * grid.getNumY() );
        ones.setAll( 1. );
        _phasefield.current().setBlock( _phasefield.current().size() - ones.size(), ones );

        // define energy and derivative
        aol::MultiVector<RealType> surfaceLoad, displacement;
        _surfaceLoad.appendReferencesTo( surfaceLoad );
        _displacement.appendReferencesTo( displacement );
        ElasticProblemOp<ConfiguratorType> elasticEnergyMinimizer( grid, domain, dirichletNodes, surfaceLoad, _lambda, _mu, _delta, true );
        elasticEnergyMinimizer.setDisplacement( displacement );
        int k = level - _coarsestLevel;
        ShapeOptimizationEnergy<ConfiguratorType> e( grid, domain, _lambda, _mu, _alpha[k], _nu[k], _eta[k], _delta, _epsilon, elasticEnergyMinimizer );
        ShapeOptimizationDerivative<ConfiguratorType> de( grid, domain, _lambda, _mu, _alpha[k], _nu[k], _eta[k], _delta, _epsilon, elasticEnergyMinimizer, dirichletNodes, _maxIterations );

        // define minimization method
        aol::QuasiNewtonIteration<RealType,aol::Vector<RealType>,aol::Vector<RealType> > minimizer( e, de, _maxSteps, 1.e-6, 5, false );
        minimizer.setTauMin( .1 ); // experience shows: either tau=1 or tau=2 or tau=0
        minimizer.setBeta( .9 );
        minimizer.setTimestepController( aol::NewtonInfo<RealType>::WOLFE ); // otherwise the Quasi-Newton method yields no descent directions
        //aol::GradientDescent<ConfiguratorType,aol::Vector<RealType> > minimizer( grid, e, de, _maxSteps, 1., 1.e-6 );
        //minimizer.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::Vector<RealType> >::USE_NONLINEAR_CG | aol::GradientDescent<ConfiguratorType,aol::Vector<RealType> >::DO_NOT_SMOOTH_DESCENT_DIRECTION );
        aol::Vector<RealType> initialGuess( _phasefield.current(), aol::DEEP_COPY );

        // perform minimization
        minimizer.apply( initialGuess, _phasefield.current() );

        // if the computational domain is only a subset of the unit square, smoothly extend the displacement onto the unit square
        elasticEnergyMinimizer.getDisplacement( displacement );
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          aol::Vector<RealType> dispComp( displacement[j], aol::DEEP_COPY );
          qc::SmoothlyExtendImage<ConfiguratorType>( dispComp, grid, displacement[j], domain );
        }

        // save results
        if ( ConfiguratorType::Dim == qc::QC_2D ) {
          saveAsImage<ConfiguratorType,PGM>( _phasefield.current(), grid, _destDirectory, "phasefield", level * 10 + i, -1., 1., true );
          saveDisplacement<ConfiguratorType>( displacement, grid, _destDirectory, "displacement", level * 10 + i, true );
        } else {
          save<ConfiguratorType>( _phasefield.current(), grid, _destDirectory, "phasefield", level * 10 + i, true );
          aol::Vector<RealType> deformedPhasefield( _phasefield.current(), aol::STRUCT_COPY );
          deformedPhasefield.setAll( -1. );
          qc::DeformImage<ConfiguratorType>( _phasefield.current(), grid, deformedPhasefield, displacement, false );
          save<ConfiguratorType>( deformedPhasefield, grid, _destDirectory, "defPhasefield", level * 10 + i, true );
          saveDisplacement<ConfiguratorType>( displacement, grid, _destDirectory, "displacement", level * 10 + i, true );
        }
        aol::Vec<3,RealType> energyComps;
        e.computeEnergyComponents( _phasefield.current(), energyComps );
        cerr<<energyComps[0]<<" "<<energyComps[1]<<" "<<energyComps[2]<<" "<<endl;
      }

      _displacement.levProlongate();
      _phasefield.levProlongate();
    }
  }
};

template <typename ConfiguratorType>
class BucklingTest {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef typename ConfiguratorType::InitType InitType;

  const InitType _grid;
  const RealType _delta, _lambda, _mu;
  const RealType _minThickness, _maxThickness;
  const RealType _load;
  ArrayType _structure;
  qc::MultiArray<RealType,ConfiguratorType::Dim,ConfiguratorType::Dim> _buckledDisplacement, _symmetricDisplacement;
  qc::MultiArray<RealType,ConfiguratorType::Dim,ConfiguratorType::Dim> _surfaceLoad;
  qc::BitArray<ConfiguratorType::Dim> _dirichletNodes;

public:
  BucklingTest( const int Level = 6, const RealType MinThickness = .05, const RealType MaxThickness = .6, const RealType Delta = .01, const RealType Lambda = 1., const RealType Mu = 1., const RealType Load = .05 ) :
    _grid( Level, ConfiguratorType::Dim ),
    _delta( Delta ),
    _lambda( Lambda ),
    _mu( Mu ),
    _minThickness( MinThickness ),
    _maxThickness( MaxThickness ),
    _load( Load ),
    _structure( _grid ),
    _buckledDisplacement( _grid ),
    _symmetricDisplacement( _grid ),
    _surfaceLoad( _grid ),
    _dirichletNodes( qc::GridSize<ConfiguratorType::Dim>( _grid ) ) {
    // generate a Dirichlet mask at z = 0
    for ( int y = 0; y < ( ( ConfiguratorType::Dim == qc::QC_2D ) ? 1 : _grid.getNumY() ); y++ )
      for ( int x = 0; x < _grid.getNumX(); x++ ) {
        qc::CoordType node( x, y, 0 );
        _dirichletNodes.set( node, true );
      }
  }

  void execute() {
    const int minPixelWidth = ceil( _minThickness * _grid.getNumX() );
    const int maxPixelWidth = floor( _maxThickness * _grid.getNumX() );
    aol::MultiVector<RealType> energies( 5, maxPixelWidth - minPixelWidth + 1 );

    // try to find a hysteresis by first finding (initially symmetric) displacements with decreasing bar thickness,
    // then the buckling displacements with increasing bar thickness
    // if only the symmetric solution is to be found, exchange trust region method in ElasticProblemOp by Newton-method (iterative solver)

    // find symmetric displacements
    for ( int pixelWidth = maxPixelWidth, index = 0; pixelWidth >= minPixelWidth; pixelWidth--, index++ ) {
      // generate a vertical bar with prescribed thickness
      _structure.setAll( -1. );
      for ( int z = 0; z < _grid.getNumZ(); z++ )
        for ( int y = 0; y < _grid.getNumY(); y++ )
          for ( int x = _grid.getNumX() / 2 - ( pixelWidth / 2 ); x <= _grid.getNumX() / 2 - ( pixelWidth / 2 ) + pixelWidth; x++ ) {
            qc::CoordType node( x, y, z );
            _structure.set( node, 1. );
          }
      // generate a compression load at y = 1
      RealType stress = _load / ( pixelWidth + 1 ) / _grid.H();
      _surfaceLoad[ConfiguratorType::Dim - 1] = _structure;
      _surfaceLoad[ConfiguratorType::Dim - 1] *= - stress;
      _surfaceLoad[ConfiguratorType::Dim - 1].clamp( - stress, 0. ); // load is only integrated on the boundary z=1, anyway
      // find symmetric displacement by a homotopy method
      ElasticProblemOp<ConfiguratorType> elasticEnergyMinimizer( _grid, _surfaceLoad, _dirichletNodes, _lambda, _mu, _delta, false );
      elasticEnergyMinimizer.setDisplacement( _symmetricDisplacement );
      elasticEnergyMinimizer.apply( _structure, _symmetricDisplacement );
      // compute the energy contributions and save them
      energies[0][index] = ( pixelWidth + 1 ) * _grid.H(); // approximate bar thickness
      aol::Vec<3,RealType> lenVolElast;
      ShapeOptimizationEnergy<ConfiguratorType>( _grid, _lambda, _mu, 1, 1, 1, _delta, _grid.H(), elasticEnergyMinimizer ).computeEnergyComponents( _structure,  lenVolElast );
      for ( int i = 0; i < 3; i++ )
        energies[i + 1][index] = lenVolElast[i];
      elasticEnergyMinimizer.getDisplacement( _symmetricDisplacement );
      energies[4][index] = elasticEnergyMinimizer.computeEnergy( _structure, _symmetricDisplacement );
      saveAsImage<ConfiguratorType,PGM>( _structure, _grid, "./", "structure", pixelWidth, -1, 1, true );
      saveDisplacement<ConfiguratorType>( _symmetricDisplacement, _grid, "./", "symmetricDispl", pixelWidth, true );
    }

    // save results
    aol::Plotter<RealType> plotter;
    plotter.set_outfile_base_name( "symmetricBranch" );
    aol::PlotDataFileHandler<RealType> plotHandler;
    for ( int i = 1; i < 5; i++ )
      plotHandler.generateFunctionPlot( energies[0], energies[i] );
    plotter.addPlotCommandsFromHandler( plotHandler );
    plotter.setLabels( "bar thickness", "energies" );
    plotter.genPlot( aol::GNUPLOT_PNG );
    cerr<<energies<<endl;

    // generate a coarse approximation to a bending displacement
    qc::DataGenerator<ConfiguratorType>( _grid ).generateLineLevelset( _buckledDisplacement[0], 1, 0 );
    for ( int i = 0; i < _buckledDisplacement[0].size(); i++ )
      _buckledDisplacement[0][i] = cos( aol::NumberTrait<RealType>::pi * _buckledDisplacement[0][i] ) - 1.; // buckling mode 1
      //_buckledDisplacement[0][i] = 4 * aol::Sqr( _buckledDisplacement[0][i] ) * ( _buckledDisplacement[0][i] - 1. ); // buckling mode 2
    // alternatively, start from solution
    _buckledDisplacement = _symmetricDisplacement;

    // find buckling displacements
    for ( int pixelWidth = minPixelWidth, index = 0; pixelWidth <= maxPixelWidth; pixelWidth++, index++ ) {
      // generate a vertical bar with prescribed thickness
      _structure.setAll( -1. );
      for ( int z = 0; z < _grid.getNumZ(); z++ )
        for ( int y = 0; y < _grid.getNumY(); y++ )
          for ( int x = _grid.getNumX() / 2 - ( pixelWidth / 2 ); x <= _grid.getNumX() / 2 - ( pixelWidth / 2 ) + pixelWidth; x++ ) {
            qc::CoordType node( x, y, z );
            _structure.set( node, 1. );
          }
      // generate a compression load at y = 1
      RealType stress = _load / ( pixelWidth + 1 ) / _grid.H();
      _surfaceLoad[ConfiguratorType::Dim - 1] = _structure;
      _surfaceLoad[ConfiguratorType::Dim - 1] *= - stress;
      _surfaceLoad[ConfiguratorType::Dim - 1].clamp( - stress, 0. ); // load is only integrated on the boundary z=1, anyway
      // find buckling displacement by a homotopy method
      ElasticProblemOp<ConfiguratorType> elasticEnergyMinimizer( _grid, _surfaceLoad, _dirichletNodes, _lambda, _mu, _delta, false );
      elasticEnergyMinimizer.setDisplacement( _buckledDisplacement );
      elasticEnergyMinimizer.apply( _structure, _buckledDisplacement );
      // compute the energy contributions and save them
      energies[0][index] = ( pixelWidth + 1 ) * _grid.H(); // approximate bar thickness
      aol::Vec<3,RealType> lenVolElast;
      ShapeOptimizationEnergy<ConfiguratorType>( _grid, _lambda, _mu, 1, 1, 1, _delta, _grid.H(), elasticEnergyMinimizer ).computeEnergyComponents( _structure,  lenVolElast );
      for ( int i = 0; i < 3; i++ )
        energies[i + 1][index] = lenVolElast[i];
      elasticEnergyMinimizer.getDisplacement( _buckledDisplacement );
      energies[4][index] = elasticEnergyMinimizer.computeEnergy( _structure, _buckledDisplacement );
      saveDisplacement<ConfiguratorType>( _buckledDisplacement, _grid, "./", "bucklingDispl", pixelWidth, true );
    }

    // save results
    aol::Plotter<RealType> plotter2;
    plotter2.set_outfile_base_name( "buckledBranch" );
    aol::PlotDataFileHandler<RealType> plotHandler2;
    for ( int i = 1; i < 5; i++ )
      plotHandler2.generateFunctionPlot( energies[0], energies[i] );
    plotter2.addPlotCommandsFromHandler( plotHandler2 );
    plotter2.setLabels( "bar thickness", "energies" );
    plotter2.genPlot( aol::GNUPLOT_PNG );
    cerr<<energies<<endl;
  }
};

class TestEnergy :
  public aol::Op<aol::Vector<double>, aol::Scalar<double> > {
public:
  virtual void applyAdd( const aol::Vector<double> &Arg, aol::Scalar<double> &Dest ) const {
    Dest.v += -10 * aol::Sqr( Arg[0] ) + 10 * aol::Sqr( Arg[1] ) + 4 * sin( Arg[0] * Arg[1] ) - 2 * Arg[0] + aol::Sqr( aol::Sqr( Arg[0] ) );
  }
};

class TestDerivative :
  public aol::Op<aol::Vector<double>, aol::Vector<double> > {
public:
  virtual void applyAdd( const aol::Vector<double> &Arg, aol::Vector<double> &Dest ) const {
    Dest[0] += -20 * Arg[0] + 4 * cos( Arg[0] * Arg[1] ) * Arg[1] - 2 + 4 * aol::Cub( Arg[0] );
    Dest[1] += 20 * Arg[1] + 4 * cos( Arg[0] * Arg[1] ) * Arg[0];
  }
};

class TestSecondDerivative :
  public aol::Op<aol::Vector<double>, aol::FullMatrix<double> > {
public:
  virtual void applyAdd( const aol::Vector<double> &Arg, aol::FullMatrix<double> &Dest ) const {
    Dest.add( 0, 0, -20 - 4 * sin( Arg[0] * Arg[1] ) * aol::Sqr( Arg[1] ) + 12 * aol::Sqr( Arg[0] ) );
    Dest.add( 0, 1, 4 * cos( Arg[0] * Arg[1] ) - 4 * sin( Arg[0] * Arg[1] ) * Arg[0] * Arg[1] );
    Dest.add( 1, 0, 4 * cos( Arg[0] * Arg[1] ) - 4 * sin( Arg[0] * Arg[1] ) * Arg[0] * Arg[1] );
    Dest.add( 1, 1, 20 - 4 * sin( Arg[0] * Arg[1] ) * aol::Sqr( Arg[0] ) );
  }
};



// define the settings/configuration of the program, i.e. computation accuracy, dimension, gridtype,
// finite element type, used quadrature rules
typedef double RealType;
const qc::Dimension Dim = qc::QC_2D;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, Dim, aol::GaussQuadrature<RealType, Dim, 3> > ConfiguratorType;

int main( int argc, char *argv[] ) {

/*  // testing the trust region method
  TestEnergy e;
  TestDerivative de;
  TestSecondDerivative d2e;
  aol::Vector<double> start( 2 ), end( 2 );
  start[0] = .7067;
  start[1] = -3.2672;
  aol::TrustRegionMethod<RealType,aol::Vector<double>,aol::FullMatrix<double> >(e, de, d2e, new aol::FullMatrix<double>( 2, 2 ) ).apply( start, end );
  abort();
*/

/*  // testing the elastic optimization; parameters: level, minThickness, maxThickness, delta, lambda, mu, load
  BucklingTest<ConfiguratorType>( atoi( argv[1] ), atof( argv[2] ), atof( argv[3] ), atof( argv[4] ), atof( argv[5] ), atof( argv[6] ), atof( argv[7] ) ).execute(); abort();
*/

  try {
    char parameterFileName[1024];

    // read in file names of parameter files
    if ( argc > 50 ) {             // too many input arguments specified; explain correct syntax
      cerr << "Too many input files. Use n<=49: averagePhaseField <parameter filename 1> ... <parameter filename n>" << endl;
      return EXIT_FAILURE;
    }
    else {
      for ( int i=1; i<argc; i++ ) {
        sprintf( parameterFileName, "%s",  argv[i] );

        // read in parameters from specified parameter file
        cerr << endl << "Reading parameters from '" << parameterFileName << "'." << endl;
        aol::ParameterParser parameterParser( parameterFileName );

        // execute the algorithm
        ShapeOptimizationOp<ConfiguratorType>( parameterParser ).execute();
      }
    }
  }
  catch ( aol::Exception &el ) {
    // print out error message
    el.dump();
  }

  // pause before closing the program (to let the user first read the output)
  aol::callSystemPauseIfNecessaryOnPlatform();
  return EXIT_SUCCESS;
}
