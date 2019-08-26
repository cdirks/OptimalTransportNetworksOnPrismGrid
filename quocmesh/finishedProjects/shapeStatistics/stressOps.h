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

#ifndef __STRESSOPS_H
#define __STRESSOPS_H

#include <aol.h>
#include <deformations.h>
#include <hyperelastic.h>
#include <FELevelsetOpInterface.h>

/**************************************************************
 * Classes for stress evaluation and transformation.
 **************************************************************/

template <typename ConfiguratorType>
class ResolveInterpenetrationEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  const RealType _threshold;
  const RealType _factor;
  const std::vector<qc::CoordType> &_problematicNodes;
  std::vector<qc::CoordType> _adjacentElements;
  const typename ConfiguratorType::InitType &_grid;

  inline RealType evaluateNonlinearity( const RealType Det ) const {
    return ( Det > _threshold ) ? 0 : _factor * aol::Sqr( Det - _threshold );
  }

public:
  ResolveInterpenetrationEnergy( const typename ConfiguratorType::InitType &Grid, const std::vector<qc::CoordType> &ProblematicNodes, const RealType Threshold, const RealType Factor ) :
    _threshold( Threshold * ( ConfiguratorType::Dim == qc::QC_3D ? aol::Cub( Grid.H() ) : aol::Sqr( Grid.H() ) ) ),
    _factor( Factor ),
    _problematicNodes( ProblematicNodes ),
    _grid( Grid ) {
    const short elementsPerNode = 1 << ConfiguratorType::Dim;
    qc::BitArray< ConfiguratorType::Dim > elementHandled( qc::GridSize<ConfiguratorType::Dim>::createFrom( Grid ) );
    for ( unsigned int k = 0; k < _problematicNodes.size(); k++ )
      for ( int el = 0; el < elementsPerNode; el++ ) {
        qc::CoordType elem( _problematicNodes[k] );
        for ( int dir = 0; dir < ConfiguratorType::Dim; dir++ )
          elem[dir] = aol::Max( 0, aol::Min( _grid.getSize()[dir] - 2, elem[dir] - ( ( el >> dir ) % 2 ) ) );
        if ( !elementHandled.get( elem ) )
          _adjacentElements.push_back( elem );
        elementHandled.set( elem, true );
      }
  }

  void applyAdd( const aol::MultiVector<typename ConfiguratorType::RealType> &Arg, aol::Scalar<typename ConfiguratorType::RealType> &Dest ) const {
    qc::MultiArray<typename ConfiguratorType::RealType,ConfiguratorType::Dim,ConfiguratorType::Dim> arg( _grid, Arg, aol::FLAT_COPY );
    const short nodesPerElement = 1 << ConfiguratorType::Dim;
    for ( unsigned int k = 0; k < _adjacentElements.size(); k++ )
      for ( int node = 0; node < nodesPerElement; node++ ) {
        qc::CoordType coord( _adjacentElements[k][0] + node % 2, _adjacentElements[k][1] + ( node >> 1 ) % 2, _adjacentElements[k][2] + ( node >> 2 ) % 2 );
        typename ConfiguratorType::MatType dPhi;
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          qc::CoordType coord1( coord ), coord2( coord );
          coord1[j] = _adjacentElements[k][j];
          coord2[j] = _adjacentElements[k][j] + 1;
          for ( int i = 0; i < ConfiguratorType::Dim; i++ )
            dPhi[i][j] = arg[i].get( coord2 ) - arg[i].get( coord1 );
        }
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          dPhi[i][i] += _grid.H();
        Dest += evaluateNonlinearity( dPhi.det() );
      }
  }
};

template <typename ConfiguratorType>
class ResolveInterpenetrationGrad :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  const RealType _threshold;
  const RealType _factor;
  const std::vector<qc::CoordType> &_problematicNodes;
  std::vector<qc::CoordType> _adjacentElements;
  const typename ConfiguratorType::InitType &_grid;

  inline RealType evaluateNonlinearity( const RealType Det ) const {
    return ( Det > _threshold ) ? 0 : _factor * 2 * ( Det - _threshold );
  }

public:
  ResolveInterpenetrationGrad( const typename ConfiguratorType::InitType &Grid, const std::vector<qc::CoordType> &ProblematicNodes, const RealType Threshold, const RealType Factor ) :
    _threshold( Threshold * ( ConfiguratorType::Dim == qc::QC_3D ? aol::Cub( Grid.H() ) : aol::Sqr( Grid.H() ) ) ),
    _factor( Factor ),
    _problematicNodes( ProblematicNodes ),
    _grid( Grid ) {
    const short elementsPerNode = 1 << ConfiguratorType::Dim;
    qc::BitArray< ConfiguratorType::Dim > elementHandled( qc::GridSize<ConfiguratorType::Dim>::createFrom( Grid ) );
    for ( unsigned int k = 0; k < _problematicNodes.size(); k++ )
      for ( int el = 0; el < elementsPerNode; el++ ) {
        qc::CoordType elem( _problematicNodes[k] );
        for ( int dir = 0; dir < ConfiguratorType::Dim; dir++ )
          elem[dir] = aol::Max( 0, aol::Min( _grid.getSize()[dir] - 2, elem[dir] - ( ( el >> dir ) % 2 ) ) );
        if ( !elementHandled.get( elem ) )
          _adjacentElements.push_back( elem );
        elementHandled.set( elem, true );
      }
  }

  void applyAdd( const aol::MultiVector<typename ConfiguratorType::RealType> &Arg, aol::MultiVector<typename ConfiguratorType::RealType> &Dest ) const {
    qc::MultiArray<typename ConfiguratorType::RealType,ConfiguratorType::Dim,ConfiguratorType::Dim> arg( _grid, Arg, aol::FLAT_COPY ), grad( arg, aol::STRUCT_COPY ), dest( _grid, Dest, aol::FLAT_COPY );
    const short nodesPerElement = 1 << ConfiguratorType::Dim;
    for ( unsigned int k = 0; k < _adjacentElements.size(); k++ )
      for ( int node = 0; node < nodesPerElement; node++ ) {
        qc::CoordType coord( _adjacentElements[k][0] + node % 2, _adjacentElements[k][1] + ( node >> 1 ) % 2, _adjacentElements[k][2] + ( node >> 2 ) % 2 );
        typename ConfiguratorType::MatType dPhi;
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          qc::CoordType coord1( coord ), coord2( coord );
          coord1[j] = _adjacentElements[k][j];
          coord2[j] = _adjacentElements[k][j] + 1;
          for ( int i = 0; i < ConfiguratorType::Dim; i++ )
            dPhi[i][j] = arg[i].get( coord2 ) - arg[i].get( coord1 );
        }
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          dPhi[i][i] += _grid.H();
        typename ConfiguratorType::MatType deriv;
        deriv.makeCofactorMatrix( dPhi );
        deriv *= evaluateNonlinearity( dPhi.det() );
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          qc::CoordType coord1( coord ), coord2( coord );
          coord1[j] = _adjacentElements[k][j];
          coord2[j] = _adjacentElements[k][j] + 1;
          for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
            grad[i].add( coord2, deriv[i][j] );
            grad[i].add( coord1, - deriv[i][j] );
          }
        }
      }

    for ( unsigned int k = 0; k < _problematicNodes.size(); k++ )
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        dest[i].add( _problematicNodes[k], grad[i].get( _problematicNodes[k] ) );
  }
};

/**
 * Formally hyperelastic energy density used to remedy local material interpenetration.
 */
template <typename ConfiguratorType>
class AntiInterpenetrationHyperelasticEnergyDensity {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  const RealType _threshold;
  const RealType _factor;

public:
  AntiInterpenetrationHyperelasticEnergyDensity( const RealType Threshold, const RealType Factor ) :
    _threshold( Threshold ),
    _factor( Factor ) {}

  inline RealType evaluate ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                             const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    return ( I3 <= _threshold ? _factor * aol::Sqr( _threshold - I3 ) : 0. );
  }

  inline RealType evaluateAtQuadPoint ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                        const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return ( I3 <= _threshold ? _factor * aol::Sqr( _threshold - I3 ) : 0. );
  }

  inline aol::Vec3<RealType> evaluateDerivative ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    return aol::Vec3<RealType>( 0., 0., ( I3 <= _threshold ? 2 * _factor * ( I3 - _threshold ) : 0 ) );
  }

  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return aol::Vec3<RealType>( 0., 0., ( I3 <= _threshold ? 2 * _factor * ( I3 - _threshold ) : 0 ) );
  }
};

template <typename ConfiguratorType>
class CompressionAllowingEnergyDensity {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  // weights of the different energy components
  const RealType _weightLengthEnergy;
  const RealType _weightSurfaceEnergy;
  const RealType _weightVolumeEnergy;
  // Maximum compression allowed. Beyond, the energy is continued linearly.
  const RealType _compressionThreshold;
public:
  CompressionAllowingEnergyDensity( const RealType WeightLengthEnergy, const RealType WeightSurfaceEnergy, const RealType WeightVolumeEnergy ) :
    _weightLengthEnergy ( WeightLengthEnergy ),
    _weightSurfaceEnergy ( WeightSurfaceEnergy ),
    _weightVolumeEnergy ( WeightVolumeEnergy ),
    _compressionThreshold ( 1.e-2 ) {}

  inline RealType evaluate ( const RealType I1, const RealType /*I2*/, const RealType I3,
                             const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return _weightLengthEnergy * I1 / 2 + _weightVolumeEnergy * aol::Sqr( I3 ) / 4 - ( ConfiguratorType::Dim / 2. ) * _weightLengthEnergy - _weightVolumeEnergy / 4
            - ( _weightLengthEnergy + _weightVolumeEnergy / 2 )
              * ( I3 > _compressionThreshold
                ? log( I3 )
                : log( _compressionThreshold ) + ( I3 - _compressionThreshold ) / _compressionThreshold );
  }

  inline RealType evaluateAtQuadPoint ( const RealType I1, const RealType /*I2*/, const RealType I3,
                                        const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return _weightLengthEnergy * I1 / 2 + _weightVolumeEnergy * aol::Sqr( I3 ) / 4 - ( ConfiguratorType::Dim / 2. ) * _weightLengthEnergy - _weightVolumeEnergy / 4
            - ( _weightLengthEnergy + _weightVolumeEnergy / 2 )
              * ( I3 > _compressionThreshold
                ? log( I3 )
                : log( _compressionThreshold ) + ( I3 - _compressionThreshold ) / _compressionThreshold );
  }

  inline aol::Vec3<RealType> evaluateDerivative ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return aol::Vec3<RealType>( _weightLengthEnergy / 2,
                                0,
                                _weightVolumeEnergy * I3 / 2 - ( _weightLengthEnergy + _weightVolumeEnergy / 2 ) / aol::Max( I3, _compressionThreshold ) );
  }

  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return aol::Vec3<RealType>( _weightLengthEnergy / 2,
                                0,
                                _weightVolumeEnergy * I3 / 2 - ( _weightLengthEnergy + _weightVolumeEnergy / 2 ) / aol::Max( I3, _compressionThreshold ) );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivative ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                            const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return aol::Matrix33<RealType>( 0, 0, 0, 0, 0, 0, 0, 0, _weightVolumeEnergy / 2
                  + ( I3 > _compressionThreshold ? ( _weightLengthEnergy + _weightVolumeEnergy / 2 ) / aol::Sqr( I3 ) : 0. ) );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivativeAtQuadPoint ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                                       const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return aol::Matrix33<RealType>( 0, 0, 0, 0, 0, 0, 0, 0, _weightVolumeEnergy / 2
                  + ( I3 > _compressionThreshold ? ( _weightLengthEnergy + _weightVolumeEnergy / 2 ) / aol::Sqr( I3 ) : 0. ) );
  }
};

template <typename ConfiguratorType,typename HyperelasticEnergyDensityType>
class FirstPiolaKirchhoffStressOp :
  public aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim,FirstPiolaKirchhoffStressOp<ConfiguratorType,HyperelasticEnergyDensityType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;

public:
  FirstPiolaKirchhoffStressOp( const typename ConfiguratorType::InitType &Grid, const HyperelasticEnergyDensityType &EnergyDensity ) :
    aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim,FirstPiolaKirchhoffStressOp<ConfiguratorType,HyperelasticEnergyDensityType> >( Grid ),
    _hyperelasticEnergyDensity( EnergyDensity ) {}

  void getNonlinearity ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         typename aol::Vec<ConfiguratorType::Dim*ConfiguratorType::Dim,RealType> &NL ) const {

    // the deformation gradient
    typename ConfiguratorType::MatType dPhi;
    DiscFunc.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1;

    // compute stress
    typename ConfiguratorType::MatType firstPiolaKirchhoffStress;
    qc::HyperelasticGradient<ConfiguratorType,HyperelasticEnergyDensityType>::firstPiolaKirchhoffStress( dPhi, _hyperelasticEnergyDensity, El, QuadPoint, firstPiolaKirchhoffStress );

    // return result
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        NL[i*ConfiguratorType::Dim+j] = firstPiolaKirchhoffStress[i][j];
  }
};


//!Operator for evaluating Cauchy stress (i.e. force per unit area in deformed configuration)
/**Computes \f$ \left(\frac{\int_\Omega W_{,A}(D\phi) cof D\phi^{-1} \varphi_i(x)dx}{\int_\Omega\varphi_i(x)dx}\right)_i \f$,
 * where \f$\varphi_i\f$ is the \f$i\f$th finite element basis function and \f$\phi\f$ is passed to "apply" as argument. 
 * Note that \f$ W(D\phi) = W(I_1, I_2, I_3) \f$ with  \f$ I_1 = |D\phi|^2 \f$, \f$ I_2 = |cof D\phi|^2 \f$ and \f$ I_3 = \det D\phi \f$.
 *
 * \author Wirth
 */
template <typename ConfiguratorType,typename HyperelasticEnergyDensityType>
class CauchyStressOp :
  public aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim,CauchyStressOp<ConfiguratorType,HyperelasticEnergyDensityType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;

public:
  CauchyStressOp( const typename ConfiguratorType::InitType &Grid, const HyperelasticEnergyDensityType &EnergyDensity ) :
    aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim,CauchyStressOp<ConfiguratorType,HyperelasticEnergyDensityType> >( Grid ),
    _hyperelasticEnergyDensity( EnergyDensity ) {}

  void getNonlinearity ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         typename aol::Vec<ConfiguratorType::Dim*ConfiguratorType::Dim,RealType> &NL ) const {

    // the deformation gradient and its cofactor matrix
    typename ConfiguratorType::MatType dPhi, cof;
    DiscFunc.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1;
    cof.makeCofactorMatrix( dPhi );

    // compute the hyperelastic invariants $||\nabla\phi||_F^2,||cof(\nabla\phi)||_F^2,det(\nabla\phi)$
    RealType
      I1 = dPhi.normSqr(),
      I2 = cof.normSqr(),
      I3 = dPhi.det();

    // compute the outer derivative of the hyperelastic energy
    aol::Vec3<RealType> outerDeriv ( _hyperelasticEnergyDensity.evaluateDerivativeAtQuadPoint ( I1, I2, I3, El, QuadPoint ) );

    // compute the stress tensor
    typename ConfiguratorType::MatType cauchyStress, aux;
    // the length term
    aux = dPhi;
    aux.transpose();
    cauchyStress.makeProduct( dPhi, aux );
    cauchyStress *= 2 * outerDeriv[0] / I3;
    // the surface term
    aux = cof;
    aux.transpose();
    dPhi.makeProduct( cof, aux );
    RealType trace = dPhi.tr();
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] -= trace;
    dPhi *= -2 / I3 * outerDeriv[1];
    cauchyStress += dPhi;
    // the volume term
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      cauchyStress[i][i] += outerDeriv[2];

    // return result
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        NL[i*ConfiguratorType::Dim+j] = cauchyStress[i][j];
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // compute the Cauchy stress in the reference configuration
    aol::MultiVector<RealType> undefDest( Dest, aol::STRUCT_COPY ), dest( Dest, aol::STRUCT_COPY );
    aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim,CauchyStressOp<ConfiguratorType,HyperelasticEnergyDensityType> >::applyAdd( Arg, undefDest );

    // deform the Cauchy stress to the deformed configuration
    qc::InvDeformImage<ConfiguratorType,aol::MultiVector<RealType> >( undefDest, this->_initializer, dest, Arg );
    Dest += dest;
  }
};

template<typename ConfiguratorType>
class Cauchy2FirstPiolaKirchhoff :
  public aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim,Cauchy2FirstPiolaKirchhoff<ConfiguratorType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;

public:
  Cauchy2FirstPiolaKirchhoff( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement ) :
    aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim,Cauchy2FirstPiolaKirchhoff<ConfiguratorType> >( Grid ),
    _displacement( Grid, Displacement ) {}

  inline void getNonlinearity( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim> &DiscFunc,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                               aol::Vec<ConfiguratorType::Dim*ConfiguratorType::Dim,RealType> &NL ) const {
    NL.setZero();

    // compute the factor relating first Piola Kirchhoff and Cauchy stress
    typename ConfiguratorType::MatType defGrad, cof;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, defGrad );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      defGrad[i][i] += 1.;
    cof.makeCofactorMatrix( defGrad );

    // get Cauchy stress in the reference configuration
    aol::Vec<ConfiguratorType::Dim*ConfiguratorType::Dim,RealType> cauchyStress;
    typename ConfiguratorType::VecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    if ( !qc::transformCoord<ConfiguratorType>( this->_initializer, _displacement, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord ) )
      cauchyStress.setZero();
    else
      DiscFunc.evaluate( transformedEl, transformedLocalCoord, cauchyStress );

    // rescale and reorientate to obtain first Piola Kirchhoff stress
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        for ( int k = 0; k < ConfiguratorType::Dim; k++ )
          NL[i*ConfiguratorType::Dim+j] += cauchyStress[i*ConfiguratorType::Dim+k] * cof[k][j];
  }
};

template<typename ConfiguratorType>
class StressMultWeightedNormal :
  public aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim,ConfiguratorType::Dim,StressMultWeightedNormal<ConfiguratorType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const aol::DiscreteFunctionDefault<ConfiguratorType> _signedDistFunc;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _weight;

public:
  StressMultWeightedNormal( const typename ConfiguratorType::InitType &Grid, const aol::Vector<RealType> &SignedDistFunc, const aol::Vector<RealType> &Weight ) :
    aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim,ConfiguratorType::Dim,StressMultWeightedNormal<ConfiguratorType> >( Grid ),
    _signedDistFunc( Grid, SignedDistFunc ),
    _weight( Grid, Weight ) {}

  void getNonlinearity ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         typename aol::Vec<ConfiguratorType::Dim,RealType> &NL ) const {
    NL.setZero();
    aol::Vec<ConfiguratorType::Dim,RealType> normal;
    _signedDistFunc.evaluateGradientAtQuadPoint( El, QuadPoint, normal );
    aol::Vec<ConfiguratorType::Dim*ConfiguratorType::Dim,RealType> stress;
    DiscFunc.evaluateAtQuadPoint( El, QuadPoint, stress );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        NL[i] += stress[i*ConfiguratorType::Dim+j] * normal[j];
    NL *= _weight.evaluateAtQuadPoint( El, QuadPoint );
  }
};

template<typename ConfiguratorType>
class DiffusedBoundaryL2Product :
  public aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim, DiffusedBoundaryL2Product<ConfiguratorType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  const typename ConfiguratorType::InitType &_grid;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _weight;
  const ArrayType &_kernel;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _signedDist;
  ArrayType _scaling;
  const qc::LevelsetMassOp<ConfiguratorType> _l2Product;
  const StressMultWeightedNormal<ConfiguratorType> _normalOp;
  const qc::Convolution<ConfiguratorType::Dim> _convolution;


  inline typename aol::VecDimTrait<int,ConfiguratorType::Dim>::VecType fitVec( const aol::Vec3<int> &Vector ){
    typename aol::VecDimTrait<int,ConfiguratorType::Dim>::VecType result;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      result[i] = Vector[i];
    return result;
  }

public:
  DiffusedBoundaryL2Product( const typename ConfiguratorType::InitType &Grid,
                             const ArrayType &SignedDist,
                             const aol::Vector<RealType> &Weight,
                             const ArrayType &Kernel ) :
    aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim, DiffusedBoundaryL2Product<ConfiguratorType> >( Grid ),
    _grid( Grid ),
    _weight( Grid, Weight ),
    _kernel( Kernel ),
    _signedDist( Grid, SignedDist ),
    _scaling( Grid ),
    _l2Product( Grid, SignedDist ),
    _normalOp( Grid, SignedDist, Weight ),
    _convolution( fitVec( Grid.getSize() ) ) {
    const ArrayType weight( Weight, Grid );
    _convolution.convolve( weight, _kernel, _scaling );save<ConfiguratorType>(_scaling,_grid,"Yscaling");
    /*_scaling *= 1. / Grid.getNumberOfElements(); // not needed, if this scaling is also left out in applyAdd*/
    for ( int j = 0; j < _scaling.size(); j++ )
      if ( aol::Abs( _scaling[j] ) < 1.e-10 )
        _scaling[j] = 1.;
  }

  void getNonlinearity ( aol::auto_container<ConfiguratorType::Dim*ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         aol::Vec<ConfiguratorType::Dim*ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {
    typename ConfiguratorType::VecType normal;
    _signedDist.evaluateGradientAtQuadPoint( El, QuadPoint, normal );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        NL[i*ConfiguratorType::Dim+j] = normal[j] * DiscFuncs[i*ConfiguratorType::Dim+j].evaluateAtQuadPoint( El, QuadPoint );
    NL *= _weight.evaluateAtQuadPoint( El, QuadPoint );
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> normalArg( _grid ), normalConvArg( _grid ),
                               auxMultiVec, auxMultiVecComponents( _grid );
    _normalOp.apply( Arg, normalArg );save<ConfiguratorType>(normalArg[0],_grid,"YnormalArg0");save<ConfiguratorType>(normalArg[1],_grid,"YnormalArg1");
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      ArrayType arg( normalArg[i], _grid, aol::FLAT_COPY ), res( normalConvArg[i], _grid, aol::FLAT_COPY );
      _convolution.convolve( arg, _kernel, res );save<ConfiguratorType>(normalConvArg[0],_grid,"YnormalConvArg0");save<ConfiguratorType>(normalConvArg[1],_grid,"YnormalConvArg1");
      /*res *= 1. / _grid.getNumberOfElements(); // not needed, if this scaling is also left out in constructor*/
      for ( int j = 0; j < normalConvArg[i].size(); j++ )
        normalConvArg[i][j] /= /*aol::Sqr( */_scaling[j]/* )*/;save<ConfiguratorType>(normalConvArg[0],_grid,"YnormalConvScArg0");save<ConfiguratorType>(normalConvArg[1],_grid,"YnormalConvScArg1");
      ArrayType aux( auxMultiVecComponents[i], _grid, aol::FLAT_COPY ), auxMultiVecComp( _grid );
      _l2Product.apply( normalConvArg[i], aux );save<ConfiguratorType>(auxMultiVecComponents[0],_grid,"Yl20");save<ConfiguratorType>(auxMultiVecComponents[1],_grid,"Yl21");
      _convolution.convolve( aux, _kernel, auxMultiVecComp );save<ConfiguratorType>(auxMultiVecComp,_grid,"Yconvl21");
      auxMultiVecComp *= _grid.getNumberOfElements(); // if the scaling is left out everywhere else, this has to be scaled
      /*auxMultiVecComp *= 1. / _grid.getNumberOfElements(); // not needed, if this scaling is also left out in constructor*/
      for ( int x = 0; x < aux.getNumX(); x++ )
        for ( int y = 0; y < aux.getNumY(); y++ )
          aux.set( x, y, auxMultiVecComp.get( /*aux.getNumX() -*/ x, /*aux.getNumY() -*/ y ) );
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        auxMultiVec.appendReference( auxMultiVecComponents[i] );
    }

    aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim, DiffusedBoundaryL2Product<ConfiguratorType> >::applyAdd( auxMultiVec, Dest );
  }
};



#endif
