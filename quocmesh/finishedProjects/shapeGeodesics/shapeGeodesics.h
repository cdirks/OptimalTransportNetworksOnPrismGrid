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

#ifndef __SHAPEGEODESICS_H
#define __SHAPEGEODESICS_H

#include "../shapeStatistics/testResultSaving.h"
#include "../shapeStatistics/auxiliaryOps.h"
#include "../shapeStatistics/stressOps.h"
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

//#define DEBUGGING_OUTPUT
//#define USE_LANDMARKS
#define NUM_LEVELSETS 1

bool TEST = false;
bool TEST1 = false;
int TEST_INDEX = 1000;


/**************************************************************
 * The averaging energy of one single deformation and its variation.
 **************************************************************/

inline const qc::GridDefinition &cubicGrid( const qc::GridDefinition &Grid ) {
  return Grid;
}

template <qc::Dimension Dim>
inline const qc::GridDefinition &cubicGrid( const qc::simplex::GridStructure<qc::GridDefinition,Dim> &Grid ) {
  return Grid.getCubicGrid();
}

template <typename ConfiguratorType, typename HeavisideFunctionType>
class VolumeOfZeroSubLevelSetOp:
  public  aol::FENonlinIntegrationScalarInterface<ConfiguratorType,VolumeOfZeroSubLevelSetOp<ConfiguratorType,HeavisideFunctionType> > {

private:
  const HeavisideFunctionType _heavisideFunction;

public:
  VolumeOfZeroSubLevelSetOp( const typename ConfiguratorType::InitType &Grid,
                             const HeavisideFunctionType HeavisideFunction ) :
    aol::FENonlinIntegrationScalarInterface<ConfiguratorType, VolumeOfZeroSubLevelSetOp<ConfiguratorType,HeavisideFunctionType> >( Grid ),
    _heavisideFunction( HeavisideFunction ) {}

  typename ConfiguratorType::RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                                                          const typename ConfiguratorType::ElementType &El,
                                                          int QuadPoint,
                                                          const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return _heavisideFunction.evaluate( DiscFunc.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class VolumeOfZeroSubLevelSetVariation:
  public  aol::FENonlinOpInterface<ConfiguratorType,VolumeOfZeroSubLevelSetVariation<ConfiguratorType,HeavisideFunctionType> > {

private:
  const HeavisideFunctionType _heavisideFunction;

public:
  VolumeOfZeroSubLevelSetVariation( const typename ConfiguratorType::InitType &Grid,
                                    const HeavisideFunctionType HeavisideFunction ) :
    aol::FENonlinOpInterface<ConfiguratorType,VolumeOfZeroSubLevelSetVariation<ConfiguratorType,HeavisideFunctionType> >( Grid ),
    _heavisideFunction( HeavisideFunction ) {}

  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                         typename ConfiguratorType::RealType &NL ) const {
    NL = _heavisideFunction.evaluateDerivative( DiscFunc.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

template <typename ConfiguratorType>
class LandmarkMismatchEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const ConfiguratorType _config;
  const RealType _lambda;

public:
  LandmarkMismatchEnergy( const typename ConfiguratorType::InitType &Grid, const RealType Lambda ) :
    _grid( Grid ),
    _config( Grid ),
    _lambda( Lambda ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector<RealType> displacement, landmarksLeft, landmarksRight;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      displacement.appendReference( Arg[i] );
      landmarksLeft.appendReference( Arg[i+ConfiguratorType::Dim] );
      landmarksRight.appendReference( Arg[i+2*ConfiguratorType::Dim] );
    }
    aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> displ( _grid, displacement );

    for ( int j = 0; j < landmarksLeft[0].size(); j++ ) {
      bool insideFlag = true;
      typename ConfiguratorType::VecType landmarkLeft, landmarkRight, offset;
      for ( int k = 0; k < ConfiguratorType::Dim; k++ ) {
        landmarkLeft[k] = landmarksLeft[k][j];
        landmarkRight[k] = landmarksRight[k][j];
        insideFlag = insideFlag && ( landmarkLeft[k] >= 0. ) && ( landmarkLeft[k] <= 1. );
      }
      if ( insideFlag ) {
        typename ConfiguratorType::ElementType el;
        typename ConfiguratorType::DomVecType localCoord;
        _config.getLocalCoords( landmarkLeft, el, localCoord );
        displ.evaluate( el, localCoord, offset );
        landmarkLeft += offset;
        landmarkLeft -= landmarkRight;
        Dest.v += landmarkLeft.normSqr();
      } else
        Dest.v += aol::NumberTrait<RealType>::Inf;
    }
  }
};

template <typename ConfiguratorType>
class LandmarkMismatchGradient :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const ConfiguratorType _config;
  const RealType _lambda;

public:
  LandmarkMismatchGradient( const typename ConfiguratorType::InitType &Grid, const RealType Lambda ) :
    _grid( Grid ),
    _config( Grid ),
    _lambda( Lambda ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> displacement, landmarksLeft, landmarksRight, displacementGradient, landmarksLeftGradient, landmarksRightGradient;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      displacement.appendReference( Arg[i] );
      landmarksLeft.appendReference( Arg[i+ConfiguratorType::Dim] );
      landmarksRight.appendReference( Arg[i+2*ConfiguratorType::Dim] );
      displacementGradient.appendReference( Dest[i] );
      landmarksLeftGradient.appendReference( Dest[i+ConfiguratorType::Dim] );
      landmarksRightGradient.appendReference( Dest[i+2*ConfiguratorType::Dim] );
    }
    aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> displ( _grid, displacement );

    for ( int j = 0; j < landmarksLeft[0].size(); j++ ) {
      // obtain the landmark deviation
      typename ConfiguratorType::VecType landmarkLeft, landmarkRight, offset;
      for ( int k = 0; k < ConfiguratorType::Dim; k++ ) {
        landmarkLeft[k] = landmarksLeft[k][j];
        landmarkRight[k] = landmarksRight[k][j];
        if ( ( landmarkLeft[k] < 0. ) || ( landmarkLeft[k] > 1. ) )
          throw aol::Exception( "landmarks have to lie within [0,1]^d" );;
      }
      typename ConfiguratorType::ElementType el;
      typename ConfiguratorType::DomVecType localCoord;
      _config.getLocalCoords( landmarkLeft, el, localCoord );
      displ.evaluate( el, localCoord, offset );
      landmarkLeft += offset;
      landmarkLeft -= landmarkRight;
      // obtain the derivative wrt right landmarks
      landmarkLeft *= 2 * _lambda;
      for ( int k = 0; k < ConfiguratorType::Dim; k++ )
        landmarksRightGradient[k][j] -= landmarkLeft[k];
      // obtain the derivative wrt left landmarks
      typename ConfiguratorType::MatType defGrad;
      displ.evaluateGradient( el, localCoord, defGrad );
      for ( int k = 0; k < ConfiguratorType::Dim; k++ )
        defGrad[k][k] += 1;
      defGrad.transpose();
      defGrad.mult( landmarkLeft, offset );
      for ( int k = 0; k < ConfiguratorType::Dim; k++ )
        landmarksLeftGradient[k][j] = offset[k];
      // obtain the derivative wrt displacement
      for ( int dof = 0; dof < _config.getNumLocalDofs( el ); dof++ ) {
        const RealType bfsVal = _config.getBaseFunctionSet( el ).evaluate( dof, localCoord );
        for ( int d = 0; d < ConfiguratorType::Dim; d++ )
          displacementGradient[d][_config.localToGlobal( el, dof )] += landmarkLeft[d] * bfsVal;
      }
    }
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class LevelsetMismatchEnergyOfDisp :
  public aol::FENonlinIntegrationVectorInterface<ConfiguratorType,LevelsetMismatchEnergyOfDisp<ConfiguratorType,HeavisideFunctionType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  // the fixed and deformable levelset function to be matched
  const aol::DiscreteFunctionDefault<ConfiguratorType> _uFixed, _uTemplate;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // indicates, whether the fixed shape may be partially occluded or not
  const bool _occlusion;

public:
  LevelsetMismatchEnergyOfDisp( const typename ConfiguratorType::InitType &Grid,
                                const aol::Vector<RealType> &UFixed,
                                const aol::Vector<RealType> &UTemplate,
                                const HeavisideFunctionType &HeavisideFunction,
                                const bool Occlusion = false ) :
    aol::FENonlinIntegrationVectorInterface<ConfiguratorType,LevelsetMismatchEnergyOfDisp<ConfiguratorType,HeavisideFunctionType> >( Grid ),
    _uFixed( Grid, UFixed ),
    _uTemplate( Grid, UTemplate ),
    _heavisideFunction( HeavisideFunction ),
    _occlusion( Occlusion ) {}

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El, int QuadPoint,
                              const typename ConfiguratorType::DomVecType &RefCoord ) const {
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    // attention! Clipping is important to extend the image "continuously" outside the computational domain and thus prevent energy discontinuities wrt displacements.
    qc::transformAndClipCoord<ConfiguratorType> ( *this->_config, DiscFuncs, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord );

    RealType result = aol::Sqr( _heavisideFunction.evaluate( _uTemplate.evaluate( transformedEl, transformedLocalCoord ) )
                      - _heavisideFunction.evaluate( _uFixed.evaluateAtQuadPoint( El, QuadPoint ) ) );

    if ( _occlusion )
      // multiply result with a heaviside approximation that is zero for negative _uFixed and relatively low for _uFixed small
      //result *= aol::Sqr( _heavisideFunction.evaluate( _uFixed.evaluateAtQuadPoint( El, QuadPoint ) ) );
      result *= aol::Min( 1., aol::Max( 0., _uFixed.evaluateAtQuadPoint( El, QuadPoint ) / this->_initializer.H() - .5 ) );

    return result;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class LevelsetMismatchGradientWRTDisp :
  public aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim,LevelsetMismatchGradientWRTDisp<ConfiguratorType,HeavisideFunctionType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  // the fixed and deformable levelset function to be matched
  const aol::DiscreteFunctionDefault<ConfiguratorType> _uFixed, _uTemplate;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // indicates, whether the fixed shape may be partially occluded or not
  const bool _occlusion;

public:
  LevelsetMismatchGradientWRTDisp( const typename ConfiguratorType::InitType &Grid,
                                   const aol::Vector<RealType> &UFixed,
                                   const aol::Vector<RealType> &UTemplate,
                                   const HeavisideFunctionType &HeavisideFunction,
                                   const bool Occlusion = false ) :
    aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim,LevelsetMismatchGradientWRTDisp<ConfiguratorType,HeavisideFunctionType> >( Grid ),
    _uFixed( Grid, UFixed ),
    _uTemplate( Grid, UTemplate ),
    _heavisideFunction( HeavisideFunction ),
    _occlusion( Occlusion ) {}

  void getNonlinearity( aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El, int QuadPoint,
                        const typename ConfiguratorType::DomVecType &RefCoord,
                        aol::Vec<ConfiguratorType::Dim,RealType > &NL) const {
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    aol::Vec<ConfiguratorType::Dim,bool> coordinateWithinLimits;
    qc::transformAndClipCoord<ConfiguratorType> ( *this->_config, DiscFuncs, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord, coordinateWithinLimits );

    RealType uTemplate = _uTemplate.evaluate( transformedEl, transformedLocalCoord );
    _uTemplate.evaluateGradient( transformedEl, transformedLocalCoord, NL );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      if ( coordinateWithinLimits[i] == false )
        NL[i] = 0.;

    NL *= 2. * ( _heavisideFunction.evaluate( uTemplate ) - _heavisideFunction.evaluate( _uFixed.evaluateAtQuadPoint( El, QuadPoint ) ) )
              * _heavisideFunction.evaluateDerivative( uTemplate );

    if ( _occlusion )
      // multiply result with a heaviside approximation that is zero for negative _uFixed and relatively low for _uFixed small
      //NL *= aol::Sqr( _heavisideFunction.evaluate( _uFixed.evaluateAtQuadPoint( El, QuadPoint ) ) );
      NL *= aol::Min( 1., aol::Max( 0., _uFixed.evaluateAtQuadPoint( El, QuadPoint ) / this->_initializer.H() - .5 ) );
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class LevelsetMismatchGradientWRTSDF :
  public aol::FENonlinOpInterface<ConfiguratorType,LevelsetMismatchGradientWRTSDF<ConfiguratorType,HeavisideFunctionType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  // the levelset function to be matched with
  const aol::DiscreteFunctionDefault<ConfiguratorType> _uTemplate;
  // the matching displacement
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // indicates, whether the fixed shape may be partially occluded or not
  const bool _occlusion;

public:
  LevelsetMismatchGradientWRTSDF( const typename ConfiguratorType::InitType &Grid,
                                  const aol::Vector<RealType> &LevelsetFunction,
                                  const aol::MultiVector<RealType> &Displacement,
                                  const HeavisideFunctionType &HeavisideFunction,
                                  const bool Occlusion = false ) :
    aol::FENonlinOpInterface<ConfiguratorType,LevelsetMismatchGradientWRTSDF<ConfiguratorType,HeavisideFunctionType> >( Grid ),
    _uTemplate( Grid, LevelsetFunction ),
    _displacement( Grid, Displacement ),
    _heavisideFunction( HeavisideFunction ),
    _occlusion( Occlusion ) {}

  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         const int QuadPoint,
                         const typename ConfiguratorType::DomVecType &LocalCoord, RealType &NL ) const {
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    // attention! Clipping is important to extend the image "continuously" outside the computational domain and thus prevent energy discontinuities wrt displacements.
    qc::transformAndClipCoord<ConfiguratorType> ( *this->_config, _displacement, El, QuadPoint, LocalCoord, transformedEl, transformedLocalCoord );

    const RealType uFixed = DiscFunc.evaluateAtQuadPoint( El, QuadPoint );
    const RealType uTemplate = _uTemplate.evaluate( transformedEl, transformedLocalCoord );

    NL = 2. * ( _heavisideFunction.evaluate( uFixed ) - _heavisideFunction.evaluate( uTemplate ) ) * _heavisideFunction.evaluateDerivative( uFixed );

    if ( _occlusion ){
      // use a heaviside approximation that is zero for negative _uFixed and relatively low for _uFixed small
      //NL *= _heavisideFunction.evaluate( uFixed ) * (2. * _heavisideFunction.evaluate( uFixed ) - _heavisideFunction.evaluate( uTemplate ) );
      RealType aux = uFixed / this->_initializer.H() - .5;
      NL *= aol::Min( aol::ZOTrait<RealType>::one, aol::Max( aol::ZOTrait<RealType>::zero, aux ) );
      NL += aol::Sqr( _heavisideFunction.evaluate( uFixed ) - _heavisideFunction.evaluate( uTemplate ) )
            * ( aux < 0. || aux > 1. ? 0. : 1 / this->_initializer.H() );
    }
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class LevelsetMismatchGradientWRTDeformedSDF :
  public qc::FENonlinDeformOpInterface<ConfiguratorType,LevelsetMismatchGradientWRTDeformedSDF<ConfiguratorType,HeavisideFunctionType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  // the levelset function to be matched with
  const aol::DiscreteFunctionDefault<ConfiguratorType> _uFixed;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // indicates, whether the fixed shape may be partially occluded or not
  const bool _occlusion;

public:
  LevelsetMismatchGradientWRTDeformedSDF( const typename ConfiguratorType::InitType &Grid,
                                          const aol::Vector<RealType> &LevelsetFunction,
                                          const aol::MultiVector<RealType> &Displacement,
                                          const HeavisideFunctionType &HeavisideFunction,
                                          const bool Occlusion = false ) :
    qc::FENonlinDeformOpInterface<ConfiguratorType,LevelsetMismatchGradientWRTDeformedSDF<ConfiguratorType,HeavisideFunctionType> >( Grid, Displacement ),
    _uFixed( Grid, LevelsetFunction ),
    _heavisideFunction( HeavisideFunction ),
    _occlusion( Occlusion ) {}

  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                         const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                         typename ConfiguratorType::RealType &NL ) const {
    const RealType uFixed = _uFixed.evaluateAtQuadPoint( El, QuadPoint );
    const RealType uTemplate = DiscFunc.evaluate( TransformedEl, TransformedLocalCoord );

    NL = 2. * ( _heavisideFunction.evaluate( uTemplate ) - _heavisideFunction.evaluate( uFixed ) ) * _heavisideFunction.evaluateDerivative( uTemplate );

    if ( _occlusion )
      // multiply result with a heaviside approximation that is zero for negative _uFixed and relatively low for _uFixed small
      //NL *= aol::Sqr( _heavisideFunction.evaluate( uFixed ) );
      NL *= aol::Min( 1., aol::Max( 0., uFixed / this->_grid.H() - .5 ) );
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, typename MaterialLawType>
class LevelsetGeodesicEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,MaterialLawType> WeightedEnergyDensityType;

  const typename ConfiguratorType::InitType &_grid;
  const int _numPolygonLines, _numShapes;
  const RealType _gamma, _mu, _nu, _lambda;
  const HeavisideFunctionType _heavisideFunction;
  const MaterialLawType _elasticEnergyDensity;
  const aol::MultiVector<RealType> &_weight;
  const aol::HeavisideLevelsetLengthEnergy<ConfiguratorType,HeavisideFunctionType,1> _lengthEnergyOp;
  const VolumeOfZeroSubLevelSetOp<ConfiguratorType,HeavisideFunctionType> _volumeOp;
  // indicates, whether the first shape may be partially occluded or not
  const bool _occlusion;

public:
  LevelsetGeodesicEnergy( const typename ConfiguratorType::InitType &Grid,
                          const int NumberOfDeformations,
                          const RealType Gamma,
                          const RealType Mu,
                          const RealType Nu,
                          const RealType Lambda,
                          const RealType NormRegularization,
                          const RealType Epsilon,
                          const RealType LengthEnergyWeight,
                          const RealType VolumeEnergyWeight,
                          const aol::MultiVector<RealType> &Weight,
                          const bool Occlusion = false ) :
    _grid( Grid ),
    _numPolygonLines( NumberOfDeformations ),
    _numShapes( _numPolygonLines + 1 ),
    _gamma( Gamma ),
    _mu( Mu ),
    _nu( Nu ),
    _lambda( Lambda ),
    _heavisideFunction( Epsilon ),
    _elasticEnergyDensity( LengthEnergyWeight, 0, VolumeEnergyWeight ),
    _weight( Weight ),
    _lengthEnergyOp( _grid, _heavisideFunction, NormRegularization ),
    _volumeOp( _grid, _heavisideFunction ),
    _occlusion( Occlusion ) {}

  // Arg has (2*dim+1)*n+dim+1 components, ie. n+1 levelsets, dim*n displacement components, dim*(n+1) landmark component vectors
  void computeEnergyComponents( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &EnergyComponents ) const {
    aol::MultiVector<RealType> FidelityEnergies, LengthRegularizationEnergies, VolumeDeviationEnergies;
    for ( int i = 0; i < NUM_LEVELSETS; i++ ) {
      FidelityEnergies.appendReference( EnergyComponents[i], false );
      LengthRegularizationEnergies.appendReference( EnergyComponents[NUM_LEVELSETS+1+i], false );
      VolumeDeviationEnergies.appendReference( EnergyComponents[2*NUM_LEVELSETS+1+i], false );
    }
    aol::Vector<RealType> DeformationEnergies( EnergyComponents[NUM_LEVELSETS], aol::FLAT_COPY );
#ifdef USE_LANDMARKS
    aol::Vector<RealType> LandmarkDeviationEnergies( EnergyComponents[NUM_LEVELSETS+3], aol::FLAT_COPY );
#endif

    // compute the volumes of the start and end shapes
    RealType volumeOfEndShapes[NUM_LEVELSETS][2];
    if ( _nu > 0. )
      for ( int i = 0; i < NUM_LEVELSETS; i++ ) {
        aol::Scalar<RealType> volumeOfEndShape;
        _volumeOp.apply( Arg[i*_numShapes+0], volumeOfEndShape );
        volumeOfEndShapes[i][0] = volumeOfEndShape;
        _volumeOp.apply( Arg[(i+1)*_numShapes-1], volumeOfEndShape );
        volumeOfEndShapes[i][1] = volumeOfEndShape;
      }
    else
      for ( int i = 0; i < NUM_LEVELSETS; i++ )
        for ( int j = 0; j < 2; j++ )
          volumeOfEndShapes[i][j] = aol::NumberTrait<RealType>::NaN;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numPolygonLines; i++ ) {
      aol::MultiVector<RealType> sdfLeft, sdfRight;
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        sdfLeft.appendReference( Arg[j*_numShapes+i], false );
        sdfRight.appendReference( Arg[j*_numShapes+i+1], false );
      }
      aol::MultiVector<RealType> displacement, landmarksLeft, landmarksRight;
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        displacement.appendReference( Arg[NUM_LEVELSETS*_numShapes+i*ConfiguratorType::Dim+j], false );
#ifdef USE_LANDMARKS
        landmarksLeft.appendReference( Arg[NUM_LEVELSETS*_numShapes+_numPolygonLines*ConfiguratorType::Dim+i*ConfiguratorType::Dim+j], false );
        landmarksRight.appendReference( Arg[NUM_LEVELSETS*_numShapes+_numPolygonLines*ConfiguratorType::Dim+(i+1)*ConfiguratorType::Dim+j], false );
#endif
      }

      aol::Scalar<RealType> destComp;

      // compute deformation energy
      WeightedEnergyDensityType weightedDensity( _elasticEnergyDensity, _grid, _weight[i] );
      qc::HyperelasticEnergy<ConfiguratorType,WeightedEnergyDensityType>( _grid, weightedDensity ).apply( displacement, destComp );
      DeformationEnergies[i] = destComp.v;

      // compute mismatch energies
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        LevelsetMismatchEnergyOfDisp<ConfiguratorType,HeavisideFunctionType>( _grid, sdfLeft[j], sdfRight[j], _heavisideFunction, _occlusion && i == 0 ).apply( displacement, destComp );
        FidelityEnergies[j][i] = destComp.v;
      }

      // compute perimeter regularization
      if ( _mu > 0. && i > 0 )
        for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
          aol::MultiVector<RealType> arg;
          arg.appendReference( sdfLeft[j], false );
          _lengthEnergyOp.apply( arg, destComp );
          LengthRegularizationEnergies[j][i] = destComp.v;
        }

      // compute volume deviation energy
      if ( _nu > 0. && i > 0 )
        for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
          _volumeOp.apply( sdfLeft[j], destComp );
          VolumeDeviationEnergies[j][i] = aol::Sqr( destComp.v - ( ( ( volumeOfEndShapes[j][1] - volumeOfEndShapes[j][0] ) * i ) / _numPolygonLines + volumeOfEndShapes[j][0] ) );
        }

#ifdef USE_LANDMARKS
      // compute landmark energy
      displacement.appendReference( landmarksLeft );
      displacement.appendReference( landmarksRight );
      LandmarkMismatchEnergy<ConfiguratorType>( _grid, _lambda ).apply( displacement, destComp );
      LandmarkDeviationEnergies[i] = destComp.v;
#endif
    }
  }

  // Arg has (2*dim+1)*n+dim+1 components, ie. n+1 levelsets, dim*n displacement components, dim*(n+1) landmark component vectors
  void computeEnergyDensities( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &FidelityEnergy, aol::MultiVector<RealType> &DeformationEnergy, aol::MultiVector<RealType> &LengthRegularizationEnergy, aol::MultiVector<RealType> &VolumeDeviationEnergy ) const {
    FidelityEnergy.reallocate( NUM_LEVELSETS*_numPolygonLines, Arg[0].size() );
    DeformationEnergy.reallocate( _numPolygonLines, Arg[0].size() );
    LengthRegularizationEnergy.reallocate( NUM_LEVELSETS*_numPolygonLines, Arg[0].size() );
    VolumeDeviationEnergy.reallocate( NUM_LEVELSETS*_numPolygonLines, Arg[0].size() );

    // compute the volumes of the start and end shapes
    RealType volumeOfEndShapes[NUM_LEVELSETS][2];
    if ( _nu > 0. )
      for ( int i = 0; i < NUM_LEVELSETS; i++ ) {
        aol::Scalar<RealType> volumeOfEndShape;
        _volumeOp.apply( Arg[i*_numShapes+0], volumeOfEndShape );
        volumeOfEndShapes[i][0] = volumeOfEndShape;
        _volumeOp.apply( Arg[(i+1)*_numShapes-1], volumeOfEndShape );
        volumeOfEndShapes[i][1] = volumeOfEndShape;
      }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numPolygonLines; i++ ) {
      aol::MultiVector<RealType> sdfLeft, sdfRight;
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        sdfLeft.appendReference( Arg[j*_numShapes+i], false );
        sdfRight.appendReference( Arg[j*_numShapes+i+1], false );
      }
      aol::MultiVector<RealType> displacement;
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        displacement.appendReference( Arg[NUM_LEVELSETS*_numShapes+i*ConfiguratorType::Dim+j], false );

      // compute deformation energy
      WeightedEnergyDensityType weightedDensity( _elasticEnergyDensity, _grid, _weight[i] );
      qc::HyperelasticEnergy<ConfiguratorType,WeightedEnergyDensityType>( _grid, weightedDensity ).applyAddIntegrand( displacement, DeformationEnergy[i] );

      // compute mismatch energies
      for ( int j = 0; j < NUM_LEVELSETS; j++ )
        LevelsetMismatchEnergyOfDisp<ConfiguratorType,HeavisideFunctionType>( _grid, sdfLeft[j], sdfRight[j], _heavisideFunction, _occlusion && i == 0 ).applyAddIntegrand( displacement, FidelityEnergy[j*_numPolygonLines+i] );

      // compute perimeter regularization
      if ( _mu > 0. && i > 0 )
        for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
          aol::MultiVector<RealType> arg;
          arg.appendReference( sdfLeft[j], false );
          _lengthEnergyOp.applyAddIntegrand( arg, LengthRegularizationEnergy[j*_numPolygonLines+i] );
        }

      // compute volume deviation energy
      if ( _nu > 0. && i > 0 )
        for ( int j = 0; j < NUM_LEVELSETS; j++ )
          _volumeOp.applyAddIntegrand( sdfLeft[j], VolumeDeviationEnergy[j*_numPolygonLines+i] );
    }
  }

  // Arg has (dim+1)*n+1 components, ie. n+1 levelsets and dim*n displacement components
  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector<RealType> EnergyComponents( 3 * NUM_LEVELSETS + 2, _numPolygonLines );
    computeEnergyComponents( Arg, EnergyComponents );
    Dest += EnergyComponents[NUM_LEVELSETS].sum();
    for ( int i = 0; i < NUM_LEVELSETS; i++ )
      Dest += _gamma * EnergyComponents[i].sum() + _mu * EnergyComponents[NUM_LEVELSETS+1+i].sum() + _nu * EnergyComponents[2*NUM_LEVELSETS+1+i].sum();
#ifdef USE_LANDMARKS
    Dest += _lambda * EnergyComponents[3*NUM_LEVELSETS+1].sum();
#endif
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, typename MaterialLawType>
class LevelsetGeodesicGradient :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,MaterialLawType> WeightedEnergyDensityType;

  const typename ConfiguratorType::InitType &_grid;
  const int _numPolygonLines, _numShapes;
  const RealType _gamma, _mu, _nu, _lambda;
  const HeavisideFunctionType _heavisideFunction;
  const MaterialLawType _elasticEnergyDensity;
  const aol::MultiVector<RealType> &_weight;
  const aol::FullVariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType,HeavisideFunctionType,1> _lengthEnergyVariation;
  const VolumeOfZeroSubLevelSetOp<ConfiguratorType,HeavisideFunctionType> _volumeOp;
  const VolumeOfZeroSubLevelSetVariation<ConfiguratorType,HeavisideFunctionType> _volumeVariation;
  // indicates, whether the first shape may be partially occluded or not
  const bool _occlusion;

public:
  LevelsetGeodesicGradient( const typename ConfiguratorType::InitType &Grid,
                            const int NumberOfDeformations,
                            const RealType Gamma,
                            const RealType Mu,
                            const RealType Nu,
                            const RealType Lambda,
                            const RealType NormRegularization,
                            const RealType Epsilon,
                            const RealType LengthEnergyWeight,
                            const RealType VolumeEnergyWeight,
                            const aol::MultiVector<RealType> &Weight,
                            const bool Occlusion = false ) :
    _grid( Grid ),
    _numPolygonLines( NumberOfDeformations ),
    _numShapes( _numPolygonLines + 1 ),
    _gamma( Gamma ),
    _mu( Mu ),
    _nu( Nu ),
    _lambda( Lambda ),
    _heavisideFunction( Epsilon ),
    _elasticEnergyDensity( LengthEnergyWeight, 0, VolumeEnergyWeight ),
    _weight( Weight ),
    _lengthEnergyVariation( _grid, _heavisideFunction, NormRegularization ),
    _volumeOp( _grid, _heavisideFunction ),
    _volumeVariation( _grid, _heavisideFunction ),
    _occlusion( Occlusion ) {}

  // Arg has (2*dim+1)*n+dim+1 components, ie. n+1 levelsets, dim*n displacement components, dim*(n+1) landmark component vectors
  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {

    // compute the volumes of the start and end shapes
    RealType volumeOfEndShapes[NUM_LEVELSETS][2];
    if ( _nu > 0. )
      for ( int i = 0; i < NUM_LEVELSETS; i++ ) {
        aol::Scalar<RealType> volumeOfEndShape;
        _volumeOp.apply( Arg[i*_numShapes+0], volumeOfEndShape );
        volumeOfEndShapes[i][0] = volumeOfEndShape;
        _volumeOp.apply( Arg[(i+1)*_numShapes-1], volumeOfEndShape );
        volumeOfEndShapes[i][1] = volumeOfEndShape;
      }
    else
      for ( int i = 0; i < NUM_LEVELSETS; i++ )
        for ( int j = 0; j < 2; j++ )
          volumeOfEndShapes[i][j] = aol::NumberTrait<RealType>::NaN;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numPolygonLines; i++ ) {
      aol::MultiVector<RealType> sdfLeft, sdfRight, sdfGradient( NUM_LEVELSETS, Arg[0].size() );
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        sdfLeft.appendReference( Arg[j*_numShapes+i], false );
        sdfRight.appendReference( Arg[j*_numShapes+i+1], false );
      }
      aol::MultiVector<RealType> displacement, displacementGradient, landmarksLeft, landmarksRight, landmarksGradientLeft, landmarksGradientRight;
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        displacement.appendReference( Arg[NUM_LEVELSETS*_numShapes+i*ConfiguratorType::Dim+j], false );
        displacementGradient.appendReference( Dest[NUM_LEVELSETS*_numShapes+i*ConfiguratorType::Dim+j], false );
#ifdef USE_LANDMARKS
        landmarksLeft.appendReference( Arg[NUM_LEVELSETS*_numShapes+_numPolygonLines*ConfiguratorType::Dim+i*ConfiguratorType::Dim+j], false );
        landmarksRight.appendReference( Arg[NUM_LEVELSETS*_numShapes+_numPolygonLines*ConfiguratorType::Dim+(i+1)*ConfiguratorType::Dim+j], false );
        landmarksGradientLeft.appendReference( Dest[NUM_LEVELSETS*_numShapes+_numPolygonLines*ConfiguratorType::Dim+i*ConfiguratorType::Dim+j], false );
        landmarksGradientRight.appendReference( Dest[NUM_LEVELSETS*_numShapes+_numPolygonLines*ConfiguratorType::Dim+(i+1)*ConfiguratorType::Dim+j], false );
#endif
      }

#ifdef USE_LANDMARKS
      // derivatives of landmark energy
      aol::MultiVector<RealType> landmarksGrad( landmarksGradientLeft, aol::STRUCT_COPY ), dispLand, dispLandGrad;
      dispLand.appendReference( displacement );
      dispLand.appendReference( landmarksLeft );
      dispLand.appendReference( landmarksRight );
      dispLandGrad.appendReference( displacementGradient );
      dispLandGrad.appendReference( landmarksGrad );
      dispLandGrad.appendReference( landmarksGradientRight );
      LandmarkMismatchGradient<ConfiguratorType>( _grid, _lambda ).applyAdd( dispLand, dispLandGrad );
#endif

      // signed distance function right term
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        LevelsetMismatchGradientWRTDeformedSDF<ConfiguratorType,HeavisideFunctionType>( _grid, sdfLeft[j], displacement, _heavisideFunction, _occlusion && i == 0 ).apply( sdfRight[j], sdfGradient[j] );
        Dest[j*_numShapes+i+1].addMultiple( sdfGradient[j], _gamma );
      }

      // signed distance function left term
      for ( int j = 0; j < NUM_LEVELSETS; j++ )
        LevelsetMismatchGradientWRTSDF<ConfiguratorType,HeavisideFunctionType>( _grid, sdfRight[j], displacement, _heavisideFunction, _occlusion && i == 0 ).apply( sdfLeft[j], sdfGradient[j] );
      sdfGradient *= _gamma;

      // displacement term
      aol::MultiVector<RealType> aux( displacement, aol::STRUCT_COPY );
      WeightedEnergyDensityType weightedDensity( _elasticEnergyDensity, _grid, _weight[i] );
      qc::HyperelasticGradient<ConfiguratorType,WeightedEnergyDensityType>( _grid, weightedDensity ).applyAdd( displacement, displacementGradient );
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        LevelsetMismatchGradientWRTDisp<ConfiguratorType,HeavisideFunctionType>( _grid, sdfLeft[j], sdfRight[j], _heavisideFunction, _occlusion && i == 0 ).apply( displacement, aux );
        displacementGradient.addMultiple( aux, _gamma );
      }

      // length regularization energy of signed distance function left
      if ( _mu > 0. && i > 0 ) {
        aol::MultiVector<RealType> dest( sdfLeft, aol::STRUCT_COPY );
        _lengthEnergyVariation.applyAdd( sdfLeft, dest );
        sdfGradient.addMultiple( dest, _mu );
      }

      // volume deviation energy of signed distance function left
      if ( _nu > 0. && i > 0 )
        for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
          aol::Vector<RealType> dest( sdfLeft[j], aol::STRUCT_COPY );
          aol::Scalar<RealType> volume;
          _volumeOp.apply( sdfLeft[j], volume );
          _volumeVariation.apply( sdfLeft[j], dest );
          dest *= 2 * ( volume.v - ( ( ( volumeOfEndShapes[j][1] - volumeOfEndShapes[j][0] ) * i ) / _numPolygonLines + volumeOfEndShapes[j][0] ) );
          sdfGradient[j].addMultiple( dest, _nu );
        }

#ifdef _OPENMP
//#pragma omp barrier // prevent collision with above writing in Dest; no performance loss since it is last statement of loop
#endif
      for ( int j = 0; j < NUM_LEVELSETS; j++ )
        Dest[j*_numShapes+i] += sdfGradient[j];
#ifdef USE_LANDMARKS
      landmarksGradientLeft += landmarksGrad;
#endif
    }
    for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
      Dest[j*_numShapes].setZero();
      Dest[j*_numShapes+_numPolygonLines].setZero();
    }
#ifdef USE_LANDMARKS
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      Dest[NUM_LEVELSETS*_numShapes+_numPolygonLines*ConfiguratorType::Dim+i].setZero();
      Dest[NUM_LEVELSETS*_numShapes+2*_numPolygonLines*ConfiguratorType::Dim+i].setZero();
    }
#endif
  }
};

template <typename ConfiguratorType>
class MeanCurvatureOp :
  public aol::FENonlinDiffOpInterface<ConfiguratorType,MeanCurvatureOp<ConfiguratorType> > {
  // for mean curvature on a signed distance function, a stiffness operator would be sufficient
  // assumes negative values inside and positive outside the shape

private:
  typedef typename ConfiguratorType::RealType RealType;
  const RealType _epsSqr;
  typename ConfiguratorType::MatrixType _massMatrix;

public:
  MeanCurvatureOp( const typename ConfiguratorType::InitType &Grid, const RealType Epsilon ) :
    aol::FENonlinDiffOpInterface<ConfiguratorType,MeanCurvatureOp<ConfiguratorType> >( Grid ),
    _epsSqr( aol::Sqr( Epsilon ) ),
    _massMatrix( Grid ) {
    aol::MassOp<ConfiguratorType>( Grid ).assembleAddMatrix( _massMatrix );
  }

  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         typename ConfiguratorType::VecType &NL ) const {
    DiscFunc.evaluateGradientAtQuadPoint( El, QuadPoint, NL );
    NL *= -1. / ( _epsSqr + NL.normSqr() );
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // compute massMatrix * meanCurvature
    aol::Vector<RealType> dest( Dest, aol::STRUCT_COPY );
    aol::FENonlinDiffOpInterface<ConfiguratorType,MeanCurvatureOp<ConfiguratorType> >::applyAdd( Arg, dest );
    // compute meanCurvature
    aol::DiagonalPreconditioner< aol::Vector<RealType> > preCond( _massMatrix );
    aol::PCGInverse<aol::Vector<RealType> > inv( _massMatrix, preCond, 1e-16, 1000, aol::STOPPING_ABSOLUTE );
    inv.setQuietMode( true );
    inv.apply( dest, Dest );
  }
};

template <typename ConfiguratorType>
class FeatureEnhancingOp :
  public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  // local erosion and dilation via an epsilon-ball, whose radius depends on local curvature
  const RealType _epsFactor, _maxEpsilon;
  const MeanCurvatureOp<ConfiguratorType> _meanCurvatureOp;

public:
  FeatureEnhancingOp( const typename ConfiguratorType::InitType &Grid,
                      const RealType EpsFactor,
                      const RealType MaxEpsilon ) :
    _grid( Grid ),
    _epsFactor( EpsFactor ),
    _maxEpsilon( MaxEpsilon ),
    _meanCurvatureOp( Grid, 1.e-4 ) {}

  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // compute local curvature
    aol::Vector<RealType> curvature( Arg, aol::STRUCT_COPY );
    _meanCurvatureOp.apply( Arg, curvature );

    // blur local curvature
    qc::LinearSmoothOp<RealType> linSmooth;
    linSmooth.setCurrentGrid( _grid );
    linSmooth.setSigma( sqrt( 1.2 * _grid.H() ) );
    linSmooth.apply( curvature, curvature );

    // locally dilate or erose (assumes negative levelset function values inside and positive outside the shape)
    typename ConfiguratorType::ArrayType dilated( Arg, _grid, aol::DEEP_COPY ), erosed( Arg, _grid, aol::DEEP_COPY );
    dilated.addToAll( -20 * _grid.H() );
    erosed.addToAll( 20 * _grid.H() );
    /* alternatively:
    // Dilation and Erosion assume the signs of the levelset function in the opposite way, hence they are exchanged...
    typename ConfiguratorType::ArrayType arg( Arg, _grid, aol::FLAT_COPY );
    morph::Dilation<RealType,ConfiguratorType::Dim> erosion( _grid, 1 );
    morph::Erosion<RealType,ConfiguratorType::Dim> dilation( _grid, 1 );
    dilation.setRadius( 40 * _grid.H() );
    dilation.setTau( 40 * _grid.H() );
    erosion.setRadius( 40 * _grid.H() );
    erosion.setTau( 40 * _grid.H() );
    dilation.apply( arg, dilated );
    erosion.apply( arg, erosed );*/
    for ( int i = 0; i < Arg.size(); i++ ) {
      RealType locCurv = aol::Max( -_maxEpsilon, aol::Min( _maxEpsilon, curvature[i] * _epsFactor ) );
      Dest[i] = dilated[i] + ( _maxEpsilon - locCurv ) / ( 2. * _maxEpsilon ) * ( erosed[i] - dilated[i] );
    }
  }
};


/**************************************************************
 * Classes for averaging of given shapes and related problems.
 **************************************************************/

// code for linear morphing; reads in displacements and objects from a geodesic and linearly morphs between those
/*
int maxLev = 3;
int interpLev = 3;
typename ConfiguratorType::InitType init( _finestLevel, ConfiguratorType::Dim );
aol::RandomAccessContainer<aol::MultiVector<RealType> > displacement;
aol::RandomAccessContainer<ArrayType> image;
for ( int i = 0; i < 1<<maxLev; i++ ) {
  aol::MultiVector<RealType> disp( init );
  for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
    char filename[1024];
    sprintf( filename, "geodesicBending/geodesic2/displacement%d%d.pgm", i+1, j );
    ArrayType( disp[j], init, aol::FLAT_COPY ).load( filename );
  }
  displacement.pushBack( disp );
  char filename[1024];
  sprintf( filename, "geodesicBending/geodesic2/object%d.pgm", i );
  ArrayType img( filename );
  img *= 256 / init.H();
  sprintf( filename, "imageInit%d.png", i );
  img.savePNG( filename );
  image.pushBack( img );
}

for ( int level = 0; level <= maxLev; level++ ){
  for ( int i = 0; i < (1<<maxLev); i += 1<<(maxLev-level) ){
    aol::MultiVector<RealType> disp( init ), aux( init );
    for ( int j = i; j < i + (1<<(maxLev-level)); j++ ){
      qc::ConcatAndSmoothlyExtendDeformations<ConfiguratorType>( displacement[j], disp, init, aux );
      disp = aux;
    }
    for ( int j = 0; j < (1<<(interpLev-level)); j++ ) {
      aol::MultiVector<RealType> displ( disp, aol::DEEP_COPY );
      displ *= ( j * 1. ) / (1<<(interpLev-level));
      ArrayType defImg( init );
      qc::InvDeformImage<ConfiguratorType,aol::Vector<RealType> >( image[i], init, defImg, displ );
      //qc::DeformAndSmoothlyExtendImage<ConfiguratorType,aol::Vector<RealType> >( image[i], init, defImg, displ );
      char filename[1024];
      sprintf( filename, "image%d%d.png", level, i * (1<<(interpLev-maxLev)) + j );
      defImg.savePNG( filename );
    }
  }
}
abort();
*/

// code for color mapping
/*int size = 15;
RealType data[15] = {0.0378,0.1213,0.2587,0.2465,0.3767,0.1007,0.2554,0.2438,0.3668,0.2352,0.2363,0.3232,0.1362,0.3123,0.4037};
aol::Vector<RealType> values( data, size );
qc::MultiArray<unsigned char,2,3> colors( size, 1 );
convertToColor( values, colors, values.getMinValue(), values.getMaxValue() );
for ( int i = 0; i < size; i++ )
  cerr<<"{\\definecolor{localColor}{rgb}{"<<colors[0][i]/255.<<","<<colors[1][i]/255.<<","<<colors[2][i]/255.<<"}\\multicolumn{1}{>{\\columncolor{localColor}}c}{"<<values[i]<<"}}"<<endl;
abort();*/

template <typename ConfiguratorType>
class ShapeGeodesicOp {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef qc::RestrictOp<RealType,qc::STD_QUOC_RESTRICT/*qc::THROW_AWAY_RESTRICT*/> RestrictOpType;

  // if you want to allow compression to zero volume, choose the second
  //typedef qc::HyperelasticEnergyDensityDefault<ConfiguratorType> MaterialLawType;
  typedef CompressionAllowingEnergyDensity<ConfiguratorType> MaterialLawType;

  // in case of using multilinear FE
  //typedef qc::ProlongOp<RealType> ProlongOpType;
  //typedef ConfiguratorType ConfiguratorTypeMultiLin

  // in case of using simplicial FE
  typedef qc::simplex::ProlongOp<RealType,ConfiguratorType::Dim> ProlongOpType;
  typedef qc::QuocConfiguratorTraitMultiLin<RealType, ConfiguratorType::Dim, aol::GaussQuadrature<RealType, ConfiguratorType::Dim, 3> > ConfiguratorTypeMultiLin;


  // choose the identity function, if geodesic is computed between greyscale images instead of shapes
  typedef aol::ArcTanHeavisideFunction<RealType> HeavisideFunctionType;
  //typedef aol::IdentityFunction<RealType> HeavisideFunctionType;

  // coarsest and finest level for the multiscale method as well as time step
  const int _coarsestLevel, _finestLevel;
  // coarse to fine grids over the examined domain
  std::vector<typename ConfiguratorType::InitType*> _grids;
  // along the geodesic, we have _numberOfShapes shapes
  const int _numberOfShapes;
  const RealType _tau;
  // the objects (as signed distance functions) along the geodesic
  std::vector<qc::MultilevelArray<RealType,qc::Array<RealType>,ProlongOpType,RestrictOpType>*> _signedDistFuncs[NUM_LEVELSETS];
#ifdef USE_LANDMARKS
  // landmarks
  std::vector<aol::MultiVector<RealType>*> _landmarks;
#endif
  // the directory for the results
  char _destDirectory[1024];
  // weights of the hyperelastic deformation energy
  const RealType _lengthEnergyWeight, _volumeEnergyWeight;
  const bool _constantElasticParameters;
  aol::MultiVector<RealType> _weights;
  const bool _occlusion, _blackWhiteInputOutput, _inputRedistancing;
  // levelset parameters
  const RealType _epsFactor, _gamma, _mu, _nu, _lambda, _normReg;
  aol::Vector<RealType> _delta, _deltaInit;
  // iteration number of alternate minimisations and gradient descent steps for the deformations per minimisation
  const int _maxIterations, _maxSteps;
  // displacements
  std::vector<qc::MultiDimMultilevelArray<RealType,qc::Array<RealType>,ProlongOpType,RestrictOpType>*> _displacements;

  // energy values for testing
  aol::Vector<RealType> _totalEnergy;
  aol::RandomAccessContainer<aol::MultiVector<RealType> > _energyComponents; // fidelity, elastic, perimeter, volume, landmark energy

  // Pow = 1 means final, Pow = 0 initial delta
  inline void updateWeights( const int Level, const RealType Pow ){
    _weights.reallocate( _numberOfShapes, _grids[Level]->getNumberOfNodes() );
    if ( _constantElasticParameters )
      _weights.setAll( 1. );
    else {
      // compute weight values to be used
      aol::Vector<RealType> weight( 1 << NUM_LEVELSETS );
      weight[0] = 1.;
      for ( int j = 0; j < _delta.size(); j++ )
        weight[j+1] = _deltaInit[j] * pow( _delta[j] / _deltaInit[j], Pow );
      // update the material properties for each shape
      for ( int i = 0; i < _numberOfShapes; i++ )
        for ( int k = 0; k < _weights[i].size(); k++ ) {
          int phase = 0;
          for ( int j = 0; j < NUM_LEVELSETS; j++ )
            if ( _signedDistFuncs[j][i]->operator[]( Level )[k] < 0 )
              phase += 1 << j;
          _weights[i][k] = weight[phase];
        }
    }
  }

public:
  ShapeGeodesicOp( aol::ParameterParser &Parser ) :
    // prepare the grids
    _coarsestLevel( Parser.getInt( "coarsestLevel" ) ),
    _finestLevel( Parser.getInt( "GridDepth" ) ),
    _grids( _finestLevel + 1 ),
    // load number of objects and create the landmark-vector
    _numberOfShapes( Parser.getInt( "numberOfPoints" ) ),
    _tau( 1. / ( _numberOfShapes - 1 ) ),
#ifdef USE_LANDMARKS
    _landmarks( _numberOfShapes ),
#endif
    // load hyperelastic and levelset parameters
    _lengthEnergyWeight( Parser.getDouble( "lengthEnergyWeight" ) ),
    _volumeEnergyWeight( Parser.getDouble( "volumeEnergyWeight" ) ),
    _constantElasticParameters( Parser.getInt( "constantElasticParameters" ) ),
    _occlusion( Parser.getInt( "occlusion" ) == 1 ),
    _blackWhiteInputOutput( Parser.getInt( "blackWhiteInputOutput" ) == 1 ),
    _inputRedistancing( Parser.getInt( "inputRedistancing" ) == 1 ),
    _epsFactor( Parser.getDouble( "epsFactor" ) ),
    _gamma( Parser.getDouble( "gamma" ) ),
    _mu( Parser.getDouble( "mu" ) ),
    _nu( Parser.getDouble( "nu" ) ),
    _lambda( Parser.hasVariable( "lambda" ) ? Parser.getDouble( "lambda" ) : 0. ),
    _normReg( Parser.getDouble( "normReg" ) ),
    // load number of optimisation steps per level and gradient descent steps per optimisation step
    _maxIterations( Parser.getInt( "maxIterations" ) ),
    _maxSteps( Parser.getInt( "maxSteps" ) ),
    // create vector of displacements
    _displacements( _numberOfShapes ),
    // energy value vectors for testing
    _totalEnergy( _maxIterations + 1 ),
    _energyComponents( _maxIterations + 1 ) {
    for ( int j = 0; j <= _maxIterations; j++ )
      _energyComponents[j].reallocate( 3 * NUM_LEVELSETS + 2, _numberOfShapes - 1 );
    Parser.getRealVec( "delta", _delta );
    Parser.getRealVec( "deltaInitial", _deltaInit );

    // read in the directory, where results are to be saved
    Parser.getString( "destDirectory", _destDirectory );

    // create the grids
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int level = 0; level <= _finestLevel; level++ )
      _grids[level] = new typename ConfiguratorType::InitType( qc::GridSize<ConfiguratorType::Dim>( static_cast<qc::CoordType::DataType>( ( 1 << level ) + 1 ) ) );

    // create the objects and displacements
    for ( int j = 0; j < NUM_LEVELSETS; j++ )
      _signedDistFuncs[j].resize( _numberOfShapes );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      for ( int j = 0; j < NUM_LEVELSETS; j++ )
        _signedDistFuncs[j][i] = new qc::MultilevelArray<RealType,qc::Array<RealType>,ProlongOpType,RestrictOpType>( _finestLevel, ConfiguratorType::Dim );
      _displacements[i] = new qc::MultiDimMultilevelArray<RealType,qc::Array<RealType>,ProlongOpType,RestrictOpType>( _finestLevel, ConfiguratorType::Dim, ConfiguratorType::Dim );
    }

    // load the input objects and make them levelset functions to the mean value level set
    for ( int i = 0; i < _numberOfShapes; i += _numberOfShapes - 1 )
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        char variable[1024], fileName[1024];
        sprintf( variable, ( i == 0 ) ? "startImageFile%d" : "endImageFile%d", j );
        strncpy ( fileName, Parser.getString( variable ).c_str(), 1024 );
        cerr<<"load " + string( fileName ) + "...\n";
        ArrayType tmpArray( fileName );
        if ( _blackWhiteInputOutput )
          tmpArray.addToAll( -.5 * ( tmpArray.getMinValue() + tmpArray.getMaxValue() ) );
        _signedDistFuncs[j][i]->current() = tmpArray;
      }

#ifdef USE_LANDMARKS
    // load the landmark positions
    _landmarks[0] = new aol::MultiVector<RealType>( ConfiguratorType::Dim, 0 );
    _landmarks[_numberOfShapes-1] = new aol::MultiVector<RealType>( ConfiguratorType::Dim, 0 );
    if ( Parser.hasVariable( "landmarksStart" ) ) {
      Parser.getRealMultiVec( "landmarksStart", *_landmarks[0] );
      Parser.getRealMultiVec( "landmarksEnd", *_landmarks[_numberOfShapes-1] );
    }
    for ( int i = 1; i < _numberOfShapes - 1; i++ )
      _landmarks[i] = new aol::MultiVector<RealType>( *_landmarks[0], aol::DEEP_COPY );
#endif

    // initialization with earlier results
    const typename ConfiguratorType::InitType &grid( *_grids[_finestLevel] );
    if ( Parser.hasVariable( "InitDeformFileNameTemplate" ) ) {
      char deformFileTempl[1024], levelFileTempl[1024];
      Parser.getString( "InitDeformFileNameTemplate", deformFileTempl );
      Parser.getString( "InitLevelsetFileNameTemplate", levelFileTempl );
      int kNew = _numberOfShapes - 1, kOld = Parser.getInt( "InitNumberOfShapes" ) - 1;

      // load earlier displacements and level set functions
      aol::RandomAccessContainer<aol::MultiVector<RealType> > displacementOld( kOld + 1 ),
                                                              displacement( kNew + 1 );
      aol::RandomAccessContainer<aol::MultiVector<RealType> > levelsetFuncOld( NUM_LEVELSETS );
      for ( int j = 0; j < NUM_LEVELSETS; j++ )
        levelsetFuncOld[j].reallocate( kOld + 1, _signedDistFuncs[0][0]->current().size() );
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 1; i <= kOld; i++ ) {
        displacementOld[i].reallocate( ConfiguratorType::Dim, grid.getNumberOfNodes() );
        char filename[1024];
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          sprintf( filename, deformFileTempl, i, j );
          cerr<<"load " + string( filename ) + "...\n";
          ArrayType( displacementOld[i][j], grid, aol::FLAT_COPY ).load( filename );
        }
        for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
          sprintf( filename, levelFileTempl, j, i - 1 );
          cerr<<"load " + string( filename ) + "...\n";
          ArrayType( levelsetFuncOld[j][i-1], grid, aol::FLAT_COPY ).load( filename );
        }
      }
      if ( ConfiguratorType::Dim == qc::QC_2D ) // black-white images instead signed distance functions
        for ( int j = 0; j < NUM_LEVELSETS; j++ )
          levelsetFuncOld[j].addToAll( - .5 * ( levelsetFuncOld[j].getMaxValue() - levelsetFuncOld[j].getMinValue() ) );
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 1; i <= kNew; i++ )
        _displacements[i]->appendReferencesTo( displacement[i] );

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 1; i <= kNew; i++ ){
        int i1 = ( ( i * kOld - 1 ) / kNew ) + 1,     // ceil( i*kOld/kNew )
            i2 = ( ( ( i - 1 ) * kOld ) / kNew ) + 1; // floor( (i-1)*kOld/kNew )+1
        aol::MultiVector<RealType> displ1( displacementOld[i1], aol::DEEP_COPY ),
                                   displ2( displacementOld[i2], aol::DEEP_COPY );
        // compute displacement
        if ( i1 == i2 ){
          displ1 *= ( 1. * kOld ) / kNew;
          displ2 *= ( i - 1 ) * ( 1. * kOld ) / kNew - i1 + 1;
          qc::InvDeformAndSmoothlyExtendImage<ConfiguratorTypeMultiLin,aol::MultiVector<RealType> >( displ1, cubicGrid( grid ), displacement[i], displ2 );
        } else {
          aol::MultiVector<RealType> displ2Copy( displacementOld[i2], aol::DEEP_COPY ),
                                     displAux( displacementOld[i2], aol::STRUCT_COPY );
          displ2 *= i2 - ( i - 1 ) * ( 1. * kOld ) / kNew;
          displ2Copy *= 1 - i2 + ( i - 1 ) * ( 1. * kOld ) / kNew;
          qc::InvDeformAndSmoothlyExtendImage<ConfiguratorTypeMultiLin,aol::MultiVector<RealType> >( displ2, cubicGrid( grid ), displAux, displ2Copy );
          for ( int j = i2 + 1; j < i1; j++ ) {
            qc::InvDeformAndSmoothlyExtendImage<ConfiguratorTypeMultiLin,aol::MultiVector<RealType> >( displacementOld[j], cubicGrid( grid ), displ2Copy, displAux );
            displAux += displ2Copy;
          }
          displ1 *= i * ( 1. * kOld ) / kNew - i1 + 1;
          qc::InvDeformAndSmoothlyExtendImage<ConfiguratorTypeMultiLin,aol::MultiVector<RealType> >( displ1, cubicGrid( grid ), displacement[i], displAux );
          displacement[i] += displAux;
        }
        // compute levelset function
        displ1 = displacementOld[i1];
        displ1 *= 1 - i1 + i * ( 1. * kOld ) / kNew;
        if ( i < kNew )
          for ( int j = 0; j < NUM_LEVELSETS; j++ )
            qc::InvDeformAndSmoothlyExtendImage<ConfiguratorTypeMultiLin,aol::Vector<RealType> >( levelsetFuncOld[j][i1-1], cubicGrid( grid ), _signedDistFuncs[j][i]->current(), displ1 );
      }
    } else
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        ArrayType aux( grid );
        char variable[1024];
        sprintf( variable, "InitLevelsetFileName%d", j );
        if ( Parser.hasVariable( variable ) ) {
          aux.load( Parser.getString( variable ).c_str() );
          if ( _blackWhiteInputOutput )
            aux.addToAll( -.5 * ( aux.getMinValue() + aux.getMaxValue() ) );
        } else {
          typename ConfiguratorType::VecType center( .5 );
          qc::DataGenerator<ConfiguratorTypeMultiLin>( cubicGrid( grid ) ).generateSphereLevelset( aux, center, 0.35 );
          aux *= -1.;
        }
        for ( int i = 1; i < _numberOfShapes - 1; i++ )
          _signedDistFuncs[j][i]->current() = aux;
      }

    // generate the needed coarser scales and make all levelset functions signed distance functions of the zero levelset
#ifdef _OPENMP
//#pragma omp parallel for // SignedDistOp not thread-safe?!
#endif
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      if ( _coarsestLevel != _finestLevel ) {
        for ( int j = 0; j < NUM_LEVELSETS; j++ )
          _signedDistFuncs[j][i]->levRestrict( _coarsestLevel, _finestLevel );
        _displacements[i]->levRestrict( _coarsestLevel, _finestLevel );
      }
      if ( _inputRedistancing )
        for ( int level = _coarsestLevel; level <= _finestLevel; level++ ) // on all levels to prevent plateaus from restriction
          for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
            aol::Vector<RealType> aux( _signedDistFuncs[j][i]->operator[]( level ), aol::DEEP_COPY );
            ( typename qc::SignedDistanceOpTrait<ConfiguratorTypeMultiLin,ConfiguratorType::Dim>::OpType( cubicGrid( *_grids[level] ) ) ).apply( aux, _signedDistFuncs[j][i]->operator[]( level ) );
          }
    }

    // notify the user
    cerr << "Parameters loaded." << endl;
  }

  ~ShapeGeodesicOp() {
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      for ( int j = 0; j < NUM_LEVELSETS; j++ )
        delete _signedDistFuncs[j][i];
      delete _displacements[i];
#ifdef USE_LANDMARKS
      delete _landmarks[i];
#endif
    }
  }

  inline void putGeodesicIntoMultiVector( const int Level, aol::MultiVector<RealType> &Geodesic ){
    for ( int j = 0; j < NUM_LEVELSETS; j++ )
      for ( int i = 0; i < _numberOfShapes; i++ )
        Geodesic.appendReference( _signedDistFuncs[j][i]->operator[]( Level ), false );
    for ( int i = 0; i < _numberOfShapes - 1; i++ ) {
      _displacements[i+1]->setCurLevel( Level );
      _displacements[i+1]->appendReferencesTo( Geodesic );
    }
#ifdef USE_LANDMARKS
    for ( int i = 0; i < _numberOfShapes; i++ )
      Geodesic.appendReference( *_landmarks[i], false );
#endif
  }

  inline void debuggingOutput( const int Level, const int Index ) {
    typename ConfiguratorType::InitType &grid( *_grids[ Level ] );
    RealType epsilon = grid.H() * _epsFactor;

    // append all level set functions and displacements and landmarks to one MultiVector
    aol::MultiVector<RealType> geodesic;
    putGeodesicIntoMultiVector( Level, geodesic );

    // compute the different energy component values
    LevelsetGeodesicEnergy<ConfiguratorType,HeavisideFunctionType,MaterialLawType> energy( grid, _numberOfShapes - 1, _gamma, _mu, _nu, _lambda, _normReg, epsilon, _lengthEnergyWeight, _volumeEnergyWeight, _weights, _occlusion );
    energy.computeEnergyComponents( geodesic, _energyComponents[Index] );
    _totalEnergy[Index] = _energyComponents[Index][NUM_LEVELSETS].sum() + _lambda * _energyComponents[Index][3*NUM_LEVELSETS+1].sum();
    for ( int j = 0; j < NUM_LEVELSETS; j++ )
      _totalEnergy[Index] += _gamma * _energyComponents[Index][j].sum() + _mu * _energyComponents[Index][NUM_LEVELSETS+j+1].sum() + _nu * _energyComponents[Index][2*NUM_LEVELSETS+j+1].sum();

#ifdef DEBUGGING_OUTPUT
    // compute the different energy component densities
    aol::MultiVector<RealType> fidelityEnergy, elasticEnergy, perimeterEnergy, volumeEnergy;
    energy.computeEnergyDensities( geodesic, fidelityEnergy, elasticEnergy, perimeterEnergy, volumeEnergy );
    for ( int i = 0; i < _numberOfShapes - 1; i++ ) {
      saveAsColorImage<ConfiguratorType>( elasticEnergy[i], grid, _destDirectory, "elasticEnergy", Level * 10000 + Index * 100 + i );
      char fileName[1024];
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        sprintf( fileName, "mismatch%dEnergy", j );
        saveAsColorImage<ConfiguratorType>( fidelityEnergy[j*(_numberOfShapes-1)+i], grid, _destDirectory, fileName, Level * 10000 + Index * 100 + i );
      }
    }

    // compute the energy variations
    aol::MultiVector<RealType> geodesicGrad( geodesic, aol::STRUCT_COPY );
    LevelsetGeodesicGradient<ConfiguratorType,HeavisideFunctionType,MaterialLawType>( grid, _numberOfShapes - 1, _gamma, _mu, _nu, _lambda, _normReg, epsilon, _lengthEnergyWeight, _volumeEnergyWeight, _weights, _occlusion ).apply( geodesic, geodesicGrad );

    for ( int i = 0; i < _numberOfShapes - 1; i++ ) {
      // save variation wrt level set function
      char fileName[1024];
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        sprintf( fileName, "levelset%dGrad", j );
        saveAsColorImage<ConfiguratorType>( geodesicGrad[j*_numberOfShapes+i], grid, _destDirectory, fileName, Level * 10000 + Index * 100 + i, -.0004, .0004 );
      }

      // save variation wrt displacement
      aol::MultiVector<RealType> dispGrad;
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        dispGrad.appendReference( geodesicGrad[NUM_LEVELSETS*_numberOfShapes+i*ConfiguratorType::Dim+j], false );
      dispGrad *= 500;
      saveDisplacementAsChessPattern<ConfiguratorTypeMultiLin,PNG>( dispGrad, cubicGrid( grid ), _destDirectory, "displGrad", Level * 10000 + Index * 100 + i );
      //saveAsColorImage<ConfiguratorType>( dispGrad[0], grid, _destDirectory, "dispxGrad", Level * 10000 + Index * 100 + i, -.05, .05 );
      //saveAsColorImage<ConfiguratorType>( dispGrad[1], grid, _destDirectory, "dispyGrad", Level * 10000 + Index * 100 + i, -.05, .05 );
    }
#endif

    saveEnergyValues( Level, Index );
  }

  inline void saveEnergyValues( const int Level, const int /*Index*/ ) {
    char fileName[1024];
    sprintf( fileName, "%s/energyValuesOnLevel%d.txt", _destDirectory, Level );
    ofstream energyFile( fileName );
    energyFile << "elastic energy:" << endl;
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      for ( int j = 0; j < _totalEnergy.size(); j++ )
        energyFile << _energyComponents[j][NUM_LEVELSETS][i] << " ";
      energyFile << endl;
    }
    energyFile << endl << "mismatch energy" << endl;
    for ( int k = 0; k < NUM_LEVELSETS; k++ )
      for ( int i = 0; i < _numberOfShapes; i++ ) {
        for ( int j = 0; j < _totalEnergy.size(); j++ )
          energyFile << _energyComponents[j][k][i] << " ";
        energyFile << endl;
      }
    energyFile << endl << "total energy" << endl;
    for ( int i = 0; i < _totalEnergy.size(); i++ )
      energyFile << _totalEnergy[i] << " ";
    energyFile << endl;
    energyFile << endl << "perimeter energy:" << endl;
    for ( int k = 0; k < NUM_LEVELSETS; k++ )
      for ( int i = 0; i < _numberOfShapes; i++ ) {
        for ( int j = 0; j < _totalEnergy.size(); j++ )
          energyFile << _energyComponents[j][NUM_LEVELSETS+1+k][i] << " ";
        energyFile << endl;
      }
    energyFile << endl << "volume energy" << endl;
    for ( int k = 0; k < NUM_LEVELSETS; k++ )
      for ( int i = 0; i < _numberOfShapes; i++ ) {
        for ( int j = 0; j < _totalEnergy.size(); j++ )
          energyFile << _energyComponents[j][2*NUM_LEVELSETS+1+k][i] << " ";
        energyFile << endl;
      }
    energyFile << endl << "landmark energy" << endl;
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      for ( int j = 0; j < _totalEnergy.size(); j++ )
        energyFile << _energyComponents[j][3*NUM_LEVELSETS+1][i] << " ";
      energyFile << endl;
    }
    energyFile.close();
  }

  void saveResults( const int Index, const int Level ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numberOfShapes; i++ ) {
      // save the objects
      char fileName[1024];
      ArrayType image( _signedDistFuncs[0][0]->operator[]( Level ), aol::STRUCT_COPY );
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        sprintf( fileName, "object%d_", j );
        image = _signedDistFuncs[j][i]->operator[]( Level );
        if ( _blackWhiteInputOutput )
          image.clamp( - _grids[Level]->H() / 2, _grids[Level]->H() / 2 );
        //save<ConfiguratorType>( image, *_grids[Level], _destDirectory, fileName, i );
        saveAsImage<ConfiguratorType,PGM>( image, *_grids[Level], _destDirectory, fileName, i );
#ifdef DEBUGGING_OUTPUT
        saveAsImage<ConfiguratorType,PNG>( image, *_grids[Level], _destDirectory, fileName, Level * 10000 + Index * 100 + i );
        //sprintf( fileName, "levelset%dFunc", j );
        //saveAsColorImage<ConfiguratorType>( _signedDistFuncs[j][i]->operator[]( Level ), *_grids[Level], _destDirectory, fileName, Level * 10000 + Index * 100 + i, -.6, .2 );
#endif
      }

      if ( i > 0 ) {
        // save the displacements
        aol::MultiVector<RealType> displacement;
        _displacements[i]->setCurLevel( Level );
        _displacements[i]->appendReferencesTo( displacement );
        saveDisplacement<ConfiguratorType>( displacement, *_grids[Level], _destDirectory, "displacement", i );
#ifdef DEBUGGING_OUTPUT
        saveDisplacementAsChessPattern<ConfiguratorTypeMultiLin,PNG>( displacement, cubicGrid( *_grids[Level] ), _destDirectory, "dispChessPat", Level * 10000 + Index * 100 + i );

        // save the back-deformed objects
        for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
          qc::DeformImage<ConfiguratorTypeMultiLin>( _signedDistFuncs[j][i]->operator[]( Level ), cubicGrid( *_grids[Level] ), image, displacement );
          if ( ConfiguratorType::Dim == qc::QC_2D && _inputRedistancing )
            image.clamp( 0., _grids[Level]->H() );
          char fileName[1024];
          sprintf( fileName, "backDeformed%dImage", j );
          //save<ConfiguratorType>( image, *_grids[Level], _destDirectory, fileName, i );
          saveAsImage<ConfiguratorType,PNG>( image, *_grids[Level], _destDirectory, fileName, Level * 10000 + Index * 100 + i );
        }
#endif
      }
    }

    debuggingOutput( Level, Index );
  }

  void relax( const int Level ) {
    aol::MultiVector<RealType> geodesic;
    putGeodesicIntoMultiVector( Level, geodesic );
    aol::MultiVector<RealType> initialGeodesic( geodesic, aol::DEEP_COPY );

    LevelsetGeodesicEnergy<ConfiguratorType,HeavisideFunctionType,MaterialLawType> e( *_grids[Level], _numberOfShapes - 1, _gamma, _mu, _nu, _lambda, _normReg, _grids[Level]->H() * _epsFactor, _lengthEnergyWeight, _volumeEnergyWeight, _weights, _occlusion );
    LevelsetGeodesicGradient<ConfiguratorType,HeavisideFunctionType,MaterialLawType> de( *_grids[Level], _numberOfShapes - 1, _gamma, _mu, _nu, _lambda, _normReg, _grids[Level]->H() * _epsFactor, _lengthEnergyWeight, _volumeEnergyWeight, _weights, _occlusion );

    // perform a few steps of gradient descent
#ifdef USE_LANDMARKS
    aol::Vector<int> tauBlocks( 3 );
    tauBlocks[2] = _numberOfShapes * ConfiguratorType::Dim;
#else
    aol::Vector<int> tauBlocks( 2 );
#endif
    tauBlocks[0] = NUM_LEVELSETS * _numberOfShapes;
    tauBlocks[1] = ( _numberOfShapes - 1 ) * ConfiguratorType::Dim;
    aol::GradientDescentComponentWiseTimestepControlled<ConfiguratorType> gradientDescent( *_grids[Level], e, de, tauBlocks, _maxSteps/*, aol::ZOTrait<RealType>::one, aol::ZOTrait<RealType>::zero, geodesic.numComponents()*/ );
    // attention! If smoothing is not deactivated, the smoothed descent direction for the level set functions is so bad that the gradient descent fails!
    gradientDescent.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::USE_NONLINEAR_CG | aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::DO_NOT_SMOOTH_DESCENT_DIRECTION );
    gradientDescent.setTauMin( 1.e-12 ); // close to complete compression there might be very strong changes of the energy!
    //gradientDescent.setTimestepController( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::GRADIENT_SIMPLE );

    //aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > gradientDescent( *_grids[Level], e, de, _maxSteps );
    //aol::QuasiNewtonIterationComponentWiseTimestepControlled<RealType> gradientDescent( e, de, tauBlocks, _maxSteps, 1.e-6, 50 );

    gradientDescent.apply( initialGeodesic, geodesic );
  }

  void execute() {
#ifdef _OPENMP
    omp_set_nested( 1 );
#endif
    cerr<<"relax the geodesic..."<<endl;
    //int fac = 1 << ( _finestLevel - _coarsestLevel );
    //_gamma = _gamma / fac;
    for ( int level = _coarsestLevel; level <= _finestLevel; level++ ) {
      cerr << aol::color::green << "Relaxing on level " << level << " of " << _finestLevel << aol::color::reset << endl;
      for ( int iteration = 1; iteration <= _maxIterations; iteration++ ) {
        updateWeights( level, ( level == _coarsestLevel ) ? iteration * 1. / _maxIterations : 1 ); // make material differences smaller to distribute elastic energy better among the shapes in the first iteration
        if ( iteration == 1 )
          saveResults( iteration - 1, level );
        relax( level );
        saveResults( iteration, level );
      }
      if ( level < _finestLevel )
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int i = 0; i < _numberOfShapes; i++ ) {
          if ( ( i > 0 ) && ( i < _numberOfShapes - 1 ) )
            for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
              _signedDistFuncs[j][i]->setCurLevel( level );
              _signedDistFuncs[j][i]->levProlongate();
            }
          _displacements[i]->setCurLevel( level );
          _displacements[i]->levProlongate();
          // remedy local interpenetrations of matter, which were overlooked on coarse grid
          aol::MultiVector<RealType> displacement;
          _displacements[i]->appendReferencesTo( displacement );
          aol::MultiVector<RealType> oldDisplacement( displacement, aol::DEEP_COPY );
          AntiInterpenetrationHyperelasticEnergyDensity<ConfiguratorType> antiInterpenetrationDensity( 1.e-5, 1.e0 ); // threshold has to be larger than or equal to the threshold employed in the used hyperelastic law
          qc::HyperelasticEnergy<ConfiguratorType,AntiInterpenetrationHyperelasticEnergyDensity<ConfiguratorType> > e( *_grids[level+1], antiInterpenetrationDensity );
          qc::HyperelasticGradient<ConfiguratorType,AntiInterpenetrationHyperelasticEnergyDensity<ConfiguratorType> > de( *_grids[level+1], antiInterpenetrationDensity );
          aol::QuasiNewtonIteration<RealType,aol::MultiVector<RealType>,aol::MultiVector<RealType> >( e, de ).apply( oldDisplacement, displacement );
        }
      //_gamma *= 2;
    }
  }
};

template <typename ConfiguratorType>
class ShapeGeodesicMultiGridOp {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef qc::RestrictOp<RealType,qc::STD_QUOC_RESTRICT/*qc::THROW_AWAY_RESTRICT*/> RestrictOpType;

  // if you want to allow compression to zero volume, choose the second
  //typedef qc::HyperelasticEnergyDensityDefault<ConfiguratorType> MaterialLawType;
  typedef CompressionAllowingEnergyDensity<ConfiguratorType> MaterialLawType;

  // in case of using multilinear FE
  //typedef qc::ProlongOp<RealType> ProlongOpType;
  //typedef ConfiguratorType ConfiguratorTypeMultiLin

  // in case of using simplicial FE
  typedef qc::simplex::ProlongOp<RealType,ConfiguratorType::Dim> ProlongOpType;
  typedef qc::QuocConfiguratorTraitMultiLin<RealType, ConfiguratorType::Dim, aol::GaussQuadrature<RealType, ConfiguratorType::Dim, 3> > ConfiguratorTypeMultiLin;

  // choose the identity function, if geodesic is computed between greyscale images instead of shapes
  typedef aol::ArcTanHeavisideFunction<RealType> HeavisideFunctionType;
  //typedef aol::IdentityFunction<RealType> HeavisideFunctionType;

  typedef typename qc::MultiDimMultilevelArray<RealType,qc::Array<RealType>,ProlongOpType,RestrictOpType> DispMultiLevelType;

  // coarsest and finest level L such that space and time are discretized using 2^L+1 points in each direction
  const int _coarsestSpaceLevel, _finestSpaceLevel;
  const int _coarsestTimeLevel, _finestTimeLevel;
  const int _maxNumberOfTimePoints;
  // coarse to fine grids over the examined domain (will be indexed by the space level)
  std::vector<typename ConfiguratorType::InitType*> _grids;
  // the objects (as levelset functions) along the geodesic (indexed by their number on the finest geodesic in time, starting from 0)
  std::vector<qc::MultilevelArray<RealType,qc::Array<RealType>,ProlongOpType,RestrictOpType>*> _signedDistFuncs[NUM_LEVELSETS];
#ifdef USE_LANDMARKS
  // the landmarks (indexed by the number of the shape on the finest geodesic in time, to which they belong)
  std::vector<aol::MultiVector<RealType>*> _landmarks;
#endif
  // the displacements (indexed by the number of the shape on the finest geodesic in time, to which they map)
  std::vector<DispMultiLevelType*> _displacements;
  // the directory for the results
  char _destDirectory[1024];
  // parameters of the hyperelastic deformation energy
  const RealType _defEnergyParam1, _defEnergyParam2;
  const bool _constantElasticParameters;
  const bool _occlusion, _blackWhiteInputOutput, _inputRedistancing;
  // parameters of the fitting energy
  const RealType _epsFactor, _gamma, _mu, _nu, _lambda, _normReg;
  aol::Vector<RealType> _delta, _deltaInit;
  // iteration numbers
  const int _maxTimeVCycles, _maxSpaceVCycles, _maxSmoothingSteps;

  //! Afterwards, "DefEnergyWeights[i]" is one inside the ith object and close to zero outside, but only on the current time and space level.
  //! Pow = 1 means final, Pow = 0 initial weight distribution
  inline void computeDeformationEnergyWeights( const int TimeLevel, const int SpaceLevel, const RealType Pow, aol::MultiVector<RealType> &DefEnergyWeights ){
    DefEnergyWeights.reallocate( ( 1 << TimeLevel ) + 1, _grids[SpaceLevel]->getNumberOfNodes() );
    if ( _constantElasticParameters )
      DefEnergyWeights.setAll( 1. );
    else {
      aol::Vector<RealType> weight( 1 << NUM_LEVELSETS );
      weight[0] = 1.;
      for ( int j = 0; j < _delta.size(); j++ )
        weight[j+1] = _deltaInit[j] * pow( _delta[j] / _deltaInit[j], Pow );
      const int timePointInc = 1 << ( _finestTimeLevel - TimeLevel );
      for ( int i = 0; i < DefEnergyWeights.numComponents(); i++ )
        for ( int k = 0; k < DefEnergyWeights[i].size(); k++ ) {
          int phase = 0;
          for ( int j = 0; j < NUM_LEVELSETS; j++ )
            if ( _signedDistFuncs[j][i*timePointInc]->operator[]( SpaceLevel )[k] < 0 )
              phase += 1 << j;
          DefEnergyWeights[i][k] = weight[phase];
        }
    }
  }

  //! Afterwards, the first 2^"TimeLevel"+1 components of "Geodesic" contain the objects
  //! and the following dimension*2^"TimeLevel" components the displacements.
  inline void putGeodesicIntoMultiVector( const int TimeLevel, const int SpaceLevel, aol::MultiVector<RealType> &Geodesic ){
    Geodesic.reallocate( 0, 0 );
    const int timePointInc = 1 << ( _finestTimeLevel - TimeLevel );
    for ( int j = 0; j < NUM_LEVELSETS; j++ )
      for ( int i = 0; i < _maxNumberOfTimePoints; i += timePointInc )
        Geodesic.appendReference( _signedDistFuncs[j][i]->operator[]( SpaceLevel ), false );
    for ( int i = timePointInc; i < _maxNumberOfTimePoints; i += timePointInc ) {
      _displacements[i]->setCurLevel( SpaceLevel );
      _displacements[i]->appendReferencesTo( Geodesic );
    }
#ifdef USE_LANDMARKS
    for ( int i = 0; i < _maxNumberOfTimePoints; i += timePointInc )
      Geodesic.appendReference( *_landmarks[i], false );
#endif
  }

  //! Returns a list with all nodes at which Displacement produces self-interpenetrations and true if there are none.
  inline bool checkInterpenetrations( const aol::MultiVector<RealType> &Displacement, const typename ConfiguratorType::InitType &Grid, std::vector<qc::CoordType> &CoordinatesOfProblematicNodes, std::vector<qc::Element> &ProblematicElements, const RealType Threshold ){

    CoordinatesOfProblematicNodes.clear();
    ProblematicElements.clear();

    // compute finite differences of the displacement components in x, y, and z direction
    aol::auto_container<ConfiguratorType::Dim,qc::MultiArray<RealType,ConfiguratorType::Dim,ConfiguratorType::Dim> > finiteDifferences;
    qc::MultiArray<RealType,ConfiguratorType::Dim,ConfiguratorType::Dim> displacement( Grid, Displacement, aol::FLAT_COPY );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      qc::MultiArray<RealType,ConfiguratorType::Dim,ConfiguratorType::Dim> disp( displacement, aol::STRUCT_COPY );
      qc::CoordType offset( Grid.getSize() );
      offset[i] -= 1;
      disp.shiftByOffsetFrom( offset, displacement );
      disp -= displacement;
      finiteDifferences.set_copy( i, disp );
    }
    // compute determinants for each element and node
    RealType elementVolume = ConfiguratorType::Dim == qc::QC_2D ? aol::Sqr( Grid.H() ) : aol::Cub( Grid.H() );
    RealType threshold = Threshold * elementVolume;
    qc::BitArray< ConfiguratorType::Dim > problematicNode( qc::GridSize<ConfiguratorType::Dim>::createFrom( Grid ) );
    const short nodesPerElement = 1 << ConfiguratorType::Dim;
    for ( int z = 0; z < Grid.getNumZ() - 1; z++ )
      for ( int y = 0; y < Grid.getNumY() - 1; y++ )
        for ( int x = 0; x < Grid.getNumX() - 1; x++ ) {
          qc::Element elem( x, y, z );
          bool elementProblematic = false;
          for ( int node = 0; node < nodesPerElement; node++ ) {
            qc::CoordType coord( x + node % 2, y + ( node >> 1 ) % 2, z + ( node >> 2 ) % 2 );
            typename ConfiguratorType::MatType dPhi;
            for ( int i = 0; i < ConfiguratorType::Dim; i++ )
              for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
                qc::CoordType diffCoord( coord );
                diffCoord[j] = elem[j];
                dPhi[i][j] = finiteDifferences[j][i].get( diffCoord );
              }
            for ( int i = 0; i < ConfiguratorType::Dim; i++ )
              dPhi[i][i] += Grid.H();
            if ( dPhi.det() <= threshold ) {
              if ( !problematicNode.get( coord ) ) {
                problematicNode.set( coord, true );
                CoordinatesOfProblematicNodes.push_back( coord );
              }
              if ( !elementProblematic ) {
                elementProblematic = true;
                ProblematicElements.push_back( elem );
              }
            }
          }
        }

    return ( CoordinatesOfProblematicNodes.size() == 0 );
  }

  //! Remedies local interpenetrations of matter, which were e.g. overlooked on a coarse resolution
  inline void removeInterpenetrations( aol::MultiVector<RealType> &Displacement, const typename ConfiguratorType::InitType &Grid ){
    const aol::MultiVector<RealType> oldDisplacement( Displacement, aol::DEEP_COPY );
    /* // alternative approach
    AntiInterpenetrationHyperelasticEnergyDensity<ConfiguratorType> antiInterpenetrationDensity( 3.e-5, 2.e0 ); // threshold has to be larger than or equal to the threshold employed in the used hyperelastic law
    qc::HyperelasticEnergy<ConfiguratorType,AntiInterpenetrationHyperelasticEnergyDensity<ConfiguratorType> > e( Grid, antiInterpenetrationDensity );
    qc::HyperelasticGradient<ConfiguratorType,AntiInterpenetrationHyperelasticEnergyDensity<ConfiguratorType> > de( Grid, antiInterpenetrationDensity );
    //aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > minimization( Grid, e, de, 50 );
    //minimization.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::USE_NONLINEAR_CG | aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::DO_NOT_SMOOTH_DESCENT_DIRECTION );
    aol::QuasiNewtonIteration<RealType,aol::MultiVector<RealType>,aol::MultiVector<RealType> > minimization( e, de, 20, 0. );
    minimization.apply( oldDisplacement, Displacement );*/

    RealType elementVolume = ConfiguratorType::Dim == qc::QC_2D ? aol::Sqr( Grid.H() ) : aol::Cub( Grid.H() );
    RealType threshold = 3e-5;
    std::vector<qc::CoordType> coordinatesOfProblematicNodes;
    std::vector<qc::Element> problematicElements;
    checkInterpenetrations( Displacement, Grid, coordinatesOfProblematicNodes, problematicElements, threshold );
    // resolve interpenetrations
    ResolveInterpenetrationEnergy<ConfiguratorType> e( Grid, coordinatesOfProblematicNodes, threshold, 1. / elementVolume );
    ResolveInterpenetrationGrad<ConfiguratorType> de( Grid, coordinatesOfProblematicNodes, threshold, 1. / elementVolume );
    aol::QuasiNewtonIteration<RealType,aol::MultiVector<RealType>,aol::MultiVector<RealType> > minimization( e, de, 1000, 0. );
    minimization.apply( oldDisplacement, Displacement );
  }

public:
  //! Constructor reads in all data from a parameter file.
  ShapeGeodesicMultiGridOp( aol::ParameterParser &Parser ) :
    // prepare space and time discretization
    _coarsestSpaceLevel( Parser.getInt( "coarsestSpaceLevel" ) ),
    _finestSpaceLevel( Parser.getInt( "finestSpaceLevel" ) ),
    _coarsestTimeLevel( Parser.getInt( "coarsestTimeLevel" ) ),
    _finestTimeLevel( Parser.getInt( "finestTimeLevel" ) ),
    _maxNumberOfTimePoints( ( 1 << _finestTimeLevel ) + 1 ),
    _grids( _finestSpaceLevel + 1 ),
    // create vector of landmarks and displacements
#ifdef USE_LANDMARKS
    _landmarks( _maxNumberOfTimePoints ),
#endif
    _displacements( _maxNumberOfTimePoints ),
    // load energy parameters
    _defEnergyParam1( Parser.getDouble( "defEnergyParam1" ) ),
    _defEnergyParam2( Parser.getDouble( "defEnergyParam2" ) ),
    _constantElasticParameters( Parser.getInt( "constantElasticParameters" ) ),
    _occlusion( Parser.getInt( "occlusion" ) == 1 ),
    _blackWhiteInputOutput( Parser.getInt( "blackWhiteInputOutput" ) == 1 ),
    _inputRedistancing( Parser.getInt( "inputRedistancing" ) == 1 ),
    _epsFactor( Parser.getDouble( "epsFactor" ) ),
    _gamma( Parser.getDouble( "gamma" ) ),
    _mu( Parser.getDouble( "mu" ) ),
    _nu( Parser.getDouble( "nu" ) ),
    _lambda( Parser.hasVariable( "lambda" ) ? Parser.getDouble( "lambda" ) : 0. ),
    _normReg( Parser.getDouble( "normReg" ) ),
    // load iteration numbers
    _maxTimeVCycles( Parser.getInt( "maxTimeVCycles" ) ),
    _maxSpaceVCycles( Parser.getInt( "maxSpaceVCycles" ) ),
    _maxSmoothingSteps( Parser.getInt( "maxSmoothingSteps" ) ) {
    // read in the directory, where results are to be saved
    Parser.getString( "destDirectory", _destDirectory );
    Parser.getRealVec( "delta", _delta );
    Parser.getRealVec( "deltaInitial", _deltaInit );


    // create the grids
    cerr << "Create the space discretization...\n";
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int level = 0; level <= _finestSpaceLevel; level++ )
      _grids[level] = new typename ConfiguratorType::InitType( qc::GridSize<ConfiguratorType::Dim>( static_cast<qc::CoordType::DataType>( ( 1 << level ) + 1 ) ) );

    // create the objects and displacements
    cerr << "Create the time discretization...\n";
    for ( int j = 0; j < NUM_LEVELSETS; j++ )
      _signedDistFuncs[j].resize( _maxNumberOfTimePoints );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _maxNumberOfTimePoints; i++ ) {
      for ( int j = 0; j < NUM_LEVELSETS; j++ )
        _signedDistFuncs[j][i] = new qc::MultilevelArray<RealType,qc::Array<RealType>,ProlongOpType,RestrictOpType>( _finestSpaceLevel, ConfiguratorType::Dim );
      _displacements[i] = new qc::MultiDimMultilevelArray<RealType,qc::Array<RealType>,ProlongOpType,RestrictOpType>( _finestSpaceLevel, ConfiguratorType::Dim, ConfiguratorType::Dim );
    }

    // load the input objects (in 2D, shift the mean value to the zero level set)
    cerr << "Load the input objects...\n";
    for ( int i = 0; i < _maxNumberOfTimePoints; i += _maxNumberOfTimePoints - 1 )
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        ArrayType tmpArray( _signedDistFuncs[j][i]->current(), aol::FLAT_COPY );
        char variable[1024], fileName[1024];
        sprintf( variable, ( i == 0 ) ? "startImageFile%d" : "endImageFile%d", j );
        strncpy ( fileName, Parser.getString( variable ).c_str(), 1024 );
        cerr<<"load " + string( fileName ) + "...\n";
        tmpArray.load( fileName );
        if ( _blackWhiteInputOutput ) // black-white images instead of signed distance functions
          tmpArray.addToAll( -.5 * ( tmpArray.getMinValue() + tmpArray.getMaxValue() ) );
      }

#ifdef USE_LANDMARKS
    // load the landmark positions
    _landmarks[0] = new aol::MultiVector<RealType>( ConfiguratorType::Dim, 0 );
    _landmarks[_maxNumberOfTimePoints-1] = new aol::MultiVector<RealType>( ConfiguratorType::Dim, 0 );
    if ( Parser.hasVariable( "landmarksStart" ) ) {
      Parser.getRealMultiVec( "landmarksStart", *_landmarks[0] );
      Parser.getRealMultiVec( "landmarksEnd", *_landmarks[_maxNumberOfTimePoints-1] );
    }
    for ( int i = 1; i < _maxNumberOfTimePoints - 1; i++ )
      _landmarks[i] = new aol::MultiVector<RealType>( *_landmarks[0], aol::DEEP_COPY );
#endif

    // initialization of the objects
    cerr << "Initialize objects on the geodesic...\n";
    ArrayType initObject( *_grids[_finestSpaceLevel] );
    for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
      char variable[1024];
      sprintf( variable, "InitLevelsetFileName%d", j );
      if ( Parser.hasVariable( variable ) ) {
        initObject.load( Parser.getString( variable ).c_str() );
        if ( _blackWhiteInputOutput )
          initObject.addToAll( -.5 * ( initObject.getMinValue() + initObject.getMaxValue() ) );
      } else {
        typename ConfiguratorType::VecType center( .5 );
        qc::DataGenerator<ConfiguratorTypeMultiLin>( cubicGrid( *_grids[_finestSpaceLevel] ) ).generateSphereLevelset( initObject, center, 0.35 );
        initObject *= -1.;
      }
      for ( int i = 1; i < _maxNumberOfTimePoints - 1; i++ )
        _signedDistFuncs[j][i]->current() = initObject;
    }

    // generate the needed coarser scales and make levelset functions signed distance functions
    cerr << "Prepare space coarsening and convert objects into signed distance functions...\n";
#ifdef _OPENMP
//#pragma omp parallel for // SignedDistOp not thread-safe?!
#endif
    for ( int i = 0; i < _maxNumberOfTimePoints; i++ ) {
      if ( _coarsestSpaceLevel != _finestSpaceLevel ) {
        for ( int j = 0; j < NUM_LEVELSETS; j++ )
          _signedDistFuncs[j][i]->levRestrict( _coarsestSpaceLevel, _finestSpaceLevel );
        _displacements[i]->levRestrict( _coarsestSpaceLevel, _finestSpaceLevel );
      }
      if ( _inputRedistancing )
        for ( int level = _coarsestSpaceLevel; level <= _finestSpaceLevel; level++ ) // on all levels to prevent plateaus from restriction
          for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
            aol::Vector<RealType> aux( _signedDistFuncs[j][i]->operator[]( level ), aol::DEEP_COPY );
            ( typename qc::SignedDistanceOpTrait<ConfiguratorTypeMultiLin,ConfiguratorType::Dim>::OpType( cubicGrid( *_grids[level] ) ) ).apply( aux, _signedDistFuncs[j][i]->operator[]( level ) );
          }
    }

    // notify the user
    cerr << "Initialization finished.\n";
  }

  //! Destructor only has to delete objects and displacements
  ~ShapeGeodesicMultiGridOp() {
    for ( int i = 0; i < _maxNumberOfTimePoints; i++ ) {
      for ( int j = 0; j < NUM_LEVELSETS; j++ )
        delete _signedDistFuncs[j][i];
      delete _displacements[i];
#ifdef USE_LANDMARKS
      delete _landmarks[i];
#endif
    }
  }

  //! Computes the energy components and the energy gradients and saves them
  inline void saveEnergiesAndVariations( const int TimeLevel, const int SpaceLevel, const int Index ) {
    const typename ConfiguratorType::InitType &grid( *_grids[SpaceLevel] );
    aol::MultiVector<RealType> geodesic, defEnergyWeights;
    computeDeformationEnergyWeights( TimeLevel, SpaceLevel, 1, defEnergyWeights );
    putGeodesicIntoMultiVector( TimeLevel, SpaceLevel, geodesic );

    LevelsetGeodesicEnergy<ConfiguratorType,HeavisideFunctionType,MaterialLawType> e( grid, 1 << TimeLevel, _gamma, _mu, _nu, _lambda, _normReg, grid.H() * _epsFactor, _defEnergyParam1, _defEnergyParam2, defEnergyWeights, _occlusion );
    LevelsetGeodesicGradient<ConfiguratorType,HeavisideFunctionType,MaterialLawType> de( grid, 1 << TimeLevel, _gamma, _mu, _nu, _lambda, _normReg, grid.H() * _epsFactor, _defEnergyParam1, _defEnergyParam2, defEnergyWeights, _occlusion );

    // save the energies
    aol::MultiVector<RealType> energyComponents( 3 * NUM_LEVELSETS + 1, 1 << TimeLevel );
#ifdef USE_LANDMARKS
    energyComponents.resize( 3 * NUM_LEVELSETS + 2, 1 << TimeLevel );
#endif
    e.computeEnergyComponents( geodesic, energyComponents );
    fstream energyFile( strcat( _destDirectory, "/energyValues.txt" ) );
    energyFile.seekp( 0, ios_base::end );
    energyFile << endl << endl << Index << endl << "deformation energies: " << energyComponents[NUM_LEVELSETS] << endl;
    energyFile << "fidelity energies: ";
    for ( int j = 0; j < NUM_LEVELSETS; j++ )
      energyFile << energyComponents[j];
    energyFile << "length regularization energies: ";
    for ( int j = NUM_LEVELSETS + 1; j <= 2*NUM_LEVELSETS; j++ )
      energyFile << energyComponents[j];
    energyFile << "volume deviation energies: ";
    for ( int j = 2 * NUM_LEVELSETS + 1; j <= 3 * NUM_LEVELSETS; j++ )
      energyFile << energyComponents[j];
#ifdef USE_LANDMARKS
    energyFile << endl << "landmark deviation energies: " << energyComponents[3*NUM_LEVELSETS+1];
#endif
    energyFile.close();

    // save the variations
    aol::MultiVector<RealType> geodesicGrad( geodesic, aol::STRUCT_COPY );
    de.apply( geodesic, geodesicGrad );
    const int timePointInc = 1 << ( _finestTimeLevel - TimeLevel );
    for ( int i = 0; i < _maxNumberOfTimePoints; i += timePointInc ) {
      aol::MultiVector<RealType> dispGrad;
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        dispGrad.appendReference( geodesicGrad[NUM_LEVELSETS*((1<<TimeLevel)+1)+i*ConfiguratorType::Dim/timePointInc+j], false );
      dispGrad *= 500;
      //saveDisplacementAsChessPattern<ConfiguratorType,PNG>( dispGrad, grid, _destDirectory, "displGrad", Index * 100 + i );
      char fileName[1024];
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        sprintf( fileName, "levelset%dGrad", j );
        //saveAsColorImage<ConfiguratorType>( geodesicGrad[j*NUM_LEVELSETS+i], grid, _destDirectory, fileName, Index * 100 + i, -.0004, .0004 );
      }
    }
  }

  //! Saves the geodesic
  void saveResults( const int TimeLevel, const int SpaceLevel, const int /*Index*/ ) {
    const int timePointInc = 1 << ( _finestTimeLevel - TimeLevel );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _maxNumberOfTimePoints; i += timePointInc ) {
      // save the objects
      char fileName[1024];
      ArrayType image( _signedDistFuncs[0][0]->operator[]( SpaceLevel ), aol::STRUCT_COPY );
      for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
        sprintf( fileName, "object%d_", j );
        image = _signedDistFuncs[j][i]->operator[]( SpaceLevel );
        if ( ConfiguratorType::Dim == qc::QC_2D && _inputRedistancing )
          image.clamp( 0., _grids[SpaceLevel]->H() );
        //save<ConfiguratorType>( image, *_grids[SpaceLevel], _destDirectory, fileName, i );
        saveAsImage<ConfiguratorType,PGM>( image, *_grids[SpaceLevel], _destDirectory, fileName, i );
        //saveAsImage<ConfiguratorType,PNG>( image, *_grids[SpaceLevel], _destDirectory, fileName, Index * 100 + i );
        sprintf( fileName, "levelset%dFunc", j );
        //saveAsColorImage<ConfiguratorType>( _signedDistFuncs[j][i]->operator[]( SpaceLevel ), *_grids[SpaceLevel], _destDirectory, fileName, Index * 100 + i, -.6, .2 );
      }

      if ( i > 0 ) {
        // save the displacements
        aol::MultiVector<RealType> displacement;
        _displacements[i]->setCurLevel( SpaceLevel );
        _displacements[i]->appendReferencesTo( displacement );
        saveDisplacement<ConfiguratorType>( displacement, *_grids[SpaceLevel], _destDirectory, "displacement", i );
        //saveDisplacementAsChessPattern<ConfiguratorType,PNG>( displacement, *_grids[SpaceLevel], _destDirectory, "dispChessPat", Index * 100 + i );

        // save the back-deformed objects
        for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
          qc::DeformImage<ConfiguratorTypeMultiLin>( _signedDistFuncs[j][i]->operator[]( SpaceLevel ), cubicGrid( *_grids[SpaceLevel] ), image, displacement );
          if ( ConfiguratorType::Dim == qc::QC_2D && _inputRedistancing )
            image.clamp( 0., _grids[SpaceLevel]->H() );
          char fileName[1024];
          sprintf( fileName, "backDeformed%dImage", j );
          //save<ConfiguratorType>( image, *_grids[Level], _destDirectory, fileName, i );
          //saveAsImage<ConfiguratorType,PNG>( image, *_grids[Level], _destDirectory, fileName, Level * 10000 + Index * 100 + i );
        }
      }
    }
  }

  //! Does "NumSteps" gradient descent steps for the geodesic energy.
  void smoothGeodesic( const int TimeLevel, const int SpaceLevel, const int NumSteps ) {
    cerr << "Smoothing on time-like level " << TimeLevel << " and space-like level " << SpaceLevel << endl;
    const typename ConfiguratorType::InitType &grid( *_grids[SpaceLevel] );
    aol::MultiVector<RealType> geodesic, defEnergyWeights;
    computeDeformationEnergyWeights( TimeLevel, SpaceLevel, 1, defEnergyWeights );
    putGeodesicIntoMultiVector( TimeLevel, SpaceLevel, geodesic );
    aol::MultiVector<RealType> initialGeodesic( geodesic, aol::DEEP_COPY );

    LevelsetGeodesicEnergy<ConfiguratorType,HeavisideFunctionType,MaterialLawType> e( grid, 1 << TimeLevel, _gamma, _mu, _nu, _lambda, _normReg, grid.H() * _epsFactor, _defEnergyParam1, _defEnergyParam2, defEnergyWeights, _occlusion );
    LevelsetGeodesicGradient<ConfiguratorType,HeavisideFunctionType,MaterialLawType> de( grid, 1 << TimeLevel, _gamma, _mu, _nu, _lambda, _normReg, grid.H() * _epsFactor, _defEnergyParam1, _defEnergyParam2, defEnergyWeights, _occlusion );

#ifdef USE_LANDMARKS
    aol::Vector<int> tauBlocks( 3 );
    tauBlocks[2] = NUM_LEVELSETS * ( ( 1 << TimeLevel ) + 1 ) * ConfiguratorType::Dim;
#else
    aol::Vector<int> tauBlocks( 2 );
#endif
    tauBlocks[0] = NUM_LEVELSETS * ( ( 1 << TimeLevel ) + 1 );
    tauBlocks[1] = ( ( 1 << TimeLevel ) + 1 ) * ConfiguratorType::Dim;
    aol::GradientDescentComponentWiseTimestepControlled<ConfiguratorType> gradientDescent( grid, e, de, tauBlocks, NumSteps );
    // attention! If smoothing is not deactivated, the smoothed descent direction for the level set functions is so bad that the gradient descent fails!
    gradientDescent.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::USE_NONLINEAR_CG | aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::DO_NOT_SMOOTH_DESCENT_DIRECTION );
    gradientDescent.setTauMin( 1.e-20 ); // close to complete compression there might be very strong changes of the energy!
    //gradientDescent.setTimestepController( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::GRADIENT_SIMPLE );

    //aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > gradientDescent( grid, e, de, NumSteps );
    //aol::QuasiNewtonIterationComponentWiseTimestepControlled<RealType> gradientDescent( e, de, tauBlocks, NumSteps, 1.e-6, 50 );

    gradientDescent.apply( initialGeodesic, geodesic );
  }

  //! Restricts the spatial discretization by one level
  void spaceRestrict( const int TimeLevel, const int SpaceLevelOld ) {
    const int timePointInc = 1 << ( _finestTimeLevel - TimeLevel );
    const int SpaceLevelNew = SpaceLevelOld - 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = timePointInc; i < _maxNumberOfTimePoints; i += timePointInc ) {
      // restrict the levelset functions
      if ( i < _maxNumberOfTimePoints - 1 ) // the objects at the end (and beginning) of the geodesic stay unchanged
        for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
          _signedDistFuncs[j][i]->setCurLevel( SpaceLevelOld );
          ArrayType highlyResolvedFunction( _signedDistFuncs[j][i]->current(), aol::DEEP_COPY );
          _signedDistFuncs[j][i]->levRestrict();
          // now compute the residuum and save it on the fine level
          _signedDistFuncs[j][i]->levProlongate();
          _signedDistFuncs[j][i]->current() -= highlyResolvedFunction;
          _signedDistFuncs[j][i]->setCurLevel( SpaceLevelNew );
        }
      // restrict the displacements
      _displacements[i]->setCurLevel( SpaceLevelOld );
      aol::MultiVector<RealType> displacementAux;
      _displacements[i]->appendReferencesTo( displacementAux );
      aol::MultiVector<RealType> highlyResolvedDisplacement( displacementAux, aol::DEEP_COPY );
      _displacements[i]->levRestrict();
      // now compute the residuum and save it on the fine level
      _displacements[i]->levProlongate();
      displacementAux -= highlyResolvedDisplacement;
      _displacements[i]->setCurLevel( SpaceLevelNew );
    }
  }

  //! Prolongates the spatial discretization by one level
  void spaceProlongate( const int TimeLevel, const int SpaceLevelOld, const bool DoNotAddOldResiduums = false ) {
    const int timePointInc = 1 << ( _finestTimeLevel - TimeLevel );
    const int spaceLevelNew = SpaceLevelOld + 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = timePointInc; i < _maxNumberOfTimePoints; i += timePointInc ) {
      // prolongate the levelset functions
      if ( i < _maxNumberOfTimePoints - 1 ) // the objects at the end (and beginning) of the geodesic stay unchanged
        for ( int j = 0; j < NUM_LEVELSETS; j++ ) {
          _signedDistFuncs[j][i]->setCurLevel( spaceLevelNew );
          ArrayType residuum( _signedDistFuncs[j][i]->current(), aol::DEEP_COPY );
          _signedDistFuncs[j][i]->setCurLevel( SpaceLevelOld );
          _signedDistFuncs[j][i]->levProlongate();
          if ( !DoNotAddOldResiduums )
            // correct prolongation with the residuum which was saved on the fine level
            _signedDistFuncs[j][i]->current() -= residuum;
        }
      // prolongate the displacements
      _displacements[i]->setCurLevel( spaceLevelNew );
      aol::MultiVector<RealType> displacementAux;
      _displacements[i]->appendReferencesTo( displacementAux );
      aol::MultiVector<RealType> residuum( displacementAux, aol::DEEP_COPY );
      _displacements[i]->setCurLevel( SpaceLevelOld );
      _displacements[i]->levProlongate();
      if ( !DoNotAddOldResiduums ) {
        // correct prolongation with the residuum which was saved on the fine level
        displacementAux -= residuum;
        RealType factor = 1.;
        std::vector<qc::CoordType> nodes;
        std::vector<qc::Element> elems;
        while ( !checkInterpenetrations( displacementAux, *_grids[spaceLevelNew], nodes, elems, 3e-5 ) && ( factor > 0.02 ) ) {
          factor /= 2.;
          displacementAux.addMultiple( residuum, factor );
        }
        if ( factor < 0.02 )
          displacementAux.addMultiple( residuum, factor );
      }

      // remedy local interpenetrations of matter, which were overlooked on coarse grid
      aol::MultiVector<RealType> displacement;
      _displacements[i]->appendReferencesTo( displacement );
      // removeInterpenetrations( displacement, *_grids[SpaceLevelOld+1] );
    }
  }

  //! Restricts the time discretization by one level
  void timeRestrict( const int TimeLevelOld, const int SpaceLevel ) {
    const int timePointInc = 1 << ( _finestTimeLevel - TimeLevelOld + 1 );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = timePointInc; i < _maxNumberOfTimePoints; i += timePointInc ) {
      aol::MultiVector<RealType> displacementFirstHalf, displacementSecondHalf;
      _displacements[i]->setCurLevel( SpaceLevel );
      _displacements[i-(timePointInc/2)]->setCurLevel( SpaceLevel );
      _displacements[i]->appendReferencesTo( displacementSecondHalf );
      _displacements[i-(timePointInc/2)]->appendReferencesTo( displacementFirstHalf );
      aol::MultiVector<RealType> displacementSecondHalfOld( displacementSecondHalf, aol::DEEP_COPY );

      qc::ConcatAndSmoothlyExtendDeformations<ConfiguratorType>( displacementSecondHalfOld, displacementFirstHalf, *_grids[SpaceLevel], displacementSecondHalf );

      // remedy local interpenetrations of matter
      // (multilinear approximation of concatenation of two noninterpenetrating multilinear deformations can result in interpenetrations)
      //removeInterpenetrations( displacementSecondHalf, *_grids[SpaceLevel] );
    }
 }

  //! Prolongates the time discretization by one level
  void timeProlongate( const int TimeLevelOld, const int SpaceLevel, const bool DoNotUseOldResiduums = false ) {
    const int timePointInc = 1 << ( _finestTimeLevel - TimeLevelOld );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = timePointInc; i < _maxNumberOfTimePoints; i += timePointInc ) {
      aol::MultiVector<RealType> displacementFirstHalf, displacementSecondHalf;
      _displacements[i]->setCurLevel( SpaceLevel );
      _displacements[i-(timePointInc/2)]->setCurLevel( SpaceLevel );
      _displacements[i]->appendReferencesTo( displacementSecondHalf );
      _displacements[i-(timePointInc/2)]->appendReferencesTo( displacementFirstHalf );
      aol::MultiVector<RealType> displacementSecondHalfOld( displacementSecondHalf, aol::DEEP_COPY );

      if ( DoNotUseOldResiduums ) {
        displacementFirstHalf = displacementSecondHalfOld;
        displacementFirstHalf *= .5;
        removeInterpenetrations( displacementFirstHalf, *_grids[SpaceLevel] );
        for ( int j = 0; j < NUM_LEVELSETS; j++ )
          qc::InvDeformAndSmoothlyExtendImage<ConfiguratorTypeMultiLin>( _signedDistFuncs[j][i-timePointInc]->operator[]( SpaceLevel ), cubicGrid( *_grids[SpaceLevel] ), _signedDistFuncs[j][i-(timePointInc/2)]->operator[]( SpaceLevel ), displacementFirstHalf );
      }

      qc::InvConcatAndSmoothlyExtendDeformations<ConfiguratorTypeMultiLin>( displacementSecondHalfOld, displacementFirstHalf, cubicGrid( *_grids[SpaceLevel] ), displacementSecondHalf );

      // remedy local interpenetrations of matter, which resulted from the inverse deformation
      removeInterpenetrations( displacementSecondHalf, *_grids[SpaceLevel] );
    }
  }

  //! Performs a space-like V-cycle
  void spaceVCycle( const int TimeLevel, const int SpaceLevel, const int NumSteps, const bool OnlyRefiningBranch = false ) {
    if ( !OnlyRefiningBranch ) {
      for ( int level = SpaceLevel; level > _coarsestSpaceLevel; level-- ) {
        smoothGeodesic( TimeLevel, level, NumSteps );
        spaceRestrict( TimeLevel, level );
      }
      smoothGeodesic( TimeLevel, _coarsestSpaceLevel, NumSteps ); // so that on the coarsest grid level there is also double smoothing...
    }
    for ( int level = _coarsestSpaceLevel; level < SpaceLevel; level++ ) {
      smoothGeodesic( TimeLevel, level, NumSteps );
      spaceProlongate( TimeLevel, level, OnlyRefiningBranch );
    }
    smoothGeodesic( TimeLevel, SpaceLevel, NumSteps );
  }

  //! Performs a time-like V-cycle, whose smoothing represents a space-like V-cycle
  void timeVCycle( const int TimeLevel, const int SpaceLevel, const int NumSpaceVCycles, const int NumSteps, const bool OnlyRefiningBranch = false ) {
    if ( !OnlyRefiningBranch ) {
      for ( int level = TimeLevel; level > _coarsestTimeLevel; level-- ) {
        for ( int i = 0; i < NumSpaceVCycles; i++ ) {
          cerr << aol::color::blue << "Do space-like V-cycle number " << i + 1 << " of " << NumSpaceVCycles << aol::color::reset << endl;
          spaceVCycle( level, SpaceLevel, NumSteps );
        }
        timeRestrict( level, SpaceLevel );
      }
      for ( int i = 0; i < NumSpaceVCycles; i++ ) {
        cerr << aol::color::blue << "Do space-like V-cycle number " << i + 1 << " of " << NumSpaceVCycles << aol::color::reset << endl;
        spaceVCycle( _coarsestTimeLevel, SpaceLevel, NumSteps ); // so that on the coarsest time level there is also double smoothing...
      }
    }
    for ( int level = _coarsestTimeLevel; level < TimeLevel; level++ ) {
      for ( int i = 0; i < NumSpaceVCycles; i++ ) {
        cerr << aol::color::blue << "Do space-like V-cycle number " << i + 1 << " of " << NumSpaceVCycles << aol::color::reset << endl;
        spaceVCycle( level, SpaceLevel, NumSteps, OnlyRefiningBranch && level == _coarsestTimeLevel );
      }
      timeProlongate( level, SpaceLevel, OnlyRefiningBranch );
    }
    for ( int i = 0; i < NumSpaceVCycles; i++ ) {
      cerr << aol::color::blue << "Do space-like V-cycle number " << i + 1 << " of " << NumSpaceVCycles << aol::color::reset << endl;
      spaceVCycle( TimeLevel, SpaceLevel, NumSteps );
    }
  }

  //! Computes the geodesic
  void execute() {
#ifdef _OPENMP
    omp_set_nested( 1 );
#endif
    cerr << "relax the geodesic..." << endl;
    for ( int i = 0; i < _maxTimeVCycles; i++ ) {
      cerr << aol::color::red << "Do time-like V-cycle number " << i + 1 << " of " << _maxTimeVCycles << aol::color::reset << endl;
      timeVCycle( _finestTimeLevel, _finestSpaceLevel, _maxSpaceVCycles, _maxSmoothingSteps, i == 0 );
      saveResults( _finestTimeLevel, _finestSpaceLevel, i );
    }
  }
};


#endif
