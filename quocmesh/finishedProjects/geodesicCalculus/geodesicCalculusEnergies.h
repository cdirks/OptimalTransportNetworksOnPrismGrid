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

#ifndef __GEODESICCALCULUSENERGIES_H
#define __GEODESICCALCULUSENERGIES_H

#include <aol.h>
#include <convolution.h>
#include <RudinOsherFatemi.h>
#include <hyperelastic.h>
#include <deformations.h>
#include <parameterParser.h>
#include <gradientDescent.h>
#include <trustRegionMethod.h>
#include <FEOpInterface.h>

#include <gradientflow.h>

//#define TV_OF_SMOOTHED_CHAR_FUN

//! represents the hyperelastic energy density \f$ W(I_1,I_2,I_3) = \log(\sqrt{1+\exp(2a(b-I_3))}) \f$,
/** where a and b are passed to the constructor as Slope and Shift and \f$ I_1,I_2,I_3 \f$ stand for invariants
  * \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$, \f$ I_3=det(\nabla\phi) \f$
  * of a deformation gradient \f$ \nabla\phi \f$.
  *
  * \f$ W \f$ is a smoothed version of \f$ a(b-I_3) \f$ for \f$ I_3 \leq b \f$ and 0 else.
  *
  * \author wirth
  */
template <typename ConfiguratorType>
class AntiSelfPenetrationEnergyDensity {
private:
  typedef typename ConfiguratorType::RealType RealType;
  RealType _slope, _shift;
  RealType _cutOff, _cutOffError;
public:
  AntiSelfPenetrationEnergyDensity( const RealType Slope, const RealType Shift ) :
    _slope( Slope ),
    _shift( Shift ),
    _cutOff( 5. ),
    _cutOffError( .5 * log( 1. + exp( 2. * _cutOff ) ) - _cutOff ) {}

  inline RealType evaluate ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                             const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    RealType x = _slope * ( _shift - I3 );
    return x < _cutOff ? .5 * log( 1. + exp( 2. * x ) ) : x + _cutOffError;
  }

  inline RealType evaluateAtQuadPoint ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                        const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    RealType x = _slope * ( _shift - I3 );
    return x < _cutOff ? .5 * log( 1. + exp( 2. * x ) ) : x + _cutOffError;
  }

  inline aol::Vec3<RealType> evaluateDerivative ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    RealType x = _slope * ( I3 - _shift );
    return aol::Vec3<RealType>( 0., 0., x > - _cutOff ? - _slope / ( 1 + exp( 2. * x ) ) : - _slope );
  }

  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    RealType x = _slope * ( I3 - _shift );
    return aol::Vec3<RealType>( 0., 0., x > - _cutOff ? - _slope / ( 1 + exp( 2. * x ) ) : - _slope );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivative ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                            const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    RealType x = _slope * ( _shift - I3 ), aux = exp( 2. * x );
    return aol::Matrix33<RealType>( 0., 0., 0., 0., 0., 0., 0., 0., x < _cutOff ? 2. * aol::Sqr( _slope ) * aux / aol::Sqr( 1 + aux ) : 0. );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivativeAtQuadPoint ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                                       const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    RealType x = _slope * ( _shift - I3 ), aux = exp( 2. * x );
    return aol::Matrix33<RealType>( 0., 0., 0., 0., 0., 0., 0., 0., x < _cutOff ? 2. * aol::Sqr( _slope ) * aux / aol::Sqr( 1 + aux ) : 0. );
  }
};

//! ===========================================================================
//! ResolveSelfPenetrationEnergy

template <typename ConfiguratorType>
class ResolveSelfPenetrationEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const aol::MultiVector<RealType> _originalDisplacement;
  const AntiSelfPenetrationEnergyDensity<ConfiguratorType> _antiInterpenetrationEnergyDensity;
  const qc::HyperelasticEnergy<ConfiguratorType,AntiSelfPenetrationEnergyDensity<ConfiguratorType> > _antiInterpenetrationEnergy;
  const aol::MassOp<ConfiguratorType> _fidelityOp;

public:
  // note: the first MorphingDisplacement is a dummy
  ResolveSelfPenetrationEnergy ( const typename ConfiguratorType::InitType &Grid,
                                  const aol::MultiVector<RealType> &OriginalDisplacement,
                                  const RealType Slope,
                                  const RealType Shift ) :
    _originalDisplacement( OriginalDisplacement ),
    _antiInterpenetrationEnergyDensity( Slope, Shift ),
    _antiInterpenetrationEnergy( Grid, _antiInterpenetrationEnergyDensity ),
    _fidelityOp( Grid, aol::ASSEMBLED ){}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    _antiInterpenetrationEnergy.applyAdd( Arg, Dest );
    aol::MultiVector<RealType> aux( Arg, aol::DEEP_COPY ), aux2( Arg, aol::STRUCT_COPY );
    aux -= _originalDisplacement;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      _fidelityOp.apply( aux[i], aux2[i] );
    Dest += .5 * ( aux * aux2 );
  }
};

template <typename ConfiguratorType>
class ResolveSelfPenetrationDeriv :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const aol::MultiVector<RealType> _originalDisplacement;
  const AntiSelfPenetrationEnergyDensity<ConfiguratorType> _antiInterpenetrationEnergyDensity;
  const qc::HyperelasticGradient<ConfiguratorType,AntiSelfPenetrationEnergyDensity<ConfiguratorType> > _antiInterpenetrationDeriv;
  const aol::MassOp<ConfiguratorType> _fidelityOp;

public:
  // note: the first MorphingDisplacement is a dummy
  ResolveSelfPenetrationDeriv( const typename ConfiguratorType::InitType &Grid,
                                const aol::MultiVector<RealType> &OriginalDisplacement,
                                const RealType Slope,
                                const RealType Shift ) :
    _originalDisplacement( OriginalDisplacement ),
    _antiInterpenetrationEnergyDensity( Slope, Shift ),
    _antiInterpenetrationDeriv( Grid, _antiInterpenetrationEnergyDensity ),
    _fidelityOp( Grid, aol::ASSEMBLED ){}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    _antiInterpenetrationDeriv.applyAdd( Arg, Dest );
    aol::MultiVector<RealType> aux( Arg, aol::DEEP_COPY );
    aux -= _originalDisplacement;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      _fidelityOp.applyAdd( aux[i], Dest[i] );
  }
};

template <typename ConfiguratorType>
class ResolveSelfPenetrationHessian :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::BlockOpBase<typename ConfiguratorType::RealType,typename ConfiguratorType::MatrixType> >{

private:
  typedef typename ConfiguratorType::RealType RealType;

  const AntiSelfPenetrationEnergyDensity<ConfiguratorType> _antiInterpenetrationEnergyDensity;
  const qc::HyperelasticHessian<ConfiguratorType,AntiSelfPenetrationEnergyDensity<ConfiguratorType> > _antiInterpenetrationHessian;
  const aol::MassOp<ConfiguratorType> _fidelityOp;

public:
  // note: the first MorphingDisplacement is a dummy
  ResolveSelfPenetrationHessian( const typename ConfiguratorType::InitType &Grid,
                                  const aol::MultiVector<RealType> &/*OriginalDisplacement*/,
                                  const RealType Slope,
                                  const RealType Shift ) :
    _antiInterpenetrationEnergyDensity( Slope, Shift ),
    _antiInterpenetrationHessian( Grid, _antiInterpenetrationEnergyDensity ),
    _fidelityOp( Grid, aol::ASSEMBLED ){}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::BlockOpBase<RealType,typename ConfiguratorType::MatrixType> &Dest ) const {
    _antiInterpenetrationHessian.applyAdd( Arg, Dest );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      _fidelityOp.assembleAddMatrix( Dest.getReference( i, i ) );
  }
};

template <typename ConfiguratorType>
void resolveSelfPenetration( aol::MultiVector<typename ConfiguratorType::RealType> &Displacement,
                              const typename ConfiguratorType::InitType &Grid,
                              const typename ConfiguratorType::RealType Slope,
                              const typename ConfiguratorType::RealType Shift,
                              const bool QuietMode = true ) {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef typename ConfiguratorType::MatrixType MatrixType;
  // perform minimization from coarse scale to fine scale
  qc::MultiDimMultilevelArray<RealType,ArrayType> origDisp( Grid, Displacement.numComponents() ), newDisp( Grid, Displacement.numComponents() );
  for ( int i = 0; i < Displacement.numComponents(); i++ )
    origDisp.getArray( i ) = Displacement[i];
  if ( Grid.getGridDepth() > 3 )
    origDisp.levRestrict( 3, Grid.getGridDepth() );
  for ( int level = aol::Min( 3, Grid.getGridDepth() ); level <= Grid.getGridDepth(); level++ ) {
    if ( !QuietMode )
      cerr<<"resolve self-interpenetration on grid level "<<level<<" ...\n";
    // set the current scale
    typename ConfiguratorType::InitType grid( level, ConfiguratorType::Dim );
    origDisp.setCurLevel( level );
    newDisp.setCurLevel( level );
    aol::MultiVector<RealType> disp, dispNew;
    origDisp.appendReferencesTo( disp );
    newDisp.appendReferencesTo( dispNew );
    // minimize on this scale
    ResolveSelfPenetrationEnergy<ConfiguratorType>    e( grid, disp, Slope, Shift );
    ResolveSelfPenetrationDeriv<ConfiguratorType>    de( grid, disp, Slope, Shift );
    ResolveSelfPenetrationHessian<ConfiguratorType> d2e( grid, disp, Slope, Shift );
    aol::TrustRegionMethod<RealType, aol::MultiVector<RealType>, aol::BlockOpBase<RealType, MatrixType>, aol::DiagonalPreconditioner<aol::MultiVector<RealType> >, MatrixType> minimizer( e, de, d2e, new aol::BlockMatrix<MatrixType>( grid ), 3000, 1.e-6 );
    minimizer.setIterativeMode( ( level > 6 ) || ( ConfiguratorType::Dim > qc::QC_2D ) );
    minimizer.setQuietMode( QuietMode );
    aol::MultiVector<RealType> initialGuess( dispNew, aol::DEEP_COPY );
    minimizer.apply( initialGuess, dispNew );
    newDisp.levProlongate();
  }
  for ( int i = 0; i < Displacement.numComponents(); i++ )
    Displacement[i] = newDisp.getArray( i );
}

//! ===========================================================================
template <typename ConfiguratorType>
class MismatchGradientWRTnondefCharFun :
  public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType _grid;
  // the smoothed characteristic function to be matched with
  const aol::Vector<RealType> _uTemplate;
  // the matching displacement
  const aol::MultiVector<RealType> _displacement;
  const aol::MassOp<ConfiguratorType> _massOp;
  qc::LinearConvolutionSmoothOp<RealType,ConfiguratorType::Dim> _gaussKernelSmoother;

public:
  MismatchGradientWRTnondefCharFun( const typename ConfiguratorType::InitType &Grid,
                                    const aol::Vector<RealType> &UTemplate,
                                    const aol::MultiVector<RealType> &Displacement,
                                    const RealType Epsilon ) :
    _grid( Grid ),
    _uTemplate( UTemplate ),
    _displacement( Displacement ),
    _massOp( Grid ),
    _gaussKernelSmoother( Grid.getNumX(), Grid.getNumY() ) {
    _gaussKernelSmoother.setSigma( Epsilon );
  }

  // Arg is the smoothed nondeformed characteristic function
  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // compute -2(G*u_k\circ\phi_k-G*u_{k-1})
    aol::Vector<RealType> deformedImage( Arg, aol::STRUCT_COPY );
    qc::DeformImage<ConfiguratorType>( _uTemplate, _grid, deformedImage, _displacement, true );
    aol::Vector<RealType> difference( Arg, aol::DEEP_COPY );
    difference -= deformedImage;
    difference *= 2.;
    // convolve result with G(-\cdot); use that Gauss kernel is point symmetric
    aol::Vector<RealType> smoothedDifference( Arg, aol::STRUCT_COPY );
    _gaussKernelSmoother.apply( difference, smoothedDifference );
    // result is the application of a mass op
    _massOp.applyAdd( smoothedDifference, Dest );
  }
};

template <typename ConfiguratorType>
class MismatchGradientWRTdefCharFun :
  public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType _grid;
  // the smoothed characteristic function to be matched with
  const aol::Vector<RealType> &_uFixed;
  // the matching displacement
  const aol::MultiVector<RealType> &_displacement;
  const aol::MassOp<ConfiguratorType> _massOp;
  qc::LinearConvolutionSmoothOp<RealType,ConfiguratorType::Dim> _gaussKernelSmoother;

public:
  MismatchGradientWRTdefCharFun( const typename ConfiguratorType::InitType &Grid,
                                 const aol::Vector<RealType> &UFixed,
                                 const aol::MultiVector<RealType> &Displacement,
                                 const RealType Epsilon ) :
    _grid( Grid ),
    _uFixed( UFixed ),
    _displacement( Displacement ),
    _massOp( Grid ),
    _gaussKernelSmoother( Grid.getNumX(), Grid.getNumY() ) {
    _gaussKernelSmoother.setSigma( Epsilon );
  }

  // Arg is the smoothed characteristic function to be deformed
  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // compute -2(G*u_k-G*u_{k-1}\circ\phi_k^{-1})|\det\nabla\phi_k^{-1}|
    aol::Vector<RealType> deformedImage( Arg, aol::STRUCT_COPY );
    qc::InvDeformImage<ConfiguratorType,aol::Vector<RealType> >( _uFixed, _grid, deformedImage, _displacement );
    aol::Vector<RealType> difference( Arg, aol::DEEP_COPY );
    difference -= deformedImage;
    difference *= 2.;
    aol::MultiVector<RealType> identity( _grid ), invDisp( _grid );
    qc::DataGenerator<ConfiguratorType>( _grid ).generateIdentity( identity );
    // compute \phi_k^{-1} in "invDisp"
    qc::InvConcatAndSmoothlyExtendDeformations<ConfiguratorType>( identity, _displacement, _grid, invDisp );
    invDisp -= identity;
    // compute \det\nabla\phi_k^{-1} pointwise in "detInvDef"
    aol::Vector<RealType> detInvDef( Arg, aol::STRUCT_COPY );
    qc::DeformedVolumeIntegrator<ConfiguratorType>( _grid ).applyAddIntegrand( invDisp, detInvDef );
    // compute -2(G*u_k-G*u_{k-1}\circ\phi_k^{-1})|\det\nabla\phi_k^{-1}|
    for ( int i = 0; i < Arg.size(); i++ )
      difference[i] *= detInvDef[i];
    // convolve result with G(-\cdot); use that Gauss kernel is point symmetric
    aol::Vector<RealType> smoothedDifference( Arg, aol::STRUCT_COPY );
    _gaussKernelSmoother.apply( difference, smoothedDifference );
    // result is the application of a mass op
    _massOp.applyAdd( smoothedDifference, Dest );
  }
};
//! ===========================================================================
template <typename ConfiguratorType>
class RegLengthEnergy :
  public aol::FENonlinIntegrationVectorInterface<ConfiguratorType,RegLengthEnergy<ConfiguratorType>,ConfiguratorType::Dim> {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const RealType _delta;
  qc::LinearConvolutionSmoothGradientOp<RealType,ConfiguratorType::Dim> _gaussKernelSmootherGrad;

public:
  RegLengthEnergy( const typename ConfiguratorType::InitType &Grid, const RealType Epsilon, const RealType Delta ) :
    aol::FENonlinIntegrationVectorInterface<ConfiguratorType,RegLengthEnergy<ConfiguratorType>,ConfiguratorType::Dim>( Grid ),
    _delta( Delta ),
    _gaussKernelSmootherGrad( Grid.getNumX(), Grid.getNumY() ) {
    _gaussKernelSmootherGrad.setSigma( Epsilon );
  }

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El, int QuadPoint,
                              const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    RealType result = 0;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      result += aol::Sqr( DiscFuncs[i].evaluateAtQuadPoint( El, QuadPoint ) );
    result += aol::Sqr( _delta );
    return sqrt( result );
  }

  void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/ ) const {}

  void applyAdd( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector<RealType> aux( ConfiguratorType::Dim, Arg.size() );
    _gaussKernelSmootherGrad.apply( Arg, aux );

    aol::FENonlinIntegrationVectorInterface<ConfiguratorType,RegLengthEnergy<ConfiguratorType>,ConfiguratorType::Dim>::applyAdd( aux, Dest );
  }
};

template <typename ConfiguratorType>
class RegLengthEnergyVariation :
  public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const RealType _delta;
  const aol::MassOp<ConfiguratorType> _massOp;
  qc::LinearConvolutionSmoothGradientOp<RealType,ConfiguratorType::Dim> _gaussKernelSmootherGrad;

public:
  RegLengthEnergyVariation( const typename ConfiguratorType::InitType &Grid, const RealType Epsilon, const RealType Delta ) :
    _delta( Delta ),
    _massOp( Grid ),
    _gaussKernelSmootherGrad( Grid.getNumX(), Grid.getNumY() ) {
    _gaussKernelSmootherGrad.setSigma( Epsilon );
  }

  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // compute (\nabla G*u/\sqrt{|\nabla G*u|^2+\delta^2})
    aol::MultiVector<RealType> smoothedGradNormalized( ConfiguratorType::Dim, Arg.size() );
    _gaussKernelSmootherGrad.apply( Arg, smoothedGradNormalized );
    for ( int i = 0; i < smoothedGradNormalized[0].size(); i++ ) {
      RealType normSqr = aol::Sqr( _delta );
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        normSqr += aol::Sqr( smoothedGradNormalized[j][i] );
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        smoothedGradNormalized[j][i] /= sqrt( normSqr );
    }
    // convolve result with \nabla G(-\cdot)=-\nabla G(\cdot)
    aol::Vector<RealType> localNormalField( Arg, aol::STRUCT_COPY );
    aol::MultiVector<RealType> aux( ConfiguratorType::Dim, Arg.size() );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
      _gaussKernelSmootherGrad.apply( smoothedGradNormalized[j], aux );
      localNormalField -= aux[j];
    }
    // result is the application of a mass op
    _massOp.applyAdd( localNormalField, Dest );
  }
};
//! ===========================================================================
template <typename ConfiguratorType, typename MaterialLawType>
class GeodesicDeformationCorrectionEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const int _numPolygonLines, _numShapes;
  const MaterialLawType _elasticEnergyDensity;
  const aol::RandomAccessContainer< aol::MultiVector<typename ConfiguratorType::RealType> > &_displacements;

public:
  GeodesicDeformationCorrectionEnergy( const typename ConfiguratorType::InitType &Grid,
                                       const aol::RandomAccessContainer<aol::MultiVector<typename ConfiguratorType::RealType> > &Displacements,
                                       const MaterialLawType &ElasticEnergyDensity ) :
    _grid( Grid ),
    _numPolygonLines( Displacements.size() ),
    _numShapes( _numPolygonLines + 1 ),
    _elasticEnergyDensity( ElasticEnergyDensity ),
    _displacements( Displacements ){}

  // Arg contains (K-1) displacements, ie. dim*(K-1) displacement components
  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 0; k < _numPolygonLines; k++ ) {
      aol::MultiVector<RealType> preCorrection, postCorrection;
      // extract the pre-correction
      if ( k > 0 )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          preCorrection.appendReference( Arg[ConfiguratorType::Dim*(k-1)+j] );
      else
        preCorrection.reallocate( _grid );
      // extract the post-correction
      if ( k < _numPolygonLines - 1 )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          postCorrection.appendReference( Arg[ConfiguratorType::Dim*k+j] );
      else
        postCorrection.reallocate( _grid );
      // compute deformation energy of kth deformation
      aol::Scalar<RealType> dest;
      postCorrection.appendReference( preCorrection );
      qc::HyperelasticPrePostDeformEnergy<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity, _displacements[k] ).applyAdd( postCorrection, dest );
#ifdef _OPENMP
#pragma omp critical (GeodesicDeformationCorrectionEnergyApplyAdd)
#endif
      Dest += dest;
    }
    if ( aol::isNaN( Dest.v ) )
      Dest = aol::NumberTrait<RealType>::Inf;
  }
};

template <typename ConfiguratorType, typename MaterialLawType>
class GeodesicDeformationCorrectionGradient :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const int _numPolygonLines, _numShapes;
  const MaterialLawType _elasticEnergyDensity;
  const aol::RandomAccessContainer< aol::MultiVector<typename ConfiguratorType::RealType> > &_displacements;

public:
  GeodesicDeformationCorrectionGradient( const typename ConfiguratorType::InitType &Grid,
                                         const aol::RandomAccessContainer<aol::MultiVector<typename ConfiguratorType::RealType> > &Displacements,
                                         const MaterialLawType &ElasticEnergyDensity ) :
    _grid( Grid ),
    _numPolygonLines( Displacements.size() ),
    _numShapes( _numPolygonLines + 1 ),
    _elasticEnergyDensity( ElasticEnergyDensity ),
    _displacements( Displacements ){}

  // Arg contains (K-1) displacements, ie. dim*(K-1) displacement components
  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 0; k < _numPolygonLines; k++ ) {
      aol::MultiVector<RealType> preCorrection, postCorrection, preCorrGrad, postCorrGrad;
      // extract the pre-correction
      if ( k > 0 )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          preCorrection.appendReference( Arg[ConfiguratorType::Dim*(k-1)+j] );
          preCorrGrad.appendReference( Dest[ConfiguratorType::Dim*(k-1)+j] );
        }
      else
        preCorrection.reallocate( _grid );
      // extract the post-correction
      if ( k < _numPolygonLines - 1 )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ){
          postCorrection.appendReference( Arg[ConfiguratorType::Dim*k+j] );
          postCorrGrad.appendReference( Dest[ConfiguratorType::Dim*k+j] );
        }
      else
        postCorrection.reallocate( _grid );
      // compute gradient of deformation energy of kth deformation wrt pre- and post-correction
      postCorrection.appendReference( preCorrection );
      if ( k > 0 )
        qc::HyperelasticPrePostDeformGradientPreComp<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity, _displacements[k] ).applyAdd( postCorrection, preCorrGrad );
#ifdef _OPENMP
#pragma omp barrier
#endif
      if ( k < _numPolygonLines - 1 )
        qc::HyperelasticPrePostDeformGradientPostComp<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity, _displacements[k] ).applyAdd( postCorrection, postCorrGrad );
    }
  }
};

template <typename ConfiguratorType, typename MaterialLawType, typename SubMatrixType = typename ConfiguratorType::MatrixType>
class GeodesicDeformationCorrectionHessian
      : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::BlockOpBase<typename ConfiguratorType::RealType,SubMatrixType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const int _numPolygonLines, _numShapes;
  const MaterialLawType _elasticEnergyDensity;
  const aol::RandomAccessContainer< aol::MultiVector<typename ConfiguratorType::RealType> > &_displacements;

public:
  GeodesicDeformationCorrectionHessian( const typename ConfiguratorType::InitType &Grid,
                                        const aol::RandomAccessContainer<aol::MultiVector<typename ConfiguratorType::RealType> > &Displacements,
                                        const MaterialLawType &ElasticEnergyDensity ) :
    _grid( Grid ),
    _numPolygonLines( Displacements.size() ),
    _numShapes( _numPolygonLines + 1 ),
    _elasticEnergyDensity( ElasticEnergyDensity ),
    _displacements( Displacements ){}

  // Arg contains (K-1) displacements, ie. dim*(K-1) displacement components
  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::BlockOpBase<RealType,SubMatrixType> &Dest ) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 0; k < _numPolygonLines; k++ ) {
      aol::MultiVector<RealType> preCorrection, postCorrection;
      // extract the pre-correction
      if ( k > 0 )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          preCorrection.appendReference( Arg[ConfiguratorType::Dim*(k-1)+j] );
      else
        preCorrection.reallocate( _grid );
      // extract the post-correction
      if ( k < _numPolygonLines - 1 )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          postCorrection.appendReference( Arg[ConfiguratorType::Dim*k+j] );
      else
        postCorrection.reallocate( _grid );
      // assemble the different components of the Hessian
      if ( k > 0 )
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          for ( int j = 0; j < ConfiguratorType::Dim; j++ )
            qc::HyperelasticPrePostDeformSubHessianPreComp<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity, _displacements[k], preCorrection, postCorrection, i, j ).assembleAddMatrix( Dest.getReference( ConfiguratorType::Dim*(k-1)+i, ConfiguratorType::Dim*(k-1)+j ) );
      if ( ( k > 0 ) && ( k < _numPolygonLines - 1 ) )
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
            qc::HyperelasticPrePostDeformSubHessianMixedPart<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity, _displacements[k], preCorrection, postCorrection, i, j, true ).assembleAddMatrix( Dest.getReference( ConfiguratorType::Dim*k+i, ConfiguratorType::Dim*(k-1)+j ) );
            qc::HyperelasticPrePostDeformSubHessianMixedPart<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity, _displacements[k], preCorrection, postCorrection, i, j, false ).assembleAddMatrix( Dest.getReference( ConfiguratorType::Dim*(k-1)+i, ConfiguratorType::Dim*k+j ) );
          }
#ifdef _OPENMP
#pragma omp barrier
#endif
      if ( k < _numPolygonLines - 1 )
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          for ( int j = 0; j < ConfiguratorType::Dim; j++ )
            qc::HyperelasticPrePostDeformSubHessianPostComp<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity, _displacements[k], preCorrection, postCorrection, i, j ).assembleAddMatrix( Dest.getReference( ConfiguratorType::Dim*k+i, ConfiguratorType::Dim*k+j ) );
    }
  }
};
//! ===========================================================================

/**
 * Given a vector of shapes and deformations which connect these shapes (i.e. a path in shape space),
 * this method redistributes the deformations along the path so that a new path is generated.
 * Shapes contains the zeroth to the kth shape, Displacements[k] connects the kth shape to the (k+1)th.
 * CurrentPercentagesOfDeformation[k] contains the percentage of Displacements[k] wrt the total deformation along the path.
 * WantedPercentagesOfDeformation[k] contains the desired percentage of Displacements[k] wrt the total deformation along the path.
 */
template <typename ConfiguratorType>
void distributeDeformationsAndShapesAlongPath( const typename ConfiguratorType::InitType &Grid,
                                               aol::MultiVector<typename ConfiguratorType::RealType> &Shapes,
                                               aol::RandomAccessContainer< aol::MultiVector<typename ConfiguratorType::RealType> > &Displacements,
                                               aol::Vector<typename ConfiguratorType::RealType> &CurrentPercentagesOfDeformation,
                                               aol::Vector<typename ConfiguratorType::RealType> &WantedPercentagesOfDeformation ) {
  typedef typename ConfiguratorType::RealType RealType;
  int numShapes = Shapes.numComponents();
  aol::RandomAccessContainer<aol::MultiVector<typename ConfiguratorType::RealType> > oldDisplacements( Displacements );
  aol::MultiVector<RealType> oldShapes( Shapes, aol::DEEP_COPY );
  aol::Vector<RealType> currentAccumPerc( numShapes + 1 ), wantedAccumPerc( numShapes + 1 );
  for ( int i = 1; i < numShapes; i++ ) {
    currentAccumPerc[i] = currentAccumPerc[i-1] + CurrentPercentagesOfDeformation[i-1];
    wantedAccumPerc[i] = wantedAccumPerc[i-1] + WantedPercentagesOfDeformation[i-1];
  }
  currentAccumPerc[numShapes] = aol::NumberTrait<RealType>::Inf;
  wantedAccumPerc[numShapes] = aol::NumberTrait<RealType>::Inf;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int i = 0; i < numShapes - 1; i++ ){
    int i1 = 0, i2 = 0;
    while ( currentAccumPerc[i1+1] <= wantedAccumPerc[i] ) i1++;
    while ( currentAccumPerc[i2+1] < wantedAccumPerc[i+1] ) i2++;
    i2 = aol::Min( i2, numShapes - 2 );

    // compute displacement i connecting shape i with shape i+1
    aol::MultiVector<RealType> displ1( oldDisplacements[i1], aol::DEEP_COPY ),
                               displ2( oldDisplacements[i2], aol::DEEP_COPY ),
                               displAux( displ1, aol::STRUCT_COPY );
    displ1 *= ( wantedAccumPerc[i] - currentAccumPerc[i1] ) / CurrentPercentagesOfDeformation[i1];
    displ2 *= ( wantedAccumPerc[i+1] - currentAccumPerc[i2] ) / CurrentPercentagesOfDeformation[i2];
    for ( int j = i2 - 1; j >= i1; j-- ) {
      qc::ConcatAndSmoothlyExtendDeformations<ConfiguratorType>( displ2, oldDisplacements[j], Grid, displAux );
      displ2 = displAux;
    }
    qc::InvConcatAndSmoothlyExtendDeformations<ConfiguratorType>( displ2, displ1, Grid, Displacements[i] );

    // compute shape i+1
    if ( i < numShapes - 2 ) {
      displ2 = oldDisplacements[i2];
      displ2 *= ( wantedAccumPerc[i+1] - currentAccumPerc[i2] ) / CurrentPercentagesOfDeformation[i2];
      qc::InvDeformAndSmoothlyExtendImage<ConfiguratorType,aol::Vector<RealType> >( oldShapes[i2], Grid, Shapes[i+1], displ2 );
    }
  }
}

template <typename ConfiguratorType, typename MaterialLawType>
void correctDeformationsAndShapes ( const typename ConfiguratorType::InitType &Grid,
                                    aol::MultiVector<typename ConfiguratorType::RealType> &Shapes,
                                    aol::RandomAccessContainer< aol::MultiVector<typename ConfiguratorType::RealType> > &Displacements,
                                    MaterialLawType &ElasticEnergyDensity, 
                                    int steps = 50 ) {
  // for the correction we have to prevent non-invertibility!
  typedef qc::AntiCompressionHyperelasticEnergyDensity<ConfiguratorType,MaterialLawType> MatLawType;
  MatLawType elasticEnergyDensity( ElasticEnergyDensity );

  // compute the correction displacements
  // attention! the elastic energy must tend to infinity for total compression, otherwise the method does not work!
  typedef typename ConfiguratorType::RealType RealType;
  typedef aol::SparseMatrix<RealType> SubMatrixType;
  aol::MultiVector<RealType> zeros( ConfiguratorType::Dim * ( Displacements.size() - 1 ), Grid.getNumberOfNodes() ), corrections( zeros );
  GeodesicDeformationCorrectionEnergy<ConfiguratorType,MatLawType> e( Grid, Displacements, elasticEnergyDensity );
  GeodesicDeformationCorrectionGradient<ConfiguratorType,MatLawType> de( Grid, Displacements, elasticEnergyDensity );
  GeodesicDeformationCorrectionHessian<ConfiguratorType,MatLawType,SubMatrixType> d2e( Grid, Displacements, elasticEnergyDensity );

  //alternative minimization approaches
  /*// Quasi-Newton method
  aol::QuasiNewtonIteration<RealType,aol::MultiVector<RealType>,aol::MultiVector<RealType> > minimizer( e, de, 200, 1.e-5, 50, false );
  minimizer.setTauMin( .001 ); // in case the Armijo condition is used
  minimizer.setBeta( .9 );
  minimizer.setTimestepController( aol::NewtonInfo<RealType>::WOLFE );*/
  /*// gradient descent
  aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > minimizer( Grid, e, de, 50, 1., 1.e-5 );
  minimizer.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::USE_NONLINEAR_CG | aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::DO_NOT_SMOOTH_DESCENT_DIRECTION );*/
  // trust region method
  aol::TrustRegionMethod<RealType,aol::MultiVector<RealType>,aol::BlockOpBase<RealType,SubMatrixType>,aol::DiagonalPreconditioner<aol::MultiVector<RealType> >,SubMatrixType> minimizer( e, de, d2e, new aol::BlockMatrix<SubMatrixType>( corrections.numComponents(), corrections.numComponents(), Grid.getNumberOfNodes(), Grid.getNumberOfNodes() ), steps, 1.e-5 );
  minimizer.setIterativeMode();

  minimizer.apply( zeros, corrections );

  // compute the corrected shapes
  for ( int k = 1; k < Displacements.size(); k++ ) {
    aol::MultiVector<RealType> correction;
    aol::Vector<RealType> deformedImage( Grid );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      correction.appendReference( corrections[ConfiguratorType::Dim*(k-1)+j] );
    qc::InvDeformImage<ConfiguratorType, aol::Vector<RealType> > ( Shapes[k], Grid, deformedImage, correction );
    Shapes[k] = deformedImage;
  }

  // compute the corrected displacements
  for ( int k = 0; k < Displacements.size(); k++ ) {
    aol::MultiVector<RealType> correctionL, correctionR, auxDisp( Grid );
    if ( k > 0 )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        correctionR.appendReference( corrections[ConfiguratorType::Dim*(k-1)+j] );
    else
      correctionR.reallocate( Grid );
    if ( k < Displacements.size() - 1 )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        correctionL.appendReference( corrections[ConfiguratorType::Dim*k+j] );
    else
      correctionL.reallocate( Grid );

    qc::ConcatAndSmoothlyExtendDeformations<ConfiguratorType>( correctionL, Displacements[k], Grid, auxDisp );
    qc::InvConcatAndSmoothlyExtendDeformations<ConfiguratorType>( auxDisp, correctionR, Grid, Displacements[k] );
  }
}

template <typename ConfiguratorType, typename MaterialLawType>
void correctDeformationsAndShapes ( const typename ConfiguratorType::InitType &Grid,
                                    aol::MultiVector<typename ConfiguratorType::RealType> &Shapes,
                                    aol::MultiVector<typename ConfiguratorType::RealType>  &Displacements,
                                    MaterialLawType &ElasticEnergyDensity, 
                                    int steps = 50 ) {
  aol::RandomAccessContainer< aol::MultiVector<typename ConfiguratorType::RealType> > displacements;
  int N = Displacements.numComponents() / ConfiguratorType::Dim; 
  // write to RandomAccessContainer
  for ( int k = 0; k <  N; k++ ){
    aol::MultiVector<typename ConfiguratorType::RealType> aux;
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      aux.appendReference( Displacements[ k * ConfiguratorType::Dim + j ] );
    displacements.pushBack( aux );
  }

  correctDeformationsAndShapes<ConfiguratorType, MaterialLawType> ( Grid, Shapes, displacements, ElasticEnergyDensity, steps );

  for ( int k = 0; k <  N; k++ )
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      Displacements[k *ConfiguratorType::Dim + j ] = displacements[k][j];

}

#endif
