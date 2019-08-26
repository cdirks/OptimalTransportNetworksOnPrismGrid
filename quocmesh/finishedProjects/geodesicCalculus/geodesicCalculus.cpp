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

#include "geodesicCalculusEnergies.h"
#include <trustRegionMethod.h>
#include <smoother.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gradientflow.h>

template <typename ArrayType>
void loadCharFun( const char* FileName, ArrayType &CharFun ) {
  // load characteristic function
  cerr<<"load " + string( FileName ) + "...\n";
  CharFun.load( FileName );
  // intensify the contrast
  CharFun.addToAll( -CharFun.getMinValue() );
  CharFun /= CharFun.getMaxValue();
}

template <typename ConfiguratorType>
void saveMultiArray( const aol::MultiVector<typename ConfiguratorType::RealType> &Arg,
                     const typename ConfiguratorType::InitType &Grid,
                     const char* FileNameTemplate,
                     const char* Directory,
                     const int Number = -1 ) {
  for ( int i = 0; i < Arg.numComponents(); i++ ) {
    char fileName[1024];
    if ( Number < 0 )
      sprintf( fileName, FileNameTemplate, Directory, i );
    else
      sprintf( fileName, FileNameTemplate, Directory, Number, i );
    typename ConfiguratorType::ArrayType( Arg[i], Grid, aol::FLAT_COPY ).save( fileName, qc::PGM_DOUBLE_BINARY );
  }
}

template <typename RealType>
class MaskingOp :
  public aol::Op<aol::MultiVector<RealType> > {

private:
  const aol::BitVector &_mask;
  
public:
  MaskingOp( const aol::BitVector &Mask ) :
    _mask( Mask ) {}
    
  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    for ( int i = 0; i < Arg.numComponents(); i++ )
      for ( int j = 0; j < Arg[i].size(); j++ )
        if ( !_mask[j] )
          Dest[i][j] += Arg[i][j];
  }
};

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType, typename SubMatrixType = typename ConfiguratorType::MatrixType>
class MaskedHyperelasticHessian
      : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::BlockOpBase<typename ConfiguratorType::RealType,SubMatrixType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the underlying FE grid
  const typename ConfiguratorType::InitType &_grid;
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the Dirichlet nodes
  const aol::BitVector &_mask;

public:
  MaskedHyperelasticHessian ( const typename ConfiguratorType::InitType &Grid,
                              const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                              const aol::BitVector &Mask )
      : _grid ( Grid ),
        _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
        _mask( Mask ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::BlockOpBase<RealType,SubMatrixType> &Dest ) const {
    for ( int k = 0; k < ConfiguratorType::Dim; k++ )
      for ( int l = 0; l < ConfiguratorType::Dim; l++ )
        qc::HyperelasticSubHessian<ConfiguratorType,HyperelasticEnergyDensityType>( _grid, _hyperelasticEnergyDensity, Arg, k, l ).assembleAddMatrix( Dest.getReference( k, l ), &_mask, k == l );
  }
};

template <typename ConfiguratorType, typename MaterialLawType>
void elasticallySmoothDeformation ( const aol::MultiVector<typename ConfiguratorType::RealType> &Displacement,
                                    const typename ConfiguratorType::InitType &Grid,
                                    const MaterialLawType &MaterialLaw,
                                    aol::MultiVector<typename ConfiguratorType::RealType> &ResultDisplacement,
                                    const typename qc::BitArrayTrait<ConfiguratorType::Dim>::ArrayType &ValueFixed,
                                    const bool QuietMode = true ) {
  typedef typename ConfiguratorType::MatrixType MatrixType;
  typedef typename ConfiguratorType::RealType RealType;

  if ( !QuietMode )
    cerr<<"extend deformation elastically optimally\n";

  qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>     e( Grid, MaterialLaw );
  qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>  de( Grid, MaterialLaw );
  MaskingOp<RealType> maskingOp( ValueFixed );
  aol::CompositeOp<aol::MultiVector<RealType> > mde;
  mde.appendReference( de );
  mde.appendReference( maskingOp );
  MaskedHyperelasticHessian<ConfiguratorType,MaterialLawType> d2e( Grid, MaterialLaw, ValueFixed );

  aol::TrustRegionMethod<RealType, aol::MultiVector<RealType>, aol::BlockOpBase<RealType,MatrixType>, aol::DiagonalPreconditioner<aol::MultiVector<RealType> >, MatrixType> minimizer( e, mde, d2e, new aol::BlockMatrix<MatrixType>( Grid ), 2000, 1.e-6 );
  minimizer.setIterativeMode( false );
  minimizer.setQuietMode( QuietMode );
  minimizer.apply( Displacement, ResultDisplacement );
}

template <typename ConfiguratorType, typename MaterialLawType>
void invConcatAndSmoothlyExtendDeformations ( const aol::MultiVector<typename ConfiguratorType::RealType> &Displacement1,
                                              const aol::MultiVector<typename ConfiguratorType::RealType> &Displacement2,
                                              const typename ConfiguratorType::InitType &Grid,
                                              const MaterialLawType &MaterialLaw,
                                              aol::MultiVector<typename ConfiguratorType::RealType> &ResultDisplacement,
                                              const bool QuietMode = true ) {
  typedef typename ConfiguratorType::MatrixType MatrixType;
  typedef typename ConfiguratorType::RealType RealType;

  // compute the deformation \phi_1
  aol::MultiVector<RealType> identity( Grid ), tmp( Displacement1 );
  qc::DataGenerator<ConfiguratorType>( Grid ).generateIdentity( identity );
  tmp += identity;

  // concatenate \phi_1 with \phi_2^{-1}
  typename qc::BitArray<ConfiguratorType::Dim> invDispDefined;
  qc::TransformFunction<RealType,ConfiguratorType::Dim> transformFunction( Grid );
  transformFunction.setDeformation( Displacement2 );
  transformFunction.transform( tmp, ResultDisplacement, invDispDefined );

  // compute the corresponding displacement
  ResultDisplacement -= identity;

  // harmonically extend the displacement where it is not yet defined
  qc::SmoothlyExtendImage<ConfiguratorType,aol::MultiVector<RealType> >( ResultDisplacement, Grid, tmp, invDispDefined );

  // resolve any potential self-interpenetrations
  resolveSelfPenetration<ConfiguratorType>( tmp, Grid, 100, 0, QuietMode );

  // turn harmonic extension into elastically optimal extension
  elasticallySmoothDeformation<ConfiguratorType,MaterialLawType>( tmp, Grid, MaterialLaw, ResultDisplacement, invDispDefined, QuietMode );
}

template <typename ConfiguratorType, typename MaterialLawType>
void elasticallySmoothDeformationAtBoundary ( const aol::MultiVector<typename ConfiguratorType::RealType> &Displacement,
                                              const typename ConfiguratorType::InitType &Grid,
                                              const MaterialLawType &MaterialLaw,
                                              aol::MultiVector<typename ConfiguratorType::RealType> &ResultDisplacement,
                                              const int NumBoundaryPixels = 4,
                                              const bool QuietMode = true ) {
  typedef typename ConfiguratorType::MatrixType MatrixType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType VecType;

  if ( !QuietMode )
    cerr << "elastically smooth displacement in " << NumBoundaryPixels << " pixels wide boundary region ...\n";
  typename qc::BitArrayTrait<ConfiguratorType::Dim>::ArrayType valueFixed( Grid );
  aol::Vector<RealType> aux( Grid );
  qc::DataGenerator<ConfiguratorType>( Grid ).generateScaledLInfinitySphereLevelset( aux, VecType( .5 ), VecType( 1. / Grid.H() ), 0.5 / Grid.H() - NumBoundaryPixels );
  aux.clamp( -1., 0. );
  static_cast<aol::BitVector*>( &valueFixed )->setNonzeroPatternFrom<RealType,aol::Vector<RealType> >( aux );
  elasticallySmoothDeformation<ConfiguratorType,MaterialLawType>( Displacement, Grid, MaterialLaw, ResultDisplacement, valueFixed, QuietMode );
}

//! Energy as described in Rumpf and Wirth, "DISCRETE GEODESIC CALCULUS IN THE SPACE OF VISCOUS FLUIDIC OBJECTS", chapter 4.1
//! \author Wirth
//! In particular, see figure 4.1 for reference. The (\phi_k)_k are the genuine degrees of freedom, denoted as "correction displacements" here.
//! The reference displacements (\bar \psi_k)_k are assumed to be given and denoted by "morphing displacements" here.
//! Note that the main advantage of this approach is having only K deformations as DOFs instead of K deformations AND K+1 shapes!
//! The drawback is having concatenated arguments and the assumption of suitbale reference shapes and reference deformations to be given.
template <typename ConfiguratorType, typename MaterialLawType>
class MorphingCorrectionEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,MaterialLawType> WeightElastDens;

  const typename ConfiguratorType::InitType &_grid;
  const int _numTimePoints;
  const aol::RandomAccessContainer< qc::MultiArray<RealType, ConfiguratorType::Dim> > &_morphingDisplacements;
  const aol::Vector<RealType> &_initialImage, &_finalImage, &_initialImageApprox, &_finalImageApprox;
  const MaterialLawType _elasticEnergyDensity;
  const aol::MultiVector<RealType> &_elasticEnergyDensityWeight;
  const RealType _matchingWeight;
  const bool _withMatching;
  
public:
  // note: the first MorphingDisplacement is a dummy
  MorphingCorrectionEnergy( const typename ConfiguratorType::InitType &Grid,
                            const aol::RandomAccessContainer< qc::MultiArray<RealType, ConfiguratorType::Dim> > &MorphingDisplacements,
                            const aol::Vector<RealType> &InitialImage,
                            const aol::Vector<RealType> &FinalImage,
                            const aol::Vector<RealType> &InitialImageApprox,
                            const aol::Vector<RealType> &FinalImageApprox,
                            const MaterialLawType &ElasticEnergyDensity,
                            const aol::MultiVector<RealType> &ElasticEnergyDensityWeight,
                            const RealType MatchingWeight,
                            const bool WithMatching = true ) :
    _grid( Grid ),
    _numTimePoints( MorphingDisplacements.size() ),
    _morphingDisplacements( MorphingDisplacements ),
    _initialImage( InitialImage ),
    _finalImage( FinalImage ),
    _initialImageApprox( InitialImageApprox ),
    _finalImageApprox( FinalImageApprox ),
    _elasticEnergyDensity( ElasticEnergyDensity ),
    _elasticEnergyDensityWeight( ElasticEnergyDensityWeight ),
    _matchingWeight( MatchingWeight ),
    _withMatching( WithMatching ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    // bring argument into reasonable format
    aol::RandomAccessContainer<aol::MultiVector<RealType> > arg( _numTimePoints );
    for ( int k = 0; k < _numTimePoints; k++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        arg[k].appendReference( Arg[k*ConfiguratorType::Dim+j] );
    // elastic deformation energy
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 1; k < _numTimePoints; k++ ) {
      aol::Scalar<RealType> dest;
      aol::MultiVector<RealType> correction;
      correction.appendReference( arg[k] );
      correction.appendReference( arg[k-1] );
      WeightElastDens elastDens( _elasticEnergyDensity, _grid, _elasticEnergyDensityWeight[k-1] );
      qc::HyperelasticPrePostDeformEnergy<ConfiguratorType,WeightElastDens>( _grid, elastDens, _morphingDisplacements[k] ).applyAdd( correction, dest );
#ifdef _OPENMP
#pragma omp critical (MorphingCorrectionEnergyApplyAdd)
#endif
      Dest += dest;
    }
    // small regularization of all displacements
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 0; k < _numTimePoints; k++ ) {
      aol::Scalar<RealType> dest;
      qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity ).applyAdd( arg[_numTimePoints-1], dest );
#ifdef _OPENMP
#pragma omp critical (MorphingCorrectionEnergyApplyAdd)
#endif
      Dest.addMultiple( dest, 1.e-3 );
    }
    // matching energy
    if ( _withMatching ) {
      aol::Scalar<RealType> dest;
      qc::MismatchEnergyWRTDisp<ConfiguratorType>( _grid, _initialImageApprox, _initialImage ).apply( arg[0], dest );
      qc::MismatchEnergyWRTDisp<ConfiguratorType>( _grid, _finalImageApprox, _finalImage ).applyAdd( arg[_numTimePoints-1], dest );
      Dest.addMultiple( dest, _matchingWeight );
    }
    // to make step size control work, replace NaN by Inf
    if ( aol::isNaN( Dest.v ) )
      Dest = aol::NumberTrait<RealType>::Inf;
  }
};

template <typename ConfiguratorType, typename MaterialLawType>
class MorphingCorrectionGradient :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,MaterialLawType> WeightElastDens;

  const typename ConfiguratorType::InitType &_grid;
  const int _numTimePoints;
  const aol::RandomAccessContainer< qc::MultiArray<RealType, ConfiguratorType::Dim> > &_morphingDisplacements;
  const aol::Vector<RealType> &_initialImage, &_finalImage, &_initialImageApprox, &_finalImageApprox;
  const MaterialLawType _elasticEnergyDensity;
  const aol::MultiVector<RealType> &_elasticEnergyDensityWeight;
  const RealType _matchingWeight;
  const bool _withMatching;

public:
  // note: the first MorphingDisplacement is a dummy
  MorphingCorrectionGradient( const typename ConfiguratorType::InitType &Grid,
                              const aol::RandomAccessContainer< qc::MultiArray<RealType, ConfiguratorType::Dim> > &MorphingDisplacements,
                              const aol::Vector<RealType> &InitialImage,
                              const aol::Vector<RealType> &FinalImage,
                              const aol::Vector<RealType> &InitialImageApprox,
                              const aol::Vector<RealType> &FinalImageApprox,
                              const MaterialLawType &ElasticEnergyDensity,
                              const aol::MultiVector<RealType> &ElasticEnergyDensityWeight,
                              const RealType MatchingWeight,
                              const bool WithMatching = true ) :
    _grid( Grid ),
    _numTimePoints( MorphingDisplacements.size() ),
    _morphingDisplacements( MorphingDisplacements ),
    _initialImage( InitialImage ),
    _finalImage( FinalImage ),
    _initialImageApprox( InitialImageApprox ),
    _finalImageApprox( FinalImageApprox ),
    _elasticEnergyDensity( ElasticEnergyDensity ),
    _elasticEnergyDensityWeight( ElasticEnergyDensityWeight ),
    _matchingWeight( MatchingWeight ),
    _withMatching( WithMatching ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // bring argument and result into reasonable format
    aol::RandomAccessContainer<aol::MultiVector<RealType> > arg( _numTimePoints ), dest( _numTimePoints );
    for ( int k = 0; k < _numTimePoints; k++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        arg[k].appendReference( Arg[k*ConfiguratorType::Dim+j] );
        dest[k].appendReference( Dest[k*ConfiguratorType::Dim+j] );
      }
    // elastic deformation energy
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 1 + ( _withMatching ? 0 : 1 ); k < _numTimePoints; k++ ) {
      aol::MultiVector<RealType> correction;
      correction.appendReference( arg[k] );
      correction.appendReference( arg[k-1] );
      WeightElastDens elastDens( _elasticEnergyDensity, _grid, _elasticEnergyDensityWeight[k-1] );
      qc::HyperelasticPrePostDeformGradientPreComp<ConfiguratorType,WeightElastDens>( _grid, elastDens, _morphingDisplacements[k] ).applyAdd( correction, dest[k-1] );
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 1; k < _numTimePoints - ( _withMatching ? 0 : 1 ); k++ ) {
      aol::MultiVector<RealType> correction;
      correction.appendReference( arg[k] );
      correction.appendReference( arg[k-1] );
      WeightElastDens elastDens( _elasticEnergyDensity, _grid, _elasticEnergyDensityWeight[k-1] );
      qc::HyperelasticPrePostDeformGradientPostComp<ConfiguratorType,WeightElastDens>( _grid, elastDens, _morphingDisplacements[k] ).applyAdd( correction, dest[k] );
    }
    // small regularization of all displacements
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = ( _withMatching ? 0 : 1 ); k < _numTimePoints - ( _withMatching ? 0 : 1 ); k++ ) {
      aol::MultiVector<RealType> destComp( dest[0], aol::STRUCT_COPY );
      qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity ).applyAdd( arg[_numTimePoints-1], destComp );
      dest[k].addMultiple( destComp, 1.e-3 );
    }
    // matching energy
    if ( _withMatching ) {
      aol::MultiVector<RealType> destComp( dest[0], aol::STRUCT_COPY );
      qc::MismatchEnergyVariationWRTDisp<ConfiguratorType>( _grid, _initialImageApprox, _initialImage ).apply( arg[0], destComp );
      dest[0].addMultiple( destComp, _matchingWeight );
      qc::MismatchEnergyVariationWRTDisp<ConfiguratorType>( _grid, _finalImageApprox, _finalImage ).apply( arg[_numTimePoints-1], destComp );
      dest[_numTimePoints-1].addMultiple( destComp, _matchingWeight );
    }
  }
};

template <typename ConfiguratorType, typename MaterialLawType, typename SubMatrixType = typename ConfiguratorType::MatrixType>
class MorphingCorrectionHessian :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::BlockOpBase<typename ConfiguratorType::RealType,SubMatrixType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,MaterialLawType> WeightElastDens;

  const typename ConfiguratorType::InitType &_grid;
  const int _numTimePoints;
  const aol::RandomAccessContainer< qc::MultiArray<RealType, ConfiguratorType::Dim> > &_morphingDisplacements;
  const aol::Vector<RealType> &_initialImage, &_finalImage, &_initialImageApprox, &_finalImageApprox;
  const aol::MultiVector<RealType> &_gradInitialImage, &_gradFinalImage;
  const MaterialLawType _elasticEnergyDensity;
  const aol::MultiVector<RealType> &_elasticEnergyDensityWeight;
  const RealType _matchingWeight;
  const bool _withMatching;

public:
  // note: the first MorphingDisplacement is a dummy
  MorphingCorrectionHessian( const typename ConfiguratorType::InitType &Grid,
                             const aol::RandomAccessContainer< qc::MultiArray<RealType, ConfiguratorType::Dim> > &MorphingDisplacements,
                             const aol::Vector<RealType> &InitialImage,
                             const aol::Vector<RealType> &FinalImage,
                             const aol::Vector<RealType> &InitialImageApprox,
                             const aol::Vector<RealType> &FinalImageApprox,
                             const aol::MultiVector<RealType> &GradInitialImage,
                             const aol::MultiVector<RealType> &GradFinalImage,
                             const MaterialLawType &ElasticEnergyDensity,
                             const aol::MultiVector<RealType> &ElasticEnergyDensityWeight,
                             const RealType MatchingWeight,
                             const bool WithMatching = true ) :
    _grid( Grid ),
    _numTimePoints( MorphingDisplacements.size() ),
    _morphingDisplacements( MorphingDisplacements ),
    _initialImage( InitialImage ),
    _finalImage( FinalImage ),
    _initialImageApprox( InitialImageApprox ),
    _finalImageApprox( FinalImageApprox ),
    _gradInitialImage( GradInitialImage ),
    _gradFinalImage( GradFinalImage ),
    _elasticEnergyDensity( ElasticEnergyDensity ),
    _elasticEnergyDensityWeight( ElasticEnergyDensityWeight ),
    _matchingWeight( MatchingWeight ),
    _withMatching( WithMatching ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::BlockOpBase<RealType,SubMatrixType> &Dest ) const {
    // bring argument into reasonable format
    aol::RandomAccessContainer<aol::MultiVector<RealType> > arg( _numTimePoints );
    for ( int k = 0; k < _numTimePoints; k++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        arg[k].appendReference( Arg[k*ConfiguratorType::Dim+j] );
    // elastic deformation energy
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 1 + ( _withMatching ? 0 : 1 ); k < _numTimePoints; k++ )
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          WeightElastDens elastDens( _elasticEnergyDensity, _grid, _elasticEnergyDensityWeight[k-1] );
          qc::HyperelasticPrePostDeformSubHessianPreComp<ConfiguratorType,WeightElastDens>( _grid, elastDens, _morphingDisplacements[k], arg[k-1], arg[k], i, j ).assembleAddMatrix( Dest.getReference( (k-1)*ConfiguratorType::Dim+i, (k-1)*ConfiguratorType::Dim+j ) );
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 1; k < _numTimePoints - ( _withMatching ? 0 : 1 ); k++ )
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          WeightElastDens elastDens( _elasticEnergyDensity, _grid, _elasticEnergyDensityWeight[k-1] );
          qc::HyperelasticPrePostDeformSubHessianPostComp<ConfiguratorType,WeightElastDens>( _grid, elastDens, _morphingDisplacements[k], arg[k-1], arg[k], i, j ).assembleAddMatrix( Dest.getReference( k*ConfiguratorType::Dim+i, k*ConfiguratorType::Dim+j ) );
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 1 + ( _withMatching ? 0 : 1 ); k < _numTimePoints - ( _withMatching ? 0 : 1 ); k++ )
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          WeightElastDens elastDens( _elasticEnergyDensity, _grid, _elasticEnergyDensityWeight[k-1] );
          qc::HyperelasticPrePostDeformSubHessianMixedPart<ConfiguratorType,WeightElastDens>( _grid, elastDens, _morphingDisplacements[k], arg[k-1], arg[k], i, j, true ).assembleAddMatrix( Dest.getReference( k*ConfiguratorType::Dim+i, (k-1)*ConfiguratorType::Dim+j ) );
          qc::HyperelasticPrePostDeformSubHessianMixedPart<ConfiguratorType,WeightElastDens>( _grid, elastDens, _morphingDisplacements[k], arg[k-1], arg[k], i, j, false ).assembleAddMatrix( Dest.getReference( (k-1)*ConfiguratorType::Dim+i, k*ConfiguratorType::Dim+j ) );
        }
    // small regularization of all displacements
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = ( _withMatching ? 0 : 1 ); k < _numTimePoints - ( _withMatching ? 0 : 1 ); k++ )
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          qc::HyperelasticSubHessian<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity, arg[k], i, j ).assembleAddMatrix( Dest.getReference( k*ConfiguratorType::Dim+i, k*ConfiguratorType::Dim+j ), 1.e-3 );
    // matching energy
    if ( _withMatching ) {
      aol::BlockOp<RealType,SubMatrixType> destComp( ConfiguratorType::Dim, ConfiguratorType::Dim );
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          destComp.setReference( i, j, Dest.getReference( i, j ) );
      qc::MismatchEnergySecondVariationWRTDisp<ConfiguratorType,aol::BlockOp<RealType,SubMatrixType> >
        ( _grid, _initialImageApprox, _initialImage, _gradInitialImage ).applyAddMultiple( arg[0], destComp, _matchingWeight );
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          destComp.setReference( i, j, Dest.getReference( (_numTimePoints-1)*ConfiguratorType::Dim+i, (_numTimePoints-1)*ConfiguratorType::Dim+j ) );
      qc::MismatchEnergySecondVariationWRTDisp<ConfiguratorType,aol::BlockOp<RealType,SubMatrixType> >
        ( _grid, _finalImageApprox, _finalImage, _gradFinalImage ).applyAddMultiple( arg[_numTimePoints-1], destComp, _matchingWeight );
    } else
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int row = 0; row < Dest.getReference( i, i ).getNumRows(); row++ ) {
          Dest.getReference( i, i ).add( row, row, 1. );
          Dest.getReference( (_numTimePoints-1)*ConfiguratorType::Dim+i, (_numTimePoints-1)*ConfiguratorType::Dim+i ).add( row, row, 1. );
        }
  }
};

template<typename RealType, typename VectorType, typename MatrixType, bool DirectSolver = false, typename PreconditionType = aol::IdentityOp<VectorType>, typename SubMatrixType = MatrixType >
class NewtonMinimizer :
  public aol::NewtonMinimizationBase< RealType, VectorType, VectorType, MatrixType > {

  mutable PreconditionType *_pPrecond;

public:
  NewtonMinimizer( const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                   const aol::Op<VectorType> &DE,
                   const aol::Op<VectorType,MatrixType> &D2E,
                   MatrixType *PMatD2E,
                   const int MaxIterations = 50,
                   const RealType StopEpsilon = 1.e-6,
                   const bool WriteTimeSteps = false,
                   const char *BaseSaveName = NULL ) :
    aol::NewtonMinimizationBase<RealType,VectorType,VectorType,MatrixType>( E, DE, D2E, PMatD2E, MaxIterations, StopEpsilon, WriteTimeSteps, BaseSaveName ),
    _pPrecond( NULL ) {}

  ~NewtonMinimizer() {
    delete _pPrecond;
  }

  virtual void prepareSolver() const {
    if ( DirectSolver ) {
      if ( this->_pSolver == NULL )
        this->_pSolver = new aol::CholeskyBlockInverseOp<RealType,SubMatrixType>( 243 );
      static_cast<aol::CholeskyBlockInverseOp<RealType,SubMatrixType>*>( this->_pSolver )->setMatrix( *this->_pMatDF );
    } else {
      delete this->_pSolver;
      delete _pPrecond;
      _pPrecond = new PreconditionType( *this->_pMatDF );
      this->_pSolver = new aol::PCGInverse<VectorType>( *this->_pMatDF, *_pPrecond, this->_pInfo->getSolverInfo() );
      // the solver-accuracy is automatically set to fit to the Newton-minimizer accuracy
      static_cast<aol::PCGInverse<VectorType>*>( this->_pSolver )->setStopping( aol::STOPPING_ABSOLUTE );
    }
  }
};

template <typename ConfiguratorType, typename VectorType, typename SecondDerivativeType, typename SubMatrixType>
class NewtonIteration :
  public aol::NewtonIterationBase<typename ConfiguratorType::RealType,VectorType,VectorType,SecondDerivativeType> {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;

  void prepareSolver () const {
    // aol::CholeskyBiBlockInverseOp<ConfiguratorType,SubMatrixType,ConfiguratorType::Dim>* pSolver = new aol::CholeskyBiBlockInverseOp<ConfiguratorType,SubMatrixType,ConfiguratorType::Dim>( _grid, *(this->_pMatDF) );
    aol::UMFPACKBlockInverseOp<RealType,SubMatrixType>* pSolver = new aol::UMFPACKBlockInverseOp<RealType,SubMatrixType>( *(this->_pMatDF) );
    this->_pSolver = pSolver;
  }
  
public:
  NewtonIteration( const typename ConfiguratorType::InitType &Grid,
                   const aol::Op<VectorType> &F,
                   const aol::Op<VectorType, SecondDerivativeType> &DF,
                   const int MaxIterations = 50,
                   const RealType StopEpsilon = 1.e-6,
                   const bool WriteTimeSteps = false,
                   const char *BaseSaveName = NULL) :
    aol::NewtonIterationBase<RealType,VectorType,VectorType,SecondDerivativeType>( new SecondDerivativeType( ConfiguratorType::Dim, ConfiguratorType::Dim ), F, DF, MaxIterations, StopEpsilon, WriteTimeSteps, BaseSaveName ),
    _grid( Grid ) {
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        this->_pMatDF->allocateMatrix ( i, j, Grid );
  }
};

/** Computes the energy \f$ \int_{\phi_l(\Omega)} W(D\phi) dx + \int_{\phi(\phi_l(\Omega))} W(D(\phi_r\circ\phi_l^{-1}\circ\phi^{-1})) dx \f$,
 *  where \f$ \phi-id \f$ is passed to apply and \f$ \phi_l-id,\phi_r-id \f$ are passed to the constructor as "DisplacementL" and "DisplacementR".
 *  The computation is done on \f$ \Omega \f$ (after a change of variables), which is passed to the constructor as "Grid".
 */
template <typename ConfiguratorType, typename MaterialLawType>
class MidPointEnergy :
  public qc::FENonlinDeformIntegrationVectorInterface < ConfiguratorType, MidPointEnergy<ConfiguratorType,MaterialLawType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_l and \phi_r
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementL, _displacementR;
  
public:
  MidPointEnergy ( const typename ConfiguratorType::InitType &Grid,
                   const MaterialLawType &HyperelasticEnergyDensity,
                   const aol::MultiVector<RealType> &DisplacementL,
                   const aol::MultiVector<RealType> &DisplacementR ) :
    qc::FENonlinDeformIntegrationVectorInterface < ConfiguratorType, MidPointEnergy<ConfiguratorType,MaterialLawType> > ( Grid, DisplacementL ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementL( Grid, DisplacementL ),
    _displacementR( Grid, DisplacementR ) {}

  RealType evaluateIntegrand ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord ) const {
    // compute the deformation gradients $\nabla\phi$, $\nabla\phi_l$, $\nabla\phi_r$
    typename ConfiguratorType::MatType dPhi,dPhiL,dPhiR;
    DiscFuncs.evaluateGradient( TransformedEl, TransformedLocalCoord, dPhi );
    _displacementL.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiL );
    _displacementR.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiR );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      dPhi[i][i]  += 1.;
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute the deformation energy
    dPhiR *= dPhiL.inverse();
    dPhiR *= dPhi.inverse();
    return dPhiL.det() * ( qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhi, _hyperelasticEnergyDensity, El, QuadPoint )
         + dPhi.det() * qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhiR, _hyperelasticEnergyDensity, El, QuadPoint ) );
  }
};

//! Computes the Gateaux derivative of MidPointEnergy wrt \f$ \phi \f$.
template <typename ConfiguratorType, typename MaterialLawType>
class MidPointDeriv :
  public qc::FENonlinDeformVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, MidPointDeriv<ConfiguratorType,MaterialLawType>, false> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_l and \phi_r
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementL, _displacementR;

public:
  MidPointDeriv ( const typename ConfiguratorType::InitType &Grid,
                  const MaterialLawType &HyperelasticEnergyDensity,
                  const aol::MultiVector<RealType> &DisplacementL,
                  const aol::MultiVector<RealType> &DisplacementR ) :
    qc::FENonlinDeformVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, MidPointDeriv<ConfiguratorType,MaterialLawType>, false> ( Grid, DisplacementL ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementL( Grid, DisplacementL ),
    _displacementR( Grid, DisplacementR ) {}

  void getNonlinearity ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                         const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                         aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, RealType> &NL ) const {
    // compute the deformation gradients $\nabla\phi$, $\nabla\phi_l$, $\nabla\phi_r$
    typename ConfiguratorType::MatType dPhi,dPhiL,dPhiR;
    DiscFuncs.evaluateGradient( TransformedEl, TransformedLocalCoord, dPhi );
    _displacementL.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiL );
    _displacementR.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiR );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      dPhi[i][i]  += 1.;
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute \partial_{,A}W(\nabla\phi), \partial_{,A}W(\nabla\phi_r\nabla\phi_l^{-1}\nabla\phi^{-1})
    typename ConfiguratorType::MatType cofDPhi, dPhiRL( dPhiR ), dWdPhi, dWdPhiRL;
    dPhiRL *= dPhiL.inverse();
    dPhiRL *= dPhi.inverse();
    cofDPhi.makeCofactorMatrix( dPhi );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhi, _hyperelasticEnergyDensity, El, QuadPoint, dWdPhi );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiRL, _hyperelasticEnergyDensity, El, QuadPoint, dWdPhiRL );

    // compute NL = det(\nabla\phi_l)\partial_{,A}W(\nabla\phi)+det(\nabla\phi_l)W(\nabla\phi_r\nabla\phi_l^{-1}\nabla\phi^{-1})cof(\nabla\phi)-...
    // cof(\nabla\phi)cof(\nabla\phi_l)\nabla\phi_r^T\partial_{,A}W(\nabla\phi_r\nabla\phi_l^{-1}\nabla\phi^{-1})\nabla\phi^{-T}
    NL.makeProductAtransposedB( dPhiRL, dWdPhiRL );
    NL *= cofDPhi;
    NL.addMultiple( cofDPhi, -qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhiRL, _hyperelasticEnergyDensity, El, QuadPoint ) );
    NL -= dWdPhi;
    NL *= -dPhiL.det();
  }
};

//! Computes (k,l)-block of the Hessian of MidPointEnergy wrt \f$ \phi \f$.
template <typename ConfiguratorType, typename MaterialLawType>
class MidPointSubHessian :
  public qc::FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType,MidPointSubHessian<ConfiguratorType,MaterialLawType>,false> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_l and \phi_r and \phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementL, _displacementR, _displacement;
  // the row and column index of the sub-Hessian
  const int _k, _l;
  
public:
  MidPointSubHessian ( const typename ConfiguratorType::InitType &Grid,
                       const MaterialLawType &HyperelasticEnergyDensity,
                       const aol::MultiVector<RealType> &DisplacementL,
                       const aol::MultiVector<RealType> &DisplacementR,
                       const aol::MultiVector<RealType> &Displacement,
                       const int K,
                       const int L ) :
    qc::FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType,MidPointSubHessian<ConfiguratorType,MaterialLawType>,false> ( Grid, DisplacementL ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementL( Grid, DisplacementL ),
    _displacementR( Grid, DisplacementR ),
    _displacement( Grid, Displacement ),
    _k( K ),
    _l( L ) {}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*LocalCoord*/,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                               typename ConfiguratorType::MatType &Matrix ) const {
    // compute the deformation gradients $\nabla\phi$, $\nabla\phi_l$, $\nabla\phi_r$
    typename ConfiguratorType::MatType dPhi,dPhiL,dPhiR;
    _displacement.evaluateGradient( TransformedEl, TransformedLocalCoord, dPhi );
    _displacementL.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiL );
    _displacementR.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiR );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      dPhi[i][i]  += 1.;
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute D\phi_r D\phi_l^{-1} D\phi^{-1} and W, \partial_{,A}W
    typename ConfiguratorType::MatType dPhiRInvDPhiLInvDPhi( dPhiR ), dPhiInv( dPhi.inverse() ), dW;
    dPhiRInvDPhiLInvDPhi *= dPhiL.inverse();
    dPhiRInvDPhiLInvDPhi *= dPhiInv;
    RealType W = qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhiRInvDPhiLInvDPhi, _hyperelasticEnergyDensity, El, QuadPoint );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiRInvDPhiLInvDPhi, _hyperelasticEnergyDensity, El, QuadPoint, dW );

    // compute the Hessian matrix; denote variations of \phi by \psi and \chi; they only act on the component "_k" and "_l"
    // \psi shall belong to vectors multiplied from the left, \chi from the right
    typename ConfiguratorType::MatType cofDPhi, aux, aux2;
    typename ConfiguratorType::VecType auxVec1, auxVec2;
    cofDPhi.makeCofactorMatrix( dPhi );
    
    qc::HyperelasticSubHessian<ConfiguratorType,MaterialLawType>::elasticityTensor( dPhi, _hyperelasticEnergyDensity, El, QuadPoint, _k, _l, Matrix );

    if ( _k != _l ) {
      aux.makeProductABtransposed( dPhiInv, dW );
      aux *= dPhiRInvDPhiLInvDPhi;
      aux.getColumn( _l, auxVec1 );
      aux.getColumn( _k, auxVec2 );
      aux.makeTensorProduct( cofDPhi[_k], auxVec1 );
      Matrix += aux;

      aux.transpose();
      Matrix -= aux;

      aux.makeTensorProduct( cofDPhi[_l], auxVec2 );
      Matrix -= aux;

      aux.transpose();
      Matrix += aux;

      // M s.t. D\psi_k M D\chi_l = W(D\phi_r D\phi_l^{-1} D\phi^{-1}) \partial_{,A}^2 det(D\phi)(D\psi,D\chi)
      aux.makeTensorProduct( cofDPhi[_k], cofDPhi[_l] );
      Matrix.addMultiple( aux, W / dPhi.det() );
      aux.transpose();
      Matrix.addMultiple( aux, -W / dPhi.det() );
    }

    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        qc::HyperelasticSubHessian<ConfiguratorType,MaterialLawType>::elasticityTensor( dPhiRInvDPhiLInvDPhi, _hyperelasticEnergyDensity, El, QuadPoint, i, j, aux );
        aux *= cofDPhi;
        aux2.makeProduct( dPhiInv, aux );
        Matrix.addMultiple( aux2, dPhiRInvDPhiLInvDPhi[i][_k] * dPhiRInvDPhiLInvDPhi[j][_l] );
      }
    
    Matrix *= dPhiL.det();
  }
};

/** Computes the energy \f$ \int_{\phi_l(\Omega)} W(D(\phi\circ\phi_l^{-1})) dx + \int_{\phi(\Omega)} W(D(\phi_r\circ\phi^{-1})) dx \f$,
 *  where \f$ \phi-id \f$ is passed to apply and \f$ \phi_l-id,\phi_r-id \f$ are passed to the constructor as "DisplacementL" and "DisplacementR".
 *  The computation is done on \f$ \Omega \f$ (after a change of variables), which is passed to the constructor as "Grid".
 */
template <typename ConfiguratorType, typename MaterialLawType>
class HalfTransportEnergy :
  public aol::FENonlinIntegrationVectorInterface< ConfiguratorType, HalfTransportEnergy<ConfiguratorType,MaterialLawType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_l and \phi_r
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementL, _displacementR;

public:
  HalfTransportEnergy ( const typename ConfiguratorType::InitType &Grid,
                        const MaterialLawType &HyperelasticEnergyDensity,
                        const aol::MultiVector<RealType> &DisplacementL,
                        const aol::MultiVector<RealType> &DisplacementR ) :
    aol::FENonlinIntegrationVectorInterface < ConfiguratorType, HalfTransportEnergy<ConfiguratorType,MaterialLawType> > ( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementL( Grid, DisplacementL ),
    _displacementR( Grid, DisplacementR ) {}

  RealType evaluateIntegrand ( const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/ ) const {
    // compute the deformation gradients $\nabla\phi$, $\nabla\phi_l$, $\nabla\phi_r$
    typename ConfiguratorType::MatType dPhi, dPhiL, dPhiR, auxMat;
    _displacementL.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiL );
    _displacementR.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiR );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      DiscFuncs[i].evaluateGradientAtQuadPoint( El, QuadPoint, dPhi[i] );
      dPhi[i][i]  += 1.;
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute the deformation energy
    dPhiR *= dPhi.inverse();
    auxMat.makeProduct( dPhi, dPhiL.inverse() );
    return dPhiL.det() * qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( auxMat, _hyperelasticEnergyDensity, El, QuadPoint )
         + dPhi.det() * qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhiR, _hyperelasticEnergyDensity, El, QuadPoint );
  }

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::FENonlinIntegrationVectorInterface< ConfiguratorType, HalfTransportEnergy<ConfiguratorType,MaterialLawType> >::applyAdd( Arg, Dest );
    // fix translation
    int halfWidth = ( this->_initializer ).getWidth() / 8 * 7;
    typename ConfiguratorType::ElementType midNode( halfWidth, halfWidth, halfWidth );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      Dest.add( aol::Sqr( Arg[j][(this->_config)->localToGlobal( midNode, 0 )] ) );
    // fix rotation
    if ( ConfiguratorType::Dim == qc::QC_2D ) {
      midNode[0] += 1;
      Dest.add( aol::Sqr( Arg[1][(this->_config)->localToGlobal( midNode, 0 )] ) );
    } else
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        midNode[j] += 1;
        Dest.add( aol::Sqr( Arg[(j+1)%3][(this->_config)->localToGlobal( midNode, 0 )] ) );
        midNode[j] -= 1;
      }
  }
};

//! Computes the Gateaux derivative of HalfTransportEnergy wrt \f$ \phi \f$.
template <typename ConfiguratorType, typename MaterialLawType>
class HalfTransportDeriv :
  public aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, HalfTransportDeriv<ConfiguratorType,MaterialLawType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_l and \phi_r
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementL, _displacementR;

public:
  HalfTransportDeriv( const typename ConfiguratorType::InitType &Grid,
                      const MaterialLawType &HyperelasticEnergyDensity,
                      const aol::MultiVector<RealType> &DisplacementL,
                      const aol::MultiVector<RealType> &DisplacementR ) :
    aol::FENonlinVectorDiffOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim,HalfTransportDeriv<ConfiguratorType,MaterialLawType> > ( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementL( Grid, DisplacementL ),
    _displacementR( Grid, DisplacementR ) {}

  void getNonlinearity ( const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                         aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, RealType> &NL ) const {
    // compute the deformation gradients $D\phi$, $D\phi_l$, $D\phi_r$
    typename ConfiguratorType::MatType dPhi,dPhiL,dPhiR;
    _displacementL.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiL );
    _displacementR.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiR );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      DiscFuncs[i].evaluateGradientAtQuadPoint( El, QuadPoint, dPhi[i] );
      dPhi[i][i]  += 1.;
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute \partial_{,A}W(D\phi D\phi_l^{-1}), \partial_{,A}W(D\phi_r D\phi^{-1})
    typename ConfiguratorType::MatType cofDPhi, cofDPhiL, dPhiInvDPhiL, dPhiRInvDPhi, dWdPhi, dWdPhiR;
    cofDPhi.makeCofactorMatrix( dPhi );
    cofDPhiL.makeCofactorMatrix( dPhiL );
    dPhiInvDPhiL.makeProduct( dPhi, dPhiL.inverse() );
    dPhiRInvDPhi.makeProduct( dPhiR, dPhi.inverse() );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiInvDPhiL, _hyperelasticEnergyDensity, El, QuadPoint, dWdPhi );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiRInvDPhi, _hyperelasticEnergyDensity, El, QuadPoint, dWdPhiR );

    // compute NL = \partial_{,A}W(D\phi D\phi_l^{-1})cof(D\phi_l)+W(D\phi_r D\phi^{-1})cof(D\phi)-...
    // (D\phi_r D\phi^{-1})^T\partial_{,A}W(D\phi_r D\phi^{-1})cof(D\phi)
    NL.makeProductAtransposedB( dPhiRInvDPhi, dWdPhiR );
    NL *= cofDPhi;
    NL *= -1.;
    NL.addMultiple( cofDPhi, qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhiRInvDPhi, _hyperelasticEnergyDensity, El, QuadPoint ) );
    dWdPhi *= cofDPhiL;
    NL += dWdPhi;
  }

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::FENonlinVectorDiffOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim,HalfTransportDeriv<ConfiguratorType,MaterialLawType> >::applyAdd( Arg, Dest );
    // fix translation
    int halfWidth = ( this->_initializer ).getWidth() / 8 * 7;
    typename ConfiguratorType::ElementType midNode( halfWidth, halfWidth, halfWidth );
    int index = (this->_config)->localToGlobal( midNode, 0 );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      Dest[j][index] += 2. * Arg[j][index];
    // fix rotation
    if ( ConfiguratorType::Dim == qc::QC_2D ) {
      midNode[0] += 1;
      index = (this->_config)->localToGlobal( midNode, 0 );
      Dest[1][index] += 2. * Arg[1][index];
    } else
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        midNode[j] += 1;
        index = (this->_config)->localToGlobal( midNode, 0 );
        Dest[(j+1)%3][index] += 2. * Arg[(j+1)%3][index];
        midNode[j] -= 1;
      }
  }
};

//! Computes (k,l)-block of the Hessian of HalfTransportEnergy wrt \f$ \phi \f$.
template <typename ConfiguratorType, typename MaterialLawType>
class HalfTransportSubHessian :
  public aol::FELinAsymMatrixWeightedStiffInterface<ConfiguratorType,HalfTransportSubHessian<ConfiguratorType,MaterialLawType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_l and \phi_r and \phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementL, _displacementR, _displacement;
  // the row and column index of the sub-Hessian
  const int _k, _l;

public:
  HalfTransportSubHessian ( const typename ConfiguratorType::InitType &Grid,
                       const MaterialLawType &HyperelasticEnergyDensity,
                       const aol::MultiVector<RealType> &DisplacementL,
                       const aol::MultiVector<RealType> &DisplacementR,
                       const aol::MultiVector<RealType> &Displacement,
                       const int K,
                       const int L ) :
    aol::FELinAsymMatrixWeightedStiffInterface<ConfiguratorType,HalfTransportSubHessian<ConfiguratorType,MaterialLawType> > ( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementL( Grid, DisplacementL ),
    _displacementR( Grid, DisplacementR ),
    _displacement( Grid, Displacement ),
    _k( K ),
    _l( L ) {}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*LocalCoord*/,
                               typename ConfiguratorType::MatType &Matrix ) const {
    // compute the deformation gradients $D\phi$, $D\phi_l$, $D\phi_r$
    typename ConfiguratorType::MatType dPhi,dPhiL,dPhiR;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    _displacementL.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiL );
    _displacementR.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiR );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      dPhi[i][i]  += 1.;
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute D\phi D\phi_l^{-1}, D\phi_r D\phi^{-1} and W, \partial_{,A}W
    typename ConfiguratorType::MatType dPhiInvDPhiL, dPhiRInvDPhi, dPhiInv( dPhi.inverse() ), dPhiLInv( dPhiL.inverse() ), dW;
    dPhiInvDPhiL.makeProduct( dPhi, dPhiLInv );
    dPhiRInvDPhi.makeProduct( dPhiR, dPhiInv );
    RealType W = qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhiRInvDPhi, _hyperelasticEnergyDensity, El, QuadPoint );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiRInvDPhi, _hyperelasticEnergyDensity, El, QuadPoint, dW );

    // compute the Hessian matrix; denote variations of \phi by \psi and \chi; they only act on the component "_k" and "_l"
    // \psi shall belong to vectors multiplied from the left, \chi from the right
    typename ConfiguratorType::MatType cofDPhi, cofDPhiL, aux, aux2;
    typename ConfiguratorType::VecType auxVec1, auxVec2;
    cofDPhi.makeCofactorMatrix( dPhi );
    cofDPhiL.makeCofactorMatrix( dPhiL );

    qc::HyperelasticSubHessian<ConfiguratorType,MaterialLawType>::elasticityTensor( dPhiInvDPhiL, _hyperelasticEnergyDensity, El, QuadPoint, _k, _l, aux );
    aux *= cofDPhiL;
    Matrix.makeProduct( dPhiLInv, aux );

    if ( _k != _l ) {
      aux.makeProductABtransposed( dPhiInv, dW );
      aux *= dPhiRInvDPhi;
      aux.getColumn( _l, auxVec1 );
      aux.getColumn( _k, auxVec2 );
      aux.makeTensorProduct( cofDPhi[_k], auxVec1 );
      Matrix -= aux;

      aux.transpose();
      Matrix += aux;

      aux.makeTensorProduct( cofDPhi[_l], auxVec2 );
      Matrix += aux;

      aux.transpose();
      Matrix -= aux;

      // M s.t. D\psi_k M D\chi_l = W(D\phi_r D\phi^{-1}) \partial_{,A}^2 det(D\phi)(D\psi,D\chi)
      aux.makeTensorProduct( cofDPhi[_k], cofDPhi[_l] );
      Matrix.addMultiple( aux, W / dPhi.det() );
      aux.transpose();
      Matrix.addMultiple( aux, -W / dPhi.det() );
    }

    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        qc::HyperelasticSubHessian<ConfiguratorType,MaterialLawType>::elasticityTensor( dPhiRInvDPhi, _hyperelasticEnergyDensity, El, QuadPoint, i, j, aux );
        aux *= cofDPhi;
        aux2.makeProduct( dPhiInv, aux );
        Matrix.addMultiple( aux2, dPhiRInvDPhi[i][_k] * dPhiRInvDPhi[j][_l] );
      }
  }
};

//! Computes the Hessian of HalfTransportEnergy wrt \f$ \phi \f$.
template <typename ConfiguratorType, typename MaterialLawType>
class HalfTransportHessian :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::BlockOpBase<typename ConfiguratorType::RealType,typename ConfiguratorType::MatrixType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  const ConfiguratorType _config;
  const typename ConfiguratorType::InitType &_grid;
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_l and \phi_r
  const aol::MultiVector<RealType> &_displacementL, &_displacementR;

public:
  HalfTransportHessian( const typename ConfiguratorType::InitType &Grid,
                        const MaterialLawType &HyperelasticEnergyDensity,
                        const aol::MultiVector<RealType> &DisplacementL,
                        const aol::MultiVector<RealType> &DisplacementR ) :
    _config( Grid ),
    _grid( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementL( DisplacementL ),
    _displacementR( DisplacementR ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::BlockOpBase<RealType,typename ConfiguratorType::MatrixType> &Dest ) const {
    // assemble the different components of the Derivative
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        HalfTransportSubHessian<ConfiguratorType,MaterialLawType>( _grid, _hyperelasticEnergyDensity, _displacementL, _displacementR, Arg, i, j ).assembleAddMatrix( Dest.getReference( i, j ) );
    // fix translation
    int halfWidth = _grid.getWidth()  / 8 * 7;
    typename ConfiguratorType::ElementType midNode( halfWidth, halfWidth, halfWidth );
    int index = _config.localToGlobal( midNode, 0 );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      Dest.getReference( j, j ).add( index, index, 2. );
    // fix rotation
    if ( ConfiguratorType::Dim == qc::QC_2D ) {
      midNode[0] += 1;
      index = _config.localToGlobal( midNode, 0 );
      Dest.getReference( 1, 1 ).add( index, index, 2. );
    } else
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        midNode[j] += 1;
        index = _config.localToGlobal( midNode, 0 );
        Dest.getReference( (j+1)%3, (j+1)%3 ).add( index, index, 2. );
        midNode[j] -= 1;
      }
  }
};

/** Computes the energy \f$ \int_{\Omega} W(D\phi_d) dx + \int_{\phi_d(\Omega)} W(D(\phi\circ\phi_s\circ\phi_d^{-1})) dx \f$,
 *  where \f$ \phi_d-id \f$ is passed to apply and \f$ \phi_s-id,\phi-id \f$ are passed to the constructor as "DisplacementSide" and "Displacement".
 *  The computation is done on \f$ \Omega \f$ (after a change of variables), which is passed to the constructor as "Grid".
 */
template <typename ConfiguratorType, typename MaterialLawType>
class DiagEnergy :
  public qc::FENonlinDeformIntegrationVectorInterface < ConfiguratorType, DiagEnergy<ConfiguratorType,MaterialLawType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_s and \phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementSide, _displacement;

public:
  DiagEnergy ( const typename ConfiguratorType::InitType &Grid,
               const MaterialLawType &HyperelasticEnergyDensity,
               const aol::MultiVector<RealType> &DisplacementSide,
               const aol::MultiVector<RealType> &Displacement ) :
    qc::FENonlinDeformIntegrationVectorInterface < ConfiguratorType, DiagEnergy<ConfiguratorType,MaterialLawType> > ( Grid, DisplacementSide ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementSide( Grid, DisplacementSide ),
    _displacement( Grid, Displacement ) {}

  RealType evaluateIntegrand ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord ) const {
    // compute the deformation gradients $\nabla\phi$, $\nabla\phi_s$, $\nabla\phi_d$
    typename ConfiguratorType::MatType dPhi,dPhiS,dPhiD;
    _displacement.evaluateGradient( TransformedEl, TransformedLocalCoord, dPhi );
    _displacementSide.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiS );
    DiscFuncs.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiD );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      dPhi[i][i]  += 1.;
      dPhiS[i][i] += 1.;
      dPhiD[i][i] += 1.;
    }

    // compute the deformation energy
    dPhi *= dPhiS;
    dPhi *= dPhiD.inverse();
    return qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhiD, _hyperelasticEnergyDensity, El, QuadPoint )
         + dPhiD.det() * qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhi, _hyperelasticEnergyDensity, El, QuadPoint );
  }
};

/** Computes the derivative of \f$ \int_{\Omega} W(D\phi_d) dx + \int_{\phi_d(\Omega)} W(D(\phi\circ\phi_s\circ\phi_d^{-1})) dx \f$ wrt \f$ \phi_d \f$,
 *  where \f$ \phi-id \f$ is passed to apply and \f$ \phi_s-id,\phi_d-id \f$ are passed to the constructor as "DisplacementSide" and "DisplacementDiag".
 *  The computation is done on \f$ \Omega \f$ (after a change of variables), which is passed to the constructor as "Grid".
 */
template <typename ConfiguratorType, typename MaterialLawType>
class EndPointConstr :
  public aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, EndPointConstr<ConfiguratorType,MaterialLawType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_s and \phi_d
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementSide, _displacementDiag;
  // mask showing which elements lie in the image of _displacementSide
  mutable qc::BitArray<ConfiguratorType::Dim> _mask;
  
public:
  EndPointConstr ( const typename ConfiguratorType::InitType &Grid,
                   const MaterialLawType &HyperelasticEnergyDensity,
                   const aol::MultiVector<RealType> &DisplacementSide,
                   const aol::MultiVector<RealType> &DisplacementDiag ) :
    aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, EndPointConstr<ConfiguratorType,MaterialLawType> > ( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementSide( Grid, DisplacementSide ),
    _displacementDiag( Grid, DisplacementDiag ),
    _mask( Grid ) {}

  const qc::BitArray<ConfiguratorType::Dim>& getMaskReference() const {
    return _mask;
  }

  void getNonlinearity ( const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &LocalCoord,
                         aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, RealType> &NL ) const {
    // compute $\phi_s(x)$
    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType  transformedCoord;
    qc::transformAndClipCoord<ConfiguratorType>( *this->_config, _displacementSide, El, QuadPoint, LocalCoord, transformedEl, transformedCoord );
    _mask.set( transformedEl, true );

    // compute the deformation gradients $\nabla\phi$, $\nabla\phi_s$, $\nabla\phi_d$
    typename ConfiguratorType::MatType dPhi,dPhiS,dPhiD;
    _displacementSide.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiS );
    _displacementDiag.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiD );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      DiscFuncs[i].evaluateGradient( transformedEl, transformedCoord, dPhi[i] );
      dPhi[i][i]  += 1.;
      dPhiS[i][i] += 1.;
      dPhiD[i][i] += 1.;
    }

    // compute \partial_{,A}W(\nabla\phi_d), \partial_{,A}W(\nabla\phi\nabla\phi_s\nabla\phi_d^{-1})
    typename ConfiguratorType::MatType cofDPhiD, dPhiSD( dPhi ), dWdPhiSD, dW;
    dPhiSD *= dPhiS;
    dPhiSD *= dPhiD.inverse();
    cofDPhiD.makeCofactorMatrix( dPhiD );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiD, _hyperelasticEnergyDensity, El, QuadPoint, dW );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiSD, _hyperelasticEnergyDensity, El, QuadPoint, dWdPhiSD );

    // compute NL = \partial_{,A}W(\nabla\phi_d)+W(\nabla\phi\nabla\phi_s\nabla\phi_d^{-1})cof(\nabla\phi_d)-...
    // cof(\nabla\phi_d)\nabla\phi_s^T\nabla\phi^T\partial_{,A}W(\nabla\phi\nabla\phi_s\nabla\phi_d^{-1})\nabla\phi_d^{-T}
    NL.makeProductAtransposedB( dPhiSD, dWdPhiSD );
    NL *= cofDPhiD;
    NL *= -1.;
    NL.addMultiple( cofDPhiD, qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhiSD, _hyperelasticEnergyDensity, El, QuadPoint ) );
    NL += dW;
  }
};

//! Computes the Gateaux derivative of DiagEnergy wrt \f$ \phi_d \f$.
template <typename ConfiguratorType, typename MaterialLawType>
class DiagDeriv :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_s and \phi
  const aol::MultiVector<RealType> _displacementSide, _displacement;

public:
  DiagDeriv ( const typename ConfiguratorType::InitType &Grid,
              const MaterialLawType &HyperelasticEnergyDensity,
              const aol::MultiVector<RealType> &DisplacementSide,
              const aol::MultiVector<RealType> &Displacement ) :
    _grid( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementSide( DisplacementSide ),
    _displacement( Displacement ) {}

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    EndPointConstr<ConfiguratorType,MaterialLawType>( _grid, _hyperelasticEnergyDensity, _displacementSide, Arg ).applyAdd( _displacement, Dest );
  }
};

//! Computes (k,l)-block of the derivative of EndPointConstr wrt \f$ \phi \f$.
template <typename ConfiguratorType, typename MaterialLawType>
class EndPointSubDeriv :
  public qc::FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType,EndPointSubDeriv<ConfiguratorType,MaterialLawType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_l and \phi_r and \phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementSide, _displacementDiag, _displacement;
  // the row and column index of the sub-Derivative
  const int _k, _l;
  // mask showing which elements lie in the image of _displacementSide
  mutable qc::BitArray<ConfiguratorType::Dim> _mask;
  
public:
  EndPointSubDeriv ( const typename ConfiguratorType::InitType &Grid,
                     const MaterialLawType &HyperelasticEnergyDensity,
                     const aol::MultiVector<RealType> &DisplacementSide,
                     const aol::MultiVector<RealType> &DisplacementDiag,
                     const aol::MultiVector<RealType> &Displacement,
                     const int K,
                     const int L ) :
    qc::FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType,EndPointSubDeriv<ConfiguratorType,MaterialLawType> > ( Grid, DisplacementSide, aol::ONTHEFLY, qc::RIGHT ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementSide( Grid, DisplacementSide ),
    _displacementDiag( Grid, DisplacementDiag ),
    _displacement( Grid, Displacement ),
    _k( K ),
    _l( L ),
    _mask( Grid ) {}

  const qc::BitArray<ConfiguratorType::Dim>& getMaskReference() const {
    return _mask;
  }

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*LocalCoord*/,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                               typename ConfiguratorType::MatType &Matrix ) const {
    // compute the deformation gradients $\nabla\phi$, $\nabla\phi_s$, $\nabla\phi_d$
    typename ConfiguratorType::MatType dPhi,dPhiS,dPhiD;
    _displacement.evaluateGradient( TransformedEl, TransformedLocalCoord, dPhi );
    _displacementSide.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiS );
    _displacementDiag.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiD );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      dPhi[i][i]  += 1.;
      dPhiS[i][i] += 1.;
      dPhiD[i][i] += 1.;
    }

    // compute \partial_{,A}W(\nabla\phi\nabla\phi_s\nabla\phi_d^{-1})
    typename ConfiguratorType::MatType dPhiSD( dPhi ), invDPhiD( dPhiD.inverse() ), cofDPhiD, dW;
    dPhiSD *= dPhiS;
    dPhiSD *= invDPhiD;
    cofDPhiD.makeCofactorMatrix( dPhiD );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiSD, _hyperelasticEnergyDensity, El, QuadPoint, dW );

    // compute the derivative matrix; denote variations of \phi_d and \phi\circ\phi_s by \psi and \chi; they only act on the component "_k" and "_l"
    // \psi shall belong to vectors multiplied from the left, \chi from the right
    typename ConfiguratorType::MatType aux1, aux2;
    // M s.t. D\psi_k^T M D\chi_l = [\partial_{,A}W(D\phi D\phi_s D\phi_d^{-1}) D\phi_d^{-T} : D\chi][cof D\phi_d : D\psi]
    Matrix.makeTensorProduct( cofDPhiD[_k], invDPhiD * dW[_l] );
    // M s.t. D\psi_k^T M D\chi_l = -cof D\phi_d D\chi^T \partial_{,A}W(D\phi D\phi_s D\phi_d^{-1}) D\phi_d^{-T} : D\psi
    Matrix -= Matrix.transposed();
    // M s.t. D\psi_k^T M D\chi_l = -\partial_{,AA}W(D\phi D\phi_s D\phi_d^{-1})(D\phi D\phi_s D\phi_d^{-1} D\psi cof D\phi_d^T)(D\chi D\phi_d^{-1})
    for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
      qc::HyperelasticSubHessian<ConfiguratorType,MaterialLawType>::elasticityTensor( dPhiSD, _hyperelasticEnergyDensity, El, QuadPoint, j, _l, aux1 );
      aux1 *= cofDPhiD;
      aux2.makeProduct( invDPhiD, aux1 );
      Matrix.addMultiple( aux2, -dPhiSD[j][_k] );
    }
  }
};

//! Computes the derivative of EndPointConstr wrt \f$ \phi \f$.
template <typename ConfiguratorType, typename MaterialLawType, typename MatrixType>
class EndPointDeriv :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, MatrixType > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  const ConfiguratorType _config;
  const typename ConfiguratorType::InitType &_grid;
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_s and \phi_d
  const aol::MultiVector<RealType> &_displacementSide, &_displacementDiag;

public:
  EndPointDeriv( const typename ConfiguratorType::InitType &Grid,
                 const MaterialLawType &HyperelasticEnergyDensity,
                 const aol::MultiVector<RealType> &DisplacementSide,
                 const aol::MultiVector<RealType> &DisplacementDiag ) :
    _config( Grid ),
    _grid( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementSide( DisplacementSide ),
    _displacementDiag( DisplacementDiag ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, MatrixType &Dest ) const {
    // assemble the different components of the Derivative
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        EndPointSubDeriv<ConfiguratorType,MaterialLawType>( _grid, _hyperelasticEnergyDensity, _displacementSide, _displacementDiag, Arg, i, j ).assembleAddMatrix( Dest.getReference( i, j ) );
  }
};

/** Computes the derivative of \f$ \int_{\Omega} W(D\phi_d) dx + \int_{\phi_h(\Omega)} W(D(\phi\circ\phi_h^{-1})) dx \f$ wrt \f$ \phi_h \f$,
 *  where \f$ \phi-id \f$ is passed to apply and \f$ \phi_h-id \f$ is passed to the constructor as "DisplacementHalf".
 *  The computation is done on \f$ \Omega \f$ (after a change of variables), which is passed to the constructor as "Grid".
 */
template <typename ConfiguratorType, typename MaterialLawType>
class DoubleTransportConstr :
  public aol::FENonlinVectorDiffOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim,DoubleTransportConstr<ConfiguratorType,MaterialLawType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformation \phi_h
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementHalf;

public:
  DoubleTransportConstr ( const typename ConfiguratorType::InitType &Grid,
                          const MaterialLawType &HyperelasticEnergyDensity,
                          const aol::MultiVector<RealType> &DisplacementHalf ) :
    aol::FENonlinVectorDiffOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim,DoubleTransportConstr<ConfiguratorType,MaterialLawType> > ( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementHalf( Grid, DisplacementHalf ) {}

  void getNonlinearity ( const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                         aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, RealType> &NL ) const {
    // compute the deformation gradients $\nabla\phi$, $\nabla\phi_h$
    typename ConfiguratorType::MatType dPhi,dPhiH;
    _displacementHalf.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiH );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      DiscFuncs[i].evaluateGradientAtQuadPoint( El, QuadPoint, dPhi[i] );
      dPhi[i][i]  += 1.;
      dPhiH[i][i] += 1.;
    }

    // compute \partial_{,A}W(\nabla\phi_h), \partial_{,A}W(\nabla\phi\nabla\phi_h^{-1})
    typename ConfiguratorType::MatType cofDPhiH, dPhiInvDPhiH( dPhi ), dWdPhiInvDPhiH, dW;
    dPhiInvDPhiH *= dPhiH.inverse();
    cofDPhiH.makeCofactorMatrix( dPhiH );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiH, _hyperelasticEnergyDensity, El, QuadPoint, dW );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiInvDPhiH, _hyperelasticEnergyDensity, El, QuadPoint, dWdPhiInvDPhiH );

    // compute NL = \partial_{,A}W(\nabla\phi_h)+W(\nabla\phi\nabla\phi_h^{-1})cof(\nabla\phi_h)-...
    // cof(\nabla\phi_h)\nabla\phi^T\partial_{,A}W(\nabla\phi\nabla\phi_h^{-1})\nabla\phi_h^{-T}
    NL.makeProductAtransposedB( dPhiInvDPhiH, dWdPhiInvDPhiH );
    NL *= cofDPhiH;
    NL *= -1.;
    NL.addMultiple( cofDPhiH, qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>::energyDensity( dPhiInvDPhiH, _hyperelasticEnergyDensity, El, QuadPoint ) );
    NL += dW;
  }
  
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::FENonlinVectorDiffOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim,DoubleTransportConstr<ConfiguratorType,MaterialLawType> >::applyAdd( Arg, Dest );
    // fix translation
    int halfWidth2 = ( this->_initializer ).getWidth() / 8 * 5;
    int halfWidth = ( this->_initializer ).getWidth() / 8 * 7;
    typename ConfiguratorType::ElementType midNode( halfWidth2, halfWidth, halfWidth );
    int index = (this->_config)->localToGlobal( midNode, 0 );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      Dest[j][index] += 2. * Arg[j][index];
    // fix rotation
    if ( ConfiguratorType::Dim == qc::QC_2D ) {
      midNode[0] += 1;
      index = (this->_config)->localToGlobal( midNode, 0 );
      Dest[1][index] += 2. * Arg[1][index];
    } else
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        midNode[j] += 1;
        index = (this->_config)->localToGlobal( midNode, 0 );
        Dest[(j+1)%3][index] += 2. * Arg[(j+1)%3][index];
        midNode[j] -= 1;
      }
  }
};

//! Computes (k,l)-block of the derivative of DoubleTransportConstr wrt \f$ \phi \f$.
template <typename ConfiguratorType, typename MaterialLawType>
class DoubleTransportSubDeriv :
  public aol::FELinAsymMatrixWeightedStiffInterface<ConfiguratorType,DoubleTransportSubDeriv<ConfiguratorType,MaterialLawType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformations \phi_h and \phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacementHalf, _displacement;
  // the row and column index of the sub-Derivative
  const int _k, _l;

public:
  DoubleTransportSubDeriv ( const typename ConfiguratorType::InitType &Grid,
                            const MaterialLawType &HyperelasticEnergyDensity,
                            const aol::MultiVector<RealType> &DisplacementHalf,
                            const aol::MultiVector<RealType> &Displacement,
                            const int K,
                            const int L ) :
    aol::FELinAsymMatrixWeightedStiffInterface<ConfiguratorType,DoubleTransportSubDeriv<ConfiguratorType,MaterialLawType> > ( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementHalf( Grid, DisplacementHalf ),
    _displacement( Grid, Displacement ),
    _k( K ),
    _l( L ) {}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*LocalCoord*/,
                               typename ConfiguratorType::MatType &Matrix ) const {
    // compute the deformation gradients $\nabla\phi$, $\nabla\phi_h$
    typename ConfiguratorType::MatType dPhi,dPhiH;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    _displacementHalf.evaluateGradientAtQuadPoint( El, QuadPoint, dPhiH );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      dPhi[i][i]  += 1.;
      dPhiH[i][i] += 1.;
    }

    // compute \partial_{,A}W(\nabla\phi\nabla\phi_h^{-1})
    typename ConfiguratorType::MatType cofDPhiH, dPhiInvDPhiH( dPhi ), invDPhiH( dPhiH.inverse() ), dW;
    dPhiInvDPhiH *= invDPhiH;
    cofDPhiH.makeCofactorMatrix( dPhiH );
    qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>::firstPiolaKirchhoffStress( dPhiInvDPhiH, _hyperelasticEnergyDensity, El, QuadPoint, dW );

    // compute the derivative matrix; denote variations of \phi_h and \phi by \psi and \chi; they only act on the component "_k" and "_l"
    // \psi shall belong to vectors multiplied from the left, \chi from the right
    typename ConfiguratorType::MatType aux;
    // M s.t. D\psi_k^T M D\chi_l = [\partial_{,A}W(D\phi D\phi_h^{-1}) D\phi_h^{-T} : D\chi][cof D\phi_h : D\psi]
    aux[_k] += dW[_l];
    // M s.t. D\psi_k^T M D\chi_l = -cof D\phi_h D\chi^T \partial_{,A}W(D\phi D\phi_h^{-1}) D\phi_h^{-T} : D\psi
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      aux[i][_k] -= dW[_l][i];
    // M s.t. D\psi_k^T M D\chi_l = -\partial_{,AA}W(D\phi D\phi_h^{-1})(D\phi D\phi_h^{-1} D\psi cof D\phi_h^T)(D\chi D\phi_h^{-1})
    for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
      qc::HyperelasticSubHessian<ConfiguratorType,MaterialLawType>::elasticityTensor( dPhiInvDPhiH, _hyperelasticEnergyDensity, El, QuadPoint, j, _l, Matrix );
      aux.addMultiple( Matrix, -dPhiInvDPhiH[j][_k] );
    }

    aux *= cofDPhiH;
    Matrix.makeProduct( invDPhiH, aux );
  }
};

//! Computes the derivative of DoubleTransportConstr wrt \f$ \phi \f$.
template <typename ConfiguratorType, typename MaterialLawType, typename MatrixType>
class DoubleTransportDeriv :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, MatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  const ConfiguratorType _config;
  const typename ConfiguratorType::InitType &_grid;
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const MaterialLawType &_hyperelasticEnergyDensity;
  // the deformation \phi_h
  const aol::MultiVector<RealType> &_displacementHalf;

public:
  DoubleTransportDeriv( const typename ConfiguratorType::InitType &Grid,
                        const MaterialLawType &HyperelasticEnergyDensity,
                        const aol::MultiVector<RealType> &DisplacementHalf ) :
    _config( Grid ),
    _grid( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacementHalf( DisplacementHalf ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, MatrixType &Dest ) const {
    // assemble the different components of the Derivative
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        DoubleTransportSubDeriv<ConfiguratorType,MaterialLawType>( _grid, _hyperelasticEnergyDensity, _displacementHalf, Arg, i, j ).assembleAddMatrix( Dest.getReference( i, j ) );
    // fix translation
    int halfWidth2 = _grid.getWidth() / 8 * 5;
    int halfWidth = _grid.getWidth() / 8 * 7;
    typename ConfiguratorType::ElementType midNode( halfWidth2, halfWidth, halfWidth );
    int index = _config.localToGlobal( midNode, 0 );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      Dest.getReference( j, j ).add( index, index, 2. );
    // fix rotation
    if ( ConfiguratorType::Dim == qc::QC_2D ) {
      midNode[0] += 1;
      index = _config.localToGlobal( midNode, 0 );
      Dest.getReference( 1, 1 ).add( index, index, 2. );
    } else
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        midNode[j] += 1;
        index = _config.localToGlobal( midNode, 0 );
        Dest.getReference( (j+1)%3, (j+1)%3 ).add( index, index, 2. );
        midNode[j] -= 1;
      }
  }
};

template <typename ConfiguratorType,typename MaterialLawType>
class ParallelTransportOp :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::VectorContainer<aol::MultiVector<typename ConfiguratorType::RealType> > > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  typedef typename ConfiguratorType::MatrixType SubMatrixType;

  // grid on which displacements live
  const GridType _grid;
  // displacements along the discrete geodesic
  const int _numTimePoints;
  const aol::RandomAccessContainer<aol::MultiVector<RealType> > &_geodesicDisplacements;
  // elastic deformation energy density
  const aol::RandomAccessContainer<MaterialLawType> &_elasticMaterialLaws;
  // number of Newton iterations
  const int _maxIterations;

public:
  // note: first component of GeodesicDisplacements is a dummy
  ParallelTransportOp( const GridType &Grid,
                       const aol::RandomAccessContainer<aol::MultiVector<RealType> > &GeodesicDisplacements,
                       const aol::RandomAccessContainer<MaterialLawType> &ElasticMaterialLaws,
                       const int MaxIterations = 1e3 ) :
    _grid( Grid ),
    _numTimePoints( GeodesicDisplacements.size() ),
    _geodesicDisplacements( GeodesicDisplacements ),
    _elasticMaterialLaws( ElasticMaterialLaws ),
    _maxIterations( MaxIterations ) {}

  static void computeParallelogramMidpoint( const GridType &Grid,
                                            const MaterialLawType &MaterialLaw,
                                            const aol::MultiVector<RealType> &DispL,
                                            const aol::MultiVector<RealType> &DispR,
                                            aol::MultiVector<RealType> &Disp,
                                            const int MaxIterations = 1e3,
                                            const RealType MaxAccuracy = 1e-7 ) {
    // prepare energy and derivatives
    HalfTransportEnergy<ConfiguratorType,MaterialLawType>    e( Grid, MaterialLaw, DispL, DispR );
    HalfTransportDeriv<ConfiguratorType,MaterialLawType>    de( Grid, MaterialLaw, DispL, DispR );
    HalfTransportHessian<ConfiguratorType,MaterialLawType> d2e( Grid, MaterialLaw, DispL, DispR );
    // initialize minimization
    aol::TrustRegionMethod<RealType,aol::MultiVector<RealType>,aol::BlockOpBase<RealType,SubMatrixType>,aol::DiagonalPreconditioner<aol::MultiVector<RealType> >,SubMatrixType> minimizer( e, de, d2e, new aol::BlockMatrix<SubMatrixType>( Grid ), MaxIterations, MaxAccuracy );
    minimizer.setIterativeMode( false );
    // perform minimization
    aol::MultiVector<RealType> initialGuess( Disp, aol::DEEP_COPY );
    minimizer.apply( initialGuess, Disp );
  }

  static void computeParallelogramEndpoint( const GridType &Grid,
                                            const MaterialLawType &MaterialLaw,
                                            const aol::MultiVector<RealType> &DispSide,
                                            const aol::MultiVector<RealType> &DispDiag,
                                            aol::MultiVector<RealType> &Disp,
                                            const int MaxIterations = 1e3 ) {
    // prepare function of which to find a zero and its derivative
    EndPointConstr<ConfiguratorType,MaterialLawType> f( Grid, MaterialLaw, DispSide, DispDiag );
    EndPointDeriv<ConfiguratorType,MaterialLawType,aol::SparseBlockMatrix<aol::SparseMatrix<RealType> > > df( Grid, MaterialLaw, DispSide, DispDiag );
    // initialize Newton's method
    NewtonIteration<ConfiguratorType,aol::MultiVector<RealType>,aol::SparseBlockMatrix<aol::SparseMatrix<RealType> >,aol::SparseMatrix<RealType> > zeroFinder( Grid, f, df, MaxIterations, 1.e-7 );
    // perform zero finding
    aol::MultiVector<RealType> initialGuess( Disp, aol::DEEP_COPY );
    zeroFinder.apply( initialGuess, Disp );
  }

  static void shootGeodesic ( const GridType &Grid,
                              const MaterialLawType &MaterialLaw,
                              const aol::MultiVector<RealType> &DispHalf,
                              aol::MultiVector<RealType> &Disp,
                              const int MaxIterations = 1e3,
                              const RealType MaxAccuracy = 1e-7 ) {
    // prepare function of which to find a zero and its derivative
    DoubleTransportConstr<ConfiguratorType,MaterialLawType> f( Grid, MaterialLaw, DispHalf );
    DoubleTransportDeriv<ConfiguratorType,MaterialLawType,aol::SparseBlockMatrix<aol::SparseMatrix<RealType> > > df( Grid, MaterialLaw, DispHalf );
    // initialize Newton's method
    NewtonIteration<ConfiguratorType,aol::MultiVector<RealType>,aol::SparseBlockMatrix<aol::SparseMatrix<RealType> >,aol::SparseMatrix<RealType> > zeroFinder( Grid, f, df, MaxIterations, MaxAccuracy );
    zeroFinder.setTimestepController( aol::NewtonInfo<RealType>::ARMIJO );
    zeroFinder.setSigma( .1 );
    // perform zero finding
    aol::MultiVector<RealType> initialGuess( Disp, aol::DEEP_COPY );
    zeroFinder.apply( initialGuess, Disp );
  }

  void computeParallelTransport( const aol::MultiVector<RealType> &TangentVector,
                                 aol::RandomAccessContainer<aol::MultiVector<RealType> > &MidPointDisps,
                                 aol::RandomAccessContainer<aol::MultiVector<RealType> > &TransportedVector ) const {
    TransportedVector[0] = TangentVector;
    for( int k = 1; k < _numTimePoints; k++ ) {
      //! transport vector along kth geodesic segment
      // compute parallelogram midpoint
      computeParallelogramMidpoint( _grid, _elasticMaterialLaws[k-1], TransportedVector[k-1], _geodesicDisplacements[k], MidPointDisps[k], _maxIterations );
      // compute transported vector
      aol::MultiVector<RealType> continuedDisp( _grid );
      qc::ConcatAndSmoothlyExtendDeformations<ConfiguratorType>( TransportedVector[k], _geodesicDisplacements[k], _grid, continuedDisp );
      shootGeodesic( _grid, _elasticMaterialLaws[k-1], MidPointDisps[k], continuedDisp, _maxIterations );
      qc::InvConcatAndSmoothlyExtendDeformations<ConfiguratorType>( continuedDisp, _geodesicDisplacements[k], _grid, TransportedVector[k] );
    }
  }

  // Dest[i] will contain the displacement Arg transported to the ith node of the discrete geodesic
  virtual void applyAdd( const aol::MultiVector<RealType> &Arg, aol::VectorContainer<aol::MultiVector<RealType> > &Dest ) const {
    aol::VectorContainer<aol::MultiVector<RealType> > dest( Dest ), midPointDisps( Dest );
    computeParallelTransport( Arg, midPointDisps, dest );
    Dest += dest;
  }
};

//! Performs computation of a discrete geodesic as described in Rumpf and Wirth, "DISCRETE GEODESIC CALCULUS IN THE SPACE OF VISCOUS FLUIDIC OBJECTS", chapter 4.1
//! \author Wirth
//! In particular, see figure 4.1 for reference. The (\phi_k)_k are the genuine degrees of freedom, denoted as "correction displacements" here.
//! The reference displacements (\bar \psi_k)_k are assumed to be given and denoted by "morphing displacements" here.
//! Note that the main advantage of this approach is having only K deformations as DOFs instead of K deformations AND K+1 shapes!
//! The drawback is having concatenated arguments and the assumption of suitbale reference shapes and reference deformations to be given.
template <typename ConfiguratorType>
class GeodesicComputationViaFastMorphing {

private:
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef aol::SparseMatrix<RealType> MatrixType;
  typedef qc::HyperelasticEnergyDensity2DLinear<ConfiguratorType> MaterialLawType;
  typedef qc::AntiCompressionHyperelasticEnergyDensity<ConfiguratorType,MaterialLawType> MatLawType;

  // finest grid over the examined domain
  const qc::GridSize<ConfiguratorType::Dim> _gridSize;
  const typename ConfiguratorType::InitType _grid;
  // number of time points along the morphing path
  const int _numberOfTimePoints;
  // given initial and final object as well as objects along the morphing path
  qc::MultiDimMultilevelArray< RealType, ArrayType > _endObjects, _objects;
  // the directory for the results
  char _destDirectory[1024];
  // the hyperelastic deformation energy density
  const MaterialLawType _linearElasticEnergyDensity;
  const MatLawType _elasticEnergyDensity;
  const RealType _weakMaterialStiffness;
  // further parameters
  const bool _withMatching;
  aol::Vector<RealType> _matchingWeight, _smoothingPixelWidth;
  // levels for the multiscale method
  const int _finestLevel;
  aol::Vector<int> _level;
  // minimization method, descent steps and accuracy
  aol::Vector<int> _method, _mode, _maxDescentSteps;
  aol::Vector<RealType> _maxAccuracy;
  // 0: no shape reinitialisation, 1: last shape is transported backwards, 2: first shape is transported forwards
  aol::Vector<int> _shapeReinitMode;
  // flag indicating whether the initial morphing path shall be updated during the algorithm
  aol::BitVector _doPathReinitialization;
  // matching displacements along morphing path and correction displacements of each object along the path
  qc::MultiDimMultilevelArray<RealType, qc::MultiArray<RealType, ConfiguratorType::Dim>, qc::ProlongOp<RealType>,qc::RestrictOp<RealType,qc::THROW_AWAY_RESTRICT> > _morphingDisplacements, _correctionDisplacements;

  void generateDeformationMask( const aol::MultiVector<RealType> &Displacement, qc::BitArray<ConfiguratorType::Dim> &Mask ) {
    const qc::GridSize<ConfiguratorType::Dim> gridSize( Mask );
    const typename ConfiguratorType::InitType grid( gridSize );
    const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> displacement( grid, Displacement );
    Mask.setZero();

    typename ConfiguratorType::InitType::OldAllNodeIterator fnit;
    for ( fnit = grid._nBeginIt; fnit != grid._nEndIt; ++fnit ) {
      qc::Element transformedEl;
      typename ConfiguratorType::DomVecType localCoord, transformedLocalCoord;
      if ( qc::transformCoord<ConfiguratorType>( grid, displacement, *fnit, 0, localCoord, transformedEl, transformedLocalCoord ) )
        Mask.set( transformedEl, true );
    }
  }
  
public:
  GeodesicComputationViaFastMorphing( aol::ParameterParser &Parser ) :
    // define the finest grid over the domain
    _gridSize( aol::Vec<ConfiguratorType::Dim,short>( ( 1 << Parser.getInt( "GridDepth" ) ) + 1 ) ),
    _grid( _gridSize ),
    // load number of time points and create vector of morphing objects
    _numberOfTimePoints( Parser.getInt( "numberOfTimePoints" ) ),
    _endObjects( Parser.getInt( "GridDepth" ), ConfiguratorType::Dim, 2 ),
    _objects( Parser.getInt( "GridDepth" ), ConfiguratorType::Dim, _numberOfTimePoints ),
    // load hyperelastic and other parameters
    _linearElasticEnergyDensity( Parser.getDouble( "stiffness" ) ),
    _elasticEnergyDensity( _linearElasticEnergyDensity ),
    _weakMaterialStiffness( Parser.getDouble( "weakMaterialStiffness" ) ),
    _withMatching( Parser.getInt( "withMatching" ) ),
    _finestLevel( Parser.getInt( "GridDepth" ) ),
    // create the matching displacements (the zeroth is just a dummy) and the correction displacements
    _morphingDisplacements( _finestLevel, ConfiguratorType::Dim, _numberOfTimePoints ),
    _correctionDisplacements( _finestLevel, ConfiguratorType::Dim, _numberOfTimePoints ) {

    // read in the directory, where results are to be saved
    Parser.getString( "destDirectory", _destDirectory );

    // read in parameters
    aol::Vector<int> doPathReinitialization;
    Parser.getRealVec( "matchingWeight", _matchingWeight );
    Parser.getRealVec( "smoothingPixelWidth", _smoothingPixelWidth );
    Parser.getIntVec( "level", _level );
    Parser.getIntVec( "method", _method );
    Parser.getIntVec( "mode", _mode );
    Parser.getIntVec( "maxDescentSteps", _maxDescentSteps );
    Parser.getRealVec( "maxAccuracy", _maxAccuracy );
    Parser.getIntVec( "shapeReinitMode", _shapeReinitMode );
    Parser.getIntVec( "doPathReinitialization", doPathReinitialization );
    _doPathReinitialization.resize( doPathReinitialization.size() );
    for ( int i = 0; i < doPathReinitialization.size(); i++ )
      _doPathReinitialization.set( i, doPathReinitialization[i] != 0 );

    // load the end objects and initialization objects and scale them to the interval [0,1]
    loadCharFun( Parser.getString( "initialObject" ).c_str(), _endObjects.getArray( 0 ) );
    loadCharFun( Parser.getString( "finalObject" ).c_str(), _endObjects.getArray( 1 ) );
    if ( Parser.hasVariable( "initializationObjectFileNameTemplate" ) )
      for ( int i = 0; i < _numberOfTimePoints; i++ ) {
        char fileName[1024];
        sprintf( fileName, Parser.getString( "initializationObjectFileNameTemplate" ).c_str(), i );
        loadCharFun( fileName, _objects.getArray( i ) );
      }
    else
      for ( int i = 0; i < _numberOfTimePoints; i++ )
        _objects.getArray( i ) = _endObjects.getArray( 1 );

    // initialise matching displacements
    if ( Parser.hasVariable( "initializationDisplacementFileNameTemplate" ) )
      for ( int i = 1; i < _numberOfTimePoints; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          char fileName[1024];
          sprintf( fileName, Parser.getString( "initializationDisplacementFileNameTemplate" ).c_str(), i, j );
          cerr<<"load " + string( fileName ) + "...\n";
          _morphingDisplacements.getArray(i)[j].load( fileName );
        }

    // generate the needed coarser scales of objects and displacements
    if ( _level.getMinValue() < _finestLevel ) {
      _endObjects.levRestrict( _level.getMinValue(), _finestLevel );
      _objects.levRestrict( _level.getMinValue(), _finestLevel );
      _morphingDisplacements.levRestrict( _level.getMinValue(), _finestLevel );
      _correctionDisplacements.levRestrict( _level.getMinValue(), _finestLevel );
    }

    // fix number of needed threads
    omp_set_num_threads( _numberOfTimePoints );

    // notify the user
    cerr << "Parameters loaded." << endl;
  }

  template<typename ConfType, typename VecType, typename MatType, typename SubMatType>
  void minimize( const typename ConfType::InitType &Grid,
                 const aol::Op<VecType,aol::Scalar<typename ConfType::RealType> > &E,
                 const aol::Op<VecType> &DE,
                 const aol::Op<VecType,MatType> &D2E,
                 VecType &Minimizer,
                 const short Method,
                 const int MaxDescentSteps,
                 const typename ConfType::RealType MaxAccuracy,
                 const int Mode,
                 MatType* PMatrix /* will be deleted in this method */ ) const {
    typedef typename ConfType::RealType RealType;
    typedef typename aol::DiagonalPreconditioner<VecType> PreconditionType;
    // typedef typename aol::SSORPreconditioner<VecType,SubMatType> PreconditionType;

    VecType initialGuess( Minimizer );

    switch ( Method ) { //alternative minimization approaches
    case 1: { // gradient descent
              // fastest and good matching, but (for same accuracy) higher energy
      typedef aol::GradientDescent<ConfiguratorType,VecType> GradDescType;
      GradDescType minimizer( Grid, E, DE, MaxDescentSteps, 1., MaxAccuracy );
      minimizer.setConfigurationFlags( GradDescType::USE_NONLINEAR_CG | GradDescType::DO_NOT_SMOOTH_DESCENT_DIRECTION | GradDescType::USE_GRADIENT_BASED_STOPPING );
      minimizer.apply( initialGuess, Minimizer );
      delete PMatrix;
    } break;
    case 2: { // DFP Quasi-Newton method (Mode = #iterations before reset)
      aol::QuasiNewtonIteration<RealType,aol::MultiVector<RealType>,aol::MultiVector<RealType> > minimizer( E, DE, MaxDescentSteps, MaxAccuracy, Mode, false );
      minimizer.setTauMin( 1.e-5 );
      minimizer.setBeta( .9 );
      minimizer.setTimestepController( aol::NewtonInfo<RealType>::WOLFE );
      minimizer.apply( initialGuess, Minimizer );
      delete PMatrix;
    } break;
    case 3: { // Newton iteration
      NewtonMinimizer<RealType,VecType,MatType,false,PreconditionType,SubMatType> minimizer( E, DE, D2E, PMatrix, MaxDescentSteps, MaxAccuracy );
      minimizer.setTimestepController( aol::NewtonInfo<RealType>::ARMIJO );
      minimizer.apply( initialGuess, Minimizer );
    } break;
    case 4: { // trust region method (Mode = iterative or not)
              // noniterative works up to level 6, iterative is faster (despite 100 times more iterations) and yields less energy but slightly worse matching
      aol::TrustRegionMethod<RealType,aol::MultiVector<RealType>,aol::BlockOpBase<RealType,MatrixType>,PreconditionType,MatrixType> minimizer( E, DE, D2E, PMatrix, MaxDescentSteps, MaxAccuracy );
      minimizer.setIterativeMode( Mode );
      minimizer.apply( initialGuess, Minimizer );
    } break;
    case 5: { // BFGS Quasi-Newton method (Mode = #iterations before reset)
      aol::QuasiNewtonBFGS<RealType,aol::MultiVector<RealType>,aol::MultiVector<RealType> > minimizer( E, DE, MaxDescentSteps, MaxAccuracy, Mode, false );
      minimizer.setTimestepController( aol::NewtonInfo<RealType>::WOLFE );
      minimizer.apply( initialGuess, Minimizer );
      delete PMatrix;
    } break;
    default: throw aol::Exception( "Minimization method not provided" );
    }
  }

  void execute() {
    for ( int iter = 0; iter < _level.size(); iter++ ) {
      cerr << endl << aol::color::red << "iteration " << iter + 1 << " of " << _level.size() << aol::color::reset << endl;
      //! print parameters
      cerr << " level        " << _level[iter] << " of " << _finestLevel << endl;
      cerr << " method       " << _method[iter] << endl;
      cerr << " min mode     " << _mode[iter] << endl;
      cerr << " reinit       " << _doPathReinitialization[iter] << endl;
      cerr << " shape transp " << _shapeReinitMode[iter] << endl;

      //! define mesh
      const qc::GridSize<ConfiguratorType::Dim> gridSize( aol::Vec<ConfiguratorType::Dim,short>( ( 1 << _level[iter] ) + 1 ) );
      const typename ConfiguratorType::InitType grid( gridSize );
      _endObjects.setCurLevel( _level[iter] );
      _objects.setCurLevel( _level[iter] );
      _morphingDisplacements.setCurLevel( _level[iter] );
      _correctionDisplacements.setCurLevel( _level[iter] );

      //! put variables into a more convenient form
      aol::RandomAccessContainer< qc::MultiArray<RealType, ConfiguratorType::Dim> > morphingDisplacements;
      _morphingDisplacements.getArrayReferences( morphingDisplacements );

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 1; i < _numberOfTimePoints; i++ )
        resolveSelfPenetration<ConfiguratorType>( morphingDisplacements[i], grid, 100, 0, true );

      //! find optimal correction displacements
      // smooth the end objects
      qc::LinearConvolutionSmoothOp<RealType, ConfiguratorType::Dim> gaussKernelSmoother( grid.getNumX(), grid.getNumY() );
      qc::LinearConvolutionSmoothGradientOp<RealType, ConfiguratorType::Dim> gaussKernelSmootherGrad( grid.getNumX(), grid.getNumY() );
      gaussKernelSmoother.setSigma( _smoothingPixelWidth[iter] * grid.H() );
      gaussKernelSmootherGrad.setSigma( _smoothingPixelWidth[iter] * grid.H() );
      aol::Vector<RealType> initialObject( grid ), finalObject( grid ), initialObjectApprox( grid ), finalObjectApprox( grid );
      aol::MultiVector<RealType> initialObjectGrad( grid ), finalObjectGrad( grid );
      gaussKernelSmoother.apply( _endObjects.getArray( 0, _level[iter] ), initialObject );
      gaussKernelSmoother.apply( _endObjects.getArray( 1, _level[iter] ), finalObject );
      gaussKernelSmoother.apply( _objects.getArray( 0, _level[iter] ), initialObjectApprox );
      gaussKernelSmoother.apply( _objects.getArray( _numberOfTimePoints - 1, _level[iter] ), finalObjectApprox );
      gaussKernelSmootherGrad.apply( _endObjects.getArray( 0, _level[iter] ), initialObjectGrad );
      gaussKernelSmootherGrad.apply( _endObjects.getArray( 1, _level[iter] ), finalObjectGrad );
      // declare the energy and its variation
      aol::MultiVector<RealType> aux;
      _objects.appendReferencesTo( aux );
      aol::MultiVector<RealType> elastWeight( aux, aol::DEEP_COPY );
      elastWeight *= 1. - _weakMaterialStiffness;
      elastWeight.addToAll( _weakMaterialStiffness );
      MorphingCorrectionEnergy<ConfiguratorType,MatLawType>
        e( grid, morphingDisplacements, initialObject, finalObject, initialObjectApprox, finalObjectApprox, _elasticEnergyDensity, elastWeight, _matchingWeight[iter], _withMatching );
      MorphingCorrectionGradient<ConfiguratorType,MatLawType>
        de( grid, morphingDisplacements, initialObject, finalObject, initialObjectApprox, finalObjectApprox, _elasticEnergyDensity, elastWeight, _matchingWeight[iter], _withMatching );
      MorphingCorrectionHessian<ConfiguratorType,MatLawType,MatrixType>
        d2e( grid, morphingDisplacements, initialObject, finalObject, initialObjectApprox, finalObjectApprox, initialObjectGrad, finalObjectGrad, _elasticEnergyDensity, elastWeight, _matchingWeight[iter], _withMatching );
      // perform minimization
      aol::MultiVector<RealType> corrections;
      _correctionDisplacements.appendReferencesTo( corrections );
      aol::BlockMatrix<MatrixType>* pHessian = new aol::BlockMatrix<MatrixType>( _numberOfTimePoints * ConfiguratorType::Dim, _numberOfTimePoints * ConfiguratorType::Dim, grid.getNumberOfNodes(), grid.getNumberOfNodes() );
      minimize<ConfiguratorType,aol::MultiVector<RealType>,aol::BlockOpBase<RealType,MatrixType>,MatrixType>( grid, e, de, d2e, corrections, _method[iter], _maxDescentSteps[iter], _maxAccuracy[iter], _mode[iter], pHessian );

      //! change the corrections s. t. they have small average translation and rotation and remove boundary artifacts
      aol::MultiVector<RealType> identity( grid );
      qc::DataGenerator<ConfiguratorType>( grid ).generateIdentity( identity );
      for ( int i = 0; i < _numberOfTimePoints; i++ ) {
        aol::MultiVector<RealType> correction;
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          correction.appendReference( corrections[ConfiguratorType::Dim*i+j] );
        if ( i > 0 && i < _numberOfTimePoints - 1 ) {
          typename ConfiguratorType::VecType center(.5), translation;
          typename ConfiguratorType::MatType rotation;
          qc::projectOntoRigidDisplacement<ConfiguratorType>( correction, elastWeight[i], grid, center, translation, rotation, true );
          rotation = rotation.inverse();
          for ( int k = 0; k < correction[0].size(); k++ ){
            typename ConfiguratorType::VecType pos;
            for ( int j = 0; j < ConfiguratorType::Dim; j++ )
              pos[j] = identity[j][k] + correction[j][k];
            pos -= center;
            pos -= translation;
            pos = rotation * pos;
            for ( int j = 0; j < ConfiguratorType::Dim; j++ )
              correction[j][k] = pos[j] + center[j] - identity[j][k];
          }
        }
        aol::MultiVector<RealType> origDisp( correction, aol::DEEP_COPY );
        if ( _level[iter] > 3 )
          elasticallySmoothDeformationAtBoundary<ConfiguratorType,MatLawType>( origDisp, grid, _elasticEnergyDensity, correction, aol::Min( 4, 1 << ( _level[iter] - 4 ) ) );
      }

      //! prolongate the intermediate results onto finest grid
      _correctionDisplacements.levProlongate( _level[iter], _finestLevel );
      _morphingDisplacements.setCurLevel( _finestLevel );
      _objects.setCurLevel( _finestLevel );
      aol::RandomAccessContainer< qc::MultiArray<RealType, ConfiguratorType::Dim> > correctionDisplacementsFine, morphingDisplacementsFine;
      _correctionDisplacements.getArrayReferences( correctionDisplacementsFine );
      _morphingDisplacements.getArrayReferences( morphingDisplacementsFine );

      //! compute corrected morphing objects
      aol::MultiVector<RealType> corrObj( _numberOfTimePoints, _grid.getNumberOfNodes() );
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < _numberOfTimePoints; i++ )
        qc::InvDeformImage<ConfiguratorType, aol::Vector<RealType> > ( _objects.getArray( i, _finestLevel ), _grid, corrObj[i], correctionDisplacementsFine[i] );

      //! compute corrected morphing displacements
      aol::RandomAccessContainer< qc::MultiArray<RealType, ConfiguratorType::Dim> > corrMorphDisp( _numberOfTimePoints, morphingDisplacementsFine[0] );
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 1; i < _numberOfTimePoints; i++ ) {
        aol::MultiVector<RealType> auxDisp( _grid );
        qc::ConcatAndSmoothlyExtendDeformations<ConfiguratorType>( correctionDisplacementsFine[i], morphingDisplacementsFine[i], _grid, auxDisp );
        invConcatAndSmoothlyExtendDeformations<ConfiguratorType>( auxDisp, correctionDisplacementsFine[i-1], _grid, _elasticEnergyDensity, corrMorphDisp[i], false );
      }

      //! if desired, update morphing objects and displacements (after prolongation)
      switch ( _shapeReinitMode[iter] ) {
        case 2: 
          corrObj[0] = _endObjects.getArray( 0, _finestLevel );
          for ( int i = 1; i < _numberOfTimePoints; i++ )
            qc::InvDeformImage<ConfiguratorType, aol::Vector<RealType> > ( corrObj[i-1], _grid, corrObj[i], corrMorphDisp[i] );
          break;
        case 1:
          corrObj[_numberOfTimePoints-1] = _endObjects.getArray( 1, _finestLevel );
          for ( int i = _numberOfTimePoints - 1; i > 0; i-- )
            qc::DeformImage<ConfiguratorType> ( corrObj[i], _grid, corrObj[i-1], corrMorphDisp[i] );
          break;
        default:
          break;
      }
      if ( _doPathReinitialization[iter] )
        for ( int i = 0; i < _numberOfTimePoints; i++ ) {
          _objects.getArray( i, _finestLevel ) = corrObj[i];
          morphingDisplacementsFine[i] = corrMorphDisp[i];
          correctionDisplacementsFine[i].setZero();
        }
      if ( _level.getMinValue() < _finestLevel ) {
        _morphingDisplacements.levRestrict( _level.getMinValue(), _finestLevel );
        _objects.levRestrict( _level.getMinValue(), _finestLevel );
        _correctionDisplacements.levRestrict( _level.getMinValue(), _finestLevel );
      }

      //! save data
      for ( int i = 0; i < _numberOfTimePoints; i++ ) {
        char fileName[1024];
        // save morphing object
        ArrayType object( corrObj[i], _grid, aol::DEEP_COPY );
        object *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
        sprintf( fileName, "%s/object_%d.%s", _destDirectory, i, ConfiguratorType::Dim == 2 ? "pgm" : "bz2" );
        object.save( fileName, ConfiguratorType::Dim == 2 ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );
        // save morphing displacement
        saveMultiArray<ConfiguratorType>( corrMorphDisp[i], _grid, "%s/displacement_%d_%d.bz2", _destDirectory, i );
      }
    }
  }
};

template <typename ConfiguratorType>
class GeodesicShooting {

private:
  typedef typename ConfiguratorType::VecType   VecType;
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef aol::SparseMatrix<RealType> MatrixType;
  typedef qc::HyperelasticEnergyDensity2DLinear<ConfiguratorType> MaterialLawType;
  typedef qc::AntiCompressionHyperelasticEnergyDensity<ConfiguratorType,MaterialLawType> MatLawType;
  typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,MatLawType> WeightMatLawType;

  // finest grid over the examined domain
  const qc::GridSize<ConfiguratorType::Dim> _gridSize;
  const typename ConfiguratorType::InitType _grid;
  // number of time points along the morphing path
  const int _numberOfTimePoints, _numberOfSegments;
  // objects along the geodesic path
  aol::MultiVector<RealType> _objects;
  // matching displacements along geodesic path
  aol::RandomAccessContainer<aol::MultiVector<RealType> > _morphingDisplacements;
  // the directory for the results
  char _destDirectory[1024];
  // the hyperelastic deformation energy density
  const MaterialLawType _linearElasticEnergyDensity;
  const MatLawType _elasticEnergyDensity;
  const RealType _weakMaterialStiffness;
  // coarsest and finest level for the multiscale method
  const int _coarsestLevel, _finestLevel;
  // iteration number for Newton's method and desired accuracy
  const int _maxIterations;
  const RealType _maxAccuracy;

public:
  GeodesicShooting( aol::ParameterParser &Parser ) :
    // define the finest grid over the domain
    _gridSize( aol::Vec<ConfiguratorType::Dim,short>( ( 1 << Parser.getInt( "GridDepth" ) ) + 1 ) ),
    _grid( _gridSize ),
    // load number of time points and create vector of morphing objects and all needed displacements
    _numberOfTimePoints( Parser.getInt( "numberOfTimePoints" ) ),
    _numberOfSegments( _numberOfTimePoints - 1 ),
    _objects( _numberOfTimePoints, _grid.getNumberOfNodes() ),
    _morphingDisplacements( _numberOfTimePoints, aol::MultiVector<RealType>( _grid ) ),
    // load hyperelastic and other parameters
    _linearElasticEnergyDensity( Parser.getDouble( "stiffness" ) ),
    _elasticEnergyDensity( _linearElasticEnergyDensity ),
    _weakMaterialStiffness( Parser.getDouble( "weakMaterialStiffness" ) ),
    // load multiscale levels and number of Newton iterations
    _coarsestLevel( Parser.getInt( "coarsestLevel" ) ),
    _finestLevel( Parser.getInt( "GridDepth" ) ),
    _maxIterations( Parser.getInt( "maxIterations" ) ),
    _maxAccuracy( Parser.getDouble( "maxAccuracy" ) ) {

    // read in the directory, where results are to be saved
    Parser.getString( "destDirectory", _destDirectory );

    // load the first geodesic object and scale it to the interval [0,1]
    ArrayType object( _objects[0], _grid, aol::FLAT_COPY );
    loadCharFun( Parser.getString( "initObjectFileName" ).c_str(), object );
    char fileName[1024];
    ArrayType aux( object, aol::DEEP_COPY );
    aux *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
    sprintf( fileName, "%s/object_0.%s", _destDirectory, ConfiguratorType::Dim == 2 ? "pgm" : "bz2" );
    aux.save( fileName, ConfiguratorType::Dim == 2 ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );

    // load the first matching displacement
    aol::MultiVector<RealType> initDisp( _grid );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
      sprintf( fileName, Parser.getString( "initDisplacementFileNameTemplate" ).c_str(), j );
      cerr<<"load " + string( fileName ) + "...\n";
      ArrayType( initDisp[j], _grid, aol::FLAT_COPY ).load( fileName );
    }

    // remove possible boundary artifacts
    elasticallySmoothDeformationAtBoundary<ConfiguratorType,MatLawType>( initDisp, _grid, _elasticEnergyDensity, _morphingDisplacements[1] );
#ifdef _OPENMP
    // fix number of needed threads
    omp_set_num_threads( ConfiguratorType::Dim );
#endif
    // notify the user
    cerr << "Parameters loaded." << endl;
  }

  static void shootGeodesicOneTimestep( const typename ConfiguratorType::InitType &Grid,
                                        const MatLawType &MaterialLaw,
                                        const aol::Vector<RealType> &Weight,
                                        const aol::MultiVector<RealType> &DispHalf,
                                        aol::MultiVector<RealType> &DispFull,
                                        const int CoarseLevel = 3,
                                        const int MaxIterations = 1e3,
                                        const RealType MaxAccuracy = 1e-7 ) {
    //! initialize different levels of the multilevel method
    const int FineLevel = Grid.getGridDepth();
    qc::MultiDimMultilevelArray<RealType,ArrayType> tanVecDisp( Grid, DispFull.numComponents() ), twoStepDisp( Grid, DispFull.numComponents() );
    for ( int i = 0; i < DispHalf.numComponents(); i++ )
      tanVecDisp.getArray( i ) = DispHalf[i];
    qc::MultilevelArray<RealType,ArrayType> weight( Grid );
    weight.current() = Weight;
    if ( FineLevel > CoarseLevel ) {
      tanVecDisp.levRestrict( CoarseLevel, FineLevel );
      weight.levRestrict( CoarseLevel, FineLevel );
    }

    //! on each level of the multilevel method do...
    for ( int level = CoarseLevel; level <= FineLevel; level++ ) {
      cerr << endl << aol::color::red << "level " << level << " of " << FineLevel << aol::color::reset << endl;

      //! refine mesh
      typename ConfiguratorType::InitType grid( level, ConfiguratorType::Dim );
      tanVecDisp.setCurLevel( level );
      twoStepDisp.setCurLevel( level );

      //! put variables into a more convenient form
      aol::MultiVector<RealType> tanVec, twoStep;
      tanVecDisp.appendReferencesTo( tanVec );
      twoStepDisp.appendReferencesTo( twoStep );

      //! initialize the geodesic step
      if ( level == CoarseLevel ) {
        twoStep = tanVec;
        twoStep *= 2.;
        resolveSelfPenetration<ConfiguratorType>( twoStep, grid, 100, 0, true );
      }

      //! shoot one geodesic step on this scale
      WeightMatLawType materialLaw( MaterialLaw, grid, weight[level] );
      ParallelTransportOp<ConfiguratorType,WeightMatLawType>::shootGeodesic( grid, materialLaw, tanVec, twoStep, MaxIterations, MaxAccuracy );

      //! prolongate the result
      twoStepDisp.levProlongate();
    }

    //! return result
    for ( int i = 0; i < DispFull.numComponents(); i++ )
      DispFull[i] = twoStepDisp.getArray( i );
  }

  void execute() {
    cerr << "transport 0. shape onto 1. position..." << endl;
    qc::InvDeformAndSmoothlyExtendImage<ConfiguratorType>( _objects[0], _grid, _objects[1], _morphingDisplacements[1] );
    char fileName[1024];
    ArrayType object( _objects[1], _grid, aol::DEEP_COPY );
    object *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
    sprintf( fileName, "%s/object_1.%s", _destDirectory, ConfiguratorType::Dim == 2 ? "pgm" : "bz2" );
    object.save( fileName, ConfiguratorType::Dim == 2 ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );
    saveMultiArray<ConfiguratorType>( _morphingDisplacements[1], _grid, "%s/displacement_%d_%d.bz2", _destDirectory, 1 );

    for ( int k = 1; k < _numberOfSegments; k++ ) {
      cerr << "compute displacement from time step " << k-1 << " to " << k+1 << endl;
      aol::Vector<RealType> weight( _objects[k-1], aol::DEEP_COPY );
      weight.addToAll( _weakMaterialStiffness );
      weight.clamp( 0., 1. );
      aol::MultiVector<RealType> contDisp( _grid );
      shootGeodesicOneTimestep( _grid, _elasticEnergyDensity, weight, _morphingDisplacements[k], contDisp, _coarsestLevel, _maxIterations, _maxAccuracy );

      cerr << "compute displacement from time step " << k << " to " << k+1 << " (and resolve interpenetrations)" << endl;
      aol::MultiVector<RealType> auxDisp( _grid );
      invConcatAndSmoothlyExtendDeformations<ConfiguratorType>( contDisp, _morphingDisplacements[k], _grid, _elasticEnergyDensity, auxDisp, false );
      elasticallySmoothDeformationAtBoundary<ConfiguratorType,MatLawType>( auxDisp, _grid, _elasticEnergyDensity, _morphingDisplacements[k+1], 8 );

      cerr << "transport " << k << ". shape onto " << k+1 << ". position..." << endl;
      qc::InvDeformAndSmoothlyExtendImage<ConfiguratorType>( _objects[k], _grid, _objects[k+1], _morphingDisplacements[k+1] );

      cerr << "save result" << endl;
      // save morphing object
      ArrayType object( _objects[k+1], _grid, aol::DEEP_COPY );
      object *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
      sprintf( fileName, "%s/object_%d.%s", _destDirectory, k+1, ConfiguratorType::Dim == 2 ? "pgm" : "bz2" );
      object.save( fileName, ConfiguratorType::Dim == 2 ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );
      // save matching displacement
      saveMultiArray<ConfiguratorType>( _morphingDisplacements[k+1], _grid, "%s/displacement_%d_%d.bz2", _destDirectory, k+1 );
    }
  }
};

template <typename ConfiguratorType>
class ParallelTransport {

private:
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef aol::SparseMatrix<RealType> MatrixType;
  typedef qc::HyperelasticEnergyDensity2DLinear<ConfiguratorType> MaterialLawType;
  typedef qc::AntiCompressionHyperelasticEnergyDensity<ConfiguratorType,MaterialLawType> MatLawType;
  typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,MatLawType> WeightMatLawType;

  // finest grid over the examined domain
  const qc::GridSize<ConfiguratorType::Dim> _gridSize;
  const typename ConfiguratorType::InitType _grid;
  // number of time points along the morphing path
  const int _numberOfTimePoints;
  // given objects along the geodesic path
  aol::MultiVector<RealType> _objects;
  // matching displacements along geodesic path and parallel transported displacements of each object along the path
  aol::RandomAccessContainer<aol::MultiVector<RealType> > _morphingDisplacements, _transportedDisplacements;
  // the directory for the results
  char _destDirectory[1024];
  // the hyperelastic deformation energy density
  const MaterialLawType _linearElasticEnergyDensity;
  const MatLawType _elasticEnergyDensity;
  const RealType _weakMaterialStiffness;
  // coarsest and finest level for the multiscale method
  const int _coarsestLevel, _finestLevel;
  // iteration number for Newton's method and desired accuracy
  const int _maxIterations;
  const RealType _maxAccuracy;

public:
  ParallelTransport( aol::ParameterParser &Parser ) :
    // define the finest grid over the domain
    _gridSize( aol::Vec<ConfiguratorType::Dim,short>( ( 1 << Parser.getInt( "GridDepth" ) ) + 1 ) ),
    _grid( _gridSize ),
    // load number of time points and create vector of morphing objects and all needed displacements
    _numberOfTimePoints( Parser.getInt( "numberOfTimePoints" ) ),
    _objects( _numberOfTimePoints, _grid.getNumberOfNodes() ),
    _morphingDisplacements( _numberOfTimePoints, aol::MultiVector<RealType>( _grid ) ),
    _transportedDisplacements( _numberOfTimePoints, aol::MultiVector<RealType>( _grid ) ),
    // load hyperelastic and other parameters
    _linearElasticEnergyDensity( Parser.getDouble( "stiffness" ) ),
    _elasticEnergyDensity( _linearElasticEnergyDensity ),
    _weakMaterialStiffness( Parser.getDouble( "weakMaterialStiffness" ) ),
    // load multiscale levels and number of Newton iterations
    _coarsestLevel( Parser.getInt( "coarsestLevel" ) ),
    _finestLevel( Parser.getInt( "GridDepth" ) ),
    _maxIterations( Parser.getInt( "maxIterations" ) ),
    _maxAccuracy( Parser.getDouble( "maxAccuracy" ) ) {

    // read in the directory, where results are to be saved
    Parser.getString( "destDirectory", _destDirectory );

    // load the geodesic objects and scale them to the interval [0,1]
    for ( int i = 0; i < _numberOfTimePoints; i++ ) {
      char fileName[1024];
      sprintf( fileName, Parser.getString( "geodesicObjectFileNameTemplate" ).c_str(), i );
      ArrayType object( _objects[i], _grid, aol::FLAT_COPY );
      loadCharFun( fileName, object );
    }

    // load matching displacements
    for ( int i = 1; i < _numberOfTimePoints; i++ ) {
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        char fileName[1024];
        sprintf( fileName, Parser.getString( "geodesicDisplacementFileNameTemplate" ).c_str(), i, j );
        cerr<<"load " + string( fileName ) + "...\n";
        ArrayType( _morphingDisplacements[i][j], _grid, aol::FLAT_COPY ).load( fileName );
      }
    }

    // load displacement to be transported
    for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
      char fileName[1024];
      sprintf( fileName, Parser.getString( "tangentVectorDisplacement" ).c_str(), j );
      cerr<<"load " + string( fileName ) + "...\n";
      ArrayType( _transportedDisplacements[0][j], _grid, aol::FLAT_COPY ).load( fileName );
    }

    // save initial perturbed object
    saveMultiArray<ConfiguratorType>( _transportedDisplacements[0], _grid, "%s/tanVecDisp_%d_%d.bz2", _destDirectory, 0 );
    char fileName[1024];
    ArrayType pertObject( _grid );
    qc::InvDeformImage<ConfiguratorType,aol::Vector<RealType> > ( _objects[0], _grid, pertObject, _transportedDisplacements[0] );
    pertObject *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
    sprintf( fileName, "%s/tanVecObject_%d.%s", _destDirectory, 0, ConfiguratorType::Dim == 2 ? "pgm" : "bz2" );
    pertObject.save( fileName, ConfiguratorType::Dim == 2 ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );
#ifdef _OPENMP
    // fix number of needed threads
    omp_set_num_threads( ConfiguratorType::Dim );
#endif
    // notify the user
    cerr << "Parameters loaded." << endl;
  }

  static void computeParallelogramMidpoint( const typename ConfiguratorType::InitType &Grid,
                                            const MatLawType &MaterialLaw,
                                            const aol::Vector<RealType> &Weight,
                                            const aol::MultiVector<RealType> &TanVecDisp,
                                            const aol::MultiVector<RealType> &GeodesicDisp,
                                            aol::MultiVector<RealType> &MidPointDisp,
                                            const int CoarseLevel = 3,
                                            const int MaxIterations = 1e3,
                                            const RealType MaxAccuracy = 1e-7 ) {
    //! initialize different levels of the multilevel method
    const int FineLevel = Grid.getGridDepth();
    qc::MultiDimMultilevelArray<RealType,ArrayType> tanVecDisp( Grid, TanVecDisp.numComponents() ),
                                                    geodesicDisp( Grid, GeodesicDisp.numComponents() ),
                                                    midPointDisp( Grid, MidPointDisp.numComponents() );
    for ( int i = 0; i < TanVecDisp.numComponents(); i++ ) {
      tanVecDisp.getArray( i ) = TanVecDisp[i];
      geodesicDisp.getArray( i ) = GeodesicDisp[i];
    }
    qc::MultilevelArray<RealType,ArrayType> weight( Grid );
    weight.current() = Weight;
    if ( FineLevel > CoarseLevel ) {
      tanVecDisp.levRestrict( CoarseLevel, FineLevel );
      geodesicDisp.levRestrict( CoarseLevel, FineLevel );
      weight.levRestrict( CoarseLevel, FineLevel );
    }

    //! on each level of the multilevel method do...
    for ( int level = CoarseLevel; level <= FineLevel; level++ ) {
      cerr << endl << aol::color::red << "level " << level << " of " << FineLevel << aol::color::reset << endl;

      //! refine mesh
      typename ConfiguratorType::InitType grid( level, ConfiguratorType::Dim );
      tanVecDisp.setCurLevel( level );
      geodesicDisp.setCurLevel( level );
      midPointDisp.setCurLevel( level );

      //! put variables into a more convenient form
      aol::MultiVector<RealType> tanVec, geodesic, midPoint;
      tanVecDisp.appendReferencesTo( tanVec );
      geodesicDisp.appendReferencesTo( geodesic );
      midPointDisp.appendReferencesTo( midPoint );
      resolveSelfPenetration<ConfiguratorType>( geodesic, grid, 100, 0, true );

      //! initialize the mid point displacement
      if ( level == CoarseLevel ) {
        midPoint = tanVec;
        midPoint += geodesic;
        midPoint /= 2.;
        resolveSelfPenetration<ConfiguratorType>( midPoint, grid, 100, 0, true );
      }

      //! shoot one geodesic step on this scale
      WeightMatLawType materialLaw( MaterialLaw, grid, weight[level] );
      ParallelTransportOp<ConfiguratorType,WeightMatLawType>::computeParallelogramMidpoint( grid, materialLaw, tanVec, geodesic, midPoint, MaxIterations, MaxAccuracy );

      //! prolongate the result
      midPointDisp.levProlongate();
    }

    //! return result
    for ( int i = 0; i < MidPointDisp.numComponents(); i++ )
      MidPointDisp[i] = midPointDisp.getArray( i );
  }

  void execute() {
    for ( int k = 1; k < _numberOfTimePoints; k++ ) {
      cerr << "compute parallel transport from time step " << k-1 << " to " << k << ":" << endl;
      aol::Vector<RealType> weight( _objects[k-1], aol::DEEP_COPY );
      weight.addToAll( _weakMaterialStiffness );
      weight.clamp( 0., 1. );
      aol::MultiVector<RealType> midDisp( _grid ), diagDisp( _grid );

      cerr << "compute parallelogram midpoint" << endl;
      computeParallelogramMidpoint( _grid, _elasticEnergyDensity, weight, _transportedDisplacements[k-1], _morphingDisplacements[k], midDisp, _coarsestLevel, _maxIterations, _maxAccuracy );

      cerr << "compute parallelogram diagonal" << endl;
      GeodesicShooting<ConfiguratorType>::shootGeodesicOneTimestep( _grid, _elasticEnergyDensity, weight, midDisp, diagDisp, _coarsestLevel, _maxIterations, _maxAccuracy );

      cerr << "compute transported vector (and resolve interpenetrations)" << endl;
      invConcatAndSmoothlyExtendDeformations<ConfiguratorType>( diagDisp, _morphingDisplacements[k], _grid, _elasticEnergyDensity, _transportedDisplacements[k], false );
      
      cerr << "save result" << endl;
      // save tangent vector displacement
      saveMultiArray<ConfiguratorType>( _transportedDisplacements[k], _grid, "%s/tanVecDisp_%d_%d.bz2", _destDirectory, k );
      // save perturbed object
      char fileName[1024];
      ArrayType pertObject( _grid );
      qc::InvDeformImage<ConfiguratorType,aol::Vector<RealType> > ( _objects[k], _grid, pertObject, _transportedDisplacements[k] );
      pertObject *= ( ConfiguratorType::Dim == 2 ) ? 255 : 1;
      sprintf( fileName, "%s/tanVecObject_%d.%s", _destDirectory, k, ConfiguratorType::Dim == 2 ? "pgm" : "bz2" );
      pertObject.save( fileName, ConfiguratorType::Dim == 2 ? qc::PGM_UNSIGNED_CHAR_BINARY : qc::PGM_DOUBLE_BINARY );
    }
  }
};

typedef double RealType;
const qc::Dimension Dim = qc::QC_2D;
//typedef qc::QuocConfiguratorTraitMultiLin<RealType, Dim, aol::GaussQuadrature<RealType, Dim, 3> > ConfiguratorType;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, Dim, aol::Quadrature2D<RealType,aol::SimpsonQuadrature<RealType> > > ConfiguratorType;

int main( int argc, char *argv[] ) {

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
        switch ( parameterParser.getInt( "programType" ) ) {
        // computing shape geodesics between two given shapes
        case 1: ( GeodesicComputationViaFastMorphing<ConfiguratorType>( parameterParser ) ).execute(); break;
        // computing parallel transport of shapes
        case 2: ( ParallelTransport<ConfiguratorType>( parameterParser ) ).execute(); break;
        // computing shape geodesics via shooting
        case 3: ( GeodesicShooting<ConfiguratorType>( parameterParser ) ).execute(); break;
        }
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
