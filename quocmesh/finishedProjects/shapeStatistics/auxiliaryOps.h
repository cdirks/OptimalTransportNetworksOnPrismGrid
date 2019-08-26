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

#ifndef __AUXILIARYOPS_H
#define __AUXILIARYOPS_H

#include <aol.h>
#include <multilevelArray.h>
#include <vectorExtensions.h>
#include <preconditioner.h>
#include <convolution.h>
#include <generator.h>
#include <deformations.h>
#include <signedDistanceOp.h>

/**************************************************************
 * General auxiliary classes.
 **************************************************************/

template <typename ConfiguratorType>
inline void dilateAndClamp( const aol::Vector<typename ConfiguratorType::RealType> &CharacteristicFunc,
                            const typename ConfiguratorType::InitType &Grid,
                            aol::Vector<typename ConfiguratorType::RealType> &Result,
                            const typename ConfiguratorType::RealType NumPixelDilation = 3.5,
                            const typename ConfiguratorType::RealType ClampLower = 1e-3,
                            const typename ConfiguratorType::RealType ClampUpper = 1. ) {
  // singed distance ops work via marching/sweeping methods and are thus only appropriate for rectangular meshes
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::QuocConfiguratorTraitMultiLin<RealType, ConfiguratorType::Dim, aol::GaussQuadrature<RealType, ConfiguratorType::Dim, 3> > MultiLinConfType;
  aol::Vector<RealType> levelsetFunc( CharacteristicFunc, aol::DEEP_COPY );
  levelsetFunc.addToAll( -.5 );
  typename MultiLinConfType::InitType grid( qc::GridSize<MultiLinConfType::Dim>::createFrom( Grid ) );
  ( typename qc::SignedDistanceOpTrait<MultiLinConfType,MultiLinConfType::Dim>::OpType( grid ) ).apply( levelsetFunc, Result );
  Result /= grid.H(); // now each pixel corresponds to an increase by 1
  Result.addToAll( NumPixelDilation );
  Result.clamp( ClampLower, ClampUpper );
}

template <typename ConfiguratorType, class ProlongOpType = qc::ProlongOp<typename ConfiguratorType::RealType>, class RestrictOpType = qc::RestrictOp<typename ConfiguratorType::RealType, qc::STD_QUOC_RESTRICT> >
class MultiLevelMultiArray :
  public qc::MultiDimMultilevelArray<typename ConfiguratorType::RealType, typename ConfiguratorType::ArrayType, ProlongOpType, RestrictOpType > {

private:
  const int _numComp, _numSubComp;

public:
  MultiLevelMultiArray( const typename ConfiguratorType::InitType &Grid, const int NumComp, const int NumSubComp ) :
    qc::MultiDimMultilevelArray<typename ConfiguratorType::RealType, typename ConfiguratorType::ArrayType, ProlongOpType, RestrictOpType >( Grid, NumComp * NumSubComp ),
    _numComp( NumComp ),
    _numSubComp( NumSubComp ) {}

  MultiLevelMultiArray( const int FineDepth, const qc::Dimension Dim, const int NumComp, const int NumSubComp ) :
    qc::MultiDimMultilevelArray<typename ConfiguratorType::RealType, typename ConfiguratorType::ArrayType, ProlongOpType, RestrictOpType >( FineDepth, Dim, NumComp * NumSubComp ),
    _numComp( NumComp ),
    _numSubComp( NumSubComp ) {}

  void restrictMultiArray( const int I, const int CoarseLevel, const int FineLevel ) {
    for ( int j = 0; j < _numSubComp; j++ )
      this->getMultilevelArray( I * _numSubComp + j ).levRestrict( CoarseLevel, FineLevel );
  }

  void prolongateMultiArray( const int I ) {
    for ( int j = 0; j < _numSubComp; j++ )
      this->getMultilevelArray( I * _numSubComp + j ).levProlongate();
  }

  void setCurLevelMultiArray( const int I, const int Level ) {
    for ( int j = 0; j < _numSubComp; j++ )
      this->getMultilevelArray( I * _numSubComp + j ).setCurLevel( Level );
  }

  void appendMultiArrayReferenceTo( const int I, aol::MultiVector<typename ConfiguratorType::RealType> &Result ) {
    for ( int j = 0; j < _numSubComp; j++ )
      Result.appendReference( this->getArray( I * _numSubComp + j ), false );
  }

  void appendMultiArrayReferenceTo( const int I, const int level, aol::MultiVector<typename ConfiguratorType::RealType> &Result ) {
    for ( int j = 0; j < _numSubComp; j++ )
      Result.appendReference( this->getArray( I * _numSubComp + j, level ), false );
  }

  typename ConfiguratorType::ArrayType& getMultiArrayComp( const int I, const int J ) {
    return this->getArray( I * _numSubComp + J );
  }

  typename ConfiguratorType::ArrayType& getMultiArrayComp( const int I, const int J, const int Level ) {
    return this->getArray( I * _numSubComp + J, Level );
  }

  void getMultiArrays( aol::RandomAccessContainer<aol::MultiVector<typename ConfiguratorType::RealType> > &Result ) {
    if ( Result.size() != 0 )
      throw aol::Exception( "Non-empty result-RandomAccessContainer", __FILE__, __LINE__ );
    else {
      Result.reallocate( _numComp );
      for ( int i = 0; i < _numComp; i++ )
        for ( int j = 0; j < _numSubComp; j++ )
          Result[i].appendReference( this->getArray( i * _numSubComp + j ), false );
    }
  }
};

template <typename ConfiguratorType>
class AmbrosioTortorelliPhasefieldOp :
  public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const RealType _beta, _nu, _epsilon;
  typename ConfiguratorType::MatrixType* _regularizationMat;
  typename ConfiguratorType::ArrayType* _regularizationRHS;
  const aol::OperatorType _opType;
  const bool _quiet;

public:
  AmbrosioTortorelliPhasefieldOp( const typename ConfiguratorType::InitType &Grid,
                                  const RealType Beta,
                                  const RealType Nu,
                                  const RealType Epsilon,
                                  const aol::OperatorType OpType = aol::ONTHEFLY,
                                  const bool Quiet = true ) :
    _grid( Grid ),
    _beta( Beta ),
    _nu( Nu ),
    _epsilon( Epsilon ),
    _regularizationMat( NULL ),
    _regularizationRHS( NULL ),
    _opType( OpType ),
    _quiet( Quiet ) {
    if ( OpType == aol::ASSEMBLED ) {
      _regularizationMat = new typename ConfiguratorType::MatrixType( _grid );
      aol::MassOp<ConfiguratorType>( _grid ).assembleAddMatrix( *_regularizationMat, _nu / ( 2. * _epsilon ) );
      aol::StiffOp<ConfiguratorType>( _grid ).assembleAddMatrix( *_regularizationMat, 2. * _nu * _epsilon );
      _regularizationRHS = new typename ConfiguratorType::ArrayType( _grid );
      aol::Vector<RealType> aux( _grid );
      aux.setAll( _nu / ( 2. * _epsilon ) );
      aol::MassOp<ConfiguratorType>( _grid ).apply( aux , *_regularizationRHS );
    }
  }

  ~AmbrosioTortorelliPhasefieldOp() {
    delete _regularizationMat;
    delete _regularizationRHS;
  }

  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    //! assemble the system matrix
    typename ConfiguratorType::MatrixType systemMatrix( _grid );
    ( aol::SquaredDiffWeightMassOp<ConfiguratorType>( _grid, Arg ) ).assembleAddMatrix( systemMatrix, _beta );
    if ( _opType == aol::ASSEMBLED )
      systemMatrix += *_regularizationMat;
    else {
      ( aol::MassOp<ConfiguratorType>( _grid ) ).assembleAddMatrix( systemMatrix, _nu / ( 2. * _epsilon ) );
      ( aol::StiffOp<ConfiguratorType>( _grid ) ).assembleAddMatrix( systemMatrix, 2. * _nu * _epsilon );
    }

    //! assemble the rhs
    aol::Vector<RealType> rhs( _grid );
    if ( _opType == aol::ASSEMBLED )
      rhs = *_regularizationRHS;
    else {
      aol::Vector<RealType> aux( _grid );
      aux.setAll( _nu / ( 2. * _epsilon ) );
      ( aol::MassOp<ConfiguratorType>( _grid ) ).apply( aux , rhs );
    }

    //! solve for the phase field
    aol::DiagonalPreconditioner< aol::Vector<RealType> > preCond( systemMatrix );
    aol::PCGInverse<aol::Vector<RealType> > inv( systemMatrix, preCond, 1e-16, 1000 );
    //aol::CGInverse<aol::Vector<RealType> > inv( systemMatrix, info );
    //inv.getSolverInfoReference().setQuietMode( true );
    inv.setStopping( aol::STOPPING_ABSOLUTE );
    inv.setQuietMode( _quiet );
    inv.applyAdd( rhs, Dest );
  }
};

template <typename ConfiguratorType>
class AmbrosioTortorelliImageSmootherOp :
  public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_phasefield;
  const RealType _alpha, _beta;

public:
  AmbrosioTortorelliImageSmootherOp( const typename ConfiguratorType::InitType &Grid,
                                     const aol::Vector<RealType> &Phasefield,
                                     const RealType Alpha,
                                     const RealType Beta ) :
    _grid( Grid ),
    _phasefield( Phasefield ),
    _alpha( Alpha ),
    _beta( Beta ) {}

  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    //! assemble the system matrix
    typename ConfiguratorType::MatrixType systemMatrix( _grid );
    aol::MassOp<ConfiguratorType>( _grid ).assembleAddMatrix( systemMatrix );
    aol::SquaredWeightStiffOp<ConfiguratorType>( _grid, _phasefield ).assembleAddMatrix( systemMatrix, _beta / _alpha );

    //! assemble the rhs
    aol::Vector<RealType> rhs( _grid );
    aol::MassOp<ConfiguratorType>( _grid ).apply( Arg, rhs );

    //! solve for the smoothed image
    aol::DiagonalPreconditioner< aol::Vector<RealType> > preCond( systemMatrix );
    aol::PCGInverse<aol::Vector<RealType> > inv( systemMatrix, preCond, 1e-16, 1000 );
    //aol::CGInverse<aol::Vector<RealType> > inv( systemMatrix );
    inv.setStopping( aol::STOPPING_ABSOLUTE );
    inv.applyAdd( rhs, Dest );
  }
};

template <typename ConfiguratorType>
class WeightedConvolutionOp :
  public aol::Op<typename ConfiguratorType::ArrayType> {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  const ArrayType &_weight;
  const ArrayType &_kernel;
  const qc::Convolution<ConfiguratorType::Dim,RealType> _convolution;
  ArrayType _convolvedWeight;

  inline typename aol::VecDimTrait<int,ConfiguratorType::Dim>::VecType size( const ArrayType &Vector ) {
    typename aol::VecDimTrait<int,ConfiguratorType::Dim>::VecType result;
    result[0] = Vector.getNumX();
    result[1] = Vector.getNumY();
    if ( ConfiguratorType::Dim == qc::QC_3D )
      result[2] = Vector.getNumZ();
    return result;
  }

public:
  WeightedConvolutionOp( const ArrayType &Weight, const ArrayType &Kernel ) :
    _weight( Weight ),
    _kernel( Kernel ),
    _convolution( size( Kernel ) ),
    _convolvedWeight( Weight, aol::STRUCT_COPY ) {
    _convolution.convolve( _weight, _kernel, _convolvedWeight );
    for ( int i = 0; i < _convolvedWeight.size(); i++ )
      if ( aol::Abs( _convolvedWeight[i] ) < 1.e-10 )
        _convolvedWeight[i] = 1.;
  }

  void apply( const ArrayType &Arg, ArrayType &Dest ) const {
    // multiply the argument pointwise with the weight
    ArrayType arg( Arg, aol::STRUCT_COPY );
    for ( int i = 0; i < arg.size(); i++ )
      arg[i] = Arg[i] * _weight[i];
    // convolve the result with the kernel and scale it with the convolved weight
    _convolution.convolve( arg, _kernel, Dest );
    for ( int i = 0; i < Dest.size(); i++ )
      Dest[i] /= _convolvedWeight[i];

    /* alternatively:
    // use qc::generate...Kernel(...) (but where is kernel center?)
    // 3d and 2d work with same code, since 2d-Array is just 3d-Array with z-coordinate 1
    for ( int z = 0; z < Arg.getNumZ(); z++ )
      for ( int y = 0; y < Arg.getNumY(); y++ )
        for ( int x = 0; x < Arg.getNumX(); x++ )
          qc::Array<RealType>( Dest, Arg.getNumX(), Arg.getNumY(), Arg.getNumZ(), aol::FLAT_COPY ).set( x, y, z, Arg.getWeightedConvolveValue( x, y, z, _kernel, _weight ) );
    */
  }

  void applyAdd( const ArrayType &Arg, ArrayType &Dest ) const {
    ArrayType dest( Dest, aol::STRUCT_COPY );
    apply( Arg, dest );
    Dest += dest;
  }
};

/**
 * \brief This class represents the weighted mass matrix \f$ \left(\int_\Omega w(\phi^{-1}(x))\varphi_i(x)\varphi_j(x)|det(\nabla\phi^{-1}(x))|dx\right)_{ij} \f$.
 *
 * Here, the \f$ \varphi \f$ represent the finite element base functions, and the domain \f$ \Omega \f$ is represented by the grid
 * which is passed to the constructor. Furthermore, \f$ d \f$ is passed to the constructor as aol::MultiVector, where \f$ d \f$ is
 * the displacement such that the deformation \f$ \phi \f$ is given by \f$ d \f$+identity. The weight \f$ w \f$ has to be implemented
 * in the derived class as \c getWeight().
 *
 * Same usage as aol::MassOp.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename Imp, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class InvDeformedWeightMassInterface :
  public aol::FELinScalarWeightedMassInterface<ConfiguratorType, InvDeformedWeightMassInterface<ConfiguratorType, Imp, IndexMode>, IndexMode> {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the displacement $d$ in all Cartesian directions
  const aol::DiscreteVectorFunctionDefault< ConfiguratorType, ConfiguratorType::Dim > _d;
  // the inverse displacement, $\phi^{-1}$-identity, in all Cartesian directions as vector of the dofs and as multilinear interpolation
  aol::MultiVector<RealType> _invDispDOFs;
  const aol::DiscreteVectorFunctionDefault< ConfiguratorType, ConfiguratorType::Dim > _invDisp;
  // array encoding whether $\phi^{-1}$ is defined at the corresponding node
  typename qc::BitArray<ConfiguratorType::Dim> _invPhiDefined;

public:
  InvDeformedWeightMassInterface( const typename ConfiguratorType::InitType &Grid,
                                  const aol::MultiVector<RealType> &D,
                                  aol::OperatorType OpType = aol::ONTHEFLY ) :
    // initialise the grid
    aol::FELinScalarWeightedMassInterface<ConfiguratorType, InvDeformedWeightMassInterface<ConfiguratorType, Imp, IndexMode>, IndexMode>( Grid, OpType ),
    _d( Grid, D ),
    _invDispDOFs( D, aol::STRUCT_COPY ),
    _invDisp( Grid, _invDispDOFs ),
    _invPhiDefined( Grid ) {
    // generate the identity function
    aol::MultiVector<RealType> identity( D, aol::STRUCT_COPY );
    qc::DataGenerator<ConfiguratorType>( Grid ).generateIdentity( identity );
    // compute $\phi^{-1}$-identity
    qc::TransformFunction<RealType, ConfiguratorType::Dim> transformOp( Grid );
    transformOp.setDeformation( D );
    transformOp.transform( identity, _invDispDOFs, _invPhiDefined );
    _invDispDOFs -= identity;
  }

  //! Returns \f$ (w\circ\phi^{-1})|det(\nabla\phi^{-1})| \f$ evaluated at the point specified by the element "El" and quadrature point "QuadPoint".
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::VecType& RefCoord ) const {
    // if $\phi^{-1}$ is not defined at this position, return 0
    if ( !_invPhiDefined.elementTrue( El ) )
      return 0.;

    // compute $\phi^{-1}(x)$
    typename ConfiguratorType::VecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    if ( !qc::transformCoord<ConfiguratorType> ( this->_grid, _invDisp, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord ) )
      return 0.;

    // compute $D\phi(\phi^{-1}(x))$
    typename ConfiguratorType::MatType dPhi;
    _d.evaluateGradient( transformedEl, transformedLocalCoord, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute $w(\phi^{-1}(x))|det(D\phi^{-1}(x))|$
    return this->asImp().getWeight( transformedEl, transformedLocalCoord ) / aol::Abs( dPhi.det() );
  }

  //! Returns \f$ w \f$, evaluated at specified point. Has to be implemented in the derived classes of the interface.
  inline RealType getWeight ( const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::VecType& RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().getWeight ( El, RefCoord );
  }

protected:
  //! Barton-Nackman trick
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  //! Barton-Nackman trick
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }
};

template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticEnergyEvaluator :
  public aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim,1,HyperelasticEnergyEvaluator<ConfiguratorType,HyperelasticEnergyDensityType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
public:
  HyperelasticEnergyEvaluator ( const qc::GridDefinition &Grid, const HyperelasticEnergyDensityType &HyperelasticEnergyDensity ) :
    aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim,1,HyperelasticEnergyEvaluator<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ) {};

  void getNonlinearity( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint,
                        const typename ConfiguratorType::VecType &/*RefCoord*/,
                        aol::Vec<1,RealType> &NL ) const {
    // the deformation gradient and its cofactor matrix
    typename ConfiguratorType::MatType dphi;
    typename ConfiguratorType::MatType cof;

    // compute the deformation gradient $\nabla\phi$
    DiscFuncs.evaluateGradientAtQuadPoint ( El, QuadPoint, dphi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dphi[i][i] += 1.;

    // compute the hyperelastic invariants $||\nabla\phi||_F^2,||cof(\nabla\phi)||_F^2,det(\nabla\phi)$
    RealType
      I1 = dphi.normSqr(),
      I2,
      I3 = dphi.det();
    if ( ConfiguratorType::Dim != 2 ) {
      cof.makeCofactorMatrix( dphi );
      I2 = cof.normSqr();
    } else
      // in 2d, $||cof(\nabla\phi)||_F$ equals $||\nabla\phi||_F$
      I2 = I1;
    NL[0] = _hyperelasticEnergyDensity.evaluateAtQuadPoint ( I1, I2, I3, El, QuadPoint );
  }
};



#endif
