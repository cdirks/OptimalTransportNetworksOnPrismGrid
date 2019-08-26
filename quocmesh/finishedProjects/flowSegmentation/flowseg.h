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

#ifndef __FLOWSEG_H
#define __FLOWSEG_H

#include <signedDistanceOp.h>
#include <ChanVese.h>
#include <ArmijoSearch.h>
#include <gradientflow.h>
#include <hyperelastic.h>
#include <deformations.h>
#include <quocTimestepSaver.h>
#include <chanVeseDescent.h>
#include <quocDescent.h>
#include <Newton.h>

/**
 * Class that collects everything the FlowField_Ab_Energy and its
 * variations have in common.
 *
 * \author Berkels
 *
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int FlowFieldDimension>
class FlowFieldCommons {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
  aol::auto_container< ConfiguratorType::Dim, const aol::DiscreteFunctionDefault<ConfiguratorType> > _discrFlowField;
public:
  FlowFieldCommons( const typename ConfiguratorType::InitType &Initializer,
                    const aol::MultiVector<RealType> &FlowField )
    : _grid( Initializer )
  {
    for ( int c=0; c<FlowFieldDimension; c++ ) {
      aol::DiscreteFunctionDefault<ConfiguratorType> temp( Initializer, FlowField[c] );
      _discrFlowField.set_copy( c, temp );
    }
  }
  ~FlowFieldCommons () {}

  RealType evaluateIndicator ( const aol::Vector<RealType> &IndicatorParameter,
                               const int /*IndicatorParameterIndex*/,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &RefCoord ) const {
    RealType indicator = 0.;
    aol::Vec<FlowFieldDimension, RealType> x;
    for ( int i = 0; i < FlowFieldDimension; i++ ){
      x[i]= (El[i]+RefCoord[i])*(_grid.H());
    }
    aol::Vec<FlowFieldDimension, RealType> Ax;
    for ( int i = 0; i < FlowFieldDimension; i++ ){ 
      for ( int j = 0; j < FlowFieldDimension; j++ ){ 
        Ax[i] += IndicatorParameter[FlowFieldDimension*i+j]  * x[j];
      }
    }
    for ( int i = 0; i < FlowFieldDimension; i++ ){  
      indicator += aol::Sqr( _discrFlowField[i].evaluateAtQuadPoint( El, QuadPoint)
                             - Ax[i] 
                             - IndicatorParameter[aol::Sqr(FlowFieldDimension)+i] );
    }
    return indicator;
  }
};


/*
 * Calulates the classical Chan Vese energy of vector valued images for multi levelset gray value segmentation
 *
 * \author Berkels, Olischlaeger
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int FlowFieldDimension>
class FlowField_Ab_Energy
: public aol::ChanVeseEnergyInterface< ConfiguratorType,
                                       HeavisideFunctionType,
                                       NumberOfLevelsetFunctions,
                                       FlowFieldDimension,
                                       FlowField_Ab_Energy<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension> >,
  public FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension> {
public:
  typedef typename ConfiguratorType::RealType RealType;

  FlowField_Ab_Energy( const typename ConfiguratorType::InitType &Initializer,
                       const HeavisideFunctionType &HeavisideFunction,
                       const aol::MultiVector<RealType> &FlowField )
  : aol::ChanVeseEnergyInterface< ConfiguratorType,
                                  HeavisideFunctionType,
                                  NumberOfLevelsetFunctions,
                                  FlowFieldDimension,
                                  FlowField_Ab_Energy<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension> >
                                ( Initializer, HeavisideFunction ),
    FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension> ( Initializer, FlowField ) {}

  ~FlowField_Ab_Energy () {}

  using FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension>::evaluateIndicator;
};

/*
 * Calulates the variation of FlowField_Ab_Energy with respect to the levelset functions.
 *
 * \author Berkels, Olischlaeger
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int FlowFieldDimension>
class FlowField_Ab_EnergyMultiPhiVariation
: public aol::ChanVeseEnergyLevelsetDerivativeInterface< ConfiguratorType,
                                                         HeavisideFunctionType,
                                                         NumberOfLevelsetFunctions,
                                                         FlowField_Ab_EnergyMultiPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension> >,
  public FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension> {
public:
  typedef typename ConfiguratorType::RealType RealType;

  FlowField_Ab_EnergyMultiPhiVariation( const typename ConfiguratorType::InitType &Initializer,
                                        const HeavisideFunctionType &HeavisideFunction,
                                        const aol::MultiVector<RealType> &FlowField,
                                        const aol::MultiVector<RealType> &Parameters )
  : aol::ChanVeseEnergyLevelsetDerivativeInterface< ConfiguratorType,
                                                    HeavisideFunctionType,
                                                    NumberOfLevelsetFunctions,
                                                    FlowField_Ab_EnergyMultiPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension> >
                                                   ( Initializer, Parameters, HeavisideFunction ),
    FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension> ( Initializer, FlowField ) {}

  ~FlowField_Ab_EnergyMultiPhiVariation () {}

  using FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension>::evaluateIndicator;
};


/*
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int FlowFieldDimension>
class FlowField_Ab_EnergyMultiSecondPhiVariation
: public aol::QuadraticCEMEnergySecondLevelsetDerivativeInterface< ConfiguratorType,
                                                                   HeavisideFunctionType,
                                                                   NumberOfLevelsetFunctions,
                                                                   FlowField_Ab_EnergyMultiSecondPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension> >,
  public FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension> {
public:
  typedef typename ConfiguratorType::RealType RealType;

  FlowField_Ab_EnergyMultiSecondPhiVariation ( const typename ConfiguratorType::InitType &Initializer,
                                               const HeavisideFunctionType &HeavisideFunction,
                                               const aol::MultiVector<RealType> &FlowField,
                                               const aol::MultiVector<RealType> &LevelsetFunctions,
                                               const aol::MultiVector<RealType> &Parameters,
                                               const int LevelsetFunctionNumber )
  : aol::QuadraticCEMEnergySecondLevelsetDerivativeInterface< ConfiguratorType,
                                                              HeavisideFunctionType,
                                                              NumberOfLevelsetFunctions,
                                                              FlowField_Ab_EnergyMultiSecondPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension>
                                                            > ( Initializer, LevelsetFunctions, Parameters, HeavisideFunction, LevelsetFunctionNumber ),
    FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension> ( Initializer, FlowField ) {}

  ~FlowField_Ab_EnergyMultiSecondPhiVariation () {}

  using FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension>::evaluateIndicator;
};

/*
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int FlowFieldDimension>
class FlowField_Ab_EnergyMultiMixedSecondPhiVariation
: public aol::QuadraticCEMEnergyMixedSecondLevelsetDerivativeInterface< ConfiguratorType,
                                                                        HeavisideFunctionType,
                                                                        NumberOfLevelsetFunctions,
                                                                        FlowField_Ab_EnergyMultiMixedSecondPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension> >,
  public FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension> {
public:
  typedef typename ConfiguratorType::RealType RealType;

  FlowField_Ab_EnergyMultiMixedSecondPhiVariation ( const typename ConfiguratorType::InitType &Initializer,
                                                    const HeavisideFunctionType &HeavisideFunction,
                                                    const aol::MultiVector<RealType> &FlowField,
                                                    const aol::MultiVector<RealType> &LevelsetFunctions,
                                                    const aol::MultiVector<RealType> &Parameters,
                                                    const int LevelsetFunctionNumberOne,
                                                    const int LevelsetFunctionNumberTwo )
  : aol::QuadraticCEMEnergyMixedSecondLevelsetDerivativeInterface< ConfiguratorType,
                                                                   HeavisideFunctionType,
                                                                   NumberOfLevelsetFunctions,
                                                                   FlowField_Ab_EnergyMultiMixedSecondPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension>
                                                                 > ( Initializer, LevelsetFunctions, Parameters, HeavisideFunction, LevelsetFunctionNumberOne, LevelsetFunctionNumberTwo ),
    FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension> ( Initializer, FlowField ) {}

  ~FlowField_Ab_EnergyMultiMixedSecondPhiVariation () {}

  using FlowFieldCommons<ConfiguratorType, HeavisideFunctionType, FlowFieldDimension>::evaluateIndicator;
};

/*
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int FlowFieldDimension, typename MatrixType>
class SecondPhiVaritionOfFullFlowField_Ab_Energy
: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::SparseBlockMatrix<MatrixType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::MultiVector<RealType> &_flowField;
  const aol::MultiVector<RealType> &_parameters;
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _deltaMCM, _lengthParameter;
public:

  SecondPhiVaritionOfFullFlowField_Ab_Energy ( const typename ConfiguratorType::InitType &Initializer,
                                               const HeavisideFunctionType &HeavisideFunction,
                                               const aol::MultiVector<RealType> &FlowField,
                                               const aol::MultiVector<RealType> &Parameters,
                                               const RealType DeltaMCM,
                                               const RealType LengthParameter )
    : _grid ( Initializer ),
      _flowField ( FlowField ),
      _parameters ( Parameters ),
      _heavisideFunction ( HeavisideFunction ),
      _deltaMCM ( DeltaMCM ),
      _lengthParameter ( LengthParameter )
  {}

  ~SecondPhiVaritionOfFullFlowField_Ab_Energy () {}

  virtual void apply( const aol::MultiVector<RealType> &MArg,  aol::SparseBlockMatrix<MatrixType> &MDest ) const{
    MDest.deleteBlockEntries();

    for ( int i = 0; i < NumberOfLevelsetFunctions; i++ ) {
      FlowField_Ab_EnergyMultiSecondPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension>
        DDEfidelity( _grid, _heavisideFunction, _flowField, MArg, _parameters, i );

      MatrixType &mat = MDest.allocateMatrix( i, i, _grid );
      DDEfidelity.assembleAddMatrix ( mat );

      aol::IsotropicMatrixStiffOp<ConfiguratorType> isotropicMatrixStiffOp ( _grid, MArg[i], _deltaMCM );
      isotropicMatrixStiffOp.assembleAddMatrix ( mat, _lengthParameter );

      for ( int j = 0; j < i; j++ ) {
        FlowField_Ab_EnergyMultiMixedSecondPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension>
          DDEfidelityMixed( _grid, _heavisideFunction, _flowField, MArg, _parameters, i, j );

        MatrixType &matOne = MDest.allocateMatrix( i, j, _grid );
        DDEfidelityMixed.assembleAddMatrix ( matOne );

        MatrixType &matTwo = MDest.allocateMatrix( j, i, _grid );
        matTwo = matOne;
      }
    }
  }

  virtual void applyAdd( const aol::MultiVector<RealType> &/*MArg*/, aol::SparseBlockMatrix<MatrixType> &/*MDest*/ ) const{
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }
};


/*
 * We need to know the number of degrees of freedom in a symmetric matrix at compile time.
 * Further we need to know numVolExvDofs = SymmMatDofsTrait<dim>::num + 2*dim + aol::Sqr(dim).
 *
 * I didn't find a better way to handle this than this trait.
 *
 * \author Berkels
 */
template <int Dim>
class SymmMatDofsTrait {};

template <>
class SymmMatDofsTrait<2> {
public:
  static const int num = 3;
  static const int numVolExvDofs = 3 + 2*2 + 2*2;
};

template <>
class SymmMatDofsTrait<3> {
public:
  static const int num = 6;
  static const int numVolExvDofs = 6 + 2*3 + 3*3;
};


enum VFIELDMODE {
  AFFINE,
  AFFINE_WITH_SKEW_PENALTY
};

template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int ImageDimension>
class VectorFieldAffineCoefficientsUpdater
{
private:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const HeavisideFunctionType &_heavisideFunction;
  aol::MultiVector<RealType> &_flowFieldAb;
  const aol::MultiVector<RealType> &_imageMVec;
  const VFIELDMODE _vFieldMode;
public:
  VectorFieldAffineCoefficientsUpdater ( const typename ConfiguratorType::InitType &Initializer,
                                         const HeavisideFunctionType &HeavisideFunction,
                                         const aol::MultiVector<RealType> &ImageMVec,
                                         aol::MultiVector<RealType> &FlowFieldAb,
                                         const VFIELDMODE VFieldMode = AFFINE )
    : _grid ( Initializer ),
      _heavisideFunction ( HeavisideFunction ),
      _flowFieldAb( FlowFieldAb ), // (A,b)
      _imageMVec( ImageMVec ), //VectorField V
      _vFieldMode ( VFieldMode )
  {
  }

  void update ( aol::MultiVector<RealType> &Position ) const {
    const int dim = ConfiguratorType::Dim;
    const int numberOfNodes = this->_grid.getNumberOfNodes();

    // Calculate the new (A,b) values based on the current levelset functions (stored in Position).

    // Calculate H
    aol::MultiLevelsetVolumes<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions>
      volE( this->_grid, Position, this->_heavisideFunction );

    aol::MultiVector<RealType> H_Vec( _flowFieldAb.numComponents(), 1 );
    volE.apply( H_Vec, H_Vec );

    // Calculate H_ij, H_i, V_i, V_ij
    aol::MultiVector<RealType> Identity( ImageDimension, numberOfNodes );
    qc::DataGenerator<ConfiguratorType> generator( this->_grid );
    generator.generateIdentity( Identity );

    aol::MultiVector<RealType> volExv( SymmMatDofsTrait<dim>::numVolExvDofs, numberOfNodes );

    for ( int i = 0; i < numberOfNodes; i++ ){
      // Prepare H_ij
      int indexOffset = 0;
      for ( int k = 0; k < dim; k++ ) {
        for ( int l = k; l >= 0; l-- ) {
          volExv[indexOffset][i]=Identity[k][i]*Identity[l][i];
          indexOffset++;
        }
      }

      // Prepare H_i
      for ( int k = 0; k < dim; k++ )
        volExv[SymmMatDofsTrait<dim>::num+k][i]=Identity[k][i];

      // Prepare V_i
      for ( int k = 0; k < dim; k++ )
        volExv[SymmMatDofsTrait<dim>::num+dim+k][i]=_imageMVec[k][i];

      // Prepare V_ij
      for ( int k = 0; k < dim; k++ )
        for ( int l = 0; l < dim; l++ )
          volExv[SymmMatDofsTrait<dim>::num+2*dim+k*dim+l][i]=Identity[l][i]*_imageMVec[k][i];
    }

    aol::MultiLevelsetWeightedVolumes<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, SymmMatDofsTrait<dim>::numVolExvDofs>
      volEweight( this->_grid, Position, this->_heavisideFunction, volExv );

    aol::MultiVector<RealType> H_index_Vec( _flowFieldAb.numComponents(), SymmMatDofsTrait<dim>::numVolExvDofs );
    volEweight.apply( H_index_Vec, H_index_Vec );

    /* 2D situation:
     * H_index_Vec[][0] ->H_00         ->x[0]*x[0]
     * H_index_Vec[][1] ->H_11         ->x[1]*x[1]
     * H_index_Vec[][2] ->H_01 = H_10  ->x[1]*x[0]
     * H_index_Vec[][3] ->H_0          ->x[0]
     * H_index_Vec[][4] ->H_1          ->x[1]
     * H_index_Vec[][5] ->V_0          ->v[0]
     * H_index_Vec[][6] ->V_1          ->v[1]
     * H_index_Vec[][7] ->V_00         ->x[0]*v[0]
     * H_index_Vec[][8] ->V_10         ->x[1]*v[0]
     * H_index_Vec[][9] ->V_01         ->x[0]*v[1]
     * H_index_Vec[][10] ->V_11        ->x[1]*v[1]
     */

    // solve for each component:
    for ( int i = 0; i < _flowFieldAb.numComponents(); i++ ){
      // This segment has zero volume, therefore we can't solve for it.
      // Attention: A comparision of a floating point value to zero is made here,
      // perhaps it's better only to do an approximate comparision.
      if ( H_Vec[i][0] == aol::NumberTrait<RealType>::zero )
        continue;

      aol::Mat<dim+1, dim+1, RealType> H_mat;
      
      /* System matrix in 2D:
       * H_0,  H_1,  H
       * H_00, H_01, H_0
       * H_01, H_11; H_1
       */

      // Set the H entry (top right).
      H_mat.set( 0, dim, H_Vec[i][0] );

      // Set the H_i vecs (first row and last column).
      for ( int k = 0; k < dim; k++ ){
        H_mat.set( 0,   k,   H_index_Vec[i][SymmMatDofsTrait<dim>::num+k] );
        H_mat.set( k+1, dim, H_index_Vec[i][SymmMatDofsTrait<dim>::num+k] );
      }

      // Set the H_ij matrix.
      int indexOffset = 0;
      // First set the lower left part of the matrix.
      for ( int k = 0; k < dim; k++ ) {
        for ( int l = k; l >= 0; l-- ) {
          H_mat.set( k+1, l, H_index_Vec[i][indexOffset] );
          indexOffset++;
        }
      }
      // Now copy the lower left non diag entries to the upper right (H_ij is symmetric).
      for ( int k = 0; k < dim; k++ ) {
        for ( int l = k+1; l < dim; l++ ) {
          H_mat.set( k+1, l, H_mat.get( l+1, k ) );
        }
      }

      std::vector<aol::Vec<dim+1, RealType> > argV ( dim );

      for ( int k = 0; k < dim; k++ )
      {
        argV[k][0] = H_index_Vec[i][SymmMatDofsTrait<dim>::num+dim+k];
        for ( int l = 0; l < dim; l++ )
          argV[k][l+1] = H_index_Vec[i][SymmMatDofsTrait<dim>::num+2*dim+l+k*dim];
      }

      std::vector<aol::Vec<dim+1, RealType> > destV ( dim );

     switch( _vFieldMode ) {
      case AFFINE:
        {
          aol::Mat<dim+1, dim+1, RealType> H_mat_inverse;
          H_mat_inverse.makeInverse( H_mat );

          for ( int k = 0; k < dim; k++ )
            H_mat_inverse.mult(argV[k], destV[k]);
        }
        break;

      case AFFINE_WITH_SKEW_PENALTY:
        {
          aol::Mat<dim*(dim+1), dim*(dim+1), RealType> systemMat;
          // Copy the H_mat matrix on the three diagonal block positions of systemMat.
          for ( int i = 0; i < dim; i++ )
            for ( int k = 0; k < dim+1; k++ )
              for ( int l = 0; l < dim+1; l++ )
                systemMat.set( k + i*(dim+1), l + i*(dim+1), H_mat.get ( k, l ) );

          for ( int k = 0; k < dim; k++ ) {
            for ( int l = k+1; l < dim; l++ ) {
              systemMat.add( k*(dim+1) + l+1, k*(dim+1) + l, 1 );
              systemMat.add( k*(dim+1) + l+1, l*(dim+1) + k, 1 );
              systemMat.add( l*(dim+1) + k+1, k*(dim+1) + l, 1 );
              systemMat.add( l*(dim+1) + k+1, l*(dim+1) + k, 1 );
            }
          }

          aol::Mat<dim*(dim+1), dim*(dim+1), RealType> systemMatInverse;
          systemMatInverse.makeInverse( systemMat );

          // Copy the argV vecs into rhs
          aol::Vec<dim*(dim+1), RealType> rhs;
          for ( int i = 0; i < dim; i++ )
            for ( int k = 0; k < dim+1; k++ )
              rhs[k + i*(dim+1)] = argV[i][k];

          aol::Vec<dim*(dim+1), RealType> dest;
          systemMatInverse.mult ( rhs, dest );

          // Copy the dest into the destV vecs
          for ( int i = 0; i < dim; i++ )
            for ( int k = 0; k < dim+1; k++ )
              destV[i][k] = dest[k + i*(dim+1)];

        }
        break;

      default:
        {
          throw aol::UnimplementedCodeException( "mode not implemented", __FILE__, __LINE__ );
        }
        break;
     }

      /* 2D situation:
       * destV[0][0] <=> _flowFieldAb[][0] <=> A_00
       * destV[0][1] <=> _flowFieldAb[][1] <=> A_01
       * destV[0][2] <=> _flowFieldAb[][4] <=> b_0
       *
       * destV[1][0] <=> _flowFieldAb[][2] <=> A_10
       * destV[1][1] <=> _flowFieldAb[][3] <=> A_11
       * destV[1][2] <=> _flowFieldAb[][5] <=> b_1
       */

      // Get A from destV
      for ( int k = 0; k < dim; k++ )
        for ( int l = 0; l < dim; l++ )
          _flowFieldAb[i][k*dim+l]= destV[k][l];

      // Get b from destV
      for ( int k = 0; k < dim; k++ )
        _flowFieldAb[i][aol::Sqr(dim)+k]= destV[k][dim];

    }
  }

  const aol::MultiVector<RealType>& getFlowFieldAb ( ) const {
    return _flowFieldAb;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int ImageDimension>
class VectorFieldGradientDescent
: public qc::ChanVeseGradientDescent<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> {
public:
private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::ChanVeseGradientDescent<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> SuperClass;
  const VectorFieldAffineCoefficientsUpdater<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ImageDimension> _coeffUpdater;
public:
  VectorFieldGradientDescent( const typename ConfiguratorType::InitType &Initializer,
                              const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                              const aol::Op<aol::MultiVector<RealType> > &DE,
                              const HeavisideFunctionType &HeavisideFunction,
                              const aol::MultiVector<RealType> &ImageMVec,
                              aol::MultiVector<RealType> &FlowFieldAb,
                              const aol::Vector<RealType> &VisualizationBackgroundImage,
                              const int MaxIterations = 50,
                              const RealType StartTau = aol::ZOTrait<RealType>::one,
                              const RealType StopEpsilon = aol::ZOTrait<RealType>::zero,
                              const VFIELDMODE VFieldMode = AFFINE )
    : SuperClass( Initializer, E, DE, HeavisideFunction, VisualizationBackgroundImage, MaxIterations, StartTau, StopEpsilon, SuperClass::MINIMIZE_WRT_LEVELSETFUNCTIONS ),
      _coeffUpdater ( Initializer, HeavisideFunction, ImageMVec, FlowFieldAb, VFieldMode )
  {
  }
  virtual bool postProcess( aol::MultiVector<RealType> &Position, const int Iteration ) const {
    bool postProcessed = SuperClass::postProcess( Position, Iteration );

    _coeffUpdater.update ( Position );

    // Output the timestep.
    this->writeMultiVectorToTextFile ( _coeffUpdater.getFlowFieldAb(), "flowFieldAb.txt", Iteration, this->_maxIterations );

    if ( this->checkSaveConditions ( Iteration ) ) {
      aol::MultiVector<RealType> currentVecField( ImageDimension, this->_grid.getNumberOfNodes() );
      qc::DataGenerator<ConfiguratorType> generator( this->_grid );
      generator.generatePiecewiseAffineFlowField ( currentVecField, Position, _coeffUpdater.getFlowFieldAb(), this->_displayedIsovalue );

      std::vector<aol::PlotOutFileType> outTypeVec;
      outTypeVec.push_back(aol::GNUPLOT_PNG);
      //outTypeVec.push_back(aol::GNUPLOT_PS);
      //outTypeVec.push_back(aol::GNUPLOT_EPS);
      this->plotVectorFieldSegmented ( currentVecField, Position, this->_grid, outTypeVec, Position.numComponents()+1, Iteration, true, this->_displayedIsovalue );
    }

    // This bool is only used for the stopping criterion to
    // prevent energy increases caused by the post processing
    // from stopping the iteration. Since the explicit minimization
    // with respect to the gray values never increases the energy,
    // we don't need to return a true if we just updated the
    // gray values.
    return postProcessed;
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int ImageDimension, typename MatrixType, typename SolverType>
class VectorFieldNewtonIteration
  : public aol::NewtonIterationSparseBlockMatrix<MatrixType, SolverType, NumberOfLevelsetFunctions>,
    public qc::QuocTimestepSaver<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;
  const VectorFieldAffineCoefficientsUpdater<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ImageDimension> _coeffUpdater;
  const typename ConfiguratorType::InitType &_grid;
public:
  VectorFieldNewtonIteration ( const typename ConfiguratorType::InitType &Initializer,
                               const aol::Op<aol::MultiVector<RealType> > &F,
                               const aol::Op<aol::MultiVector<RealType>, aol::SparseBlockMatrix<MatrixType> > &DF,
                               const HeavisideFunctionType &HeavisideFunction,
                               const aol::MultiVector<RealType> &ImageMVec,
                               aol::MultiVector<RealType> &FlowFieldAb,
                               const int MaxIterations = 50,
                               const RealType StopEpsilon = 1.e-6,
                               const VFIELDMODE VFieldMode = AFFINE )
    : aol::NewtonIterationSparseBlockMatrix<MatrixType, SolverType, NumberOfLevelsetFunctions>
        ( F, DF, MaxIterations, StopEpsilon ),
      qc::QuocTimestepSaver<ConfiguratorType>( 1, 1, "Me", true ),
      _coeffUpdater ( Initializer, HeavisideFunction, ImageMVec, FlowFieldAb, VFieldMode ),
      _grid ( Initializer )
  {
    this->_pInfo->getSolverInfo().setQuietMode ( true );
  }
  virtual void postProcess( aol::MultiVector<RealType> &Position, const int Iteration ) const {
    _coeffUpdater.update ( Position );

    // Output the timestep.
    writeMultiVectorToTextFile ( _coeffUpdater.getFlowFieldAb(), "flowFieldAb.txt", Iteration, this->_pInfo->getMaxIterations() );

    for ( int i = 0; i < NumberOfLevelsetFunctions; i++ ){
      this->saveTimestepBZ2( i, Iteration, Position[i], this->_grid );
      this->savePNG( Position[i], this->_grid, Iteration, i );
    }
  }

  virtual void writeResult( const aol::MultiVector<RealType> &/*Dest*/, const char* /*FileName*/ ) const {
  }
};

template <typename RealType>
class FlowParams{
public:
  //! gamma weights the perimiter length term
  const RealType gamma;
  //! epsilon is the regularization parameter of the Heaviside function
  const RealType epsilon;
  //! epsilon is the regularization parameter of the MCMStiffOp used in the variation of the length term.
  const RealType epsilonMCM;
  const int numSteps;
  //! outerSteps is the number of times epsilon will be halved
  const int outerSteps;
  FlowParams( aol::ParameterParser &Parser )
  : gamma( Parser.getDouble( "gamma" ) ),
    epsilon( Parser.getDouble( "epsilon" ) ),
    epsilonMCM( Parser.getDouble( "epsilonMCM" ) ),
    numSteps( Parser.getInt( "numSteps" ) ),
    outerSteps( Parser.getInt( "outerSteps" ) )
  {}
};

/*
 * Class to generate example segmentations. Currently VERY limited.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, int numberOfLevelsetFunctions>
class ExampleSegmentationGenerator : public qc::DataGenerator <ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType; 
  ExampleSegmentationGenerator ( const typename ConfiguratorType::InitType &Initializer )
    : qc::DataGenerator <ConfiguratorType> ( Initializer ) {}
  
  void generateExample ( aol::MultiVector<RealType> &LevelsetFunctions, const int ExampleNumber ) {
    switch ( ExampleNumber ) {
    case 0:
      {
        typename ConfiguratorType::VecType normal;
        normal[0] = 1.;
        normal[1] = 0.5;
        if ( ConfiguratorType::Dim == 3 )
          normal[2] = 0.5;
        normal /= normal.norm();
        typename ConfiguratorType::VecType point;
        point.setAll( 0.5 );
        this->generatePlaneLevelset( LevelsetFunctions[0], normal, point );
        typename ConfiguratorType::VecType center;
        center.setAll ( 0.5 );
        if ( numberOfLevelsetFunctions > 1 )
          this->generateSphereLevelset ( LevelsetFunctions[1], center, 0.3 );
      }
      break;

    default:
      throw aol::UnimplementedCodeException( "Unknown example number", __FILE__, __LINE__ );
      break;
    };
  }
};

//! The operator that controls all the necessary calls and data structures.
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int FlowFieldDimension = ConfiguratorType::Dim>
class VectorFieldSegmentation : public qc::QuocTimestepSaver<ConfiguratorType>
{
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  // ******** some more parameters ************
    
  typename ConfiguratorType::InitType _grid;
  aol::ParameterParser &_parser;
  const FlowParams<RealType> _params;
  const bool _usingThresholdModel;
  qc::DataGenerator<ConfiguratorType> _generator;
public:
  //! If GridDepth > 0, it overrides the grid depth setting in Parser.
  VectorFieldSegmentation( aol::ParameterParser &Parser, const int GridDepth = 0 ) :
    aol::TimestepSaver<RealType>( Parser.getInt("saveOffset"), Parser.getInt("numSaveFirst"), "PleaseRenameMe", true ),
    qc::QuocTimestepSaver<ConfiguratorType>( Parser.getInt("saveOffset"), Parser.getInt("numSaveFirst"), "PleaseRenameMe", true ),
    _grid( (GridDepth > 0) ? GridDepth : Parser.getInt( "depth" ) , ConfiguratorType::Dim ),
    _parser( Parser ),
    _params( Parser ),
    // Somewhat hackish way to find out, if we are using the Chan Vese or a threshold model.
    _usingThresholdModel( HeavisideFunctionType::ScalingFunctionType::evaluate == aol::SquareFunction<RealType>::evaluate ),
    _generator( _grid )
  { 
    //cerr<<"parser num save first"<<(Parser.getInt("numSaveFirst"))<<endl;
    string saveDirectory = aol::strprintf( "%s_%d_%f_%.2f_%d", Parser.getString( "saveDirectory" ).c_str(), qc::logBaseTwo ( _grid.getNumX() ), _params.gamma, _params.epsilon, _params.numSteps);
    aol::makeDirectory( saveDirectory.c_str() );
    saveDirectory += "/";
    this->setSaveDirectory( saveDirectory.c_str() );
    Parser.dumpToFile( "parameter-dump.txt", saveDirectory.c_str() );
  }

  virtual ~VectorFieldSegmentation() {
  }

  bool isUsingThresholdModel () const {
    return _usingThresholdModel;
  }

  //! Reads a vertor field from the file named FileName to VectorField (accepts different formats).
  //! The MultiArray is resized to the size of the read data set. 
  void readVectorField ( qc::MultiArray<RealType, ConfiguratorType::Dim, ConfiguratorType::Dim> &VectorField, const string &FileName ) {
    cerr << "Reading vectorfield from file " << FileName << endl;

    ifstream vfFile ( FileName.c_str(), ios::binary );

    if ( !vfFile.good() )
      throw aol::Exception ( "Couldn't open file for reading.", __FILE__, __LINE__ );

    // Assume data to be in vtk format.
    if ( ConfiguratorType::Dim == qc::QC_3D ) {
      string temp;

      aol::checkNextLineOrString ( vfFile, "# vtk DataFile Version 1.0" );
      aol::checkNextLineOrString ( vfFile, "vtk output" );
      aol::checkNextLineOrString ( vfFile, "ASCII" );
      aol::checkNextLineOrString ( vfFile, "DATASET STRUCTURED_POINTS" );
      aol::checkNextLineOrString ( vfFile, "DIMENSIONS", false );

      int inNumX, inNumY, inNumZ;
      vfFile >> inNumX >> inNumY >> inNumZ;
      vfFile.get();

      aol::checkNextLineOrString ( vfFile, "ASPECT_RATIO 1.0 1.0 1.0" );
      aol::checkNextLineOrString ( vfFile, "ORIGIN 0.0 0.0 0.0" );
      aol::checkNextLineOrString ( vfFile, "POINT_DATA", false );

      int numberOfPoints;
      vfFile >> numberOfPoints;
      if ( numberOfPoints != inNumX * inNumY * inNumZ )
        throw aol::Exception ( "Wrong entry in field POINT_DATA", __FILE__, __LINE__ );
      vfFile.get();

      aol::checkNextLineOrString ( vfFile, "" );
      aol::checkNextLineOrString ( vfFile, "VECTORS vectors float" );
      aol::checkNextLineOrString ( vfFile, "" );

      RealType dummy;

      VectorField.reallocate ( inNumX, inNumY, inNumZ );

      aol::ProgressBar<> pb ( "Progress" );
      pb.start ( numberOfPoints );

      for ( int j = 0; j < numberOfPoints; ++j ) {
        for ( int i = 0; i < 3; ++i ) {
          vfFile >> dummy;
          VectorField[i][j] = dummy;
        }
        pb++;
      }
      cerr << endl;
    }
    // Data in "benard" format.
    else if ( vfFile.peek() == 'V' ) {

      string temp;
      getline ( vfFile, temp );

      if ( temp.compare ( "V2" ) )
        throw aol::Exception ( "Unknown header.", __FILE__, __LINE__ );

      int inWidth, inHeight;

      vfFile >> inWidth;
      vfFile >> inHeight;

      cerr << "Input width = " << inWidth << ", height = " << inHeight << endl;

      RealType dummy;

      VectorField.reallocate ( inWidth, inHeight );

      for ( int j = 0; j < inWidth*inHeight; ++j ) {
        for ( int i = 0; i < 2; ++i ) {
          vfFile >> dummy;
          VectorField[i][j] = dummy;
        }
      }
    }
    else if ( _parser.hasVariable ( "grape2DVecInput" ) ) {
      const int inWidth = _parser.getInt ( "inNumX" );
      const int inHeight = _parser.getInt ( "inNumY" );
      VectorField.reallocate ( inWidth, inHeight );

      for ( int j = 0; j < inWidth*inHeight; ++j ) {
        for ( int i = 0; i < 2; ++i ) {
          VectorField[i][j] = aol::readBinaryData<double, RealType>( vfFile );
        }
      }
    }
    // Assume data to be in "rayleigh" format.
    else {
      /*const int xlength = static_cast<int>*/(aol::readBinaryData<float, float>( vfFile ));
      /*const int ylength = static_cast<int>*/(aol::readBinaryData<float, float>( vfFile ));
      const int inWidth = static_cast<int>(aol::readBinaryData<float, float>( vfFile ));
      const int inHeight = static_cast<int>(aol::readBinaryData<float, float>( vfFile ));
      cerr << "Input width = " << inWidth << ", height = " << inHeight << endl;
      if ( (inWidth <= 0) || (inHeight <= 0) )
        throw aol::Exception ( "Error reading \"rayleigh\" format, non-positive width or height\n", __FILE__, __LINE__ );
      VectorField.reallocate ( inWidth, inHeight );
      aol::readBinaryData<float, RealType>( vfFile, VectorField[0].getData(), inWidth*inHeight );
      aol::readBinaryData<float, RealType>( vfFile, VectorField[1].getData(), inWidth*inHeight );
    }

    if ( vfFile.good() )
      cerr << "Success.\n";
    else
      throw aol::Exception ( "Error reading from file", __FILE__, __LINE__ );

    vfFile.close();

    if ( ConfiguratorType::Dim == qc::QC_2D ) {
      string filename = aol::strprintf ( "%sinput.png", this->getSaveDirectory() );
      qc::writeColorField ( VectorField[0], VectorField[1], filename.c_str() );
    }
  }


  void initializeTargetVField( aol::MultiVector<RealType> &U0MVec ) {
    aol::MultiVector<RealType> flowFieldAb( 1<<NumberOfLevelsetFunctions, aol::Sqr(FlowFieldDimension)+FlowFieldDimension );
    if( _parser.hasVariable ( "Ab" ) ){
      aol::MultiVector<RealType> temp;
      _parser.getRealMultiVec<RealType> ( "Ab", temp );
      cerr << "Parsed Ab:\n" << temp << endl;
      flowFieldAb = temp;
    }
    else {
      flowFieldAb[0][1] = 0.1;
      flowFieldAb[0][2] = -0.1;
      flowFieldAb[0][4] = -0.1;
      flowFieldAb[0][5] = 0.1;

      //flowFieldAb[0][1] = -0.1;
      //flowFieldAb[0][2] = 0.1;
      flowFieldAb[1][4] = -0.1;
      flowFieldAb[1][5] = 0.1;
      /*
      flowFieldAb[2][1]=0.3;
      flowFieldAb[2][2]=-0.3;
      flowFieldAb[2][4]=-0.3;
      flowFieldAb[2][5]=0.3;

      flowFieldAb[3][4]=0.1;
      flowFieldAb[3][5]=-0.1;
      */
    }

    aol::MultiVector<RealType> levelsetFunctionsRef( NumberOfLevelsetFunctions, _grid.getNumberOfNodes() );
    //_generator.generateLineLevelset( levelsetFunctionsRef[0], 1, 0.5 );
    //_generator.generateChanVeseInitialization( levelsetFunctionsRef );
    _generator.generateCircleLevelset( levelsetFunctionsRef[0], 0.3, 0.5,0.5 );
    //_generator.generateLineLevelset( levelsetFunctionsRef[0], 0, 0.5 );
    //_generator.generateDiagonalLevelsetPlusCirclesInTheCorners(levelsetFunctionsRef[0]);
    //_generator.generateDiagonalLevelset(levelsetFunctionsRef[1]);
    //ExampleSegmentationGenerator<ConfiguratorType, NumberOfLevelsetFunctions> exampleSegmentationGenerator ( _grid );
    //exampleSegmentationGenerator.generateExample ( levelsetFunctionsRef, 0 );
    this->setSaveName( "phi0" );
    this->saveTimestepBZ2( -1, levelsetFunctionsRef, _grid );

    _generator.generatePiecewiseAffineFlowField ( U0MVec, levelsetFunctionsRef, flowFieldAb ); 

    bool segmentationKnown = true;
    if( _parser.hasVariable ( "VFFileName" ) ) {
      qc::MultiArray<RealType, ConfiguratorType::Dim, ConfiguratorType::Dim> inputVectorField;
      readVectorField ( inputVectorField, _parser.getString ( "VFFileName" ) );
      qc::MultiArray<RealType, ConfiguratorType::Dim, ConfiguratorType::Dim> _v0Array ( _grid, U0MVec );
      _v0Array.resampleFrom( inputVectorField );
      segmentationKnown = false;
    }

    cerr << "Print example ....";
    this->setSaveName( "v0" );
    if ( segmentationKnown )
      this->plotVectorFieldSegmented ( U0MVec, levelsetFunctionsRef, _grid, aol::GNUPLOT_PNG, -1, -1, true );
    else
      this->plotVectorField ( U0MVec, _grid, aol::GNUPLOT_PNG );
    cerr << "done." <<endl;
  }

  void analyzeFidelityError( )
  {
    aol::MultiVector<RealType> u0MVec( FlowFieldDimension, _grid.getNumberOfNodes() );
    initializeTargetVField( u0MVec );
    qc::MultiArray<RealType, ConfiguratorType::Dim, NumberOfLevelsetFunctions> levelsetFunctions ( _grid );

    int dimSize = _parser.getDimSize ( "segmentation", 0 );
    char filename[1024];
    if ( NumberOfLevelsetFunctions != dimSize )
      throw aol::Exception( "NumberOfLevelsetFunctions != dimSize", __FILE__, __LINE__);

    for ( int i = 0; i < dimSize; ++i ){
      _parser.getString ( "segmentation", filename, i );
      levelsetFunctions[i].load( filename );
    }

    aol::MultiVector<RealType> flowFieldAb( 1<<NumberOfLevelsetFunctions, aol::Sqr(FlowFieldDimension)+FlowFieldDimension );

    // Attention: Make sure that the same epsilon is used here that was used when caluculating the segmentation.
    HeavisideFunctionType heavisideFunction( _params.epsilon );

    VectorFieldAffineCoefficientsUpdater<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ConfiguratorType::Dim>
      coeffUpdater ( _grid, heavisideFunction, u0MVec, flowFieldAb, AFFINE );
    coeffUpdater.update ( levelsetFunctions );

    aol::MultiVector<RealType> estimatedFlowField ( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
    qc::DataGenerator<ConfiguratorType> generator( _grid );
    generator.generatePiecewiseAffineFlowField ( estimatedFlowField, levelsetFunctions, coeffUpdater.getFlowFieldAb(), isUsingThresholdModel() ? 0.5 : 0. );

    ArrayType error ( _grid );
    for ( int i = 0; i < error.size(); ++i ) {
      aol::Vec<ConfiguratorType::Dim, RealType> tempVec;
      for ( int j = 0; j < ConfiguratorType::Dim; ++j )
        tempVec[j] = ( estimatedFlowField[j][i] - u0MVec[j][i] );
      error[i] = tempVec.norm();
    }

    this->setSaveName ( "error" );
    this->saveTimestepBZ2( -1, error, _grid );
    this->savePNG( error, _grid );

    cerr << "L_infty error " << error.getMaxValue() << endl;
    cerr << "L_2 error " << error.norm() / error.size() << endl;
  }

  void computeTimesteps( )
  {
    cerr << "NumberOfLevelsetFunctions: " << NumberOfLevelsetFunctions << endl;
    aol::MultiVector<RealType> u0MVec( FlowFieldDimension, _grid.getNumberOfNodes() );

    initializeTargetVField( u0MVec );

    // Startvektoren A,b = 0
    aol::MultiVector<RealType> flowFieldAb( 1<<NumberOfLevelsetFunctions, aol::Sqr(FlowFieldDimension)+FlowFieldDimension );
    flowFieldAb.setZero();

    aol::MultiVector<RealType> levelsetFunctions( NumberOfLevelsetFunctions, _grid.getNumberOfNodes() );
    _generator.generateChanVeseInitialization( levelsetFunctions );
    //_generator.generateCircleLevelset( levelsetFunctions[0], 0.2, 0.25,0.25 );
    //_generator.generateLineLevelset( levelsetFunctions[0], 1, 0.5 );
    //_generator.generateDiagonalLevelsetPlusCirclesInTheCorners(levelsetFunctions[0]);
    //_generator.generateDiagonalLevelset(levelsetFunctions[0]);
    //_generator.generateCircleLevelset ( levelsetFunctions[0], 0.4, 0.3, 0.3 );
    //_generator.generateCircleLevelset ( levelsetFunctions[1], 0.4, 0.7, 0.7 );

    // This thresholding is only appropriate when converting the zero levelline, so
    // that it can be used in a Esedoglu like model.
    if ( _usingThresholdModel )
      levelsetFunctions.threshold ( 0, 0, 1 );

    segment ( u0MVec, levelsetFunctions, flowFieldAb );
  }

  void segment ( const aol::MultiVector<RealType> &U0MVec,
                 aol::MultiVector<RealType> &LevelsetFunctions,
                 aol::MultiVector<RealType> &FlowFieldAb ) {
    for ( int outerStep = 0; outerStep < _params.outerSteps; outerStep++ ) {
      HeavisideFunctionType heavisideFunction( _params.epsilon / aol::Pow( 2, outerStep ) );

      // Energy mit Ab
      FlowField_Ab_Energy<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension>
        tempE( _grid, heavisideFunction, U0MVec );

      aol::HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions>
        Ereg( _grid, heavisideFunction );

      //Variation der Energie
      aol::VariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, NumberOfLevelsetFunctions>
        DEreg( _grid, _params.epsilonMCM, _params.gamma );

      aol::ChanVeseEnergyLevelsetPart<RealType> Efidelity( tempE, FlowFieldAb );

      FlowField_Ab_EnergyMultiPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension>
        DEfidelity( _grid, heavisideFunction, U0MVec, FlowFieldAb );

      aol::LinCombOp<aol::MultiVector<RealType>, aol::Scalar<RealType> > E;
      E.appendReference( Efidelity );
      E.appendReference( Ereg, _params.gamma );

      aol::LinCombOp<aol::MultiVector<RealType> > DE;
      DE.appendReference( DEfidelity );
      DE.appendReference( DEreg );

      ArrayType backgroundImage ( _grid );
      backgroundImage.setAll ( 1. );
      // VectorFieldGradientDescent solves for the levelset functions and
      // for the mean gray values (see class descriptions).
      typedef VectorFieldGradientDescent<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension> GradientDescentType;
      GradientDescentType gradientDescent_solver( _grid, E, DE, heavisideFunction, U0MVec, FlowFieldAb, backgroundImage, _params.numSteps, 
                                                  aol::ZOTrait<RealType>::one, aol::ZOTrait<RealType>::zero, AFFINE );
      gradientDescent_solver.setConfigurationFlags ( GradientDescentType::USE_NONLINEAR_CG );
      gradientDescent_solver.setFlags( GradientDescentType::DO_NOT_SAVE_TIMESTEPS_AS_ISOLINE_PNG );
      if ( _usingThresholdModel )
        gradientDescent_solver.configureForQuadraticEsedogluModel();
      gradientDescent_solver.activateSaving();
      gradientDescent_solver.setNumberSaveFirstPics( 2 );
      gradientDescent_solver.setSaveTimestepOffset( _parser.getInt( "saveOffset") );
      gradientDescent_solver.setSaveName( "phi", (_params.outerSteps > 1 ) ? outerStep : -1 );
      gradientDescent_solver.setSaveDirectory( this->_saveDirectory );
      //gradientDescent_solver.setTauMax( 1.0 );
      //gradientDescent_solver.setStopFilterWidth( 0.01 );
      aol::MultiVector<RealType> mtmp( LevelsetFunctions, aol::STRUCT_COPY );

      gradientDescent_solver.apply( LevelsetFunctions, mtmp );
      LevelsetFunctions = mtmp;

      //_generator.generateTexturePNG(LevelsetFunctions);

      this->setSaveName ( "u" );
      this->setWriteAllTimeSteps ( true );
      this->saveTimestepBZ2( outerStep, LevelsetFunctions, _grid );
      if ( _usingThresholdModel ) {
        for ( int j = 0; j < NumberOfLevelsetFunctions; ++j )
          this->savePNG( LevelsetFunctions[j], this->_grid, j, outerStep );
      }
      this->setWriteAllTimeSteps ( false );

      for (int i = 0; i < FlowFieldAb.numComponents(); i++){
        for (int j = 0; j < FlowFieldAb[i].size(); j++){
          cerr << "FlowFieldAb[" << i << "][" << j << "] = " << FlowFieldAb[i][j] << endl;
        }
        cerr << endl;
      }
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int FlowFieldDimension = ConfiguratorType::Dim>
class VectorFieldSegmentationMultilevelDescent : public qc::MultilevelDescentInterface<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  typedef qc::MultilevelDescentInterface<ConfiguratorType> SuperClass;
  aol::ParameterParser &_parser;
  qc::MultiDimMultilevelArray<RealType> _v0, _levelsetFunctions;
  aol::MultiVector<RealType> _flowFieldAb;

public:
  VectorFieldSegmentationMultilevelDescent ( aol::ParameterParser &Parser )
   : SuperClass ( Parser.getInt( "depth" ) ),
     _parser ( Parser ),
     _v0 ( this->_grid, ConfiguratorType::Dim ),
     _levelsetFunctions ( this->_grid, NumberOfLevelsetFunctions ),
     _flowFieldAb ( 1<<NumberOfLevelsetFunctions, aol::Sqr(FlowFieldDimension)+FlowFieldDimension ) {

    aol::MultiVector<RealType> v0MVec ( 0, 0 );
    _v0.appendReferencesTo ( v0MVec );

    VectorFieldSegmentation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension> vectorFieldSegmentation( _parser );
    vectorFieldSegmentation.initializeTargetVField ( v0MVec );
    _v0.levRestrict( 0, this->getLevel() );

    aol::MultiVector<RealType> _levelsetFunctionsMVec ( 0, 0 );
    _levelsetFunctions.appendReferencesTo ( _levelsetFunctionsMVec );

    qc::DataGenerator<ConfiguratorType>( this->_grid ).generateChanVeseInitialization( _levelsetFunctionsMVec );
    if ( vectorFieldSegmentation.isUsingThresholdModel() )
      _levelsetFunctionsMVec.threshold ( 0, 0, 1 );
    _levelsetFunctions.levRestrict( 0, this->getLevel() );
  }

  virtual void setLevel ( const int Level ) {
    SuperClass::setLevel ( Level );
    _v0.setCurLevel( Level );
    _levelsetFunctions.setCurLevel( Level );
  }

  virtual void prolongate( ) {
    if ( this->getLevel() < this->_grid.getGridDepth() ) {
      _levelsetFunctions.prolongate( );
      setLevel( this->getLevel() + 1 );
    }
    else
      throw aol::Exception( "Can't prolongate if getLevel() >= _grid.getGridDepth()!", __FILE__, __LINE__ );
  }
  
  virtual void descentOnCurrentGrid( ) {
    VectorFieldSegmentation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, FlowFieldDimension> vectorFieldSegmentation( _parser, this->getLevel() );

    aol::MultiVector<RealType> v0MVec ( 0, 0 );
    _v0.appendReferencesTo ( v0MVec );
    aol::MultiVector<RealType> _levelsetFunctionsMVec ( 0, 0 );
    _levelsetFunctions.appendReferencesTo ( _levelsetFunctionsMVec );

    vectorFieldSegmentation.segment ( v0MVec, _levelsetFunctionsMVec, _flowFieldAb );
  }

};

#endif
