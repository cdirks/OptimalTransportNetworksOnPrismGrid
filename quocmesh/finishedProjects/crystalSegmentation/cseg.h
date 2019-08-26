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

#ifndef __CSEG_H
#define __CSEG_H

#include <signedDistanceOp.h>
#include <ChanVese.h>
#include <ArmijoSearch.h>
#include <gradientflow.h>
#include <hyperelastic.h>
#include <deformations.h>
#include <quocTimestepSaver.h>
#include <chanVeseDescent.h>
#include <finiteDifferences.h>


template <typename RealType>
class GrainParams{
public:
  const RealType atomDistance;
  //! gamma weights the perimiter length term
  const RealType gamma;
  //! delta weights the phase fitting term
  const RealType delta;
  //! lambda weights the deformation regularization term
  const RealType lambda;
  const RealType thresholdLower;
  const RealType thresholdUpper;
  const RealType thresholdGradient;
  //! epsilon is the regularization parameter of the Heaviside function
  const RealType epsilon;
  //! epsilon is the regularization parameter of the MCMStiffOp used in the variation of the length term.
  const RealType epsilonMCM;
  //! numStepsis the maximal number of allowed steps in the gradient descent.
  const int numSteps;
  string dirName;
  GrainParams( aol::ParameterParser &Parser )
  : atomDistance( Parser.getDouble( "atomDistance" ) ),
    gamma( Parser.getDouble( "gamma" ) ),
    delta( Parser.getDouble( "delta" ) ),
    lambda( Parser.getDouble( "lambda" ) ),
    thresholdLower( Parser.getDouble( "thresholdLower" ) ),
    thresholdUpper( Parser.getDouble( "thresholdUpper" ) ),
    thresholdGradient( Parser.getDouble( "thresholdGradient" ) ),
    epsilon( Parser.getDouble( "epsilon" ) ),
    epsilonMCM( Parser.getDouble( "epsilonMCM" ) ),
    numSteps( Parser.getInt( "numSteps" ) )
  {
    char fn[1024];
    sprintf( fn, "%f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_%f_%d", gamma, epsilon, lambda, delta, thresholdLower, thresholdUpper, thresholdGradient, atomDistance, numSteps);
    dirName = fn;
  }
  const char* getDirName() const{
    return dirName.c_str();
  }
};

/*
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class RotationEnergyMultiWithDeformationAlphaSegmentationGradientDescent
: public qc::ChanVeseGradientDescent<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions>{
private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::ChanVeseGradientDescent<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> SuperClass;
public:
  RotationEnergyMultiWithDeformationAlphaSegmentationGradientDescent( const typename ConfiguratorType::InitType &Initializer,
                                                                      const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                                                      const aol::Op<aol::MultiVector<RealType> > &DE,
                                                                      const HeavisideFunctionType &HeavisideFunction,
                                                                      const aol::Vector<RealType> &Image,
                                                                      const int MaxIterations = 50,
                                                                      const RealType StartTau = aol::ZOTrait<RealType>::one,
                                                                      const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
    : SuperClass( Initializer, E, DE, HeavisideFunction, Image, MaxIterations, StartTau, StopEpsilon )
  {
    this->_customTauBlocks.resize(NumberOfLevelsetFunctions+2);
    this->_customTauBlocks.setAll( 1 );
    this->_customTauBlocks[NumberOfLevelsetFunctions] = 1<<NumberOfLevelsetFunctions;
    this->_customTauBlocks[NumberOfLevelsetFunctions] = ConfiguratorType::Dim;
  }

  virtual void writeTimeStep ( const aol::MultiVector<RealType> &Data, const int Iteration ) const 
  {
    SuperClass::writeTimeStep( Data, Iteration );
    aol::MultiVector<RealType> deformation( 0, 0 );
    for( int i = 0; i < ConfiguratorType::Dim; i++ )
      deformation.appendReference( Data[i+NumberOfLevelsetFunctions+(1<<NumberOfLevelsetFunctions)] );

    aol::Vector<RealType> tmp( this->_image, aol::STRUCT_COPY );
    qc::DeformImage<ConfiguratorType>( this->_image, this->_grid, tmp, deformation);
    this->savePNG( tmp, this->_grid, Iteration, NumberOfLevelsetFunctions+1 );
    for( int i = 0; i < ConfiguratorType::Dim; i++ )
      this->saveTimestepBZ2 ( NumberOfLevelsetFunctions+2+i, Iteration, Data[i+NumberOfLevelsetFunctions+(1<<NumberOfLevelsetFunctions)], this->_grid );
  }
  virtual void smoothDirection( const aol::MultiVector<RealType> &Position, const RealType Sigma, aol::MultiVector<RealType> &Direction ) const {
    SuperClass::smoothDirection( Position, Sigma, Direction );

    // Smooth the psi-variation
    qc::LinearSmoothOp<RealType> linSmooth;
    linSmooth.setCurrentGrid( this->_grid );
    linSmooth.setSigma( Sigma );

    const int indexOffset = NumberOfLevelsetFunctions+(1<<NumberOfLevelsetFunctions);
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
      linSmooth.apply( Direction[i+indexOffset], Direction[i+indexOffset] );
    }
  }

  virtual bool postProcess( aol::MultiVector<RealType> &Position, const int Iteration ) const {
    SuperClass::postProcess( Position, Iteration );
    qc::SymmetricProjector<ConfiguratorType> symmetricProjector( this->_grid );
    aol::MultiVector<RealType> deformation( 0, 0 );
    for( int i = 0; i < ConfiguratorType::Dim; i++ )
      deformation.appendReference( Position[i+NumberOfLevelsetFunctions+(1<<NumberOfLevelsetFunctions)] );
    symmetricProjector.project( deformation );
    return true;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class SymmetricTransformationGradientDescent
: public aol::GradientDescentWithAutomaticFilterWidth<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> >,
  public qc::QuocTimestepSaver<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;
  const HeavisideFunctionType &_heavisideFunction;
  const aol::Vector<RealType> &_image;
public:
  SymmetricTransformationGradientDescent( const typename ConfiguratorType::InitType &Initializer,
                                          const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                          const aol::Op<aol::MultiVector<RealType> > &DE,
                                          const HeavisideFunctionType &HeavisideFunction,
                                          const aol::Vector<RealType> &Image,
                                          const int MaxIterations = 50,
                                          const RealType StartTau = aol::ZOTrait<RealType>::one,
                                          const RealType StopEpsilon = aol::ZOTrait<RealType>::zero)
  : aol::GradientDescentWithAutomaticFilterWidth<ConfiguratorType, aol::MultiVector<RealType> >(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon),
    _heavisideFunction( HeavisideFunction ),
    _image( Image )
  {
    this->setMaximumNumIterationsPerFilterWidth( 50 );
  }

  virtual bool postProcess( aol::MultiVector<RealType> &Position, const int Iteration ) const {
    if( (Position.numComponents() == ConfiguratorType::Dim + 1) || (Position.numComponents() == ConfiguratorType::Dim + 2) ){
      qc::SymmetricProjector<ConfiguratorType> symmetricProjector( this->_grid );
      symmetricProjector.project( Position );
      if( ( Position.numComponents() == ConfiguratorType::Dim + 2 ) && ( Iteration % 5 == 0 ) ){
        qc::SignedDistanceOp<ConfiguratorType> signedDistOp( this->_grid );
        signedDistOp.apply( Position[ConfiguratorType::Dim+1], Position[ConfiguratorType::Dim+1] );
      }
      // If the zero levelset does not intersect with [0,1]^2, the levelset function is reinitialized.
      if ( Position.numComponents() == ConfiguratorType::Dim + 2 ){
        const RealType minVal = Position[ConfiguratorType::Dim+1].getMinValue();
        const RealType maxVal = Position[ConfiguratorType::Dim+1].getMaxValue();
        if( minVal * maxVal > 0. ){
          qc::DataGenerator<ConfiguratorType> generator( this->_grid );
          generator.generateLineLevelset( Position[ConfiguratorType::Dim+1], 0, 0.5 );
        }
      }
      return true;
    }
    else
      return false;
  }
  virtual void writeTimeStep ( const aol::MultiVector<RealType> &Data, const int Iteration ) const 
  {
    aol::Vector<RealType> tmp( _image, aol::STRUCT_COPY );
    qc::DeformImage<ConfiguratorType>( _image, this->_grid, tmp, Data);
    if( Data.numComponents() == ConfiguratorType::Dim + 2 )
      this->savePNG( tmp, this->_grid, Iteration, 0 );
    else
      this->savePNG( tmp, this->_grid, Iteration );
    if( (Data.numComponents() == ConfiguratorType::Dim + 1) || (Data.numComponents() == ConfiguratorType::Dim + 2) ){
      char fn[1024];
      sprintf( fn, "%salpha.txt", this->_saveDirectory);
      ofstream outfile(fn, ofstream::out | ofstream::app);
      outfile << this->_iterationsFormat(Iteration) << " " << Data[ConfiguratorType::Dim][0];
      if( Data.numComponents() == ConfiguratorType::Dim + 2 ){
        this->saveTimestepPNGDrawRedAndBlueIsoline( 1, Iteration, 0., 0.1, Data[ConfiguratorType::Dim+1], tmp, this->_grid );
        outfile << " " << Data[ConfiguratorType::Dim][1];
        //this->saveTimestepPNGDrawRedAndBlueIsoline( 1, Iteration, 0., 0.1, Data[ConfiguratorType::Dim+1], _image, this->_grid );
      }
      outfile << endl;
      outfile.close();
    }
    if( Data.numComponents() == 2*ConfiguratorType::Dim + 1){
      qc::DeformImageFromMVecPart<ConfiguratorType>( _image, this->_grid, tmp, Data, true, ConfiguratorType::Dim );
      this->savePNG( tmp, this->_grid, Iteration, 1 );
      this->saveTimestepPNGDrawRedAndBlueIsoline( 2, Iteration, 0., 0.1, Data[2*ConfiguratorType::Dim], _image, this->_grid );
    }
  }
  virtual void smoothDirection( const aol::MultiVector<RealType> &Position, const RealType Sigma, aol::MultiVector<RealType> &Direction ) const {
    aol::Vector<RealType> tmp( this->_grid );

    int levelsetFunctionComponentNumber = -1;
    if( Direction.numComponents() == 2*ConfiguratorType::Dim+1 )
      levelsetFunctionComponentNumber = 2*ConfiguratorType::Dim;
    if( Direction.numComponents() == ConfiguratorType::Dim+2 )
      levelsetFunctionComponentNumber = ConfiguratorType::Dim+1;

    qc::LinearSmoothOp<RealType> linSmooth;
    linSmooth.setCurrentGrid( this->_grid );
    linSmooth.setSigma( Sigma );

    if( levelsetFunctionComponentNumber != -1 ){
      aol::HeavisideFunctionLumpedMassOp<ConfiguratorType, HeavisideFunctionType >
        heavisideFunctionLumpedMassInvOpPhaseLevelset( this->_grid, Position[levelsetFunctionComponentNumber], _heavisideFunction, aol::INVERT );
      heavisideFunctionLumpedMassInvOpPhaseLevelset.apply( Direction[levelsetFunctionComponentNumber], tmp );
      Direction[levelsetFunctionComponentNumber] = tmp;
      linSmooth.apply( Direction[levelsetFunctionComponentNumber], Direction[levelsetFunctionComponentNumber] );
    }

    int numberOfComponentsToSmooth = ConfiguratorType::Dim;
    if( Position.numComponents() == 2*ConfiguratorType::Dim + 1 )
      numberOfComponentsToSmooth = 2*ConfiguratorType::Dim;
    for( int i = 0; i < numberOfComponentsToSmooth; i++ )
      linSmooth.apply( Direction[i], Direction[i] );
  }


  virtual void getTauAndUpdateDescentDir( aol::MultiVector<RealType> &DescentDir,
                                          const aol::MultiVector<RealType> &CurrentPosition,
                                          aol::Vector<RealType> &Tau ) const {
    aol::Vector<int> tauBlocks;
    if( DescentDir.numComponents() == ConfiguratorType::Dim ){
      tauBlocks.resize( 1 );
      tauBlocks[0] = ConfiguratorType::Dim;
    }
    if( DescentDir.numComponents() == ConfiguratorType::Dim + 1 ){
      tauBlocks.resize( 2 );
      tauBlocks[0] = ConfiguratorType::Dim;
      tauBlocks[1] = 1;
    }
    if( DescentDir.numComponents() == ConfiguratorType::Dim + 2 ){
      tauBlocks.resize( 3 );
      tauBlocks[0] = ConfiguratorType::Dim;
      tauBlocks[1] = 1;
      tauBlocks[2] = 1;
    }
    if( DescentDir.numComponents() == 2*ConfiguratorType::Dim + 1 ){
      tauBlocks.resize( 2 );
      tauBlocks[0] = 2*ConfiguratorType::Dim;
      tauBlocks[1] = 1;
    }
    if( tauBlocks.sum() != DescentDir.numComponents() )
      throw aol::Exception( "tauBlocks.sum() != DescentDir.numComponents()", __FILE__, __LINE__ );

    const int numComponents = tauBlocks.size();
    if( numComponents != Tau.size() ){
      Tau.resize( numComponents );
      Tau.setAll( 1. );
    }
    aol::MultiVector<RealType> descentDirComponent( DescentDir, aol::STRUCT_COPY );

    int handledDDComponents = 0;
    for( int component = 0; component < numComponents; component++){
      for( int i = 0; i < tauBlocks[component]; i++)
        descentDirComponent[i+handledDDComponents] = DescentDir[i+handledDDComponents];
      Tau[component] = this->getTimestepWidthWithArmijoLineSearch(descentDirComponent, CurrentPosition, Tau[component]);
      for( int i = 0; i < tauBlocks[component]; i++){
        DescentDir[i+handledDDComponents] *= Tau[component];
        descentDirComponent[i+handledDDComponents].setZero();
      }
      handledDDComponents += tauBlocks[component];
    }
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, typename VectorType>
class LevelsetGradientDescent{
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class LevelsetGradientDescent<ConfiguratorType, HeavisideFunctionType, aol::Vector<typename ConfiguratorType::RealType> >
: public aol::GradientDescent<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> >,
  public qc::QuocTimestepSaver<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;
  const HeavisideFunctionType &_heavisideFunction;
  const aol::Vector<RealType> &_image;
public:
  LevelsetGradientDescent( const typename ConfiguratorType::InitType &Initializer,
                           const aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > &E,
                           const aol::Op<aol::Vector<RealType> > &DE,
                           const HeavisideFunctionType &HeavisideFunction,
                           const aol::Vector<RealType> &Image,
                           const int MaxIterations = 50,
                           const RealType StartTau = aol::ZOTrait<RealType>::one,
                           const RealType StopEpsilon = aol::ZOTrait<RealType>::zero)
  : aol::GradientDescent<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> >(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon),
   _heavisideFunction( HeavisideFunction ),
   _image(Image)
  {
    this->setUpdateFilterWidth( true );
  }
  virtual bool postProcess( aol::Vector<RealType> &Position, const int /*Iteration*/ ) const {
    qc::SignedDistanceOp<ConfiguratorType> signedDistOp( this->_grid );
    signedDistOp.apply( Position, Position );
    return true;
  }
  virtual void smoothDirection( const aol::Vector<RealType> &Position, const RealType Sigma, aol::Vector<RealType> &Direction ) const {
    aol::Vector<RealType> tmp( this->_grid );

    aol::HeavisideFunctionLumpedMassOp<ConfiguratorType, HeavisideFunctionType >
      heavisideFunctionLumpedMassInvOp( this->_grid, Position, _heavisideFunction, aol::INVERT );
    heavisideFunctionLumpedMassInvOp.apply( Direction, tmp );

    Direction = tmp;

    qc::LinearSmoothOp<RealType> linSmooth;
    linSmooth.setCurrentGrid( this->_grid );
    linSmooth.setSigma( Sigma );
//    linSmooth.setTau( 5 * this->_grid.H() );
    linSmooth.apply( Direction, Direction );

  }
  virtual void writeTimeStep( const aol::Vector<RealType> &Dest, const int Iterations ) const {
    this->saveTimestepPNGDrawRedAndBlueIsoline( 1, Iterations, 0., 0.1, Dest, _image, this->_grid );
    this->saveTimestepPNGDrawRedIsoline( 2, Iterations, 0., Dest, _image, this->_grid );
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class LevelsetGradientDescent<ConfiguratorType, HeavisideFunctionType, aol::MultiVector<typename ConfiguratorType::RealType> >
//: public aol::GradientDescent<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> >,
: public aol::GradientDescentComponentWiseTimestepControlled<ConfiguratorType>,
  public qc::QuocTimestepSaver<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;
  const HeavisideFunctionType &_heavisideFunction;
  const aol::Vector<RealType> &_image;
public:
  LevelsetGradientDescent( const typename ConfiguratorType::InitType &Initializer,
                           const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                           const aol::Op<aol::MultiVector<RealType> > &DE,
                           const HeavisideFunctionType &HeavisideFunction,
                           const aol::Vector<RealType> &Image,
                           const int MaxIterations = 50,
                           const RealType StartTau = aol::ZOTrait<RealType>::one,
                           const RealType StopEpsilon = aol::ZOTrait<RealType>::zero)
//  : aol::GradientDescent<ConfiguratorType, aol::MultiVector<RealType> >(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon),
  : aol::GradientDescentComponentWiseTimestepControlled<ConfiguratorType>(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon),
   _heavisideFunction( HeavisideFunction ),
   _image(Image)
  {
    this->setUpdateFilterWidth( true );
    this->setMaximumNumIterationsPerFilterWidth( 20 );
  }

  virtual bool postProcess( aol::MultiVector<RealType> &Position, const int Iteration ) const {
    qc::SignedDistanceOp<ConfiguratorType> signedDistOp( this->_grid );
    if( Iteration % 5 == 0 ){
      signedDistOp.apply( Position[0], Position[0] );
      if( Position.numComponents() == 3 ){
        signedDistOp.apply( Position[2], Position[2] );
      }
      return true;
    }
    else
      return false;
  }

  virtual void smoothDirection( const aol::MultiVector<RealType> &Position, const RealType Sigma, aol::MultiVector<RealType> &Direction ) const {
    aol::Vector<RealType> tmp( this->_grid );

    qc::LinearSmoothOp<RealType> linSmooth;
    linSmooth.setCurrentGrid( this->_grid );
    linSmooth.setSigma( Sigma );
//    linSmooth.setTau( 5 * this->_grid.H() );

    aol::HeavisideFunctionLumpedMassOp<ConfiguratorType, HeavisideFunctionType >
      heavisideFunctionLumpedMassInvOpAlphaLevelset( this->_grid, Position[0], _heavisideFunction, aol::INVERT );
    heavisideFunctionLumpedMassInvOpAlphaLevelset.apply( Direction[0], tmp );

    Direction[0] = tmp;

    linSmooth.apply( Direction[0], Direction[0] );

    if( Direction.numComponents() == 3 ){
      aol::HeavisideFunctionLumpedMassOp<ConfiguratorType, HeavisideFunctionType >
        heavisideFunctionLumpedMassInvOpPhaseLevelset( this->_grid, Position[2], _heavisideFunction, aol::INVERT );
      heavisideFunctionLumpedMassInvOpPhaseLevelset.apply( Direction[2], tmp );

      Direction[2] = tmp;

      linSmooth.apply( Direction[2], Direction[2] );
      Direction[2] *= 0.1;
    }

    Direction[1] *= 0.5;

  }
  virtual void writeTimeStep( const aol::MultiVector<RealType> &MDest, const int Iterations ) const {
    this->saveTimestepPNGDrawRedAndBlueIsoline( 1, Iterations, 0., 0.1, MDest[0], _image, this->_grid );
    this->saveTimestepPNGDrawRedIsoline( 3, Iterations, 0., MDest[0], _image, this->_grid );
    this->saveTimestepBZ2( 5, Iterations, MDest[0], this->_grid );
    if( MDest.numComponents() == 3 ){
      this->saveTimestepPNGDrawRedAndBlueIsoline( 2, Iterations, 0., 0.1, MDest[2], _image, this->_grid );
      this->saveTimestepPNGDrawRedIsoline( 4, Iterations, 0., MDest[2], _image, this->_grid );
      this->saveTimestepBZ2( 6, Iterations, MDest[2], this->_grid );
    }
    char fn[1024];
    sprintf( fn, "%salpha.txt", this->_saveDirectory);
    ofstream outfile(fn, ofstream::out | ofstream::app);
    outfile << this->_iterationsFormat(Iterations) << " " << MDest[1][0] << " " << MDest[1][1] << endl;
    outfile.close();
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
inline typename ConfiguratorType::RealType AlphaIndicatorSum
  ( const qc::GridDefinition &Grid,
    const qc::Element &El,
    int QuadPoint, 
    const typename ConfiguratorType::VecType& RefCoord,
    const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
    const typename ConfiguratorType::RealType Threshold,
    const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
    const HeavisideFunctionType &HeavisideFunction )
{
  typedef typename ConfiguratorType::RealType RealType;
  RealType temp = 0.;
  for( int j = 0; j < 6; j++ ){
    typename ConfiguratorType::VecType transformed_local_coord;
    qc::Element transformed_el;
    if( qc::transformCoord<ConfiguratorType> ( Grid, El, RefCoord, Q[j], transformed_el, transformed_local_coord ) )
      temp += (1 - HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold) );
  }
  return (HeavisideFunction.evaluate( DiscrU0.evaluateAtQuadPoint( El, QuadPoint ) - Threshold)*temp);
}

template <typename ConfiguratorType, typename HeavisideFunctionType>
inline typename ConfiguratorType::RealType AlphaPsiIndicatorSum
  ( const qc::GridDefinition &Grid,
    const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
    const qc::Element &El,
    int QuadPoint, 
    const typename ConfiguratorType::VecType& RefCoord,
    const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
    const typename ConfiguratorType::RealType Threshold,
    const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
    const HeavisideFunctionType &HeavisideFunction )
{
  typedef typename ConfiguratorType::RealType RealType;
  RealType temp = 0.;
  typename ConfiguratorType::VecType transformed_local_coord;
  qc::Element transformed_el;
  for( int j = 0; j < 6; j++ ){
    if( qc::translateAndTransformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, RefCoord, Q[j], transformed_el, transformed_local_coord) ){
      temp += (1 - HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold) );
    }
  }
  if( qc::transformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord ) )
    return (HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold)*temp);
  else
    return 0.;
}

//! The first part is not affected by the integral transformation, the second part is.
//! This has to be taken into account, if the integrand is multiplied by another function,
//! for example the H(\phi(x))
template <typename ConfiguratorType, typename HeavisideFunctionType>
inline void AlphaPsiIndicatorSumDerivativeWRTPsi
  ( const qc::GridDefinition &Grid,
    const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
    const qc::Element &El,
    int QuadPoint,
    const typename ConfiguratorType::VecType& RefCoord,
    const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
    const typename ConfiguratorType::RealType Threshold,
    const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
    const HeavisideFunctionType &HeavisideFunction,
    aol::Vec<ConfiguratorType::Dim,typename ConfiguratorType::RealType> &NL )
{
  //first the part, which is not affected by the integral transformation
  typedef typename ConfiguratorType::RealType RealType;
  RealType part1 = 0.;
  RealType part2 = 0.;
  typename ConfiguratorType::VecType transformed_local_coord;
  qc::Element transformed_el;
  for( int j = 0; j < 6; j++ ){
    if( qc::translateAndTransformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, RefCoord, Q[j], transformed_el, transformed_local_coord) ){
      part1 += (1 - HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold) );
    }
  }
  //second the part, which is affected by the integral transformation
  for( int j = 0; j < 6; j++ ){
    if( qc::translateAndTransformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, RefCoord, Q[j+12], transformed_el, transformed_local_coord) ){
      part2 -= HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold);
    }
  }
  if( qc::transformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord ) ){
    DiscrU0.evaluateGradient( transformed_el, transformed_local_coord, NL );
    NL *= (HeavisideFunction.evaluateDerivative( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold)*(part1+part2));
  }
  else{
    NL.setZero();
    return;
  }
}

template <typename ConfiguratorType, typename HeavisideFunctionType>
inline typename ConfiguratorType::RealType AlphaPsiIndicatorSumDerivativeWRTAlpha
  ( const qc::GridDefinition &Grid,
    const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
    const qc::Element &El,
    int QuadPoint, 
    const typename ConfiguratorType::VecType& RefCoord,
    const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
    const typename ConfiguratorType::RealType Threshold,
    const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
    const HeavisideFunctionType &HeavisideFunction )
{
  typename ConfiguratorType::MatType dPsi;
  typename ConfiguratorType::VecType tempVec;
  typedef typename ConfiguratorType::RealType RealType;
  RealType temp = 0.;
  typename ConfiguratorType::VecType transformed_local_coord;
  qc::Element transformed_el;
  for( int j = 0; j < 6; j++ ){
    if( qc::transformCoord<ConfiguratorType> ( Grid, El, RefCoord, Q[j], transformed_el, transformed_local_coord ) ){
      for ( int c = 0; c < ConfiguratorType::Dim; c++ ) {
        DiscrTransformation[c].evaluateGradient ( transformed_el, transformed_local_coord, tempVec );
        dPsi.setRow ( c, tempVec );
        dPsi[c][c] += 1.;
      }
      if( qc::translateAndTransformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, RefCoord, Q[j], transformed_el, transformed_local_coord) ){
        typename ConfiguratorType::VecType gradU;
        DiscrU0.evaluateGradient( transformed_el, transformed_local_coord, gradU );
        dPsi.mult( Q[j+6], tempVec );
        temp += HeavisideFunction.evaluateDerivative( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold)
                *(gradU*tempVec);
      }
    }
  }
  if( qc::transformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord ) )
    return (-1.*HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold)*temp);
  else
    return 0.;
}

//! The first part is not affected by the integral transformation, the second part is.
//! This has to be taken into account, if the integrand is multiplied by another function,
//! for example the H(\phi(x))
template <typename ConfiguratorType, typename HeavisideFunctionType>
inline void AlphaPsiPhiIndicatorSumDerivativeWRTPsi
  ( const qc::GridDefinition &Grid,
    const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
    const qc::Element &El,
    int QuadPoint,
    const typename ConfiguratorType::VecType& RefCoord,
    const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
    const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrLevelsetFunction,
    const typename ConfiguratorType::RealType Threshold,
    const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
    const HeavisideFunctionType &HeavisideFunction,
    aol::Vec<ConfiguratorType::Dim,typename ConfiguratorType::RealType> &NL,
    const bool OneMinusLevelsetHeaviside )
{
  //first the part, which is not affected by the integral transformation
  typedef typename ConfiguratorType::RealType RealType;
  RealType part1 = 0.;
  RealType part2 = 0.;
  RealType temp = 0.;
  typename ConfiguratorType::VecType transformed_local_coord;
  qc::Element transformed_el;
  for( int j = 0; j < 6; j++ ){
    if( qc::translateAndTransformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, RefCoord, Q[j], transformed_el, transformed_local_coord) ){
      part1 += (1 - HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold) );
    }
  }
  const RealType h = HeavisideFunction.evaluate( DiscrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint ) );
  if( OneMinusLevelsetHeaviside )
    part1 *= (1.-h);
  else
    part1 *= h;
  //second the part, which is affected by the integral transformation
  for( int j = 0; j < 6; j++ ){
    if( qc::transformCoord<ConfiguratorType> ( Grid, El, RefCoord, Q[j+12], transformed_el, transformed_local_coord ) ){
      const RealType h = HeavisideFunction.evaluate( DiscrLevelsetFunction.evaluate( transformed_el, transformed_local_coord ) );
      if( OneMinusLevelsetHeaviside )
        temp = (1.-h);
      else
        temp = h;
      if( qc::translateAndTransformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, RefCoord, Q[j+12], transformed_el, transformed_local_coord) ){
        part2 -= temp * HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold);
      }
    }
  }
  if( qc::transformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord ) ){
    DiscrU0.evaluateGradient( transformed_el, transformed_local_coord, NL );
    NL *= (HeavisideFunction.evaluateDerivative( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold)*(part1+part2));
  }
  else{
    NL.setZero();
    return;
  }
}

//! The first part is not affected by the integral transformation, the second part is.
//! This has to be taken into account, if the integrand is multiplied by another function,
//! for example the H(\phi(x))
template <typename ConfiguratorType, typename HeavisideFunctionType>
inline void AlphaPsiMultiPhiIndicatorSumDerivativeWRTPsi
  ( const qc::GridDefinition &Grid,
    const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
    const qc::Element &El,
    int QuadPoint,
    const typename ConfiguratorType::VecType& RefCoord,
    const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
    const aol::MultiVector<typename ConfiguratorType::RealType> &LevelsetFunctions,
    const typename ConfiguratorType::RealType Threshold,
    const std::vector<std::vector<aol::Vec2<typename ConfiguratorType::RealType> >*> QVec,
    const HeavisideFunctionType &HeavisideFunction,
    aol::Vec<ConfiguratorType::Dim,typename ConfiguratorType::RealType> &NL )
{
  //first the part, which is not affected by the integral transformation
  typedef typename ConfiguratorType::RealType RealType;
  RealType part1 = 0.;
  RealType part2 = 0.;
  RealType temp = 0.;
  typename ConfiguratorType::VecType transformed_local_coord;
  qc::Element transformed_el;

  // Construct HeavisideFunctionProduct at the point, where we want to evaluate it.
  aol::HeavisideFunctionProduct<ConfiguratorType, HeavisideFunctionType>
    heavisideFunctionProduct( Grid, HeavisideFunction, LevelsetFunctions, El, QuadPoint );
  {
    // Loop over all combinations of Levelset function signs, each represents one of the indicator parameters.
    aol::PowerSetIterator iterator( LevelsetFunctions.numComponents() );
    do
    {
      temp = 0;
      const int i = iterator.getCurrentPosition();
      for( int j = 0; j < 6; j++ ){
        if( qc::translateAndTransformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, RefCoord, (*(QVec[i]))[j], transformed_el, transformed_local_coord) ){
          temp += (1 - HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold) );
        }
      }
      part1 += temp*heavisideFunctionProduct.evaluate( iterator );
      iterator.increment();
    } while ( !iterator.end() );
  }
  //second the part, which is affected by the integral transformation
  {
    // Loop over all combinations of Levelset function signs, each represents one of the indicator parameters.
    aol::PowerSetIterator iterator( LevelsetFunctions.numComponents() );
    do
    {
      const int i = iterator.getCurrentPosition();
      for( int j = 0; j < 6; j++ ){
        if( qc::transformCoord<ConfiguratorType> ( Grid, El, RefCoord, (*(QVec[i]))[j+12], transformed_el, transformed_local_coord ) ){
          // Construct HeavisideFunctionProduct at the point, where we want to evaluate it.
          aol::HeavisideFunctionProduct<ConfiguratorType, HeavisideFunctionType>
            heavisideFunctionProductTranslated( Grid, HeavisideFunction, LevelsetFunctions, transformed_el, transformed_local_coord );
          temp = heavisideFunctionProductTranslated.evaluate( iterator );
          if( qc::translateAndTransformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, RefCoord, (*(QVec[i]))[j+12], transformed_el, transformed_local_coord) ){
            part2 -= temp * HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold);
          }
        }
      }
      iterator.increment();
    } while ( !iterator.end() );
  }
  // Moving this check up to the beginning could speed up things a little.
  if( qc::transformCoord<ConfiguratorType> ( Grid, DiscrTransformation, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord ) ){
    DiscrU0.evaluateGradient( transformed_el, transformed_local_coord, NL );
    NL *= (HeavisideFunction.evaluateDerivative( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold)*(part1+part2));
  }
  else{
    NL.setZero();
    return;
  }
}

//! This function is not properly tested and not needed at the moment!
template <typename ConfiguratorType, typename HeavisideFunctionType>
inline typename ConfiguratorType::RealType AlphaIndicatorProduct
  ( const qc::GridDefinition &Grid,
    const qc::Element &El,
    int QuadPoint, 
    const typename ConfiguratorType::VecType& RefCoord,
    const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
    const typename ConfiguratorType::RealType Threshold,
    const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
    const HeavisideFunctionType &HeavisideFunction )
{
  typedef typename ConfiguratorType::RealType RealType;
  RealType temp = 1.;
  for( int j = 0; j < 6; j++ ){
    typename ConfiguratorType::VecType transformed_local_coord;
    qc::Element transformed_el;
    if( qc::transformCoord<ConfiguratorType> ( Grid, El, RefCoord, Q[j], transformed_el, transformed_local_coord ) )
      temp *= (1 - HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold) );
  }
  return (HeavisideFunction.evaluate( DiscrU0.evaluateAtQuadPoint( El, QuadPoint ) - Threshold)*temp);
}

template <typename ConfiguratorType, typename HeavisideFunctionType>
inline typename ConfiguratorType::RealType AlphaIndicatorSumDerivative
  ( const qc::GridDefinition &Grid,
    const qc::Element &El,
    int QuadPoint, 
    const typename ConfiguratorType::VecType& RefCoord,
    const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
    const typename ConfiguratorType::RealType Threshold,
    const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
    const HeavisideFunctionType &HeavisideFunction )
{
  typedef typename ConfiguratorType::RealType RealType;
  RealType temp = 0.;
  for( int j = 0; j < 6; j++ ){
    typename ConfiguratorType::VecType transformed_local_coord;
    qc::Element transformed_el;
    if( qc::transformCoord<ConfiguratorType> ( Grid, El, RefCoord, Q[j], transformed_el, transformed_local_coord ) ){
      typename ConfiguratorType::VecType gradU;
      DiscrU0.evaluateGradient( transformed_el, transformed_local_coord, gradU );
      temp += HeavisideFunction.evaluateDerivative( DiscrU0.evaluate( transformed_el, transformed_local_coord ) - Threshold)
              *(gradU*Q[j+6]);
    }
  }
  return (-1.*HeavisideFunction.evaluate( DiscrU0.evaluateAtQuadPoint( El, QuadPoint ) - Threshold)*temp);
}

//! This function is not properly tested and not needed at the moment!
template <typename ConfiguratorType, typename HeavisideFunctionType>
inline typename ConfiguratorType::RealType AlphaIndicatorProductDerivative
  ( const qc::GridDefinition &Grid,
    const qc::Element &El,
    int QuadPoint, 
    const typename ConfiguratorType::VecType& RefCoord,
    const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
    const typename ConfiguratorType::RealType Threshold,
    const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
    const HeavisideFunctionType &HeavisideFunction )
{
  typedef typename ConfiguratorType::RealType RealType;
  RealType temp = 0.;
  std::vector<typename ConfiguratorType::VecType> transformed_local_coordVec(6);
  std::vector<qc::Element> transformed_elVec(6);
  std::vector<bool> transformedPositionStillInDomainVec(6);
  for( int j = 0; j < 6; j++ ){
    transformedPositionStillInDomainVec[j] = qc::transformCoord<ConfiguratorType> ( Grid, El, RefCoord, Q[j], transformed_elVec[j], transformed_local_coordVec[j] );
  }
  for( int j = 0; j < 6; j++ ){
    if( transformedPositionStillInDomainVec[j] ){
      RealType productValue = 1.;
      for( int k = 0; k < 6; k++ ){
        if( transformedPositionStillInDomainVec[k] && k != j ){
          productValue *= (1 - HeavisideFunction.evaluate( DiscrU0.evaluate( transformed_elVec[k], transformed_local_coordVec[k] ) - Threshold) );
        }
      }
      typename ConfiguratorType::VecType gradU;
      DiscrU0.evaluateGradient( transformed_elVec[j], transformed_local_coordVec[j], gradU );
      temp += HeavisideFunction.evaluateDerivative( DiscrU0.evaluate( transformed_elVec[j], transformed_local_coordVec[j] ) - Threshold)
              *(gradU*Q[j+6])*productValue;
    }
  }
  return (-1.*HeavisideFunction.evaluate( DiscrU0.evaluateAtQuadPoint( El, QuadPoint ) - Threshold)*temp);
}

template <typename ConfiguratorType, typename HeavisideFunctionType>
class PhaseEnergy
: public aol::FENonlinIntegrationVectorInterface<ConfiguratorType,
                                                   PhaseEnergy<ConfiguratorType, HeavisideFunctionType>,
                                                   1 > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU0;
  const RealType _delta;
  const RealType _thresholdLower;
  const RealType _thresholdUpper;
  const RealType _thresholdGradient;
public:
  PhaseEnergy( const typename ConfiguratorType::InitType &Initializer,
               const HeavisideFunctionType &HeavisideFunction,
               const aol::Vector<RealType> &U0,
               const RealType Delta,
               const RealType ThresholdLower,
               const RealType ThresholdUpper,
               const RealType ThresholdGradient )
  : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                               PhaseEnergy<ConfiguratorType, HeavisideFunctionType>,
                                               1 > ( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _discrU0( Initializer, U0 ),
      _delta( Delta ),
      _thresholdLower( ThresholdLower ),
      _thresholdUpper( ThresholdUpper ),
      _thresholdGradient( ThresholdGradient ) {
  }

  RealType evaluateIntegrand( const aol::auto_container<1,aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    const RealType phi = discrFuncs[0].evaluateAtQuadPoint( El, QuadPoint);
    const RealType u0 = _discrU0.evaluateAtQuadPoint( El, QuadPoint);
    const RealType HPhi = _heavisideFunction.evaluate( phi );
    typename ConfiguratorType::VecType gradU;
    _discrU0.evaluateGradientAtQuadPoint( El, QuadPoint, gradU );
    const RealType HU0MinusThresholdLower = _heavisideFunction.evaluate( u0 - _thresholdLower );
    const RealType HThresholdUpperMinusU0 = _heavisideFunction.evaluate( _thresholdUpper - u0 );
    const RealType HEpsilonMinusNormGradU = _heavisideFunction.evaluate( _thresholdGradient - gradU.norm() );
    const RealType temp = 1. - HU0MinusThresholdLower*HThresholdUpperMinusU0*HEpsilonMinusNormGradU;
    return ( HPhi*temp + _delta*(1.-HPhi)*(1.-temp) );
  }
};

//! Computes \f$ \int ( 1 - 2H(u-\theta) )\phi_i \f$, where arg=u.
template <typename ConfiguratorType, typename HeavisideFunctionType>
class PhaseGradient : public aol::FENonlinOpInterface< ConfiguratorType, PhaseGradient<ConfiguratorType, HeavisideFunctionType> > {
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _delta;
  const RealType _thresholdLower;
  const RealType _thresholdUpper;
  const RealType _thresholdGradient;
public:
  PhaseGradient( const typename ConfiguratorType::InitType &Initializer,
                 const HeavisideFunctionType &HeavisideFunction,
                 const RealType Delta,
                 const RealType ThresholdLower,
                 const RealType ThresholdUpper,
                 const RealType ThresholdGradient )
    : aol::FENonlinOpInterface< ConfiguratorType,
                                 PhaseGradient<ConfiguratorType, HeavisideFunctionType> > ( Initializer ),
    _heavisideFunction( HeavisideFunction ),
    _delta( Delta ),
    _thresholdLower( ThresholdLower ),
    _thresholdUpper( ThresholdUpper ),
    _thresholdGradient( ThresholdGradient ) {
  }

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc, 
			const typename ConfiguratorType::ElementType &El, 
			int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
			typename ConfiguratorType::RealType &NL ) const {
    const RealType u0 = DiscFunc.evaluateAtQuadPoint( El, QuadPoint );
    typename ConfiguratorType::VecType gradU;
    DiscFunc.evaluateGradientAtQuadPoint( El, QuadPoint, gradU );
    const RealType HU0MinusThresholdLower = _heavisideFunction.evaluate( u0 - _thresholdLower );
    const RealType HThresholdUpperMinusU0 = _heavisideFunction.evaluate( _thresholdUpper - u0 );
    const RealType HEpsilonMinusNormGradU = _heavisideFunction.evaluate( _thresholdGradient - gradU.norm() );
    const RealType temp = 1. - HU0MinusThresholdLower*HThresholdUpperMinusU0*HEpsilonMinusNormGradU;

    NL = ((1+_delta)*temp-_delta);
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class SolidLiquidPhaseMSSegmentor : public qc::TwoPhaseMSSegmentor<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::ArrayType &_image;
  mutable typename ConfiguratorType::ArrayType _indicatorBase;
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _delta;
  const RealType _thresholdLower;
  const RealType _thresholdUpper;
  const RealType _thresholdGradient;
  const bool _normalizeIndicatorParts;
public:
  SolidLiquidPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                                const typename ConfiguratorType::ArrayType &Image,
                                const HeavisideFunctionType &HeavisideFunction,
                                const RealType Gamma,
                                const RealType Delta,
                                const RealType ThresholdLower,
                                const RealType ThresholdUpper,
                                const RealType ThresholdGradient,
                                const bool NormalizeIndicatorParts )
    : qc::TwoPhaseMSSegmentor<ConfiguratorType> ( Initializer, Gamma ),
      _image ( Image ),
      _indicatorBase ( Initializer ),
      _heavisideFunction ( HeavisideFunction ),
      _delta( Delta ),
      _thresholdLower( ThresholdLower ),
      _thresholdUpper( ThresholdUpper ),
      _thresholdGradient( ThresholdGradient ),
      _normalizeIndicatorParts ( NormalizeIndicatorParts ) {}

protected:
  void prepareIndicatorFunctionGeneration ( ) const {
    typename ConfiguratorType::ArrayType imageGradientNorm ( this->_grid );
    qc::DataGenerator<ConfiguratorType> generator ( this->_grid );
    generator.generateGradientNormFD ( _image, imageGradientNorm );

    const RealType lowerNormalization = _normalizeIndicatorParts ? ( _thresholdLower - _image.getMinValue() ) : 1;
    const RealType upperNormalization = _normalizeIndicatorParts ? ( _image.getMaxValue() - _thresholdUpper ) : 1;
    const RealType gradNormalization = _normalizeIndicatorParts ? ( imageGradientNorm.getMaxValue() - _thresholdGradient ) : 1;

    for ( int i = 0; i < _indicatorBase.size(); ++i) {
      const RealType u0 = _image[i];
      const RealType HU0MinusThresholdLower = _heavisideFunction.evaluate( ( u0 - _thresholdLower ) / lowerNormalization );
      const RealType HThresholdUpperMinusU0 = _heavisideFunction.evaluate( ( _thresholdUpper - u0 ) / upperNormalization );
      const RealType HEpsilonMinusNormGradU = _heavisideFunction.evaluate( (_thresholdGradient - imageGradientNorm[i]) / gradNormalization );
      _indicatorBase[i] = 1. - HU0MinusThresholdLower*HThresholdUpperMinusU0*HEpsilonMinusNormGradU;
    }
  }
  virtual void generateIndicatorFunction ( const int IndicatorNumber, qc::ScalarArray<RealType, qc::QC_2D> &IndicatorFunction ) const {
    const RealType shift = 0.5;
    for ( int i = 0; i < IndicatorFunction.size(); ++i )
      IndicatorFunction[i] = ( ( IndicatorNumber == 0 ) ? _indicatorBase[i] : ( 1 - _indicatorBase[i] ) ) + shift;
  }
};

template <typename RealType>
void initQVectors( const RealType Alpha, const RealType AtomDistance, std::vector<aol::Vec2<RealType> > &Q ){
  const RealType cosAlpha = cos( Alpha );
  const RealType sinAlpha = sin( Alpha );
  //const RealType scaling = 4*M_PI/(sqrt(3)*100);
  //const RealType scaling = 4*M_PI/(sqrt(3)*200);
  //const RealType scaling = sqrt( aol::Sqr( (71.-86.)/129. ) + aol::Sqr( (97.-98.)/129. ) );
  //const RealType scaling = 5./129.;
  //const RealType scaling = 9./129.;
  //const RealType scaling = 7./257.;
  //const RealType scaling = 11./513.;
  //const RealType scaling = 11./129.;
  //const RealType scaling = 18.5/257.;
  //const RealType scaling = 12./257.;
  //const RealType scaling = 9./257.;
  //const RealType scaling = 6./129.;
  for( int i = 0; i < 6; i ++ ){
    const RealType angle = 2.*aol::NumberTrait<RealType>::pi/6.*i;
    //! \f$ M(\alpha)q_i \f$
    Q[i][0] = AtomDistance*(cosAlpha*cos( angle ) - sinAlpha*sin( angle ));
    Q[i][1] = AtomDistance*(sinAlpha*cos( angle ) + cosAlpha*sin( angle ));
    //! \f$ M^\prime(\alpha)q_i \f$
    Q[i+6][0] = AtomDistance*( -sinAlpha*cos( angle ) - cosAlpha*sin( angle ) );
    Q[i+6][1] = AtomDistance*( cosAlpha*cos( angle ) - sinAlpha*sin( angle ) );
    //! \f$ M(\alpha)q_i \f$
    Q[i+12][0] = -Q[i][0];
    Q[i+12][1] = -Q[i][1];
  }
}

/**
 * Helper function to calculate the data term if the grain segmentation is formulated as continuous labeling
 * problem and lifted to get a convex formulation.
 *
 * \author Berkels
 */
template <typename RealType>
void grainSegConvexDataTerm ( qc::ScalarArray<RealType, qc::QC_3D> &G, const qc::ScalarArray<RealType, qc::QC_2D> &Image, const RealType Epsilon, const RealType AtomDistance, const RealType Threshold ) {
  std::vector<aol::Vec2<RealType> > qVec ( 18 );
  aol::ArcTanHeavisideFunction<RealType> heavisideFunc ( Epsilon );

  const int m = G.getNumX();
  const int n = G.getNumY();
  const int l = G.getNumZ();
  for ( int k = 0; k < l; ++k ) {
    const RealType alpha = k*aol::NumberTrait<RealType>::pi/l;
    initQVectors( alpha, AtomDistance, qVec );
    for ( int j = 0; j < n; ++j ) {
      for ( int i = 0; i < m; ++i ) {
        RealType tmp = 0;
        for( int c = 0; c < 6; c ++ ) {
          const RealType x = i+qVec[c][0];
          const RealType y = j+qVec[c][1];

          if ( ( x >= 0 ) && ( x < m ) && ( y >= 0 ) && ( y < n ) )
            tmp += ( 1 - heavisideFunc.evaluate ( Image.interpolate ( x, y ) - Threshold ) );
        }
        G.set ( i, j, k, ( heavisideFunc.evaluate ( Image.get ( i, j )  - Threshold ) * tmp ) );
      }
    }
  }
}

/**
 * Calculates the rotated base vectors of the atom grid for all
 * angles alpha contained in the first entry of all components
 * of Parameters.
 *
 * \author Berkels
 */
template <typename RealType>
void initQVectorsForAllParameters( const aol::MultiVector<RealType> &Parameters, const RealType AtomDistance, const int NumberOfLevelsetFunctions, const std::vector<std::vector<aol::Vec2<RealType> >*> &QVec ){
#ifdef BOUNDS_CHECK
  if( Parameters.numComponents() != static_cast<int>(QVec.size()) )
    throw aol::Exception( "Parameters.size() != QVec.size()", __FILE__, __LINE__ );
#endif
  // The PowerSetIterator worries about how the indicator parameters are numbered.
  aol::PowerSetIterator iterator( NumberOfLevelsetFunctions );
  do
  {
    const int i = iterator.getCurrentPosition();
    initQVectors<RealType>( Parameters[i][0], AtomDistance, *(QVec[i]) );
    iterator.increment();
  } while ( !iterator.end() );
}

//! Computes \f$ \int H(u-\theta) \sum_j(1-H(u(x+M(\alpha)q_j)-\theta) \f$, where arg[0]=u.
template <typename ConfiguratorType, typename HeavisideFunctionType>
class RotationEnergy
: public aol::FENonlinIntegrationVectorInterface<ConfiguratorType,
                                                   RotationEnergy<ConfiguratorType, HeavisideFunctionType>,
                                                   1 > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _threshold;
  std::vector<aol::Vec2<RealType> > _q;
  RealType (*_pIndicatorType)( const qc::GridDefinition &Grid,
                               const qc::Element &El,
                               int QuadPoint, 
                               const typename ConfiguratorType::VecType& RefCoord,
                               const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
                               const typename ConfiguratorType::RealType Threshold,
                               const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
                               const HeavisideFunctionType &HeavisideFunction );
public:
  RotationEnergy( const typename ConfiguratorType::InitType &Initializer,
                  const HeavisideFunctionType &HeavisideFunction,
                  const RealType Alpha,
                  const RealType AtomDistance,
                  const RealType Threshold,
                  RealType (*PIndicatorType)( const qc::GridDefinition &Grid,
                                               const qc::Element &El,
                                               int QuadPoint, 
                                               const typename ConfiguratorType::VecType& RefCoord,
                                               const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
                                               const typename ConfiguratorType::RealType Threshold,
                                               const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
                                               const HeavisideFunctionType &HeavisideFunction ) )
  : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                               RotationEnergy<ConfiguratorType, HeavisideFunctionType>,
                                               1 > ( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _alpha( Alpha ),
      _atomDistance( AtomDistance ),
      _threshold( Threshold ),
      _q(18)
  {
    initQVectors<RealType>( _alpha, _atomDistance, _q );
    _pIndicatorType = PIndicatorType;
  }

  RealType evaluateIntegrand( const aol::auto_container<1,aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &RefCoord ) const {
    return (*_pIndicatorType)( this->_initializer, El, QuadPoint, RefCoord, discrFuncs[0], _threshold, _q, _heavisideFunction );
  }
};

//! Computes \f$ \int H(u-\theta) \sum_j(1-H(u(x+M(\alpha)q_j)-\theta)\phi_i \f$, where arg=u.
template <typename ConfiguratorType, typename HeavisideFunctionType>
class RotationGradient : public aol::FENonlinOpInterface< ConfiguratorType, RotationGradient<ConfiguratorType, HeavisideFunctionType> > {
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _threshold;
  std::vector<aol::Vec2<RealType> > _q;
  RealType (*_pIndicatorType)( const qc::GridDefinition &Grid,
                               const qc::Element &El,
                               int QuadPoint, 
                               const typename ConfiguratorType::VecType& RefCoord,
                               const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
                               const typename ConfiguratorType::RealType Threshold,
                               const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
                               const HeavisideFunctionType &HeavisideFunction );
public:
  RotationGradient( const typename ConfiguratorType::InitType &Initializer,
                    const HeavisideFunctionType &HeavisideFunction,
                    const RealType Alpha,
                    const RealType AtomDistance,
                    const RealType Threshold,
                    RealType (*PIndicatorType)( const qc::GridDefinition &Grid,
                                                const qc::Element &El,
                                                int QuadPoint, 
                                                const typename ConfiguratorType::VecType& RefCoord,
                                                const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
                                                const typename ConfiguratorType::RealType Threshold,
                                                const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
                                                const HeavisideFunctionType &HeavisideFunction ) )
    : aol::FENonlinOpInterface< ConfiguratorType,
                                 RotationGradient<ConfiguratorType, HeavisideFunctionType> > ( Initializer ),
    _heavisideFunction( HeavisideFunction ),
    _alpha( Alpha ),
    _atomDistance( AtomDistance ),
    _threshold( Threshold ),
    _q(18),
    _pIndicatorType( PIndicatorType )
  {
    initQVectors<RealType>( _alpha, _atomDistance, _q );
  }

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc, 
			const typename ConfiguratorType::ElementType &El, 
			int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
			typename ConfiguratorType::RealType &NL ) const {
    NL = (*_pIndicatorType)( this->_initializer, El, QuadPoint, RefCoord, DiscFunc, _threshold, _q, _heavisideFunction );
  }
};

//! Computes \f$ \int H(u(psi)-\theta) \sum_j(1-H(u(psi(x+M(\alpha)q_j))-\theta) \f$, where MArg=psi.
template <typename ConfiguratorType, typename HeavisideFunctionType>
class RotationEnergyWithDeformation
: public aol::FENonlinIntegrationVectorInterface<ConfiguratorType,
                                                   RotationEnergyWithDeformation<ConfiguratorType, HeavisideFunctionType>,
                                                   ConfiguratorType::Dim > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _threshold;
  std::vector<aol::Vec2<RealType> > _q;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU0;
  RealType (*_pIndicatorType)( const qc::GridDefinition &Grid,
                               const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
                               const qc::Element &El,
                               int QuadPoint, 
                               const typename ConfiguratorType::VecType& RefCoord,
                               const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
                               const typename ConfiguratorType::RealType Threshold,
                               const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
                               const HeavisideFunctionType &HeavisideFunction );
public:
  RotationEnergyWithDeformation( const typename ConfiguratorType::InitType &Initializer,
                  const HeavisideFunctionType &HeavisideFunction,
                  const RealType Alpha,
                  const RealType AtomDistance,
                  const RealType Threshold,
                  const aol::Vector<RealType> &U0Dofs,
                  RealType (*PIndicatorType)( const qc::GridDefinition &Grid,
                                              const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
                                              const qc::Element &El,
                                              int QuadPoint, 
                                              const typename ConfiguratorType::VecType& RefCoord,
                                              const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
                                              const typename ConfiguratorType::RealType Threshold,
                                              const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
                                              const HeavisideFunctionType &HeavisideFunction ) )
  : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                               RotationEnergyWithDeformation<ConfiguratorType, HeavisideFunctionType>,
                                               ConfiguratorType::Dim > ( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _alpha( Alpha ),
      _atomDistance( AtomDistance ),
      _threshold( Threshold ),
      _q(18),
      _discrU0( Initializer, U0Dofs )
  {
    initQVectors<RealType>( _alpha, _atomDistance, _q );
    _pIndicatorType = PIndicatorType;
  }

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &RefCoord ) const {
    return (*_pIndicatorType)( this->_initializer, discrFuncs, El, QuadPoint, RefCoord, _discrU0, _threshold, _q, _heavisideFunction );
  }
};


template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfRotationEnergyWithDeformationWRTPsi
: public aol::FENonlinVectorOpInterface< ConfiguratorType,
                                          ConfiguratorType::Dim,
                                          ConfiguratorType::Dim,
                                          VariationOfRotationEnergyWithDeformationWRTPsi<ConfiguratorType, HeavisideFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _threshold;
  std::vector<aol::Vec2<RealType> > _q;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU0;
public:
  VariationOfRotationEnergyWithDeformationWRTPsi( const typename ConfiguratorType::InitType &Initializer,
                                                  const HeavisideFunctionType &HeavisideFunction,
                                                  const RealType Alpha,
                                                  const RealType AtomDistance,
                                                  const RealType Threshold,
                                                  const aol::Vector<RealType> &U0Dofs )
  : aol::FENonlinVectorOpInterface< ConfiguratorType,
                                     ConfiguratorType::Dim,
                                     ConfiguratorType::Dim,
                                     VariationOfRotationEnergyWithDeformationWRTPsi<ConfiguratorType, HeavisideFunctionType> > ( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _alpha( Alpha ),
      _atomDistance( AtomDistance ),
      _threshold( Threshold ),
      _q(18),
      _discrU0( Initializer, U0Dofs )
  {
    initQVectors<RealType>( _alpha, _atomDistance, _q );
  }

  void getNonlinearity( aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                        aol::Vec<ConfiguratorType::Dim,typename ConfiguratorType::RealType> &NL ) const {
    AlphaPsiIndicatorSumDerivativeWRTPsi<ConfiguratorType, HeavisideFunctionType>( this->_initializer, DiscFuncs, El, QuadPoint, RefCoord, _discrU0, _threshold, _q, _heavisideFunction, NL );
  }
};


template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfRotationEnergyWithDeformationAndLevelsetWRTPsi
: public aol::FENonlinVectorOpInterface< ConfiguratorType,
                                          ConfiguratorType::Dim,
                                          ConfiguratorType::Dim,
                                          VariationOfRotationEnergyWithDeformationAndLevelsetWRTPsi<ConfiguratorType, HeavisideFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _threshold;
  std::vector<aol::Vec2<RealType> > _q;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU0;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const bool _oneMinusLevelsetHeaviside;
public:
  VariationOfRotationEnergyWithDeformationAndLevelsetWRTPsi( const typename ConfiguratorType::InitType &Initializer,
                                                             const HeavisideFunctionType &HeavisideFunction,
                                                             const RealType Alpha,
                                                             const RealType AtomDistance,
                                                             const RealType Threshold,
                                                             const aol::Vector<RealType> &U0Dofs,
                                                             const aol::Vector<RealType> &LevelsetFunctionDofs,
                                                             const bool OneMinusLevelsetHeaviside )
  : aol::FENonlinVectorOpInterface< ConfiguratorType,
                                     ConfiguratorType::Dim,
                                     ConfiguratorType::Dim,
                                     VariationOfRotationEnergyWithDeformationAndLevelsetWRTPsi<ConfiguratorType, HeavisideFunctionType> > ( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _alpha( Alpha ),
      _atomDistance( AtomDistance ),
      _threshold( Threshold ),
      _q(18),
      _discrU0( Initializer, U0Dofs ),
      _discrLevelsetFunction( Initializer, LevelsetFunctionDofs ),
      _oneMinusLevelsetHeaviside( OneMinusLevelsetHeaviside )
  {
    initQVectors<RealType>( _alpha, _atomDistance, _q );
  }

  void getNonlinearity( aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                        aol::Vec<ConfiguratorType::Dim,typename ConfiguratorType::RealType> &NL ) const {
    AlphaPsiPhiIndicatorSumDerivativeWRTPsi<ConfiguratorType, HeavisideFunctionType>( this->_initializer, DiscFuncs, El, QuadPoint, RefCoord, _discrU0, _discrLevelsetFunction, _threshold, _q, _heavisideFunction, NL, _oneMinusLevelsetHeaviside );
  }
};


//! arg = u0
template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfRotationEnergyWithDeformationAndLevelsetWRTPhi
: public aol::FENonlinOpInterface< ConfiguratorType,
                                    VariationOfRotationEnergyWithDeformationAndLevelsetWRTPhi<ConfiguratorType, HeavisideFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _threshold;
  std::vector<aol::Vec2<RealType> > _q;
  aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > _discrTrans;
public:
  VariationOfRotationEnergyWithDeformationAndLevelsetWRTPhi
  ( const typename ConfiguratorType::InitType &Initializer,
    const HeavisideFunctionType &HeavisideFunction,
    const RealType Alpha,
    const RealType AtomDistance,
    const RealType Threshold,
    const aol::MultiVector<RealType> &TransformationDofs )
  : aol::FENonlinOpInterface< ConfiguratorType,
                               VariationOfRotationEnergyWithDeformationAndLevelsetWRTPhi<ConfiguratorType, HeavisideFunctionType>
                               > ( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _alpha( Alpha ),
      _atomDistance( AtomDistance ),
      _threshold( Threshold ),
      _q(18)
  {
    for ( int c = 0; c < TransformationDofs.numComponents(); c++ ) {
      aol::DiscreteFunctionDefault<ConfiguratorType> temp( Initializer, TransformationDofs[c] );
      _discrTrans.set_copy( c, temp );
    }
    initQVectors<RealType>( _alpha, _atomDistance, _q );
  }

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
			const typename ConfiguratorType::ElementType &El,
			int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
			typename ConfiguratorType::RealType &NL ) const {
    NL = AlphaPsiIndicatorSum<ConfiguratorType, HeavisideFunctionType>( this->_initializer, _discrTrans, El, QuadPoint, RefCoord, DiscFunc, _threshold, _q, _heavisideFunction );
  }
};


//! Computes \f$ -\int H(u(\psi)-\theta) \sum_{i=1,\cdots,m} H^\prime(u(\psi(x+M(\alpha)q_i))-\theta)\nabla u(\psi(x+M(\alpha)q_i))\cdot D\psi(x) M^\prime(\alpha)q_i \f$, where marg[0]=\psi.
template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfRotationEnergyWithDeformationWRTAlpha
: public aol::FENonlinIntegrationVectorInterface<ConfiguratorType,
                                                   VariationOfRotationEnergyWithDeformationWRTAlpha<ConfiguratorType, HeavisideFunctionType>,
                                                   ConfiguratorType::Dim > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _threshold;
  std::vector<aol::Vec2<RealType> > _q;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU0;
public:
  VariationOfRotationEnergyWithDeformationWRTAlpha( const typename ConfiguratorType::InitType &Initializer,
                                                    const HeavisideFunctionType &HeavisideFunction,
                                                    const RealType Alpha,
                                                    const RealType AtomDistance,
                                                    const RealType Threshold,
                                                    const aol::Vector<RealType> &U0Dofs )
  : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                               VariationOfRotationEnergyWithDeformationWRTAlpha<ConfiguratorType, HeavisideFunctionType>,
                                               ConfiguratorType::Dim > ( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _alpha( Alpha ),
      _atomDistance( AtomDistance ),
      _threshold( Threshold ),
      _q(18),
      _discrU0( Initializer, U0Dofs )
  {
    initQVectors<RealType>( _alpha, _atomDistance, _q );
  }

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &RefCoord ) const {
    return AlphaPsiIndicatorSumDerivativeWRTAlpha<ConfiguratorType, HeavisideFunctionType>( this->_initializer, DiscFuncs, El, QuadPoint, RefCoord, _discrU0, _threshold, _q, _heavisideFunction );
  }
};


//! Computes \f$ -\int H(u-\theta) \sum_j H^\prime(u(x+M(\alpha)q_j)-\theta)\nabla u(x+M(\alpha)q_j)\cdot M^\prime(\alpha)q_i \f$, where arg[0]=u.
template <typename ConfiguratorType, typename HeavisideFunctionType>
class AlphaDerivative
: public aol::FENonlinIntegrationVectorInterface<ConfiguratorType,
                                                   AlphaDerivative<ConfiguratorType, HeavisideFunctionType>,
                                                   1 > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _threshold;
  std::vector<aol::Vec2<RealType> > _q;
  RealType (*_pIndicatorDerivativeType)( const qc::GridDefinition &Grid,
                                         const qc::Element &El,
                                         int QuadPoint, 
                                         const typename ConfiguratorType::VecType& RefCoord,
                                         const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
                                         const typename ConfiguratorType::RealType Threshold,
                                         const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
                                         const HeavisideFunctionType &HeavisideFunction );
public:
  AlphaDerivative( const typename ConfiguratorType::InitType &Initializer,
                   const HeavisideFunctionType &HeavisideFunction,
                   const RealType Alpha,
                   const RealType AtomDistance,
                   const RealType Threshold,
                   RealType (*PIndicatorDerivativeType)( const qc::GridDefinition &Grid,
                                                         const qc::Element &El,
                                                         int QuadPoint, 
                                                         const typename ConfiguratorType::VecType& RefCoord,
                                                         const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrU0,
                                                         const typename ConfiguratorType::RealType Threshold,
                                                         const std::vector<aol::Vec2<typename ConfiguratorType::RealType> > &Q,
                                                         const HeavisideFunctionType &HeavisideFunction ) )
  : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                               AlphaDerivative<ConfiguratorType, HeavisideFunctionType>,
                                               1 > ( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _alpha( Alpha ),
      _atomDistance( AtomDistance ),
      _threshold( Threshold ),
      _q(18),
      _pIndicatorDerivativeType( PIndicatorDerivativeType )
  {
    initQVectors<RealType>( _alpha, _atomDistance, _q );
  }

  RealType evaluateIntegrand( const aol::auto_container<1,aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &RefCoord ) const {
    return (*_pIndicatorDerivativeType)( this->_initializer, El, QuadPoint, RefCoord, discrFuncs[0], _threshold, _q, _heavisideFunction );
  }
};


/**
 * Class that collects everything the RotationEnergyMulti and its
 * variations have in common.
 *
 * \author Berkels
 *
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class RotationEnergyCommons {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType _grid;
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _atomDistance;
  const RealType _threshold;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU0;
  std::vector<std::vector<aol::Vec2<RealType> >*> _qVec;
public:
  RotationEnergyCommons( const typename ConfiguratorType::InitType &Initializer,
                         const HeavisideFunctionType &HeavisideFunction,
                         const RealType AtomDistance,
                         const RealType Threshold,
                         const aol::Vector<RealType> &U0,
                         const int NumberOfLevelsetFunctions )
    : _grid( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _atomDistance( AtomDistance ),
      _threshold( Threshold ),
      _discrU0( Initializer, U0 ),
      _qVec(0)
  {
    for( int i = 0; i < 1<<NumberOfLevelsetFunctions; i++ )
      _qVec.push_back( new std::vector<aol::Vec2<RealType> >(18) );
  }
  ~RotationEnergyCommons () {
    for( unsigned int i = 0; i < _qVec.size(); i++ )
      delete _qVec[i];
  }

  RealType evaluateIndicator ( const aol::Vector<RealType> &/*IndicatorParameter*/,
                               const int IndicatorParameterIndex,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &RefCoord ) const {
    return AlphaIndicatorSum<ConfiguratorType, HeavisideFunctionType>
             ( _grid, El, QuadPoint, RefCoord, _discrU0, _threshold, *(_qVec[IndicatorParameterIndex]), _heavisideFunction );
  }

  void evaluateIndicatorDerivative ( const aol::Vector<RealType> &/*IndicatorParameter*/,
                                     const int IndicatorParameterIndex,
                                     const typename ConfiguratorType::ElementType &El,
                                     int QuadPoint,
                                     const typename ConfiguratorType::VecType &RefCoord,
                                     aol::Vec<1, RealType> &Derivative ) const {
    Derivative[0] = AlphaIndicatorSumDerivative<ConfiguratorType, HeavisideFunctionType>
             ( _grid, El, QuadPoint, RefCoord, _discrU0, _threshold, *(_qVec[IndicatorParameterIndex]), _heavisideFunction );
  }
};

/**
 * Calulates the Rotation Energy. Alternative to calculate it with the classes 
 * RotationEnergy and HeavisideAndOneMinusHFENonlinIntegrationVectorInterface,
 * which is currently used in CompleteRotationEnergy
 *
 * \author Berkels
 *
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class RotationEnergyMulti
: public aol::ChanVeseEnergyInterface< ConfiguratorType,
                                       HeavisideFunctionType,
                                       NumberOfLevelsetFunctions,
                                       1,
                                       RotationEnergyMulti<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >,
  public RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  RotationEnergyMulti( const typename ConfiguratorType::InitType &Initializer,
                       const HeavisideFunctionType &HeavisideFunction,
                       const RealType AtomDistance,
                       const RealType Threshold,
                       const aol::Vector<RealType> &U0 )
    : aol::ChanVeseEnergyInterface< ConfiguratorType,
                                    HeavisideFunctionType,
                                    NumberOfLevelsetFunctions,
                                    1,
                                    RotationEnergyMulti<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >
                                  ( Initializer, HeavisideFunction ),
      RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType>( Initializer, HeavisideFunction, AtomDistance, Threshold, U0, NumberOfLevelsetFunctions )
  {}
  ~RotationEnergyMulti () {}

  virtual void cacheIndicatorParameters( const aol::MultiVector<RealType> &Parameters ) const{
    initQVectorsForAllParameters<RealType>( Parameters, this->_atomDistance, NumberOfLevelsetFunctions, this->_qVec );
  };

  using RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType>::evaluateIndicator;
};

/**
 * Calulates the alpha derivative of the Rotation Energy. Alternative to calculate
 * it with the classes AlphaDerivative and HeavisideFENonlinIntegrationVectorInterface,
 * which is currently used in VariationOfCompleteRotationEnergy
 *
 * \author Berkels
 *
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class RotationEnergyAlphaVariation
: public aol::ChanVeseEnergyParameterDerivativeInterface< ConfiguratorType,
                                                          HeavisideFunctionType,
                                                          NumberOfLevelsetFunctions,
                                                          1,
                                                          RotationEnergyAlphaVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >,
  public RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  RotationEnergyAlphaVariation( const typename ConfiguratorType::InitType &Initializer,
                                const HeavisideFunctionType &HeavisideFunction,
                                const RealType AtomDistance,
                                const RealType Threshold,
                                const aol::Vector<RealType> &U0,
                                const aol::MultiVector<RealType> &LevelsetFunctions )
    : aol::ChanVeseEnergyParameterDerivativeInterface< ConfiguratorType,
                                                       HeavisideFunctionType,
                                                       NumberOfLevelsetFunctions,
                                                       1,
                                                       RotationEnergyAlphaVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >
                                                     ( Initializer, LevelsetFunctions, HeavisideFunction ),
      RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType>( Initializer, HeavisideFunction, AtomDistance, Threshold, U0, LevelsetFunctions.numComponents() )
  {}
  ~RotationEnergyAlphaVariation () {}
  virtual void cacheIndicatorParameters( const aol::MultiVector<RealType> &Arg ) const{
    initQVectorsForAllParameters<RealType>( Arg, this->_atomDistance, this->_levelsetFunctions.numComponents(), this->_qVec );
  };

  using RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType>::evaluateIndicatorDerivative;
};

/**
 * Calulates the phi derivative of the Rotation Energy. Alternative to calculate
 * it with two instances of the class RotationGradient,
 * which is currently used in VariationOfCompleteRotationEnergy
 *
 * \author Berkels
 *
 */

template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class RotationEnergyPhiVariation
: public aol::ChanVeseEnergyLevelsetDerivativeInterface< ConfiguratorType,
                                                         HeavisideFunctionType,
                                                         NumberOfLevelsetFunctions,
                                                         RotationEnergyPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >,
  public RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  RotationEnergyPhiVariation( const typename ConfiguratorType::InitType &Initializer,
                              const HeavisideFunctionType &HeavisideFunction,
                              const RealType AtomDistance,
                              const RealType Threshold,
                              const aol::Vector<RealType> &U0,
                              const aol::MultiVector<RealType> &Parameters )
    : aol::ChanVeseEnergyLevelsetDerivativeInterface< ConfiguratorType,
                                                      HeavisideFunctionType,
                                                      NumberOfLevelsetFunctions,
                                                      RotationEnergyPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >
                                                     ( Initializer, Parameters, HeavisideFunction ),
      RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType>( Initializer, HeavisideFunction, AtomDistance, Threshold, U0, NumberOfLevelsetFunctions )
  {
    initQVectorsForAllParameters<RealType>( Parameters, this->_atomDistance, NumberOfLevelsetFunctions, this->_qVec );
  }
  ~RotationEnergyPhiVariation () {}

  using RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType>::evaluateIndicator;
};

/**
 * Class that collects everything the RotationEnergyMultiWithDeformation and its
 * variations have in common.
 *
 * \author Berkels
 *
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class RotationEnergyWithDeformationCommons 
: public RotationEnergyCommons <ConfiguratorType, HeavisideFunctionType>{
public:
  typedef typename ConfiguratorType::RealType RealType;
  aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > _discrTransformation;
  RotationEnergyWithDeformationCommons( const typename ConfiguratorType::InitType &Initializer,
                                        const HeavisideFunctionType &HeavisideFunction,
                                        const RealType AtomDistance,
                                        const RealType Threshold,
                                        const aol::MultiVector<RealType> &TransformationDofs,
                                        const aol::Vector<RealType> &U0,
                                        const int NumberOfLevelsetFunctions )
    : RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType>( Initializer, HeavisideFunction, AtomDistance, Threshold, U0, NumberOfLevelsetFunctions )
  {
    for ( int c = 0; c < TransformationDofs.numComponents(); c++ ) {
      aol::DiscreteFunctionDefault<ConfiguratorType> temp( Initializer, TransformationDofs[c] );
      _discrTransformation.set_copy( c, temp );
    }
  }

  RealType evaluateIndicator ( const aol::Vector<RealType> &/*IndicatorParameter*/,
                               const int IndicatorParameterIndex,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &RefCoord ) const {
    return AlphaPsiIndicatorSum<ConfiguratorType, HeavisideFunctionType>
             ( this->_grid, _discrTransformation, El, QuadPoint, RefCoord, this->_discrU0, this->_threshold, *(this->_qVec[IndicatorParameterIndex]), this->_heavisideFunction );
  }

  void evaluateIndicatorDerivative ( const aol::Vector<RealType> &/*IndicatorParameter*/,
                                     const int IndicatorParameterIndex,
                                     const typename ConfiguratorType::ElementType &El,
                                     int QuadPoint,
                                     const typename ConfiguratorType::VecType &RefCoord,
                                     aol::Vec<1, RealType> &Derivative ) const {
    Derivative[0] = AlphaPsiIndicatorSumDerivativeWRTAlpha<ConfiguratorType, HeavisideFunctionType>
                      ( this->_grid, _discrTransformation, El, QuadPoint, RefCoord, this->_discrU0, this->_threshold, *(this->_qVec[IndicatorParameterIndex]), this->_heavisideFunction );
  }
};

/*
 * Calulates the Rotation Energy with a deformation. Alternative to calculate it with the classes 
 * RotationEnergyWithDeformation and HeavisideAndOneMinusHFENonlinIntegrationVectorInterface,
 * which is currently used in CompleteElasticEnergyWithAlphaSegmentation
 *
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class RotationEnergyMultiWithDeformation
: public aol::ChanVeseEnergyInterface< ConfiguratorType,
                                       HeavisideFunctionType,
                                       NumberOfLevelsetFunctions,
                                       1,
                                       RotationEnergyMultiWithDeformation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >,
  public RotationEnergyWithDeformationCommons<ConfiguratorType, HeavisideFunctionType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  RotationEnergyMultiWithDeformation( const typename ConfiguratorType::InitType &Initializer,
                                      const HeavisideFunctionType &HeavisideFunction,
                                      const RealType AtomDistance,
                                      const RealType Threshold,
                                      const aol::MultiVector<RealType> &TransformationDofs,
                                      const aol::Vector<RealType> &U0 )
    : aol::ChanVeseEnergyInterface< ConfiguratorType,
                                    HeavisideFunctionType,
                                    NumberOfLevelsetFunctions,
                                    1,
                                    RotationEnergyMultiWithDeformation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >
                                  ( Initializer, HeavisideFunction ),
      RotationEnergyWithDeformationCommons<ConfiguratorType, HeavisideFunctionType>( Initializer, HeavisideFunction, AtomDistance, Threshold, TransformationDofs, U0, NumberOfLevelsetFunctions )
  {}
  ~RotationEnergyMultiWithDeformation () {}

  virtual void cacheIndicatorParameters( const aol::MultiVector<RealType> &Parameters ) const{
    initQVectorsForAllParameters<RealType>( Parameters, this->_atomDistance, NumberOfLevelsetFunctions, this->_qVec );
  };

  using RotationEnergyWithDeformationCommons<ConfiguratorType, HeavisideFunctionType>::evaluateIndicator;
};

/**
 * Calulates the alpha derivative of the Rotation Energy with a deformation. Alternative to calculate
 * it with two instances of VariationOfRotationEnergyWithDeformationWRTAlpha and HeavisideFENonlinIntegrationVectorInterface
 * which is currently used in VariationOfCompleteElasticEnergyWithAlphaSegmentation.
 *
 * \author Berkels
 *
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class RotationEnergyMultiWithDeformationAlphaVariation
: public aol::ChanVeseEnergyParameterDerivativeInterface< ConfiguratorType,
                                                          HeavisideFunctionType,
                                                          NumberOfLevelsetFunctions,
                                                          1,
                                                          RotationEnergyMultiWithDeformationAlphaVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >,
  public RotationEnergyWithDeformationCommons<ConfiguratorType, HeavisideFunctionType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  RotationEnergyMultiWithDeformationAlphaVariation( const typename ConfiguratorType::InitType &Initializer,
                                                    const HeavisideFunctionType &HeavisideFunction,
                                                    const RealType AtomDistance,
                                                    const RealType Threshold,
                                                    const aol::MultiVector<RealType> &TransformationDofs,
                                                    const aol::Vector<RealType> &U0,
                                                    const aol::MultiVector<RealType> &LevelsetFunctions )
    : aol::ChanVeseEnergyParameterDerivativeInterface< ConfiguratorType,
                                                       HeavisideFunctionType,
                                                       NumberOfLevelsetFunctions,
                                                       1,
                                                       RotationEnergyMultiWithDeformationAlphaVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >
                                                     ( Initializer, LevelsetFunctions, HeavisideFunction ),
      RotationEnergyWithDeformationCommons<ConfiguratorType, HeavisideFunctionType>( Initializer, HeavisideFunction, AtomDistance, Threshold, TransformationDofs, U0, LevelsetFunctions.numComponents() )
  {}
  ~RotationEnergyMultiWithDeformationAlphaVariation () {}
  virtual void cacheIndicatorParameters( const aol::MultiVector<RealType> &Arg ) const{
    initQVectorsForAllParameters<RealType>( Arg, this->_atomDistance, this->_levelsetFunctions.numComponents(), this->_qVec );
  };

  using RotationEnergyWithDeformationCommons<ConfiguratorType, HeavisideFunctionType>::evaluateIndicatorDerivative;
};

/**
 * Calulates the phi derivative of the Rotation Energy with a deformation. Alternative to calculate
 * it with two instances of the class VariationOfRotationEnergyWithDeformationAndLevelsetWRTPhi,
 * which is currently used in VariationOfCompleteElasticEnergyWithAlphaSegmentation
 *
 * \author Berkels
 *
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class RotationEnergyMultiWithDeformationPhiVariation
: public aol::ChanVeseEnergyLevelsetDerivativeInterface< ConfiguratorType,
                                                         HeavisideFunctionType,
                                                         NumberOfLevelsetFunctions,
                                                         RotationEnergyMultiWithDeformationPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >,
  public RotationEnergyWithDeformationCommons<ConfiguratorType, HeavisideFunctionType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  RotationEnergyMultiWithDeformationPhiVariation( const typename ConfiguratorType::InitType &Initializer,
                                                  const HeavisideFunctionType &HeavisideFunction,
                                                  const RealType AtomDistance,
                                                  const RealType Threshold,
                                                  const aol::MultiVector<RealType> &TransformationDofs,
                                                  const aol::Vector<RealType> &U0,
                                                  const aol::MultiVector<RealType> &Parameters )
    : aol::ChanVeseEnergyLevelsetDerivativeInterface< ConfiguratorType,
                                                      HeavisideFunctionType,
                                                      NumberOfLevelsetFunctions,
                                                      RotationEnergyMultiWithDeformationPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >
                                                     ( Initializer, Parameters, HeavisideFunction ),
      RotationEnergyWithDeformationCommons<ConfiguratorType, HeavisideFunctionType>( Initializer, HeavisideFunction, AtomDistance, Threshold, TransformationDofs, U0, NumberOfLevelsetFunctions )
  {
    initQVectorsForAllParameters<RealType>( Parameters, this->_atomDistance, NumberOfLevelsetFunctions, this->_qVec );
  }
  ~RotationEnergyMultiWithDeformationPhiVariation () {}

  using RotationEnergyWithDeformationCommons<ConfiguratorType, HeavisideFunctionType>::evaluateIndicator;
};

/**
 * Calulates the psi derivative of the Rotation Energy with a deformation. Alternative to calculate
 * it with two instances of the class VariationOfRotationEnergyWithDeformationAndLevelsetWRTPsi,
 * which is currently used in VariationOfCompleteElasticEnergyWithAlphaSegmentation
 *
 * \author Berkels
 *
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class RotationEnergyMultiWithDeformationPsiVariation
: public aol::FENonlinVectorOpInterface< ConfiguratorType,
                                          ConfiguratorType::Dim,
                                          ConfiguratorType::Dim,
                                          RotationEnergyMultiWithDeformationPsiVariation<ConfiguratorType, HeavisideFunctionType> >,
  public RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::MultiVector<RealType> &_levelsetFunctions;
public:
  RotationEnergyMultiWithDeformationPsiVariation( const typename ConfiguratorType::InitType &Initializer,
                                                  const HeavisideFunctionType &HeavisideFunction,
                                                  const RealType AtomDistance,
                                                  const RealType Threshold,
                                                  const aol::Vector<RealType> &U0,
                                                  const aol::MultiVector<RealType> &LevelsetFunctions,
                                                  const aol::MultiVector<RealType> &Parameters )
  : aol::FENonlinVectorOpInterface< ConfiguratorType,
                                     ConfiguratorType::Dim,
                                     ConfiguratorType::Dim,
                                     RotationEnergyMultiWithDeformationPsiVariation<ConfiguratorType, HeavisideFunctionType> > ( Initializer ),
      RotationEnergyCommons<ConfiguratorType, HeavisideFunctionType>( Initializer, HeavisideFunction, AtomDistance, Threshold, U0, LevelsetFunctions.numComponents() ),
      _levelsetFunctions( LevelsetFunctions )
  {
    initQVectorsForAllParameters<RealType>( Parameters, this->_atomDistance, LevelsetFunctions.numComponents(), this->_qVec );
  }

  void getNonlinearity( aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                        aol::Vec<ConfiguratorType::Dim,typename ConfiguratorType::RealType> &NL ) const {
    AlphaPsiMultiPhiIndicatorSumDerivativeWRTPsi<ConfiguratorType, HeavisideFunctionType>
      ( this->_grid, DiscFuncs, El, QuadPoint, RefCoord, this->_discrU0, _levelsetFunctions, this->_threshold, this->_qVec, this->_heavisideFunction, NL );
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class CompletePhaseEnergy : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _gamma;
  const RealType _delta;
  const RealType _thresholdLower;
  const RealType _thresholdUpper;
  const RealType _thresholdGradient;
  const HeavisideFunctionType &_heavisideFunction;
public:
  CompletePhaseEnergy( const typename ConfiguratorType::InitType &Initializer,
                       const aol::Vector<RealType> &U0,
                       const RealType Gamma,
                       const RealType Delta,
                       const RealType ThresholdLower,
                       const RealType ThresholdUpper,
                       const RealType ThresholdGradient,
                       const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _gamma( Gamma ),
  _delta( Delta ),
  _thresholdLower( ThresholdLower ),
  _thresholdUpper( ThresholdUpper ),
  _thresholdGradient( ThresholdGradient ),
  _heavisideFunction( HeavisideFunction )
  {
  }
  virtual void apply( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const{
    //! Perimiter length term
    aol::HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType >
      heavisideLevelsetLengthEnergy( _grid, _heavisideFunction );
    aol::MultiVector<RealType> mPhi( 0, 0 );
    mPhi.appendReference( Arg );
    heavisideLevelsetLengthEnergy.apply( mPhi, Dest);
    Dest *= _gamma;

    PhaseEnergy<ConfiguratorType, HeavisideFunctionType> phaseEnergy( _grid, _heavisideFunction, _u0, _delta, _thresholdLower, _thresholdUpper, _thresholdGradient );
    phaseEnergy.applyAdd( mPhi, Dest );
  }
  virtual void applyAdd( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( Arg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfCompletePhaseEnergy : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _gamma;
  const RealType _delta;
  const RealType _thresholdLower;
  const RealType _thresholdUpper;
  const RealType _thresholdGradient;
  const HeavisideFunctionType &_heavisideFunction;
public:
  VariationOfCompletePhaseEnergy( const typename ConfiguratorType::InitType &Initializer,
                                  const aol::Vector<RealType> &U0,
                                  const RealType Gamma,
                                  const RealType Delta,
                                  const RealType ThresholdLower,
                                  const RealType ThresholdUpper,
                                  const RealType ThresholdGradient,
                                  const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _gamma( Gamma ),
  _delta( Delta ),
  _thresholdLower( ThresholdLower ),
  _thresholdUpper( ThresholdUpper ),
  _thresholdGradient( ThresholdGradient ),
  _heavisideFunction( HeavisideFunction )
  {
  }
  virtual void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const{
    //! phi-derivative of perimiter length term
    qc::MCMStiffOp<ConfiguratorType> mcmStiffOp( _grid, aol::ONTHEFLY, 0.1 );
    mcmStiffOp.setImageReference( Arg );
    mcmStiffOp.apply( Arg, Dest );
    Dest *= _gamma;

    PhaseGradient<ConfiguratorType, HeavisideFunctionType> phaseGradient( _grid, _heavisideFunction, _delta, _thresholdLower, _thresholdUpper, _thresholdGradient );
    phaseGradient.applyAdd( _u0, Dest );
  }
  virtual void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const{
    aol::Vector<RealType> tmp( Arg, aol::STRUCT_COPY );
    apply( Arg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class CompleteCombinedEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _atomDistance;
  const RealType _gamma;
  const RealType _delta;
  const RealType _thresholdLower;
  const RealType _thresholdUpper;
  const RealType _thresholdGradient;
  const HeavisideFunctionType &_heavisideFunction;
  const CompletePhaseEnergy<ConfiguratorType, HeavisideFunctionType> _completePhaseEnergy;
public:
  CompleteCombinedEnergy( const typename ConfiguratorType::InitType &Initializer,
                          const aol::Vector<RealType> &U0,
                          const RealType AtomDistance,
                          const RealType Gamma,
                          const RealType Delta,
                          const RealType ThresholdLower,
                          const RealType ThresholdUpper,
                          const RealType ThresholdGradient,
                          const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _atomDistance( AtomDistance ),
  _gamma( Gamma ),
  _delta( Delta ),
  _thresholdLower( ThresholdLower ),
  _thresholdUpper( ThresholdUpper ),
  _thresholdGradient( ThresholdGradient ),
  _heavisideFunction( HeavisideFunction ),
  _completePhaseEnergy( _grid, _u0, _gamma, 1., _thresholdLower, _thresholdUpper, _thresholdGradient, _heavisideFunction )  {
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    //! Perimiter length term
    aol::HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType>
      heavisideLevelsetLengthEnergy( _grid, _heavisideFunction );
    aol::MultiVector<RealType> mPhi( 0, 0 );
    mPhi.appendReference( MArg[0] );
    heavisideLevelsetLengthEnergy.apply( mPhi, Dest);
    Dest *= _gamma;

    //! Alpha fitting term
    aol::MultiVector<RealType> mU0( 0, 0 );
    mU0.appendReference( _u0 );
    RotationEnergy<ConfiguratorType, HeavisideFunctionType>
      rotarionEnergySumAlpha1( _grid, _heavisideFunction, MArg[1][0], _atomDistance, _thresholdUpper, &AlphaIndicatorSum<ConfiguratorType, HeavisideFunctionType> );
    RotationEnergy<ConfiguratorType, HeavisideFunctionType>
      rotarionEnergySumAlpha2( _grid, _heavisideFunction, MArg[1][1], _atomDistance, _thresholdUpper, &AlphaIndicatorSum<ConfiguratorType, HeavisideFunctionType> );

    aol::HeavisideAndOneMinusHFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, RotationEnergy<ConfiguratorType, HeavisideFunctionType> >
      heavisideRotationEnergySum( _grid, MArg[0], _heavisideFunction, rotarionEnergySumAlpha1, rotarionEnergySumAlpha2 );

    //heavisideRotationEnergySum.applyAdd( mU0, Dest );

    aol::OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, aol::HeavisideAndOneMinusHFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, RotationEnergy<ConfiguratorType, HeavisideFunctionType> > >
      oneMinusHeavisidePhaseHeavisideRotationEnergy( _grid, MArg[2], _heavisideFunction, heavisideRotationEnergySum );
    oneMinusHeavisidePhaseHeavisideRotationEnergy.applyAdd( mU0, Dest );

    //! Phase fitting term
    Dest *= 1./_delta;
    _completePhaseEnergy.applyAdd( MArg[2], Dest);
    Dest *= _delta;
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( MArg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfCompleteCombinedEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _atomDistance;
  const RealType _gamma;
  const RealType _delta;
  const RealType _thresholdLower;
  const RealType _thresholdUpper;
  const RealType _thresholdGradient;
  const HeavisideFunctionType &_heavisideFunction;
  const VariationOfCompletePhaseEnergy<ConfiguratorType, HeavisideFunctionType> _variationOfCompletePhaseEnergy;
public:
  VariationOfCompleteCombinedEnergy( const typename ConfiguratorType::InitType &Initializer,
                                     const aol::Vector<RealType> &U0,
                                     const RealType AtomDistance,
                                     const RealType Gamma,
                                     const RealType Delta,
                                     const RealType ThresholdLower,
                                     const RealType ThresholdUpper,
                                     const RealType ThresholdGradient,
                                     const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _atomDistance( AtomDistance ),
  _gamma( Gamma ),
  _delta( Delta ),
  _thresholdLower( ThresholdLower ),
  _thresholdUpper( ThresholdUpper ),
  _thresholdGradient( ThresholdGradient ),
  _heavisideFunction( HeavisideFunction ),
  _variationOfCompletePhaseEnergy( _grid, _u0, _gamma, 1., _thresholdLower, _thresholdUpper, _thresholdGradient, _heavisideFunction )
  {
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{

    //! phi_alpha-derivative
    //! phi_alpha-derivative of perimiter length term
    qc::MCMStiffOp<ConfiguratorType> mcmStiffOp( _grid, aol::ONTHEFLY, 0.1 );
    mcmStiffOp.setImageReference( MArg[0] );
    mcmStiffOp.apply( MArg[0], MDest[0] );
    MDest[0] *= _gamma;

    aol::LinCombOp<aol::Vector<RealType> > phiAlphaGradient;
    //! phi_alpha-derivative of alpha fitting term
    RotationGradient<ConfiguratorType, HeavisideFunctionType>
      rotationGradientSumAlpha1( _grid, _heavisideFunction, MArg[1][0], _atomDistance, _thresholdUpper, &AlphaIndicatorSum<ConfiguratorType, HeavisideFunctionType> );
    RotationGradient<ConfiguratorType, HeavisideFunctionType>
      rotationGradientSumAlpha2( _grid, _heavisideFunction, MArg[1][1], _atomDistance, _thresholdUpper, &AlphaIndicatorSum<ConfiguratorType, HeavisideFunctionType> );

    aol::OneMinusHeavisideFENonlinOpInterface<ConfiguratorType, HeavisideFunctionType, RotationGradient<ConfiguratorType, HeavisideFunctionType> >
      phiPhaseRotationGradientSumAlpha1 ( _grid, MArg[2], _heavisideFunction, rotationGradientSumAlpha1 );
    aol::OneMinusHeavisideFENonlinOpInterface<ConfiguratorType, HeavisideFunctionType, RotationGradient<ConfiguratorType, HeavisideFunctionType> >
      phiPhaseRotationGradientSumAlpha2 ( _grid, MArg[2], _heavisideFunction, rotationGradientSumAlpha2 );

    //phiAlphaGradient.appendReference(rotationGradientSumAlpha1);
    //phiAlphaGradient.appendReference(rotationGradientSumAlpha2, -1.);
    phiAlphaGradient.appendReference(phiPhaseRotationGradientSumAlpha1);
    phiAlphaGradient.appendReference(phiPhaseRotationGradientSumAlpha2, -1.);
    phiAlphaGradient.applyAdd( _u0, MDest[0] );
    //! end of phi_alpha-derivative

    //! phi_phase-derivative
    //! phi_phase-derivative of alpha fitting term
    aol::HeavySideOneMinusHFENonlinOpInterface<ConfiguratorType, HeavisideFunctionType, RotationGradient<ConfiguratorType, HeavisideFunctionType> >
      alphaHeavisideRotationGradientSumAlpha( _grid, MArg[0], _heavisideFunction, rotationGradientSumAlpha1, rotationGradientSumAlpha2 );
    alphaHeavisideRotationGradientSumAlpha.apply( _u0, MDest[2] );
    MDest[2] *= -1./_delta;

    //! phi_phase-derivative of phase fitting term
    _variationOfCompletePhaseEnergy.applyAdd( MArg[2], MDest[2] );
    MDest[2] *= _delta;
    //! end of phi_phase-derivative

    //! alpha-derivative of alpha fitting term
    aol::Scalar<RealType> tempScalar;
    aol::MultiVector<RealType> mU0( 0, 0 );
    mU0.appendReference( _u0 );

    AlphaDerivative<ConfiguratorType, HeavisideFunctionType>
      alpha1SumDerivative( _grid, _heavisideFunction,  MArg[1][0], _atomDistance, _thresholdUpper, AlphaIndicatorSumDerivative<ConfiguratorType, HeavisideFunctionType> );
    AlphaDerivative<ConfiguratorType, HeavisideFunctionType>
      alpha2SumDerivative( _grid, _heavisideFunction,  MArg[1][1], _atomDistance, _thresholdUpper, AlphaIndicatorSumDerivative<ConfiguratorType, HeavisideFunctionType> );

    aol::HeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, AlphaDerivative<ConfiguratorType, HeavisideFunctionType> >
      heavisideSumAlpha1Derivative( _grid, MArg[0], _heavisideFunction, alpha1SumDerivative );
    aol::OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, aol::HeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, AlphaDerivative<ConfiguratorType, HeavisideFunctionType> > >
      phiPhaseheavisideSumAlpha1Derivative( _grid, MArg[2], _heavisideFunction, heavisideSumAlpha1Derivative );

    aol::OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, AlphaDerivative<ConfiguratorType, HeavisideFunctionType> >
      oneMinusHeavisideSumAlpha2Derivative( _grid, MArg[0], _heavisideFunction, alpha2SumDerivative );
    aol::OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, aol::OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, AlphaDerivative<ConfiguratorType, HeavisideFunctionType> > >
      phiPhaseoneMinusHeavisideSumAlpha2Derivative( _grid, MArg[2], _heavisideFunction, oneMinusHeavisideSumAlpha2Derivative );

    //heavisideSumAlpha1Derivative.apply( mU0, tempScalar);
    phiPhaseheavisideSumAlpha1Derivative.apply( mU0, tempScalar);
    MDest[1][0] = tempScalar[0];
    //oneMinusHeavisideSumAlpha2Derivative.apply( mU0, tempScalar);
    phiPhaseoneMinusHeavisideSumAlpha2Derivative.apply( mU0, tempScalar);
    MDest[1][1] = tempScalar[0];

  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
    apply( MArg, tmp );
    MDest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class CompleteRotationEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _atomDistance;
  const RealType _gamma;
  const RealType _threshold;
  const HeavisideFunctionType &_heavisideFunction;
public:
  CompleteRotationEnergy( const typename ConfiguratorType::InitType &Initializer,
                          const aol::Vector<RealType> &U0,
                          const RealType AtomDistance,
                          const RealType Gamma,
                          const RealType Threshold,
                          const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _atomDistance( AtomDistance ),
  _gamma( Gamma ),
  _threshold( Threshold ),
  _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    //! Perimiter length term
    aol::HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType>
      heavisideLevelsetLengthEnergy( _grid, _heavisideFunction );
    aol::MultiVector<RealType> mPhi( 0, 0 );
    mPhi.appendReference( MArg[0] );
    heavisideLevelsetLengthEnergy.apply( mPhi, Dest);
    Dest *= _gamma;

    //! Alpha fitting term
    aol::MultiVector<RealType> mU0( 0, 0 );
    mU0.appendReference( _u0 );
    RotationEnergy<ConfiguratorType, HeavisideFunctionType>
      rotarionEnergySumAlpha1( _grid, _heavisideFunction, MArg[1][0], _atomDistance, _threshold, &AlphaIndicatorSum<ConfiguratorType, HeavisideFunctionType> );
    RotationEnergy<ConfiguratorType, HeavisideFunctionType>
      rotarionEnergySumAlpha2( _grid, _heavisideFunction, MArg[1][1], _atomDistance, _threshold, &AlphaIndicatorSum<ConfiguratorType, HeavisideFunctionType> );

    aol::HeavisideAndOneMinusHFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, RotationEnergy<ConfiguratorType, HeavisideFunctionType> >
      heavisideRotationEnergySum( _grid, MArg[0], _heavisideFunction, rotarionEnergySumAlpha1, rotarionEnergySumAlpha2 );

    heavisideRotationEnergySum.applyAdd( mU0, Dest );
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( MArg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfCompleteRotationEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _atomDistance;
  const RealType _gamma;
  const RealType _threshold;
  const HeavisideFunctionType &_heavisideFunction;
public:
  VariationOfCompleteRotationEnergy( const typename ConfiguratorType::InitType &Initializer,
                                     const aol::Vector<RealType> &U0,
                                     const RealType AtomDistance,
                                     const RealType Gamma,
                                     const RealType Threshold,
                                     const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _atomDistance( AtomDistance ),
  _gamma( Gamma ),
  _threshold( Threshold ),
  _heavisideFunction( HeavisideFunction )
  {
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{

    //! phi_alpha-derivative
    //! phi_alpha-derivative of perimiter length term
    qc::MCMStiffOp<ConfiguratorType> mcmStiffOp( _grid, aol::ONTHEFLY, 0.1 );
    mcmStiffOp.setImageReference( MArg[0] );
    mcmStiffOp.apply( MArg[0], MDest[0] );
    MDest[0] *= _gamma;

    aol::LinCombOp<aol::Vector<RealType> > phiAlphaGradient;
    //! phi_alpha-derivative of alpha fitting term
    RotationGradient<ConfiguratorType, HeavisideFunctionType>
      rotationGradientSumAlpha1( _grid, _heavisideFunction, MArg[1][0], _atomDistance, _threshold, &AlphaIndicatorSum<ConfiguratorType, HeavisideFunctionType> );
    RotationGradient<ConfiguratorType, HeavisideFunctionType>
      rotationGradientSumAlpha2( _grid, _heavisideFunction, MArg[1][1], _atomDistance, _threshold, &AlphaIndicatorSum<ConfiguratorType, HeavisideFunctionType> );

    phiAlphaGradient.appendReference(rotationGradientSumAlpha1);
    phiAlphaGradient.appendReference(rotationGradientSumAlpha2, -1.);
    phiAlphaGradient.applyAdd( _u0, MDest[0] );
    //! end of phi_alpha-derivative

    //! alpha-derivative of alpha fitting term
    aol::Scalar<RealType> tempScalar;
    aol::MultiVector<RealType> mU0( 0, 0 );
    mU0.appendReference( _u0 );

    AlphaDerivative<ConfiguratorType, HeavisideFunctionType>
      alpha1SumDerivative( _grid, _heavisideFunction,  MArg[1][0], _atomDistance, _threshold, AlphaIndicatorSumDerivative<ConfiguratorType, HeavisideFunctionType> );
    AlphaDerivative<ConfiguratorType, HeavisideFunctionType>
      alpha2SumDerivative( _grid, _heavisideFunction,  MArg[1][1], _atomDistance, _threshold, AlphaIndicatorSumDerivative<ConfiguratorType, HeavisideFunctionType> );

    aol::HeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, AlphaDerivative<ConfiguratorType, HeavisideFunctionType> >
      heavisideSumAlpha1Derivative( _grid, MArg[0], _heavisideFunction, alpha1SumDerivative );

    aol::OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, AlphaDerivative<ConfiguratorType, HeavisideFunctionType> >
      oneMinusHeavisideSumAlpha2Derivative( _grid, MArg[0], _heavisideFunction, alpha2SumDerivative );

    heavisideSumAlpha1Derivative.apply( mU0, tempScalar);
    MDest[1][0] = tempScalar[0];
    oneMinusHeavisideSumAlpha2Derivative.apply( mU0, tempScalar);
    MDest[1][1] = tempScalar[0];
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
    apply( MArg, tmp );
    MDest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, int numberOfLevelsetFunctions>
class CompleteRotationMultiEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const GrainParams<RealType> &_params;
  const HeavisideFunctionType &_heavisideFunction;
public:
  CompleteRotationMultiEnergy( const typename ConfiguratorType::InitType &Initializer,
                               const aol::Vector<RealType> &U0,
                               const GrainParams<RealType> &Params,
                               const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0( U0 ),
  _params( Params ),
  _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::checkMultiChanVeseStructure<RealType>(MArg, numberOfLevelsetFunctions);

    //! Perimiter length term
    aol::HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      heavisideLevelsetLengthEnergy( _grid, _heavisideFunction );
    heavisideLevelsetLengthEnergy.apply( MArg, Dest);
    Dest *= _params.gamma;

    RotationEnergyMulti<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      rotEMulti( _grid, _heavisideFunction, _params.atomDistance, _params.thresholdUpper, _u0 );
    rotEMulti.applyAdd( MArg, Dest );

  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( MArg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, int numberOfLevelsetFunctions>
class VariationOfCompleteRotationMultiEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const GrainParams<RealType> &_params;
  const HeavisideFunctionType &_heavisideFunction;
public:
  VariationOfCompleteRotationMultiEnergy( const typename ConfiguratorType::InitType &Initializer,
                                          const aol::Vector<RealType> &U0,
                                          const GrainParams<RealType> &Params,
                                          const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _params( Params ),
  _heavisideFunction( HeavisideFunction )
  {
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::checkMultiChanVeseStructure<RealType>(MArg, numberOfLevelsetFunctions);

    aol::MultiVector<RealType> levelsetFunctionsArg( 0, 0 );
    aol::MultiVector<RealType> levelsetFunctionsDest( 0, 0 );
    for( int i = 0; i < numberOfLevelsetFunctions; i++ ){
      levelsetFunctionsArg.appendReference( MArg[i] );
      levelsetFunctionsDest.appendReference( MDest[i] );
    }
    aol::MultiVector<RealType> alphaArg( 0, 0 );
    aol::MultiVector<RealType> alphaDest( 0, 0 );
    for( int i = 0; i < (1<<numberOfLevelsetFunctions); i++ )
    {
      alphaArg.appendReference( MArg[i+numberOfLevelsetFunctions] );
      alphaDest.appendReference( MDest[i+numberOfLevelsetFunctions] );
    }

    //! phi_alpha-derivative
    //! phi_alpha-derivative of perimiter length term
    aol::VariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, numberOfLevelsetFunctions>
      variationOfHeavisideLevelsetLengthEnergy( _grid, _params.epsilonMCM, _params.gamma );
    variationOfHeavisideLevelsetLengthEnergy.apply( levelsetFunctionsArg, levelsetFunctionsDest );

    RotationEnergyPhiVariation<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      phiVar( _grid, _heavisideFunction, _params.atomDistance, _params.thresholdUpper, _u0, alphaArg );

    phiVar.applyAdd( levelsetFunctionsArg, levelsetFunctionsDest );
    //! end of phi_alpha-derivative

    //! alpha-derivative of alpha fitting term
    RotationEnergyAlphaVariation<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      alphaVar( _grid, _heavisideFunction, _params.atomDistance, _params.thresholdUpper, _u0, levelsetFunctionsArg );
    alphaVar.apply( alphaArg, alphaDest );
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
    apply( MArg, tmp );
    MDest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, int numberOfLevelsetFunctions>
class CompleteElasticEnergyMultiWithAlphaSegmentation : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const GrainParams<RealType> &_params;
  const HeavisideFunctionType &_heavisideFunction;
public:
  CompleteElasticEnergyMultiWithAlphaSegmentation( const typename ConfiguratorType::InitType &Initializer,
                                                   const aol::Vector<RealType> &U0,
                                                   const GrainParams<RealType> &Params,
                                                   const HeavisideFunctionType &HeavisideFunction )
  : _grid( Initializer ),
    _u0( U0 ),
    _params( Params ),
    _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::checkMultiChanVeseStructure<RealType>(MArg, numberOfLevelsetFunctions,ConfiguratorType::Dim);

    // Get the deformation components from MArg
    aol::MultiVector<RealType> deformation(0, 0);
    for( int i = 0; i < ConfiguratorType::Dim; i++ )
      deformation.appendReference( MArg[i+numberOfLevelsetFunctions+(1<<numberOfLevelsetFunctions)] );

    qc::SymmetricLengthEnergy<ConfiguratorType> symmetricLengthEnergy( _grid );
    symmetricLengthEnergy.apply( deformation, Dest );
    Dest *= _params.lambda;

    //! Perimiter length term
    aol::HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      heavisideLevelsetLengthEnergy( _grid, _heavisideFunction );
    Dest /= _params.gamma;
    heavisideLevelsetLengthEnergy.applyAdd( MArg, Dest);
    Dest *= _params.gamma;

    RotationEnergyMultiWithDeformation<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      rotEMulti( _grid, _heavisideFunction, _params.atomDistance, _params.thresholdUpper, deformation, _u0 );
    rotEMulti.setNumOfAdditionalVectorComponentsInArgument( ConfiguratorType::Dim );
    rotEMulti.applyAdd( MArg, Dest );

  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( MArg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, int numberOfLevelsetFunctions>
class VariationOfCompleteElasticEnergyMultiWithAlphaSegmentation : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const GrainParams<RealType> &_params;
  const HeavisideFunctionType &_heavisideFunction;
public:
  VariationOfCompleteElasticEnergyMultiWithAlphaSegmentation( const typename ConfiguratorType::InitType &Initializer,
                                                              const aol::Vector<RealType> &U0,
                                                              const GrainParams<RealType> &Params,
                                                              const HeavisideFunctionType &HeavisideFunction )
  : _grid( Initializer ),
    _u0(U0),
    _params( Params ),
    _heavisideFunction( HeavisideFunction )
  {
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::checkMultiChanVeseStructure<RealType>(MArg, numberOfLevelsetFunctions,ConfiguratorType::Dim);

    aol::MultiVector<RealType> levelsetFunctionsArg( 0, 0 );
    aol::MultiVector<RealType> levelsetFunctionsDest( 0, 0 );
    for( int i = 0; i < numberOfLevelsetFunctions; i++ ){
      levelsetFunctionsArg.appendReference( MArg[i] );
      levelsetFunctionsDest.appendReference( MDest[i] );
    }
    aol::MultiVector<RealType> alphaArg( 0, 0 );
    aol::MultiVector<RealType> alphaDest( 0, 0 );
    for( int i = 0; i < (1<<numberOfLevelsetFunctions); i++ )
    {
      alphaArg.appendReference( MArg[i+numberOfLevelsetFunctions] );
      alphaDest.appendReference( MDest[i+numberOfLevelsetFunctions] );
    }
    aol::MultiVector<RealType> deformationArg( 0, 0 );
    aol::MultiVector<RealType> deformationDest( 0, 0 );
    for( int i = 0; i < ConfiguratorType::Dim; i++ )
    {
      deformationArg.appendReference( MArg[i+numberOfLevelsetFunctions+(1<<numberOfLevelsetFunctions)] );
      deformationDest.appendReference( MDest[i+numberOfLevelsetFunctions+(1<<numberOfLevelsetFunctions)] );
    }

    //! phi_alpha-derivative
    //! phi_alpha-derivative of perimiter length term
    aol::VariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, numberOfLevelsetFunctions>
      variationOfHeavisideLevelsetLengthEnergy( _grid, _params.epsilonMCM, _params.gamma );
    variationOfHeavisideLevelsetLengthEnergy.apply( levelsetFunctionsArg, levelsetFunctionsDest );

    RotationEnergyMultiWithDeformationPhiVariation<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      phiVar( _grid, _heavisideFunction, _params.atomDistance, _params.thresholdUpper, deformationArg, _u0, alphaArg );

    phiVar.applyAdd( levelsetFunctionsArg, levelsetFunctionsDest );
    //! end of phi_alpha-derivative

    //! alpha-derivative of alpha fitting term
    RotationEnergyMultiWithDeformationAlphaVariation<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
      alphaVar( _grid, _heavisideFunction, _params.atomDistance, _params.thresholdUpper, deformationArg, _u0, levelsetFunctionsArg );
    alphaVar.apply( alphaArg, alphaDest );

    //! psi-derivative
    qc::VariationOfSymmetricLengthEnergy<ConfiguratorType> variationOfSymmetricLengthEnergy( _grid );
    variationOfSymmetricLengthEnergy.apply( deformationArg, deformationDest );
    deformationDest *= _params.lambda;

    RotationEnergyMultiWithDeformationPsiVariation<ConfiguratorType, HeavisideFunctionType>
      psiVar( _grid, _heavisideFunction, _params.atomDistance, _params.thresholdUpper, _u0, levelsetFunctionsArg, alphaArg );
    psiVar.applyAdd( deformationArg, deformationDest );
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
    apply( MArg, tmp );
    MDest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class CompleteElasticEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _lambda;
  const RealType _threshold;
  const HeavisideFunctionType &_heavisideFunction;
public:
  CompleteElasticEnergy( const typename ConfiguratorType::InitType &Initializer,
                         const aol::Vector<RealType> &U0,
                         const RealType Alpha,
                         const RealType AtomDistance,
                         const RealType Lambda,
                         const RealType Threshold,
                         const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _alpha( Alpha ),
  _atomDistance( AtomDistance ),
  _lambda( Lambda ),
  _threshold( Threshold ),
  _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    qc::SymmetricLengthEnergy<ConfiguratorType> symmetricLengthEnergy( _grid );
    symmetricLengthEnergy.apply( MArg, Dest);
    Dest *= _lambda;

    RotationEnergyWithDeformation<ConfiguratorType, HeavisideFunctionType>
      rotationEnergyWithDeformation( _grid, _heavisideFunction, _alpha, _atomDistance, _threshold, _u0, &AlphaPsiIndicatorSum<ConfiguratorType, HeavisideFunctionType> );
    rotationEnergyWithDeformation.applyAdd( MArg, Dest );
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( MArg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfCompleteElasticEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _lambda;
  const RealType _threshold;
  const HeavisideFunctionType &_heavisideFunction;
public:
  VariationOfCompleteElasticEnergy( const typename ConfiguratorType::InitType &Initializer,
                                    const aol::Vector<RealType> &U0,
                                    const RealType Alpha,
                                    const RealType AtomDistance,
                                    const RealType Lambda,
                                    const RealType Threshold,
                                    const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _alpha( Alpha ),
  _atomDistance( AtomDistance ),
  _lambda( Lambda ),
  _threshold( Threshold ),
  _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    qc::VariationOfSymmetricLengthEnergy<ConfiguratorType> variationOfSymmetricLengthEnergy( _grid );
    variationOfSymmetricLengthEnergy.apply( MArg, MDest);
    MDest *= _lambda;

    VariationOfRotationEnergyWithDeformationWRTPsi<ConfiguratorType, HeavisideFunctionType>
      variationOfRotationEnergyWithDeformationWRTPsi( _grid, _heavisideFunction, _alpha, _atomDistance, _threshold, _u0 );
    variationOfRotationEnergyWithDeformationWRTPsi.applyAdd( MArg, MDest );
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
    apply( MArg, tmp );
    MDest += tmp;
  }
};


template <typename ConfiguratorType, typename HeavisideFunctionType>
class CompleteAlphaPsiElasticEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _atomDistance;
  const RealType _lambda;
  const RealType _threshold;
  const HeavisideFunctionType &_heavisideFunction;
public:
  CompleteAlphaPsiElasticEnergy( const typename ConfiguratorType::InitType &Initializer,
                                 const aol::Vector<RealType> &U0,
                                 const RealType AtomDistance,
                                 const RealType Lambda,
                                 const RealType Threshold,
                                 const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _atomDistance( AtomDistance ),
  _lambda( Lambda ),
  _threshold( Threshold ),
  _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::MultiVector<RealType> psi( 0, 0 );
    for( int i = 0; i < ConfiguratorType::Dim; i++ )
      psi.appendReference( MArg[i] );
    CompleteElasticEnergy<ConfiguratorType, HeavisideFunctionType>
      E( this->_grid, _u0, MArg[ConfiguratorType::Dim][0], _atomDistance, _lambda, _threshold, _heavisideFunction );
    E.apply( psi, Dest );
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( MArg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfCompleteAlphaPsiElasticEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _atomDistance;
  const RealType _lambda;
  const RealType _threshold;
  const HeavisideFunctionType &_heavisideFunction;
public:
  VariationOfCompleteAlphaPsiElasticEnergy( const typename ConfiguratorType::InitType &Initializer,
                                            const aol::Vector<RealType> &U0,
                                            const RealType AtomDistance,
                                            const RealType Lambda,
                                            const RealType Threshold,
                                            const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _atomDistance( AtomDistance ),
  _lambda( Lambda ),
  _threshold( Threshold ),
  _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    // psi-variation
    aol::MultiVector<RealType> psi( 0, 0 );
    aol::MultiVector<RealType> variationOfPsi( 0, 0 );
    for( int i = 0; i < ConfiguratorType::Dim; i++ ){
      psi.appendReference( MArg[i] );
      variationOfPsi.appendReference( MDest[i] );
    }
    VariationOfCompleteElasticEnergy<ConfiguratorType, HeavisideFunctionType>
      DE( _grid, _u0, MArg[ConfiguratorType::Dim][0], _atomDistance, _lambda, _threshold, _heavisideFunction );
    DE.apply( psi, variationOfPsi );

    // alpha-variation
    aol::Scalar<RealType> tmp;
    VariationOfRotationEnergyWithDeformationWRTAlpha<ConfiguratorType, HeavisideFunctionType>
      varWRTAlpha( _grid, _heavisideFunction, MArg[ConfiguratorType::Dim][0], _atomDistance, _threshold, _u0 );
    varWRTAlpha.apply( psi, tmp );
    MDest[ConfiguratorType::Dim][0] = tmp[0];
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
    apply( MArg, tmp );
    MDest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class CompleteElasticEnergyWithPsiSegmentation : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _gamma;
  const RealType _lambda;
  const RealType _threshold;
  const HeavisideFunctionType &_heavisideFunction;
public:
  CompleteElasticEnergyWithPsiSegmentation( const typename ConfiguratorType::InitType &Initializer,
                                            const aol::Vector<RealType> &U0,
                                            const RealType Alpha,
                                            const RealType AtomDistance,
                                            const RealType Gamma,
                                            const RealType Lambda,
                                            const RealType Threshold,
                                            const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _alpha( Alpha ),
  _atomDistance( AtomDistance ),
  _gamma( Gamma ),
  _lambda( Lambda ),
  _threshold( Threshold ),
  _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::MultiVector<RealType> psi1( 0, 0);
    aol::MultiVector<RealType> psi2( 0, 0);
    for( int i = 0; i < ConfiguratorType::Dim; i++ ){
      psi1.appendReference( MArg[i] );
      psi2.appendReference( MArg[i+ConfiguratorType::Dim] );
    }
    qc::SymmetricLengthEnergy<ConfiguratorType> symmetricLengthEnergy( _grid );
    symmetricLengthEnergy.apply( psi1, Dest );
    symmetricLengthEnergy.applyAdd( psi2, Dest );
    Dest *= _lambda;

    RotationEnergyWithDeformation<ConfiguratorType, HeavisideFunctionType>
      rotationEnergyWithDeformation( _grid, _heavisideFunction, _alpha, _atomDistance, _threshold, _u0, &AlphaPsiIndicatorSum<ConfiguratorType, HeavisideFunctionType> );
      
    aol::HeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, RotationEnergyWithDeformation<ConfiguratorType, HeavisideFunctionType> >
      heavisideRotationEnergyWithDeformation( _grid, MArg[2*ConfiguratorType::Dim],_heavisideFunction, rotationEnergyWithDeformation);
    aol::OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, RotationEnergyWithDeformation<ConfiguratorType, HeavisideFunctionType> >
      oneMinusHeavisideRotationEnergyWithDeformation( _grid, MArg[2*ConfiguratorType::Dim],_heavisideFunction, rotationEnergyWithDeformation);

    heavisideRotationEnergyWithDeformation.applyAdd( psi1, Dest );
    oneMinusHeavisideRotationEnergyWithDeformation.applyAdd( psi2, Dest );
    
    //! Perimiter length term
    aol::HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType >
      heavisideLevelsetLengthEnergy( _grid, _heavisideFunction );
    aol::MultiVector<RealType> mPhi( 0, 0 );
    mPhi.appendReference( MArg[2*ConfiguratorType::Dim] );
    Dest /= _gamma;
    heavisideLevelsetLengthEnergy.applyAdd( mPhi, Dest);
    Dest *= _gamma;

  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( MArg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfCompleteElasticEnergyWithPsiSegmentation : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _alpha;
  const RealType _atomDistance;
  const RealType _gamma;
  const RealType _lambda;
  const RealType _threshold;
  const HeavisideFunctionType &_heavisideFunction;
public:
  VariationOfCompleteElasticEnergyWithPsiSegmentation( const typename ConfiguratorType::InitType &Initializer,
                                                       const aol::Vector<RealType> &U0,
                                                       const RealType Alpha,
                                                       const RealType AtomDistance,
                                                       const RealType Gamma,
                                                       const RealType Lambda,
                                                       const RealType Threshold,
                                                       const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _alpha( Alpha ),
  _atomDistance( AtomDistance ),
  _gamma( Gamma ),
  _lambda( Lambda ),
  _threshold( Threshold ),
  _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> psi1( 0, 0);
    aol::MultiVector<RealType> psi2( 0, 0);
    aol::MultiVector<RealType> varOfPsi1( 0, 0);
    aol::MultiVector<RealType> varOfPsi2( 0, 0);
    for( int i = 0; i < ConfiguratorType::Dim; i++ ){
      psi1.appendReference( MArg[i] );
      psi2.appendReference( MArg[i+ConfiguratorType::Dim] );
      varOfPsi1.appendReference( MDest[i] );
      varOfPsi2.appendReference( MDest[i+ConfiguratorType::Dim] );
    }

    // psi1+psi2-variation
    qc::VariationOfSymmetricLengthEnergy<ConfiguratorType> variationOfSymmetricLengthEnergy( _grid );
    variationOfSymmetricLengthEnergy.apply( psi1, varOfPsi1 );
    variationOfSymmetricLengthEnergy.apply( psi2, varOfPsi2 );
    varOfPsi1 *= _lambda;
    varOfPsi2 *= _lambda;

    VariationOfRotationEnergyWithDeformationAndLevelsetWRTPsi<ConfiguratorType, HeavisideFunctionType>
      varOfRotEnergyWRTPsi1( _grid, _heavisideFunction, _alpha, _atomDistance, _threshold, _u0, MArg[2*ConfiguratorType::Dim], false );
    varOfRotEnergyWRTPsi1.applyAdd( psi1, varOfPsi1 );
    VariationOfRotationEnergyWithDeformationAndLevelsetWRTPsi<ConfiguratorType, HeavisideFunctionType>
      varOfRotEnergyWRTPsi2( _grid, _heavisideFunction, _alpha, _atomDistance, _threshold, _u0, MArg[2*ConfiguratorType::Dim], true );
    varOfRotEnergyWRTPsi2.applyAdd( psi2, varOfPsi2 );
    
    // phi-variation
    qc::MCMStiffOp<ConfiguratorType> mcmStiffOp( _grid, aol::ONTHEFLY, 0.1 );
    mcmStiffOp.setImageReference( MArg[2*ConfiguratorType::Dim] );
    mcmStiffOp.apply( MArg[2*ConfiguratorType::Dim], MDest[2*ConfiguratorType::Dim] );
    MDest[2*ConfiguratorType::Dim] *= _gamma;

    VariationOfRotationEnergyWithDeformationAndLevelsetWRTPhi<ConfiguratorType, HeavisideFunctionType>
      psi1PartOfVarOfRotEnergyWRTPhi( _grid, _heavisideFunction, _alpha, _atomDistance, _threshold, psi1 );
    psi1PartOfVarOfRotEnergyWRTPhi.applyAdd( _u0, MDest[2*ConfiguratorType::Dim] );
    VariationOfRotationEnergyWithDeformationAndLevelsetWRTPhi<ConfiguratorType, HeavisideFunctionType>
      psi2PartOfVarOfRotEnergyWRTPhi( _grid, _heavisideFunction, _alpha, _atomDistance, _threshold, psi2 );
    MDest[2*ConfiguratorType::Dim] *= -1;
    psi2PartOfVarOfRotEnergyWRTPhi.applyAdd( _u0, MDest[2*ConfiguratorType::Dim] );
    MDest[2*ConfiguratorType::Dim] *= -1;
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
    apply( MArg, tmp );
    MDest += tmp;
  }
};

// marg[0]-marg[d-1] = psi, marg[d] = alpha, marg[d+1] = phi
template <typename ConfiguratorType, typename HeavisideFunctionType>
class CompleteElasticEnergyWithAlphaSegmentation : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _atomDistance;
  const RealType _gamma;
  const RealType _lambda;
  const RealType _threshold;
  const HeavisideFunctionType &_heavisideFunction;
public:
  CompleteElasticEnergyWithAlphaSegmentation( const typename ConfiguratorType::InitType &Initializer,
                                              const aol::Vector<RealType> &U0,
                                              const RealType AtomDistance,
                                              const RealType Gamma,
                                              const RealType Lambda,
                                              const RealType Threshold,
                                              const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _atomDistance( AtomDistance ),
  _gamma( Gamma ),
  _lambda( Lambda ),
  _threshold( Threshold ),
  _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::MultiVector<RealType> psi( 0, 0);
    for( int i = 0; i < ConfiguratorType::Dim; i++ ){
      psi.appendReference( MArg[i] );
    }
    qc::SymmetricLengthEnergy<ConfiguratorType> symmetricLengthEnergy( _grid );
    symmetricLengthEnergy.apply( psi, Dest );
    Dest *= _lambda;

    RotationEnergyWithDeformation<ConfiguratorType, HeavisideFunctionType>
      rotationEnergyWithDeformationAlpha1( _grid, _heavisideFunction, MArg[ConfiguratorType::Dim][0], _atomDistance, _threshold, _u0, &AlphaPsiIndicatorSum<ConfiguratorType, HeavisideFunctionType> );
    RotationEnergyWithDeformation<ConfiguratorType, HeavisideFunctionType>
      rotationEnergyWithDeformationAlpha2( _grid, _heavisideFunction, MArg[ConfiguratorType::Dim][1], _atomDistance, _threshold, _u0, &AlphaPsiIndicatorSum<ConfiguratorType, HeavisideFunctionType> );

    aol::HeavisideAndOneMinusHFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, RotationEnergyWithDeformation<ConfiguratorType, HeavisideFunctionType> >
      heavisideAndOneMinusHRotationEnergyWithDeformation( _grid, MArg[ConfiguratorType::Dim+1],_heavisideFunction, rotationEnergyWithDeformationAlpha1, rotationEnergyWithDeformationAlpha2 );

    heavisideAndOneMinusHRotationEnergyWithDeformation.applyAdd( psi, Dest );

    //! Perimiter length term
    aol::HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType >
      heavisideLevelsetLengthEnergy( _grid, _heavisideFunction );
    aol::MultiVector<RealType> mPhi( 0, 0 );
    mPhi.appendReference( MArg[ConfiguratorType::Dim+1] );
    Dest /= _gamma;
    heavisideLevelsetLengthEnergy.applyAdd( mPhi, Dest);
    Dest *= _gamma;

  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( MArg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfCompleteElasticEnergyWithAlphaSegmentation : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_u0;
  const RealType _atomDistance;
  const RealType _gamma;
  const RealType _lambda;
  const RealType _threshold;
  const HeavisideFunctionType &_heavisideFunction;
public:
  VariationOfCompleteElasticEnergyWithAlphaSegmentation( const typename ConfiguratorType::InitType &Initializer,
                                                         const aol::Vector<RealType> &U0,
                                                         const RealType AtomDistance,
                                                         const RealType Gamma,
                                                         const RealType Lambda,
                                                         const RealType Threshold,
                                                         const HeavisideFunctionType &HeavisideFunction ):
  _grid( Initializer ),
  _u0(U0),
  _atomDistance( AtomDistance ),
  _gamma( Gamma ),
  _lambda( Lambda ),
  _threshold( Threshold ),
  _heavisideFunction( HeavisideFunction ){
  }
  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> psi( 0, 0);
    aol::MultiVector<RealType> varOfPsi( 0, 0);
    for( int i = 0; i < ConfiguratorType::Dim; i++ ){
      psi.appendReference( MArg[i] );
      varOfPsi.appendReference( MDest[i] );
    }

    // psi-variation
    qc::VariationOfSymmetricLengthEnergy<ConfiguratorType> variationOfSymmetricLengthEnergy( _grid );
    variationOfSymmetricLengthEnergy.apply( psi, varOfPsi );
    varOfPsi *= _lambda;

    VariationOfRotationEnergyWithDeformationAndLevelsetWRTPsi<ConfiguratorType, HeavisideFunctionType>
      varOfRotEnergyAlpha1WRTPsi( _grid, _heavisideFunction, MArg[ConfiguratorType::Dim][0], _atomDistance, _threshold, _u0, MArg[ConfiguratorType::Dim+1], false );
    varOfRotEnergyAlpha1WRTPsi.applyAdd( psi, varOfPsi );
    VariationOfRotationEnergyWithDeformationAndLevelsetWRTPsi<ConfiguratorType, HeavisideFunctionType>
      varOfRotEnergyAlpha2WRTPsi( _grid, _heavisideFunction, MArg[ConfiguratorType::Dim][1], _atomDistance, _threshold, _u0, MArg[ConfiguratorType::Dim+1], true );
    varOfRotEnergyAlpha2WRTPsi.applyAdd( psi, varOfPsi );

    // alpha-variation
    aol::Scalar<RealType> tmp;
    VariationOfRotationEnergyWithDeformationWRTAlpha<ConfiguratorType, HeavisideFunctionType>
      varWRTAlpha1( _grid, _heavisideFunction, MArg[ConfiguratorType::Dim][0], _atomDistance, _threshold, _u0 );
    aol::HeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, VariationOfRotationEnergyWithDeformationWRTAlpha<ConfiguratorType, HeavisideFunctionType> >
      heavisideVarWRTAlpha1( _grid, MArg[ConfiguratorType::Dim+1],_heavisideFunction, varWRTAlpha1);
    heavisideVarWRTAlpha1.apply( psi, tmp );
    MDest[ConfiguratorType::Dim][0] = tmp[0];

    VariationOfRotationEnergyWithDeformationWRTAlpha<ConfiguratorType, HeavisideFunctionType>
      varWRTAlpha2( _grid, _heavisideFunction, MArg[ConfiguratorType::Dim][1], _atomDistance, _threshold, _u0 );
    aol::OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, VariationOfRotationEnergyWithDeformationWRTAlpha<ConfiguratorType, HeavisideFunctionType> >
      oneMinusHeavisideVarWRTAlpha2( _grid, MArg[ConfiguratorType::Dim+1],_heavisideFunction, varWRTAlpha2);
    oneMinusHeavisideVarWRTAlpha2.apply( psi, tmp );
    MDest[ConfiguratorType::Dim][1] = tmp[0];

    // phi-variation
    qc::MCMStiffOp<ConfiguratorType> mcmStiffOp( _grid, aol::ONTHEFLY, 0.1 );
    mcmStiffOp.setImageReference( MArg[ConfiguratorType::Dim+1] );
    mcmStiffOp.apply( MArg[ConfiguratorType::Dim+1], MDest[ConfiguratorType::Dim+1] );
    MDest[ConfiguratorType::Dim+1] *= _gamma;

    VariationOfRotationEnergyWithDeformationAndLevelsetWRTPhi<ConfiguratorType, HeavisideFunctionType>
      alpha1PartOfVarOfRotEnergyWRTPhi( _grid, _heavisideFunction, MArg[ConfiguratorType::Dim][0], _atomDistance, _threshold, psi );
    alpha1PartOfVarOfRotEnergyWRTPhi.applyAdd( _u0, MDest[ConfiguratorType::Dim+1] );
    VariationOfRotationEnergyWithDeformationAndLevelsetWRTPhi<ConfiguratorType, HeavisideFunctionType>
      alpha2PartOfVarOfRotEnergyWRTPhi( _grid, _heavisideFunction, MArg[ConfiguratorType::Dim][1], _atomDistance, _threshold, psi );
    MDest[ConfiguratorType::Dim+1] *= -1;
    alpha2PartOfVarOfRotEnergyWRTPhi.applyAdd( _u0, MDest[ConfiguratorType::Dim+1] );
    MDest[ConfiguratorType::Dim+1] *= -1;
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
    apply( MArg, tmp );
    MDest += tmp;
  }
};

/*
 * Class for various types of grain segmentation.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class CrystalSegmentation : public qc::QuocTimestepSaver<ConfiguratorType>
{
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  // ******** some more parameters ************
    
  //! Image at atomic scale we are processing.
  ArrayType _u0;
  const typename ConfiguratorType::InitType _grid;
  aol::ParameterParser &_parser;
  const GrainParams<RealType> _params;
    
public:
  CrystalSegmentation( aol::ParameterParser &Parser ) :
    qc::QuocTimestepSaver<ConfiguratorType>( Parser.getInt("saveOffset"), Parser.getInt("numSaveFirst"), "PleaseRenameMe", true ),
    _u0( Parser.getString( "reference" ).c_str() ), 
    _grid( _u0.getSize () ),
    _parser( Parser ),
    _params( Parser )
  {    
    char fn_reference[1024];
    Parser.getString( "reference", fn_reference );
    char filename[1024];
    strncpy(filename, aol::getFileName(fn_reference), 1023);

    _u0 /= _u0.getMaxValue();

    Parser.getString( "saveDirectory", fn_reference );
    char fn[1024];

    sprintf( fn, "%s_%s_%s", fn_reference, filename, _params.getDirName() );
    sprintf( fn_reference, "%s/", fn );
    this->setSaveDirectory( fn_reference );
    aol::makeDirectory( fn );
    Parser.dumpToFile( "parameter-dump.txt", fn_reference );
  }

  virtual ~CrystalSegmentation() {
  }

  /*
   * Loads the node values of the levelset functions from disk. The parameter
   * ParameterName in the parameter file has to specify the filenames of the
   * components.
   *
   * Returns true if loading was sucessful, false otherwise.
   *
   */
  bool loadLevelsetFunctions( aol::MultiVector<RealType> &LevelsetFunctions, const char *ParameterName ) const {
    if( _parser.hasVariable ( ParameterName ) ){
      int dimSize = _parser.getDimSize (ParameterName, 0 );
      char filename[1024];
      if ( LevelsetFunctions.numComponents() != dimSize ){
        cerr << "Can't load levelsetfunctions: LevelsetFunctions.numComponents() != dimSize!\n";
        return false;
      }
      for ( int i = 0; i < dimSize; i ++ ){
        _parser.getString ( ParameterName, filename, i );
        ArrayType tempArray( LevelsetFunctions[i], _grid );
        tempArray.load( filename );
      }
      return true;
    }
    else
      return false;
  }

  void computeTimesteps( )
  {
    enum MODE {
        SEGMENT_PHASE,
        SEGMENT_PHASE_UNCONSTRAINED_CONVEX,
        SEGMENT_GRAINS,
        SEGMENT_GRAINS_MULTI,
        SEGMENT_GRAINS_MULTI_WITH_ALPHA_AND_GET_DEFORMATION,
        SEGMENT_GRAY_VALUES_MULTI,
        SEGMENT_GRAINS_AND_PHASE,
        GET_DEFORMATION,
        SEGMENT_GRAINS_WITH_DEFORMATION,
        SEGMENT_GRAINS_WITH_ALPHA_AND_GET_DEFORMATION
    };
    MODE mode = SEGMENT_GRAINS_MULTI_WITH_ALPHA_AND_GET_DEFORMATION;

    ArrayType phaseLevelsetFunction( _grid );
    ArrayType alphaLevelsetFunction( _grid );
    aol::Vector<RealType> alpha( 2 );
    alpha[0] = 0.;
    alpha[1] = -1.*aol::NumberTrait<RealType>::pi/4.;
//    alpha[0] = -1.*M_PI/4.;
//    alpha[1] = 0.;
//    alpha[0] = -1.*M_PI/8.;
//    alpha[1] = -1.*M_PI/8.;
//    alpha[0] = 3.*M_PI/180.;
//    alpha[1] = -37.*M_PI/180.;
//    alpha[0] = 27.*M_PI/180.;
//    alpha[1] = -26.5*M_PI/180.;
//    alpha[1] = -11.*M_PI/180.;
//    alpha[0] = 4.*M_PI/180.;
//    alpha[1] = 4.6*M_PI/180.;
//    alpha[0] = -28.*M_PI/180.;
//    alpha[0] = -0.0413194;
//    alpha[1] = -0.293298;



    // Initialize levelset function
    {
      qc::DataGenerator<ConfiguratorType> generator( _grid );

      //generator.generateLineLevelset ( phaseLevelsetFunction, 0, 0.5 );
      //phaseLevelsetFunction *= -1;
      //generator.generateCircleLevelset ( alphaLevelsetFunction, 0.5, 0.5, 0.5 );
      generator.generateCircleLevelset ( phaseLevelsetFunction, 0.5, 0.5, 0.5 );
      phaseLevelsetFunction *= -1;

      HeavisideFunctionType heavisideFunction( _params.epsilon );
      generator.generateLineLevelset( alphaLevelsetFunction, 0, 0.5 );
/*
      aol::MultiVector<RealType> testMV( 0, 0 );
      testMV.appendReference( alphaLevelsetFunction );
      const typename ConfiguratorType::ElementType element;
      aol::HeavisideFunctionProduct<ConfiguratorType, HeavisideFunctionType> temp( _grid, heavisideFunction, testMV, element, 0);
      aol::PowerSetIterator myit(1);
      cerr << temp.evaluate( myit ) << endl;
      cerr << temp.evaluateDerivative( 0, myit );
      aol::callSystemPauseIfNecessaryOnPlatform();
*/
      switch( mode ) {
      case SEGMENT_PHASE:
        {
          CompletePhaseEnergy<ConfiguratorType, HeavisideFunctionType>
            E( _grid, _u0, _params.gamma, _params.delta, _params.thresholdLower, _params.thresholdUpper, _params.thresholdGradient, heavisideFunction );
        
          VariationOfCompletePhaseEnergy<ConfiguratorType, HeavisideFunctionType>
            DE( _grid, _u0, _params.gamma, _params.delta, _params.thresholdLower, _params.thresholdUpper, _params.thresholdGradient, heavisideFunction );

          LevelsetGradientDescent<ConfiguratorType, HeavisideFunctionType, aol::Vector<RealType> >
            gradientDescent_solver( _grid, E, DE, heavisideFunction, _u0, _params.numSteps, aol::ZOTrait<RealType>::one );
          gradientDescent_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "phi", this->_saveDirectory );

          aol::Vector<RealType> tmp( _grid.getNumberOfNodes() );
    
          gradientDescent_solver.apply( phaseLevelsetFunction, tmp );
          phaseLevelsetFunction = tmp;
        }
        break;
      case SEGMENT_PHASE_UNCONSTRAINED_CONVEX:
        {
          SolidLiquidPhaseMSSegmentor<ConfiguratorType, HeavisideFunctionType>
            solidLiquidPhaseMSSegmentor ( _grid, _u0, heavisideFunction, _params.gamma, _params.delta, _params.thresholdLower, _params.thresholdUpper, _params.thresholdGradient, true );

          solidLiquidPhaseMSSegmentor.segment ( phaseLevelsetFunction );

          qc::QuocTimestepSaver<ConfiguratorType> timestepSaver( _parser.getInt( "saveOffset"), _parser.getInt( "numSaveFirst"), "phi", true );
          timestepSaver.setSaveDirectory ( this->_saveDirectory );
          timestepSaver.setSaveName ( "u" );
          timestepSaver.savePNG ( phaseLevelsetFunction, _grid, 0 );
          timestepSaver.saveTimestepBZ2 ( 0, phaseLevelsetFunction, _grid );
          timestepSaver.setSaveName ( "seg" );
          timestepSaver.saveTimestepPNGDrawRedIsoline( -1, 1, 0.5, phaseLevelsetFunction, _u0, _grid );
        }
        break;
      case SEGMENT_GRAINS_AND_PHASE:
        {
          CompleteCombinedEnergy<ConfiguratorType, HeavisideFunctionType>
            E( _grid, _u0, _params.atomDistance, _params.gamma, _params.delta, _params.thresholdLower, _params.thresholdUpper, _params.thresholdGradient, heavisideFunction );

          VariationOfCompleteCombinedEnergy<ConfiguratorType, HeavisideFunctionType>
            DE( _grid, _u0, _params.atomDistance, _params.gamma, _params.delta, _params.thresholdLower, _params.thresholdUpper, _params.thresholdGradient, heavisideFunction );

          LevelsetGradientDescent<ConfiguratorType, HeavisideFunctionType, aol::MultiVector<RealType> >
            gradientDescent_solver( _grid, E, DE, heavisideFunction, _u0, _params.numSteps, aol::ZOTrait<RealType>::one );
          gradientDescent_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "phi", this->_saveDirectory );

          aol::MultiVector<RealType> levelsetFunctionsAndAlpha( 0, 0 );
          levelsetFunctionsAndAlpha.appendReference( alphaLevelsetFunction );
          levelsetFunctionsAndAlpha.appendReference( alpha );
          levelsetFunctionsAndAlpha.appendReference( phaseLevelsetFunction );

          aol::MultiVector<RealType> mtmp( levelsetFunctionsAndAlpha, aol::STRUCT_COPY );
    
          gradientDescent_solver.apply( levelsetFunctionsAndAlpha, mtmp );
          levelsetFunctionsAndAlpha = mtmp;

        }
        break;
      case SEGMENT_GRAINS:
        {
          //CompleteRotationMultiEnergy<ConfiguratorType, HeavisideFunctionType>
          CompleteRotationEnergy<ConfiguratorType, HeavisideFunctionType>
            E( _grid, _u0, _params.atomDistance, _params.gamma, _params.thresholdUpper, heavisideFunction );

          //VariationOfCompleteRotationMultiEnergy<ConfiguratorType, HeavisideFunctionType>
          VariationOfCompleteRotationEnergy<ConfiguratorType, HeavisideFunctionType>
            DE( _grid, _u0, _params.atomDistance, _params.gamma, _params.thresholdUpper, heavisideFunction );

          LevelsetGradientDescent<ConfiguratorType, HeavisideFunctionType, aol::MultiVector<RealType> >
            gradientDescent_solver( _grid, E, DE, heavisideFunction, _u0, _params.numSteps, aol::ZOTrait<RealType>::one );
          gradientDescent_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "phi", this->_saveDirectory );
          aol::EnergyDescentPlotter<LevelsetGradientDescent<ConfiguratorType, HeavisideFunctionType, aol::MultiVector<RealType> > >
            energyDescentPlotter( "energy.txt", "sigma.txt", this->_saveDirectory, gradientDescent_solver );

          aol::MultiVector<RealType> levelsetFunctionAndAlpha( 0, 0 );
          levelsetFunctionAndAlpha.appendReference( alphaLevelsetFunction );
          levelsetFunctionAndAlpha.appendReference( alpha );

          aol::MultiVector<RealType> mtmp( levelsetFunctionAndAlpha, aol::STRUCT_COPY );
    
          gradientDescent_solver.apply( levelsetFunctionAndAlpha, mtmp );
          levelsetFunctionAndAlpha = mtmp;

          energyDescentPlotter.plot( "EnergyPlot" );
        }
        break;
      case SEGMENT_GRAINS_MULTI:
        {
          const int numberOfLevelsetFunctions = 2;
          CompleteRotationMultiEnergy<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
            E( _grid, _u0, _params, heavisideFunction );

          VariationOfCompleteRotationMultiEnergy<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
            DE( _grid, _u0, _params, heavisideFunction );

          qc::ChanVeseGradientDescent<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
            gradientDescent_solver( _grid, E, DE, heavisideFunction, _u0, _params.numSteps, aol::ZOTrait<RealType>::one );
          gradientDescent_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "phi", this->_saveDirectory );
          //gradientDescent_solver.setTauMin( pow(2.,-10.) );

          aol::MultiVector<RealType> levelsetFunctionsAndParameters( numberOfLevelsetFunctions, _grid.getNumberOfNodes() );

          if( !loadLevelsetFunctions( levelsetFunctionsAndParameters, "phi" ) ){
            generator.generateChanVeseInitialization( levelsetFunctionsAndParameters );
            //generator.generateLineLevelset( levelsetFunctionsAndParameters[0], 1, 0.7 );
            //generator.generateLineLevelset( levelsetFunctionsAndParameters[1], 1, 0.5 );
          }
          aol::MultiVector<RealType> parameters( 1<<numberOfLevelsetFunctions, 1 );
          levelsetFunctionsAndParameters.appendReference( parameters );
          if ( numberOfLevelsetFunctions > 1 ){
            parameters[0][0] = 0.;
            parameters[1][0] = -0.5*aol::NumberTrait<RealType>::pi/4.;
            parameters[2][0] = 1.*aol::NumberTrait<RealType>::pi/4.;
            parameters[3][0] = -1.*aol::NumberTrait<RealType>::pi/4.;
          }

          aol::MultiVector<RealType> mtmp( levelsetFunctionsAndParameters, aol::STRUCT_COPY );
    
          gradientDescent_solver.apply( levelsetFunctionsAndParameters, mtmp );
          levelsetFunctionsAndParameters = mtmp;

        }
        break;
      case SEGMENT_GRAINS_MULTI_WITH_ALPHA_AND_GET_DEFORMATION:
        {
          const int numberOfLevelsetFunctions = 1;
          CompleteElasticEnergyMultiWithAlphaSegmentation<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
            E( _grid, _u0, _params, heavisideFunction );

          VariationOfCompleteElasticEnergyMultiWithAlphaSegmentation<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
            DE( _grid, _u0, _params, heavisideFunction );

          RotationEnergyMultiWithDeformationAlphaSegmentationGradientDescent<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
            gradientDescent_solver( _grid, E, DE, heavisideFunction, _u0, _params.numSteps, aol::ZOTrait<RealType>::one );
          gradientDescent_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "phi", this->_saveDirectory );

          aol::MultiVector<RealType> levelsetFunctionsParametersAndDeformation( numberOfLevelsetFunctions, _grid.getNumberOfNodes() );

          if( !loadLevelsetFunctions( levelsetFunctionsParametersAndDeformation, "phi" ) ){
            generator.generateChanVeseInitialization( levelsetFunctionsParametersAndDeformation );
          }
          aol::MultiVector<RealType> parameters( 1<<numberOfLevelsetFunctions, 1 );
          levelsetFunctionsParametersAndDeformation.appendReference( parameters );
          aol::MultiVector<RealType> psi( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
          levelsetFunctionsParametersAndDeformation.appendReference( psi );

          aol::MultiVector<RealType> mtmp( levelsetFunctionsParametersAndDeformation, aol::STRUCT_COPY );
    
          gradientDescent_solver.apply( levelsetFunctionsParametersAndDeformation, mtmp );
          levelsetFunctionsParametersAndDeformation = mtmp;

        }
        break;
      case SEGMENT_GRAY_VALUES_MULTI:
        {
          const int numberOfLevelsetFunctions = 2;

          const int imageDimension = 1;
          aol::MultiVector<RealType> u0MVec( 0, 0 );
          for ( int i = 0; i < imageDimension; i++ )
            u0MVec.appendReference( _u0 );
/*
          const int imageDimension = 2;
          aol::MultiVector<RealType> u0MVec( imageDimension, _grid.getNumberOfNodes() );
          generator.generateRectangular ( u0MVec[0] , 1.0, 0.25, 0.25, 1.0, 1.0 );
          generator.generateRectangular ( u0MVec[1] , 1.0, 0., 0., 0.75, 0.75 );
*/
          aol::MultiVector<RealType> meanGrayValues( 1<<numberOfLevelsetFunctions, imageDimension );

          aol::ClassicalChanVeseVectorEnergyMulti<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions, imageDimension>
            tempE( _grid, heavisideFunction, u0MVec );
/*                                      
          aol::ClassicalChanVeseEnergyMulti<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
            tempE( _grid, heavisideFunction, _u0 );
*/
          aol::HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>
            Ereg( _grid, heavisideFunction, _params.epsilonMCM );
          aol::VariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, numberOfLevelsetFunctions>
            DEreg( _grid, _params.epsilonMCM, _params.gamma );

          aol::ChanVeseEnergyLevelsetPart<RealType> Efidelity( tempE, meanGrayValues );
          aol::ClassicalChanVeseVectorEnergyMultiPhiVariation<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions, imageDimension>
            DEfidelity( _grid, heavisideFunction, u0MVec, meanGrayValues );
/*
          aol::ClassicalChanVeseEnergyMultiPhiVariation<ConfiguratorType, HeavisideFunctionType>
            DEfidelity( _grid, heavisideFunction, _u0, meanGrayValues );
*/

          aol::LinCombOp<aol::MultiVector<RealType>, aol::Scalar<RealType> > E;
          E.appendReference( Efidelity );
          E.appendReference( Ereg, _params.gamma );
          aol::LinCombOp<aol::MultiVector<RealType> > DE;
          DE.appendReference( DEfidelity );
          DE.appendReference( DEreg );

          aol::MultiVector<RealType> levelsetFunctions( numberOfLevelsetFunctions, _grid.getNumberOfNodes() );
          generator.generateChanVeseInitialization( levelsetFunctions );
          //generator.generateLineLevelset( levelsetFunctions[0], 0, 0.4 );
          //generator.generateLineLevelset( levelsetFunctions[1], 1, 0.6 );
          //generator.generateCircleLevelset ( levelsetFunctions[0], 0.4, 0.3, 0.3 );
          //generator.generateCircleLevelset ( levelsetFunctions[1], 0.4, 0.7, 0.7 );

          // ClassicalChanVeseGradientDescent solves for the levelset functions and
          // for the mean gray values (see class descriptions).
          qc::ClassicalChanVeseGradientDescent<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions, imageDimension>
            gradientDescent_solver( _grid, E, DE, heavisideFunction, u0MVec, meanGrayValues, _params.numSteps );
          gradientDescent_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "phi", this->_saveDirectory );
          //gradientDescent_solver.setTauMax( 1.0 );
          aol::MultiVector<RealType> mtmp( levelsetFunctions, aol::STRUCT_COPY );

          gradientDescent_solver.apply( levelsetFunctions, mtmp );
          levelsetFunctions = mtmp;
        }
        break;
      case GET_DEFORMATION:
        {
/*
          {
            aol::MultiVector<RealType> def1( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
            aol::MultiVector<RealType> def2( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
            generator.generateSkewSymmetricDeformation( -0.3, def1 );
            //generator.generateShearDeformation( 0.3, def1, 0.5 );
            generator.generateSkewSymmetricDeformation( 0.3, def2 );
            this->setSaveName( "def1" );
            this->plotVectorField ( def1, _grid, aol::GNUPLOT_PNG );
            this->setSaveName( "def2" );
            this->plotVectorField ( def2, _grid, aol::GNUPLOT_PNG );
            aol::MultiVector<RealType> dest( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
            qc::concatenateTwoDeformations<ConfiguratorType> ( _grid, def1, def2, dest );
            this->setSaveName( "conc" );
            this->plotVectorField ( dest, _grid, aol::GNUPLOT_PNG );
            this->plotVectorField ( dest, _grid, aol::GNUPLOT_PS );
            this->plotVectorField ( dest, _grid, aol::GNUPLOT_EPS );
          }
*/
          // Generate deformation psi and calculate u0 \circ psi
          aol::MultiVector<RealType> psi( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
          //generator.generateNonLinearDeformation( 0.3, psi );
          //generator.generateShearDeformation( 0.3, psi, 0.5 );
          //generator.generateSkewSymmetricDeformation( 0.3, psi );
          //generator.generateSymmetricDeformation( 0.3, psi, 0., 0. );
          this->setSaveName( "def" );
          this->plotVectorField ( psi, _grid, aol::GNUPLOT_PNG );
          this->plotVectorField ( psi, _grid, aol::GNUPLOT_PS );
          this->plotVectorField ( psi, _grid, aol::GNUPLOT_EPS );
          this->setSaveName( "col" );
          this->writeColorField( psi, _grid );

          ArrayType u0Deformed( _grid );
          qc::TransformFunction<RealType, ConfiguratorType::Dim> transform( _grid );
          transform.setDeformation( psi );
          transform.apply( _u0, u0Deformed );

          qc::SkewSymmetricMean<ConfiguratorType> skewSymmetricMean01( _grid, 0, 1 );
          qc::SymmetricProjector<ConfiguratorType> symmetricProjector( _grid );

          psi.setZero();

          this->setSaveName( "deformed" );
          this->savePNG( u0Deformed, _grid );
/*
          {
            qc::ScalarArray<RealType, qc::QC_2D> defX( "def_found_000.bz2" );
            qc::ScalarArray<RealType, qc::QC_2D> defY( "def_found_001.bz2");
            aol::MultiVector<RealType> deformation( 0, 0 );
            deformation.appendReference( defX );
            deformation.appendReference( defY );

            aol::MultiVector<RealType> rotation( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
            generator.generateSkewSymmetricDeformation( -0.147087, rotation );
            this->setSaveName( "all_onlydef" );
            this->saveDeformedImage( u0Deformed, _grid, deformation );
            this->setSaveName( "all" );
            this->saveDeformedImage( u0Deformed, _grid, deformation, rotation );
            aol::MultiVector<RealType> dest( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
            qc::concatenateTwoDeformations<ConfiguratorType> ( _grid, rotation, deformation, dest );
            this->setSaveName( "all_conc" );
            this->saveDeformedImage( u0Deformed, _grid, dest );
            system( "PAUSE" );
          }
*/
          CompleteAlphaPsiElasticEnergy<ConfiguratorType, HeavisideFunctionType>
            E( _grid, u0Deformed, _params.atomDistance, _params.lambda, _params.thresholdUpper, heavisideFunction );
          VariationOfCompleteAlphaPsiElasticEnergy<ConfiguratorType, HeavisideFunctionType>
            DE( _grid, u0Deformed, _params.atomDistance, _params.lambda, _params.thresholdUpper, heavisideFunction );

          aol::MultiVector<RealType> psiAndAlpha( 0, 0 );
          psiAndAlpha.appendReference( psi );
          psiAndAlpha.appendReference( alpha );

/*
          CompleteElasticEnergy<ConfiguratorType, HeavisideFunctionType>
            E( _grid, u0Deformed, 0., _params.atomDistance, _params.lambda, _params.thresholdUpper, heavisideFunction );
          VariationOfCompleteElasticEnergy<ConfiguratorType, HeavisideFunctionType>
            DE( _grid, u0Deformed, 0., _params.atomDistance, _params.lambda, _params.thresholdUpper, heavisideFunction );
*/
          SymmetricTransformationGradientDescent<ConfiguratorType, HeavisideFunctionType>
            gradientDescent_solver( _grid, E, DE, heavisideFunction, u0Deformed, 1000, 1., 0.0000001 );
          gradientDescent_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "deform", this->_saveDirectory );
          aol::EnergyDescentPlotter<SymmetricTransformationGradientDescent<ConfiguratorType, HeavisideFunctionType> >
            energyDescentPlotter( "energy.txt", "sigma.txt", this->_saveDirectory, gradientDescent_solver );

          aol::MultiVector<RealType> mtmp( psiAndAlpha, aol::STRUCT_COPY );
          gradientDescent_solver.apply( psiAndAlpha, mtmp );
/*
          aol::MultiVector<RealType> mtmp( psi, aol::STRUCT_COPY );
          gradientDescent_solver.apply( psi, mtmp );
*/
          energyDescentPlotter.plot( "EnergyPlot" );

          this->setSaveName( "deform" );
          this->setSaveDirectory( this->_saveDirectory );
          this->saveDeformedImage( u0Deformed, _grid, mtmp );
          aol::MultiVector<RealType> rotation( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
          generator.generateSkewSymmetricDeformation( mtmp[ConfiguratorType::Dim][0], rotation );
          this->setSaveName( "deform_rot" );
          this->saveDeformedImage( u0Deformed, _grid, mtmp, rotation );


          this->setSaveName( "def_found" );
          this->plotVectorField ( mtmp, _grid, aol::GNUPLOT_PNG );
          this->plotVectorField ( mtmp, _grid, aol::GNUPLOT_PS );
          this->plotVectorField ( mtmp, _grid, aol::GNUPLOT_EPS );
          this->setSaveName( "col_found" );
          this->writeColorField( mtmp, _grid );

          this->setSaveName( "def_found" );
          this->saveTimestepBZ2 ( 0, mtmp[0], _grid );
          this->setSaveName( "def_found" );
          this->saveTimestepBZ2 ( 1, mtmp[1], _grid );

          aol::MultiVector<RealType> dest( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
          qc::concatenateTwoDeformations<ConfiguratorType> ( _grid, mtmp, rotation, dest );
          this->setSaveName( "conc_res" );
          this->plotVectorField ( dest, _grid, aol::GNUPLOT_PNG );
          this->plotVectorField ( dest, _grid, aol::GNUPLOT_PS );
          this->plotVectorField ( dest, _grid, aol::GNUPLOT_EPS );

          //qc::SymmetricLengthEnergy<ConfiguratorType> symmetricLengthEnergy( _grid );
          //qc::VariationOfSymmetricLengthEnergy<ConfiguratorType> variationOfSymmetricLengthEnergy( _grid );
          //aol::testFirstDerivativeMultiVectorAllDirections<ConfiguratorType>( _grid, psi, symmetricLengthEnergy, variationOfSymmetricLengthEnergy, "t", 0.01 );
        }
        break;
      case SEGMENT_GRAINS_WITH_DEFORMATION:
        {
          aol::MultiVector<RealType> psi1( ConfiguratorType::Dim, _grid.getNumberOfNodes() );
          aol::MultiVector<RealType> psi2( ConfiguratorType::Dim, _grid.getNumberOfNodes() );

          generator.generateSkewSymmetricDeformation( -0.0764884, psi1 );
          generator.generateSkewSymmetricDeformation( -0.34 , psi2 );

          ArrayType psiLevelsetFunction( _grid );
          generator.generateLineLevelset( psiLevelsetFunction, 0, 0.4 );
          //generator.generateCircleLevelset( psiLevelsetFunction, 0.4 );
          
          aol::MultiVector<RealType> psisAndLevelset( 0, 0 );
          psisAndLevelset.appendReference( psi1 );
          psisAndLevelset.appendReference( psi2 );
          psisAndLevelset.appendReference( psiLevelsetFunction );
          
          aol::MultiVector<RealType> mtmp( psisAndLevelset, aol::STRUCT_COPY );

          CompleteElasticEnergyWithPsiSegmentation<ConfiguratorType, HeavisideFunctionType>
            E( _grid,_u0, 0., _params.atomDistance, _params.gamma, _params.lambda, _params.thresholdUpper, heavisideFunction );
          VariationOfCompleteElasticEnergyWithPsiSegmentation<ConfiguratorType, HeavisideFunctionType>
            DE( _grid,_u0, 0., _params.atomDistance, _params.gamma, _params.lambda, _params.thresholdUpper, heavisideFunction );

          //aol::testFirstDerivativeMultiVectorAllDirections<ConfiguratorType>( _grid, psisAndLevelset, E, DE, "t", 0.01 );
          //aol::testFirstDerivativeMultiVectorSingleComponentAllDirections<ConfiguratorType>( _grid, psisAndLevelset, E, DE, 4, "t", 0.001 );

          SymmetricTransformationGradientDescent<ConfiguratorType, HeavisideFunctionType>
            gradientDescent_solver( _grid, E, DE, heavisideFunction, _u0, 1000, 1., 0.0000001 );
          gradientDescent_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "deform", this->_saveDirectory );

          gradientDescent_solver.apply( psisAndLevelset, mtmp );
          psisAndLevelset = mtmp;

          this->setSaveName( "phi" );
          this->saveTimestepPNGDrawRedAndBlueIsoline( -1, -1, 0., 0.1, psiLevelsetFunction, _u0, _grid );
          this->setSaveName( "def1_found" );
          this->plotVectorField ( psi1, _grid, aol::GNUPLOT_PNG );
          this->plotVectorField ( psi1, _grid, aol::GNUPLOT_PS );
          this->plotVectorField ( psi1, _grid, aol::GNUPLOT_EPS );
          this->setSaveName( "def2_found" );
          this->plotVectorField ( psi2, _grid, aol::GNUPLOT_PNG );
          this->plotVectorField ( psi2, _grid, aol::GNUPLOT_PS );
          this->plotVectorField ( psi2, _grid, aol::GNUPLOT_EPS );
        }
        break;
      case SEGMENT_GRAINS_WITH_ALPHA_AND_GET_DEFORMATION:
        {
          aol::MultiVector<RealType> psi( ConfiguratorType::Dim, _grid.getNumberOfNodes() );

          generator.generateNonLinearDeformation( 0.025, psi );
          //generator.generateShearDeformation( 0.3, psi, 0.5 );
          //generator.generateSkewSymmetricDeformation( 0.3, psi );
          //generator.generateSymmetricDeformation( 0.1, psi, 0.5, 0.5 );
          //psi += 0.1;
          ArrayType u0Deformed( _grid );
          qc::TransformFunction<RealType, ConfiguratorType::Dim> transform( _grid );
          transform.setDeformation( psi );
          transform.apply( _u0, u0Deformed );
          _u0 = u0Deformed;
          //psi += -0.1;
          psi.setZero();

          ArrayType psiLevelsetFunction( _grid );
          generator.generateLineLevelset( psiLevelsetFunction, 0, 0.5 );
/*
          psiLevelsetFunction.load( "phi_found.bz2" );
          ArrayType psi0Array( psi[0], _grid );
          psi0Array.load( "def_found_000.bz2" );
          ArrayType psi1Array( psi[1], _grid );
          psi1Array.load( "def_found_001.bz2" );
          //generator.generateCircleLevelset( psiLevelsetFunction, 0.4 );
          this->setSaveName( "def" );
          ofstream vecoff( "vec.dat" );
          qc::WriteVectorFieldAsGnuplotFile<RealType>( vecoff, psi0Array, psi1Array, 0.1 );
*/

          aol::MultiVector<RealType> psiAlphaAndLevelset( 0, 0 );
          psiAlphaAndLevelset.appendReference( psi );
          psiAlphaAndLevelset.appendReference( alpha );
          psiAlphaAndLevelset.appendReference( psiLevelsetFunction );

          aol::MultiVector<RealType> mtmp( psiAlphaAndLevelset, aol::STRUCT_COPY );

          CompleteElasticEnergyWithAlphaSegmentation<ConfiguratorType, HeavisideFunctionType>
            E( _grid,_u0, _params.atomDistance, _params.gamma, _params.lambda, _params.thresholdUpper, heavisideFunction );
          VariationOfCompleteElasticEnergyWithAlphaSegmentation<ConfiguratorType, HeavisideFunctionType>
            DE( _grid,_u0, _params.atomDistance, _params.gamma, _params.lambda, _params.thresholdUpper, heavisideFunction );

          //aol::testFirstDerivativeMultiVectorSingleComponentAllDirections<ConfiguratorType>( _grid, psiAlphaAndLevelset, E, DE, 2, "t", 0.001 );

          SymmetricTransformationGradientDescent<ConfiguratorType, HeavisideFunctionType>
            gradientDescent_solver( _grid, E, DE, heavisideFunction, _u0, 1000, 1., 0.0000001 );
          gradientDescent_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "deform", this->_saveDirectory );
          aol::EnergyDescentPlotter<SymmetricTransformationGradientDescent<ConfiguratorType, HeavisideFunctionType> >
            energyDescentPlotter( "energy.txt", "sigma.txt", this->_saveDirectory, gradientDescent_solver );

          gradientDescent_solver.apply( psiAlphaAndLevelset, mtmp );
          psiAlphaAndLevelset = mtmp;

          energyDescentPlotter.plot( "EnergyPlot" );

          this->setSaveName( "phi_found" );
          this->saveTimestepPNGDrawRedAndBlueIsoline( -1, -1, 0., 0.1, psiLevelsetFunction, _u0, _grid );
          this->setSaveName( "phi_found" );
          this->saveTimestepBZ2( -1, -1, psiLevelsetFunction, _grid );
          this->setSaveName( "col_found" );
          this->writeColorField( psi, _grid );
          this->setSaveName( "def_found" );
          this->plotVectorField ( psi, _grid, aol::GNUPLOT_PNG );
          this->plotVectorField ( psi, _grid, aol::GNUPLOT_PS );
          this->plotVectorField ( psi, _grid, aol::GNUPLOT_EPS );
          this->setSaveName( "def_found" );
          this->saveTimestepBZ2 ( -1, psi, _grid );
        }
        break;
        default:
          throw aol::UnimplementedCodeException ( "Unsupported mode", __FILE__, __LINE__ );
      }

/*
      SimpleGradientDescent<ConfiguratorType, aol::Vector<RealType> > simple(_grid, DE, _params.numSteps );
      simple.apply( levelsetFunction, tmp );
      levelsetFunction = tmp;

      this->setSaveName( "phi" );
      saveTimestepPNGDrawRedIsoline( -1, 0., levelsetFunction, _u0, _grid);
*/
    }
  }   // end computeTimesteps
};


#endif
