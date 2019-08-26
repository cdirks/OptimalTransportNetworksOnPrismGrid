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

#ifndef __CHANVESEDESCENT_H
#define __CHANVESEDESCENT_H

#include <signedDistanceOp.h>
#include <ChanVese.h>
#include <gradientDescent.h>
#include <quocTimestepSaver.h>

namespace qc {

/*
 * Gradient descent to minimize energies of type ChanVeseEnergyInterface
 *
 * Needs SignedDistanceOp, therefore part of module eikonal.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class ChanVeseGradientDescent
: public aol::GradientDescentComponentWiseTimestepControlled<ConfiguratorType>,
  public qc::QuocTimestepSaver<ConfiguratorType> {
public:
  enum MODE {
    FULL_MINIMIZATION,
    MINIMIZE_WRT_LEVELSETFUNCTIONS
  };
  /*
   * With these flags the behavior of this class can be configured.
   */
  enum FLAGS {
    //! Never do redistancing in postProcess.
    NO_REDISTANCING = 1,
    //! Don't apply LinearSmoothOp in smoothDirection.
    NO_SMOOTHING = 2,
    //! If the zero levelset does not intersect with [0,1]^d, don' push it in postProcess.
    DO_NOT_FORCE_ZERO_LEVELLINE_INSIDE_DOMAIN = 4,
    //! Don't save the timesteps with saveTimestepPNGDrawRedAndBlueIsoline.
    DO_NOT_SAVE_TIMESTEPS_AS_ISOLINE_PNG = 8,
    //! Clamp the levelset functions to [0,1], only meant to be used in combination with aol::ClampedIdentityHeavisideFunction.
    CLAMP_TO_ZERO_ONE = 16,
    //! Save the timesteps with savePNG (only done if DO_NOT_SAVE_TIMESTEPS_AS_ISOLINE_PNG is also set).
    SAVE_LEVELSET_FUNCTIONS_AS_PNG = 32,
    //! Instead of aol::HeavisideFunctionLumpedMassOp use aol::LumpedMassOp as metric.
    USE_STANDARD_LUMPED_MASS_MATRIX_AS_METRIC = 64
  };
private:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const aol::Vector<RealType> &_image;
  const MODE _mode;
  unsigned int _flags;
  RealType _displayedIsovalue;
  const aol::LumpedMassOp<ConfiguratorType> _lumpedMassInv;
public:
  ChanVeseGradientDescent( const typename ConfiguratorType::InitType &Initializer,
                           const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                           const aol::Op<aol::MultiVector<RealType> > &DE,
                           const HeavisideFunctionType &HeavisideFunction,
                           const aol::Vector<RealType> &Image,
                           const int MaxIterations = 50,
                           const RealType StartTau = aol::ZOTrait<RealType>::one,
                           const RealType StopEpsilon = aol::ZOTrait<RealType>::zero,
                           const MODE Mode = FULL_MINIMIZATION )
  : aol::GradientDescentComponentWiseTimestepControlled<ConfiguratorType>(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon),
   _heavisideFunction( HeavisideFunction ),
   _image(Image),
   _mode(Mode),
   _flags(0),
   _displayedIsovalue(aol::NumberTrait<RealType>::zero),
   _lumpedMassInv ( Initializer, aol::INVERT )
  {
    this->setFilterWidth( 1. );
    this->setUpdateFilterWidth( true );
    this->setMaximumNumIterationsPerFilterWidth( 20 );
    if( _mode == FULL_MINIMIZATION ){
      this->_customTauBlocks.resize(NumberOfLevelsetFunctions+1);
      this->_customTauBlocks.setAll( 1 );
      this->_customTauBlocks[NumberOfLevelsetFunctions] = 1<<NumberOfLevelsetFunctions;
    }
  }
  void setFlags ( const unsigned int Flags ){
    _flags = Flags;
    // NO_SMOOTHING is not compatible with automatic filter width control.
    // As long as the flags are not moved to GradientDescentBase we have
    // to make this workaround.
    if( !(_flags & NO_SMOOTHING) )
      this->setUpdateFilterWidth( true );
    else
      this->setUpdateFilterWidth( false );
  }

  void setDisplayedIsovalue ( const RealType DisplayedIsovalue ) {
    _displayedIsovalue = DisplayedIsovalue;
  }

  RealType getDisplayedIsovalue ( ) const {
    return _displayedIsovalue;
  }

  void configureForEsedogluModel () {
    setFlags ( _flags | NO_REDISTANCING | NO_SMOOTHING | DO_NOT_FORCE_ZERO_LEVELLINE_INSIDE_DOMAIN | CLAMP_TO_ZERO_ONE | USE_STANDARD_LUMPED_MASS_MATRIX_AS_METRIC );
    setDisplayedIsovalue ( 0.5 );
    this->setTimestepController( this->GRADIENT_SIMPLE );
  }

  void configureForQuadraticEsedogluModel () {
    setFlags ( _flags | NO_REDISTANCING | NO_SMOOTHING | DO_NOT_FORCE_ZERO_LEVELLINE_INSIDE_DOMAIN | SAVE_LEVELSET_FUNCTIONS_AS_PNG | USE_STANDARD_LUMPED_MASS_MATRIX_AS_METRIC );
    setDisplayedIsovalue ( 0.5 );
  }

  virtual bool postProcess( aol::MultiVector<RealType> &Position, const int Iteration ) const {
    bool postProcessed = false;
    // If the zero levelset does not intersect with [0,1]^d, push it back so that it does.
    if( !(_flags & DO_NOT_FORCE_ZERO_LEVELLINE_INSIDE_DOMAIN) ){
      for ( int i = 0; i < NumberOfLevelsetFunctions; i++ ){
        const RealType minVal = Position[i].getMinValue();
        const RealType maxVal = Position[i].getMaxValue();
        if ( minVal > 0. ){
          Position[i].addToAll( - (minVal + 0.1*(maxVal-minVal)) );
          postProcessed = true;
        }
        else if ( maxVal < 0. ){
          Position[i].addToAll( - (maxVal - 0.1*(maxVal-minVal)) );
          postProcessed = true;
        }
      }
    }
    // Redistance all levelset functions every five iterations.
    if( !(_flags & NO_REDISTANCING) && (Iteration % 5 == 0) ){
      typename qc::SignedDistanceOpTrait<ConfiguratorType>::OpType signedDistOp( this->_grid );
      for ( int i = 0; i < NumberOfLevelsetFunctions; i++ ){
        signedDistOp.apply( Position[i], Position[i] );
      }
      postProcessed = true;
    }
    // Set values in the levelset functions that are smaller than zero to zero and those bigger than one to one.
    if( (_flags & CLAMP_TO_ZERO_ONE) ){
      for ( int i = 0; i < NumberOfLevelsetFunctions; i++ ){
        Position[i].clamp ( aol::NumberTrait<RealType>::zero, aol::NumberTrait<RealType>::one );
      }
      // postProcessed is only used for the stopping criterion to
      // prevent energy increases caused by the post processing
      // from stopping the iteration. Since clamping is only
      // meant to be used with aol::ClampedIdentityHeavisideFunction
      // we don't need to set postProcessed to true here.
    }

    return postProcessed;
  }

  virtual void smoothDirection( const aol::MultiVector<RealType> &Position, const RealType Sigma, aol::MultiVector<RealType> &Direction ) const {
    aol::Vector<RealType> tmp( this->_grid.getNumberOfNodes() );

    qc::LinearSmoothOp<RealType> linSmooth;
    if( !(_flags & NO_SMOOTHING) ) {
      linSmooth.setCurrentGrid( this->_grid );
      linSmooth.setSigma( Sigma );
    }

    for ( int i = 0; i < NumberOfLevelsetFunctions; i++ ){
      if ( _flags & USE_STANDARD_LUMPED_MASS_MATRIX_AS_METRIC ) {
        _lumpedMassInv.apply( Direction[i], tmp );
      }
      else {
        aol::HeavisideFunctionLumpedMassOp<ConfiguratorType, HeavisideFunctionType >
          heavisideFunctionLumpedMassInvOpAlphaLevelset( this->_grid, Position[i], _heavisideFunction, aol::INVERT );
        heavisideFunctionLumpedMassInvOpAlphaLevelset.apply( Direction[i], tmp );
      }
      Direction[i] = tmp;
      if( !(_flags & NO_SMOOTHING) ){
        linSmooth.apply( Direction[i], Direction[i] );
      }
    }
  }
  virtual void writeTimeStep( const aol::MultiVector<RealType> &MDest, const int Iterations ) const {
    for ( int i = 0; i < NumberOfLevelsetFunctions; i++ ){
      this->saveTimestepBZ2( i, Iterations, MDest[i], this->_grid );
      if( !(_flags & DO_NOT_SAVE_TIMESTEPS_AS_ISOLINE_PNG) )
        this->saveTimestepPNGDrawRedAndBlueIsoline( i, Iterations, _displayedIsovalue, _displayedIsovalue + 0.1, MDest[i], _image, this->_grid );
      else if ( (_flags & SAVE_LEVELSET_FUNCTIONS_AS_PNG) )
       this->savePNG( MDest[i], this->_grid, Iterations, i );
    }

    // The following visualization is mainly intended for 2D and handles 3D data as 2D slices.
    if ( this->_grid.getDimOfWorld() == 2 ) {
      // Generate and write an image showing the current segmentation color coded.
      qc::MultiArray<RealType, 2, 3> vImageArray ( this->_grid );
      qc::DataGenerator<ConfiguratorType> generator( this->_grid );
      generator.generateColoredSegmentation ( NumberOfLevelsetFunctions, MDest, _image, vImageArray, _displayedIsovalue );
      this->savePNG( vImageArray, Iterations, NumberOfLevelsetFunctions );
    }
    // Assume 3D:
    else {
      SliceExtractor<ConfiguratorType> _sliceExtractor( this->_grid, MDest );
      qc::ScalarArray<RealType, qc::QC_2D> whiteImage ( _sliceExtractor.getNumX(), _sliceExtractor.getNumY() );
      whiteImage.setAll( 1. );
      aol::MultiVector<RealType> levelsetFunctionsSlice ( NumberOfLevelsetFunctions, _sliceExtractor.getNumX()*_sliceExtractor.getNumY() );

      for ( int z = 0; z < _sliceExtractor.getNumOfSlices(); z++ ) {
        _sliceExtractor.getSlice ( z, levelsetFunctionsSlice );
        qc::MultiArray<RealType, 2, 3> vImageArray ( _sliceExtractor.get2DGridReference() );
        qc::DataGenerator<typename SliceExtractor<ConfiguratorType>::ConfiguratorType2D> generator( _sliceExtractor.get2DGridReference() );
        generator.generateColoredSegmentation ( NumberOfLevelsetFunctions, levelsetFunctionsSlice, whiteImage, vImageArray, _displayedIsovalue );
        this->savePNG( vImageArray, Iterations, NumberOfLevelsetFunctions, z );
      }
    }

    if( _mode == FULL_MINIMIZATION ){
      this->writeMultiVectorToTextFile ( MDest, "parameters.txt", Iterations, this->_maxIterations, NumberOfLevelsetFunctions, NumberOfLevelsetFunctions+(1<<NumberOfLevelsetFunctions));
    }
  }
};

/**
 * Gradient descent to minimize the classical Chan Vese energy for
 * multi levelset gray value segmentation. The minimazation with
 * respect to the gray values is done explicitly in the function
 * postProcess. We don't need to do a gradient descent for this.
 *
 * The gray values must not be included in the apply arguments,
 * but a reference to the vector containing them has to be supplied
 * in the constructor. postProcess alters the values behind the
 * reference. Therefore the energy and the variation must use the
 * same gray value reference. Otherwise they won't get the updated
 * values calculated by postProcess.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int ImageDimension>
class ClassicalChanVeseGradientDescent
: public qc::ChanVeseGradientDescent<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> {
private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::ChanVeseGradientDescent<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> SuperClass;
  aol::MultiVector<RealType> &_meanGrayValues;
  const aol::MultiVector<RealType> &_imageMVec;
public:
  ClassicalChanVeseGradientDescent( const typename ConfiguratorType::InitType &Initializer,
                                    const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                    const aol::Op<aol::MultiVector<RealType> > &DE,
                                    const HeavisideFunctionType &HeavisideFunction,
                                    const aol::MultiVector<RealType> &ImageMVec,
                                    aol::MultiVector<RealType> &MeanGrayValues,
                                    const int MaxIterations = 50,
                                    const RealType StartTau = aol::ZOTrait<RealType>::one,
                                    const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
    : SuperClass( Initializer, E, DE, HeavisideFunction, ImageMVec[0], MaxIterations, StartTau, StopEpsilon, SuperClass::MINIMIZE_WRT_LEVELSETFUNCTIONS ),
      _meanGrayValues( MeanGrayValues ),
      _imageMVec( ImageMVec )
  {
  }
  virtual bool postProcess( aol::MultiVector<RealType> &Position, const int Iteration ) const {
    bool postProcessed = SuperClass::postProcess( Position, Iteration );

    // Calculate the new gray values based on the current levelset functions (stored in Position).
    aol::MultiLevelsetVolumes<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions>
      volE( this->_grid, Position, this->_heavisideFunction );

    aol::MultiLevelsetWeightedVolumes<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ImageDimension>
      volEweight( this->_grid, Position, this->_heavisideFunction, _imageMVec );

    aol::MultiVector<RealType> tempMVec( _meanGrayValues.numComponents(), 1 );
    volE.apply( tempMVec, tempMVec );
    // MultiLevelsetWeightedVolumes only uses the structure of the first argument, not the actual content.
    volEweight.apply( _meanGrayValues, _meanGrayValues );
    for( int i = 0; i < _meanGrayValues.numComponents(); i++ )
      _meanGrayValues[i] /= tempMVec[i][0];
    this->writeMultiVectorToTextFile ( _meanGrayValues, "grayvals.txt", Iteration, this->_maxIterations );

    if ( ( ImageDimension == 1 ) && ( ConfiguratorType::Dim == 2 ) ) {
      this->savePiecewiseConstantImageSegmented ( _meanGrayValues, Position, this->_grid, this->_displayedIsovalue, Iteration, Position.numComponents()+1 );
    }
    else if ( ( ImageDimension == 3 ) && ( ConfiguratorType::Dim == 2 ) ) {
      this->savePiecewiseConstantColorImageSegmented ( _meanGrayValues, Position, this->_grid, this->_displayedIsovalue, Iteration, Position.numComponents()+1 );
    }
    else
    {
      std::vector<aol::PlotOutFileType> outTypeVec;
      outTypeVec.push_back(aol::GNUPLOT_PNG);
      outTypeVec.push_back(aol::GNUPLOT_PS);
      outTypeVec.push_back(aol::GNUPLOT_EPS);
      this->plotVectorFieldSegmented ( _imageMVec, Position, this->_grid, outTypeVec, Position.numComponents()+1, Iteration, true, this->_displayedIsovalue );
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

} // end namespace qc

#endif //__CHANVESEDESCENT_H
