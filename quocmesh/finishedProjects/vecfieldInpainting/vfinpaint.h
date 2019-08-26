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

#ifndef __VFINPAINT_H
#define __VFINPAINT_H

#include <ArmijoSearch.h>
#include <gradientflow.h>
#include <quocTimestepSaver.h>
#include <chanVeseDescent.h>
#include <narrow.h>
#include <narrowBandConfigurators.h>
#include <anisoStiffOps.h>

template <typename RealType>
class DiscreteVecInpaintParams{
  string dirName;
public:
  //! epsilon is the regularization parameter of the absolute value
  const RealType epsilon;
  //!
  const RealType mu;
  //!
  const RealType nu;
  //!
  const RealType lambda;
  //! delta is the regularization parameter of the Heaviside function
  const RealType delta;
  //! theta is the threshold for the mask, i.e. the value of the levelline that is used to mark the destroyed regions
  const RealType theta;
  //! numSteps is the maximal number of allowed steps in the gradient descent.
  const int numSteps;
  DiscreteVecInpaintParams( aol::ParameterParser &Parser )
  : epsilon( Parser.getDouble( "epsilon" ) ),
    mu( Parser.getDouble( "mu" ) ),
    nu( Parser.getDouble( "nu" ) ),
    lambda( Parser.getDouble( "lambda" ) ),
    delta( Parser.getDouble( "delta" ) ),
    theta( Parser.getDouble( "theta" ) ),
    numSteps( Parser.getInt( "numSteps" ) )
  {
    char fn[1024];
    sprintf( fn, "%.2f_%.2f_%.2f_%.2f_%.2f_%d", epsilon, mu, nu, lambda, delta, numSteps);
    dirName = fn;
  }
  const char* getDirName() const{
    return dirName.c_str();
  }
};

namespace qc {

}

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename PeronaMalikFunctionType, typename HeavisideFunctionType>
class VFieldInpaintingEnergy : public aol::Op< aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_image;
  const aol::MultiVector<RealType> &_v0;
  const PeronaMalikFunctionType &_peronaMalikFunction;
  const DiscreteVecInpaintParams<RealType> &_params;
  const aol::DiagonalBlockOp<RealType> &_heavisideFunctionWeightedMassMatBlockOp;
  const aol::Vector<RealType> &_levelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
  const bool _useFidelityTerm;
  const bool _weightRegTermWithOneMinusH;
public:
  VFieldInpaintingEnergy ( const typename ConfiguratorType::InitType &Initializer,
                           const aol::Vector<RealType> &Image,
                           const aol::MultiVector<RealType> &V0,
                           const PeronaMalikFunctionType &PeronaMalikFunction,
                           const DiscreteVecInpaintParams<RealType> &Params,
                           const aol::DiagonalBlockOp<RealType> &HeavisideFunctionWeightedMassMatBlockOp,
                           const aol::Vector<RealType> &LevelsetFunction,
                           const HeavisideFunctionType &HeavisideFunction,
                           const bool UseFidelityTerm,
                           const bool WeightRegTermWithOneMinusH )
    : _grid ( Initializer ),
      _image ( Image ),
      _v0 ( V0 ),
      _peronaMalikFunction ( PeronaMalikFunction ),
      _params ( Params ),
      _heavisideFunctionWeightedMassMatBlockOp ( HeavisideFunctionWeightedMassMatBlockOp ),
      _levelsetFunction ( LevelsetFunction ),
      _heavisideFunction ( HeavisideFunction ),
      _useFidelityTerm ( UseFidelityTerm ),
      _weightRegTermWithOneMinusH ( WeightRegTermWithOneMinusH )
  {}

  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    qc::ImageDrivenAnisoTVVecNorm2D<ConfiguratorType, PeronaMalikFunctionType> imageDrivenAnisoTVVecNorm( _grid, _image, _peronaMalikFunction,  _params.lambda, _params.epsilon );

    aol::OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, qc::ImageDrivenAnisoTVVecNorm2D<ConfiguratorType, PeronaMalikFunctionType> >
      oneMinusHImageDrivenAnisoTVVecNorm( _grid, _levelsetFunction, _heavisideFunction, imageDrivenAnisoTVVecNorm );

    if ( _weightRegTermWithOneMinusH == false )
      imageDrivenAnisoTVVecNorm.apply( MArg, Dest );
    else
      oneMinusHImageDrivenAnisoTVVecNorm.apply( MArg, Dest );
    Dest *= _params.nu;

    if ( _useFidelityTerm ) {
      aol::MultiVector<RealType> temp ( MArg, aol::DEEP_COPY );
      aol::MultiVector<RealType> temp2 ( MArg, aol::STRUCT_COPY );
      temp -= _v0;
      _heavisideFunctionWeightedMassMatBlockOp.apply( temp, temp2 );
      Dest += 0.5*(temp * temp2);
    }
  }

  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    aol::Scalar<RealType> tmp;
    apply( MArg, tmp );
    Dest += tmp;
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename PeronaMalikFunctionType, typename HeavisideFunctionType>
class VariationOfVFieldInpaintingEnergy : public aol::Op< aol::MultiVector<typename ConfiguratorType::RealType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_image;
  const aol::MultiVector<RealType> &_v0;
  const PeronaMalikFunctionType &_peronaMalikFunction;
  const DiscreteVecInpaintParams<RealType> &_params;
  const aol::DiagonalBlockOp<RealType> &_heavisideFunctionWeightedMassMatBlockOp;
  const aol::Vector<RealType> &_levelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
  const bool _useFidelityTerm;
  const bool _weightRegTermWithOneMinusH;
public:
  VariationOfVFieldInpaintingEnergy ( const typename ConfiguratorType::InitType &Initializer,
                                      const aol::Vector<RealType> &Image,
                                      const aol::MultiVector<RealType> &V0,
                                      const PeronaMalikFunctionType &PeronaMalikFunction,
                                      const DiscreteVecInpaintParams<RealType> &Params,
                                      const aol::DiagonalBlockOp<RealType> &HeavisideFunctionWeightedMassMatBlockOp,
                                      const aol::Vector<RealType> &LevelsetFunction,
                                      const HeavisideFunctionType &HeavisideFunction,
                                      const bool UseFidelityTerm,
                                      const bool WeightRegTermWithOneMinusH )
    : _grid ( Initializer ),
      _image ( Image ),
      _v0 ( V0 ),
      _peronaMalikFunction ( PeronaMalikFunction ),
      _params ( Params ),
      _heavisideFunctionWeightedMassMatBlockOp ( HeavisideFunctionWeightedMassMatBlockOp ),
      _levelsetFunction ( LevelsetFunction ),
      _heavisideFunction ( HeavisideFunction ),
      _useFidelityTerm ( UseFidelityTerm ),
      _weightRegTermWithOneMinusH ( WeightRegTermWithOneMinusH )
  {}

  virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    qc::VariationOfImageDrivenAnisoTVVecNorm2D<ConfiguratorType, PeronaMalikFunctionType> regVar( _grid, _image, _peronaMalikFunction, _params.lambda, _params.epsilon );

    aol::OneMinusHeavisideFENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, HeavisideFunctionType, qc::VariationOfImageDrivenAnisoTVVecNorm2D<ConfiguratorType, PeronaMalikFunctionType> >
      oneMinusHRegVar ( _grid, _levelsetFunction, _heavisideFunction, regVar );

    if ( _weightRegTermWithOneMinusH == false )
      regVar.apply ( MArg, MDest );
    else
      oneMinusHRegVar.apply ( MArg, MDest );
    MDest *= _params.nu;

    if ( _useFidelityTerm ) {
      aol::MultiVector<RealType> temp ( MArg, aol::DEEP_COPY );
      temp -= _v0;
      _heavisideFunctionWeightedMassMatBlockOp.applyAdd( temp, MDest );
    }
  }

  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
    aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
    apply( MArg, tmp );
    MDest += tmp;
  }
};

/**
 * Class for optical flow field inpainting using image intenfsity information.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class VectorFieldInpainting : public qc::QuocTimestepSaver<ConfiguratorType>
{
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  // ******** some more parameters ************

  const typename ConfiguratorType::InitType _grid;
  aol::ParameterParser &_parser;
  ArrayType _u0;
  aol::MultiVector<RealType> _v0;
  const DiscreteVecInpaintParams<RealType> _params;

public:
  VectorFieldInpainting( aol::ParameterParser &Parser ) :
    qc::QuocTimestepSaver<ConfiguratorType>( Parser.getInt("saveOffset"), Parser.getInt("numSaveFirst"), "d", true ),
    _grid( qc::getSizeFromArrayFile( Parser.getString( "input-image" ) ) ),
    _parser( Parser ),
    _u0( _grid ),
    _v0( _grid ),
    _params( Parser )
  {
    char fn_reference[1024];
    Parser.getString( "input-image", fn_reference );
    _u0.load( fn_reference );
    _u0 /= _u0.getMaxValue();
    char filename[1024];
    strncpy(filename, aol::getFileName(fn_reference), 1023);

    Parser.getString( "saveDirectory", fn_reference );
    char fn[1024];

    sprintf( fn, "%s_%s_%s", fn_reference, filename, _params.getDirName() );
    sprintf( fn_reference, "%s/", fn );
    this->setSaveDirectory( fn_reference );
    aol::makeDirectory( fn );
    Parser.dumpToFile( "parameter-dump.txt", fn_reference );

    this->activateSaving();
  }

  virtual ~VectorFieldInpainting() {
  }

  /*
   * Loads the node values of the velocity from disk. The parameter
   * ParameterName in the parameter file has to specify the filenames of the
   * components.
   *
   * Returns true if loading was sucessful, false otherwise.
   *
   * TODO: This function is copy and paste from CrystalSegmentation::loadLevelsetFunctions()!
   * Resolve this!
   *
   */
  bool loadVelocity( aol::MultiVector<RealType> &Velocity, const char *ParameterName ) const {
    if( _parser.hasVariable ( ParameterName ) ){
      int dimSize = _parser.getDimSize (ParameterName, 0 );
      char filename[1024];
      if ( Velocity.numComponents() != dimSize ){
        cerr << "Can't load levelsetfunctions: Velocity.numComponents() != dimSize!\n";
        return false;
      }
      for ( int i = 0; i < dimSize; i ++ ){
        _parser.getString ( ParameterName, filename, i );
        ArrayType tempArray( Velocity[i], _grid );
        tempArray.load( filename );
      }
      return true;
    }
    else
      return false;
  }

  void computeTimesteps( )
  {
    typedef qc::PeronaMalikWeightingFunction<RealType> PeronaMalikFunctionType;
    PeronaMalikFunctionType peronaMalikFunction( 1, aol::ZOTrait<RealType>::one/aol::Sqr ( _params.mu ) );

    //typedef aol::IdentityFunction<RealType> HeavisideFunctionType;
    typedef aol::ArcTanHeavisideFunction<RealType> HeavisideFunctionType;
    HeavisideFunctionType heavisideFunction( _params.delta );

    typedef typename ConfiguratorType::MatrixType MatrixType;
    MatrixType heavisideFunctionWeightedMassMat( qc::GridSize<ConfiguratorType::Dim>::createFrom ( _grid ) );

    qc::DataGenerator<ConfiguratorType> generator ( _grid );

    ArrayType levelsetFunction ( _grid );

    if( _parser.hasVariable ( "input-mask" ) ){
      levelsetFunction.load( _parser.getString ( "input-mask" ).c_str() );
      levelsetFunction.addToAll( -_params.theta );
    }
    else{
      //generator.generateCircleLevelset ( levelsetFunction, 0.05, 0.6 );
      generator.generateCircleLevelset ( levelsetFunction, 0.15 );
    }
    aol::HeavisideFunctionWeightedMassOp<ConfiguratorType, HeavisideFunctionType>
      heavisideFunctionWeightedMassOp( _grid, levelsetFunction, heavisideFunction );

    aol::MultiVector<RealType> velocity ( _grid );
    generator.generateIdentity( velocity );

    this->loadVelocity( _v0, "input-velocity" );

    for ( int i = 0; i < _v0[0].size(); i++ ){
      if ( levelsetFunction[i] <= 0. ){
        _v0[0][i] = 0.;
        _v0[1][i] = 0.;
      }
    }
    if ( !loadVelocity( velocity, "initial-velocity" ) )
      velocity = _v0;

    this->setSaveName( "v_orig_col" );
    this->writeColorField( _v0, _grid );
    this->setSaveName( "v_orig" );
    this->plotVectorField ( _v0, _grid, aol::GNUPLOT_PNG );
    this->setSaveName( "u" );
    this->saveTimestepPNGDrawRedIsoline ( -1, -1, 0., levelsetFunction, _u0, _grid );

    // Save the original (i.e. unextended) mask. It will be necessary for the proper
    // choice of Dirichlet nodes if we use "true inpainting".
    ArrayType originalMask ( _grid );
    originalMask = levelsetFunction;
    originalMask.threshold ( 0, 0, 1 );

    // Now make the mask marking the destroyed part of the levelset bigger.
    // Either using dilation
    if ( _parser.hasVariable ( "useDilationToIncreaseMaskSize" ) ) {
      qc::BitArray<qc::QC_2D> mask ( qc::GridSize<qc::QC_2D>::createFrom ( _grid ) );
      for ( int i = 0; i < mask.size(); ++i ) {
        if ( levelsetFunction[i] <= 0 )
          mask.set( i, true );
      }
      string filename = aol::strprintf ( "%smask0.pgm", this->getSaveDirectory() );
      mask.save ( filename.c_str() );
      int dilateMaskBy = _parser.getInt( "dilateMaskBy" );
      for ( int i = 0; i < dilateMaskBy; ++i ) {
        mask.dilateByOne();
        filename = aol::strprintf ( "%smask%d.pgm", this->getSaveDirectory(), i+1 );
        mask.save ( filename.c_str() );
      }

      for ( int i = 0; i < mask.size(); ++i ) {
        if ( mask[i] )
          levelsetFunction[i] = 0.;
        else
          levelsetFunction[i] = 1.;
      }
    }
    // or by constructing a signed distance function.
    else {
      qc::SignedDistanceOp<ConfiguratorType> signedDistOp( _grid );
      signedDistOp.apply( levelsetFunction, levelsetFunction );
      levelsetFunction.addToAll ( - 2*_grid.H() );
      levelsetFunction.threshold ( 0, 0, 1 );
    }

    this->setSaveName( "u2" );
    this->saveTimestepPNGDrawRedIsoline ( -1, -1, 0., levelsetFunction, _u0, _grid );

    // The Dirichlet boundary model only uses the levelsetFunction to build the narrow band,
    // therefore we don't need to redistance it to make it compatible with the chosen H type.
    if ( _parser.hasVariable ( "doTrueInpainting" ) == false ) {
      levelsetFunction.addToAll ( -0.5 );
      qc::SignedDistanceOp<ConfiguratorType> signedDistOp( _grid );
      signedDistOp.apply( levelsetFunction, levelsetFunction );
    }

    heavisideFunctionWeightedMassMat.setZero();
    heavisideFunctionWeightedMassOp.assembleAddMatrix( heavisideFunctionWeightedMassMat );
    aol::DiagonalBlockOp<RealType> heavisideFunctionWeightedMassMatBlockOp ( heavisideFunctionWeightedMassMat );

    if ( _parser.hasVariable ( "doTrueInpainting" ) ) {
      typedef nb::NarrowBandGrid<typename ConfiguratorType::InitType, qc::QC_2D> NarrowBandGridType;
      typedef nb::FullGridConfiguratorHull<ConfiguratorType, NarrowBandGridType > NarrowBandConfType;

      NarrowBandGridType narrowGrid ( _grid );
      NarrowBandConfType narrowConf ( narrowGrid );

      nb::initNarrowBandFromSubLevelset<NarrowBandConfType, ConfiguratorType> ( levelsetFunction, narrowGrid );

      // Test if the narrow band is properly initialized.
      qc::BitArray<qc::QC_2D> narrowMask ( qc::GridSize<qc::QC_2D>::createFrom ( _grid ) );
      narrowConf.writeNodeExistMaskTo ( narrowMask );
      string filename = aol::strprintf ( "%smask-nb.pgm", this->getSaveDirectory() );
      narrowMask.save ( filename.c_str() );

      qc::BitArray<qc::QC_2D> dirichletMask ( qc::GridSize<qc::QC_2D>::createFrom ( _grid ) );
      narrowGrid.extractEdges();
      for ( typename NarrowBandGridType::biterator bit = narrowGrid.bbegin(); bit != narrowGrid.bend(); ++bit ) {
        // In case originalMask.get ( *bit ) == 0., the velocity value at *bit is destroyed.
        // So we can't use such nodes as Dirichlet nodes.
        if ( originalMask.get ( *bit ) > 0. )
          dirichletMask.set ( *bit, true );
      }

      filename = aol::strprintf ( "%smask-nb-boundary.pgm", this->getSaveDirectory() );
      dirichletMask.save ( filename.c_str() );

      VFieldInpaintingEnergy<NarrowBandConfType, PeronaMalikFunctionType, HeavisideFunctionType>
        E( narrowGrid, _u0, _v0, peronaMalikFunction, _params, heavisideFunctionWeightedMassMatBlockOp, levelsetFunction, heavisideFunction, false, false );
      VariationOfVFieldInpaintingEnergy<NarrowBandConfType, PeronaMalikFunctionType, HeavisideFunctionType>
        DE( narrowGrid, _u0, _v0, peronaMalikFunction, _params, heavisideFunctionWeightedMassMatBlockOp, levelsetFunction, heavisideFunction, false, false );

      aol::L2GradientDescentDirichletBCs<NarrowBandConfType, aol::MultiVector<RealType> > gradient_solver( narrowGrid, E, DE, dirichletMask, _params.numSteps );
      gradient_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "v", this->_saveDirectory );
      gradient_solver.setApplyMetric ( false );
      gradient_solver.setConfigurationFlags ( aol::L2GradientDescentDirichletBCs<NarrowBandConfType, aol::MultiVector<RealType> >::USE_NONLINEAR_CG );
      aol::MultiVector<RealType> mtmp ( velocity, aol::STRUCT_COPY );
      gradient_solver.apply( velocity, mtmp );
      velocity = mtmp;
    }
    else {
      const bool weightRegTermWithOneMinusH = _parser.hasVariable ( "weightRegTermWithOneMinusH" );

      if ( weightRegTermWithOneMinusH ) {
        levelsetFunction.addToAll ( - _parser.getDouble( "rho" ) *_grid.H() );
        ArrayType newMask ( levelsetFunction, aol::DEEP_COPY );
        newMask.threshold ( 0, 0, 1 );
        this->setSaveName( "regMask" );
        this->savePNG ( newMask, _grid );
      }

      VFieldInpaintingEnergy<ConfiguratorType, PeronaMalikFunctionType, HeavisideFunctionType>
        E( _grid, _u0, _v0, peronaMalikFunction, _params, heavisideFunctionWeightedMassMatBlockOp, levelsetFunction, heavisideFunction, true, weightRegTermWithOneMinusH );
      VariationOfVFieldInpaintingEnergy<ConfiguratorType, PeronaMalikFunctionType, HeavisideFunctionType>
        DE( _grid, _u0, _v0, peronaMalikFunction, _params, heavisideFunctionWeightedMassMatBlockOp, levelsetFunction, heavisideFunction, true, weightRegTermWithOneMinusH );

      /*
      aol::FirstDerivativeValidator<aol::MultiVector<RealType> > tester ( E, DE, _grid.H(), aol::FirstDerivativeValidator<aol::MultiVector<RealType> >::LINEAR, 0.0001 );
      tester.testDirection ( velocity, "test/test" );
      tester.setSkipPlottingWithSmallError( true );
      tester.setSkippingThreshold( 1e-8 );
      tester.testAllDirections( velocity, "test/test" );
      */

      //aol::GradientDescentWithAutomaticFilterWidth<ConfiguratorType, aol::MultiVector<RealType> >
      aol::GradientDescent<ConfiguratorType, aol::MultiVector<RealType> >
        gradientDescent_solver( _grid, E, DE, _params.numSteps );

      gradientDescent_solver.activateSavingAndConfigure( 2, _parser.getInt( "saveOffset"), "v", this->_saveDirectory );
      aol::MultiVector<RealType> mtmp ( velocity, aol::STRUCT_COPY );
      gradientDescent_solver.apply( velocity, mtmp );
      velocity = mtmp;
    }

    this->setSaveName( "v_col" );
    this->writeColorField( velocity, _grid );
    this->setSaveName( "v_vec" );
    this->plotVectorField ( velocity, _grid, aol::GNUPLOT_PNG );
    this->plotVectorField ( velocity, _grid, aol::GNUPLOT_PS );
    this->plotVectorField ( velocity, _grid, aol::GNUPLOT_EPS );
    this->setSaveName( "v_vec" );
    this->saveTimestepBZ2 ( -1, velocity, _grid );

  }   // end computeTimesteps
};


#endif // __VFINPAINT_H
