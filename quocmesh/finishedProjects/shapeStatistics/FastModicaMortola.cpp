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

#include <aol.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <configurators.h>
#include <quocTimestepSaver.h>
#include <smallVec.h>
#include <matrix.h>
#include <gradientDescent.h>
#include <Newton.h>


/**
 * \brief This class represents a scaled mass operator, assembling the matrix \f$\left(\int_\Omega w^2\psi_i\psi_j\right)_{ij}\f$,
 * where the real-valued function \f$ w \f$ is passed to the constructor. \f$\psi_{i,j}\f$ represent the \f$i\f$-th and \f$j\f$-th FE basis function.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename WeightFunctionType>
class SquaredWeightMassOp :
  public aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredWeightMassOp<ConfiguratorType,WeightFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  // the weight to be squared
  const WeightFunctionType &_w;
public:
  SquaredWeightMassOp( const typename ConfiguratorType::InitType &Grid,
                       const WeightFunctionType &W,
                       aol::OperatorType OpType = aol::ONTHEFLY ) :
    // initialise the grid
    aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredWeightMassOp<ConfiguratorType,WeightFunctionType> >( Grid, OpType ),
    // load the weight
    _w( W ) {}

  /**
   * \brief Returns \f$w^2\f$ evaluated at the point specified by the element "El" and quadrature point "QuadPoint".
   */
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::VecType& /*RefCoord*/ ) const {
    return aol::Sqr( _w.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

template <typename ConfiguratorType, typename WeightFunctionType>
class WeightedSquaredFunctionIntegration
  : public aol::FENonlinIntegrationScalarInterface<ConfiguratorType, WeightedSquaredFunctionIntegration<ConfiguratorType, WeightFunctionType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the weight
  const WeightFunctionType &_w;

public:
  WeightedSquaredFunctionIntegration ( const typename ConfiguratorType::InitType &Grid,
                                       const WeightFunctionType &W )
    : aol::FENonlinIntegrationScalarInterface<ConfiguratorType, WeightedSquaredFunctionIntegration<ConfiguratorType, WeightFunctionType> > ( Grid ),
      _w( W ) {}

  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    const RealType f = DiscFunc.evaluateAtQuadPoint( El, QuadPoint );
    const RealType w = _w.evaluateAtQuadPoint( El, QuadPoint );
    return aol::Sqr( f ) * w;
  }
};

template <typename RealType>
class Convolution1dFiniteSupport {

  void convolve( const qc::Array<RealType> &U, const aol::Vector<RealType> &Kernel, const qc::Comp Dir, qc::Array<RealType> &Dest ) {
    int dimZ = U.getNumZ(), dimY = U.getNumY(), dimX = U.getNumX();
    // int kernSupp = Kernel.size();
    for ( int z = 0; z < dimZ; z++ )
      for ( int y = 0; y < dimY; y++ )
        for ( int x = 0; x < dimX; x++ )
          Dest( x, y, z ) = convolveLine( U, Kernel, x, y, z, Dir );
  }

  RealType convolveLine( const qc::Array<RealType> &U, const aol::Vector<RealType> &Kernel, const int X, const int Y, const int Z, const qc::Comp Dir) {
    int dimZ = U.getNumZ(), dimY = U.getNumY(), dimX = U.getNumX();
    int kernSupp = Kernel.size(), supp;
    RealType result = 0;
    switch ( Dir ) {
      case qc::QC_X:
        supp = std::min( dimX - X, kernSupp );
        for ( int x = 0; x < supp; x++ )
          result += U( X + x, Y, Z ) * Kernel( x );
        break;
      case qc::QC_Y:
        supp = std::min( dimY - Y, kernSupp );
        for ( int y = 0; y < supp; y++ )
          result += U( X, Y + y, Z ) * Kernel( y );
        break;
      case qc::QC_Z:
        supp = std::min( dimZ - Z, kernSupp );
        for ( int z = 0; z < supp; z++ )
          result += U( X, Y, Z + z ) * Kernel( z );
        break;
    }
  }
};

template <typename ConfiguratorType>
class DiscreteFunction : public aol::DiscreteFunctionInterface<ConfiguratorType, DiscreteFunction<ConfiguratorType> > {
public:
  typedef ConfiguratorType ConfType;
  typedef typename ConfType::RealType RealType;
  typedef typename ConfType::VecType  VecType;
  typedef typename ConfType::DomVecType  DomVecType;
  typedef typename ConfType::ElementType ElementType;

  const aol::DiscreteFunctionDefault<ConfiguratorType> &_function;

  DiscreteFunction ( const typename ConfiguratorType::InitType &Initializer, const aol::DiscreteFunctionDefault<ConfiguratorType> &Function ) :
    aol::DiscreteFunctionInterface<ConfiguratorType, DiscreteFunction<ConfiguratorType> > ( Initializer ),
    _function ( Function ) {
  }

  //! copy constructor
  DiscreteFunction ( const DiscreteFunction<ConfiguratorType> &DiscrFunc ) :
    aol::DiscreteFunctionInterface<ConfiguratorType, DiscreteFunction<ConfiguratorType> > ( DiscrFunc._initializer ),
    _function ( DiscrFunc._function ) {
  }

  RealType evaluate ( const ElementType& El, const DomVecType& RefCoord ) const {
    ElementType newEl;
    DomVecType newLocalCoord;
    switchCoordinates( *(this->_configPtr), El, RefCoord, _function.getConfigurator(), newEl, newLocalCoord );
    return _function.evaluate( newEl, newLocalCoord );
  }

  RealType evaluateAtQuadPoint ( const ElementType &El, int QuadPoint ) const {
    ElementType newEl;
    DomVecType localCoord = (this->_configPtr)->getBaseFunctionSet( El ).getRefCoord( QuadPoint ), newLocalCoord;
    switchCoordinates( *(this->_configPtr), El, localCoord, _function.getConfigurator(), newEl, newLocalCoord );
    return _function.evaluate( newEl, newLocalCoord );
  }

  void evaluateGradient ( const ElementType& El, const DomVecType& RefCoord, VecType& Grad ) const {
    ElementType newEl;
    DomVecType newLocalCoord;
    switchCoordinates( *(this->_configPtr), El, RefCoord, _function.getConfigurator(), newEl, newLocalCoord );
    _function.evaluateGradient( newEl, newLocalCoord, Grad );
  }

  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, aol::Vec<ConfiguratorType::DomDim, RealType>& Grad ) const {
    ElementType newEl;
    DomVecType localCoord = *(this->_configPtr).getBaseFunctionSet( El ).getRefCoord( QuadPoint ), newLocalCoord;
    switchCoordinates( *(this->_configPtr), El, localCoord, _function.getConfigurator(), newEl, newLocalCoord );
    _function.evaluateGradient( newEl, newLocalCoord, Grad );
  }

private:
  void switchCoordinates ( const ConfType &Config, const ElementType &El, const DomVecType &LocalCoord,
                           const ConfType &NewConfig, ElementType &NewEl, DomVecType &NewLocalCoord) const {
    //NewEl = NewConfig.getEmptyElement();
    //NewLocalCoord.setZero();

    // find global coordinates
    VecType globalCoord;
    Config.getGlobalCoords ( El, LocalCoord, globalCoord );
    // find grid element and local coordinates for _function
    NewConfig.getLocalCoords ( globalCoord, NewEl, NewLocalCoord );
  }
};


//////////////////////////////////////////////////////////////////////////////////////
//! ---------------------------- main program --------------------------------------//
//////////////////////////////////////////////////////////////////////////////////////

// define the settings/configuration of the program, i.e. computation accuracy, dimension, gridtype, finite element type, used quadrature rules
typedef double RealType;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfiguratorType;
typedef ConfiguratorType::ArrayType ArrayType;
typedef ConfiguratorType::InitType GridType;

/**
 * \brief This main function coordinates the segmentation of an image via Modica-Mortola.
 * The parameter data is read in and the segmentation algorithm is executed.
 * If started from the command line, argc is the number of passed arguments,
 * and argv is a pointer on the passed arguments (which are of type char[]).
 * The name of the image file should be passed, as well as two parameters.
 *
 * \author Wirth
 */
int main( int argc, char *argv[] ) {

  try {
    if ( argc != 3 ) {
      // too many input arguments specified; explain correct syntax
      cerr << "Incorrect usage. Use: " << argv[0] << " <image filename> <grid depth of computation>" << endl;
      return EXIT_FAILURE;
    }
    else {
      // load the image (must be quadratic with width 2^n+1)
      ArrayType u( argv[1] );
      u *= 1. / u.getMaxValue();
      int numX = u.getNumX();
      int numY = u.getNumY();
      GridType imageGrid( aol::Vec3<int>( numX, numY, 1 ) );

      // employ the fast Modica-Mortola model on the chosen grid depth
      GridType grid( atoi( argv[2] ) );
      ArrayType phi( grid );
      phi.setAll( 1 );
      aol::Scalar<RealType> c1 = 0, c2 = 1, tmp;
      int numIter = 10;

      // perform the fixed point iteration between phi and the greyscales
      for ( int iter = 0; iter < numIter; ++iter ) {
        // assemble system matrix and rhs
        qc::FastUniformGridMatrix<RealType, ConfiguratorType::Dim> systemMatrix( grid );
        aol::Vector<RealType> rhs( grid ), ones( grid );
        ones.setAll( 1 );
        ArrayType uC1( u ), uC2( u );
        uC1.addToAll( -c1 );
        uC2.addToAll( -c2 );
        aol::DiscreteFunctionDefault<ConfiguratorType> uC1Fine( imageGrid, uC1 );
        aol::DiscreteFunctionDefault<ConfiguratorType> uC2Fine( imageGrid, uC2 );
        DiscreteFunction<ConfiguratorType> uC1Coarse( grid, uC1Fine );
        DiscreteFunction<ConfiguratorType> uC2Coarse( grid, uC2Fine );
        SquaredWeightMassOp<ConfiguratorType,DiscreteFunction<ConfiguratorType> >( grid, uC2Coarse ).assembleAddMatrix( systemMatrix, 1 );
        systemMatrix.apply( ones, rhs );
        rhs *= -2;
        SquaredWeightMassOp<ConfiguratorType,DiscreteFunction<ConfiguratorType> >( grid, uC1Coarse ).assembleAddMatrix( systemMatrix, 1 );
        systemMatrix.applyAdd( ones, rhs );

        // solve for phi
        aol::CGInverse<aol::Vector<RealType> > inv( systemMatrix );
        inv.setStopping( aol::STOPPING_ABSOLUTE );
        inv.apply( rhs, phi );

        // compute the greyscales
        ArrayType phiPlusOne( phi ), phiMinusOne( phi );
        phiPlusOne.addToAll( 1. );
        phiMinusOne.addToAll( -1. );

        aol::DiscreteFunctionDefault<ConfiguratorType> uFine( imageGrid, u );
        DiscreteFunction<ConfiguratorType> uCoarse( grid, uFine );
        aol::DiscreteFunctionDefault<ConfiguratorType> one( grid, ones );
        WeightedSquaredFunctionIntegration<ConfiguratorType,DiscreteFunction<ConfiguratorType> >( grid, uCoarse ).apply( phiPlusOne, c2 );
        WeightedSquaredFunctionIntegration<ConfiguratorType,DiscreteFunction<ConfiguratorType> >( grid, uCoarse ).apply( phiMinusOne, c1 );
        WeightedSquaredFunctionIntegration<ConfiguratorType,aol::DiscreteFunctionDefault<ConfiguratorType> >( grid, one ).apply( phiPlusOne, tmp );
        c2 /= tmp;
        WeightedSquaredFunctionIntegration<ConfiguratorType,aol::DiscreteFunctionDefault<ConfiguratorType> >( grid, one ).apply( phiMinusOne, tmp );
        c1 /= tmp;

        // save the result
        cerr<<phi.getMinValue()<<" "<<phi.getMaxValue()<<" "<<c1<<" "<<c2<<endl;
        qc::QuocTimestepSaver<ConfiguratorType> timeStepSaver( 1, 0, "fastModMor", true );
        timeStepSaver.setSaveDirectory( "ModMor/" );
        timeStepSaver.saveTimestepPNGDrawRedAndBlueIsoline( iter, 0, 0., .9, phi, phi, grid );
        timeStepSaver.setSaveName( "phaseField" );
        timeStepSaver.savePNG( phi, grid, iter, 0 );
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
