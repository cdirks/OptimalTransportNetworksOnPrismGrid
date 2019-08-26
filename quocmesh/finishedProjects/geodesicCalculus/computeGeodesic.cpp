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

//! \f$ E[s_0,...,s_K, \phi_1, ... , \phi_K ] = \sum_{k=1}^K \int_\Omega W(D\phi_k) dx + \mu \sum_{k=0}^K \int_\Omega | \nabla s_k |_{2,\delta} d x  + \gamma \sum_{k=1}^K \int_\Omega | s_k \circ \phi_k - \phi_{k-1}|^2 d x \f$
template <typename ConfiguratorType, typename MaterialLawType>
class GeodesicEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const int _numPolygonLines, _numShapes;
  const RealType _gamma, _mu, _epsilon, _delta;
  qc::LinearConvolutionSmoothOp<RealType,ConfiguratorType::Dim> _gaussKernelSmoother;
  const MaterialLawType _elasticEnergyDensity;
  const aol::IsoEnergyOp<ConfiguratorType> _lengthEnergy;

public:
  GeodesicEnergy( const typename ConfiguratorType::InitType &Grid,
                  const int NumberOfDeformations,
                  const RealType Gamma,
                  const RealType Mu,
                  const RealType Epsilon,
                  const RealType Delta,
                  const RealType LengthEnergyWeight,
                  const RealType VolumeEnergyWeight ) :
    _grid( Grid ),
    _numPolygonLines( NumberOfDeformations ),
    _numShapes( _numPolygonLines + 1 ),
    _gamma( Gamma ),
    _mu( Mu ),
    _epsilon( Epsilon ),
    _delta( Delta ),
    _gaussKernelSmoother( Grid.getNumX(), Grid.getNumY() ),
    _elasticEnergyDensity( LengthEnergyWeight, 0, VolumeEnergyWeight ),
    _lengthEnergy( Grid, Delta ) {
    _gaussKernelSmoother.setSigma( Epsilon );
  }

  // Arg has (dim+1)*K+1 components, ie. K+1 characteristic functions and dim*K displacement components
  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    // convolve characteristic functions with a Gaussian
    aol::MultiVector<RealType> smoothedCharFncs( _numShapes, _grid.getNumberOfNodes() );
    for ( int k = 0; k < _numShapes; k++ )
      _gaussKernelSmoother.apply( Arg[k], smoothedCharFncs[k] );

    // add length regularization
    aol::Scalar<RealType> aux;
    for ( int k = 1; k < _numPolygonLines; k++ )
#ifdef TV_OF_SMOOTHED_CHAR_FUN
      _lengthEnergy.applyAdd( smoothedCharFncs[k], aux ); // alternative: RegLengthEnergy<ConfiguratorType>( _grid, _epsilon, _delta ).applyAdd( Arg[k], aux );
#else
      _lengthEnergy.applyAdd( Arg[k], aux );
#endif
    Dest.addMultiple( aux, _mu );

    // add elastic energies
    for ( int k = 1; k < _numShapes; k++ ) {
      aol::MultiVector<RealType> displacement;
      displacement.appendReference( Arg[_numShapes+ConfiguratorType::Dim*(k-1)] );
      displacement.appendReference( Arg[_numShapes+ConfiguratorType::Dim*(k-1)+1] );
      qc::HyperelasticEnergy<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity ).applyAdd( displacement, Dest );
    }

    // add matching energies
    aux.setZero();
    for ( int k = 1; k < _numShapes; k++ ) {
      aol::MultiVector<RealType> arg;
      arg.appendReference( Arg[_numShapes+ConfiguratorType::Dim*(k-1)] );
      arg.appendReference( Arg[_numShapes+ConfiguratorType::Dim*(k-1)+1] );
      qc::MismatchEnergyWRTDisp<ConfiguratorType>( _grid, smoothedCharFncs[k-1], smoothedCharFncs[k] ).applyAdd( arg, aux );
    }
    Dest.addMultiple( aux, _gamma );
  }
};

//! derivative w.r.t. shapes AND displacements of GeodesicEnergy<>
template <typename ConfiguratorType, typename MaterialLawType>
class GeodesicGradient :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::MultiVector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const int _numPolygonLines, _numShapes;
  const RealType _gamma, _mu, _epsilon, _delta;
  qc::LinearConvolutionSmoothOp<RealType,ConfiguratorType::Dim> _gaussKernelSmoother;
  const MaterialLawType _elasticEnergyDensity;
  const RegLengthEnergyVariation<ConfiguratorType> _lengthEnergyVariation;

public:
  GeodesicGradient( const typename ConfiguratorType::InitType &Grid,
                    const int NumberOfDeformations,
                    const RealType Gamma,
                    const RealType Mu,
                    const RealType Epsilon,
                    const RealType Delta,
                    const RealType LengthEnergyWeight,
                    const RealType VolumeEnergyWeight ) :
    _grid( Grid ),
    _numPolygonLines( NumberOfDeformations ),
    _numShapes( _numPolygonLines + 1 ),
    _gamma( Gamma ),
    _mu( Mu ),
    _epsilon( Epsilon ),
    _delta( Delta ),
    _gaussKernelSmoother( Grid.getNumX(), Grid.getNumY() ),
    _elasticEnergyDensity( LengthEnergyWeight, 0, VolumeEnergyWeight ),
    _lengthEnergyVariation( Grid, Epsilon, Delta ) {
    _gaussKernelSmoother.setSigma( Epsilon );
  }

  // Arg has (dim+1)*K+1 components, ie. K+1 characteristic functions and dim*K displacement components
  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // convolve characteristic functions with a Gaussian
    aol::MultiVector<RealType> smoothedCharFncs( _numShapes, _grid.getNumberOfNodes() );
    for ( int k = 0; k < _numShapes; k++ )
      _gaussKernelSmoother.apply( Arg[k], smoothedCharFncs[k] );

    // add length regularization variations
    aol::Vector<RealType> aux( Dest[0], aol::STRUCT_COPY );
    for ( int k = 1; k < _numPolygonLines; k++ ) {
#ifdef TV_OF_SMOOTHED_CHAR_FUN
      _lengthEnergyVariation.apply( Arg[k], aux );
#else
      aol::VariationOfIsoEnergyOp<ConfiguratorType>( _grid, _delta ).apply( Arg[k], aux );
#endif
      Dest[k].addMultiple( aux, _mu );
    }

    // add elastic energy variation
    for ( int k = 1; k < _numShapes; k++ ) {
      aol::MultiVector<RealType> displacement, dest;
      displacement.appendReference( Arg[_numShapes+ConfiguratorType::Dim*(k-1)] );
      displacement.appendReference( Arg[_numShapes+ConfiguratorType::Dim*(k-1)+1] );
      dest.appendReference( Dest[_numShapes+ConfiguratorType::Dim*(k-1)] );
      dest.appendReference( Dest[_numShapes+ConfiguratorType::Dim*(k-1)+1] );
      qc::HyperelasticGradient<ConfiguratorType,MaterialLawType>( _grid, _elasticEnergyDensity ).applyAdd( displacement, dest );
    }

    // add matching energy variation
    aux.setZero();
    for ( int k = 1; k < _numShapes; k++ ) {
      aol::MultiVector<RealType> displacement;
      displacement.appendReference( Arg[_numShapes+ConfiguratorType::Dim*(k-1)] );
      displacement.appendReference( Arg[_numShapes+ConfiguratorType::Dim*(k-1)+1] );
      aol::MultiVector<RealType> dispGrad( displacement, aol::STRUCT_COPY );
      // variation wrt displacement
      qc::MismatchEnergyVariationWRTDisp<ConfiguratorType>( _grid, smoothedCharFncs[k-1], smoothedCharFncs[k] ).applyAdd( displacement, dispGrad );
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        Dest[_numShapes+ConfiguratorType::Dim*(k-1)+i].addMultiple( dispGrad[i], _gamma );
      // variation wrt characteristic fncs
      aol::Vector<RealType> aux( Arg[0], aol::STRUCT_COPY );
      MismatchGradientWRTnondefCharFun<ConfiguratorType>( _grid, smoothedCharFncs[k], displacement, _epsilon ).apply( smoothedCharFncs[k-1], aux );
      Dest[k-1].addMultiple( aux, _gamma );
      MismatchGradientWRTdefCharFun<ConfiguratorType>( _grid, smoothedCharFncs[k-1], displacement, _epsilon ).apply( smoothedCharFncs[k], aux );
      Dest[k].addMultiple( aux, _gamma );
    }

    // set gradient wrt u0 and uK to zero
    Dest[0].setZero();
    Dest[_numPolygonLines].setZero();
  }
};

//!
typedef double RealType;
const qc::Dimension Dim = qc::QC_2D;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, Dim, aol::GaussQuadrature<RealType, Dim, 3> > ConfiguratorType;

// main
//! \brief Computes discrete geodesic \f$ s_0, ..., S_K \f$ of length K+1 between two given input shapes \f$ s_A =: s_0 \f$ und  \f$ s_B =: S_K \f$.
//! \author Wirth (also Boerdgen, Heeren)
//!
//! Shapes are represented by characteristic functions which have been smoothed by some Gaussian kernel.
//! We assume there are K elatic deformations \f$ \phi_k: s_{k-1} \rightarrow s_k \f$ that fulfill the matching condition \f$ \phi_k( s_{k-1} ) = s_k \f$.
//! The objective energy E depends on all (intermediate) shapes and all deformations, ie
//! \f$ E[s_0,...,s_K, \phi_1, ... , \phi_K ] = \sum_{k=1}^K \int_\Omega W(D\phi_k) dx + \mu \sum_{k=0}^K \int_\Omega | \nabla s_k |_{2,\delta} d x  + \gamma \sum_{k=1}^K \int_\Omega | s_k \circ \phi_k - \phi_{k-1}|^2 d x \f$
//! As we want to fix the first and last shape s_0 and s_K, the corresponding gradient components are set to zero while minimzing.
int main( int argc, char *argv[] ) {

  if ( argc < 11 )
    cerr<<"usage: "<<argv[0]<<" <image1> <image2> <numShapes> <gamma> <mu> "
        <<"<elasticLengthModulus> <elasticVolumeModulus> <epsilon/h> <delta> <solutionWordStem>\n"
        <<"e.g. "<<argv[0]<<" image1.pgm image2.pgm 3 10 1 1 1 5 1e-10 solution1_\n";
  else
    try {

      // load given and initial data
      qc::ScalarArray<RealType,Dim> u0( argv[1] ), uK( argv[2] );
      u0.addToAll( - u0.getMinValue() );
      u0 /= u0.getMaxValue();
      uK.addToAll( - uK.getMinValue() );
      uK /= uK.getMaxValue();
      int numShapes = atoi( argv[3] );
      qc::GridSize<Dim> gridSize( u0 );
      ConfiguratorType::InitType grid( gridSize );
      RealType gamma = strtod( argv[4], NULL ), mu = strtod( argv[5], NULL );
      RealType epsilon = strtod( argv[8], NULL ) * grid.H(), delta = strtod( argv[9], NULL );
      RealType lengthEnergyWeight = strtod( argv[6], NULL ), volumeEnergyWeight = strtod( argv[7], NULL );

      // define the energy and gradient
      GeodesicEnergy<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > geodesicEnergy( grid, numShapes - 1, gamma, mu, epsilon, delta, lengthEnergyWeight, volumeEnergyWeight );
      GeodesicGradient<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > geodesicGradient( grid, numShapes - 1, gamma, mu, epsilon, delta, lengthEnergyWeight, volumeEnergyWeight );

      // do gradient descent
      aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> > gradientDescent( grid, geodesicEnergy, geodesicGradient, 500 );
      gradientDescent.setConfigurationFlags( aol::GradientDescent<ConfiguratorType,aol::MultiVector<RealType> >::USE_NONLINEAR_CG );
      aol::MultiVector<RealType> initialData( numShapes * ( Dim + 1 ) - Dim, grid.getNumberOfNodes() ), solution( initialData, aol::STRUCT_COPY );
      for ( int k = 0; k < numShapes - 1; k++ )
        initialData[k] = u0;
      initialData[numShapes-1] = uK;

      gradientDescent.apply( initialData, solution );

      // save results
      for ( int k = 0; k < numShapes; k++ ) {
        char fileName[256];
        sprintf( fileName, "%s%d.pgm", argv[10], k );
        solution[k] *= 255;
        ConfiguratorType::ArrayType( solution[k], grid ).save( fileName, qc::PGM_UNSIGNED_CHAR_BINARY );
        // also save the smoothed characteristic functions
        aol::Vector<RealType> smoothedSol( solution[k], aol::STRUCT_COPY );
        qc::LinearConvolutionSmoothOp<RealType,ConfiguratorType::Dim> smoothOp( grid.getNumX(), grid.getNumY() );
        smoothOp.setSigma( epsilon );
        smoothOp.apply( solution[k], smoothedSol );
        sprintf( fileName, "%sSmoothed_%d.pgm", argv[10], k );
        ConfiguratorType::ArrayType( smoothedSol, grid ).save( fileName, qc::PGM_UNSIGNED_CHAR_BINARY );
        // save deformations
        if ( k != 0 ) {
          sprintf( fileName, "%sdisplacement_%d_0.pgm", argv[10], k );
          ConfiguratorType::ArrayType( solution[numShapes+Dim*(k-1)], grid ).save( fileName, qc::PGM_DOUBLE_BINARY );
          sprintf( fileName, "%sdisplacement_%d_1.pgm", argv[10], k );
          ConfiguratorType::ArrayType( solution[numShapes+Dim*(k-1)+1], grid ).save( fileName, qc::PGM_DOUBLE_BINARY );
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
