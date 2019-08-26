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
template <typename ConfiguratorType>
class SquaredWeightMassOp :
  public aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredWeightMassOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  // the weight to be squared
  const aol::DiscreteFunctionDefault<ConfiguratorType> _w;
public:
  SquaredWeightMassOp( const typename ConfiguratorType::InitType &Grid,
                       const aol::Vector<RealType> &W,
                       aol::OperatorType OpType = aol::ONTHEFLY ) :
    // initialise the grid
    aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredWeightMassOp<ConfiguratorType> >( Grid, OpType ),
    // load the weight
    _w( Grid, W ) {}

  /**
   * \brief Returns \f$w^2\f$ evaluated at the point specified by the element "El" and quadrature point "QuadPoint".
   */
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::VecType& /*RefCoord*/ ) const {
    return aol::Sqr( _w.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};



/**
 * \brief This class computes the energy \f$\int_\Omega\frac12((\phi-1)^2(u-c_1)^2+(\phi+1)^2(u-c_2)^2)+\nu(\epsilon\|\nabla\phi\|^2+\frac1\epsilon(|\phi|^2-1)^2)dx\f$,
 * where the real-valued image \f$ u \f$, the greyscales \f$ c_1, c_2 \f$, and the parameters \f$ \nu, \epsilon \f$ are passed to the constructor.
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class ModicaMortolaSegmentationEnergy
  : public aol::FENonlinIntegrationScalarInterface<ConfiguratorType, ModicaMortolaSegmentationEnergy<ConfiguratorType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the image
  const aol::DiscreteFunctionDefault<ConfiguratorType> _u;
  // the greyscales
  const RealType _c1, _c2;
  // the parameters
  const RealType _epsilon, _nu;

public:
  ModicaMortolaSegmentationEnergy ( const typename ConfiguratorType::InitType &Grid,
                                    const aol::Vector<RealType> &U,
                                    const RealType C1,
                                    const RealType C2,
                                    const RealType Epsilon,
                                    const RealType Nu )
    : aol::FENonlinIntegrationScalarInterface<ConfiguratorType, ModicaMortolaSegmentationEnergy<ConfiguratorType> > ( Grid ),
      _u( Grid, U ),
      _c1( C1 ),
      _c2( C2 ),
      _epsilon( Epsilon ),
      _nu( Nu ) {}

  /**
   * \brief Returns \f$\frac12((\phi-1)^2(u-c_1)^2+(\phi+1)^2(u-c_2)^2)+\nu(\epsilon\|\nabla\phi\|^2+\frac1\epsilon(|\phi|^2-1)^2)\f$
   * evaluated at the point specified by the element "El" and quadrature point "QuadPoint".
   */
  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    const RealType phi = DiscFunc.evaluateAtQuadPoint( El, QuadPoint );
    const RealType u = _u.evaluateAtQuadPoint( El, QuadPoint );
    typename ConfiguratorType::VecType gradPhi;
    DiscFunc.evaluateGradientAtQuadPoint( El, QuadPoint, gradPhi );
    return .5 * ( aol::Sqr( phi - 1 ) * aol::Sqr( u - _c1 ) + aol::Sqr( phi + 1 ) * aol::Sqr( u - _c2 ) )
           + _nu * ( _epsilon * gradPhi.normSqr() + aol::Sqr( aol::Sqr( phi ) - 1 ) / _epsilon );
  }
};



/**
 * \brief This class computes the energy variation
 * \f$\left(\int_\Omega((\phi-1)(u-c_1)^2+(\phi+1)(u-c_2)^2+\nu\frac4\epsilon(|\phi|^2-1)\phi)\psi_i+2\nu\epsilon\nabla\phi\cdot\nabla\psi_idx\right)_i\f$,
 * where the real-valued image \f$ u \f$, the greyscales \f$ c_1, c_2 \f$, and the parameters \f$ \nu, \epsilon \f$ are passed to the constructor.
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class ModicaMortolaSegmentationEnergyVariation
  : public aol::FENonlinOpInterface<ConfiguratorType, ModicaMortolaSegmentationEnergyVariation<ConfiguratorType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the image
  const aol::DiscreteFunctionDefault<ConfiguratorType> _u;
  // the greyscales
  const RealType _c1, _c2;
  // the parameters
  const RealType _epsilon, _nu;
  // the stiffness operator for the "stiffness" term
  const aol::StiffOp<ConfiguratorType> _stiffOp;

public:
  ModicaMortolaSegmentationEnergyVariation ( const typename ConfiguratorType::InitType &Grid,
                                             const aol::Vector<RealType> &U,
                                             const RealType C1,
                                             const RealType C2,
                                             const RealType Epsilon,
                                             const RealType Nu )
    : aol::FENonlinOpInterface<ConfiguratorType, ModicaMortolaSegmentationEnergyVariation<ConfiguratorType> > ( Grid ),
      _u( Grid, U ),
      _c1( C1 ),
      _c2( C2 ),
      _epsilon( Epsilon ),
      _nu( Nu ),
      _stiffOp( Grid ) {}

  /**
   * \brief Returns \f$(\phi-1)(u-c_1)^2+(\phi+1)(u-c_2)^2+\nu\frac4\epsilon(|\phi|^2-1)\phi\f$
   * evaluated at the point specified by the element "El" and quadrature point "QuadPoint".
   */
  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                             const typename ConfiguratorType::ElementType &El,
                             int QuadPoint,
                             const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                             typename ConfiguratorType::RealType &NL ) const {
    const RealType phi = DiscFunc.evaluateAtQuadPoint( El, QuadPoint );
    const RealType u = _u.evaluateAtQuadPoint( El, QuadPoint );
    NL = ( phi - 1) * aol::Sqr( u - _c1 ) + ( phi + 1) * aol::Sqr( u - _c2 ) + 4 * _nu / _epsilon * ( aol::Sqr( phi ) - 1 ) * phi;
  }

  /**
   * \brief This method computes the energy variation.
   */
  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // compute "mass" term of the energy variation
    aol::FENonlinOpInterface<ConfiguratorType, ModicaMortolaSegmentationEnergyVariation<ConfiguratorType> >::applyAdd( Arg, Dest );

    // compute the "stiffness" term
    if ( _nu ) {
      Dest /= _epsilon * _nu * 2;
      _stiffOp.applyAdd( Arg, Dest );
      Dest *= _epsilon * _nu * 2;
    }
  }
};



/**
 * \brief This class computes the second energy variation
 * \f$ \left(\int_\Omega((u-c_1)^2+(u-c_2)^2+\nu\frac4\epsilon(3|\phi|^2-1))\psi_i\psi_j+2\nu\epsilon\nabla\psi_i\cdot\nabla\psi_jdx\right)_{ij}\f$,
 * where the real-valued image \f$ u \f$, the greyscales \f$ c_1, c_2 \f$, and the parameters \f$ \nu, \epsilon \f$ are passed to the constructor.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename MatrixType>
class ModicaMortolaSegmentationEnergySecondVariation : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, MatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  // the grid
  const typename ConfiguratorType::InitType _grid;
  // the image
  const ArrayType _u;
  // the greyscales
  const RealType _c1, _c2;
  // the parameters
  const RealType _epsilon, _nu;

public:
  ModicaMortolaSegmentationEnergySecondVariation ( const typename ConfiguratorType::InitType &Grid,
                                                   const aol::Vector<RealType> &U,
                                                   const RealType C1,
                                                   const RealType C2,
                                                   const RealType Epsilon,
                                                   const RealType Nu )
    : _grid( Grid ),
      _u( U, Grid ),
      _c1( C1 ),
      _c2( C2 ),
      _epsilon( Epsilon ),
      _nu( Nu ) {}

  /**
   * \brief This method computes \f$M[(u-c_1)^2+(u-c_2)^2+4\frac\nu\epsilon(3\phi^2-1)]+2\nu\epsilon L\f$.
   */
  virtual void applyAdd( const aol::Vector<RealType> &Arg, MatrixType &Dest ) const{
    ArrayType uC1( _u ), uC2( _u );
    uC1.addToAll( -_c1 );
    uC2.addToAll( -_c2 );
    SquaredWeightMassOp<ConfiguratorType>( _grid, uC1 ).assembleAddMatrix( Dest, 1 );
    SquaredWeightMassOp<ConfiguratorType>( _grid, uC2 ).assembleAddMatrix( Dest, 1 );
    aol::StiffOp<ConfiguratorType>( _grid ).assembleAddMatrix( Dest, 2 * _nu * _epsilon );
    aol::MassOp<ConfiguratorType>( _grid ).assembleAddMatrix( Dest, -4 * _nu / _epsilon );
    SquaredWeightMassOp<ConfiguratorType>( _grid, Arg ).assembleAddMatrix( Dest, 12 * _nu / _epsilon );
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
    if ( argc != 4 ) {
      // too many input arguments specified; explain correct syntax
      cerr << "Incorrect usage. Use: " << argv[0] << " <image filename> <parameter 1> <parameter 2>" << endl;
      return EXIT_FAILURE;
    }
    else {
      // load the parameters
      RealType nu = atof( argv[2] ),
               lambda = atof( argv[3] );

      // execute the algorithm

      // load the image (must be quadratic with width 2^n+1)
      ArrayType u( argv[1] );
      int numX = u.getNumX();
      int numY = u.getNumY();
      GridType grid( aol::Vec3<int>( numX, numY, 1 ) );
      RealType epsilon = grid.H();
      u *= 1. / u.getMaxValue();

      // initialize the double well phi, psi, and the greyscales c1, c2
      ArrayType phi( grid ),
                psi( grid );
      phi.setAll( 1 );
      psi.setAll( 1 );
      RealType c1 = 0, c2 = 1;

      // perform the fixed point iteration between phi and psi
      int numIter = 10;
      for ( int iter = 0; iter < numIter; ++iter ) {

        RealType eps = epsilon*2;// * 10 / ( iter + 1 );
        // enhance complementarity of phi and psi
        lambda = lambda * 2;

        // assemble the system matrices and rhs
        qc::FastUniformGridMatrix<RealType, ConfiguratorType::Dim> systemMatrix( grid ), massMatrix( grid ), stiffMatrix( grid );
        aol::Vector<RealType> rhs( grid ), ones( grid );

        aol::MassOp<ConfiguratorType>( grid ).assembleAddMatrix( massMatrix, nu * ( lambda - ( 2. / eps ) ) );
        aol::StiffOp<ConfiguratorType>( grid ).assembleAddMatrix( stiffMatrix, nu * eps );
        ones.setAll( 1. );

        //! compute phi
        // produce intermediate matrix
        ArrayType shiftedImage1( u );
        shiftedImage1.addToAll( -c1 );
        SquaredWeightMassOp<ConfiguratorType>( grid, shiftedImage1 ).assembleAddMatrix( systemMatrix, 1 );

        // produce rhs
        massMatrix.apply( psi , rhs );
        rhs *= -lambda / ( lambda - ( 2. / eps ) );
        systemMatrix.applyAdd( ones, rhs );

        // produce the system matrix
        systemMatrix += massMatrix;
        systemMatrix += stiffMatrix;
        SquaredWeightMassOp<ConfiguratorType>( grid, psi ).assembleAddMatrix( systemMatrix, nu * 2. / eps );

        // solve for phi
        aol::CGInverse<aol::Vector<RealType> > inv1( systemMatrix );
        inv1.setStopping( aol::STOPPING_ABSOLUTE );
        inv1.apply( rhs, phi );

        //! compute psi
        systemMatrix.setZero();
        rhs.setZero();
        // produce intermediate matrix
        ArrayType shiftedImage2( u );
        shiftedImage2.addToAll( -c2 );
        SquaredWeightMassOp<ConfiguratorType>( grid, shiftedImage2 ).assembleAddMatrix( systemMatrix, 1 );

        // produce rhs
        massMatrix.apply( phi , rhs );
        rhs *= -lambda / ( lambda - ( 2. / eps ) );
        systemMatrix.applyAdd( ones, rhs );

        // produce the system matrix
        systemMatrix += massMatrix;
        systemMatrix += stiffMatrix;
        SquaredWeightMassOp<ConfiguratorType>( grid, phi ).assembleAddMatrix( systemMatrix, nu * 2. / eps );

        // solve for phi
        aol::CGInverse<aol::Vector<RealType> > inv2( systemMatrix );
        inv2.setStopping( aol::STOPPING_ABSOLUTE );
        inv2.apply( rhs, psi );

// moderately sharpen the interface
phi *= 2;
psi *= 2;
for ( int i = 0; i < phi.size(); i++ ) {
  if ( phi[i] > 1 ) phi[i] = 1; else if ( phi[i] < -1 ) phi[i] = -1;
  if ( psi[i] > 1 ) psi[i] = 1; else if ( psi[i] < -1 ) psi[i] = -1;
}
/*
for ( int i = 0; i < phi.getSize(); i++ ) {
  if ( phi[i] > 0 ) phi[i] = 1; else phi[i] = -1;
  if ( psi[i] > 0 ) psi[i] = 1; else psi[i] = -1;
}
*/
        //! compute the grey scales
        ArrayType phi_1( phi );
        ArrayType psi_1( psi );

        phi_1.addToAll( -1 );
        psi_1.addToAll( -1 );

        /*
        for ( int i = 0; i < phi.getSize(); i++ ) {
          phi_1[i] = (phi[i]-1)*(psi[i]+1);
          psi_1[i] = (psi[i]-1)*(phi[i]+1);
        }
        */

        SquaredWeightMassOp<ConfiguratorType>( grid, phi_1 ).apply( ones, rhs );  // rhs is used as auxiliary vector
        c1 = ( rhs * u ) / ( rhs * ones );
        SquaredWeightMassOp<ConfiguratorType>( grid, psi_1 ).apply( ones, rhs );  // rhs is used as auxiliary vector
        c2 = ( rhs * u ) / ( rhs * ones );

        //! save the result
        cerr<<phi.getMinValue()<<" "<<phi.getMaxValue()<<" "<<psi.getMinValue()<<" "<<psi.getMaxValue()<<" "<<c1<<" "<<c2<<endl;
        qc::QuocTimestepSaver<ConfiguratorType> timeStepSaver( 1, 0, "doubleWell", true );
        timeStepSaver.setSaveDirectory( "ModMor/" );
        timeStepSaver.saveTimestepPNGDrawRedAndBlueIsoline( iter, 0, 0., .9, phi, u, grid );
        timeStepSaver.saveTimestepPNGDrawRedAndBlueIsoline( iter, 1, 0., .9, psi, u, grid );
      }


      //! now employ the real Modica-Mortola model
      psi.setAll( 1 );
      c1 = 0;
      c2 = 1;

      // perform the fixed point iteration between phi and the greyscales
      for ( int iter = 0; iter < numIter; ++iter ) {
        /*
        // define total energy $E$=:"E" and its variation with respect to the double well, "DE"
        ModicaMortolaSegmentationEnergy<ConfiguratorType> E( grid, u, c1, c2, 2 * epsilon, nu );
        ModicaMortolaSegmentationEnergyVariation<ConfiguratorType> DE( grid, u, c1, c2, 2 * epsilon, nu );
        // perform a few steps of gradient descent
        aol::GradientDescentWithAutomaticFilterWidth<ConfiguratorType,aol::Vector<RealType> > gradientDescent( grid, E, DE, 200 );
        gradientDescent.apply( psi, phi );
        */
        // define variation of Energy with respect to the double well, "DE", and the second derivative, "D2E"
        ModicaMortolaSegmentationEnergyVariation<ConfiguratorType> DE( grid, u, c1, c2, 2 * epsilon, nu );
        ModicaMortolaSegmentationEnergySecondVariation<ConfiguratorType, qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> > D2E( grid, u, c1, c2, 2 * epsilon, nu );
        // perform a few steps of Newton's method
        aol::NewtonIteration<ConfiguratorType> newtonDescent( grid , DE, D2E, 2 );
        newtonDescent.apply( psi, phi );

        // compute the greyscales
        ArrayType phiPlusOne( phi ), phiMinusOne( phi ), ones( grid ), tmp( grid );
        phiPlusOne.addToAll( 1. );
        phiMinusOne.addToAll( -1. );
        ones.setAll( 1. );

        SquaredWeightMassOp<ConfiguratorType>( grid, phiMinusOne ).apply( ones, tmp );
        c1 = ( tmp * u ) / ( tmp * ones );
        SquaredWeightMassOp<ConfiguratorType>( grid, phiPlusOne ).apply( ones, tmp );
        c2 = ( tmp * u ) / ( tmp * ones );

        // save the result
        psi = phi;
        cerr<<phi.getMinValue()<<" "<<phi.getMaxValue()<<" "<<c1<<" "<<c2<<endl;
        qc::QuocTimestepSaver<ConfiguratorType> timeStepSaver( 1, 0, "singleDoubleWell", true );
        timeStepSaver.setSaveDirectory( "ModMor/" );
        timeStepSaver.saveTimestepPNGDrawRedAndBlueIsoline( iter, 0, 0., .9, phi, u, grid );
        timeStepSaver.setSaveName( "phaseField" );
        timeStepSaver.savePNG( phi, grid, iter, 0 );
      }


      //! now employ the real Modica-Mortola model, but use a linearization in each iteration
      psi.setAll( 1 );
      c1 = 0;
      c2 = 1;

      // perform the fixed point iteration between phi and the greyscales
      for ( int iter = 0; iter < numIter; ++iter ) {

        RealType eps = epsilon * 2 * 40 / ( iter + 1 );

        // assemble the system matrices and rhs
        qc::FastUniformGridMatrix<RealType, ConfiguratorType::Dim> systemMatrix( grid ), auxMatrix( grid );
        aol::Vector<RealType> rhs( grid ), ones( grid );
        ones.setAll( 1. );

        //! compute phi
        // produce intermediate matrix
        ArrayType shiftedImage1( u ), shiftedImage2( u );
        shiftedImage1.addToAll( -c1 );
        shiftedImage2.addToAll( -c2 );
        SquaredWeightMassOp<ConfiguratorType>( grid, shiftedImage1 ).assembleAddMatrix( systemMatrix, 1 );
        SquaredWeightMassOp<ConfiguratorType>( grid, shiftedImage2 ).assembleAddMatrix( auxMatrix, 1 );

        // produce rhs
        auxMatrix.apply( ones, rhs );
        rhs *= -1;
        systemMatrix.applyAdd( ones, rhs );

        // produce the system matrix
        systemMatrix += auxMatrix;
        aol::StiffOp<ConfiguratorType>( grid ).assembleAddMatrix( systemMatrix, 2 * nu * eps );
        aol::MassOp<ConfiguratorType>( grid ).assembleAddMatrix( systemMatrix, -4 * nu / eps );
        SquaredWeightMassOp<ConfiguratorType>( grid, phi ).assembleAddMatrix( systemMatrix, 4 * nu / eps  );

        // solve for phi
        aol::CGInverse<aol::Vector<RealType> > inv( systemMatrix );
        inv.setStopping( aol::STOPPING_ABSOLUTE );
        inv.apply( rhs, phi );

// moderately sharpen the interface
phi *= 2;
psi *= 2;
for ( int i = 0; i < phi.size(); i++ ) {
  if ( phi[i] > 1 ) phi[i] = 1; else if ( phi[i] < -1 ) phi[i] = -1;
  if ( psi[i] > 1 ) psi[i] = 1; else if ( psi[i] < -1 ) psi[i] = -1;
}

        //! compute the greyscales
        ArrayType phiPlusOne( phi ), phiMinusOne( phi );
        phiPlusOne.addToAll( 1. );
        phiMinusOne.addToAll( -1. );

        SquaredWeightMassOp<ConfiguratorType>( grid, phiMinusOne ).apply( ones, rhs );  // rhs is used as auxiliary vector
        c1 = ( rhs * u ) / ( rhs * ones );
        SquaredWeightMassOp<ConfiguratorType>( grid, phiPlusOne ).apply( ones, rhs );  // rhs is used as auxiliary vector
        c2 = ( rhs * u ) / ( rhs * ones );

        //! save the result
        cerr<<phi.getMinValue()<<" "<<phi.getMaxValue()<<" "<<c1<<" "<<c2<<endl;
        qc::QuocTimestepSaver<ConfiguratorType> timeStepSaver( 1, 0, "linearizedDoubleWell", true );
        timeStepSaver.setSaveDirectory( "ModMor/" );
        timeStepSaver.saveTimestepPNGDrawRedAndBlueIsoline( iter, 0, 0., .9, phi, u, grid );
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
