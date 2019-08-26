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
#include <parameterParser.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <deformations.h>
#include <quocTimestepSaver.h>
#include <ChanVese.h>
#include <mcm.h>
#include <gradientDescent.h>
#include <signedDistanceOp.h>
#include <smallVec.h>
#include <gradientflow.h>
#include <hyperelastic.h>

template< typename ConfiguratorType >
class LevelsetGenerator {
public:
  /**
    * \brief Makes "Levelset" the signed distance function to the ball around ("CentreX","CentreY","CentreZ")
    * with radius "Radius".
    */
  static void generateBallLevelset ( typename qc::Array<typename ConfiguratorType::RealType> &Levelset,
                                     typename ConfiguratorType::RealType Radius = 0.3,
                                     typename ConfiguratorType::RealType CentreX = 0.5,
                                     typename ConfiguratorType::RealType CentreY = 0.5,
                                     typename ConfiguratorType::RealType CentreZ = 0.5 ) {
    // grid width
    const int width = Levelset.getNumX();
    for ( int i = 0; i < width; ++i )          // x-components
      for ( int j = 0; j < width; ++j )        // y-components
        if ( ConfiguratorType::Dim == 3 )
          for ( int k = 0; k < width; ++k ) {  // z-components
            const typename ConfiguratorType::RealType
              // compute the current coordinates
              x = ( 1.0 * i ) / ( width - 1.0 ),
              y = ( 1.0 * j ) / ( width - 1.0 ),
              z = ( 1.0 * k ) / ( width - 1.0 ),
              // compute the signed distance function value
              value = aol::Sqr ( x - CentreX ) + aol::Sqr ( y - CentreY ) + aol::Sqr ( z - CentreZ ) - aol::Sqr ( Radius );
            Levelset.set ( i, j, k, value );
          }
        else {
          const typename ConfiguratorType::RealType
              // compute the current coordinates
            x = ( 1.0 * i ) / ( width - 1.0 ),
            y = ( 1.0 * j ) / ( width - 1.0 ),
            // compute the signed distance function value
            value = aol::Sqr ( x - CentreX ) + aol::Sqr ( y - CentreY ) - aol::Sqr ( Radius );
          Levelset.set ( i, j, value );
        }
  }
};

/**
 * \brief This class encapsulates multiple static functions which are needed to compute the energy and its variation
 * and which can be used by the different classes during the process of segmentation and average shape finding.
 * The class hence only serves to avoid redundancy of code.
 *
 * \author Wirth
 */
template <typename ConfigurationType, typename HeavisideFunctionType>
class EvaluationFunctions {

public:
  typedef typename ConfigurationType::RealType RealType;

  /**
   * \brief This function computes the displaced point \f$\phi(x)\f$, where \f$x\f$ is specified by
   * "Grid", "El", "QuadPoint", and "LocalCoord". The result is represented by "DeformedEl" and "LocalDeformedCoord".
   * The deformation \f$\phi\f$ is the displacement plus identity.
   * The method returns true, if \f$x\f$ is displaced outside the grid, else false.
   */
  inline static bool getDeformedCoord( const aol::auto_container<ConfigurationType::Dim,aol::DiscreteFunctionDefault<ConfigurationType> > &Displacement,
                                       const typename ConfigurationType::InitType &Grid,
                                       const typename ConfigurationType::ElementType &El,
                                       const int &QuadPoint,
                                       const typename ConfigurationType::VecType &LocalCoord,
                                       typename ConfigurationType::ElementType &DeformedEl,
                                       typename ConfigurationType::VecType &LocalDeformedCoord ) {
    // find global coordinates
    typename ConfigurationType::VecType coord;
    Displacement[0].getConfigurator().getGlobalCoords ( El, LocalCoord, coord );
    // find deformed coordinates $\phi(x)$
    typename ConfigurationType::VecType deformedCoord( coord );
    for ( int j=0; j<ConfigurationType::Dim; j++ )
      deformedCoord[j] += Displacement[j].evaluateAtQuadPoint( El, QuadPoint );
    // compute corresponding grid element and local coordinates
    // (makes use of the fact that $\phi$ should map the grid onto itself)
    DeformedEl = Displacement[0].getConfigurator().getEmptyElement();
    Displacement[0].getConfigurator().getLocalCoords ( deformedCoord, DeformedEl, LocalDeformedCoord );
    // check, whether point was displaced outside the grid
    // (makes use of the fact that $\phi$ should map the grid onto itself)
    for ( int j=0 ; j<ConfigurationType::Dim; j++ )
      if ( LocalDeformedCoord[j] < 0. || ( DeformedEl[j] >= ( Grid.getWidth() - 1 ) ) )
        return true;
    return false;
  }

  /**
   * \brief This function computes \f$H(\zeta(\phi_i(x)))\f$ at the point \f$x\f$ specified by the grid element "El"
   * and the quadrature point "QuadPoint", which is also represented by the local element coordinates "LocalCoord".
   * The grid on which the computation takes place is passed via "Grid", the used Heaviside function via
   * "HeavisideFunction", the levelset function \f$\zeta\f$ by "LevelsetFunction", and the displacements \f$d_i\f$ by "Displacement".
   * The resulting value \f$H(\zeta(\phi_i(x)))\f$ is returned.
   */
  inline static aol::Scalar<RealType> evaluateHZetaPhi( const HeavisideFunctionType &HeavisideFunction,
                                                        const aol::DiscreteFunctionDefault<ConfigurationType> &LevelsetFunction,
                                                        const aol::auto_container<ConfigurationType::Dim,aol::DiscreteFunctionDefault<ConfigurationType> > &Displacement,
                                                        const typename ConfigurationType::InitType &Grid,
                                                        const typename ConfigurationType::ElementType &El,
                                                        const int &QuadPoint,
                                                        const typename ConfigurationType::VecType &LocalCoord,
                                                        bool &DisplacedOutsideDomain ) {
    // find deformed element and local coordinates
    typename ConfigurationType::VecType localDeformedCoord;
    typename ConfigurationType::ElementType deformedEl;
    DisplacedOutsideDomain = getDeformedCoord( Displacement, Grid, El, QuadPoint, LocalCoord, deformedEl, localDeformedCoord );
    // compute $H(\zeta(\phi(x)))$; if $x$ was displaced outside the domain, extend by zero
    if ( !DisplacedOutsideDomain )
      return aol::Scalar<RealType> (HeavisideFunction.evaluate( LevelsetFunction.evaluate( deformedEl, localDeformedCoord ) ) );
    else
      return aol::Scalar<RealType> (0);
  }

  /**
   * \brief This function computes the vector \f$H'(\zeta(\phi_i(x)))\nabla\zeta(\phi_i(x))\f$ at the point \f$x\f$ specified by the grid element "El"
   * and the quadrature point "QuadPoint", which is also represented by the local element coordinates "LocalCoord".
   * The grid on which the computation takes place is passed via "Grid", the used Heaviside function via
   * "HeavisideFunction", the levelset function \f$\zeta\f$ by "LevelsetFunction", the displacement \f$d_i\f$ by "Displacement".
   * The resulting vector is returned.
   */
  inline static aol::Vec<ConfigurationType::Dim,RealType> evaluateDHZetaPhi( const HeavisideFunctionType &HeavisideFunction,
                                                          const aol::DiscreteFunctionDefault<ConfigurationType> &LevelsetFunction,
                                                          const aol::auto_container<ConfigurationType::Dim,aol::DiscreteFunctionDefault<ConfigurationType> > &Displacement,
                                                          const typename ConfigurationType::InitType &Grid,
                                                          const typename ConfigurationType::ElementType &El,
                                                          const int &QuadPoint,
                                                          const typename ConfigurationType::VecType &LocalCoord,
                                                          bool &DisplacedOutsideDomain ) {
    aol::Vec<ConfigurationType::Dim,RealType> DHZP;
    // find deformed element and local coordinates
    typename ConfigurationType::VecType localDeformedCoord;
    typename ConfigurationType::ElementType deformedEl;
    DisplacedOutsideDomain = getDeformedCoord( Displacement, Grid, El, QuadPoint, LocalCoord, deformedEl, localDeformedCoord );
    // compute $H(\zeta(\phi(x)))$; if $x$ was displaced outside the domain, extend by zero (DHZP already initialised with zero)
    if ( !DisplacedOutsideDomain ) {
      LevelsetFunction.evaluateGradient( deformedEl, localDeformedCoord, DHZP );
      DHZP *= HeavisideFunction.evaluateDerivative( LevelsetFunction.evaluate( deformedEl, localDeformedCoord ) );
    }
    return DHZP;
  }
};

/**
 * \brief Interface for classes to compute an integral of the form \f$\int_{\Omega}V(\psi(x))dx\f$,
 * where the integration domain \f$\Omega\f$ is defined by the member variable _grid, the "ResultType"-
 * valued function \f$V\f$ is implemented in the subclasses by the interface method "evaluateIntegrand(...)",
 * and the piecewise multilinear, "VariableType"-valued function \f$\psi\f$ is passed as a multivector or
 * vector to the method "apply(...)" or "applyAdd(...)", which computes the integral.
 * "ResultType" can be a scalar or a vector. "VariableType" can be a vector or a multivector.
 *
 * \author Wirth
 */
template <typename ConfigurationType, typename VariableType, int VariableDim, typename ResultType, typename Imp>
class VectorIntegratorOverDomainInterface :
  public aol::Op<VariableType, ResultType> {

protected:
  // dummy variable containing the general settings; used to find the first and last grid element
  mutable ConfigurationType _configuration;
  // grid representing the domain of integration (the function domain)
  const typename ConfigurationType::InitType &_grid;

  typedef typename ConfigurationType::RealType RealType;

public:
  VectorIntegratorOverDomainInterface( const typename ConfigurationType::InitType &Grid ) :
    // load the settings
    _configuration( Grid ),
    // load the grid, representing the integration domain
    _grid( Grid ) {
  }

  virtual ~VectorIntegratorOverDomainInterface() {}

  /**
   * \brief This method computes the integral \f$\int_{\Omega}V(\psi(x))dx\f$, where the integration domain \f$\Omega\f$ is
   * specified by the class member variable _grid, the vector-valued function \f$V\f$ is implemented in the method
   * "evaluateIntegrand(...)", and the piecewise multilinear, vector-valued function \f$\psi\f$ is passed as the
   * multivector "Arg". The vector-valued result is added to the vector "Dest".
   */
  void applyAdd ( const VariableType &Arg, ResultType &Dest ) const {
    // transform the multivector or vector "Arg" into a vector of piecewise multilinear functions
    aol::MultiVector<RealType> arg;
    arg.appendReference( Arg );
    aol::auto_container<VariableDim,aol::DiscreteFunctionDefault<ConfigurationType> > psi;
    for ( int i = 0; i<VariableDim; i++ )
      psi.set_copy( i, aol::DiscreteFunctionDefault<ConfigurationType> ( _grid, arg[i] ) );
    // for each grid element do ...
    for ( typename ConfigurationType::ElementIteratorType it = _configuration.begin(); it != _configuration.end(); ++it ) {
      ResultType integralOverElement;
      // compute integrand at the first quadrature point (later for the rest)
      this->asImp().evaluateIntegrand ( psi, *it, 0, _configuration.getBaseFunctionSet( *it ).getRefCoord( 0 ), integralOverElement );
      // weight the integrand value at the first quadrature point
      integralOverElement *= _configuration.getBaseFunctionSet ( *it ).getWeight ( 0 );
      // for each remaining quadrature point on the grid element do ...
      for ( int q = 1; q < ConfigurationType::QuadType::numQuadPoints; q++ ) {
        ResultType integrandAtQuadPoint;
        // compute integrand at the quadrature point
        this->asImp().evaluateIntegrand ( psi, *it, q, _configuration.getBaseFunctionSet( *it ).getRefCoord( q ), integrandAtQuadPoint );
        // add the weighted integrand value at the quadrature point
        integrandAtQuadPoint *= _configuration.getBaseFunctionSet ( *it ).getWeight ( q );
        integralOverElement += integrandAtQuadPoint;
      }
      // make the weighted sum of integrand values the integral over the element
      integralOverElement *= _configuration.vol ( *it );
      // add the integral over the grid element to the total integral
      Dest += integralOverElement;
    }
  }

  /**
   * \brief This interface method has to be implemented by any subclass. It implements the vector-valued
   * function \f$V(\psi(x))\f$, where the vector-valued piecewise multilinear function \f$\psi\f$ is passed
   * as the container "Psi" of piecewise multilinear scalar functions. The method computes \f$V(\psi(x))\f$
   * at the point \f$x\f$, which is specified by the grid element "El" and the quadrature point on the
   * element with number "QuadPoint" or alternatively the local coordinates on the element "LocalCoord".
   * The function value is supposed to be passed to the parameter "Integrand".
   */
  void evaluateIntegrand ( const aol::auto_container<VariableDim,aol::DiscreteFunctionDefault<ConfigurationType> > &Psi,
                           const typename ConfigurationType::ElementType &El,
                           const int QuadPoint,
                           const typename ConfigurationType::VecType &LocalCoord,
                           ResultType &Integrand ) const {
    throw aol::Exception ( "Called the interface function ""evaluateIntegrand(...)""", __FILE__, __LINE__ );
    this->asImp().evaluateIntegrand ( Psi, El, QuadPoint, LocalCoord, Integrand );
  }

protected:
  // Barton-Nackman trick
  inline Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  inline const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }
};

/**
 * \brief As a subclass of the interface "VectorIntegratorOverDomainInterface", this class represents
 * an integrator of the form \f$\int_{\Omega}V(\psi(x))dx for vector- and scalar-valued (discrete) functions
 * \f$V\f$ and \f$\psi\f$. When integrating via the method "apply(...)" or "applyAdd(...)", the discrete scalar
 * function \f$\psi\f$ is passed as an argument. Here, \f$\psi\f$ is meant to be any image,
 * and \f$V(\psi(x))=(x^T,1)^T\psi(x)\f$. Thus, by dividing the first entries of the resulting integral by the last entry,
 * one obtains the centre of gravity, assuming that \f$\psi\f$ ranges from 0 to 1, 1 being highest and 0 being zero density.
 *
 * \author Wirth
 */
template <typename ConfigurationType>
class MomentIntegrator :
  public VectorIntegratorOverDomainInterface<ConfigurationType, aol::Vector<typename ConfigurationType::RealType>, 1, aol::Vector<typename ConfigurationType::RealType>, MomentIntegrator<ConfigurationType> > {

private:
  typedef typename ConfigurationType::RealType RealType;

public:
  MomentIntegrator ( const typename ConfigurationType::InitType &Grid ) :
    // load the grid
    VectorIntegratorOverDomainInterface<ConfigurationType, aol::Vector<RealType>, 1, aol::Vector<RealType>, MomentIntegrator<ConfigurationType> >( Grid ) {}

  /**
   * \brief This method evaluates a vector-valued function \f$V(\psi(x))=(x^T,1)^T(1-\psi(x))\f$, where the piecewise multilinear
   * function \f$\psi\f$ is passed as the single-element-vector "Psi" of piecewise multilinear scalar functions.
   */
  void evaluateIntegrand ( const aol::auto_container<1,aol::DiscreteFunctionDefault<ConfigurationType> > &Psi,
                           const typename ConfigurationType::ElementType &El,
                           const int QuadPoint,
                           const typename ConfigurationType::VecType &LocalCoord,
                           aol::Vector<RealType> &Integrand ) const {
    // find global coordinates
    typename ConfigurationType::VecType coord;
    Psi[0].getConfigurator().getGlobalCoords ( El, LocalCoord, coord );

    // construct $(x^T,1)^T$="Integrand"
    Integrand.resize( ConfigurationType::Dim+1 );
    for ( int i=0; i<ConfigurationType::Dim; i++ )
      Integrand[i] = coord[i];
    Integrand[ConfigurationType::Dim] = 1;

    // compute and return the correct integrand
    Integrand *= (Psi[0].evaluateAtQuadPoint( El, QuadPoint ));
  }
};

/**
 * \brief As a subclass of the interface "VectorIntegratorOverDomainInterface", this class represents
 * an integrator of the form \f$\int_{\Omega}V(\psi(x))dx\f$ for vector-valued (discrete) functions
 * \f$V\f$ and \f$\psi\f$. When integrating via the method "apply(...)" or "applyAdd(...)", the discrete
 * function \f$\psi\f$ is passed as an argument. Here, \f$\psi\f$ is meant to be the levelset function
 * \f$\zeta\f$, and \f$V(\zeta(x))\f$ depends on the computation mode:
 * SIMPLE       refers to \f$H(\zeta(\phi(x)))\f$
 * INVERSE      refers to \f$1-H(\zeta(\phi(x)))\f$
 * GREY         refers to \f$u(x)H(\zeta(\phi(x)))\f$
 * GREY_INVERSE refers to \f$u(x)(1-H(\zeta(\phi(x))))\f$
 * where the deformation \f$\phi\f$ can be represented with the help of the displacement \f$d\f$ as identity\f$+d\f$.
 * The discrete displacement \f$d\f$ and image \f$u\f$ have to be passed in the constructor. \f$d\f$ and \f$u\f$ are
 * here meant to represent chosen elements of all displacements \f$d_i\f$ and images \f$u_i\f$.
 *
 * \author Wirth
 */
template <typename ConfigurationType, typename HeavisideFunctionType>
class HZPIntegratorOverDomain : // mnemonic for Heaviside Zeta Phi
  public VectorIntegratorOverDomainInterface<ConfigurationType, aol::Vector<typename ConfigurationType::RealType>, 1, aol::Scalar<typename ConfigurationType::RealType>, HZPIntegratorOverDomain<ConfigurationType, HeavisideFunctionType> > {

public:
  // provides the following computation modes: For an image $u$, levelset function $\zeta$ and deformation $\phi$,
  // SIMPLE       refers to $\int_{\Omega}H(\zeta(\phi(x)))dx$
  // INVERSE      refers to $\int_{\Omega}1-H(\zeta(\phi(x)))dx$
  // GREY         refers to $\int_{\Omega}u(x)H(\zeta(\phi(x)))dx$
  // GREY_INVERSE refers to $\int_{\Omega}u(x)(1-H(\zeta(\phi(x))))dx$
  enum MODE { SIMPLE, INVERSE, GREY, GREY_INVERSE };

private:
  typedef typename ConfigurationType::RealType RealType;

  // the image $u_i$
  const aol::Vector<RealType> &_image;
  aol::DiscreteFunctionDefault<ConfigurationType> _imageFunction;
  // the displacement $d_i$
  const aol::MultiVector<RealType> &_d;
  aol::auto_container<ConfigurationType::Dim,aol::DiscreteFunctionDefault<ConfigurationType> > _displacement;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // the computation mode
  const MODE _mode;

public:
  HZPIntegratorOverDomain ( const typename ConfigurationType::InitType &Grid,
                            const aol::Vector<RealType> &Image,
                            const aol::MultiVector<RealType> &D,
                            const HeavisideFunctionType &HeavisideFunction,
                            const MODE Mode ) :
    // load the grid
    VectorIntegratorOverDomainInterface<ConfigurationType, aol::Vector<RealType>, 1, aol::Scalar<RealType>, HZPIntegratorOverDomain<ConfigurationType, HeavisideFunctionType> >( Grid ),
    // load the image
    _image( Image ),
    // generate a multilinear interpolant from the discrete image
    _imageFunction( Grid, Image ),
    // load the displacement
    _d( D ),
     // initialise the displacement-interpolant
    _displacement(),
   // load the smoothed version of the Heaviside function to be used
    _heavisideFunction( HeavisideFunction ),
    // define the computation mode
    _mode( Mode ) {
    // generate a multilinear interpolant from the discrete displacement
    for ( int j=0; j<ConfigurationType::Dim; j++ )
      _displacement.set_copy( j, aol::DiscreteFunctionDefault<ConfigurationType> ( Grid, D[j] ) );
  }

  /**
   * \brief This method evaluates a vector-valued function \f$V(\psi(x))\f$, where the piecewise multilinear
   * function \f$\psi\f$ is passed as the single-element-vector "Psi" of piecewise multilinear scalar functions. Here, \f$\psi\f$ is
   * meant to be the levelset function \f$\zeta\f$, and \f$V(\zeta(x))\f$ depends on the computation mode:
   * SIMPLE       refers to \f$H(\zeta(\phi(x)))\f$
   * INVERSE      refers to \f$1-H(\zeta(\phi(x)))\f$
   * GREY         refers to \f$u(x)H(\zeta(\phi(x)))\f$
   * GREY_INVERSE refers to \f$u(x)(1-H(\zeta(\phi(x))))\f$
   * where the deformation \f$\phi\f$ can be represented with the help of the displacement \f$d\f$ as identity\f$+d\f$.
   * The discrete displacement \f$d\f$ and image \f$u\f$ are contained in the member variables "_displacement" and
   * "_imageFunction". \f$d\f$ and \f$u\f$ are here meant to represent chosen elements of all displacements \f$d_i\f$ and
   * images \f$u_i\f$.
   * The method computes \f$V(\psi(x))\f$ at the point \f$x\f$, which is specified by the grid element "El" and the
   * quadrature point on the element with number "QuadPoint" or alternatively the local coordinates on the
   * element "LocalCoord". The function value is passed to the parameter "Integrand".
   */
  void evaluateIntegrand ( const aol::auto_container<1,aol::DiscreteFunctionDefault<ConfigurationType> > &Psi,
                           const typename ConfigurationType::ElementType &El,
                           const int QuadPoint,
                           const typename ConfigurationType::VecType &LocalCoord,
                           aol::Scalar<RealType> &Integrand ) const {
    // for the image compute $H(\zeta(\phi(x)))$
    aol::Scalar<RealType> HZP;
    bool displacedOutsideDomain = false;
    HZP = EvaluationFunctions<ConfigurationType, HeavisideFunctionType>::evaluateHZetaPhi( _heavisideFunction, Psi[0], _displacement, this->_grid, El, QuadPoint, LocalCoord, displacedOutsideDomain );

    // compute and return the correct integrand
    if ( displacedOutsideDomain )
      Integrand.setZero();
    else
      switch ( _mode ) {
        case SIMPLE      : Integrand = HZP;
                           break;
        case INVERSE     : Integrand = 1-HZP;
                           break;
        case GREY        : Integrand = HZP*_imageFunction.evaluateAtQuadPoint( El, QuadPoint );
                           break;
        case GREY_INVERSE: Integrand = (1-HZP)*_imageFunction.evaluateAtQuadPoint( El, QuadPoint );
                           break;
        default          : break;
      }
  }
};

/**
 * \brief This class represents the operator which computes the following energy as a function of the levelset function \f$\zeta\f$:
 * \f$E=\int_{\Omega}\lambda_1(1-H(\zeta(\phi_i(x))))(u_i(x)-c_1^i)^2+\lambda_2H(\zeta(\phi_i(x)))(u_i(x)-c_2^i)^2dx\f$
 * \f$\zeta\f$ represents the levelset function, \f$H\f$ the Heaviside function, \f$u_i\f$ the \f$i\f$-th image, and \f$\phi_i\f$ the \f$i\f$-th deformation.
 * The evaluation of the integral is performed by the method "apply(...)" or "applyAdd(...)", which uses the method
 * "evaluateIntegrand(...)" to evaluate the integrand at a point. The integrand can also be referred to as greyscale term.
 * Since the energy shall be an operator on \f$\zeta\f$, \f$H\f$, \f$\Omega\f$, \f$u_i\f$ and \f$\phi_i\f$ have to be passed to the constructor.
 *
 * \author Wirth
 */
template <typename ConfigurationType, typename HeavisideFunctionType>
class GreyscaleEnergyOfZeta : public VectorIntegratorOverDomainInterface<ConfigurationType, aol::Vector<typename ConfigurationType::RealType>, 1, aol::Scalar<typename ConfigurationType::RealType>, GreyscaleEnergyOfZeta<ConfigurationType, HeavisideFunctionType> > {

private:
  typedef typename ConfigurationType::RealType RealType;

  // the image $u_i$
  const aol::Vector<RealType> &_image;
  aol::DiscreteFunctionDefault<ConfigurationType> _imageFunction;
  // the displacement $d_i$
  const aol::MultiVector<RealType> &_d;
  aol::auto_container<ConfigurationType::Dim,aol::DiscreteFunctionDefault<ConfigurationType> > _displacement;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // Chan-Vese-parameters
  const RealType _lambda1, _lambda2;
  // greyscales of the image
  const RealType &_c1, &_c2;

public:
  GreyscaleEnergyOfZeta( const typename ConfigurationType::InitType &Grid,
                const aol::Vector<RealType> &Image,
                const aol::MultiVector<RealType> &D,
                const HeavisideFunctionType &HeavisideFunction,
                const RealType &C1,
                const RealType &C2,
                const RealType Lambda1,
                const RealType Lambda2 ) :
    // load the grid
    VectorIntegratorOverDomainInterface<ConfigurationType,aol::Vector<RealType>,1,aol::Scalar<RealType>,GreyscaleEnergyOfZeta<ConfigurationType, HeavisideFunctionType> >( Grid ),
    // load the image
    _image( Image ),
    // generate the image-interpolant
    _imageFunction( Grid, Image ),
    // load the displacement
    _d( D ),
    // initialise the displacement-interpolant
    _displacement(),
    // load the smoothed version of the Heaviside function to be used
    _heavisideFunction( HeavisideFunction ),
    // load the Chan-Vese parameters
    _lambda1( Lambda1 ),
    _lambda2( Lambda2 ),
    // load the greyscale vectors
    _c1( C1 ),
    _c2( C2 ) {
    // generate a multilinear interpolant from the discrete displacement
    for ( int j=0; j<ConfigurationType::Dim; j++ )
      _displacement.set_copy( j, aol::DiscreteFunctionDefault<ConfigurationType> ( Grid, D[j] ) );
  }

  /**
   * \brief This method evaluates the term
   * \f$\lambda_1(1-H(\zeta(\phi_i(x))))(u_i(x)-c_1^i)^2+\lambda_2H(\zeta(\phi_i(x)))(u_i(x)-c_2^i)^2\f$
   * at the point \f$x\f$ specified by the grid element "El" and the quadrature point on the element with
   * number "QuadPoint" or alternatively the local coordinates on the element "LocalCoord".
   * \f$\zeta\f$ is passed as the first component of "Psi", while \f$H\f$, \f$u_i\f$ and \f$\phi_i\f$ as well as the
   * Chan-Vese parameters are member variables of the class and are passed by the constructor.
   */
  void evaluateIntegrand ( const aol::auto_container<1,aol::DiscreteFunctionDefault<ConfigurationType> > &Psi,
                           const typename ConfigurationType::ElementType &El,
                           const int QuadPoint,
                           const typename ConfigurationType::VecType &LocalCoord,
                           aol::Scalar<RealType> &Integrand ) const {
    // initialise integrand with zero
    Integrand.setZero();

    // computation of the term $H(\zeta(\phi_i(x)))$
    bool displacedOutsideDomain = false;
    Integrand = EvaluationFunctions<ConfigurationType, HeavisideFunctionType>::evaluateHZetaPhi( _heavisideFunction, Psi[0], _displacement, this->_grid, El, QuadPoint, LocalCoord, displacedOutsideDomain );
    // compute the greyscale term
    if ( displacedOutsideDomain )
      Integrand.setZero();
    else
      Integrand = _lambda1*(1-Integrand)*aol::Sqr(_imageFunction.evaluateAtQuadPoint( El, QuadPoint )-_c1)
                     +_lambda2*Integrand*aol::Sqr(_imageFunction.evaluateAtQuadPoint( El, QuadPoint )-_c2);
  }
};

/**
 * \brief This class represents the operator which computes the following energy as a function of the levelset function \f$\zeta\f$:
 * \f$E=\int_{\Omega}\mu|\nabla H(\zeta(x))|+\nu H(\zeta(x))+\lambda_1\sum_i((1-H(\zeta(\phi_i(x))))(u_i(x)-c_1^i)^2)+\lambda_2\sum_i(H(\zeta(\phi_i(x)))(u_i(x)-c_2^i)^2)dx\f$
 * \f$\zeta\f$ represents the levelset function, \f$H\f$ the Heaviside function, \f$u_i\f$ the \f$i\f$-th image, and \f$\phi_i\f$ the \f$i\f$-th deformation.
 * The evaluation of the integral is performed by the method "apply(...)" or "applyAdd(...)", which uses the method
 * "evaluateIntegrand(...)" to evaluate the integrand at a point. The first term of the integrand can be referred to as
 * length/surface term, the second term as area/volume term, and the sums as greyscale term.
 * Since the energy shall be an operator on \f$\zeta\f$, \f$H\f$, \f$\Omega\f$, \f$u_i\f$ and \f$\phi_i\f$ have to be passed to the constructor.
 *
 * \author Wirth
 */
template <typename ConfigurationType, typename HeavisideFunctionType>
class EnergyOfZeta : public VectorIntegratorOverDomainInterface<ConfigurationType, aol::Vector<typename ConfigurationType::RealType>, 1, aol::Scalar<typename ConfigurationType::RealType>, EnergyOfZeta<ConfigurationType, HeavisideFunctionType> > {

private:
  typedef typename ConfigurationType::RealType RealType;

  // the images $u_i$
  const aol::MultiVector<RealType> &_images;
  const int _numberOfImages;
  // the displacements $d_i$
  const aol::auto_container<ConfigurationType::Dim, aol::MultiVector<RealType> > &_d;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // Chan-Vese-parameters
  const RealType _mu, _nu, _lambda1, _lambda2;
  // greyscales of the different images
  const aol::Vector<RealType> &_c1, &_c2;

public:
  EnergyOfZeta( const typename ConfigurationType::InitType &Grid,
                const aol::MultiVector<RealType> &Images,
                const aol::auto_container<ConfigurationType::Dim, aol::MultiVector<RealType> > &D,
                const HeavisideFunctionType &HeavisideFunction,
                const aol::Vector<RealType> &C1,
                const aol::Vector<RealType> &C2,
                const RealType Mu,
                const RealType Nu,
                const RealType Lambda1,
                const RealType Lambda2 ) :
    // load the grid
    VectorIntegratorOverDomainInterface<ConfigurationType,aol::Vector<RealType>,1,aol::Scalar<RealType>,EnergyOfZeta<ConfigurationType,HeavisideFunctionType> >( Grid ),
    // load the images
    _images( Images ),
    _numberOfImages( Images.numComponents() ),
    // load the displacement
    _d( D ),
    // load the smoothed version of the Heaviside function to be used
    _heavisideFunction( HeavisideFunction ),
    // load the Chan-Vese parameters
    _mu( Mu ),
    _nu( Nu ),
    _lambda1( Lambda1 ),
    _lambda2( Lambda2 ),
    // load the greyscale vectors
    _c1( C1 ),
    _c2( C2 ) {
  }

  /**
   * \brief This method evaluates the term \f$\mu|\nabla H(\zeta(x))|+\nu H(\zeta(x))\f$
   * at the point \f$x\f$ specified by the grid element "El" and the quadrature point on the element with
   * number "QuadPoint" or alternatively the local coordinates on the element "LocalCoord".
   * \f$\zeta\f$ is passed as the first component of "Psi", while \f$H\f$, \f$u_i\f$ and \f$\phi_i\f$ as well as the
   * Chan-Vese parameters are member variables of the class and are passed by the constructor.
   */
  void evaluateIntegrand ( const aol::auto_container<1,aol::DiscreteFunctionDefault<ConfigurationType> > &Psi,
                           const typename ConfigurationType::ElementType &El,
                           const int QuadPoint,
                           const typename ConfigurationType::VecType& /*LocalCoord*/,
                           aol::Scalar<RealType> &Integrand ) const {
    // initialise integrand with zero
    Integrand.setZero();

    //! computation of the length/surface term
    typename ConfigurationType::VecType grad;
    Psi[0].evaluateGradientAtQuadPoint( El, QuadPoint, grad );
    Integrand += _mu*_heavisideFunction.evaluateDerivative( Psi[0].evaluateAtQuadPoint( El, QuadPoint ) )*grad.norm();

    //! computation of the area/volume term
    Integrand += _nu*_heavisideFunction.evaluate( Psi[0].evaluateAtQuadPoint( El, QuadPoint ) );
  }

  /**
   * \brief This method computes the following energy as a function of the levelset function \f$\zeta\f$="Arg":
   * \f$E=\int_{\Omega}\mu|\nabla H(\zeta(x))|+\nu H(\zeta(x))+\lambda_1\sum_i((1-H(\zeta(\phi_i(x))))(u_i(x)-c_1^i)^2)+\lambda_2\sum_i(H(\zeta(\phi_i(x)))(u_i(x)-c_2^i)^2)dx\f$
   * The result is put into Dest.
   * \f$\zeta\f$ represents the levelset function, \f$H\f$ the Heaviside function, \f$u_i\f$ the \f$i\f$-th image, and \f$\phi_i\f$ the \f$i\f$-th deformation.
   * The first term of the integrand can be referred to as
   * length/surface term, the second term as area/volume term, and the sums as greyscale term.
   * The method uses the method "evaluateIntegrand(...)" to evaluate the length and area term of the integrand at a point.
   * Since the energy shall be an operator on \f$\zeta\f$, \f$H\f$, \f$\Omega\f$, \f$u_i\f$ and \f$\phi_i\f$ were to be passed to the constructor.
   */
  void apply( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    Dest.setZero();

    //! compute the greyscale term of the energy
    // for each image...
    for ( int i=0; i<_numberOfImages; i++ ){
      // reformulate the displacement
      aol::MultiVector<RealType> d(0,0);
      for ( int j=0; j < ConfigurationType::Dim; j++ )
        d.appendReference( _d[j][i] );
      GreyscaleEnergyOfZeta<ConfigurationType,HeavisideFunctionType> greyscaleEnergy( this->_grid, _images[i], d, _heavisideFunction, _c1[i], _c2[i], _lambda1, _lambda2 );
      greyscaleEnergy.applyAdd( Arg, Dest );
    }
    Dest /= _numberOfImages;

    //! compute the length/surface term of the energy
    // careful: apply of the parent class would not work, since it would call the virtual applyAdd of this class
    VectorIntegratorOverDomainInterface<ConfigurationType,aol::Vector<RealType>,1,aol::Scalar<RealType>,EnergyOfZeta<ConfigurationType,HeavisideFunctionType> >::applyAdd( Arg, Dest );
  }

  /**
   * \brief This method executes "apply(...)" and adds the result to "Dest".
   */
  void applyAdd( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::Scalar<RealType> tmp;
    apply( Arg, tmp );
    Dest += tmp;
  }
};

/**
 * \brief GreyscaleEnergyVariationWRTZeta represents an operator on a vector (or rather discretised function) \f$\zeta\f$,
 * which yields the variation of the "average shape"-energy-greyscale term for one single image wrt \f$\zeta\f$, divided by \f$H'(\zeta)\f$,
 * applied to the finite element base functions. The division of the variation by \f$H'(\zeta)\f$ is only performed to simplify the
 * equations; during a gradient descent, this variation would have to be multiplied with \f$H'(\zeta)\f$ again to give the correct gradient.
 * The computation is executed by "apply(...)".
 *
 * \author Wirth
 */
template <typename ConfigurationType,typename HeavisideFunctionType>
class GreyscaleEnergyVariationWRTZeta :
  public aol::FENonlinOpInterface<ConfigurationType, GreyscaleEnergyVariationWRTZeta<ConfigurationType,HeavisideFunctionType> > {

private:
  typedef typename ConfigurationType::RealType RealType;

  // the image $u_i$
  const aol::Vector<RealType> &_image;
  // the inverse deformation
  aol::auto_container<ConfigurationType::Dim,aol::DiscreteFunctionDefault<ConfigurationType> > _invDef;
  // the inversely deformed image $u_i(\phi_i^{-1}(x))$ and the inverse deformation/displacement $\phi_i^{-1}-id$
  aol::Vector<RealType> _invDefImageVector;
  aol::MultiVector<RealType> _invDefVector;
  const aol::DiscreteFunctionDefault<ConfigurationType> _invDefImage;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // Chan-Vese-parameters
  const RealType _lambda1, _lambda2;
  // greyscales of the image
  const RealType &_c1, &_c2;

public:
  GreyscaleEnergyVariationWRTZeta( const typename ConfigurationType::InitType &Grid,
                          const aol::Vector<RealType> &Image,
                          const aol::MultiVector<RealType> &D,
                          const HeavisideFunctionType &HeavisideFunction,
                          const RealType &C1,
                          const RealType &C2,
                          const RealType Lambda1,
                          const RealType Lambda2 ) :
    // load the grid
    aol::FENonlinOpInterface<ConfigurationType, GreyscaleEnergyVariationWRTZeta<ConfigurationType,HeavisideFunctionType> >( Grid ),
    // load the image
    _image( Image ),
    // initialise the displacement-interpolant
    _invDef(),
    // initialise the discrete inversely deformed image and the inverse deformation
    _invDefImageVector( Grid ),
    _invDefVector( Grid ),
    // generate a multilinear interpolant of the discrete inversely deformed image and the inverse deformation
    // (the vector will be filled with values later in the constructor)
    _invDefImage( Grid, _invDefImageVector ),
    // load the smoothed version of the Heaviside function to be used
    _heavisideFunction( HeavisideFunction ),
    // load the Chan-Vese parameters
    _lambda1( Lambda1 ),
    _lambda2( Lambda2 ),
    // load the greyscales
    _c1( C1 ),
    _c2( C2 ) {
    // generate the inversely deformed image $u_i(\phi_i^{-1})$
    // TransformFunction works in the way that apply first triangulates the domain grid of the deformation,
    // then computes all deformed triangles, given by three points $\phi(x_{1/2/3}$, then writes all grid
    // points y inside the triangle as $y=\sum_i\lambda_i\phi(x_i)$, and finally---assuming that $\phi$
    // and $u_i$ may be approximated as linear---computes $u_i(\phi^{-1}(y))\approx\sum_iu_i(\lambda_ix_i)$
    qc::TransformFunction<RealType,ConfigurationType::Dim> inverseDeformation( Grid );
    inverseDeformation.setDeformation( D );
    aol::MultiVector<RealType> image;
    image.appendReference( Image );
    aol::MultiVector<RealType> invDefImageId;
    invDefImageId.appendReference( _invDefImageVector );
    invDefImageId.appendReference( _invDefVector );
    // compute the identity
    aol::MultiVector<RealType> identity( ConfigurationType::Dim, Image.size() );
    qc::DataGenerator<ConfigurationType>( Grid ).generateIdentity( identity );
    image.appendReference( identity );
    inverseDeformation.apply( image, invDefImageId );
    // generate a multilinear interpolant from the discrete deformation
    for ( int j=0; j<ConfigurationType::Dim; j++ )
      _invDef.set_copy( j, aol::DiscreteFunctionDefault<ConfigurationType> ( Grid, _invDefVector[j] ) );
  }

  /**
   * \brief This method evaluates \f$(-\lambda_1(u_i(\phi_i^{-1}(x))-c_1^i)^2+\lambda_2(u_i(\phi_i^{-1}(x))-c_2^i)^2)*|det(D\phi_i^{-1})|\f$,
   * which is referred to as greyscale term. The point of evaluation is specified
   * by the grid element "El" and the quadrature point "QuadPoint" on the element or alternatively the local coordinates "LocalCoord".
   * The result is put into "NL". "DiscFunc" is not needed here.
   */
  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfigurationType>& /*DiscFunc*/,
                        const typename ConfigurationType::ElementType &El,
                        int QuadPoint,
                        const typename ConfigurationType::VecType& /*LocalCoord*/,
                        RealType &NL ) const {
    // compute $u_i(\phi_i^{-1}(x))$
    RealType invDefImageValue = _invDefImage.evaluateAtQuadPoint(El,QuadPoint);
    // compute $D\phi_i$ (the deformation gradient)
    typename ConfigurationType::VecType grad;
    typename ConfigurationType::MatType defGrad;
    for ( int j=0; j<ConfigurationType::Dim; j++ ) {
      _invDef[j].evaluateGradientAtQuadPoint(El,QuadPoint,grad);
      defGrad[j] = grad;
    }
    // add greyscale term
    NL = (-_lambda1*aol::Sqr(invDefImageValue-_c1)+_lambda2*aol::Sqr(invDefImageValue-_c2)*defGrad.det());
  }
};

/**
 * \brief EnergyVariationWRTZeta represents an operator on a vector (or rather discretised function) \f$\zeta\f$,
 * which yields the variation of the "average shape"-energy wrt \f$\zeta\f$ divided by \f$H'(\zeta)\f$, applied to the finite element base
 * functions. The division of the variation by \f$H'(\zeta)\f$ is only performed to simplify the equations; during a gradient descent,
 * this variation would have to be multiplied with \f$H'(\zeta)\f$ again to give the correct gradient.
 * The computation is executed by "apply(...)".
 *
 * \author Wirth
 */
template <typename ConfigurationType,typename HeavisideFunctionType>
class EnergyVariationWRTZeta :
  public aol::FENonlinOpInterface<ConfigurationType, EnergyVariationWRTZeta<ConfigurationType,HeavisideFunctionType> > {

private:
  typedef typename ConfigurationType::RealType RealType;

  // the images $u_i$
  const aol::MultiVector<RealType> &_images;
  const int _numberOfImages;
  // the displacements $d_i$
  const aol::auto_container<ConfigurationType::Dim, aol::MultiVector<RealType> > &_d;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // Chan-Vese-parameters
  const RealType _mu, _nu, _lambda1, _lambda2;
  // greyscales of the different images
  const aol::Vector<RealType> &_c1, &_c2;

public:
  EnergyVariationWRTZeta( const typename ConfigurationType::InitType &Grid,
                          const aol::MultiVector<RealType> &Images,
                          const aol::auto_container<ConfigurationType::Dim, aol::MultiVector<RealType> > &D,
                          const HeavisideFunctionType &HeavisideFunction,
                          const aol::Vector<RealType> &C1,
                          const aol::Vector<RealType> &C2,
                          const RealType Mu,
                          const RealType Nu,
                          const RealType Lambda1,
                          const RealType Lambda2 ) :
    // load the grid
    aol::FENonlinOpInterface<ConfigurationType, EnergyVariationWRTZeta<ConfigurationType,HeavisideFunctionType> >( Grid ),
    // load the images
    _images( Images ),
    _numberOfImages( Images.numComponents() ),
    // load the displacement
    _d( D ),
    // load the smoothed version of the Heaviside function to be used
    _heavisideFunction( HeavisideFunction ),
    // load the Chan-Vese parameters
    _mu( Mu ),
    _nu( Nu ),
    _lambda1( Lambda1 ),
    _lambda2( Lambda2 ),
    // load the greyscale vectors
    _c1( C1 ),
    _c2( C2 ) {
  }

  /**
   * \brief This method simply evaluates the area/volume term of the energy variation integrand (divided by \f$H'(\zeta(x))\f$) and thus returns \f$\nu\f$
   * The point of evaluation is specified
   * by the grid element "El" and the quadrature point "QuadPoint" on the element or alternatively the local coordinates "LocalCoord".
   * The result is put into "NL". "DiscFunc" is not needed here.
   */
  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfigurationType>& /*DiscFunc*/,
                        const typename ConfigurationType::ElementType& /*El*/,
                        int /*QuadPoint*/,
                        const typename ConfigurationType::VecType& /*LocalCoord*/,
                        RealType &NL ) const {
    //! length/surface and greyscale term are computed in "apply(...)"

    //! compute area/volume term
    NL = _nu;
  }

  /**
   * \brief This method takes a discretised function \f$\zeta\f$ as "Arg" (which is meant to be the levelset function) and computes the vector
   * "Dest"=\f$(\int_\Omega\psi_j(\mu div(\frac{\nabla\zeta}{|\nabla\zeta|})+\nu
   * +\sum_i(-\lambda_1(u_i(\phi_i^{-1}(x))-c_1^i)^2+\lambda_2(u_i(\phi_i^{-1}(x))-c_2^i)^2)*|det(D\phi_i^{-1})|)dx)_j\f$,
   * where the \f$\psi_j\f$ are the finite element base functions, and deformations \f$\phi_i\f$ (or rather the displacement), Chan-Vese
   * parameters, greyscales \f$c_{1/2}^i\f$ and images \f$u_i\f$ are member variables of the class. The first term in the integral is
   * referred to as length/surface term, the second as area/volume term, and the third as greyscale term.
   */
  void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    Dest.setZero();

    //! compute the greyscale term of the energy variation
    // for each image...
    for ( int i=0; i<_numberOfImages; i++ ){
      // reformulate the displacement
      aol::MultiVector<RealType> d;
      for ( int j=0; j < ConfigurationType::Dim; j++ )
        d.appendReference( _d[j][i] );
      GreyscaleEnergyVariationWRTZeta<ConfigurationType,HeavisideFunctionType> greyscaleEnergyVariation( this->_initializer, _images[i], d, _heavisideFunction, _c1[i], _c2[i], _lambda1, _lambda2 );
      greyscaleEnergyVariation.applyAdd( Arg, Dest );
    }
    Dest /= (_numberOfImages*_mu);

    //! compute the length/surface term of the energy variation
    qc::MCMStiffOp<ConfigurationType> mcmStiffOp( this->_initializer, aol::ONTHEFLY, 0.1 );
    mcmStiffOp.setImageReference( Arg );
    mcmStiffOp.applyAdd( Arg, Dest );
    Dest *= _mu;

    //! compute the area/volume term of the energy variation
    // careful! apply of the parent class would not work, since it would call applyAdd of this class!
    aol::FENonlinOpInterface<ConfigurationType, EnergyVariationWRTZeta<ConfigurationType,HeavisideFunctionType> >::applyAdd( Arg, Dest );
  }

  /**
   * \brief This method executes "apply(...)" and adds the result to "Dest".
   */
  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp( Arg, aol::STRUCT_COPY );
    apply( Arg, tmp );
    Dest += tmp;
  }
};

/**
 * \brief This class performs a gradient descent on the levelset function \f$\zeta\f$ (which has to be passed to the
 * "apply(...)"-method), using the energy and the energy variation divided by \f$H'(\zeta(x))\f$, both passed
 * to the constructor. The metric, used to define the gradient, is provided by the method "smoothDirection(...)".
 *
 * \author Wirth
 */
template <typename ConfigurationType, typename HeavisideFunctionType>
class GradientDescentZeta :
  public aol::GradientDescentWithAutomaticFilterWidth<ConfigurationType, aol::Vector<typename ConfigurationType::RealType> > {

private:
  // Heaviside function used in the metric (which is needed to define a gradient)
  const HeavisideFunctionType &_heavisideFunction;

public:
  typedef typename ConfigurationType::RealType RealType;

  GradientDescentZeta( // grid on which the discretised objective function $\zeta$ is defined
                       const typename ConfigurationType::InitType &Grid,
                       // the total energy as an operator on $\zeta$
                       const aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > &E,
                       // the variation of the energy wrt $\zeta$ as an operator on $\zeta$
                       // for efficiency of code, "DE" represents not the correct energy variation,
                       // but the energy variation divided by $H'(\zeta)$. This is remedied by
                       // using a special metric when computing the gradient.
                       const aol::Op<aol::Vector<RealType> > &DE,
                       // the smoothed version of the Heaviside function to be used in the metric
                       const HeavisideFunctionType &HeavisideFunction,
                       // maximum number of gradient descent steps
                       const int MaxIterations = 50,
                       // numerical parameters
                       const RealType StartTau = aol::ZOTrait<RealType>::one,
                       const RealType StopEpsilon = aol::ZOTrait<RealType>::zero) :
    // load the passed parameters
    aol::GradientDescentWithAutomaticFilterWidth<ConfigurationType, aol::Vector<RealType> >( Grid, E, DE, MaxIterations, StartTau, StopEpsilon ),
    // create the specified Heaviside function
    _heavisideFunction( HeavisideFunction ) {
  }

  /**
   * \brief This method is automatically executed after each step of gradient descent and could for example
   * be used to save intermediate results (which is at the moment not implemented here).
   * The argument "Position" represents the current state of the discretised objective function,
   * i.e. the current \f$\zeta\f$ of the gradient descent. The argument "Iteration" is the number of
   * the current gradient descent step.
   * If no code is to be performed after the current gradient descent step, "false" is returned,
   * else "true".
   * In order to render the gradient descent more stable and to render \f$\zeta\f$ unique, \f$\zeta\f$
   * is by this method transformed into a signed distance function each fifth iteration.
   */
  virtual bool postProcess( aol::Vector<RealType> &Position, const int Iteration ) const {
    // each fifth iteration do ...
    if ( Iteration % 5 == 0 ) {
      // transform "Position" (i.e. $\zeta$) into a signed distance function
      if ( ConfigurationType::Dim == 2 ) {
        qc::SignedDistanceOp<ConfigurationType> signedDistOp( this->_grid );
        signedDistOp.apply( Position, Position );
      } else {
        qc::SignedDistanceOp3D<RealType> signedDistOp( this->_grid );
        signedDistOp.apply( Position, Position );
      }
      return true;
    } else
      // return false, since no code was executed
      return false;
  }

  /**
   * \brief This method implements the metric (i.e. the dot product) which is used to define the gradient from the energy variation.
   * Assume, the dot product is given by \f$(u,Av)=(Au,v)\f$ for a symmetric linear positive definite operator \f$A\f$. Then
   * the method transforms the vector (or rather discretised function) "Direction" into \f$A\f$ applied on "Direction". \f$A\f$ may depend
   * on "Position", which is the current state of the discretised objective function, i.e. the current \f$\zeta\f$ of the gradient descent.
   * "Sigma" simply is a parameter which may be used by the dot product. Depending on the computation mode, "GradientDescent::apply(...)"
   * (the method which performs the gradient descent) sets "Sigma" to one or decreases it towards zero each step.
   *
   * Here, the dot product also contains a pointwise division by \f$H'(\zeta(x))\f$, since the energy variation was also divided by
   * \f$H'(\zeta(x))\f$. Hence, when defining the gradient with the help of the dot product, both are made consistent, since gradient
   * descent can be expressed as "\f$\partial_{t}\zeta=-\partial_{\zeta}E\f$" which yields the following consistent, weak (time-discretised)
   * formulation for a testfunction \f$v\f$:
   * \f$g(\zeta_{k+1},v/H'(\zeta))=g(\zeta_k,v/H'(\zeta))-\tau<\partial_{\zeta}E/H'(\zeta),v>\f$, where \f$g\f$ is the used dot product.
   * The implemented dot product uses an operator \f$A\f$ which represents the concatenation of the above mentioned pointwise division by
   * \f$H'(\zeta(x))\f$ and a smoothing operator, that is \f$A=B\circ C\f$ for two operators \f$B\f$ and \f$C\f$ as follows:
   *
   * \f$C\f$ is the operator such that \f$(u,Cv)\f$ corresponds to \f$\int_{\Omega}u(x)v(x)/H'(\zeta(x))dx\f$, only that all row entries of this
   * special mass operator are summed up and put into the diagonal element (remember that \f$C\f$ represents a matrix in the discretised
   * case). This yields a lumped mass operator.
   *
   * The smoothing operator \f$B\f$ makes the gradient descent faster (by distributing local information globally and thus enhancing flow of
   * information). \f$B\f$ smoothes most for large "Sigma"=:\f$\sigma\f$ and reduces to the normal mass operator for \f$\sigma=0\f$. \f$B\f$ is the
   * operator corresponding to the dot product \f$\int_{\Omega}uv+\sigma/2 \nabla u \cdot \nabla v dx\f$ of \f$u\f$ and \f$v\f$.
   */
  virtual void smoothDirection( const aol::Vector<RealType> &Position, const RealType Sigma, aol::Vector<RealType> &Direction ) const {
    // apply $C$ on "Direction" and put the result into "aux"
    aol::Vector<RealType> aux( this->_grid.getNumberOfNodes() );
    aol::HeavisideFunctionLumpedMassOp<ConfigurationType, HeavisideFunctionType>
      heavisideFunctionLumpedMassInvOp( this->_grid, Position, _heavisideFunction, aol::INVERT );
    heavisideFunctionLumpedMassInvOp.apply( Direction, aux );

    // apply $B$ on "aux" and put the result into "Direction"
    qc::LinearSmoothOp<RealType> linSmooth;
    linSmooth.setCurrentGrid( this->_grid );
    linSmooth.setSigma( Sigma );
    linSmooth.apply( aux, Direction );
  }

  /**
   * \brief This method is automatically executed after each step of gradient descent and could for example
   * be used to save intermediate results (which is at the moment not implemented here).
   */
  virtual void writeTimeStep( const aol::Vector<RealType>& /*Dest*/, const int /*Iterations*/ ) const {
    /*
    // uses that GradientDescent is a subclass of QuocTimestepSaver
    this->setSaveDirectory( "intermediateData" );
    this->setSaveName( "intermediateZeta" );
    this->saveTimestepPNGDrawRedAndBlueIsoline( 1, Iterations, 0., 0.1, Dest, _image, this->_grid );
    */
  }
};

/**
 * \brief This class represents the operator which computes the following energy as a function of the displacement \f$d_i\f$:
 * \f$E=\int_\Omega\lambda_1(1-H(\zeta(\phi_i(x))))(u_i(x)-c_1^i)^2+\lambda_2H(\zeta(\phi_i(x)))(u_i(x)-c_2^i)^2dx
 * \int_\Omega W(||D\phi_i||,||Cof(D\phi_i)||,det(D\phi_i))dx\f$.
 * \f$\zeta\f$ represents the levelset function, \f$H\f$ the Heaviside function, \f$u_i\f$ the \f$i\f$-th image, and \f$\phi_i\f$ the \f$i\f$-th deformation.
 * The deformation is the displacement plus the identity.
 * The evaluation of the integral is performed by the method "apply(...)" or "applyAdd(...)", which uses the method
 * "evaluateIntegrand(...)" to evaluate the integrand at a point. The first term of the integrand can be referred to as
 * length/surface term, the second term as area/volume term, and the sums as greyscale term.
 * Since the energy shall be an operator on \f$d_i\f$, \f$H\f$, \f$\Omega\f$, \f$u_i\f$ and \f$\zeta\f$ have to be passed to the constructor.
 *
 * \author Wirth
 */
template <typename ConfigurationType, typename HeavisideFunctionType>
class EnergyOfD : public aol::Op<aol::MultiVector<typename ConfigurationType::RealType>,aol::Scalar<typename ConfigurationType::RealType> > {

private:
  typedef typename ConfigurationType::RealType RealType;

  // the domain
  const typename ConfigurationType::InitType _grid;
  // the image $u_i$
  const aol::Vector<RealType> &_image;
  // the levelset function $\zeta$
  const aol::Vector<RealType> &_zeta;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // Chan-Vese-parameters
  const RealType _lambda1, _lambda2;
  // greyscales of the image
  const RealType _c1, _c2;
  // weights of the hyperelastic deformation energy
  const RealType _volumeEnergyWeight;
  const RealType _lengthEnergyWeight;

public:
  EnergyOfD( const typename ConfigurationType::InitType &Grid,
             const aol::Vector<RealType> &Image,
             const aol::Vector<RealType> &Zeta,
             const HeavisideFunctionType &HeavisideFunction,
             const RealType C1,
             const RealType C2,
             const RealType Lambda1,
             const RealType Lambda2,
             const RealType VolumeEnergyWeight,
             const RealType LengthEnergyWeight ) :
     // load the grid
    _grid( Grid ),
    // load the image
    _image( Image ),
    // load the levelset function
    _zeta( Zeta ),
    // load the smoothed version of the Heaviside function to be used
    _heavisideFunction( HeavisideFunction ),
    // load the Chan-Vese parameters
    _lambda1( Lambda1 ),
    _lambda2( Lambda2 ),
    // load the greyscale vectors
    _c1( C1 ),
    _c2( C2 ),
    // load the weights of the hyperelastic deformation energy
    _volumeEnergyWeight( VolumeEnergyWeight ),
    _lengthEnergyWeight( LengthEnergyWeight ) {
  }

  /**
   * \brief This method computes the "average shape"-energy terms in which displacement \f$\phi_i\f$ occurs.
   * \f$\phi_i\f$ (discretised) is passed via "Arg". The result is put into "Dest".
   * The energy is the sum of a hyperelastic term and a greyscale term.
   */
  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    //! compute the hyperelastic energy
    /*
    Dest /= _volumeEnergyWeight;
    qc::VolumeEnergy<ConfigurationType,ConfigurationType::Dim> volumeEnergy( _grid );
    volumeEnergy.applyAdd( Arg, Dest );
    Dest *= _volumeEnergyWeight/_lengthEnergyWeight;  // the factors are only meant to decrease the deformation energy, if needed
    qc::lengthEnergy<ConfigurationType,ConfigurationType::Dim> lengthEnergy( _grid );
    lengthEnergy.applyAdd( Arg, Dest );
    Dest *= _lengthEnergyWeight;
    */
    qc::HyperelasticEnergyDensityDefault<ConfigurationType> hyperelasticEnergyDensity( _lengthEnergyWeight, 0, _volumeEnergyWeight );
    qc::HyperelasticEnergy<ConfigurationType,qc::HyperelasticEnergyDensityDefault<ConfigurationType> > hyperelasticEnergy( _grid, hyperelasticEnergyDensity );
    hyperelasticEnergy.applyAdd( Arg, Dest );

    //! compute the greyscale term of the energy
    GreyscaleEnergyOfZeta<ConfigurationType,HeavisideFunctionType> greyscaleEnergyOfZeta( _grid, _image, Arg, _heavisideFunction, _c1, _c2, _lambda1, _lambda2 );
    greyscaleEnergyOfZeta.applyAdd( _zeta, Dest );
  }
};

/**
 * \brief EnergyVariationWRTD represents an operator on a multi-vector (or rather discretised, multidimensional function) \f$d_i\f$,
 * which yields the variation of the "average shape"-energy wrt \f$d_i\f$, applied to the multidimensional finite element base
 * functions.
 * The computation is executed by "apply(...)".
 *
 * \author Wirth
 */
template <typename ConfigurationType,typename HeavisideFunctionType>
class EnergyVariationWRTD :
  public aol::FENonlinVectorOpInterface<ConfigurationType, ConfigurationType::Dim, ConfigurationType::Dim, EnergyVariationWRTD<ConfigurationType,HeavisideFunctionType> > {

private:
  typedef typename ConfigurationType::RealType RealType;

  // the image $u_i$
  const aol::Vector<RealType> &_image;
  aol::DiscreteFunctionDefault<ConfigurationType> _imageFunction;
  // the levelset function $\zeta$
  const aol::Vector<RealType> &_zeta;
  aol::DiscreteFunctionDefault<ConfigurationType> _levelsetFunction;
  // the Heaviside function
  const HeavisideFunctionType &_heavisideFunction;
  // Chan-Vese-parameters
  const RealType _lambda1, _lambda2;
  // greyscales of the image
  const RealType _c1, _c2;
  // weights of the hyperelastic deformation energy
  const RealType _volumeEnergyWeight;
  const RealType _lengthEnergyWeight;

public:
  EnergyVariationWRTD( const typename ConfigurationType::InitType &Grid,
                       const aol::Vector<RealType> &Image,
                       const aol::Vector<RealType> &Zeta,
                       const HeavisideFunctionType &HeavisideFunction,
                       const RealType C1,
                       const RealType C2,
                       const RealType Lambda1,
                       const RealType Lambda2,
                       const RealType VolumeEnergyWeight,
                       const RealType LengthEnergyWeight ) :
    // load the grid
    aol::FENonlinVectorOpInterface<ConfigurationType, ConfigurationType::Dim, ConfigurationType::Dim, EnergyVariationWRTD<ConfigurationType,HeavisideFunctionType> >( Grid ),
    // load the image
    _image( Image ),
    // generate the image-interpolant
    _imageFunction( Grid, Image ),
    // load the levelset function
    _zeta( Zeta ),
    // generate the levelset function interpolant
    _levelsetFunction( Grid, Zeta ),
    // load the smoothed version of the Heaviside function to be used
    _heavisideFunction( HeavisideFunction ),
    // load the Chan-Vese parameters
    _lambda1( Lambda1 ),
    _lambda2( Lambda2 ),
    // load the greyscale vectors
    _c1( C1 ),
    _c2( C2 ),
    // load the weights of the hyperelastic deformation energy
    _volumeEnergyWeight( VolumeEnergyWeight ),
    _lengthEnergyWeight( LengthEnergyWeight ) {
  }

  /**
   * \brief This method computes the vector \f$H'(\zeta(\phi_i(x)))\nabla\zeta(\phi_i(x))(-\lambda_1(u_i(x)-c_1^i)^2+\lambda_2(u_i(x)-c_2^i)^2)\f$
   * for some \f$i\f$, where \f$\phi_i\f$ is passed as "DiscFuncs", \f$x\f$ is specified by "El", "QuadPoint", and "LocalCoord", and the result
   * is put into "NL".
   */
  void getNonlinearity( aol::auto_container<ConfigurationType::Dim,aol::DiscreteFunctionDefault<ConfigurationType> > &DiscFuncs,
                        const typename ConfigurationType::ElementType &El,
                        int QuadPoint,
                        const typename ConfigurationType::VecType &LocalCoord,
                        aol::Vec<ConfigurationType::Dim,typename ConfigurationType::RealType> &NL ) const {
    //! compute greyscale term
    // compute $H'(\zeta(\phi_i(x)))\nabla\zeta(\phi_i(x))$
    bool displacedOutsideDomain = false;
    aol::Vec<ConfigurationType::Dim,typename ConfigurationType::RealType> DHZP;
    DHZP = EvaluationFunctions<ConfigurationType, HeavisideFunctionType>::evaluateDHZetaPhi( _heavisideFunction, _levelsetFunction, DiscFuncs, this->_initializer, El, QuadPoint, LocalCoord, displacedOutsideDomain );
    // add greyscale term
    if ( displacedOutsideDomain )
      NL.setZero();
    else
      NL = (-_lambda1*aol::Sqr(_imageFunction.evaluateAtQuadPoint( El, QuadPoint )-_c1)
            +_lambda2*aol::Sqr(_imageFunction.evaluateAtQuadPoint( El, QuadPoint )-_c2))*DHZP;
  }

  /**
   * \brief This method computes the MultiVector \f$(\int_\Omega H'(\zeta(\phi_i(x)))\nabla\zeta(\phi_i(x))(-\lambda_1(u_i(x)-c_1^i)^2+\lambda_2(u_i(x)-c_2^i)^2)\psi_j dx)_j\f$ plus
   * the hyperelastic energy variation, where the \f$\psi_j\f$ are multidimensional finite element base functions (hence the MultiVector).
   * The deformation \f$\phi_i\f$ is the displacement \f$d\f$, which is passed as (the discretized) "Arg", the remaining parameters and
   * functions are member variables of the class and thus passed by the constructor. The result is put into "Dest".
   */
  void apply( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    //! compute the hyperelastic energy variation
    /*qc::VolumeGradient<ConfigurationType,ConfigurationType::Dim> volumeEnergy( this->_initializer );
    volumeEnergy.apply( Arg, Dest );
    Dest *= _volumeEnergyWeight/_lengthEnergyWeight;  // the factors are only meant to decrease the deformation energy, if needed
    qc::lengthGradient<ConfigurationType,ConfigurationType::Dim> lengthEnergy( this->_initializer );
    lengthEnergy.applyAdd( Arg, Dest );
    Dest *= _lengthEnergyWeight;
    */
    qc::HyperelasticEnergyDensityDefault<ConfigurationType> hyperelasticEnergyDensity( _lengthEnergyWeight, 0, _volumeEnergyWeight );
    qc::HyperelasticGradient<ConfigurationType,qc::HyperelasticEnergyDensityDefault<ConfigurationType> > hyperelasticGradient( this->_initializer, hyperelasticEnergyDensity );
    hyperelasticGradient.apply( Arg, Dest );

    //! compute the greyscale term of the energy variation
    // careful! apply of the parent class would not work, since this would call applyAdd of this class!
    aol::FENonlinVectorOpInterface<ConfigurationType, ConfigurationType::Dim, ConfigurationType::Dim, EnergyVariationWRTD<ConfigurationType,HeavisideFunctionType> >::applyAdd( Arg, Dest );
  }

  /**
   * \brief This method executes "apply(...)" and adds the result to "Dest".
   */
  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> tmp( Arg, aol::STRUCT_COPY );
    apply( Arg, tmp );
    Dest += tmp;
  }
};

/**
 * \brief This class performs a gradient descent on the displacement \f$d_i\f$ (which has to be passed to the
 * "apply(...)"-method), using the energy and the energy variation, both passed
 * to the constructor. The metric, used to define the gradient, is provided by the method "smoothDirection(...)".
 *
 * \author Wirth
 */
template <typename ConfigurationType, typename HeavisideFunctionType>
class GradientDescentD :
  public aol::GradientDescentWithAutomaticFilterWidth<ConfigurationType, aol::MultiVector<typename ConfigurationType::RealType> > {

private:
  // Heaviside function used in the metric (which is needed to define a gradient)
  const HeavisideFunctionType &_heavisideFunction;

public:
  typedef typename ConfigurationType::RealType RealType;

  GradientDescentD( // grid on which the discretised objective function $\zeta$ is defined
                       const typename ConfigurationType::InitType &Grid,
                       // the total energy as an operator on $\zeta$
                       const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                       // the variation of the energy wrt $\zeta$ as an operator on $\zeta$
                       // for efficiency of code, "DE" represents not the correct energy variation,
                       // but the energy variation divided by $H'(\zeta)$. This is remedied by
                       // using a special metric when computing the gradient.
                       const aol::Op<aol::MultiVector<RealType> > &DE,
                       // the smoothed version of the Heaviside function to be used in the metric
                       const HeavisideFunctionType &HeavisideFunction,
                       // maximum number of gradient descent steps
                       const int MaxIterations = 50,
                       // numerical parameters
                       const RealType StartTau = aol::ZOTrait<RealType>::one,
                       const RealType StopEpsilon = aol::ZOTrait<RealType>::zero) :
    // load the passed parameters
    aol::GradientDescentWithAutomaticFilterWidth<ConfigurationType, aol::MultiVector<RealType> >( Grid, E, DE, MaxIterations, StartTau, StopEpsilon ),
    // create the specified Heaviside function
    _heavisideFunction( HeavisideFunction ) {
  }

  virtual bool postProcess( aol::MultiVector<RealType> &/*Position*/, const int /*Iteration*/ ) const {
      return false;
  }

  /**
   * \brief This method implements the metric (i.e. the dot product) which is used to define the gradient from the energy variation.
   * Assume, the dot product is given by \f$(u,Av)=(Au,v)\f$ for a symmetric linear positive definite operator \f$A\f$. Then
   * the method transforms the vector (or rather discretised function) "Direction" into \f$A\f$ applied on "Direction". \f$A\f$ may depend
   * on "Position", which is the current state of the discretised objective function, i.e. the current \f$d_i\f$ of the gradient descent.
   * "Sigma" simply is a parameter which may be used by the dot product. Depending on the computation mode, "GradientDescent::apply(...)"
   * (the method which performs the gradient descent) sets "Sigma" to one or decreases it towards zero each step.
   *
   * The smoothing operator \f$A\f$ makes the gradient descent faster (by distributing local information globally and thus enhancing flow of
   * information). \f$A\f$ smoothes most for large "Sigma"=:\f$\sigma\f$ and reduces to the normal mass operator for \f$\sigma=0\f$. \f$A\f$ is the
   * operator corresponding to the dot product \f$\int_{\Omega}u^Tv+\sigma/2 \nabla u:\nabla v dx\f$ of \f$u\f$ and \f$v\f$.
   */
  virtual void smoothDirection( const aol::MultiVector<RealType> &/*Position*/, const RealType Sigma, aol::MultiVector<RealType> &Direction ) const {
    // apply $A$ on "Direction" and put the result into "Direction"
    qc::LinearSmoothOp<RealType> linSmooth;
    linSmooth.setCurrentGrid( this->_grid );
    linSmooth.setSigma( Sigma );
    for ( int i=0; i<Direction.numComponents(); i++ )
      // in "linSmooth", Arg and Dest may be the same!
      linSmooth.apply( Direction[i], Direction[i] );
  }

  /**
   * \brief This method is automatically executed after each step of gradient descent and could for example
   * be used to save intermediate results (which is at the moment not implemented here).
   */
  virtual void writeTimeStep( const aol::MultiVector<RealType> &/*Dest*/, const int /*Iterations*/ ) const {
    /*
    // uses that GradientDescent is a subclass of QuocTimestepSaver
    this->setSaveDirectory( "intermediateData" );
    this->setSaveName( "intermediateZeta" );
    this->saveTimestepPNGDrawRedAndBlueIsoline( 1, Iterations, 0., 0.1, Dest, _image, this->_grid );
    */
  }

};

/**
 * \brief This class represents the Chan-Vese segmentation and the search for an average shape. It is
 * constructed from a parser, which reads in a file containing parameters. The average shape is
 * computed in the method "execute()". Results are saved as specified in the parameter file.
 * Is currently implemented for white objects on black background.
 *
 * \author Wirth
 */
template <typename ConfigurationType>
class ChanVeseAverageShapeFinder {

private:
  typedef typename ConfigurationType::RealType RealType;
  typedef typename ConfigurationType::ArrayType ArrayType;

  // grid over the examined domain
  const typename ConfigurationType::InitType _grid;
  // number of images, the images themselves, the directory and file name word stem for the results
  const int _numberOfImages;
  aol::MultiVector<RealType> _images;
  char _destDirectory[1024], _destFileNameStem[1024], _fileExtension[1024];
  // Chan-Vese-parameters
  const RealType _mu, _nu, _lambda1, _lambda2;
  // weights of the hyperelastic deformation energy
  const RealType _volumeEnergyWeight;
  const RealType _lengthEnergyWeight;
  // sharpness of Heaviside function
  const RealType _delta;
  // iteration numbers
  const int _outerIterations, _innerIterations;
  // greyscales of the different images
  aol::Vector<RealType> _c1, _c2;
  // level set function (\zeta)
  aol::Vector<RealType> _zeta;
  // displacements ($d_{i}$): the j-th component of _phi represents a multivector,
  // containing the j-th displacement component for each image
  // the deformation is then given by $\Phi=$identity$+d$
  aol::auto_container<ConfigurationType::Dim, aol::MultiVector<RealType> > _d;
  // coarsest level for the multiscale method
  const int _coarsestScale;

  /**
   * \brief Returns Value, casted as RealType.
   */
  inline RealType castDouble( double Value ) {
    return static_cast<RealType>(Value);
  }

  /**
   * \brief Returns total energy.
   */
  inline aol::Scalar<RealType> totalEnergy( const qc::GridDefinition &Grid,
                                       const aol::MultiVector<RealType> &Images,
                                       const aol::Vector<RealType> &Zeta,
                                       const aol::ArcTanHeavisideFunction<RealType> &HeavisideFunction,
                                       const aol::auto_container<ConfigurationType::Dim,aol::MultiVector<RealType> > &D,
                                       const aol::Vector<RealType> &C1,
                                       const aol::Vector<RealType> &C2,
                                       const RealType VolumeEnergyWeight,
                                       const RealType LengthEnergyWeight ) {
    aol::Scalar<RealType> energy;
    EnergyOfZeta<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> >
      EZ( Grid, Images, D, HeavisideFunction, C1, C2, _mu, _nu, _lambda1, _lambda2 );
    EZ.apply( Zeta, energy );
    energy *= _numberOfImages;
    for ( int i=0; i < _numberOfImages; i++ ){
      EnergyOfD<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> >
        ED( Grid, Images[i], Zeta, HeavisideFunction, C1[i], C2[i], 0, 0, VolumeEnergyWeight, LengthEnergyWeight );
      aol::MultiVector<RealType> disp;
      for ( int j=0; j < ConfigurationType::Dim; j++ )
        disp.appendReference( D[j][i] );
      ED.applyAdd( disp, energy );
    }
    // split up energies if wanted
    aol::Scalar<RealType> energyComp;
    EnergyOfZeta<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> >
      E1( Grid, Images, D, HeavisideFunction, C1, C2, _mu, _nu, 0, 0 );
    E1.apply( Zeta, energyComp );
    cerr << "Chan-Vese length/area energy: " << energyComp << endl;
    for ( int i=0; i < _numberOfImages; i++ ){
      // make displacement of image i a MultiVector and an auto_container as well as the image itself
      aol::MultiVector<RealType> disp;
      aol::auto_container<ConfigurationType::Dim,aol::MultiVector<RealType> > dispAutoCont;
      for ( int j=0; j < ConfigurationType::Dim; j++ ) {
        disp.appendReference( D[j][i] );
        aol::MultiVector<RealType> dispComp;
        dispComp.appendReference( disp[j] );
        dispAutoCont.set_copy( j, dispComp );
      }
      aol::MultiVector<RealType> image;
      image.appendReference( Images[i] );
      // compute greyscale engergy of the ith image
      aol::Vector<RealType> c1(1), c2(1);
      c1[0] = C1[i];
      c2[0] = C2[i];
      EnergyOfZeta<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> >
        E2( Grid, image, dispAutoCont, HeavisideFunction, c1,  c2, 0, 0, _lambda1, _lambda2 );
      E2.apply( Zeta, energyComp );
      cerr << "greyscales of image " << i << ": " << C1[i] << ", " << C2[i] << endl;
      cerr << "Chan-Vese greyscale energy of image " << i << ": " << energyComp/_numberOfImages << endl;
      // compute elastic engergy of the ith image
      EnergyOfD<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> >
        ED( Grid, Images[i], Zeta, HeavisideFunction, C1[i], C2[i], 0, 0, VolumeEnergyWeight, LengthEnergyWeight );
      ED.apply( disp, energyComp );
      cerr << "Elastic energy of deformation " << i << ": " << energyComp/_numberOfImages << endl;
    }
    return energy/_numberOfImages;
  }

public:
  ChanVeseAverageShapeFinder( aol::ParameterParser &Parser ) :
    // define the grid over the domain
    _grid( Parser.getInt( "log2GridWidth" ), ConfigurationType::Dim ),
    // load number of images and create the image-vector
    _numberOfImages( Parser.getInt( "numberOfImages" ) ),
    _images( _numberOfImages, _grid.getNumberOfNodes() ),
    // load Chan-Vese- and hyperelastic parameters
    _mu( castDouble( Parser.getDouble( "mu" ) ) ),
    _nu( castDouble( Parser.getDouble( "nu" ) ) ),
    _lambda1( castDouble( Parser.getDouble( "lambda1" ) ) ),
    _lambda2( castDouble( Parser.getDouble( "lambda2" ) ) ),
    _volumeEnergyWeight( castDouble( Parser.getDouble( "volumeEnergyWeight" ) ) ),
    _lengthEnergyWeight( castDouble( Parser.getDouble( "lengthEnergyWeight" ) ) ),
    _delta( castDouble( Parser.getDouble( "delta" ) ) ),
    // load number of outer and inner algorithm iterations
    _outerIterations( Parser.getInt( "outerIterations" ) ),
    _innerIterations( Parser.getInt( "innerIterations" ) ),
    // create greyscale vectors
    _c1(_numberOfImages),
    _c2(_numberOfImages),
    // create the level set function and the deformations
    _zeta( _grid ),
    _coarsestScale( Parser.getInt( "log2CoarseGridWidth" ) ) {

    // load word stem of image file names
    char fileNameWordStem[1024];
    Parser.getString( "imageFileNameStem", fileNameWordStem );
    Parser.getString( "fileExtension", _fileExtension );

    // load the images themselves and intensify their contrast
    for ( int i = 1; i <= _numberOfImages; i++ ) {
      // load i-th image
      char fileName[1024];
      sprintf( fileName, "./%s%d.%s", fileNameWordStem, i, _fileExtension );
      ArrayType tmpArray( _images[i-1], _grid );
      tmpArray.load( fileName );
      // intensify the contrast
      tmpArray.addToAll( -tmpArray.getMinValue() );
      tmpArray /= tmpArray.getMaxValue();
    }

    // read in the directory and file name word stem, where results are to be saved
    Parser.getString( "destDirectory", _destDirectory );
    Parser.getString( "resultFileNameStem", _destFileNameStem );

    // initialise the level set function
    /*qc::DataGenerator<ConfigurationType> levelSetFunctionInitialiser( _grid );
    levelSetFunctionInitialiser.generateCircleLevelset( _zeta, castDouble( Parser.getDouble( "initialLevelsetRadius")), .5, .5 );
    */
    typename ConfigurationType::ArrayType aux( _zeta, _grid, aol::FLAT_COPY );
    LevelsetGenerator<ConfigurationType>::generateBallLevelset( aux, castDouble( Parser.getDouble( "initialLevelsetRadius")), .5, .5, .5 );

    // initialise the displacements
    for ( int i = 0; i < ConfigurationType::Dim; i++ ) {
      aol::MultiVector<RealType> dComponent( _numberOfImages, _grid.getNumberOfNodes() );
      // dComponent.setZero(); // MultiVector is automatically initialised with 0
      _d.set_copy ( i, dComponent );
    }

    // notify the user
    cerr << "Parameters loaded." << endl;
  }

  void execute() {
    //! produce a single image with all images overlaid
    ArrayType overlaidImage( _grid );
    for ( int i=0; i<_numberOfImages; i++ )
      overlaidImage += _images[i];
    // intensify the contrast
    overlaidImage.addToAll( -overlaidImage.getMinValue() );
    overlaidImage /= overlaidImage.getMaxValue();

    //! create a class to save the (intermediate) results
    qc::QuocTimestepSaver<ConfigurationType> timeStepSaver( 1, 0, _destFileNameStem, true );
    timeStepSaver.setSaveDirectory( _destDirectory );

    //! perform a pre-registration/matching, initialising the displacements correctly
    // for each image shift the centre of gravity to the middle of the image
    // currently implemented for white images on black background
    MomentIntegrator<ConfigurationType> momentIntegrator( _grid );
    aol::Vector<RealType> momentResult( ConfigurationType::Dim+1 );
    for ( int i=0; i<_numberOfImages; i++ ) {
      momentIntegrator.apply( _images[i], momentResult );
      // compute the centre of gravity and define a corresponding shifting displacement
      for ( int j=0; j<ConfigurationType::Dim; j++ )
        _d[j][i].setAll( 0.5-momentResult[j]/momentResult[ConfigurationType::Dim] );
    }
    /*// upscale each image so that it will shrink to the correct shape
    qc::DataGenerator<ConfigurationType> upScaler( _grid );
    aol::Vector<RealType> scaleResult( _grid );
    for ( int j=0; j<ConfigurationType::Dim; j++ ) {
      upScaler.generateLineLevelset( scaleResult, 1-j, .5 );
      scaleResult *= .1;
      for ( int i=0; i<_numberOfImages; i++ )
        _d[j][i] += scaleResult;
    }
    */
    /*// shear second image in x-direction
    // flat copy of x-component of displacement of second image
    qc::ScalarArray<RealType, qc::QC_2D> d2x( _d[0][1], _grid );
    for ( int x=0; x<_grid.getWidth(); x++ )
      for ( int y=0; y<_grid.getWidth(); y++ )
        d2x.add( x, y, -.8*((_grid.getWidth()+1)/2-y)*_grid.H() );
    */

    //! generate the different levels of the images, level set function and displacements
    // define the corresponding multilevel arrays
    qc::MultilevelArray<RealType> zeta( _grid );
    qc::MultiDimMultilevelArray<RealType> images( _grid, _numberOfImages );
    qc::MultiDimMultilevelArray<RealType> d( _grid, _numberOfImages*ConfigurationType::Dim ); // contains all displacements in x-direction, then in y-direction,...
    qc::MultilevelArray<RealType> scaledOverlaidImage( _grid );
    // initialise the multilevel arrays (i.e. pass the initial values)
    typename ConfigurationType::ArrayType auxZeta( zeta.current(), aol::FLAT_COPY );
    auxZeta = _zeta;
    typename ConfigurationType::ArrayType auxOverlaidImage( scaledOverlaidImage.current(), aol::FLAT_COPY );
    auxOverlaidImage = overlaidImage;
    for ( int i=0; i<_numberOfImages; i++ ) {
      typename ConfigurationType::ArrayType auxU( images.getArray( i ), aol::FLAT_COPY );
      auxU = _images[i];
      for ( int j=0; j<ConfigurationType::Dim; j++ ) {
        typename ConfigurationType::ArrayType auxD( d.getArray( j*_numberOfImages+i ), aol::FLAT_COPY );
        auxD = _d[j][i];
      }
    }
    // generate the needed coarser scales
    zeta.levRestrict( _coarsestScale, _grid.getGridDepth() );
    images.levRestrict( _coarsestScale, _grid.getGridDepth() );
    d.levRestrict( _coarsestScale, _grid.getGridDepth() );
    scaledOverlaidImage.levRestrict( _coarsestScale, _grid.getGridDepth() );

    //! on each scale of the multiscale method do...
    for ( int scale = _coarsestScale; scale <= _grid.getGridDepth(); scale++ ){
      // refine mesh
      typename ConfigurationType::InitType grid( scale, ConfigurationType::Dim );
      zeta.setCurLevel( scale );
      images.setCurLevel( scale );
      d.setCurLevel( scale );
      scaledOverlaidImage.setCurLevel( scale );

      // define the image, levelset function and displacement on the current scale
      aol::Vector<RealType> zetaCur( zeta.current(), aol::FLAT_COPY );
      aol::MultiVector<RealType> imagesCur;
      for ( int i=0; i<_numberOfImages; i++ )
        imagesCur.appendReference( images.getArray( i ) );
      aol::auto_container<ConfigurationType::Dim,aol::MultiVector<RealType> > dCur;
      for ( int j=0; j<ConfigurationType::Dim; j++ ) {
        dCur.set_copy( j, aol::MultiVector<RealType>(0,0) );
        for ( int i=0; i<_numberOfImages; i++ )
          dCur[j].appendReference( d.getArray( j*_numberOfImages+i ) );
      }

      // set initial delta for this scale (should be relatively small on coarse scale, but enlarging it strongly from level to level destroys the results)
      //RealType levelFactor = (_grid.getGridDepth()/scale);
      RealType delta = _delta;

      // set weights of hyperelastic deformation energy for this scale (should be smaller on coarse scale)
      RealType volumeEnergyWeight = _volumeEnergyWeight;// /aol::Sqr(aol::Sqr(aol::Sqr(levelFactor)));
      RealType lengthEnergyWeight = _lengthEnergyWeight;// /aol::Sqr(levelFactor);

      //! each outer iteration, delta is halved, making the Heaviside function sharper
      for ( int outerIteration = 0; outerIteration < _outerIterations; outerIteration++ ){
        // create the appropriate Heaviside function
        aol::ArcTanHeavisideFunction<RealType> heavisideFunction( delta );

        //! each inner iteration represents one run of gradient descent on c1, c2, zeta, and d
        for ( int innerIteration = 0; innerIteration < _innerIterations; innerIteration++ ){
          cerr << endl << "outer iteration: " << outerIteration+1 << " of " << _outerIterations << endl;
          cerr << "inner iteration: " << innerIteration+1 << " of " << _innerIterations << endl;

          //! compute \f$c_1\f$, \f$c_2\f$ for each image
          // for each image...
          for ( int i = 0; i < _numberOfImages; i++ ) {
            // reformulate the displacement
            aol::MultiVector<RealType> disp;
            for ( int j=0; j < ConfigurationType::Dim; j++ )
              disp.appendReference( dCur[j][i] );

            // compute $\int_{\Omega}H(\zeta(\phi_i(x)))dx$ =: "aux"
            aol::Scalar<RealType> aux, aux1, aux2;
            typedef HZPIntegratorOverDomain<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> > HZPIntegrator;
            HZPIntegrator simpleHZP( grid, imagesCur[i], disp, heavisideFunction, HZPIntegrator::SIMPLE );
            simpleHZP.apply( zetaCur, aux );

            // compute $\int_{\Omega}u_i(x)H(\zeta(\phi_i(x)))dx$ and put it in "aux2"
            HZPIntegrator greyHZP( grid, imagesCur[i], disp, heavisideFunction, HZPIntegrator::GREY );
            greyHZP.apply( zetaCur, aux2 );

            // compute $\int_{\Omega}u_i(x)(1-H(\zeta(\phi_i(x))))dx$ and put it in "aux1"
            HZPIntegrator greyInverseHZP( grid, imagesCur[i], disp, heavisideFunction, HZPIntegrator::GREY_INVERSE );
            greyInverseHZP.apply( zetaCur, aux1 );

            // compute the $c_2$ vector
            _c2[i] = aux2/aux;

            // compute the $c_2$ vector
            HZPIntegrator inverseHZP( grid, imagesCur[i], disp, heavisideFunction, HZPIntegrator::INVERSE );
            inverseHZP.apply( zetaCur, aux );
            _c1[i] = aux1/aux;

            /*// compute the $c_1$ vector, using $\int_{\Omega}1-H(\zeta(\phi_i(x)))dx=1-\int_{\Omega}H(\zeta(\phi_i(x)))dx$
            // (not possible for too strong deformations of the domain)
            aux -= 1.;
            aux *= -1.;
            _c1[i] = aux1/aux;
          */
          }
          cerr << "greyscales updated" << endl;

          /*// fix the appropriate greylevels, if they are already known (here for white object on black background)
          _c1.setAll(1);
          _c2.setAll(0);
          */

          //! perform gradient descent for the levelset function \f$\zeta\f$
          // define total energy $E$=:"E" and its variation with respect to $\zeta$, "DE"
          EnergyOfZeta<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> >
            E( grid, imagesCur, dCur, heavisideFunction, _c1, _c2, _mu, _nu, _lambda1, _lambda2 );
          EnergyVariationWRTZeta<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> >
            DE( grid, imagesCur, dCur, heavisideFunction, _c1, _c2, _mu, _nu, _lambda1, _lambda2 );

          // define a corresponding instance of the gradient descent method
          GradientDescentZeta<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> > gradientDescentZeta( grid, E, DE, heavisideFunction, 6 );

          // perform the gradient descent
          aol::Vector<RealType> updatedZeta( zetaCur, aol::STRUCT_COPY );
          gradientDescentZeta.apply( zetaCur, updatedZeta );
          zetaCur = updatedZeta;

          cerr << "levelset function updated" << endl;

          //! perform gradient descent for the displacement \f$(d_i)_i\f$
          // for each image ...
          for ( int i=0; i < _numberOfImages; i++ ){
            // define the i-th displacement as multivector
            aol::MultiVector<RealType> disp;
            for ( int j=0; j < ConfigurationType::Dim; j++ )
              disp.appendReference( dCur[j][i] );

            // define total energy $E$=:"E" and its variation with respect to the displacement, "DE"
            EnergyOfD<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> >
              E( grid, imagesCur[i], zetaCur, heavisideFunction, _c1[i], _c2[i], _lambda1, _lambda2, volumeEnergyWeight, lengthEnergyWeight );
            EnergyVariationWRTD<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> >
              DE( grid, imagesCur[i], zetaCur, heavisideFunction, _c1[i], _c2[i], _lambda1, _lambda2, volumeEnergyWeight, lengthEnergyWeight );

            /*// check whether derivatives are correct
            aol::testFirstDerivativeMultiVectorAllDirections<ConfigurationType>(grid,d,E,DE,"derTest");
            */

            // define a corresponding instance of the gradient descent method
            GradientDescentD<ConfigurationType, aol::ArcTanHeavisideFunction<RealType> > gradientDescentD( grid, E, DE, heavisideFunction, 6 );
            //gradientDescentD.setFilterWidth( 10 );

            // perform the gradient descent
            aol::MultiVector<RealType> updatedD( disp, aol::STRUCT_COPY );
            gradientDescentD.apply( disp, updatedD );
            for ( int j=0; j < ConfigurationType::Dim; j++ )
              dCur[j][i] = updatedD[j];

            /*// set displacements to zero, if needed (eg. for debugging)
            for ( int j=0; j < ConfigurationType::Dim; j++ )
              dCur[j][i].setZero();
            */

            cerr << "deformation " << i+1 << " updated" << endl;
          }

          //! save data
          /*// print out total energy of this step
          cerr << "Total energy of this step: " << totalEnergy( grid, imagesCur, zetaCur, heavisideFunction, dCur, _c1, _c2, volumeEnergyWeight, lengthEnergyWeight ) << endl;
          */
          // save levelset
          if (ConfigurationType::Dim==2)  // in 2D save levelset with overlaid images as pgm
            timeStepSaver.saveTimestepPNGDrawRedAndBlueIsoline( scale, outerIteration, 0., 0.1, zetaCur, scaledOverlaidImage.current(), grid );
          else {                         // in 3D save the levelset so that it can be read e.g. by Grape
            char fileName[1024];
            sprintf( fileName, "%s%s_%d_%d.%s", _destDirectory, _destFileNameStem, scale, outerIteration, _fileExtension );
            (qc::ScalarArray<RealType, qc::QC_3D>( zeta.current(), grid )).save( fileName, qc::PGM_DOUBLE_BINARY );
          }
          // produce the inversely deformed images and save them as well as the displacements
          for ( int i=0; i<_numberOfImages; i++ ) {
            // define the i-th displacement as multivector
            aol::MultiVector<RealType> disp;
            for ( int j=0; j < ConfigurationType::Dim; j++ )
              disp.appendReference( dCur[j][i] );
            // define the i-th image as multivector
            aol::MultiVector<RealType> image;
            image.appendReference( imagesCur[i] );
            // compute inversely deformed image
            qc::TransformFunction<RealType,ConfigurationType::Dim> inverseDeformation( grid );
            inverseDeformation.setDeformation( disp );
            aol::MultiVector<RealType> invDefImage( image );
            inverseDeformation.apply( image, invDefImage );
            // save this inversely deformed image
            if ( ConfigurationType::Dim == 2 ) {  // in 2D save as pgm with red levelset
              timeStepSaver.setSaveName( "deformedImage" );
              timeStepSaver.saveTimestepPNGDrawRedAndBlueIsoline( scale*10+outerIteration, i, 0., 0.1, zetaCur, invDefImage[0], grid );
              timeStepSaver.setSaveName( _destFileNameStem );
            }
            else {                               // in 3D save the arrays so that they can be read e.g. by Grape
              char fileName[1024];
              sprintf( fileName, "%sdeformedImage_%d_%d.%s", _destDirectory, scale*10+outerIteration, i, _fileExtension );
              (qc::ScalarArray<RealType, qc::QC_3D>( invDefImage[0], grid )).save( fileName, qc::PGM_DOUBLE_BINARY );
            }
            // save the displacement in x-, y- (and z-) direction
            for ( int j=0; j<ConfigurationType::Dim; j++ ) {
              char fileName[1024];
              // due to compression, file extension should be bz2
              sprintf( fileName, "%sdisplacement_%d_%d_%d.%s", _destDirectory, scale, i, j, "bz2" );
              (typename ConfigurationType::ArrayType( dCur[j][i], grid )).save( fileName, qc::PGM_DOUBLE_BINARY );
            }
          }
        }
        // make Heaviside function sharper
        delta /= 2.;
      }
      // prolongate results
      zeta.levProlongate();
      for ( int i=0; i<ConfigurationType::Dim; i++ )
        d.levProlongate();
    }
  }
};

//////////////////////////////////////////////////////////////////////////////////////
//! ---------------------------- main program --------------------------------------//
//////////////////////////////////////////////////////////////////////////////////////

// define the settings/configuration of the program, i.e. computation accuracy, dimension, gridtype,
// finite element type, used quadrature rules
typedef double AccuracyType;
typedef qc::QuocConfiguratorTraitMultiLin<AccuracyType, qc::QC_2D, aol::GaussQuadrature<AccuracyType, qc::QC_2D, 3> > SettingsType;

/**
 * \brief This main function coordinates the (Chan-Vese-) segmentation and registration/averaging of a list of images.
 * The parameter data is read in and the segmentation and registration algorithm is executed.
 * If started from the command line, argc is the number of passed arguments,
 * and argv is a pointer on the passed arguments (which are of type char[]).
 * The name of the parameter file should be passed as the only argument.
 *
 * \author Wirth
 */
int main( int argc, char *argv[] ) {

  try {
    char parameterFileName[1024];

    // read in file names of parameter files
    if ( argc > 10 ) {              // too many input arguments specified; explain correct syntax
      cerr << "Too many input files." << endl;
      return EXIT_FAILURE;
    } else if ( argc == 1 ) {     // no input argument specified; use default value for filename of parameter file
      sprintf( parameterFileName, "./shape.par" );

      // read in parameters from specified parameter file
      cerr << "Reading parameters from '" << parameterFileName << "'." << endl;
      aol::ParameterParser parameterParser( parameterFileName );

      // execute the algorithm
      ChanVeseAverageShapeFinder<SettingsType> chanVeseAverageShapeFinder( parameterParser );
      chanVeseAverageShapeFinder.execute();
    } else {                       // read in filenames of parameter files
      for ( int i=1; i<argc; i++ ) {
        sprintf( parameterFileName, "%s",  argv[i] );

        // read in parameters from specified parameter file
        cerr << endl << "Reading parameters from '" << parameterFileName << "'." << endl;
        aol::ParameterParser parameterParser( parameterFileName );

        // execute the algorithm
        ChanVeseAverageShapeFinder<SettingsType> chanVeseAverageShapeFinder( parameterParser );
        chanVeseAverageShapeFinder.execute();
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
