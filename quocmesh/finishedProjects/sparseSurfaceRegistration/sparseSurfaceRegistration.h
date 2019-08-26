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

#include <configurators.h>
#include <op.h>
#include <projector.h>

#define USE_IMPROVED_APPROXIMATION 1

/**
 * Provides an easy interface to operators of the form \f$ \frac{1}{n}\sum_{i=1}^nf_i(u,p_i)\cdot\varphi_j(p_i)\f$,
 * where \f$u=\f$Arg, \f$f_i\f$ is supplied by getIthSummand and \f$p\f$ has to be supplied in the constructor.
 * The positions \f$p\f$ are assumed to be in [0,1]^d.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, int NumComponents, typename Imp>
class FENonlinSumInterface : public aol::FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  aol::MultiVector<RealType> _positions;
public:
  explicit FENonlinSumInterface ( const typename ConfiguratorType::InitType &Grid, aol::MultiVector<RealType> &Positions )
    : aol::FEOpInterface<ConfiguratorType, aol::MultiVector<RealType> > ( Grid ),
      _grid ( Grid ),
      _positions ( Positions ) { }

  explicit FENonlinSumInterface ( const typename ConfiguratorType::InitType &Grid )
    : aol::FEOpInterface<ConfiguratorType, aol::MultiVector<RealType> > ( Grid ),
      _grid ( Grid ) { }

  virtual ~FENonlinSumInterface( ) { }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    const aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumComponents> discrVecFunc ( this->getConfigurator(), MArg );

    const int numPositions = _positions.getEqualComponentSize();
    typename ConfiguratorType::DomVecType coords;
    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType transformedCoord;

    // Since we are going to devide by n, we first have to multiply by n in order not to scale the value that is already in MDest now.
    MDest *= static_cast<RealType> ( numPositions );

    for ( int position = 0; position < numPositions; ++position ) {
      _positions.getTo ( position, coords );

      // Get local coordinates and check if the position is inside the computational domain.
      if ( this->getConfigurator().getLocalCoords( coords, transformedEl, transformedCoord ) ) {
        aol::Vec<NumComponents, RealType> a;
        asImp().getIthSummand ( position, discrVecFunc, transformedEl, transformedCoord, a );

        // For each basis function whose support contains transformedEl, add its contribution to the sum.
        const typename ConfiguratorType::BaseFuncSetType &tbfs = this->getConfigurator().getBaseFunctionSet ( transformedEl );
        const int numTransformedLocalDofs = this->getConfigurator().getNumLocalDofs ( transformedEl );
        for ( int dof = 0; dof < numTransformedLocalDofs; dof++ ) {
          const RealType temp = tbfs.evaluate ( dof, transformedCoord );
          for ( int d = 0; d < NumComponents; ++d )
            MDest[d][ this->getConfigurator().localToGlobal ( transformedEl, dof ) ] += a[d] * temp;
        }
      }
    }

    MDest /= static_cast<RealType> ( numPositions );
  }

  //! Represents \f$f_i(u,p_i)\f$ and has to be implemented in the derived class.
  void getIthSummand ( const int I,
                       const aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumComponents> &DiscrVecFunc,
                       const typename ConfiguratorType::ElementType &TransformedEl,
                       const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                       aol::Vec<NumComponents, RealType> &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().getIthSummand ( I, DiscrVecFunc, TransformedEl, TransformedLocalCoord );
  }

protected:
  // Barton-Nackman trick
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class TestFittingEnergy : public aol::FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::MultiVector<RealType> &_positions;
  const aol::Vector<RealType> &_values;
public:
  TestFittingEnergy ( const typename ConfiguratorType::InitType &Initializer,
                      const aol::MultiVector<RealType> &Positions,
                      const aol::Vector<RealType> &Values )
    : aol::FEOpInterface<ConfiguratorType, aol::MultiVector<RealType>, aol::Scalar<RealType> > ( Initializer ),
      _grid ( Initializer ),
      _positions ( Positions ),
      _values ( Values ) { }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    const aol::DiscreteFunctionDefault<ConfiguratorType> discrFunc ( this->getConfigurator(), MArg[0] );

    const int numPositions = _positions.getEqualComponentSize();
    typename ConfiguratorType::DomVecType coords;
    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType transformedCoord;

    for ( int position = 0; position < numPositions; ++position ) {
      _positions.getTo ( position, coords );
      if ( this->getConfigurator().getLocalCoords( coords, transformedEl, transformedCoord ) )
        Dest[0] += 0.5 * aol::Sqr ( discrFunc.evaluate( transformedEl, transformedCoord ) - _values[position] );
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class FENonlinSumInterfaceTest : public FENonlinSumInterface<ConfiguratorType, 1, FENonlinSumInterfaceTest<ConfiguratorType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const aol::Vector<RealType> &_values;
public:

  FENonlinSumInterfaceTest ( const typename ConfiguratorType::InitType &Initializer,
                             const aol::MultiVector<RealType> &Positions,
                             const aol::Vector<RealType> &Values )
    : FENonlinSumInterface<ConfiguratorType, 1, FENonlinSumInterfaceTest<ConfiguratorType> > ( Initializer, Positions ),
      _values ( Values ) { }

  void getIthSummand ( const int I,
                       const aol::DiscreteVectorFunctionDefault<ConfiguratorType, 1> &DiscrVecFunc,
                       const typename ConfiguratorType::ElementType &TransformedEl,
                       const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                       aol::Vec<1, RealType> &NL ) const {
    NL[0] = ( DiscrVecFunc[0].evaluate( TransformedEl, TransformedLocalCoord ) - _values[I] );
  }
};
 
/**
 * \author Bauer
 */
template <typename ConfiguratorType3D>
class MatchingEnergy : public aol::FEOpInterface<ConfiguratorType3D, aol::MultiVector<typename ConfiguratorType3D::RealType>, aol::Scalar<typename ConfiguratorType3D::RealType> > {
public:
	
	typedef typename ConfiguratorType3D::RealType RealType;
	typedef typename ConfiguratorType3D::ElementType ElementType3D;
	typedef typename ConfiguratorType3D::DomVecType DomVecType3D;

    const typename ConfiguratorType3D::InitType & _grid3D;

	/// The discrete function of the SDF of the reference surface
	const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFuncs;
	/// The measured template data (triangulation sensor)
	const aol::MultiVector<RealType> & _measuredTemplateData;

	MatchingEnergy (
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &SDF,
		const aol::MultiVector<RealType> &MeasuredTemplateData ) 
		: aol::FEOpInterface<ConfiguratorType3D, aol::MultiVector<RealType>, aol::Scalar<RealType> > ( Grid3D ),
		_grid3D ( Grid3D ), _sdfDiscFuncs ( Grid3D, SDF ), _measuredTemplateData( MeasuredTemplateData ) {}
		
	void applyAdd (	const aol::MultiVector<RealType> &Arg, 
					aol::Scalar<RealType> &Dest ) const {
	
		ElementType3D TransformedEl;
		DomVecType3D TransformedCoord;

		// \sum\limits_{i=1}^n
		for ( int i = 0; i < _measuredTemplateData.getEqualComponentSize(); i++ ) 
		{
			//std::cout << Arg[0][i] << " " << Arg[1][i] << " " << Arg[2][i] << std::endl;

			// (y_i+w_i)
			DomVecType3D TransformedPoint( _measuredTemplateData[0][i] + Arg[0][i], _measuredTemplateData[1][i] + Arg[1][i], _measuredTemplateData[2][i] + Arg[2][i] );

			// | d(y_i+w_i) |^2
			if ( this->getConfigurator().getLocalCoords( TransformedPoint, TransformedEl, TransformedCoord ) )
			{
				Dest +=  0.5 / static_cast<RealType>(_measuredTemplateData.getEqualComponentSize()) * aol::Sqr( this->_sdfDiscFuncs.evaluate( TransformedEl, TransformedCoord ) ); 
			}			
		}	
	}
};


/**
 * \author Bauer
 */
template <typename ConfiguratorType3D>
class DerivativeOfMatchingEnergy : public aol::FEOpInterface<ConfiguratorType3D, aol::MultiVector<typename ConfiguratorType3D::RealType>, aol::MultiVector<typename ConfiguratorType3D::RealType> > {
public:
	
	typedef typename ConfiguratorType3D::RealType RealType;  		
	typedef typename ConfiguratorType3D::ElementType ElementType3D;
	typedef typename ConfiguratorType3D::DomVecType DomVecType3D;
	
	const typename ConfiguratorType3D::InitType & _grid3D;

	/// The discrete function of the SDF of the reference surface
	const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFuncs;
	/// The measured template data (triangulation sensor)
	const aol::MultiVector<RealType> & _measuredTemplateData;

	DerivativeOfMatchingEnergy (
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &SDF,
		const aol::MultiVector<RealType> &MeasuredTemplateData ) 
		: aol::FEOpInterface<ConfiguratorType3D, aol::MultiVector<RealType>, aol::MultiVector<RealType> > ( Grid3D ),
		_grid3D ( Grid3D ), _sdfDiscFuncs ( Grid3D, SDF ), _measuredTemplateData( MeasuredTemplateData ) {}
		
	void applyAdd (	const aol::MultiVector<RealType> &Arg, 
					aol::MultiVector<RealType> &Dest ) const {
	
		ElementType3D TransformedEl;
		DomVecType3D TransformedCoord;

		const int numPositions = _measuredTemplateData.getEqualComponentSize();

		for ( int i = 0; i < numPositions; i++ )
		{	
			// (y_i+w_i)
			DomVecType3D TransformedPoint( _measuredTemplateData[0][i] + Arg[0][i], _measuredTemplateData[1][i] + Arg[1][i], _measuredTemplateData[2][i] + Arg[2][i] );

			if ( this->getConfigurator().getLocalCoords( TransformedPoint, TransformedEl, TransformedCoord ) )
			{
				// \nabla d(y_i+w_i)
				aol::Vec<3, RealType> GradSDF;	
				this->_sdfDiscFuncs.evaluateGradient( TransformedEl, TransformedCoord, GradSDF );
				GradSDF.normalize();
						
				// d(y_i+w_i) \nabla d(y_i+w_i)
				GradSDF *= this->_sdfDiscFuncs.evaluate( TransformedEl, TransformedCoord ) / numPositions; 
				
				aol::Vec<3, RealType> tmp;
				Dest.getTo( i, tmp );
				tmp += GradSDF;
				Dest.set( i, tmp );
			}
		}		
	}
};

/**
 * \author Bauer
 */
template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class CorrelationEnergy : public aol::FEOpInterface<ConfiguratorType3D, aol::MultiVector<typename ConfiguratorType3D::RealType>, aol::Scalar<typename ConfiguratorType3D::RealType> > {
public:
	
	typedef typename ConfiguratorType2D::RealType RealType;
    typedef typename ConfiguratorType2D::ElementType ElementType2D;
	typedef typename ConfiguratorType2D::DomVecType	DomVecType2D;
	typedef typename ConfiguratorType3D::ElementType ElementType3D;
	typedef typename ConfiguratorType3D::DomVecType DomVecType3D;
	
	const typename ConfiguratorType2D::InitType & _grid2D;
	const typename ConfiguratorType3D::InitType & _grid3D;

	/// The discrete function of the SDF of the reference surface
	const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFuncs;
	/// The measured template data (triangulation sensor)
	const aol::MultiVector<RealType> & _measuredTemplateData;

	/// Reference to deformation estimate
	const aol::MultiVector<RealType> &_deformationEstimate;

	CorrelationEnergy (
		const typename ConfiguratorType2D::InitType &Grid2D,
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &SDF,
		const aol::MultiVector<RealType> &MeasuredTemplateData,
		const aol::MultiVector<RealType> &DeformationEstimate ) 
		: aol::FEOpInterface<ConfiguratorType3D, aol::MultiVector<RealType>, aol::Scalar<typename ConfiguratorType3D::RealType> > ( Grid3D ),
		_grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdfDiscFuncs ( Grid3D, SDF ), _measuredTemplateData( MeasuredTemplateData ),
		_deformationEstimate( DeformationEstimate ) {}
		
	void applyAdd (	const aol::MultiVector<RealType> &Arg, 
					aol::Scalar<RealType> &Dest ) const {
					
		// separate DispU, DispW (remark: MultiVectors DispU, DispW have different number of elements)
		aol::MultiVector<RealType> DispU( 0, 0 );
		aol::MultiVector<RealType> DispW( 0, 0 );
		for ( int i = 0; i < ConfiguratorType3D::Dim; i++ )
		{
		  DispU.appendReference( Arg[i] );		  
		}
		for ( int i = ConfiguratorType3D::Dim; i < ConfiguratorType3D::Dim*2; i++ )
		{
		  DispW.appendReference( Arg[i] );		  
		}

		ConfiguratorType2D Configurator2D( _grid2D );
		const aol::DiscreteVectorFunctionDefault<ConfiguratorType2D, 3 > DispUDiscFuncs ( Configurator2D, DispU );

		Projector<ConfiguratorType3D> Proj( _grid3D, _sdfDiscFuncs );

		ElementType2D PointOnOmegaEl;
		DomVecType2D PointOnOmegaCoord;
		ElementType3D TransformedPointPrevEstimateEl, TransformedPointCurrEstimateEl;
		DomVecType3D TransformedPointPrevEstimateCoord, TransformedPointCurrEstimateCoord;

		for ( int i = 0; i < _measuredTemplateData.getEqualComponentSize(); i++ )
		{
			// (y_i+w_i)
			DomVecType3D TransformedPointCurrEstimate( _measuredTemplateData[0][i] + DispW[0][i], _measuredTemplateData[1][i] + DispW[1][i], _measuredTemplateData[2][i] + DispW[2][i] );
			DomVecType3D TransformedPointPrevEstimate( _measuredTemplateData[0][i] + _deformationEstimate[3][i], _measuredTemplateData[1][i] + _deformationEstimate[4][i], _measuredTemplateData[2][i] + _deformationEstimate[5][i] );
				       
			if ( (this->getConfigurator().getLocalCoords(TransformedPointPrevEstimate, TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord)) &&
				 (this->getConfigurator().getLocalCoords(TransformedPointCurrEstimate, TransformedPointCurrEstimateEl, TransformedPointCurrEstimateCoord)) )
			{
				// P(y_i+w_i)
				DomVecType3D PointOnG;
				if(USE_IMPROVED_APPROXIMATION)
				{
					PointOnG = Proj.projectOnG( TransformedPointCurrEstimate, TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord, 
						TransformedPointCurrEstimateEl, TransformedPointCurrEstimateCoord );
				}
				else
				{
					PointOnG = Proj.projectOnG( TransformedPointCurrEstimate, TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord );
				}
				// \tilde{P}(y_i+w_i)
				DomVecType2D PointOnOmega = Proj.projectOnOmega( PointOnG );

				ConfiguratorType2D Configurator2D( _grid2D );
				if ( Configurator2D.getLocalCoords( PointOnOmega, PointOnOmegaEl, PointOnOmegaCoord ) )
				{
					// u(\tilde{P}(y_i+w_i))
					DomVecType3D DispUVec;
					DispUDiscFuncs.evaluate( PointOnOmegaEl, PointOnOmegaCoord, DispUVec );
					
					// | P(y_i+w_i) + u(\tilde{P}(y_i+w_i)) - y_i |^2
					DomVecType3D MeasuredTemplatePoint;
					_measuredTemplateData.getTo( i, MeasuredTemplatePoint );
					DomVecType3D Sum = PointOnG + DispUVec - MeasuredTemplatePoint;
					Dest += .5 / static_cast<RealType>(_measuredTemplateData.getEqualComponentSize()) * Sum*Sum;		
				}
			}
		}		
	}
};


/**
 * Wrapper class to debug CorrelationEnergy / DerivativeOfCorrelationEnergyWRTDispW.
 *
 * \author Berkels
 */
template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class CorrelationEnergyWRTDispW : public aol::Op<aol::MultiVector<typename ConfiguratorType3D::RealType>, aol::Scalar<typename ConfiguratorType3D::RealType> > {
public:
	
	typedef typename ConfiguratorType2D::RealType RealType;
  CorrelationEnergy<ConfiguratorType2D, ConfiguratorType3D> _E;
  const aol::MultiVector<RealType> &_dispU;

	CorrelationEnergyWRTDispW (
		const typename ConfiguratorType2D::InitType &Grid2D,
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &SDF,
		const aol::MultiVector<RealType> &MeasuredTemplateData,
    const aol::MultiVector<RealType> &DispU ) 
    : _E ( Grid2D, Grid3D, SDF, MeasuredTemplateData ),
      _dispU ( DispU ) {}
		
	void applyAdd (	const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector<RealType> mArg ( 0, 0 );
    mArg.appendReference ( _dispU );
    mArg.appendReference ( Arg );
    _E.applyAdd ( mArg, Dest );
  }
};

/**
 * \author Bauer
 */
template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class VariationOfCorrelationEnergyWRTDispU : public FENonlinSumInterface<ConfiguratorType2D, 3, VariationOfCorrelationEnergyWRTDispU<ConfiguratorType2D, ConfiguratorType3D> > {

	typedef typename ConfiguratorType2D::RealType RealType;	
	typedef typename ConfiguratorType2D::ElementType ElementType2D;
	typedef typename ConfiguratorType2D::DomVecType	DomVecType2D;
	typedef typename ConfiguratorType3D::ElementType ElementType3D;
	typedef typename ConfiguratorType3D::DomVecType DomVecType3D;

	const typename ConfiguratorType2D::InitType & _grid2D;	
	const typename ConfiguratorType3D::InitType & _grid3D;

	/// The discrete function of the SDF of the reference surface
	const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFuncs;
	/// The measured template data (triangulation sensor)
	const aol::MultiVector<RealType> & _measuredTemplateData;
	
	const aol::MultiVector<RealType> & _dispW;

	const ConfiguratorType3D _configurator3D;

	const Projector<ConfiguratorType3D> _proj;

	/// Reference to deformation estimate
	const aol::MultiVector<RealType> &_deformationEstimate;

public:
	VariationOfCorrelationEnergyWRTDispU (	
		const typename ConfiguratorType2D::InitType &Grid2D,
		const typename ConfiguratorType3D::InitType &Grid3D,
		const aol::Vector<RealType> &SDF,			
		aol::MultiVector<RealType> &MeasuredTemplateData,										
		aol::MultiVector<RealType> &DispW,
		const aol::MultiVector<RealType> &DeformationEstimate )
		// [BB] DispW != _positions. Looks like you don't want to pass this here. I think FENonlinSumInterface needs an extra constructor
		// that doesn't need _positions.
		// [SB] I was aware of DispW != _positions. Right, an extra constructor is a clean solution at this point. 
		: FENonlinSumInterface<ConfiguratorType2D, 3, VariationOfCorrelationEnergyWRTDispU<ConfiguratorType2D, ConfiguratorType3D> > ( Grid2D ),
		_grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdfDiscFuncs ( Grid3D, SDF ), _measuredTemplateData( MeasuredTemplateData ), _dispW( DispW ),
		_configurator3D( _grid3D ), _proj( _grid3D, _sdfDiscFuncs ), _deformationEstimate( DeformationEstimate ) { 
	

		ElementType2D PointOnOmegaEl;
		DomVecType2D  PointOnOmegaCoord;
		ElementType3D TransformedEl;
		DomVecType3D  TransformedCoord;

		// Since we didn't pass the positions to FENonlinSumInterface, we have to set the size manually.
		this->_positions.reallocate ( 2, _measuredTemplateData.getEqualComponentSize() );

		// set up _positions
		for ( int i = 0; i < _measuredTemplateData.getEqualComponentSize(); i++ )
		{
			// (y_i+w_i)
			DomVecType3D TransformedPoint( _measuredTemplateData[0][i] + DispW[0][i], _measuredTemplateData[1][i] + DispW[1][i], _measuredTemplateData[2][i] + DispW[2][i] );
			
			if ( _configurator3D.getLocalCoords( TransformedPoint, TransformedEl, TransformedCoord ) )
			{
				// P(y_i+w_i)
				DomVecType3D PointOnG = _proj.projectOnG( TransformedPoint, TransformedEl, TransformedCoord );
				// \tilde{P}(y_i+w_i)
				DomVecType2D PointOnOmega = _proj.projectOnOmega( PointOnG );

				this->_positions[0][i] = PointOnOmega[0];
				this->_positions[1][i] = PointOnOmega[1];		
			}
		}
	}
		
	void getIthSummand( const int I,
						const aol::DiscreteVectorFunctionDefault<ConfiguratorType2D, ConfiguratorType3D::Dim> &DiscrVecFunc,
						const typename ConfiguratorType2D::ElementType &TransformedEl,
						const typename ConfiguratorType2D::DomVecType &TransformedLocalCoord,
						aol::Vec<3, RealType> &NL ) const {

		// (y_i+w_i)
		DomVecType3D TransformedPointCurrEstimate( _measuredTemplateData[0][I] + _dispW[0][I], _measuredTemplateData[1][I] + _dispW[1][I], _measuredTemplateData[2][I] + _dispW[2][I] );
		DomVecType3D TransformedPointPrevEstimate( _measuredTemplateData[0][I] + _deformationEstimate[3][I], _measuredTemplateData[1][I] + _deformationEstimate[4][I], _measuredTemplateData[2][I] + _deformationEstimate[5][I] );
			
		ElementType3D TransformedPointPrevEstimateEl, TransformedPointCurrEstimateEl;
		DomVecType3D  TransformedPointPrevEstimateCoord, TransformedPointCurrEstimateCoord;

		if ( (_configurator3D.getLocalCoords(TransformedPointPrevEstimate, TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord)) &&
			 (_configurator3D.getLocalCoords(TransformedPointCurrEstimate, TransformedPointCurrEstimateEl, TransformedPointCurrEstimateCoord)) )
		{
			// P(y_i+w_i)
			DomVecType3D PointOnG;
			if(USE_IMPROVED_APPROXIMATION)
			{
				PointOnG = _proj.projectOnG( TransformedPointCurrEstimate, TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord, 
					TransformedPointCurrEstimateEl, TransformedPointCurrEstimateCoord );
			}
			else
			{
				PointOnG = _proj.projectOnG( TransformedPointCurrEstimate, TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord );
			}	             

			// u(\tilde{P}(y_i+w_i))
			DomVecType3D DispUVec;
			DiscrVecFunc.evaluate( TransformedEl, TransformedLocalCoord, DispUVec );
		
			// ( P(y_i+w_i) + u(\tilde{P}(y_i+w_i)) - y_i )
			DomVecType3D MeasuredTemplatePoint;
			_measuredTemplateData.getTo( I, MeasuredTemplatePoint );
			DomVecType3D Sum = PointOnG + DispUVec - MeasuredTemplatePoint;
			NL = Sum;
		}
		else
		{
			NL.setZero();
		}
	}
};


/**
 * \author Bauer
 */
template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class DerivativeOfCorrelationEnergyWRTDispW : public aol::FEOpInterface<ConfiguratorType3D, aol::MultiVector<typename ConfiguratorType3D::RealType>, aol::MultiVector<typename ConfiguratorType3D::RealType> > {
public:
	
	typedef typename ConfiguratorType3D::RealType RealType;
	typedef typename ConfiguratorType2D::ElementType ElementType2D;
	typedef typename ConfiguratorType2D::DomVecType	DomVecType2D;
	typedef typename ConfiguratorType3D::ElementType ElementType3D;
	typedef typename ConfiguratorType3D::DomVecType DomVecType3D;

	const typename ConfiguratorType2D::InitType & _grid2D;
    const typename ConfiguratorType3D::InitType & _grid3D;

	/// The discrete function of the SDF of the reference surface
	const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFuncs;
	/// The measured template data (triangulation sensor)
	const aol::MultiVector<RealType> & _measuredTemplateData;

	const aol::MultiVector<RealType> & _dispU;

	/// Reference to deformation estimate
	const aol::MultiVector<RealType> &_deformationEstimate;

	DerivativeOfCorrelationEnergyWRTDispW (
		const typename ConfiguratorType2D::InitType &Grid2D, 
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &SDF,
		const aol::MultiVector<RealType> &MeasuredTemplateData,
		const aol::MultiVector<RealType> &DispU,
		const aol::MultiVector<RealType> &DeformationEstimate ) 
		: aol::FEOpInterface<ConfiguratorType3D, aol::MultiVector<RealType>, aol::MultiVector<RealType> > ( Grid3D ),
		_grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdfDiscFuncs ( Grid3D, SDF ), _measuredTemplateData( MeasuredTemplateData ), 
		_dispU( DispU ), _deformationEstimate( DeformationEstimate ) {}
		
	void applyAdd (	const aol::MultiVector<RealType> &Arg, 
					aol::MultiVector<RealType> &Dest ) const {

		ConfiguratorType2D Configurator2D( _grid2D );
		const aol::DiscreteVectorFunctionDefault<ConfiguratorType2D, 3 > DispUDiscFuncs ( Configurator2D, _dispU );
		
		Projector<ConfiguratorType3D> Proj( _grid3D, _sdfDiscFuncs );

		ConfiguratorType3D Configurator3D( _grid3D );

		ElementType2D PointOnOmegaEl;
		DomVecType2D PointOnOmegaCoord;
		ElementType3D TransformedPointPrevEstimateEl, TransformedPointCurrEstimateEl;
		DomVecType3D TransformedPointPrevEstimateCoord, TransformedPointCurrEstimateCoord;

		const int numPositions = _measuredTemplateData.getEqualComponentSize();

		for ( int i = 0; i < numPositions; i++ )
		{
			// (y_i+w_i)
			DomVecType3D TransformedPointCurrEstimate( _measuredTemplateData[0][i] + Arg[0][i], _measuredTemplateData[1][i] + Arg[1][i], _measuredTemplateData[2][i] + Arg[2][i] );
			DomVecType3D TransformedPointPrevEstimate( _measuredTemplateData[0][i] + _deformationEstimate[3][i], _measuredTemplateData[1][i] + _deformationEstimate[4][i], _measuredTemplateData[2][i] + _deformationEstimate[5][i] );
			
			if ( (this->getConfigurator().getLocalCoords( TransformedPointPrevEstimate, TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord)) &&
				 (this->getConfigurator().getLocalCoords( TransformedPointCurrEstimate, TransformedPointCurrEstimateEl, TransformedPointCurrEstimateCoord)) )
			{
				// P(y_i+w_i)
				DomVecType3D PointOnG;
				if(USE_IMPROVED_APPROXIMATION)
				{
					PointOnG = Proj.projectOnG( TransformedPointCurrEstimate, TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord, 
						TransformedPointCurrEstimateEl, TransformedPointCurrEstimateCoord );
				}
				else
				{
					PointOnG = Proj.projectOnG( TransformedPointCurrEstimate, TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord );
				}
				// \tilde{P}(y_i+w_i)
				DomVecType2D PointOnOmega = Proj.projectOnOmega( PointOnG );
				
				if ( Configurator2D.getLocalCoords( PointOnOmega, PointOnOmegaEl, PointOnOmegaCoord ) )
				{
					// u(\tilde{P}(y_i+w_i))
					DomVecType3D DispUVec;
					DispUDiscFuncs.evaluate( PointOnOmegaEl, PointOnOmegaCoord, DispUVec );
						
					// ( P(y_i+w_i) + u(\tilde{P}(y_i+w_i)) - y_i )
					DomVecType3D MeasuredTemplatePoint;
					_measuredTemplateData.getTo( i, MeasuredTemplatePoint );
					aol::Vec<3, RealType> Sum = PointOnG + DispUVec - MeasuredTemplatePoint;

					// \nabla u
					aol::Mat<3,2,RealType> GradDispU;
					DispUDiscFuncs.evaluateGradient( PointOnOmegaEl, PointOnOmegaCoord, GradDispU );

					for ( int k = 0; k < 3; ++k ) 
					{
						// e_k
						DomVecType3D Ek( 0., 0., 0.);
						Ek[k] = 1;
						ElementType3D EkEl;
						DomVecType3D EkCoord;
						Configurator3D.getLocalCoords( Ek, EkEl, EkCoord );

						if(USE_IMPROVED_APPROXIMATION)
						{
							aol::Mat<3,3,RealType> IdentityMat;
							IdentityMat.setIdentity();

							// \nabla d (y_i + w_i^{m-1})
							aol::Vec3<RealType> GradSDFPrevEstimate;
							this->_sdfDiscFuncs.evaluateGradient( TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord, GradSDFPrevEstimate );

							// \nabla d (y_i + w_i)
							aol::Vec3<RealType> GradSDFCurrEstimate;
							this->_sdfDiscFuncs.evaluateGradient( TransformedPointCurrEstimateEl, TransformedPointCurrEstimateCoord, GradSDFCurrEstimate );

							// \nabla d (y_i + w_i^{m-1}) \nabla^T d (y_i + w_i)
							aol::Mat<3,3,RealType> GradSDFMatrix;
							GradSDFMatrix.makeTensorProduct(GradSDFPrevEstimate,GradSDFCurrEstimate);
							
							IdentityMat -= GradSDFMatrix;
							aol::Vec<3,RealType> TermA = IdentityMat*Ek;

							DomVecType2D QTermA = Proj.projectOnOmega(TermA);
							aol::Vec<3,RealType> TermB = (GradDispU*QTermA);

							Dest[k][i] += Sum * ( TermA + TermB ) / numPositions;
						}
						else
						{
							// Qe_k
							DomVecType3D PointOnGEk = Ek;
							DomVecType2D PointOnOmegaEk = Proj.projectOnOmega( PointOnGEk );

							// e_k + \nabla u Q e_k
							aol::Vec<3, RealType> temp = GradDispU * PointOnOmegaEk;
							temp += PointOnGEk;
							Dest[k][i] += Sum * ( temp ) / numPositions;
						}
					}
				}
			}
		}		
	}
};


/**
 * \author Bauer
 */
template <typename ConfiguratorType2D, typename ConfiguratorType3D, typename RegMatrixType>
class SparseSurfaceRegistrationEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType2D::RealType>, aol::Scalar<typename ConfiguratorType2D::RealType> > {
public:
	
	typedef typename ConfiguratorType3D::RealType RealType;
	typedef typename ConfiguratorType2D::ElementType ElementType2D;
	typedef typename ConfiguratorType2D::DomVecType	DomVecType2D;
	typedef typename ConfiguratorType3D::ElementType ElementType3D;
	typedef typename ConfiguratorType3D::DomVecType DomVecType3D;

	const typename ConfiguratorType2D::InitType & _grid2D;
    const typename ConfiguratorType3D::InitType & _grid3D;

	/// The discrete function of the SDF of the reference surface
	const aol::Vector<RealType> _sdf;
	const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFuncs;
	/// The measured template data (triangulation sensor)
	const aol::MultiVector<RealType> & _measuredTemplateData;
	const RegMatrixType &_regMat;

	/// The weighting factors
	const RealType _kappa;
	const RealType _lambda;
	const RealType _mu;

	/// Reference to deformation estimate
	const aol::MultiVector<RealType> &_deformationEstimate;

public:
	SparseSurfaceRegistrationEnergy( 	
		const typename ConfiguratorType2D::InitType &Grid2D,
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &SDF,
		const aol::MultiVector<RealType> &MeasuredTemplateData,
		const RegMatrixType &RegMat,
		const RealType Kappa,
		const RealType Lambda,
		const RealType Mu,
		const aol::MultiVector<RealType> &DeformationEstimate ) 
		: _grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdf ( SDF ), _sdfDiscFuncs ( Grid3D, SDF ), _measuredTemplateData( MeasuredTemplateData ), _regMat ( RegMat ),
		_kappa( Kappa ), _lambda( Lambda ), _mu( Mu ), _deformationEstimate( DeformationEstimate ) {}

	virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{

		// separate DispU, DispW (remark: MultiVectors DispU, DispW have different number of elements)
		aol::MultiVector<RealType> DispU( 0, 0 );
		aol::MultiVector<RealType> DispW( 0, 0 );
		for ( int i = 0; i < ConfiguratorType3D::Dim; i++ )
		{
		  DispU.appendReference( MArg[i] );		  
		}
		for ( int i = ConfiguratorType3D::Dim; i < ConfiguratorType3D::Dim*2; i++ )
		{
		  DispW.appendReference( MArg[i] );		  
		}

		// E_{match}
		MatchingEnergy<ConfiguratorType3D> EMatch( _grid3D, _sdf, _measuredTemplateData );
		EMatch.apply( DispW, Dest );
		Dest *= _kappa;

		// E_{corr}
		CorrelationEnergy<ConfiguratorType2D, ConfiguratorType3D> ECorr( _grid2D, _grid3D, _sdf, _measuredTemplateData, _deformationEstimate );
		Dest /= _lambda;
		ECorr.applyAdd( MArg, Dest );
		Dest *= _lambda;

		// E_{reg}
		aol::DiagonalBlockOp<RealType> DEReg ( _regMat );
		aol::QuadraticFormOp<aol::MultiVector<RealType> > EReg ( DEReg );
		Dest /= _mu;
		EReg.applyAdd( DispU, Dest );
		Dest *= _mu;
	}
  
	virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{
		aol::Scalar<RealType> tmp;
		apply( MArg, tmp );
		Dest += tmp;
	}
};


/**
 * \author Bauer
 */
template <typename ConfiguratorType2D, typename ConfiguratorType3D, typename RegMatrixType>
class VariationOfSparseSurfaceRegistrationEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType2D::RealType> > {
	
	typedef typename ConfiguratorType3D::RealType RealType;
	typedef typename ConfiguratorType2D::ElementType ElementType2D;
	typedef typename ConfiguratorType2D::DomVecType	DomVecType2D;
	typedef typename ConfiguratorType3D::ElementType ElementType3D;
	typedef typename ConfiguratorType3D::DomVecType DomVecType3D;

	const typename ConfiguratorType2D::InitType & _grid2D;
    const typename ConfiguratorType3D::InitType & _grid3D;

	/// The discrete function of the SDF of the reference surface
	const aol::Vector<RealType> _sdf;
	const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFuncs;
	/// The measured template data (triangulation sensor)
	aol::MultiVector<RealType> & _measuredTemplateData;
	const RegMatrixType &_regMat;

	/// The weighting factors
	const RealType _kappa;
	const RealType _lambda;
	const RealType _mu;

	/// Reference to deformation estimate
	const aol::MultiVector<RealType> &_deformationEstimate;
	
public:
	VariationOfSparseSurfaceRegistrationEnergy ( 		
		const typename ConfiguratorType2D::InitType &Grid2D,
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &SDF,
		aol::MultiVector<RealType> &MeasuredTemplateData,
    const RegMatrixType &RegMat,
		const RealType Kappa,
		const RealType Lambda,
		const RealType Mu,
		const aol::MultiVector<RealType> &DeformationEstimate ) 
		: _grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdf ( SDF ), _sdfDiscFuncs ( Grid3D, SDF ), _measuredTemplateData( MeasuredTemplateData ), _regMat ( RegMat ), 
		_kappa( Kappa ), _lambda( Lambda ), _mu( Mu ), _deformationEstimate( DeformationEstimate ) {}
	
	virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
   
		// separate DispU, DispW (remark: MultiVectors DispU, DispW have different number of elements)
		aol::MultiVector<RealType> DispU( 0, 0 );
		aol::MultiVector<RealType> VarOfDispU( 0, 0 );
		aol::MultiVector<RealType> DispW( 0, 0 );
		aol::MultiVector<RealType> VarOfDispW( 0, 0 );
		for ( int i = 0; i < ConfiguratorType3D::Dim; i++ )
		{
		  DispU.appendReference( MArg[i] );		  
		  VarOfDispU.appendReference( MDest[i] );
		}
		for ( int i = ConfiguratorType3D::Dim; i < ConfiguratorType3D::Dim*2; i++ )
		{
		  DispW.appendReference( MArg[i] );		 
		  VarOfDispW.appendReference( MDest[i] );
		}

		// E_{match}
		DerivativeOfMatchingEnergy<ConfiguratorType3D> DEMatch( _grid3D, _sdf, _measuredTemplateData );
		DEMatch.apply( DispW, VarOfDispW );
		VarOfDispW *= _kappa;

		// E_{corr}
		VariationOfCorrelationEnergyWRTDispU<ConfiguratorType2D, ConfiguratorType3D> VECorrWRTDispU( _grid2D, _grid3D, _sdf, _measuredTemplateData, DispW, _deformationEstimate );
		VECorrWRTDispU.apply( DispU, VarOfDispU );
		VarOfDispU *= _lambda;

		DerivativeOfCorrelationEnergyWRTDispW<ConfiguratorType2D, ConfiguratorType3D> DECorrWRTDispW( _grid2D, _grid3D, _sdf, _measuredTemplateData, DispU, _deformationEstimate );
		VarOfDispW /= _lambda;
		DECorrWRTDispW.applyAdd( DispW, VarOfDispW );		
		VarOfDispW *= _lambda;

		// E_{reg}
		aol::DiagonalBlockOp<RealType> DEReg ( _regMat );
		VarOfDispU /= _mu;
		DEReg.applyAdd( DispU, VarOfDispU );
		VarOfDispU *= _mu;
	}

	virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
		aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
		apply( MArg, tmp );
		MDest += tmp;
	}
};
