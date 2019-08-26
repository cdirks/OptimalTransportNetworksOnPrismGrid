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


/**
 * \author Bauer
 */
template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class SDFBase {
public:
	typedef typename ConfiguratorType2D::RealType RealType;
	/// The 2D grid (parameter domain)
	const typename ConfiguratorType2D::InitType & _grid2D;
	/// The 3D grid (SDF)
	const typename ConfiguratorType3D::InitType & _grid3D;
	/// The discrete function of the SDF
	const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFuncs;
	/// The discrete function of the orthogonal range data
	const aol::DiscreteFunctionDefault<ConfiguratorType2D> _rangeDiscFuncs;

	SDFBase (	
		const typename ConfiguratorType2D::InitType &Grid2D,
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &SDF,
		const aol::Vector<RealType> &Range)
		: _grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdfDiscFuncs ( Grid3D, SDF ), 
		_rangeDiscFuncs ( Grid2D, Range ) 
		{}

	/// Compute the deformed world coord for a quad point
	/// x(\xi)
	bool evaluateWorldCoordDeformed (	
		const aol::auto_container<3,aol::DiscreteFunctionDefault<ConfiguratorType2D> > 
			&DefDiscFuncs,
		const aol::DiscreteFunctionDefault<ConfiguratorType2D> &RangeDiscFuncs,
		const typename ConfiguratorType2D::ElementType &El, 
		int QuadPoint, 
		const typename ConfiguratorType2D::DomVecType &RefCoord,
		typename ConfiguratorType3D::ElementType &WorldCoordDeformedEl, 
		typename ConfiguratorType3D::VecType &WorldCoordDeformedLocalCoord ) const {
	
		// Coordinates of the point in the parameter domain [0,1]^2
		RealType u = El[0] + RefCoord[0];
		RealType v = El[1] + RefCoord[1];
		
		// Get range data at quad point by bilinear interpolation of the range data
		RealType r = RangeDiscFuncs.evaluateAtQuadPoint( El, QuadPoint );

		// Compute 3D world coord for this range value
		typename ConfiguratorType3D::VecType WorldCoord;
		WorldCoord[0] = u*_grid2D.H(); 
		WorldCoord[1] = v*_grid2D.H(); 
		WorldCoord[2] = r;
				
		typename ConfiguratorType3D::ElementType WorldCoordEl;
		typename ConfiguratorType3D::VecType WorldCoordLocalCoord;
		typename ConfiguratorType3D::VecType Deformation;
		for(int i = 0; i < 3; ++i)
		{
			// Get the 3D deformation vector
			Deformation[i] = DefDiscFuncs[i].evaluateAtQuadPoint( El, QuadPoint );	
			
			// Get the element IDs and local coords of the world coord
			WorldCoordEl[i] = static_cast<short> (WorldCoord[i]/_grid3D.H());
			WorldCoordLocalCoord[i] = WorldCoord[i]/_grid3D.H() - WorldCoordEl[i];		
		}

		// Compute deformed world coord
		if ( !qc::transformCoord<ConfiguratorType3D> ( _grid3D, WorldCoordEl, 
			WorldCoordLocalCoord, Deformation, WorldCoordDeformedEl, 
			WorldCoordDeformedLocalCoord )) {
			return false;
		}	
		return true;
	}
};


/**
 * \author Bauer
 */
template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class MatchingEnergy : public aol::FENonlinIntegrationVectorInterface<ConfiguratorType2D, 
	MatchingEnergy<ConfiguratorType2D, ConfiguratorType3D>, 3 >, 
	public SDFBase<ConfiguratorType2D, ConfiguratorType3D> {
public:
	typedef typename ConfiguratorType2D::RealType RealType;

	MatchingEnergy (	
		const typename ConfiguratorType2D::InitType &Grid2D,
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &Sdf,
		const aol::Vector<RealType> &Range )
		: aol::FENonlinIntegrationVectorInterface< ConfiguratorType2D, 
		MatchingEnergy<ConfiguratorType2D, ConfiguratorType3D>, 3 > ( Grid2D ),
		SDFBase<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D, Grid3D, Sdf, Range ) {}

	// d(x(\xi)+\phi(\xi))
	RealType evaluateIntegrand( 
		const aol::auto_container<3,aol::DiscreteFunctionDefault<ConfiguratorType2D> > 
			&DefDiscFuncs,
		const typename ConfiguratorType2D::ElementType &El,
		int QuadPoint, 
		const typename ConfiguratorType2D::DomVecType &RefCoord ) const {
	
		typename ConfiguratorType3D::ElementType WorldCoordDeformedEl;
		typename ConfiguratorType3D::VecType WorldCoordDeformedLocalCoord;

		// compute deformed world coord
		if(this->evaluateWorldCoordDeformed( DefDiscFuncs, this->_rangeDiscFuncs, El, QuadPoint, 
			RefCoord, WorldCoordDeformedEl, WorldCoordDeformedLocalCoord ))
		{
			// (d(x(\xi)+\phi(\xi)))^2
			RealType dSqr =  aol::Sqr( this->_sdfDiscFuncs.evaluate( WorldCoordDeformedEl, 
				WorldCoordDeformedLocalCoord ) ); 

			return 0.5 * ( dSqr );
		}
		else
		{
			return 0;
		}		
	}
};


/**
 * \author Bauer
 */
template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class VariationOfMatchingEnergy : public aol::FENonlinVectorOpInterface< ConfiguratorType2D, 3, 3, VariationOfMatchingEnergy<ConfiguratorType2D, ConfiguratorType3D> > , 
	public SDFBase<ConfiguratorType2D, ConfiguratorType3D> {
public:
	typedef typename ConfiguratorType2D::RealType RealType;

	VariationOfMatchingEnergy (  const typename ConfiguratorType2D::InitType &Grid2D,
				const typename ConfiguratorType3D::InitType &Grid3D, 
				const aol::Vector<RealType> &Sdf,			  
				const aol::Vector<RealType> &Range )
	: aol::FENonlinVectorOpInterface< ConfiguratorType2D, 3, 3, VariationOfMatchingEnergy<ConfiguratorType2D, ConfiguratorType3D> > ( Grid2D ),
		SDFBase<ConfiguratorType2D, ConfiguratorType3D> ( Grid2D, 
		Grid3D, Sdf, Range ) {}

	void getNonlinearity( 
		aol::auto_container<3,aol::DiscreteFunctionDefault<ConfiguratorType2D> > 
			&DefDiscFuncs, 
        const typename ConfiguratorType2D::ElementType &El, 
        int QuadPoint, 
		const typename ConfiguratorType2D::VecType &RefCoord,
		aol::Vec<3 , typename ConfiguratorType2D::RealType> &NL ) const {
				
		// deformed WorldCoord
		typename ConfiguratorType3D::ElementType WorldCoordDeformedEl;
		typename ConfiguratorType3D::VecType WorldCoordDeformedLocalCoord;

		if(this->evaluateWorldCoordDeformed( DefDiscFuncs, this->_rangeDiscFuncs, El, 
			QuadPoint, RefCoord, WorldCoordDeformedEl, WorldCoordDeformedLocalCoord ))
		{
			// d(x(\xi)+\phi(\xi))
			RealType d = this->_sdfDiscFuncs.evaluate( WorldCoordDeformedEl, 
				WorldCoordDeformedLocalCoord ); 
			// NL = \nabla \big( d(x(\xi)+\phi(\xi)) \big)
			this->_sdfDiscFuncs.evaluateGradient( WorldCoordDeformedEl, 
				WorldCoordDeformedLocalCoord, NL );

			NL *= d;
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
template <typename ConfiguratorType2D, typename ConfiguratorType3D, typename RegMatrixType>
class SurfaceRegistrationEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType2D::RealType>, aol::Scalar<typename ConfiguratorType2D::RealType> > {
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

public:
	SurfaceRegistrationEnergy( 	
		const typename ConfiguratorType2D::InitType &Grid2D,
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &SDF,
		const aol::MultiVector<RealType> &MeasuredTemplateData,
    const RegMatrixType &RegMat,
		const RealType Kappa,
		const RealType Lambda ) 
		: _grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdf ( SDF ), _sdfDiscFuncs ( Grid3D, SDF ), _measuredTemplateData( MeasuredTemplateData ), _regMat ( RegMat ),
		_kappa( Kappa ), _lambda( Lambda ) {}

	virtual void apply( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const{

		// E_{match}
		MatchingEnergy<ConfiguratorType2D, ConfiguratorType3D> EMatch( _grid2D, _grid3D, _sdf, _measuredTemplateData[2] );
		EMatch.apply( MArg, Dest );
		Dest *= _kappa;

		// E_{reg}, SECOND ORDER
		//aol::DiagonalBlockOp<RealType> DEReg ( _regMat );
		//aol::QuadraticFormOp<aol::MultiVector<RealType> > EReg ( DEReg );
		//Dest /= _lambda;
		//EReg.applyAdd( MArg, Dest );
		//Dest *= _lambda;

		// E_{reg}, FIRST ORDER
		// Assemble stiffness matrix L
		aol::StiffOp<ConfiguratorType2D> stiff( _grid2D, aol::ASSEMBLED );
		const qc::DisplacementLengthEnergy<ConfiguratorType2D> EReg(stiff);
		Dest /= _lambda;
		EReg.applyAdd( MArg, Dest );
		Dest *= _lambda;
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
class VariationOfSurfaceRegistrationEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType2D::RealType> > {
	
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
	
public:
	VariationOfSurfaceRegistrationEnergy ( 		
		const typename ConfiguratorType2D::InitType &Grid2D,
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::Vector<RealType> &SDF,
		aol::MultiVector<RealType> &MeasuredTemplateData,
    const RegMatrixType &RegMat,
		const RealType Kappa,
		const RealType Lambda ) 
		: _grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdf ( SDF ), _sdfDiscFuncs ( Grid3D, SDF ), _measuredTemplateData( MeasuredTemplateData ), _regMat ( RegMat ), 
		_kappa( Kappa ), _lambda( Lambda ) {}
	
	virtual void apply( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
   
		// E_{match}
		VariationOfMatchingEnergy<ConfiguratorType2D, ConfiguratorType3D> DEMatch( _grid2D, _grid3D, _sdf, _measuredTemplateData[2] );		
		DEMatch.apply( MArg, MDest );
		MDest *= _kappa;

		// E_{reg}, SECOND ORDER
		//aol::DiagonalBlockOp<RealType> DEReg ( _regMat );
		//MDest /= _lambda;
		//DEReg.applyAdd( MArg, MDest );
		//MDest *= _lambda;

		// E_{reg}, FIRST ORDER
		// Assemble stiffness matrix L
		aol::StiffOp<ConfiguratorType2D> stiff( _grid2D, aol::ASSEMBLED );
		const aol::DiagonalBlockOp<RealType> DEReg(stiff);
		MDest /= _lambda;
		DEReg.applyAdd( MArg, MDest );
		MDest *= _lambda;
	}

	virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const{
		aol::MultiVector<RealType> tmp( MArg, aol::STRUCT_COPY );
		apply( MArg, tmp );
		MDest += tmp;
	}
};
