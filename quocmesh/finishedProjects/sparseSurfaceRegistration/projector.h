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

#ifndef __PROJECTOR_H
#define __PROJECTOR_H

#include "math.h"

/**
 * \author Bauer
 */
template <typename ConfiguratorType3D>
class Projector {

	typedef typename ConfiguratorType3D::RealType RealType;
	const typename ConfiguratorType3D::InitType & _grid3D;

	/// The discrete function of the SDF of the reference surface
	const aol::DiscreteFunctionDefault<ConfiguratorType3D> &_sdfDiscFuncs;
		
public:
	Projector( 
		const typename ConfiguratorType3D::InitType &Grid3D, 
		const aol::DiscreteFunctionDefault<ConfiguratorType3D> &SdfDiscFuncs )
	: _grid3D( Grid3D), _sdfDiscFuncs ( SdfDiscFuncs ) { }

	// projection approximation as implemented for our MICCAI 2012 paper
	aol::Vec3<RealType> projectOnG ( aol::Vec3<RealType> TransformedPointCurrEstimate, typename ConfiguratorType3D::ElementType TransformedPointPrevEstimateEl, typename ConfiguratorType3D::VecType TransformedPointPrevEstimateCoord ) const {

		aol::Vec3<RealType> GradSDF;
		this->_sdfDiscFuncs.evaluateGradient( TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord, GradSDF );
		GradSDF.normalize();		
		return (TransformedPointCurrEstimate - this->_sdfDiscFuncs.evaluate( TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord ) * GradSDF);		
	}

	// improved projection approximation, where the scaling term of the nonlinear part is evaluated at the current estimate
	aol::Vec3<RealType> projectOnG ( aol::Vec3<RealType> TransformedPointCurrEstimate, typename ConfiguratorType3D::ElementType TransformedPointPrevEstimateEl, typename ConfiguratorType3D::VecType TransformedPointPrevEstimateCoord,
		typename ConfiguratorType3D::ElementType TransformedPointCurrEstimateEl, typename ConfiguratorType3D::VecType TransformedPointCurrEstimateCoord) const {
			
		aol::Vec3<RealType> GradSDF;
		this->_sdfDiscFuncs.evaluateGradient( TransformedPointPrevEstimateEl, TransformedPointPrevEstimateCoord, GradSDF );
		GradSDF.normalize();		
		return (TransformedPointCurrEstimate - this->_sdfDiscFuncs.evaluate( TransformedPointCurrEstimateEl, TransformedPointCurrEstimateCoord ) * GradSDF);		
	}
	
	aol::Vec2<RealType> projectOnOmega ( aol::Vec<3, RealType> PointOnG ) const {
				
		aol::Vec2<RealType> PointOnOmega( PointOnG[0], PointOnG[1] );
		
		return PointOnOmega;
	}
};

#endif // __PROJECTOR_H
