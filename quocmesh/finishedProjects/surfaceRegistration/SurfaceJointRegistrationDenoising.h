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

#ifndef __SURFACEJOINTREGISTRATIONDENOISING_H
#define __SURFACEJOINTREGISTRATIONDENOISING_H

#include "SurfaceDenoising.h"
#include <RudinOsherFatemi.h>

// test

template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class CompleteRegistrationDenoisingEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType2D::RealType>, aol::Scalar<typename ConfiguratorType2D::RealType> > {
  typedef typename ConfiguratorType2D::RealType RealType;

  /// The 2D grid (2D range data, parameter domain)
  const typename ConfiguratorType2D::InitType & _grid2D;
  /// The 3D grid (3D world coord data, SDF of reference surface)
  const typename ConfiguratorType3D::InitType & _grid3D;

  /// The SDF of the reference surface
  const aol::Vector<RealType> & _sdf;

  /// The discrete function of the SDF of the reference surface
  const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFunc;

  /// The measured (noisy) ToF range data
  const aol::Vector<RealType> & _templateRangeData;

  /// The converter from range data to 3D world coordinates
  const RangeToWorldCoordTransformation<RealType> & _rangeToWorldCoord;

  int _smoothingType;
  int _differenceType;

  const RealType _alpha;

  // The weighting factors
  const RealType _kappa;
  const RealType _lambda;
  const RealType _mu;

  const RealType _delta1;
  const RealType _delta2;

  /// The stiffness matrix
  aol::StiffOp<ConfiguratorType2D> _stiff;

public:
  CompleteRegistrationDenoisingEnergy (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const typename ConfiguratorType3D::InitType &Grid3D,
    //aol::StiffOp<ConfiguratorType2D> &Stiff,
    const aol::Vector<RealType> &SDF,
    const aol::Vector<RealType> &TemplateRangeData,
    const RangeToWorldCoordTransformation<RealType> & RangeToWorldCoord,
    int SmoothingType,
    int DifferenceType,
    const RealType Alpha,
    const RealType Kappa,
    const RealType Lambda,
    const RealType Mu,
    RealType Delta1,
    RealType Delta2 )
    : _grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdf ( SDF ), _sdfDiscFunc ( Grid3D, SDF ),
      _templateRangeData ( TemplateRangeData ), _rangeToWorldCoord ( RangeToWorldCoord ),
      _smoothingType ( SmoothingType ), _differenceType ( DifferenceType ), _alpha ( Alpha ), _kappa ( Kappa ), _lambda ( Lambda ), _mu ( Mu ), _delta1 ( Delta1 ), _delta2 ( Delta2 ), _stiff ( Grid2D, aol::ASSEMBLED ) {}

  virtual void apply ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {

    aol::MultiVector<RealType> deformation ( 0, 0 );
    for ( int i = 0; i < 3; i++ ) {
      deformation.appendReference ( MArg[i] );
    }
    aol::Vector<RealType> r ( MArg[3], aol::FLAT_COPY );

    // E_MA
    // similarity measure (to reference surface)
    RealType weighting = 0.;
    if ( _smoothingType == 0 ) {
      weighting = _alpha;
    } else if ( _smoothingType > 0 ) {
      weighting = 0.;
    }
    //MatchingEnergyOfPhi<ConfiguratorType2D, ConfiguratorType3D> ESim( _grid2D, _grid3D, _sdf, _rangeToWorldCoord, r, weighting );
    //ESim.apply( deformation, Dest );
    MatchingEnergyOfRange<ConfiguratorType2D, ConfiguratorType3D> ESim ( _grid2D, _grid3D, _sdf, _rangeToWorldCoord, deformation, weighting );
    ESim.apply ( r, Dest );
    Dest *= _kappa;

    // E_DEN
    const aol::IsoEnergyOp<ConfiguratorType2D> EIso ( _grid2D, _delta2 );
    if ( _smoothingType > 0 ) {
      Dest /= _alpha;
      EIso.applyAdd ( r, Dest );
      Dest *= _alpha;
    }

    // E_REG
    // Assemble stiffness matrix L
    //aol::StiffOp<ConfiguratorType2D> stiff( _grid2D, aol::ASSEMBLED );
    const qc::DisplacementLengthEnergy<ConfiguratorType2D> EReg ( _stiff );
    Dest /= _lambda;
    EReg.applyAdd ( deformation, Dest );
    Dest *= _lambda;

    // E_SIM
    L1DifferenceEnergy<ConfiguratorType2D> EL1Difference ( _grid2D, _templateRangeData, _rangeToWorldCoord, _delta1 );
    L2DifferenceEnergy<ConfiguratorType2D> EL2Difference ( _grid2D, _templateRangeData, _rangeToWorldCoord, _delta1 );
    Dest /= _mu;
    switch ( _differenceType ) {
      case 1: // L1
        EL1Difference.applyAdd ( r, Dest );
        break;

      case 2: // L2
        EL2Difference.applyAdd ( r, Dest );
        break;
      default:
        throw aol::UnimplementedCodeException ( "Unsupported _differenceType", __FILE__, __LINE__ );
    }
    Dest *= _mu;
  }

  virtual void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    aol::Scalar<RealType> tmp;
    apply ( MArg, tmp );
    Dest += tmp;
  }
};


template <typename ConfiguratorType2D, typename ConfiguratorType3D>
class VariationOfCompleteRegistrationDenoisingEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType2D::RealType> > {
  typedef typename ConfiguratorType2D::RealType RealType;

  /// The 2D grid (2D range data, parameter domain)
  const typename ConfiguratorType2D::InitType & _grid2D;
  /// The 3D grid (3D world coord data, SDF of reference surface)
  const typename ConfiguratorType3D::InitType & _grid3D;

  /// The SDF of the reference surface
  const aol::Vector<RealType> & _sdf;

  /// The discrete function of the SDF of the reference surface
  const aol::DiscreteFunctionDefault<ConfiguratorType3D> _sdfDiscFunc;

  /// The measured (noisy) ToF range data
  const aol::Vector<RealType> & _templateRangeData;

  /// The converter from range data to 3D world coordinates
  const RangeToWorldCoordTransformation<RealType> & _rangeToWorldCoord;

  int _smoothingType;
  int _differenceType;

  const RealType _alpha;

  // Scale factors
  const RealType _kappa;
  const RealType _lambda;
  const RealType _mu;

  const RealType _delta1;
  const RealType _delta2;

  /// The stiffness matrix
  aol::StiffOp<ConfiguratorType2D> _stiff;

public:
  VariationOfCompleteRegistrationDenoisingEnergy (
    const typename ConfiguratorType2D::InitType &Grid2D,
    const typename ConfiguratorType3D::InitType &Grid3D,
    //aol::StiffOp<ConfiguratorType2D> &Stiff,
    const aol::Vector<RealType> &SDF,
    const aol::Vector<RealType> &TemplateRangeData,
    const RangeToWorldCoordTransformation<RealType> & RangeToWorldCoord,
    int SmoothingType,
    int DifferenceType,
    const RealType Alpha,
    const RealType Kappa,
    const RealType Lambda,
    const RealType Mu,
    const RealType Delta1,
    const RealType Delta2 )
    : _grid2D ( Grid2D ), _grid3D ( Grid3D ), _sdf ( SDF ), _sdfDiscFunc ( Grid3D, SDF ),
      _templateRangeData ( TemplateRangeData ), _rangeToWorldCoord ( RangeToWorldCoord ),
      _smoothingType ( SmoothingType ), _differenceType ( DifferenceType ), _alpha ( Alpha ), _kappa ( Kappa ), _lambda ( Lambda ), _mu ( Mu ),
      _delta1 ( Delta1 ), _delta2 ( Delta2 ), _stiff ( Grid2D, aol::ASSEMBLED ) {}

  virtual void apply ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {

    aol::MultiVector<RealType> deformation ( 0, 0 );
    aol::MultiVector<RealType> varOfDeformation ( 0, 0 );
    for ( int i = 0; i < 3; i++ ) {
      deformation.appendReference ( MArg[i] );
      varOfDeformation.appendReference ( MDest[i] );
    }
    aol::Vector<RealType> r ( MArg[3], aol::FLAT_COPY );
    aol::Vector<RealType> varOf_r ( MDest[3], aol::FLAT_COPY );

    // E_MA
    VariationOfMatchingEnergyWRTPhi<ConfiguratorType2D, ConfiguratorType3D> ESimVarOfPhi ( _grid2D, _grid3D, _sdf, _rangeToWorldCoord, r );
    ESimVarOfPhi.apply ( deformation, varOfDeformation );
    varOfDeformation *= _kappa;

    FirstPartOfVariationOfMatchingEnergyWRTRange<ConfiguratorType2D, ConfiguratorType3D>  ESimVarOfR ( _grid2D, _grid3D, _sdf, _rangeToWorldCoord, deformation );
    RealType weighting = 0.;
    if ( _smoothingType == 0 ) {
      weighting = _alpha;
    } else if ( _smoothingType > 0 ) {
      weighting = 0.;
    }
    //SecondPartOfVariationOfMatchingEnergyWRTRange<ConfiguratorType2D, ConfiguratorType3D> ESimVarOfR1( _grid2D, _grid3D, _sdf, _rangeToWorldCoord, deformation, weighting );
    //ThirdPartOfVariationOfMatchingEnergyWRTRange<ConfiguratorType2D, ConfiguratorType3D> ESimVarOfR2( _grid2D, _grid3D, _sdf, _rangeToWorldCoord, deformation, weighting );
    ESimVarOfR.apply ( r, varOf_r );
    //ESimVarOfR1.applyAdd( r, varOf_r );
    //ESimVarOfR2.applyAdd( r, varOf_r );
    varOf_r *= _kappa;

    // EReg
    // Assemble stiffness matrix L
    //aol::StiffOp<ConfiguratorType2D> stiff( _grid2D, aol::ASSEMBLED );
    const aol::DiagonalBlockOp<RealType> ERegVarOfPhi ( _stiff );
    varOfDeformation /= _lambda;
    ERegVarOfPhi.applyAdd ( deformation, varOfDeformation );
    varOfDeformation *= _lambda;

    // EDen
    VariationOfL1DifferenceEnergy <ConfiguratorType2D> DEL1Difference ( _grid2D, _templateRangeData, _rangeToWorldCoord, _delta1 );
    VariationOfL2DifferenceEnergy <ConfiguratorType2D> DEL2Difference ( _grid2D, _templateRangeData, _rangeToWorldCoord, _delta1 );
    varOf_r /= _mu;
    switch ( _differenceType ) {
      case 1: // L1
        DEL1Difference.applyAdd ( r, varOf_r );
        break;

      case 2: // L2
        DEL2Difference.applyAdd ( r, varOf_r );
        break;
      default:
        throw aol::UnimplementedCodeException ( "Unsupported _differenceType", __FILE__, __LINE__ );
    }
    varOf_r *= _mu;

    // TV
    const aol::VariationOfIsoEnergyOp<ConfiguratorType2D> EIso ( _grid2D, _delta2 );
    if ( _smoothingType > 0 ) {
      varOf_r /= _alpha;
      EIso.applyAdd ( r, varOf_r );
      varOf_r *= _alpha;
    }
  }

  virtual void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    aol::MultiVector<RealType> tmp ( MArg, aol::STRUCT_COPY );
    apply ( MArg, tmp );
    MDest += tmp;
  }
};



#endif //__SURFACEJOINTREGISTRATIONDENOISING_H
