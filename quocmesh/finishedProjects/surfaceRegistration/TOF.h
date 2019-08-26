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

#ifndef __TOF_H
#define __TOF_H

/**
 * \brief  Transformation from range data to world coordinates
 * g[r]: \omega \longrightarrow \mathbb{R}^3
 *
 * \author Bauer
 */
template <typename RealType>
class RangeToWorldCoordTransformation {
  /// The focal length
  const RealType _focalLength;
  /// The sensor width (in elements)
  const RealType _width;
  /// Teh sensor height (in elements)
  const RealType _height;
  /// The shift in u direction
  RealType _uShift;
  /// The shift in v direction
  RealType _vShift;
  /// The scale factor
  const RealType _scaleFactor;
  /// The shift offset
  const RealType _shiftOffsetX;
  const RealType _shiftOffsetY;
  const RealType _shiftOffsetZ;
public:
  RangeToWorldCoordTransformation (
    const RealType FocalLength,
    const RealType Width,
    const RealType Height,
    const RealType ScaleFactor = 1.,
    const RealType ShiftOffsetX = 0.,
    const RealType ShiftOffsetY = 0.,
    const RealType ShiftOffsetZ = 0. )
    : _focalLength ( FocalLength ), _width ( Width ), _height ( Height ),
      _scaleFactor ( ScaleFactor ), _shiftOffsetX ( ShiftOffsetX ),
      _shiftOffsetY ( ShiftOffsetY ), _shiftOffsetZ ( ShiftOffsetZ ) {

    _uShift = ( _width - 1. ) / 2.;
    _vShift = ( _height - 1. ) / 2.;
  }

  aol::Vec3<RealType> evaluateGamma ( RealType u, RealType v ) const {

    u -= _uShift;
    v -= _vShift;

    aol::Vec3<RealType> gamma;

    //gamma[0] = _scaleFactor*u/_focalLength;
    //gamma[1] = _scaleFactor*v/_focalLength;
    //gamma[2] = _scaleFactor*1./(sqrt(1.+(pow(u,2)+pow(v,2))/(pow(_focalLength,2))));

    RealType frac_scaleFactor_focalLength = _scaleFactor / _focalLength;

    gamma[0] = frac_scaleFactor_focalLength * u;
    gamma[1] = frac_scaleFactor_focalLength * v;
    gamma[2] = _scaleFactor * 1. / ( sqrt ( 1. + ( pow ( u, 2 ) + pow ( v, 2 ) ) / ( pow ( _focalLength, 2 ) ) ) );

    return gamma;
  }

  void evaluateGammaDerivative ( RealType u, RealType v, aol::Vec3<RealType> &du, aol::Vec3<RealType> & dv ) const {

    u -= _uShift;
    v -= _vShift;

    RealType tmp = _focalLength / ( pow ( pow ( static_cast<double> ( _focalLength ), 2. ) + pow ( static_cast<double> ( u ), 2. ) + pow ( static_cast<double> ( v ), 2. ), 1.5 ) );

    //du[0] = _scaleFactor*1./_focalLength;
    //du[1] = 0.;
    //du[2] = -_scaleFactor*tmp*u;

    //dv[0] = 0.;
    //dv[1] = _scaleFactor*1./_focalLength;
    //dv[2] = -_scaleFactor*tmp*v;

    RealType frac_scaleFactor_focalLength = _scaleFactor / _focalLength;
    RealType mult_scaleFactor_tmp = _scaleFactor * tmp;

    du[0] = frac_scaleFactor_focalLength;
    du[1] = 0.;
    du[2] = -mult_scaleFactor_tmp * u;

    dv[0] = 0.;
    dv[1] = frac_scaleFactor_focalLength;
    dv[2] = -mult_scaleFactor_tmp * v;

    return;
  }

  aol::Vec3<RealType> evaluateWorldCoord ( RealType u, RealType v, const RealType r ) const {

    u -= _uShift;
    v -= _vShift;

    aol::Vec3<RealType> WorldCoord;

    //WorldCoord[0] = _scaleFactor*r*u/_focalLength;
    //WorldCoord[1] = _scaleFactor*r*v/_focalLength;
    //WorldCoord[2] = _scaleFactor*r*1./(sqrt(1.+(pow(u,2)+pow(v,2))/(pow(_focalLength,2)))) + _shiftOffset;

    RealType mult_scaleFactor_r = _scaleFactor * r;
    RealType frac_mult_scaleFactor_r_focalLength = mult_scaleFactor_r / _focalLength;

    WorldCoord[0] = frac_mult_scaleFactor_r_focalLength * u + _shiftOffsetX;
    WorldCoord[1] = frac_mult_scaleFactor_r_focalLength * v + _shiftOffsetY;
    WorldCoord[2] = mult_scaleFactor_r / ( sqrt ( 1. + ( pow ( u, 2 ) + pow ( v, 2 ) ) / ( pow ( _focalLength, 2 ) ) ) ) + _shiftOffsetZ;

    // shift to [0,1]^3
    WorldCoord[0] += _scaleFactor * 0.5;
    WorldCoord[1] += _scaleFactor * 0.5;
    WorldCoord[2] += _scaleFactor * 0.5;

    return WorldCoord;
  }
};

#endif
