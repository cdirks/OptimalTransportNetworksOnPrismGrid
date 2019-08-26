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

#ifndef __ESTIMATORTPTRAIT_H
#define __ESTIMATORTPTRAIT_H

#include "estimatorTP.h"

// forward declaration
template <typename T> class EstimatorTP2d;
template <typename T> class EstimatorTP3d;

/**
 * Trait to select the proper estimator op based on the Configurator template.
 *
 *  \author Paetz
 */
template <typename T, qc::Dimension Dim>
class EstimatorTPTrait {};

template <typename T>
class EstimatorTPTrait<T, qc::QC_2D> {
public:
  typedef EstimatorTP2d<T> EstimatorType;
};

template <typename T>
class EstimatorTPTrait<T, qc::QC_3D> {
public:
  typedef EstimatorTP3d<T> EstimatorType;
};

#endif
