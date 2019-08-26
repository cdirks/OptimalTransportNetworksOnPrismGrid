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

#ifndef __GRAPEINTERFACE1D_H
#define __GRAPEINTERFACE1D_H

/* Interface between the QuocMesh-library and adaptive projection
   in Grape. Just call the function ... with a 1D-QuocMesh-Grid and get a
   Triang-strucure.
*/

#include "grape2.h"

#ifdef USE_EXTERNAL_GRAPE

#include <scalarArray.h>
#include <multiVector.h>
#include <quoc.h>


/*!
 * functions convert given scalar or vector valued data on a Quoc-grid into an
 * Triang1d which can be handled by GRAPE. Args are the
 * grid and the data which is epxected in a \b qcscalarArray1d
 * "dataName" is the name under which the data is presented in GRAPE!
 * The function returns a TRIANG1D*.
 * \remarks As an example see tools/grape/scalar1d.cpp
 * \author Lenz
 */

TRIANG1D* quocmesh_convert_to_triang1d(qc::ScalarArray<double, qc::QC_1D>* triangData, const char* dataName);

/*! Do some preliminary stuff before starting grape, usually called from within initStartGrape */
void addMethodsAndProjects1d();

/* Init and start GRAPE, params are the triang */
void initStartGrape(TRIANG1D* triang);

#endif

#endif
