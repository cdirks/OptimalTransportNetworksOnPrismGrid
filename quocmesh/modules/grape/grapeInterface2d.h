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

#ifndef __GRAPEINTERFACE2D_H
#define __GRAPEINTERFACE2D_H

/* Interface between the QuocMesh-library and adaptive projection
   in Grape. Just call the function ... with a 2D-QuocMesh-Grid and get a
   GenMesh-strucure.
*/

#include "grape2.h"

#ifdef USE_EXTERNAL_GRAPE

#include <scalarArray.h>
#include <multiVector.h>
#include <gridBase.h>
#include <quoc.h>
#include <estimator2d.h>


/*!
 * functions convert given scalar or vector valued data on a Quoc-grid into an
 * genmesh (adaptive) which can be handled by GRAPE. Args are the
 * grid and the data which is epxected in a \b qcscalarArray2d or a \b qcMultiVector.
 * "dataName" is the name under which the data is presented in GRAPE!
 * The function returns a GENMESH2D*.
 * \remarks As an example see example.cpp (in the src/projects/grape_interface2d - directory)
 * \author Nemitz
 */


GENMESH2D* quocmesh_convert_to_gmesh2d(qc::ScalarArray<double, qc::QC_2D>* meshData, const char* dataName);

GENMESH2D* quocmesh_convert_to_gmesh2d(aol::MultiVector<double>* meshData, const char* dataName);

GENMESH2D* quocmesh_convert_to_gmesh2d(qc::ScalarArray<double, qc::QC_2D>* meshData, qc::Estimator2d<double>* ErrorArray, const char* dataName);


/*! add a scalar-data-array to the existing mesh! */
void addScalarData(GENMESH2D* gmesh, qc::ScalarArray<double, qc::QC_2D>* meshData, const char* dataName);

/*! add a scalar-data-array with a saturated (use makeSaturatedErrorArray)
    qcEstimator to the existing mesh! */
void addScalarDataWithEstimator(GENMESH2D* gmesh, qc::ScalarArray<double, qc::QC_2D>* meshData, const char* dataName,
        qc::Estimator2d<double>* errorArray);

/*! add a vector-data-array to the existing mesh! */
void addVectorData(GENMESH2D* gmesh, aol::MultiVector<double>* meshData, const char* dataName);

/*! Add another timestep to scalar valued data */
void addTimestep (GENMESH2D* gmesh, qc::ScalarArray<double, qc::QC_2D>* meshData, const char* dataName, qc::Estimator2d<double>* errorArray=NULL);

/*! Add another timestep to vector valued data */
void addTimestep (GENMESH2D* gmesh, aol::MultiVector<double>* meshData, const char* dataName);

/*! Do some preliminary stuff before starting grape, usually called from within initStartGrape */
void addMethodsAndProjects2d();

/* Init and start GRAPE, params are the mesh and a function, which is selected in GRAPE */
void initStartGrape(GENMESH2D* gmesh, const char* dataName);

/* the same, without selected function name */
void initStartGrape(GENMESH2D* gmesh);

/* start and create timescene */
void initStartGrapeTime(GENMESH2D* gmesh);

void add_reload_button( GENMESH2D *mesh, qc::ScalarArray<double, qc::QC_2D> &data, const char *filename );

#endif

#endif
