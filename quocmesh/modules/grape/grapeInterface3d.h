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

#ifndef __GRAPEINTERFACE3D_H
#define __GRAPEINTERFACE3D_H


/* Interface between the QuocMesh-library and adaptive projection
  in Grape. Just call the function ... with a 3D-QuocMesh-Grid and get a
  GenMesh-strucure.
  */

#include "grape2.h"

#ifdef USE_EXTERNAL_GRAPE

#include <scalarArray.h>
#include <multiVector.h>
#include <gridBase.h>
#include <estimator3d.h>
#include <quoc.h>

/*!
* function converts given scalar or vector-valued data on a Quoc-grid into an
* genmesh (adaptive) which can be handled by GRAPE. Args are the
* grid and the data which is epxected in a \b qcscalarArray3d.
* The function returns a GENMESH3D*.
* \remarks As an example see example3d.cpp (in the src/projects/grape_interface3d - directory)
* \author Nemitz
*/


/*! function expects only one qcScalarArray! */
template <typename REAL>
GENMESH3D* quocmesh_convert_to_gmesh3d(qc::ScalarArray<REAL, qc::QC_3D>* meshData, const char* dataName);

/*! function expects only one aol::MultiVector! KNOWN BUG: Adding vector-valued data to scalar-data will cause an
incorrect display of the scalar-data. Only vector-valued-data works fine!
Only scalar-data works fine too! */
template <typename REAL>
GENMESH3D* quocmesh_convert_to_gmesh3d(aol::MultiVector<REAL>* meshData,
                                       const char* dataName); // = "Vektordaten");

/*! The functions with the qc::Estimator3d-Structure */
/*! function expects the data and a saturated qc::Estimator3d-Structure (use makeSaturateErrorArray) */
template <typename REAL>
GENMESH3D* quocmesh_convert_to_gmesh3d(qc::ScalarArray<REAL, qc::QC_3D>* meshData,
                                       qc::Estimator3d<REAL>* errorArray, const char* dataName); // = "Daten mit Est." );

/*! add a scalar-data-array to the existing mesh! */
template <typename REAL>
    void addScalarData(GENMESH3D* gmesh, qc::ScalarArray<REAL, qc::QC_3D>* meshData, const char* dataName);

/*! add a scalar-data-array with a saturated (use makeSaturatedErrorArray)
    qcEstimator to the existing mesh! */
template <typename REAL>
    void addScalarDataWithEstimator(GENMESH3D* gmesh, qc::ScalarArray<REAL, qc::QC_3D>* meshData, const char* dataName,
                                qc::Estimator3d<REAL>* errorArray);

/*! add a vector-data-array to the existing mesh! */
template <typename REAL>
    void addVectorData(GENMESH3D* gmesh, aol::MultiVector<REAL>* meshData, const char* dataName);

/*! Add another timestep to scalar valued data */
void addTimestep (GENMESH3D* gmesh, qc::ScalarArray<double, qc::QC_3D>* meshData, const char* dataName, qc::Estimator3d<double>* errorArray=NULL);

/*! Do some preliminary stuff before starting grape, usually called from within initStartGrape */
void addMethodsAndProjects3d();

/* Init and start GRAPE, params are the mesh and a function, which is selected in GRAPE */
void initStartGrape(GENMESH3D* gmesh, const char* dataName);

/* the same, without selected function name */
void initStartGrape(GENMESH3D* gmesh);

/* start and create timescene */
void initStartGrapeTime(GENMESH3D* gmesh);

#endif

#endif
