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

//! @file
//!
//! @brief Example for using the GRAPE interface for a 2d scalar array, including error estimators.
//!
//! this file is an example for handling data of a 2D-QuocMesh with
//! the grape-interface. Scalar data has to be in a qcScalarArray<QC_2D>,
//! vector-valued data in a qcMultiVector (see example below!)
//!
//! @author Nemitz

#include <aol.h>

#include "grapeInterface2d.h"

#ifdef USE_EXTERNAL_GRAPE


int main(void)
{

  try
  {
    /* Depth of the grid */
    int d=11;
    int N = (1<<d)+1;


    /* **********************************************************
    *   Now some example data, that should be displayed
    *   with GRAPE, first some scalar, then vector-valued data
    * ********************************************************** */

    // fill a < ---- 2D-ARRAY --- > with some scalar data

    cerr<<"Getting the Array...";
    qc::ScalarArray<double, qc::QC_2D> ExampleDataSc(N,N);
    cerr<<"filling the array with life ...";
    for (int i=0; i<N; ++i)
      for (int j=0; j<N; ++j)
        ExampleDataSc.set(i,j,sin(3.*(i+j)/static_cast<double>(N))/2.);
        //ExampleDataSc.set(i,j,0.5*sqrt(qcAbs(sin(6.*(i+j)/((double)N))/2. * cos(6.*(i+j)/((double)N))/2.)));


    cerr<<"getting the estimator...";

    /* **********************************************************
    // the following entry is needed for the error-estimator */
    /* **************************************************************** */

    qc::Estimator2d<double> ErrorArray(N);
    cerr<<"saturating...";
    ErrorArray.makeSaturatedErrorArray(&ExampleDataSc);
    cerr<<"ready!"<<endl;

    // fill a < --- MULTIVECTOR --- > with some vector-valued data

    aol::MultiVector<double> ExampleDataMV(2, N*N);

    double y;

    for (int k=0; k<N*N; ++k)
    {
      // x = (k/N)/static_cast<double> (N);
      y = (k%N)/static_cast<double> (N);
      ExampleDataMV[0][k] = 1.;
      ExampleDataMV[1][k] = sin(6.*y);
    }
    // ------------ End of data definition


    /* **********************************************************
    * now the things you really need for starting GRAPE
    * ********************************************************** */

    // get this mesh back
    GENMESH2D* gmeshNew;

    // **************** call the convert-function *********************
    // this function expects the data defined on a Quocmesh
    // and provides an adaptive GMesh, which can be handled by GRAPE
    // you can call the convert function with scalar or with vector-
    // valued data:
    // furthermore you can add scalar oder vector-valued data

//     gmeshNew = quocmesh_convert_to_gmesh2d(&ExampleDataMV, "Vektordaten");

    gmeshNew = quocmesh_convert_to_gmesh2d(&ExampleDataSc, &ErrorArray, "Skalare Daten mit Estimator");
//     addVectorData(gmeshNew, &ExampleDataMV, "Vektordaten");
//     addScalarDataWithEstimator(gmeshNew, &ExampleDataSc, "Skalare Daten mit Est", &ErrorArray);

    // ************** now start GRAPE ************
    initStartGrape(gmeshNew, "Skalare Daten");



  }
  catch(aol::Exception e)
  {
    e.dump();
    return 42;
  }

  // thats it - try and enjoy it :-)
  exit(0);
}

#else

int main ( int, char** ) {
  cerr << "Without grape external, this program is useless" << endl;
  return ( 0 ) ;
}

#endif
