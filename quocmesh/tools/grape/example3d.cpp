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
//! @brief Example for using the GRAPE interface for a 3d scalar array including error estimators and a MultiVector for vector valued data.
//!
//! This file is an example for handling data of a 3D-QuocMesh with
//! the grape-interface. It handles one scalar and one
//! vector-valued function, which must be in an qcScalarArray and in
//! a qcMultiVector.
//!
//! @author Nemitz

#include "grapeInterface3d.h"
#include <qmException.h>
#include <aol.h>

#ifdef USE_EXTERNAL_GRAPE

int main(void)
{
  try
  {
    // the depth of the grid
    int d=6;
    int N = (1<<d)+1;      // 2^N + 1


    /* **********************************************************
    * now define some example data - first scalar-values
    * ********************************************************** */

    qc::ScalarArray<double, qc::QC_3D> ExampleDataScalar(N,N,N);
    cerr<<"Berechne skalare Daten...";

    for (int i=0; i<N; ++i)
    {
      cerr<<".";
      for (int j=0; j<N; ++j)
        for (int k=0; k<N; ++k)
        {
          double x = static_cast<double> (i) / static_cast<double> (N);
          double y = static_cast<double> (j) / static_cast<double> (N);
          double z = static_cast<double> (k) / static_cast<double> (N);

          ExampleDataScalar.set(i,j,k, (x*x + y*y + z*z)/3.);
        }
    }


    // *****************************************************
    // now the error estimator for the scalar data
    // ********************************************************
    qc::Estimator3d<double> ErrorArray(N);
    cerr<<"saturiere...";
    ErrorArray.makeSaturatedErrorArray(&ExampleDataScalar);


    // fill a < --- MULTIVECTOR --- > with some vector-valued data

    aol::MultiVector<double> ExampleDataMV(3,N*N*N);

    double x,y;

    cerr<<"Berechne Vektordaten...";

    for (int k=0; k<N*N*N; ++k)
    {
      if (k % (N*N) == 0) cerr<<".";
      x = (k/N) / static_cast<double> (N);
      y = (k%N) / static_cast<double> (N);
      ExampleDataMV[0][k] = 0.01*(0.5-x);
      ExampleDataMV[1][k] = 0.5*sin(3.*y);
      ExampleDataMV[2][k] = 0.2;
    }

    cerr<<"fertig! Dimension des Vektors: "<<(ExampleDataMV.numComponents())<<endl;


    // **************** call the convert-function *********************
    // this function expects the data
    // and provides an adaptive GMesh, which can be handled by GRAPE

    GENMESH3D* gmeshNew = quocmesh_convert_to_gmesh3d(&ExampleDataMV, "Vektordaten");
    addScalarDataWithEstimator(gmeshNew, &ExampleDataScalar, "Scalar with Est", &ErrorArray);

    // ************** now start GRAPE ************
    initStartGrape(gmeshNew, "Skalare Daten");

    // thats it - try and enjoy it :-)
  }
  catch(aol::Exception e)
  {
    e.dump();
    return 42;
  }

  return 0;
}

#else

int main ( int, char** ) {
  cerr << "Without grape external, this program is useless" << endl;
  return ( 0 ) ;
}

#endif
