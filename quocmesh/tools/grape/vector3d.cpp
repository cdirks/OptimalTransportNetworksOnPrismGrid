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
//! @brief Load one vector valued function with a 3d domain and range by specifying three scalar functions.
//!
//! Usage: \code ./vector3d x-component y-component z-component \endcode
//!
//! The 3 components need to be defined on the same grid, and are interpreted as three components of a vector valued function.
//! GRAPE is started with one 3d mesh object having one vector valued function.
//!
//! @author Nemitz

#include "grapeInterface3d.h"
#include <qmException.h>
#include <fstream>
using namespace std;

#ifdef USE_EXTERNAL_GRAPE

int main (int argc, char** argv)
{
  try {
    if (argc != 4) {
      cerr << "usage: " << argv [0] << " comp1 comp2 comp3" << endl;
      return 23;
    }

    // load data
    qc::ScalarArray<double, qc::QC_3D> comp1 (argv [1]), comp2 (argv [2]), comp3 (argv [3]);

    int N = comp1.getNumX ();

    aol::MultiVector<double> data (0, N * N * N);
    data.appendReference (comp1);
    data.appendReference (comp2);
    data.appendReference (comp3);

    // convert this to a genmesh (it will automatically be tested
    // whether the data is quadratic or not)
    GENMESH3D* mesh = quocmesh_convert_to_gmesh3d (&data, "Vektordaten");

    initStartGrape(mesh, "Vektordaten");

  }
  catch(aol::Exception e) {
    e.dump();
    return 42;
  }

  exit(0);
}

#else

int main ( int, char** ) {
  cerr << "Without grape external, this program is useless" << endl;
  return ( 0 ) ;
}

#endif
