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
//! @brief Load one vector valued function with a 2d domain and range by specifying two scalar functions.
//!
//! Usage: \code ./vector2d x-component y-component \endcode
//!
//! The 2 components need to be defined on the same grid, and are interpreted as two components of a vector valued function.
//! GRAPE is started with one 2d mesh object having one vector valued function.
//!
//! @author Lenz

#include "grapeInterface2d.h"
#include <qmException.h>
#include <fstream>
using namespace std;

#ifdef USE_EXTERNAL_GRAPE

int main (int argc, char** argv)
{
  try {
    if (argc != 3) {
      cerr << "usage: " << argv [0] << " comp1 comp2" << endl;
      return 23;
    }

    // load data
    qc::ScalarArray<double, qc::QC_2D> comp1 (argv [1]), comp2 (argv [2]);

    int N = comp1.getNumX ();

    aol::MultiVector<double> data (0, N * N);
    data.appendReference (comp1);
    data.appendReference (comp2);

    // convert this to a genmesh (it will automatically be tested
    // whether the data is quadratic or not)
    GENMESH2D* mesh = quocmesh_convert_to_gmesh2d (&data, "Vektordaten");

    initStartGrape(mesh);

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
