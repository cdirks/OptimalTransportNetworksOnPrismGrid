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
//! @brief Load a sequence of scalar functions with a 3d domain.
//!
//! Usage: \code ./sequence3d step0 [step1 ...] \endcode
//!
//! The time steps need to be defined on the same grid, and are interpreted as values of one scalar functions at times 0, 1, 2, ...
//! GRAPE is started with one time dependent 3d mesh object having one scalar function.
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
    if (argc < 2) {
      cerr << "usage: " << argv [0] << " step0 [step1 ...]" << endl;
      return 23;
    }

    // load data into scalar array
    qc::ScalarArray<double, qc::QC_3D> data (argv [1]);

    // convert this to a genmesh (it will automatically be tested
    // whether the data is quadratic or not)
    GENMESH3D* mesh = quocmesh_convert_to_gmesh3d (&data, "scalar");


    // For all timesteps
    for (int i = 2; i < argc; ++i) {

      // load data into scalar array and add to time sequence
      qc::ScalarArray<double, qc::QC_3D>* data = new qc::ScalarArray<double, qc::QC_3D> (argv [i]);
      addTimestep (mesh, data, "scalar");
    }

    // and then start GRAPE, thats it!
    initStartGrapeTime(mesh);
  }
  catch(aol::Exception e) {
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
