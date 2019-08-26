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
//! @brief Visualize one scalar function with 2d domain in grape.
//!
//! Usage: \code ./scalar2d datafile \endcode
//!
//! GRAPE is started with one 2d mesh object having one scalar function.
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
    if (argc != 2) {
      cerr << "usage: " << argv [0] << " data" << endl;
      return 23;
    }

    // load data into scalar array
    qc::ScalarArray<double, qc::QC_2D> data (argv [1]);

    // convert this to a genmesh (it will automatically be tested
    // whether the data is quadratic or not)
    GENMESH2D* mesh = quocmesh_convert_to_gmesh2d (&data, "Skalare Daten");


    add_reload_button( mesh, data, argv [1] );

    // and then start GRAPE, thats it!
    initStartGrape(mesh);
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
