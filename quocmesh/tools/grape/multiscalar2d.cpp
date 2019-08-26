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
//! @brief Load multiple scalar functions living on the same 2d grid in GRAPE as different functions of the same mesh.
//!
//! Usage: \code ./multiscalar2d data [data...] \endcode
//!
//! The functions need to be defined on the same grid.
//! GRAPE is started with one 2d mesh object having several scalar functions, named in GRAPE according to the file names.
//!
//! @author Lenz

#include "grapeInterface2d.h"
#include <qmException.h>
#include <fstream>
#include <vector>
using namespace std;

#ifdef USE_EXTERNAL_GRAPE


//! Read multiple scalar arrays into grape
//! You can simply call it a la <<./multiscalar *.a2d>>
//! Functions are named according to filenames
int main (int argc, char** argv)
{
  try {
    if (argc == 1) {
      cerr << "usage: " << argv [0] << " data [data...]" << endl;
      return 23;
    }

    qc::ScalarArray<double, qc::QC_2D>* data;
    data = new qc::ScalarArray<double, qc::QC_2D> (argv[1]);

    int N = data->getNumX ();

    GENMESH2D* mesh = quocmesh_convert_to_gmesh2d (data, argv [1]);

    for (int i = 2; i < argc; ++i) {

      data = new qc::ScalarArray<double, qc::QC_2D> (argv[i]);
      if (data->getNumX () != N || data->getNumY () != N) throw aol::Exception ("array size not equal", __FILE__, __LINE__);
      addScalarData (mesh, data, argv [i]);
    }

    initStartGrape(mesh);


  }
  catch (aol::Exception e) {
    e.dump ();
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
