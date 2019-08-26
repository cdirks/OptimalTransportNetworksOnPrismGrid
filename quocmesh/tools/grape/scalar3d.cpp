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
//! @brief Visualize one scalar function with 3d domain in grape.
//!
//! Usage: \code ./scalar3d datafile \endcode
//!
//! GRAPE is started with one 3d mesh object having one scalar function.
//!
//! @author Lenz

#include "grapeInterface3d.h"
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

    qc::ScalarArray<double, qc::QC_3D> data (argv [1]);

    int N = data.getNumX ();
    int d = qc::logBaseTwo (N);
    if (data.getNumY () != N || data.getNumZ () != N) throw aol::Exception ("Image not cube", __FILE__, __LINE__);
    if ((1 << d) + 1 != N) throw aol::Exception ("dimension no power of two", __FILE__, __LINE__);

    //    qc::GridDefinition grid (d, qc::QC_3D); // this is not used.
    GENMESH3D* mesh = quocmesh_convert_to_gmesh3d (&data, "Skalare Daten");

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
