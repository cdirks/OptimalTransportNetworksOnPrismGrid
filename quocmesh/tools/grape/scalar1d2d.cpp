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
//! @brief Load one scalar function with 1d domain and one scalar function with 2d domain into GRAPE as two different objects.
//!
//! Usage: \code ./scalar1d2d data1d data2d \endcode
//!
//! GRAPE is started with one 1d triang object and one 2d mesh object each having one scalar function.
//!
//! @author Lenz

#include "grapeInterface1d.h"
#include "grapeInterface2d.h"
#include <qmException.h>
#include <fstream>
using namespace std;

#ifdef USE_EXTERNAL_GRAPE

int main (int argc, char** argv)
{
  try {
    if (argc != 3) {
      cerr << "usage: " << argv [0] << " data1d data2d" << endl;
      return 23;
    }

    // load data into scalar arrays
    qc::ScalarArray<double, qc::QC_1D> data1d (argv [1]);
    qc::ScalarArray<double, qc::QC_2D> data2d (argv [2]);

    // convert these to a genmesh (it will automatically be tested
    // whether the data is quadratic or not)
    GENMESH2D* mesh2d = quocmesh_convert_to_gmesh2d (&data2d, "mesh2d");
    TRIANG1D*  triang = quocmesh_convert_to_triang1d (&data1d, "triang1d");

    // Connect objects via scenes
    SCENE* scene1d = reinterpret_cast<SCENE*> (GRAPE (Scene, "new-instance") ("scene1d"));
    SCENE* scene2d = reinterpret_cast<SCENE*> (GRAPE (Scene, "new-instance") ("scene2d"));
    ASSIGN (scene1d->object, reinterpret_cast<TREEOBJECT*> (triang));
    ASSIGN (scene2d->object, reinterpret_cast<TREEOBJECT*> (mesh2d));
    ASSIGN (scene1d->next_scene, scene2d);

    // and then start GRAPE, thats it!
    addMethodsAndProjects1d();
    addMethodsAndProjects2d();
    GRAPE (GRAPE (Manager, "get-stdmgr") (), "handle") (scene1d);
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
