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
//! @brief Load one scalar function with 2d domain and one scalar function with 3d domain into GRAPE as two different objects.
//!
//! Usage: \code ./scalar2d3d data2d data3d \endcode
//!
//! GRAPE is started with one 2d mesh object and one 3d mesh object each having one scalar function.
//!
//! @author Lenz

#include "grapeInterface2d.h"
#include "grapeInterface3d.h"
#include <qmException.h>
#include <fstream>
using namespace std;

#ifdef USE_EXTERNAL_GRAPE

int main (int argc, char** argv)
{
  try {
    if (argc != 3) {
      cerr << "usage: " << argv [0] << " data2d data3d" << endl;
      return 23;
    }

    // load data into scalar arrays
    qc::ScalarArray<double, qc::QC_2D> data2d (argv [1]);
    qc::ScalarArray<double, qc::QC_3D> data3d (argv [2]);

    // convert these to a genmesh (it will automatically be tested
    // whether the data is quadratic or not)
    GENMESH2D* mesh2d = quocmesh_convert_to_gmesh2d (&data2d, "mesh2d");
    GENMESH3D* mesh3d = quocmesh_convert_to_gmesh3d (&data3d, "mesh3d");

    // Connect objects via scenes
    SCENE* scene2d = reinterpret_cast<SCENE*> (GRAPE (Scene, "new-instance") ("scene2d"));
    SCENE* scene3d = reinterpret_cast<SCENE*> (GRAPE (Scene, "new-instance") ("scene3d"));
    ASSIGN (scene2d->object, reinterpret_cast<TREEOBJECT*> (mesh2d));
    ASSIGN (scene3d->object, reinterpret_cast<TREEOBJECT*> (mesh3d));
    ASSIGN (scene2d->next_scene, scene3d);

    // and then start GRAPE, thats it!
    g_project_add (const_cast<char*> ("uif-m2"));
    g_project_add (const_cast<char*> ("uif-m3"));
    GRAPE (GRAPE (Manager, "get-stdmgr") (), "handle") (scene2d);
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
