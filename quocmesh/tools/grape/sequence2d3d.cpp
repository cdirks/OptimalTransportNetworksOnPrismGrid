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
//! @brief Loads two sequences of scalar functions, one with 2d, one with 3d domain.
//!
//! Usage: \code ./sequence2d3d any number of 2d and 3d files \endcode
//!
//! Each file is added as the next time step either to a 2d or to a 3d mesh object,
//! depending on the file header with gives the dimension of the domain.
//! Both sequences are asumed to have times 0, 1, 2, 3, ...
//! All 2d files must share the same grid, the same is true for the 3d files.
//! GRAPE is started with two time dependent mesh object each having one scalar function.
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
    if (argc < 2) {
      cerr << "usage: " << argv [0] << " any number of 2d and 3d files" << endl;
      return 23;
    }

    addMethodsAndProjects2d ();
    addMethodsAndProjects3d ();

    GENMESH2D* mesh2d = NULL;
    GENMESH3D* mesh3d = NULL;

    // for all filenames given as parameters
    for (int i = 1; i < argc; ++i) {

      qc::ArrayHeader header;
      {
  aol::Bzipifstream file ( argv [i] );
  qc::ReadArrayHeader ( file, header );
      }

      if ( header.magic[0] == 'P' ) { // 2d

  qc::ScalarArray<double, qc::QC_2D>* data = new qc::ScalarArray<double, qc::QC_2D> (argv [i]);

  if ( mesh2d == NULL ) { // first 2d timestep
    mesh2d = quocmesh_convert_to_gmesh2d ( data, "scalar" );
  }
  else { // there are already other 2d timesteps
    addTimestep ( mesh2d, data, "scalar" );
  }
      }

      if ( header.magic[0] == 'Q' ) { // 3d

  qc::ScalarArray<double, qc::QC_3D>* data = new qc::ScalarArray<double, qc::QC_3D> (argv [i]);

  if ( mesh3d == NULL ) { // first 3d timestep
    mesh3d = quocmesh_convert_to_gmesh3d ( data, "scalar" );
  }
  else { // there are already other 3d timesteps
    addTimestep ( mesh3d, data, "scalar" );
  }
      }
    }

    // Connect objects via scenes
    TIMESCENE* firstscene = NULL;
    if (mesh2d) {
      TIMESCENE* scene2d = reinterpret_cast<TIMESCENE*> (GRAPE (TimeScene, "new-instance") ("sequence2d"));
      ASSIGN (scene2d->dynamic, reinterpret_cast<TREEOBJECT*> (mesh2d));
      ASSIGN (scene2d->object,  reinterpret_cast<TREEOBJECT*> (GRAPE (mesh2d, "softcopy") (NULL)));
      GRAPE (scene2d, "switch-sync-send") ();
      firstscene = scene2d;
    }
    if (mesh3d) {
      TIMESCENE* scene3d = reinterpret_cast<TIMESCENE*> (GRAPE (TimeScene, "new-instance") ("sequence3d"));
      ASSIGN (scene3d->dynamic, reinterpret_cast<TREEOBJECT*> (mesh3d));
      ASSIGN (scene3d->object,  reinterpret_cast<TREEOBJECT*> (GRAPE (mesh3d, "softcopy") (NULL)));
      GRAPE (scene3d, "switch-sync-send") ();
      if (firstscene) {  ASSIGN (firstscene->next_scene, reinterpret_cast<SCENE*> (scene3d)); }
      else { firstscene = scene3d; }
    }

    // and then start GRAPE, thats it!
    GRAPE (GRAPE (Manager, "get-stdmgr") (), "handle") (firstscene);
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
