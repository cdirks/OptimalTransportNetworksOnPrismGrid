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

/* *******************************************************************
 * Plot tree structure and flow velocities to eps file
 * ******************************************************************* */

#include "grow_vesseltree.h"

int main ( int argc, char** argv ) {
  try {
    if ( argc == 5 ) {
      grow_vesseltree<float> my_vesseltree;
      my_vesseltree.read_trees ( argv[1], argv[2] );

      my_vesseltree.gatree.arcs = my_vesseltree.atree.arcs;
      my_vesseltree.gatree.nodes = my_vesseltree.atree.nodes;
      my_vesseltree.gvtree.arcs = my_vesseltree.vtree.arcs;
      my_vesseltree.gvtree.nodes = my_vesseltree.vtree.nodes;

      my_vesseltree.epsout_grow_trees_radii ( argv[3] ); // this only works on grow_vesseltrees

      my_vesseltree.plot_tree_speeds ( argv[4], 1. / 20000., 1. );

    } else {
      cerr << "usage: usevesseltree_epsplot <a.tree> <v.tree> <outfile.eps> <outfile_speeds.eps>" << endl;
    }
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return ( 0 );
}
