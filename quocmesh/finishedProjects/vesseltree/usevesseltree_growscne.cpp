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
 * Generate vessel tree using CCO
 * ******************************************************************* */

#include "grow_vesseltree.h"

int main() {
  try {
    aol::ParameterParser parser ( "par/growscne.par" );

    for ( int seed = parser.getInt ( "min_seed" ); seed < parser.getInt ( "max_seed" ); seed++ ) {

      cerr << "Random seed = " << seed << endl;
      srand ( seed );
      grow_vesseltree<float> my_vesseltree;

      cerr.precision ( 10 );
      aol::StopWatch timer;
      timer.start();
      my_vesseltree.grow_randomtrees_scne ( seed );
      timer.stop();
      cerr << "Generating tree took " << timer.elapsedCpuTime() << " seconds." << endl;

      my_vesseltree.save_grow_trees ( seed );
      my_vesseltree.epsout_grow_trees_radii ( seed );

    }
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return ( 0 );
}
