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

/** \file
 *  \brief using efence for detecting some illegal memory access
 *
 *   usage: ef efenceExample (or use linker option -lefence) Program
 *   with some bugs that are not detected by efence (even though they
 *   should!), showing limitations of efence.
 *
 *  \author Schwen
 */

#include <aol.h>

int main ( int, char**) {

  {
    int values[10];
    cerr << values[-1] << endl; // we've turned off -Warray-bounds, else compiler would already produce warning here
    cerr << values[10] << endl; // efence does not detect these errors
    cerr << values[11] << endl;
    cerr << values[12] << endl;
  }


  {
    int *values = new int[10];

    cerr << values[-1] << endl; // efence does not detect these errors
    cerr << values[10] << endl;
    cerr << values[11] << endl;
    cerr << values[12] << endl; // program only crashes here.

    delete[] ( values );
  }

}
