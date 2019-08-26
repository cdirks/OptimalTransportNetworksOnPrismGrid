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
 *  \brief using valgrind
 *
 *  valgrind --leak-check=full --show-reachable=yes ./valgrindExample
 *  Program with various bugs that are detected by valgrind.
 *
 *  \author Schwen
 */

#include <aol.h>

WARNING_OFF(uninitialized) // this program uses an uninitialized variable to demonstrate debugging

int main ( int, char**) {

  int *values = new int[10];   // potential memory leak?

  for ( int i = 0; i < 10; ++i )
    values[i] = i;             // legal write access

  int valuem1 = values[-1];    // illegal read access
  cerr << valuem1 << endl;

  int value;
  values[5] = value;           // this is an illegal read access not yet complained about by valgrind ...
  cerr << value << endl;       // read access to uninitialized variable
  cerr << values[5] << endl;   // ... now valgrind complains


  values[10] = 42;             // illegal write access
  cerr << values[10] << endl;

  if ( 1 == 0 )
    delete[] ( values );       // memory leak!

}
