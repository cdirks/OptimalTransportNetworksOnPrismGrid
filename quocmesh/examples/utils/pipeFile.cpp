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
 *  \brief piping an open file through an application
 *
 *   This example illustrates how to open a file and pipe it through an application,
 *   e.g. how to unzip the file before using it.
 *
 *  \author Droske (comment by Wirth)
 */

#include <iostream>
#include <ext/stdio_filebuf.h>

using namespace std;
using namespace __gnu_cxx;

int main ( void ) {
  // open a bz2-file, unzip it and assign it to an istream
  FILE *IN = popen ( "bunzip2 -c ../testdata/volume_17.dat.bz2", "r" );
  stdio_filebuf<char> buf ( IN, ios::in );
  istream in ( &buf );

  // read in and output the first 7 characters of the file
  char c;
  for ( int i = 0; i < 7; i++ ) {
    in >> c;
    cerr << c << endl;
  }

  pclose ( IN );
}
