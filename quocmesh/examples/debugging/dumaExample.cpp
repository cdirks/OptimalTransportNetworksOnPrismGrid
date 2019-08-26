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

/** @file
 *
 *  @brief Adaption of efenceExample for duma
 *
 *  Usage:  <pre>/home/prog/duma-{32,64}/bin/duma dumaExample</pre>
 *  or: add
 *  \verbatim
     CFLAGS += -DUSE_DUMA -I$(DUMADIR)/include/
     LFLAGS += -L$(DUMADIR)/lib/ -lduma
     \endverbatim
 *  to <pre>makefile.local</pre>
 *  or: use <pre>util/configure GNUDUMADEBUG</pre>
 *  If you get an error stating <i>mprotect()</i> failed, try running program with <pre>DUMA_PROTECT_FREE=0.</pre>
 *  Further information can be found <a href="/home/prog/duma-64/share/doc/duma/README.txt">here</a>
 *
 *  @author Geihe
 */

#include <aol.h>

using namespace std;

int main ( int, char**) {

  {
    int values[10];
    cerr << values[-1] << endl; // we've turned off -Warray-bounds, else compiler would already produce warning here
    cerr << values[10] << endl; // duma does not detect these errors
    cerr << values[11] << endl;
    cerr << values[12] << endl << endl;
  }


  {
    int *values = new int[15];  // 15 = 1111

    cerr << values[-1] << endl; // duma only detects this error if you set environment variable DUMA_PROTECT_BELOW=1
    cerr << values[15] << endl; // duma detects this error

//     free ( values );         // and this (free must not be used with new)
//     delete ( values );       // and this (scalar version of delete is wrong here)
    delete [] (values);

    int *leak = new int[10];    // duma also detects memory leaks
    cerr << leak[0] << endl;
  }

  return 0;
}
