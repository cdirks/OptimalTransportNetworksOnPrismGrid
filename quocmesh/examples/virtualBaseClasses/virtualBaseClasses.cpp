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

#include <aol.h>

class A{
public:
  A () {cerr << "A Default" << endl;}
  A ( int i ) {cerr << "A "<< i << endl;}
};

class B : public virtual A{
public:
  B ()
    : A(0){cerr << "B Default" << endl;}
  B ( int i )
    : A(i) {cerr << "B "<< i << endl;}
};

class C : public B{
public:
  C ( int i )
    : B(i) {cerr << "C "<< i << endl;}
};

int main ( int, char** ) {

  try {
    C c(10);
  } catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
