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
#include <ctrlCCatcher.h>

namespace aol {

static bool ctrlC_pressed = false;

void ctrlCHandler ( int sig ) {
  signal(sig, ctrlCHandler);
  ctrlC_pressed = true;
}

bool getCtrlCState() {
  return ctrlC_pressed;
}

void resetCtrlCState() {
  ctrlC_pressed = false;
}

bool askForAction ( string action, string color ) {

  cerr << color << action << " (y/n, Ctrl-C to terminate)? " << aol::color::reset;
  string choice;
  sigfunc currentHandler = signal ( InterruptSignal, DefaultHandler );
  cin >> choice;
  signal ( InterruptSignal, currentHandler );

  return ( choice == "y" || choice == "yes" || choice == "j" || choice == "ja" );
}

} // end of namespace aol.
