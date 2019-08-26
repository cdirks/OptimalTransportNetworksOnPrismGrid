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

#include <gnuplotter.h>

namespace aol {

bool runGnuplot ( const char *GnuplotCommandFileName ) {
  string systemCommand;
  // Since the Windows version of gnuplot finally also has an exe called "gnuplot" the ifdef
  // shouldn't be necessary anymore, it won't hurt to keep it here for now though.
#ifdef WIN32
  systemCommand = "gnuplot.exe ";
#else
  systemCommand = "gnuplot ";
#endif
  systemCommand += GnuplotCommandFileName;
  const bool failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
  if ( failed )
    cerr << "aol::runGnuplot: Calling gnuplot returned an error.\n";
  return !failed;
}

}
