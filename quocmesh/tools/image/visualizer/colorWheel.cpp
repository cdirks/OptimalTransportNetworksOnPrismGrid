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

/**
 * \file
 * \brief Adds a color/gray wheel in the lower right part of an image.
 *
 * \author Berkels, Drokse
 */

#include "colorWheel.h"

typedef double RType;

int main ( int argc, char **argv ) {
  if ( argc < 4 ) {
    cerr << "USAGE: " << argv[0] << " filename angle1 angle2 [rescale]" << endl;
    return EXIT_FAILURE;
  }
  qc::ScalarArray<RType, qc::QC_2D> in ( argv[1] );
  qc::ColorWheel<RType> colorWheel ( in, atof ( argv[2] ), atof ( argv[3] ), ( argc >= 5 ) );
  colorWheel.saveColoredWheel ( "output.ppm" );
  colorWheel.saveGrayWheel ( "output.pgm" );
}
