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
 *  \brief generates and modifies 2D level arrays
 *
 *  The program is called with one or more parameter files as console argument(s).
 *  Each parameter file has to specify an "action" from the following list:
 *
 *     setConst
 *     addConst
 *     linComb
 *     maxOfTwo
 *     minOfTwo
 *     product
 *     smooth
 *     redistance
 *
 *  and an "output" file. Depending on the selected action, further parameters are
 *  read from the parameter file:
 *
 *  (a) setConst
 *
 *      requires "value" and "gridLevel" and produces a constant image with this value
 *
 *  (b) addConst
 *
 *      requires "input" and "value" and adds the value to the given image
 *
 *  (c) linComb, maxOfTwo, minOfTwo, product
 *
 *      They all require "input1" and "input2", linComb additionally searches the
 *      parameter file for "coeff1" and "coeff2". They write into the output file
 *      what you would expect from their names.
 *
 *  (d) smooth
 *
 *      needs "input" data array and smoothing parameter "sigma". With higher sigma,
 *      smoothing is stronger.
 *
 *  (e) redistance
 *
 *      writes a level array to "output" that has the same zero-level-set as
 *      "input", but is a discretized signed distance function
 *
 *  \author Stefan W. von Deylen
 */

#include <generator.h>
#include <configurators.h>
#include <parameterParser.h>
#include <linearSmoothOp.h>
#include <signedDistanceOp.h>

using namespace aol;
using namespace qc;

#include <cstdlib>
#include <iostream>

using namespace std;

//------------------------------------------------------------------------------------------------

typedef long double                                          RealType;

typedef QuocConfiguratorTraitMultiLin<RealType,
                                      QC_2D,
                                      GaussQuadrature
                                      <RealType, QC_2D, 3> > ConfType;

typedef FastUniformGridMatrix<RealType, ConfType::Dim>       MatrixType;
typedef ConfType::ArrayType                                  ArrayType;

const SaveType saveType = PGM_DOUBLE_BINARY;

//------------------------------------------------------------------------------------------------

#include "dataGenerator2D3D.h"

//------------------------------------------------------------------------------------------------

int main( int argc, char ** argv )
{
  if (argc < 2) {
    cerr << aol::color::error << "Usage: " << argv[0] << " <parameter files>" << endl << aol::color::reset;
    return -1;
  }

  try {
    for ( int i = 1; i < argc; ++i ) {
      aol::ParameterParser parser ( argv[i] );
      string action = parser.getString ( "action" );

      clog << aol::color::blue;

      if ( action == "setConst" )                     setConst ( parser );
      else if ( action == "addConst" )                addConst ( parser );
      else if ( action == "linComb" )                 linComb ( parser );
      else if ( action == "maxOfTwo" )                maxOfTwo ( parser );
      else if ( action == "minOfTwo" )                minOfTwo ( parser );
      else if ( action == "product" )                 product ( parser );
      else if ( action == "smooth" )                  smooth ( parser );
      else if ( action == "redistance" )              redistance<ConfType> ( parser );

      clog << aol::color::reset;
    }
  }
  catch (Exception & exc) {
    exc.dump();
  }
  catch (...) {
    clog << aol::color::reset;
  }

  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
//------------------------------------------------------------------------------------------------
