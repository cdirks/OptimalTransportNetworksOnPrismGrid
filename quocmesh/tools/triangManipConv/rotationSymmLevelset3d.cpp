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
 *  \brief reads a 2D level array and creates a 3D level array output
 *
 *  reads a 2D level array and creates a 3D level array output by rotating
 *  the input values around the axis x = y = 0.5.
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
                                      <RealType, QC_2D, 3> > ConfType2D;

typedef QuocConfiguratorTraitMultiLin<RealType,
                                      QC_3D,
                                      GaussQuadrature
                                      <RealType, QC_3D, 3> > ConfType3D;

typedef ConfType2D::ArrayType                                ArrayType2D;
typedef ConfType3D::ArrayType                                ArrayType3D;

const SaveType saveType = PGM_DOUBLE_BINARY;

//------------------------------------------------------------------------------------------------

RealType getX ( RealType x, RealType y ) {

  x -= 0.5;
  y -= 0.5;
  RealType r = sqrt ( static_cast<RealType> ( x * x + y * y ) );
  return aol::Min ( r + 0.5, aol::ZOTrait<RealType>::one );
}

int main( int argc, char ** argv )
{
  if (argc != 3) {
    cerr << aol::color::error << "Usage: " << argv[0] << " <input file> <output filename>"
         << endl << aol::color::reset;
    return -1;
  }

  try {
    clog << aol::color::blue;

    clog << "Loading " << argv[1] << " ... ";
    ArrayType2D inArray ( argv[1] );
    RealType in_min = inArray.getMinValue (),
             in_max = inArray.getMaxValue (),
             range  = in_max - in_min;
    // scale to [-0.5, 0.5]
    inArray.addToAll ( - in_min );
    inArray /= range;
    inArray.addToAll ( -0.5 );
    ConfType2D::InitType grid2D ( aol::Vec3<int> ( inArray.getNumX(), inArray.getNumY(), 1 ) );
    aol::DiscreteFunctionDefault<ConfType2D> in ( grid2D, inArray );
    clog << "done." << endl;

    clog << "Filling 3D array ... ";
    qc::GridSize<qc::QC_3D> size3D ( inArray.getNumX(), inArray.getNumY(), inArray.getNumY() );
    ArrayType3D out ( size3D );
    qc::RectangularIterator<qc::QC_3D> iter (Vec3<int>(0,0,0),Vec3<int>(size3D[0], size3D[1], size3D[2]));
    RealType h = 1. / out.getNumX();
    for ( ; iter.notAtEnd(); ++iter ) {
      aol::Vec2<RealType> p ( getX ( h * (*iter)[0], h * (*iter)[1] ), h * (*iter)[2] );
      out.set ( *iter, in.evaluateGlobal ( p ) );
    }
    clog << "done." << endl;

    clog << "Saving as " << argv[2] << " ... ";
    out.save ( argv[2], saveType );
    clog << "done." << endl;

    clog << aol::color::reset;
  }
  catch (Exception & exc) {
    exc.dump();
  }
  catch (...) {
  }

  clog << aol::color::reset;
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
//------------------------------------------------------------------------------------------------
