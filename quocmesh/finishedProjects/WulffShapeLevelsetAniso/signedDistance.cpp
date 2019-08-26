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

#include <configurators.h>
#include <eikonalNA.h>
#include <narrow.h>
#include <quoc.h>
#include <scalarArray.h>
#include <signedDistanceOp.h>

#define REAL float

int main( int argc, char **argv)
{
  try
    {
      if ( argc < 3 )
  {
   cerr << "Computes the signed distance function to a level-set (default is 0) in 3D.\n";
    cerr << "usage: signedDistance input output\n";
    return EXIT_FAILURE;
  }
      qc::ScalarArray<REAL, qc::QC_3D> inputimg( argv[1]);
      int N = inputimg.getNumX ();
      int d = qc::logBaseTwo (N);
      qc::GridDefinition grid( d, qc::QC_3D );
      qc::ScalarArray<REAL, qc::QC_3D> outputimg(grid);
//       inputimg -= 0.02;
      qc::SignedDistanceOp3D<REAL>( grid ).apply( inputimg, outputimg );

      //outputimg -= outputimg.getMinValue();

      outputimg.save( argv[2], qc::PGM_DOUBLE_BINARY );
    }
  catch ( aol::Exception &ex )
    {
      ex.dump();
    }
}
