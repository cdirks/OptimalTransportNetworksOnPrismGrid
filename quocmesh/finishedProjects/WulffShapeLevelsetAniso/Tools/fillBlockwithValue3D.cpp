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
#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>
#include <parameterParser.h>

#define RealType double

using namespace aol::color;

// *************************************************************************
// * fillBlockWithValue3D.cpp
// * fills a cubic block of a scalarArray3d with a fixed value.
// * Args: parameter-file.
// *************************************************************************

int main( int argc, char **argv)
{
  try
    {
      if ( argc < 2 )
  {
   cerr << "Fills a cubic block of a scalarArray3d with a fixed value.\n";
    cerr <<aol::color::red<< "usage: "<<argv[0]<<" parameterfile\n"<<aol::color::black;
    return EXIT_FAILURE;
  }

      aol::ParameterParser parser( argv[1] );

      // read coordinates of the block
      int x1 = parser.getInt("x1");
      int x2 = parser.getInt("x2");
      int y1 = parser.getInt("y1");
      int y2 = parser.getInt("y2");
      int z1 = parser.getInt("z1");
      int z2 = parser.getInt("z2");

      // read the fill-value
      double value = parser.getDouble("fillValue");

      // load the original image
      char fileName[1024];
      parser.getString( "inputImage", fileName );
      cerr <<green <<"\nLoading file '" << black << fileName << green << "'...\n\n";
      qc::ScalarArray<double, qc::QC_3D> inputImg(fileName);


      // check the coordinates
      if (x2 > inputImg.getNumX() || y2 > inputImg.getNumY() || z2 > inputImg.getNumZ() )
      {
        cerr << red << "\nERROR: Block-Coordinates larger than the ScalarArray!!\n" << reset;
        exit(1);
      }
      if (x1 < 0 || y1 < 0 || z1 < 0)
      {
        cerr << red << "\nERROR: Block-Coordinates are negative!!\n" << reset;
        exit(1);
      }

      // now traverse the block and set all nodes to the value
      cerr << blue << "\nFilling the block...";
      for (int i=x1; i<x2; i++) {
        cerr << ".";
        for (int j=y1; j<y2; j++)
          for (int k=z1; k<z2; k++)
            inputImg.set(i,j,k, value);
      }

      // finally save the modified image
      parser.getString( "outputImage", fileName );
      cerr <<green <<"\nSaving to file '" << black << fileName << green << "'...\n\n";
      inputImg.save( fileName );

      cerr << blue << "\ndone! Thanx for using fillBlockWithValue3D (Non-registered version!!)\n";
      cerr << reset;
    }
  catch ( aol::Exception &ex )
    {
      ex.dump();
    }
}
