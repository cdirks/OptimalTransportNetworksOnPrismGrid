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
#include <qmException.h>
#include <aol.h>
#include <quoc.h>
#include <parameterParser.h>
#include "grapeInterface3d.h"
// #include "narrow.h"
#include <signedDistanceOp.h>

using namespace std;


double sqr(double x) { return x*x; }


// -------------------- the generating routines ----------------------------------------

void generateTorusField(qc::ScalarArray<double, qc::QC_3D> &feldX,
                        qc::ScalarArray<double, qc::QC_3D> &feldY,
                        qc::ScalarArray<double, qc::QC_3D> &feldZ, int N )
{
  // the epicenter in the x,z-plane
  int Mx = N/2;
  int My = N/2;

  for (int i=0; i<N; ++i)
   for (int j=0; j<N; ++j)
     for (int k=0; k<N; ++k)
     {
       // all vectors are in the x,y-plane => Z-coord = 0
       feldZ.set( i,j,k, 0. );

       // the other vectors are orthogonal to the vector linking
       // the epicenter in the x,z-plane to the actual point
       // all vectors have length = 1
       double d = sqrt( sqr(i-Mx) + sqr(My-j) );

       if ( d > 3 )
       {
         feldX.set( i,j,k, (My-j) / d );
         feldY.set( i,j,k, (i-Mx) / d );
       }
       else
       {
         feldX.set( i,j,k, 1. );
         feldY.set( i,j,k, 0. );
       }
     }
}


void generateConstantField(qc::ScalarArray<double, qc::QC_3D> &feldX,
                        qc::ScalarArray<double, qc::QC_3D> &feldY,
                        qc::ScalarArray<double, qc::QC_3D> &feldZ, int /*N*/,
                        double x, double y, double z )
{
  feldX.setAll( x );
  feldY.setAll( y );
  feldZ.setAll( z );
}


// -------------------------- main program ---------------------------------------------

int main(int argc, char **argv)
{
  if ( argc != 3 ) {
    cerr << "USAGE: " << argv[0] << " depth savename" << endl;
    return EXIT_FAILURE;
  }
  try
  {
    // the depth of the grid
    int d = atoi( argv[1] );
    int N = (1<<d)+1;              // 2^N + 1

    qc::ScalarArray<double, qc::QC_3D> feldX( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> feldY( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> feldZ( N,N,N );

    // generate the vector field
    cerr<<aol::color::green<<"Calculating Vektorfield...\n"<<aol::color::black;

    // the generation calls
//     generateTorusField( feldX, feldY, feldZ, N );
    generateConstantField( feldX, feldY, feldZ, N, 0.,0.,1. );


    // save it
    char filename[ 1024 ];
    sprintf( filename, "%sX.bz2", argv[2]);
    feldX.save( filename, qc::PGM_DOUBLE_BINARY );
    sprintf( filename, "%sY.bz2", argv[2]);
    feldY.save( filename, qc::PGM_DOUBLE_BINARY );
    sprintf( filename, "%sZ.bz2", argv[2]);
    feldZ.save( filename, qc::PGM_DOUBLE_BINARY );
    cerr<<aol::color::green<<"... sucessfully finished! \n\n"<<aol::color::black;


    // and view it in Grape
    aol::MultiVector<double> VF (0, feldX.getNumX() * feldY.getNumY() * feldZ.getNumZ());
    VF.appendReference (feldX);
    VF.appendReference (feldY);
    VF.appendReference (feldZ);

#ifdef USE_EXTERNAL_GRAPE
    GENMESH3D* mesh = quocmesh_convert_to_gmesh3d (&VF, "VectorField");
    // and then start GRAPE, thats it!
    initStartGrape(mesh, "VectorField");
#else
    cerr << "cannot display in GRAPE without grape external" << endl;
#endif

    // thats it - try and enjoy it :-)
  }
  catch(aol::Exception e)
  {
    e.dump();
    return 42;
  }

  return 0;
}

