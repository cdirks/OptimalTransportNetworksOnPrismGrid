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
*   \brief Loads 3 scalarArrays, the second and third interpreted as angles for unit vector fields.
*
* Loads 3 scalarArrays: the first one is a usual image, the second
* and the third one are regarded as angles, from which two vector fields
* are generated.
*
* \author Nemitz
*
* Date: 04.04.2007
* ********************************************************************** */

#include "grapeInterface2d.h"
#include <aol.h>
#include <scalarArray.h>

#ifdef USE_EXTERNAL_GRAPE

using namespace aol::color;

int main( int argc, char** argv )
{

  try
  {
    if (argc != 4) {
      cerr << "Loads 3 scalarArrays: the first one is a usual image, the second\n";
      cerr << "and the third one are regarded as angles, from which two vector fields\n";
      cerr << "are generated.\n";
      cerr << red << "usage: " << argv [0] << " image angle1 angle2" << endl << reset;
      return 23;
    }

    // load data
    qc::ScalarArray<double, qc::QC_2D> img( argv[1] );
    qc::ScalarArray<double, qc::QC_2D> alpha( argv[2] );
    qc::ScalarArray<double, qc::QC_2D>  beta( argv[3] );

    // generate the two vector fields from the angles

    aol::MultiVector<double> alphaVF(2, alpha.size() );
    aol::MultiVector<double> betaVF(2, beta.size() );

    for ( int k=0; k<alpha.size(); k++ )
    {
      alphaVF[0][k] = cos( alpha[k] );
      alphaVF[1][k] = sin( alpha[k] );
      betaVF[0][k] = cos( beta[k] );
      betaVF[1][k] = sin( beta[k] );
    }


    /* **********************************************************
    * now the things you really need for starting GRAPE
    * ********************************************************** */

    // get this mesh back
    GENMESH2D* gmeshNew;

    gmeshNew = quocmesh_convert_to_gmesh2d(&img, "Img");
    addVectorData(gmeshNew, &alphaVF, "alpha");
    addVectorData(gmeshNew, &betaVF, "beta");

    // ************** now start GRAPE ************
    initStartGrape(gmeshNew, "Img");

  }
  catch(aol::Exception e)
  {
    e.dump();
    return 42;
  }

  // thats it - try and enjoy it :-)
  exit(0);
}

#else

int main ( int, char** ) {
  cerr << "Without grape external, this program is useless" << endl;
  return ( 0 ) ;
}

#endif
