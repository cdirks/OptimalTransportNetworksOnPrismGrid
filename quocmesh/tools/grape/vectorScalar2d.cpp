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
*   \brief Loads a number of scalar and vector fields in 2D.
*
* Loads a number of scalar and vector fields. The flag -v can be used to indicate
* that the following filenames are the x- and y-components of vector fields, the flag
* -s can be used to indicate that the following filenames are scalar fields.
*
* \author Wirth
*
* Date: 01.07.2008
*/

#include "grapeInterface2d.h"
#include <aol.h>
#include <scalarArray.h>
#include <vectorExtensions.h>

#ifdef USE_EXTERNAL_GRAPE

using namespace aol::color;

int main( int argc, char** argv )
{

  try
  {
    if ( argc < 3 ) {
      cerr << "Loads a number of scalar and vector fields. The flag ""-v"" can be used to indicate\n";
      cerr << "that the following filenames are the x- and y-components of vector fields, the\n";
      cerr << "flag ""-s"" can be used to indicate that the following filenames are scalar fields.\n";
      cerr << "Default is ""-s"".\n";
      cerr << red << "usage: " << argv [0] << " -s scalar_field_1 scalar_field_2 ... -v vector_x_1 vector_y_1 vector_x_2 vector_y_2 ..." << endl << reset;
      return 23;
    }

    // load data
    aol::RandomAccessContainer<qc::ScalarArray<double,qc::QC_2D> > images;
    aol::RandomAccessContainer<aol::MultiVector<double> > vectors;
    std::vector<int> imageIndex;
    std::vector<int> vectorIndex;
    bool scalarMode = true;
    int argi = 1;
    while ( argi < argc ) {
      if ( ( argv[argi][0] == '-' ) && ( argv[argi][1] == 'v' ) )         // flag for vector data; read in vector fields
        scalarMode = false;
      else if ( ( argv[argi][0] == '-' ) && ( argv[argi][1] == 's' ) )    // flag for scalar data; read in scalar fields
        scalarMode = true;
      else if ( scalarMode ) {                                            // read in a scalar field
        imageIndex.push_back( argi );
        images.pushBack( qc::ScalarArray<double,qc::QC_2D>( argv[argi] ) );
      } else {                                                            // read in a vector field
        vectorIndex.push_back( argi );
        qc::ScalarArray<double,qc::QC_2D> vecX( argv[argi++] );
        qc::ScalarArray<double,qc::QC_2D> vecY( argv[argi] );
        aol::MultiVector<double> vec;
        vec.appendReference( vecX );
        vec.appendReference( vecY );
        vectors.pushBack( vec );
      }
      argi++;
    }

    // check data
    if ( ( imageIndex.size() == 0 ) || ( vectorIndex.size() == 0 ) )
      throw aol::Exception( "Error: At least one scalar and one vector field is required." );
    int size = images[0].size();
    for ( int i = 0; i < images.size(); i++ )
      if ( images[i].size() != size )
        throw aol::Exception( "Error: Scalar fields not all of the same size." );
    for ( int i = 0; i < vectors.size(); i++ )
      for ( int j = 0; j < 2; j++ )
        if ( vectors[i][j].size() != size )
          throw aol::Exception( "Error: Vector fields have to have the same size as scalar fields." );


    /* **********************************************************
    * now the things you really need for starting GRAPE
    * ********************************************************** */

    // get this mesh back
    GENMESH2D* gmeshNew;

    gmeshNew = quocmesh_convert_to_gmesh2d( &(images[0]), argv[imageIndex[0]] );
    for ( int i = 1; i < images.size(); i++ )
      addScalarData( gmeshNew, &(images[i]), argv[imageIndex[i]] );
    for ( int i = 0; i < vectors.size(); i++ )
      addVectorData( gmeshNew, &(vectors[i]), argv[vectorIndex[i]] );

    // ************** now start GRAPE ************
    initStartGrape( gmeshNew, argv[1] );

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
