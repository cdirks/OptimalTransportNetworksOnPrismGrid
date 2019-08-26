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

//! @file
//!
//! @brief Load multiple scalar functions living on the same 2d or 3d grid in GRAPE as different functions of the same mesh.
//!
//! Usage: \code ./grapefruit data [data...] \endcode
//!
//! The functions need to be defined on the same grid (i.e. is also has to have the same dimension of domain).
//! GRAPE is started with one 2d or 3d mesh object having several scalar functions, named in GRAPE according to the file names.
//!
//! @author Nemitz

#include "grapeInterface2d.h"
#include "grapeInterface3d.h"
#include <qmException.h>
#include <fstream>
#include <quoc.h>
#include <auxiliary.h>

#ifdef USE_EXTERNAL_GRAPE


using namespace std;
using namespace aol::color;



// get the first character of the header
char getFirstHeaderChar( char* fileName )
{
  aol::ipfstream file ( fileName );
  qc::ArrayHeader header;
  qc::ReadArrayHeader (file, header);
  file.close();

  return header.magic[0];
}

// exit if there are 2D AND 3D data-sets
void exitMixedError()
{
  cerr << red << "\nERROR: Unfortunately it is NOT possible to load 2D and 3D data at the same time!\n" << reset;
  exit(42);
}

// exit if the image has not the same resolution in each direction.
void exitDifferentResolutionError()
{
  cerr << red << "\nERROR: Image must have the same resolution in each direction.\n" << reset;
  exit(42);
}

// exit if the image has not the same resolution in each direction.
void exitImagesWithDifferentResolutionsError()
{
  cerr << red << "\nERROR: All images must have the same resolution.\n" << reset;
  exit(42);
}

int main(int argc, char** argv)
{
  try {

    if (argc <= 1) {
      cerr <<red<< "usage: " << argv [0] << " data" << reset << endl;
      return 23;
    }

    char c;

    char *suffix = strrchr ( argv[1], '.' );
    if ( !strcmp ( suffix, ".png" ) ) c = 'P';
    else c = getFirstHeaderChar( argv[1] );


    // ------------------------- --------------------------------------------
    // ------------------------- 3D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( c == 'Q' )
    {
      // Speicher für die Daten
      qc::ScalarArray<double, qc::QC_3D> ** data = new qc::ScalarArray<double, qc::QC_3D>* [argc-1];

      // erste Datei um das mesh anzulegen
      data[0]->quietMode = true;
      cerr<<blue<<"\nLoading file '"<<green<<argv[1]<<blue<<"'...\nDimension: "<<reset;
      data[0] = new qc::ScalarArray<double, qc::QC_3D>( argv[1] );
      cerr<<"Success!\n"<<reset;

      // Zeilenzahl merken (direkt pruefen, ob Daten wuerfelfoermig sind)
      int Nx = data[0]->getNumX();
      int Ny = data[0]->getNumY();
      int Nz = data[0]->getNumZ();
      if (Nx != Ny || Nx != Nz || Ny != Nz) exitDifferentResolutionError();

      GENMESH3D* mesh = quocmesh_convert_to_gmesh3d (data[0], argv[1]);

      // alle Dateien durchgehen
      for (int count=2; count<argc; count++)
      {
          data[count-1]->quietMode = true;
          cerr<<blue<<"Loading file '"<<green<<argv[count]<<blue<<"'...\nDimension: "<<reset;
          c = getFirstHeaderChar( argv[count] );
          if ( c == 'P' ) exitMixedError();

          data[count-1] = new qc::ScalarArray<double, qc::QC_3D>( argv[count] );

          // test the resolution of the image
          int Nx2 = data[count-1]->getNumX();
          int Ny2 = data[count-1]->getNumY();
          int Nz2 = data[count-1]->getNumZ();
          if (Nx2 != Ny2 || Nx2 != Nz2 || Ny2 != Nz2) exitDifferentResolutionError();
          if (Nx != Nx2) exitImagesWithDifferentResolutionsError();

          cerr<<"Success!\n"<<reset;
          addScalarData(mesh, data[count-1], argv[count]);
      }

      cerr<<blue<<"\nStarting GRAPE...\n"<<reset;
      initStartGrape(mesh, argv[1]);
    }



    // ------------------------- --------------------------------------------
    // ------------------------- 2D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( c == 'P' )
    {
      // Speicher für die Daten
      qc::ScalarArray<double, qc::QC_2D> ** data = new qc::ScalarArray<double, qc::QC_2D>* [argc-1];

      // erste Datei um das mesh anzulegen
      data[0]->quietMode = true;
      cerr<<blue<<"\nLoading file '"<<green<<argv[1]<<endl<<reset;
      data[0] = new qc::ScalarArray<double, qc::QC_2D>( argv[1] );
      cerr<<"Success!\n"<<reset;

      int Nx = data[0]->getNumX();
      int Ny = data[0]->getNumY();
      if (Nx != Ny) exitDifferentResolutionError();

      GENMESH2D* mesh = quocmesh_convert_to_gmesh2d (data[0], argv[1]);

      // alle Dateien durchgehen
      for (int count=2; count<argc; count++)
      {
          data[count-1]->quietMode = true;
          cerr<<blue<<"Loading file '"<<green<<argv[count]<<endl<<reset;
          c = getFirstHeaderChar( argv[count] );
          if ( c == 'Q' ) exitMixedError();

          data[count-1] = new qc::ScalarArray<double, qc::QC_2D>( argv[count] );

          // test the resolution of the image
          int Nx2 = data[count-1]->getNumX();
          int Ny2 = data[count-1]->getNumY();
          if (Nx2 != Ny2) exitDifferentResolutionError();
          if (Nx != Nx2) exitImagesWithDifferentResolutionsError();

          cerr<<"Success!\n"<<reset;
          addScalarData(mesh, data[count-1], argv[count]);
      }

      cerr<<blue<<"\nStarting GRAPE...\n"<<reset;
      initStartGrape(mesh, argv[1]);
    }


    cerr << "\nThanks for using GRAPEFRUIT v.1.41421356, NON-REGISTERED version!!\n"<<reset;

  }// try
  catch(aol::Exception e) {
    e.dump();
    return 42;
  }
}

#else

int main ( int, char** ) {
  cerr << "Without grape external, this program is useless" << endl;
  return ( 0 ) ;
}

#endif
