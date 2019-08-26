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
/* ***************************************************
 * Datei: read3DData.cpp
 * Funktion: lädt die Raw-Datei von Tolga ein
 *           und konvertiert sie in ein Quoc-Format
 * Autor: Oliver Nemitz
 * Datum: 07.09.2004
 * *************************************************** */

// #include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <scalarArray.h>
#include <quoc.h>
#include <aol.h>
#include <parameterParser.h>

#define pics    1           /* the two modes */
#define block3d 2


 int main( int argc, char **argv )
 {

   if ( argc != 2 ) {
      string s = "USAGE: ";
      s += argv[0];
      s += " <parameterfile>";
      throw aol::Exception( s.c_str(), __FILE__, __LINE__  );
  }

  try {
    aol::ParameterParser parser( argv[1] );


    cerr<<"Grösse von float: "<<sizeof(float)<<endl;

    // first read the dimension of the 3d-data
    int dimX = parser.getInt("dimX");
    int dimY = parser.getInt("dimY");
    int dimZ = parser.getInt("dimZ");


    // there are two modes: 1) extract the data into several pictures (mode = 1)
    // 2) extract the data into one 3D-Block (mode = 2)
    int mode = parser.getInt("mode");
    int threshold = parser.getInt("threshold");

    // read the necessary information about the block
    int BDim   = parser.getInt("BDim");    // width, height and depth of the block
    int StartX = parser.getInt("StartX");  // left upper edge of the block in coordinates
    int StartY = parser.getInt("StartY");  // of the original data
    int StartZ = parser.getInt("StartZ");

    // load image into scalar array
    char loadName[ 1024 ];
    parser.getString( "loadName", loadName );

    fstream Quelle;
    Quelle.open(loadName, ios::in);
    if (!Quelle)
    {
      cerr<<aol::color::red<<"Datei '"<<aol::color::blue<<loadName;
      cout<<aol::color::red<<"' kann nicht geöffnet werden!\n"<<aol::color::reset;
    }
    else
    {
      union {
        char b[4];
        // HACK!! before it was
//         short Wert; // now it is
        float Wert;
      } val;

      cerr<<"allocating memory...";

      // the array where the data should be saved in
      qc::ScalarArray<float, qc::QC_2D> *Bild  = new qc::ScalarArray<float, qc::QC_2D>(dimX,dimY);
      qc::ScalarArray<float, qc::QC_3D> *Block = new qc::ScalarArray<float, qc::QC_3D>(BDim+1,BDim+1,BDim+1);
      cerr<<"done, clearing structure... ";
      Block->setAll(0.);
      cerr<<"done!\n ";
      char saveName[256];

      // the size of the picture is known
      // fill several 2d-images with data and save them
      cerr <<aol::color::red<< "Starting the extraction ...";
      for (int z = 0; z < dimZ; z++)
      {
        if ( mode == pics  )   Bild->clear();
        if ( mode == block3d ) cerr<<".";

        for (int y = 0; y < dimY; y++)
          for (int x = 0; x < dimX; x++)
          {

            if (Quelle) {
              // HACK! This was only
//               Quelle.read(reinterpret_cast<char*>(&val.b[1]),1);
//               Quelle.read(reinterpret_cast<char*>(&val.b[0]),1);
//               Quelle.read(reinterpret_cast<char*>(&val.b[0]),1);
//               Quelle.read(reinterpret_cast<char*>(&val.b[1]),1);
//               Quelle.read(reinterpret_cast<char*>(&val.b[2]),1);
//               Quelle.read(reinterpret_cast<char*>(&val.b[3]),1);

              Quelle.read(reinterpret_cast<char*>(&val.Wert),4);

              if (val.Wert != 0) val.Wert = 1;//cerr<<val.Wert<<",";

              if (threshold) {
                if (val.Wert < 225) val.Wert = 0;
                else val.Wert = 1;
              }

              // save the whole picture
              if ( mode == pics ) Bild->set(x,y, val.Wert);
              // only save a cut-out of the picture
              if ( mode == block3d ) {
                if (x>=StartX && x<=StartX+BDim && y>=StartY && y<=StartY+BDim) {
                  switch(BDim) {
                    case 256:
                      Block->set(x,y,2*z+4, val.Wert);
                      Block->set(x,y,2*z+5, val.Wert);
                      break;
                    case 128:
                      Block->set(x-StartX, y-StartY, 2*z + 4, val.Wert);
                      Block->set(x-StartX, y-StartY, 2*z + 5, val.Wert);
                      break;
                    case 64:
                      Block->set(x-StartX, y-StartY, z, val.Wert);
                      break;
                    case 32:
                      if (z>=StartZ && z<StartZ+(BDim/2)) {
                          Block->set(x-StartX, y-StartY, 2*(z-StartZ), val.Wert);
                          Block->set(x-StartX, y-StartY, 2*(z-StartZ)+1, val.Wert);
                        }
                        // the last line fill out manual
                        if (z == StartZ + BDim) Block->set(x-StartX, y-StartY, BDim, val.Wert);
                      break;
                  }
                }
              }
            }
            else {
              cerr<<"\n\nERROR: END OF FILE REACHED BEFORE FINISHING !\n\n";
              return 42;
            }
          }

          // now save the image or the block respectively
          if ( mode == pics ) {
            if (z<10)
              sprintf(saveName, "/home/NOSAVE/achill/nemitz/Tolga/TolgaData_0%d.pgm", z);
            else
              sprintf(saveName, "/home/NOSAVE/achill/nemitz/Tolga/TolgaData_%d.pgm", z);
            cerr<<aol::color::reset;
            cerr<<"Saving picture number "<<aol::color::green<<z<<aol::color::reset<<"...\n";
            Bild->save(saveName);
          }
      }

      Quelle.close();


      // if BDim is 2 times higher than dimZ, interpolate the planes between the
      // data

      if ( mode == block3d ) {
        if (BDim > 2*dimZ) {
          cerr<<"done!\nStarting interpolation...";
          for (int z=1; z<BDim-1; z+=2) {
            cerr<<".";
            // now walk through the plane
            for (int y = 0; y < BDim; y++)
              for (int x = 0; x < BDim; x++) {
//                 if (z+1>128) cerr<<"Fehler!";
                double val = 0.5 * (Block->get(x,y,z-1) + Block->get(x,y,z+1));
                Block->set(x,y,z, val);
              }
          }
        }

        cerr<<aol::color::reset;
        // save everything
        cerr<<"done!\nSaving block '"<<aol::color::blue;
        // scale to [0,1]
        (*Block) -= Block->getMinValue();
        (*Block) /= Block->getMaxValue();
        char saveName[ 1024 ];
        parser.getString( "saveBlockName", saveName );
        cerr<<saveName<<aol::color::reset<<"'\n";
        Block->save( saveName );
        cerr<<"done! Thank you for extracting some data :-)\n";
      }

    }

  }   // try
    catch(aol::Exception e) {
    e.dump();
    return 42;
  }

 }

