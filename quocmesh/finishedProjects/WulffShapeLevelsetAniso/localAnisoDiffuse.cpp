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
#include <quoc.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <preconditioner.h>

#include <mcm.h>
#include <fastUniformGridMatrix.h>

#include <anisotropies.h>
#include <Willmore.h>


#include <qmException.h>
#include <aol.h>
#include <parameterParser.h>
#include "grapeInterface3d.h"

#include "moments.h"


// now the main program
int main( int argc, char **argv ) {

  // read the parameters with an parameter-file
  if ( argc != 2 ) {
    string s = "USAGE: ";
    s += argv[0];
    s += " <parameterfile>";
    throw aol::Exception( s.c_str(), __FILE__, __LINE__  );
  }

  try {
    aol::ParameterParser parser( argv[1] );

    // load image into scalar array
    char imgname[ 1024 ];
    parser.getString( "image", imgname );

    cerr<<aol::color::green<<"Lade Bild...\n"<<aol::color::black;
    moments<double> Daten(imgname, parser.getDouble( "delta"),
                                    parser.getDouble( "epsilon") );


    aol::StopWatch watch;
    watch.start();

    // now calc the masses and view them in GRAPE (change)
    if ( ! parser.getDouble( "load") )
    {
      // use no radius information
      cerr<<aol::color::red<<"Berechne Momente...\n"<<aol::color::black;
      if ( ! parser.getInt( "UseRadii") )
      {
        cerr<<aol::color::red<<"aber benutze keine Radius-Informationen.."<<aol::color::black;
        Daten.calcMass();
        Daten.calcCenterOfMass();
        parser.getString( "saveMomentsName", imgname );
        Daten.calcMoments(imgname);
      }
      else  // use radius information
      {
        cerr<<aol::color::red<<"unter Berücksichtigung von Radius-Informationen, "<<aol::color::black;

        cerr<<aol::color::green<<"\nlade Radien..."<<aol::color::black;
        parser.getString( "RadiiName", imgname );
        qc::ScalarArray<double, qc::QC_3D> Radii( imgname );
        Daten.setRadiusFactor( parser.getDouble("radiusFactor") );

        Daten.calcMassWithAdjustedRadii( Radii );
        Daten.calcCenterOfMassWithAdjustedRadii( Radii );
        parser.getString( "saveMomentsName", imgname );
        Daten.calcMomentsWithAdjustedRadii( imgname, Radii );
      }
    }
    else   // that means load moments
    {
      cerr<<aol::color::red<<"\nLade Momente ...\n"<<aol::color::black;
      parser.getString( "loadMomentsName", imgname );
      Daten.loadMoments( imgname );
    }

    cerr<<"fertig, berechne Richtungen...";

    Daten.calcMassCentroidAxis( parser.getDouble( "barrier") );
    parser.getString( "saveAxisName", imgname );
    Daten.saveMassCentroidAxis( imgname );
    Daten.showInGrape();

    watch.stop();
    cerr<<"Finished! Took "<<watch.elapsedCpuTime()<<" seconds!\nATTENTION: NOT REGISTERED YET!\n";

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}

