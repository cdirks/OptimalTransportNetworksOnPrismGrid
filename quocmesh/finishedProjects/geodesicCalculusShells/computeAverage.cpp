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

#include <aol.h>
#include <parameterParser.h>

#include "averagingOps.h"

//!------------------------------------------------------------------------------------------------------------------------------------------------------
//!------------------------------------------------------------------------------------------------------------------------------------------------------
//!------------------------------------------------------------------------------------------------------------------------------------------------------

typedef double RType;
typedef om::TriMesh<RType> TriMeshType;
typedef aol::TriangMeshConfigurator<RType, TriMeshType, aol::CenterQuadrature <RType> > TriMeshCenterQuadConf;

//! main
int main(int argc, char **argv)
{

  try {    

    if ( argc != 3 ){  
      cerr << "Wrong number of input files." << endl;
      cerr << "Usage: ./computeAverage <averageMode> <parameterParser>" << endl;
      return EXIT_FAILURE;
    }
    
    //
    int averageMode = atoi( argv[1] );
    
    // read in parameters from specified parameter file    
    char parsername[1024];
    sprintf( parsername, "%s",  argv[2] );      
    cerr << endl << "Reading parameters from '" << parsername << "'." << endl;
    aol::ParameterParser pparser( parsername );     
    
      
    cerr<<"\n===========================================================================================================\n";
    aol::StopWatch time;
    time.start();    
 
    if( averageMode == 0 )
      ElasticAverageOp<TriMeshCenterQuadConf, MembraneDeformation<TriMeshCenterQuadConf>, SimpleBendingDeformation<TriMeshCenterQuadConf> >( pparser ).execute();
    if( averageMode == 1 )
      GeodesicAverageOp<TriMeshCenterQuadConf, MembraneDeformation<TriMeshCenterQuadConf>, SimpleBendingDeformation<TriMeshCenterQuadConf> >( pparser ).execute();
    if( averageMode > 1 )
      cerr << "ERROR. Unknown average mode." << endl;
    
    time.stop();
    cerr<<"\n===========================================================================================================\n";
    cerr<<"Overall wall clock time = "<<time.elapsedWallClockTime ()<<" seconds.\n";
    cerr<<"===========================================================================================================\n";

  } catch ( aol::Exception &el ) {
    el.dump();
  }

  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
