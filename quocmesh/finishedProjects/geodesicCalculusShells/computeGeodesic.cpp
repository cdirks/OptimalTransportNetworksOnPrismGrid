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

#include <triangMeshConfigurators.h>
#include <triMesh.h>

#include "geodCalculusOps.h"

//!------------------------------------------------------------------------------------------------------------------------------------------------------
//!------------------------------------------------------------------------------------------------------------------------------------------------------
//!------------------------------------------------------------------------------------------------------------------------------------------------------//!====================================================================
typedef double RType;
typedef om::TriMesh<RType> TriMeshType;
typedef aol::TriangMeshConfigurator<RType, TriMeshType, aol::CenterQuadrature <RType> > TriMeshCenterQuadConf;


// classical membrane energy penalizing trace and determinant of Cauchy-Green tensor (i.e. edge length and triangle volume)
typedef MembraneDeformation<TriMeshCenterQuadConf> MembraneDeformationType;
// simple benidng energy refers to Discrete Shell energy introduced first by Grinspun et al 2003
typedef SimpleBendingDeformation<TriMeshCenterQuadConf> BendingDeformationType;
    

//! \author Heeren
//! Program computes discrete geodesic according to one of four optimization mode:
//!   (0) straight forward optimization (without any hierarchichal method)
//!   (1) multilevel in space (realized by progressive meshes)
//!   (2) red-black alternating optimization of odd and even shapes (final relaxation with whole sequence)
//!   (3) hierarchichal method in time realized by cascadic optimization 
//! Note: Control via parameter file possible, e.g. use your own modification of parComputeGeodesic.txt!

//! main
int main(int argc, char **argv)
{
  try {        
    
    // check whether number of arguments is correct
    if ( argc != 3 ) { 
      cerr << "Wrong number of input files." << endl;
      cerr << "Usage: ./computeGeodesic <optimizationMode> <parameterFile>" << endl;
      cerr << "optimization modes are straight forward (0), MultiLevel (1), RedBlack (2) or Cascadic (3)." << endl;
      return EXIT_FAILURE;
    }
    
    // read in file names of parameter files
    char parsername[1024];
    sprintf( parsername, "%s",  argv[2] );
    cerr << endl << "Reading parameters from '" << parsername << "'." << endl;
    aol::ParameterParser pparser( parsername );  

    // optimization mode:
    int mode = atoi ( argv[1] );
    
    // optimization    
    cerr<<"\n===========================================================================================================\n";
    aol::StopWatch time;
    time.start();        
    
    GeodesicOp<TriMeshCenterQuadConf, MembraneDeformationType, BendingDeformationType>( pparser ).execute( mode ); 
    
    time.stop();
    cerr<<"\n===========================================================================================================\n";
    cerr<<"Overall wall clock time = "<<time.elapsedWallClockTime ()<<" seconds.\n";


  } catch ( aol::Exception &el ) {
    el.dump();
  } 

  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}