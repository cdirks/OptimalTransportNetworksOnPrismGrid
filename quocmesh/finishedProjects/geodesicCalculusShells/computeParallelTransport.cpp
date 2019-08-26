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
//! Program computes parallel transport of a variational shape along a given path.
//! Discrete parallel transport is realized via Schilds ladder and geodesic parallelogram construction.
//! Diagonal geodesics in each geodesic parallelogram are given by N-point geodesics.
//! Note: Control via parameter file possible, e.g. use your own modification of parComputeParallelTransport.txt!

//! main
int main(int argc, char **argv)
{
  try {        
    
    // check whether number of arguments is correct
    if ( argc != 2 ) { 
      cerr << "Wrong number of input files." << endl;
      cerr << "Usage: ./computeParallelTransport <parameterFile>" << endl;
      return EXIT_FAILURE;
    }
    
    // read in file names of parameter files
    char parsername[1024];
    sprintf( parsername, "%s",  argv[1] );
    cerr << endl << "Reading parameters from '" << parsername << "'." << endl;
    aol::ParameterParser pparser( parsername );  
    
    // optimization    
    cerr<<"\n===========================================================================================================\n";
    aol::StopWatch time;
    time.start();        
    
    ParallelTransportOp<TriMeshCenterQuadConf, MembraneDeformationType, BendingDeformationType>( pparser ).execute( ); 
    
    time.stop();
    cerr<<"\n===========================================================================================================\n";
    cerr<<"Overall wall clock time = "<<time.elapsedWallClockTime ()<<" seconds.\n";


  } catch ( aol::Exception &el ) {
    el.dump();
  } 

  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}