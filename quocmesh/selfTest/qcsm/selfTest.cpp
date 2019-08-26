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

#include <graphToTriangMesh.h>
#include <levelsetToTriangMesh.h>
#include <triangMeshToDistFunc.h>

#include <scalarArray.h>


typedef double RealType;

int  depth = 1;
int  width = (1<<depth)+1;
RealType r = .45; // radius
RealType h =  1./( width - 1 );

int main( int, char ** ){
  try{

    bool success = true;

    {
      const RealType RR = r*r;
      qc::ScalarArray<RealType, qc::QC_3D> levelset( width, width, width );
      for ( int z = 0; z < width; z++ ) {
        const RealType ZZ = aol::Sqr( z*h-.5 );
        for ( int y = 0; y < width; y++ ) {
          const RealType YY = aol::Sqr( y*h-.5 );
          for ( int x = 0; x < width; x++ ) {
            const RealType XX = aol::Sqr( x*h-.5 );
            levelset.set( x,y,z,  XX + YY + ZZ - RR );
          }
        }
      }

      aol::TriangMesh<RealType> smesh;
      qcsm::LevelsetToTriangMesh<RealType> convert;
      convert( levelset, smesh );

      aol::TriangMesh<RealType> referenceTriangulation;
      referenceTriangulation.loadFromUDPLY ( "referenceTriangulation.ply" );

      success &= referenceTriangulation.isApproxEqual ( smesh );
    }

    {
      qc::ScalarArray<RealType, qc::QC_2D> img ( 5, 5 );
      for ( int i = 0; i < img.size(); i++ ) img[i] = static_cast<RealType>(i);

      aol::TriangMesh<RealType> smesh;
      qcsm::GraphToTriangMesh<RealType> converter;
      converter.apply ( img, smesh );

      aol::TriangMesh<RealType> referenceGraph;
      referenceGraph.loadFromUDPLY ( "referenceGraph.ply.bz2" );

      success &= referenceGraph.isApproxEqual ( smesh );
    }

    if ( success ) {
      aol::printSelfTestSuccessMessage ( "--                       QCSM Self Test Successful                            --" );
      return ( EXIT_SUCCESS );
    }
    else {
      aol::printSelfTestFailureMessage ( "!!                       QCSM Self Test FAILED                                !!" );
      return ( EXIT_FAILURE );
    }

  }
  catch( aol::Exception e ){
    e.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}





