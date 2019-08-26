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

// determine the error between a continuous function satisfying a given kink condition and the CFE representation of that function.

// the tpCFE approximation is compared to an approximation with pcw affine FE on the same cubic grid and pcw trilinear FE.

#include <tpCFEGrid.h>
#include <parameterParser.h>

#include "osTestfunction.h"

#include <levelsetToTriangMesh.h>

static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;

typedef double RealType;

int main ( int argc, char** argv) {
  try {

    aol::ParameterParser params ( argc < 2 ? "par/os_test_approximation.par" : argv[1] );
    cerr << "Reading from " << ( argc < 2 ? "par/os_test_approximation.par" : argv[1] ) << endl;

    aol::Vector<int> types;
    params.getIntVec( "types", types );

    for ( int i = 0; i < types.size(); ++i ) {
      const int
        TYPE = types[i],
        depth = params.getInt ( "draw_depth" );

      const RealType
        DIFF_COEFF_PLUS  = 1.0,
        DIFF_COEFF_MINUS = params.getDouble ( "dc_minus_draw" );


      tpcfe::CFEGrid < RealType, AT >  grid ( depth );

      grid.setVerboseMode(true);
      grid.setAdjustLevelset( 0.0 );

      qc::ScalarArray<RealType, qc::QC_3D>
        interfc ( grid ),           // levelset function defining the interface
        function ( grid );            // values of the "affine" function

      qc::AArray<RealType, qc::QC_3D>
        coeff ( grid );               // diffusion coefficient at grid points

      tpcfe::InterfaceTestFunction<RealType> *testFct = tpcfe::selectITFByType<RealType> ( TYPE );

      testFct->setCoeffs ( DIFF_COEFF_PLUS, DIFF_COEFF_MINUS );
      testFct->setWidth ( grid.getWidth() );
      testFct->createInterfaceFunctionCoeffs ( interfc, function, coeff );

      {
        const int w = function.getNumX() - 1 ;
        cerr << function.get ( 0, 0, 0 ) << endl;
        cerr << function.get ( 0, 0, w ) << endl;
        cerr << function.get ( 0, w, 0 ) << endl;
        cerr << function.get ( 0, w, w ) << endl;
        cerr << function.get ( w, 0, 0 ) << endl;
        cerr << function.get ( w, 0, w ) << endl;
        cerr << function.get ( w, w, 0 ) << endl;
        cerr << function.get ( w, w, w ) << endl;
        cerr << function.get ( w/2, w/2, w/2 ) << endl;
      }

      aol::TriangMesh<RealType> tmesh;
      qcsm::LevelsetToTriangMesh<RealType> converter;
      converter.apply ( interfc, tmesh );

      tmesh.createVertexData();

      RealType minVal = aol::NumberTrait<RealType>::Inf, maxVal = - aol::NumberTrait<RealType>::Inf;

      for ( int pt = 0; pt < tmesh.getNumVertices(); ++pt ) {
        const aol::Vec3<RealType> coord = tmesh.getVertex ( pt ); // image coordinates
        const RealType value = testFct->valueAt ( coord );
        minVal = min ( minVal, value );
        maxVal = max ( maxVal, value );
        tmesh.getVertexData()[pt] = value;
      }

      cerr << "Minimal value was " << aol::detailedFormat ( minVal ) << ", maximal value was " << aol::detailedFormat ( maxVal ) << endl;

      char filename[1024];
      sprintf ( filename, "out/fctOnLevelset_%02d.ply", TYPE );
      tmesh.saveAsUDPLY ( filename );

    } // end type loop

  } catch ( aol::Exception &exc ) {

    exc.dump();

    return ( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );
}
