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

// test program to check computation of normals for aligned planar and spherical interfaces

#include<aol.h>
#include<tpCFEGrid.h>

#include "osTestfunction.h"

const int DEPTH = 5;

template< typename RealType >
bool test_normal_computation_OK ( tpcfe::InterfaceTestFunction<RealType> &testFct, RealType tolerance = 0 ) {
  typedef tpcfe::CFEGrid<RealType, tpcfe::CFE_TPOS> GridType;
  GridType cfegrid( DEPTH );
  qc::ScalarArray<RealType, qc::QC_3D> levelset( cfegrid ), dummy ( cfegrid );
  qc::AArray<RealType, qc::QC_3D> diffCoeff( cfegrid ); // dummy: affine function; not used here.

  testFct.setWidth( cfegrid.getWidth() );
  testFct.setCoeffs( 1.0, 1.0 );

  testFct.createInterfaceFunctionCoeffs ( levelset, dummy, diffCoeff );

  cfegrid.addStructureFrom ( levelset );

  cfegrid.detectAndInitVirtualNodes ( diffCoeff );

  const qc::GridSize<qc::QC_3D> gridSize ( cfegrid );
  for ( typename GridType::FullElementIterator it ( cfegrid ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> el ( *it, gridSize, cfegrid.getElType ( *it ) );

    el.computeAssembleData ( cfegrid );

    for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {
      for ( int i = 0; i < 4; i++ ) {
        int gIdx0, gIdx1;
        gIdx0 = el.globalIndex( ( *tit ) ( i,0 ) );
        if ( ( *tit ) ( i, 1 ) == 11 ) {
          gIdx1 = -1;
        } else {
          gIdx1 = el.globalIndex( ( *tit ) ( i,1 ) );
        }
        if ( gIdx1 != -1 ) {  // Node is interpolated
          const tpcfe::CFEVirtualNode<RealType, tpcfe::CFE_TPOS, RealType> &vn = cfegrid.getVirtualNodeRef ( gIdx0, gIdx1 );
          aol::Vec3<RealType> position, normal;
          tit->computeGlobalCoordinate ( position, el, i );
          testFct.gradientInterfaceValueAt_global( position, normal );
          normal.normalize();

          if ( (normal - vn._normal).norm() > tolerance ) {
            cerr << "analytic normal = " << normal << ";  approximate normal = " << vn._normal << "; error norm = " << (normal - vn._normal).norm() << endl;
          }

        }
      }
    }
  }

  return ( true );
}


typedef double RealType;

int main ( int, char** ) {
  try {

    tpcfe::InterfaceTestFunction<RealType>* testFct = NULL;

    cerr << "Verifying computation of interface normals:" << endl;
    {
      cerr << "Testing aligned planar interface ... ";
      testFct = tpcfe::selectITFByType<RealType> ( 0 );
      cerr << ( test_normal_computation_OK ( *testFct, 1.0e-16 ) ? "OK" : "failed" ) << endl;
      delete ( testFct );
    }

    {
      cerr << "Testing rotated planar interface ... ";
      testFct = tpcfe::selectITFByType<RealType> ( 2 );
      cerr << ( test_normal_computation_OK ( *testFct, 1.0e-14 ) ? "OK" : "failed" ) << endl;
      delete ( testFct );
    }

    {
      cerr << "Testing non-centered cylindrical interface ... ";
      testFct = tpcfe::selectITFByType<RealType> ( 4 );
      cerr << ( test_normal_computation_OK ( *testFct, 1.0e-2 ) ? "OK" : "failed" ) << endl;
      delete ( testFct );
    }

    {
      cerr << "Testing centered cylindrical interface ... ";
      testFct = tpcfe::selectITFByType<RealType> ( 6 );
      cerr << ( test_normal_computation_OK ( *testFct, 5.0e-2 ) ? "OK" : "failed" ) << endl;
      delete ( testFct );
    }

    {
      cerr << "Testing non-centered spherical interface ... ";
      testFct = tpcfe::selectITFByType<RealType> ( 8 );
      cerr << ( test_normal_computation_OK ( *testFct, 5.0e-2 ) ? "OK" : "failed" ) << endl;
      delete ( testFct );
    }

    {
      cerr << "Testing centered spherical interface ... ";
      testFct = tpcfe::selectITFByType<RealType> ( 10 );
      cerr << ( test_normal_computation_OK ( *testFct, 7.5e-2 ) ? "OK" : "failed" ) << endl;
      delete ( testFct );
    }

  } catch ( aol::Exception &exc ) {

    exc.dump();
    return ( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );
}

