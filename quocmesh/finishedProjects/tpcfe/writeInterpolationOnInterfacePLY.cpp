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

#include "parameterParser.h"
#include <tpCFEUtils.h>

template< typename GridType, typename ArrayType >
void writeTemperatureSlicePLYs ( const GridType &grid, const ArrayType &temperature, const char* fnmask, const short sliceIndex ) {
  typedef typename GridType::RealType RealType;
  for ( short diri = 0; diri < 1; ++diri ) {
    const qc::Comp dir = static_cast<qc::Comp>( diri );

    tpcfe::CFEInterfaceTriangulationWithScalarDataGenerator< GridType > itg ( grid, temperature );
    itg.determineSliceTriangulation ( dir, sliceIndex );

    char filename[1024];
    sprintf ( filename, fnmask, dir, sliceIndex );
    itg.saveToPLYFile ( filename );

  }
}


typedef double RealType;

int main ( int argc, char **argv ) {
  if ( argc != 2 ) {
    cerr << "usage: writeInterpolationOnInterfacePLY PARAMETERFILE" << endl;
    return( EXIT_FAILURE );
  }

  try {

    aol::ParameterParser params ( argv[1] );

    const int depth = params.getInt ( "depth" );
    const short sliceIndex = static_cast<short> ( params.getInt ( "sliceIndex" ) );

    const RealType
      dcPlus = params.getDouble ( "dcPlus" ),
      dcMinus = params.getDouble ( "dcMinus" );

    qc::ScalarArray<RealType, qc::QC_2D>::quietMode = 1;

    cerr << "Creating domain ... ";

    typedef tpcfe::CFEGrid< RealType, tpcfe::CFE_TPOS > GridType;
    GridType grid ( depth );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid ), values ( grid );
    qc::AArray<RealType, qc::QC_3D> coeff ( grid );

    char geometryFilename[1024];
    params.getString ( "geometryFilename", geometryFilename );
    cerr << geometryFilename << endl;
    levelset.load ( geometryFilename );

    grid.addStructureFrom ( levelset );

    for ( int i = 0; i < levelset.size(); ++i )
      coeff[i] = ( levelset[i] < 0 ? dcMinus : dcPlus );

    grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-15, 5.0e-16 );

    cerr << " done." << endl;


    cerr << "Loading scalar values ... ";
    char valuesFilename[1024];
    params.getString ( "valuesFilename", valuesFilename );
    values.load ( valuesFilename );
    cerr << " done." << endl;

    cerr << "Generating and saving triangulation ... ";
    tpcfe::CFEInterfaceTriangulationWithScalarDataGenerator<GridType> itg( grid, values ); // need to load scalar data

    itg.determineInterfaceAndBoundaryTriangulation( );

    char PLYFilename[1024];
    params.getString ( "PLYFilename", PLYFilename );

    itg.saveToPLYFile( PLYFilename );
    cerr << " done." << endl;


    cerr << "Generating and saving slices ... ";
    char PLYFNMask[1024];
    params.getString ( "PLYFNMask", PLYFNMask );
    writeTemperatureSlicePLYs ( grid, values, PLYFNMask, sliceIndex );
    cerr << " done." << endl;

  } catch ( aol::Exception e ) {
    e.dump();
  }

  return ( 0 );
}
