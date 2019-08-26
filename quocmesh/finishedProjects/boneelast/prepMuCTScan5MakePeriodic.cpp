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

#include <parameterParser.h>
#include <scalarArray.h>
#include <tpCFEGrid.h>
#include <multilevelArray.h>
#include <iterators.h>

using namespace std;


typedef double RealType;

int main ( int argc, char **argv ) {
  if( argc != 2 ) {
    cerr << "usage: prepMuCTScan5MakePeriodic PARAMETERFILE" << endl;
    return( EXIT_FAILURE );
  }

  try {

    aol::ParameterParser params( argv[1] );

    char inFilename [1024], outFilename [1024];
    params.getString ( "cell_filename", inFilename );
    params.getString ( "periodic_filename", outFilename );

    const int
      inDepth = params.getInt ( "in_depth" ),
      outDepth = params.getInt ( "out_depth" ),
      periodicDepth = params.getInt ( "periodic_depth" );

    qc::GridDefinition inGrid ( inDepth, qc::QC_3D ), outGrid ( outDepth, qc::QC_3D ), periodicGrid ( periodicDepth, qc::QC_3D );
    qc::ScalarArray<RealType, qc::QC_3D> inImage ( inGrid ), outImage ( outGrid );

    inImage.load ( inFilename );

    qc::MultilevelArray< double, qc::ScalarArray<double, qc::QC_3D> > inMLArray ( inGrid );
    inMLArray[inDepth] = inImage;
    inMLArray.levRestrict ( periodicDepth , inDepth );

    const int
      npx = params.getInt ( "num_periodic_x" ),
      npy = params.getInt ( "num_periodic_y" ),
      npz = params.getInt ( "num_periodic_z" ),
      mpx = ( outImage.getNumX() ) / npx,
      mpy = ( outImage.getNumY() ) / npy,
      mpz = ( outImage.getNumZ() ) / npz;


    const int
      xOff  = params.getInt ( "x_off_p"  ),
      yOff  = params.getInt ( "y_off_p"  ),
      zOff  = params.getInt ( "z_off_p"  ),
      xRang = params.getInt ( "x_rang_p" ),
      yRang = params.getInt ( "y_rang_p" ),
      zRang = params.getInt ( "z_rang_p" );

    const aol::Vec3<short int> off ( xOff, yOff, zOff ), rang ( xRang, yRang, zRang );

    qc::ScalarArray<double, qc::QC_3D> smallImage ( aol::Vec3<int>( off + rang ) );

    for ( qc::RectangularIterator<qc::QC_3D> bit ( smallImage ); !( bit.atEnd() ); ++bit ) {
      smallImage.set ( *bit, inMLArray[periodicDepth].get ( *bit + off ) );
    }

    qc::ScalarArray<double, qc::QC_3D> periodicImage ( npx * mpx, npy * mpy, npz * mpz );

    for ( qc::RectangularIterator<qc::QC_3D> bit ( aol::Vec3<short int>( 0, 0, 0 ), aol::Vec3<short int>( npx, npy, npz ) ); !( bit.atEnd() ); ++bit ) {
      const aol::Vec3<short int> pad( (*bit)[0] * mpx, (*bit)[1] * mpy, (*bit)[2] * mpz );

      for ( qc::RectangularIterator<qc::QC_3D> ibit ( aol::Vec3<short int>( 0, 0, 0 ), aol::Vec3<short int> ( mpx, mpy, mpz ) ); !( ibit.atEnd() ); ++ibit ) {
        aol::Vec3<double> posi ( ( 1.0 * (*ibit)[0] ) / mpx, ( 1.0 * (*ibit)[1] ) / mpy, ( 1.0 * (*ibit)[2] ) / mpz );
        periodicImage.set ( *ibit + pad, smallImage.interpolate_on01 ( posi ) );
      }

    }

    for ( qc::RectangularIterator<qc::QC_3D> bit ( periodicImage ); !( bit.atEnd() ); ++bit ) {
      outImage.set ( *bit, periodicImage.get ( *bit ) );
    }

    outImage.save ( outFilename, qc::SaveTypeTrait<RealType>::BinarySaveType );

  } catch ( aol::Exception e ) {
    e.dump();
  }

  return ( 0 );
}

