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

/**
 * \file
 * \brief Extracts and saves the segments implicitly encoded in a number of levelset functions
 *        (along the lines of the Vese-Chan model).
 *
 * Usage: extractSegmentsFromProductSegmentation <threshold> <hardThreshold> <inputFile1> ... <inputFileN>
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>
#include <colorWheel.h>
#include <ChanVese.h>
#include <configurators.h>
#include <generator.h>

typedef double RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType, qc::QC_2D, 3> > ConfType2d;

int main ( int argc, char **argv ) {

  try {
    if ( argc < 4 ) {
      cerr << "USAGE: " << argv[0] << "  <threshold> <hardThreshold> <inputFile1> ... <inputFileN>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const int numberOfLevelsetFunctions = argc - 3;
    const RType thresholdValue = atof ( argv[1] );
    const bool hardThreshold = ( atoi ( argv[2] ) != 0 );
    qc::ArrayHeader header;
    {
      // Make sure that in is destroyed by this additional scope, so that the file argv[3]
      // can be opened again later.
      aol::Bzipifstream in ( argv[3] );
      qc::ReadArrayHeader ( in, header );
    }
    aol::MultiVector<RType> levelsetFunctions ( numberOfLevelsetFunctions, header.numX * header.numY * header.numZ );
    aol::MultiVector<RType> segments ( 1 << numberOfLevelsetFunctions, header.numX * header.numY * header.numZ );

    // Read the levelsetfunctions depending on their dimension.
    if ( header.magic[0] == 'P' ) {
      for ( int i = 0; i < numberOfLevelsetFunctions; ++i ) {
        qc::ScalarArray<RType, qc::QC_2D>  array ( levelsetFunctions[i], header.numX, header.numY, aol::FLAT_COPY );
        array.load ( argv[3+i] );
      }
    } else if ( header.magic[0] == 'Q' ) {
      for ( int i = 0; i < numberOfLevelsetFunctions; ++i ) {
        qc::ScalarArray<RType, qc::QC_3D>  array ( levelsetFunctions[i], header.numX, header.numY, header.numZ, aol::FLAT_COPY );
        array.load ( argv[3+i] );
      }
    }

    if ( hardThreshold ) {
      for ( int i = 0; i < levelsetFunctions[0].size(); ++i ) {
        segments[aol::PowerSetIterator::getPositionNumberFromLevelsetFunctions ( levelsetFunctions, i, thresholdValue ) ][i] = 1;
      }
    } else {
      aol::Vector<RType> temp ( numberOfLevelsetFunctions );
      for ( int i = 0; i < levelsetFunctions[0].size(); ++i ) {
        aol::PowerSetIterator iterator ( levelsetFunctions.numComponents() );
        do {
          for ( int j = 0; j < numberOfLevelsetFunctions; ++j ) {
            const RType sign = 2 * ( iterator.getComponent ( j ) - 0.5f );
            temp[j] = ( levelsetFunctions[j][i] - thresholdValue ) * sign;
          }

          segments[iterator.getCurrentPosition() ][i] = temp.getMinValue() + thresholdValue;
          iterator.increment();
        } while ( !iterator.end() );
      }
    }

    if ( header.magic[0] == 'P' ) {
      qc::RectangularGrid<qc::QC_2D> grid ( aol::Vec3<int> ( header.numX, header.numY, 1 ) );
      qc::MultiArray<RType, 2, 3> vImageArray ( grid );
      qc::DataGenerator<ConfType2d> generator ( grid );
      aol::Vector<RType> backgroundImage ( grid );
      backgroundImage.setAll ( 255. );
      generator.generateColoredSegmentation ( numberOfLevelsetFunctions, levelsetFunctions, backgroundImage, vImageArray, thresholdValue );
      vImageArray.savePNG ( "out.png" );

      for ( int i = 0; i < segments.numComponents(); ++i ) {
        qc::ScalarArray<RType, qc::QC_2D>  array ( segments[i], header.numX, header.numY, aol::FLAT_COPY );
        string outfileName = aol::strprintf ( "output%d.png", i );
        array.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
        array.savePNG ( outfileName.c_str() );
      }
    } else if ( header.magic[0] == 'Q' ) {
      for ( int i = 0; i < segments.numComponents(); ++i ) {
        qc::ScalarArray<RType, qc::QC_3D>  array ( segments[i], header.numX, header.numY, header.numZ, aol::FLAT_COPY );
        string outfileName = aol::strprintf ( "output%d.bz2", i );
        array.save ( outfileName.c_str(), qc::PGM_DOUBLE_BINARY );
      }
    }

  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
