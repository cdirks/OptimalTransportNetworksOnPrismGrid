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
 * \brief Converts multiple 2D ScalarArrays to PNG, rescaling them uniformly to the PNG intensity range.
 *
 * Usage: convertMultipleQuoc2DToSameScalePNG InputFile1 ... InputFileN
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>
#include <rgbColorMap.h>
#include <multiArray.h>

typedef double RType;

int main ( int argc, char **argv ) {

  try {
    if ( argc < 2 ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile1> ... <InputFileN>" << endl;
      return EXIT_FAILURE;
    }

    const int numberOfQuocArrays = argc - 1;
    std::vector<std::string> inputFileNames ( numberOfQuocArrays );
    for ( int i = 0; i < numberOfQuocArrays; ++i )
      inputFileNames[i] = argv[1+i];

    qc::GridSize<qc::QC_2D> gridSize ( qc::getSizeFromArrayFile ( inputFileNames[0] ) );;

    aol::MultiVector<RType> quocVectors ( numberOfQuocArrays, gridSize.getNumX() * gridSize.getNumY() );

    for ( int i = 0; i < numberOfQuocArrays; ++i ) {
      qc::ScalarArray<RType, qc::QC_2D>  tempArray ( quocVectors[i], gridSize.getNumX(), gridSize.getNumY(), aol::FLAT_COPY );
      tempArray.load ( inputFileNames[i].c_str() );
    }

    const RType minValue = quocVectors.getMinValue();
    const RType maxValue = quocVectors.getMaxValue();
    const aol::RGBColorMap<RType> hsvMap ( 0., 1.,  aol::HSV_BLUE_TO_RED );
    qc::MultiArray<RType, 2, 3> tempMArray ( gridSize );
    tempMArray.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );

    for ( int i = 0; i < numberOfQuocArrays; ++i ) {
      qc::ScalarArray<RType, qc::QC_2D>  tempArray ( quocVectors[i], gridSize.getNumX(), gridSize.getNumY(), aol::FLAT_COPY );
      tempArray.setOverflowHandling ( aol::CLIP_THEN_SCALE, minValue, maxValue );
      tempArray.savePNG ( aol::strprintf ( "%s_%f_%f.png", aol::getBaseFileName ( inputFileNames[i] ).c_str(), minValue, maxValue ).c_str() );

      for ( int k = 0; k < tempArray.size(); ++k ) {
        aol::Vec3< RType > color;
        hsvMap.scalarToColor ( ( tempArray[k] - minValue ) / ( maxValue - minValue ), color );
        for ( int j = 0; j < 3; ++j )
          tempMArray[j][k] = color[j];
      }
      tempMArray.savePNG ( aol::strprintf ( "%s_col_%f_%f.png", aol::getBaseFileName ( inputFileNames[i] ).c_str(), minValue, maxValue ).c_str() );
    }
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
