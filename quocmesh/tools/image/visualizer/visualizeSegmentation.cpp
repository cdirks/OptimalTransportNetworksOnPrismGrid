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
 * \brief Visualizes a segmentation in two segments by coloring the segmented image according to the segmentation.
 *
 * Usage: visualizeSegmentation InputSegmentation SegmentedImage
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>
#include <rgbColorMap.h>
#include <multiArray.h>

int main ( int argc, char **argv ) {

  try {
    string inputSegmentationFileName;
    string segmentedImageFileName;

    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <InputSegmentation> <SegmentedImage>" << endl;
      return EXIT_FAILURE;
    }
    if ( argc == 3 ) {
      inputSegmentationFileName = argv[1];
      segmentedImageFileName = argv[2];
    }

    cerr << inputSegmentationFileName.c_str() << endl;
    cerr << segmentedImageFileName.c_str() << endl;
    qc::ScalarArray<double, qc::QC_2D> inputSegmentation ( inputSegmentationFileName.c_str() );
    inputSegmentation /= inputSegmentation.getMaxValue();
    qc::ScalarArray<double, qc::QC_2D> segmentedImage ( segmentedImageFileName.c_str() );
    segmentedImage /= segmentedImage.getMaxValue();

    aol::RGBColorMap<double> hsvMap ( 0., 1.,  aol::HSV_BLUE_TO_RED );

    qc::MultiArray<double, 2, 3> outputImage ( inputSegmentation.getNumX(), inputSegmentation.getNumY() );

    for ( int i = 0; i < segmentedImage.size(); ++i ) {
      aol::Vec3< double > color;
      hsvMap.scalarToColor ( inputSegmentation[i], color );
      color *= segmentedImage[i];
      for ( int j = 0; j < 3; ++j )
        outputImage[j][i] = color[j];
    }

    outputImage.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
    outputImage.savePNG ( "output.png" );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
