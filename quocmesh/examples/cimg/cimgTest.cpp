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
 * \brief Very rudimentary example of how the CImg library can be used.
 *
 * Converts a qc::ScalarArray to the corresponding CImg data structure, display the image
 * and reports the pixel intensity when the user clicks on it in the display window.
 *
 * \author Berkels
 */

#include <cimgIncludes.h>
#include <scalarArray.h>

int main() {
  qc::ScalarArray<unsigned char, qc::QC_2D> imageArray ( "../testdata/image_129.pgm.bz2" );
  cimg_library::CImg<unsigned char> image(imageArray.getData(), imageArray.getNumX(), imageArray.getNumY() );
  cimg_library::CImgDisplay display ( image, "Click a point" );
  while ( !display.is_closed() ) {
    display.wait();
    if ( display.button()  ) {
      const int x = display.mouse_x();
      const int y = display.mouse_y();
      cerr << "Value at pixel (" << x << "," << y << ") is " << static_cast<int> ( image( x, y ) ) << endl;
    }
  }
  return 0;
}
