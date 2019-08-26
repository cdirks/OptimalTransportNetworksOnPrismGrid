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

// 1st program to load muCT datasets obtained by UW:
// cut out region of interest and save as quoc data set

// #define SAVE_PGMS 1

#include <parameterParser.h>
#include <scalarArray.h>
#include <tpCFEGrid.h>

#include <quocMultigrid.h>

using namespace std;

typedef double RealType;

int main ( int argc, char **argv ) {
  try {

    aol::ParameterParser params ( argc == 2 ? argv[1] : "par/prepareMuctScan.par" );

    char scan_filenamemask[1024], roiFilename[1024];
    params.getString ( "scan_filenamemask", scan_filenamemask );
    params.getString ( "roiFilename", roiFilename );

    const int
      scan_size_x      = params.getInt ( "scan_size_x" ),          // number of pixels in the muct-scan in each slice
      scan_size_y      = params.getInt ( "scan_size_y" ),
      scan_first_slice = params.getInt ( "scan_first_slice" ),     // first slice to be read (included)
      scan_last_slice  = params.getInt ( "scan_last_slice" );      // last  slice to be read (included)

    qc::ScalarArray< RealType, qc::QC_3D > loadImage ( scan_size_x, scan_size_y, scan_first_slice + scan_last_slice + 1 );
    loadImage.setQuietMode ( true );
    cerr << "loading raw data ...";
    loadImage.loadSlices ( scan_filenamemask, qc::QC_Z, scan_first_slice, scan_last_slice );
    cerr << " done." << endl;

#ifdef SAVE_PGMS
    loadImage.saveSlices ( "out/loadImage_%03d.pgm", qc::QC_Z, qc::PGM_UNSIGNED_CHAR_BINARY, NULL, aol::CLIP_THEN_SCALE, 0, 255 );
#endif

    const int
      origin_x = params.getInt ( "origin_x" ),                     // origin of the area of interest
      origin_y = params.getInt ( "origin_y" ),
      origin_z = scan_first_slice,
      length_x = params.getInt ( "length_x" ),                     // size of the area of interest (after rotation, performed e. g. in gimp)
      length_y = params.getInt ( "length_y" ),
      length_z = scan_last_slice - scan_first_slice;

    cerr << " extracting region of interest from dataset ...";

    qc::ScalarArray< RealType, qc::QC_3D > ROI ( length_x + 1, length_y + 1, length_z + 1 );
    for ( qc::RectangularIterator< qc::QC_3D > rit ( ROI ); rit.notAtEnd(); ++rit ) {
      ROI.set ( *rit, loadImage.get ( (*rit)[0] + origin_x, (*rit)[1] + origin_y, (*rit)[2] + origin_z ) );
    }

#ifdef SAVE_PGMS
    ROI.saveSlices ( "out/ROI_%03d.pgm", qc::QC_Z, qc::PGM_UNSIGNED_CHAR_BINARY, NULL, aol::CLIP_THEN_SCALE, 0, 255 );
#endif

    cerr << "Minimal and maximal values are: " << ROI.getMinValue() << " " << ROI.getMaxValue() << endl;

#if 0

    // should call RiCaThreshold<RealType> instead
    RealType estimatedThreshold = ROI.sum() / ROI.size();
    {
      aol::Vector<int> histogram ( 256 );
      for ( int i = 0; i < ROI.size(); ++i ) {
        ++( histogram[ static_cast<int>( ROI[i]  ) ] );
      }

      ofstream histo ( "out/histogram.dat" );
      for ( int i = 0; i < histogram.size(); ++i )
        histo << i << " " << histogram[i] << endl;

      for ( int iter = 0; iter < 10; ++iter ) {
        RealType lowerNum = 0, lowerDen = 0, upperNum = 0, upperDen = 0;
        for ( int i = 0; i <= static_cast<int>( estimatedThreshold ); ++i ) { // this should probably be <=, even if written wrongly in the corresponding Paper (Tr79).
          lowerNum += i * histogram[i];
          lowerDen += 2 * histogram[i];
        }
        for ( int i = static_cast<int>( estimatedThreshold ) + 1; i < 256 ; ++i ) {
          upperNum += i * histogram[i];
          upperDen += 2 * histogram[i];
        }

        estimatedThreshold = lowerNum / lowerDen + upperNum / upperDen;

        cerr << estimatedThreshold << endl;
      }
    }

    cerr << "Estimated threshold for segmentation is " << estimatedThreshold << "; thresholding ... " << endl;

    ROI.addToAll ( -estimatedThreshold );
    ROI *= -1.0;
#endif

    cerr << "Save dataset ...";
    ROI.save ( roiFilename, qc::SaveTypeTrait<RealType>::BinarySaveType );
    cerr << " done." << endl;

  } catch ( aol::Exception e ) {
    e.dump();
  }

  return ( EXIT_SUCCESS );
}

