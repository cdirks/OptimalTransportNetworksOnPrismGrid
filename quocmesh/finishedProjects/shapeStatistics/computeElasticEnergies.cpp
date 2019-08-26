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

#include <aol.h>
#include <configurators.h>
#include <smallVec.h>
#include <multiArray.h>
#include <deformations.h>
#include <hyperelastic.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * \brief This function converts the 2d-image "Image" (with ideally values in [0,1]) into a colored image
 * "ColorImage", which can be saved as png-file.
 */
template <typename RealType>
void convertToColor( const qc::ScalarArray<RealType, qc::QC_2D> &Image, qc::MultiArray<unsigned char,2,3> &ColorImage, const RealType Min = 0, const RealType Max = 1 ){
  const int numX = Image.getNumX();
  const int numY = Image.getNumY();

  // specify colors at certain values
  const int partitions = 4;
  RealType intervals[partitions+1] = {0, .25, .5, .75, 1};
  const RealType colorsRGB[partitions+1][3] = {{0,0,.5}, {0,.35,.15}, {.35,.35,.35}, {1,.35,.1}, {1,1,0}};

  // expand color intervals over the desired range
  for ( int i = 0; i <= partitions; i++ )
    intervals[i] = Min + intervals[i] * ( Max - Min );

  RealType rgb[3] = {0, 0, 0};

  // for each pixel do...
  for ( int y = 0; y < numY; y++ ) {
    for ( int x = 0; x < numX; x++ ) {
      // get pixel value
      RealType b = Image.get ( x, y );
      // assign the pixel the correct color
      if ( b <= Min ) {
        rgb[0] = colorsRGB[0][0];
        rgb[1] = colorsRGB[0][1];
        rgb[2] = colorsRGB[0][2];
      } else if ( b >= Max ) {
        rgb[0] = colorsRGB[partitions][0];
        rgb[1] = colorsRGB[partitions][1];
        rgb[2] = colorsRGB[partitions][2];
      } else {
        // find the correct color interval
        int partition = 1;
        while ( ( intervals[partition] < b ) && ( partition < partitions ) )
          partition++;
        rgb[0] = colorsRGB[partition-1][0] + ( colorsRGB[partition][0] - colorsRGB[partition-1][0] ) / ( intervals[partition] - intervals[partition-1] ) * ( b - intervals[partition-1] );
        rgb[1] = colorsRGB[partition-1][1] + ( colorsRGB[partition][1] - colorsRGB[partition-1][1] ) / ( intervals[partition] - intervals[partition-1] ) * ( b - intervals[partition-1] );
        rgb[2] = colorsRGB[partition-1][2] + ( colorsRGB[partition][2] - colorsRGB[partition-1][2] ) / ( intervals[partition] - intervals[partition-1] ) * ( b - intervals[partition-1] );
      }
      ColorImage[0][y*numX+x] = static_cast< unsigned char> ( rgb[0] * 255 );
      ColorImage[1][y*numX+x] = static_cast< unsigned char> ( rgb[1] * 255 );
      ColorImage[2][y*numX+x] = static_cast< unsigned char> ( rgb[2] * 255 );
    }
  }
}

/**
 * \brief This function produces a 2d-image "ChessImage" with a chess pattern of "NumberOfStripes" horizontal stripes.
 */
template <typename PixelValueType>
void generateChessPattern( const int NumberOfStripes, qc::ScalarArray<PixelValueType, qc::QC_2D> &ChessImage ){
  const int numX = ChessImage.getNumX();
  const int numY = ChessImage.getNumY();
  // compute stripe-width in pixels
  const int stripeWidth = static_cast<int>( numY/NumberOfStripes );
  // make image black
  ChessImage.setZero();
  // for each pixel do...
  for ( int y = 0; y < numY; y++ )
    for ( int x = 0; x < numX; x++ )
      // make white field
      if ( (y%(2*stripeWidth)<stripeWidth) == (x%(2*stripeWidth)<stripeWidth) )
        ChessImage.set( x, y, 255 );
}

//////////////////////////////////////////////////////////////////////////////////////
//! ---------------------------- main program --------------------------------------//
//////////////////////////////////////////////////////////////////////////////////////

// define the settings/configuration of the program, i.e. computation accuracy, dimension, gridtype,
// finite element type, used quadrature rules
typedef double RealType;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfigurationType;

/**
 * \brief This program computes \f$det(D\phi)\f$, \f$||D\phi||^2\f$ and \f$||cof(D\phi)||^2\f$ for a deformation \f$\phi\f$.
 * (Ie. determinant, squared Frobenius norm and squared Frobenius norm of the cofactor matrix of the deformation gradient are computed.)
 * The deformation is given as \f$\phi=identity+d\f$ for a displacement d (which can eg. be interpreted as displacement function of an image).
 * The input parameters of the program are the number of deformations (or images), the grid depth n
 * (where the grid on which the displacement is defined has \f$2^n+1\f$ nodes in each direction), the
 * directory in which the input and output files lie, the filename extension of the input files, and
 * the wordstem of the input files. The input files all represent the displacements and their file names
 * are given by "wordstem_i_j.fileextension", where i is the number of the displacement and j the direction.
 * The output files are "volumeEnergy.ext", "lengthEnergy.ext", "surfaceEnergy.ext", where "ext" stands
 * for "png" in 2D and "bz2" in 3D. They display the above mentioned three terms at all grid points.
 *
 * \author Wirth
 */
int main( int argc, char *argv[] ) {

  try {
    // read in number of images, deformation filename extension and stem of file names
    if ( argc < 6 ) { // wrong number of input arguments specified; explain correct syntax
      cerr << "Incorrect syntax. Use: ""computeElasticEnergies <number of images> <grid depth> <operating directory> <filename extension> <filename wordstem> [<weight filename wordstem>]""" << endl;
      return EXIT_FAILURE;
    } else {           // compute the elastic energies and save them
      // generate a corresponding grid
      qc::GridDefinition grid( atoi(argv[2]), ConfigurationType::Dim );

      int numberOfImages = atoi( argv[1] );
      // for each image...
      for ( int i=0; i<numberOfImages; i++ ) {
        //! read in the displacement "d"
        aol::MultiVector<RealType> d( ConfigurationType::Dim, grid.getNumberOfNodes() );
        for ( int j=0; j<ConfigurationType::Dim; j++ ) {
          // read in jth component of displacement "dComponent"
          ConfigurationType::ArrayType dComponent( grid );
          char fileName[1024];
          sprintf( fileName, "%s%s_%d_%d.%s", argv[3], argv[5], i, j, argv[4] );
          dComponent.load( fileName );
          d[j] = dComponent;
        }
        //! compute the hyperelastic energies (aux is just needed as a non-used input to "apply(...)")
        ConfigurationType::ArrayType volumeInvariant( grid ), lengthInvariant( grid ), surfaceInvariant( grid );
        aol::MultiVector<RealType> invariants;
        invariants.appendReference( lengthInvariant );
        invariants.appendReference( surfaceInvariant );
        invariants.appendReference( volumeInvariant );
        qc::DeformationGradientInvariants<ConfigurationType> invariantsEvaluator( grid );
        invariantsEvaluator.apply( d, invariants );
        //! inversely deform the corresponding energy-images and a chess pattern
        // define the energies and chess pattern as multivector
        qc::ScalarArray<RealType, qc::QC_2D> chessImage( volumeInvariant.getNumX(), volumeInvariant.getNumY() );
        generateChessPattern<RealType>( 64, chessImage );
        if ( ConfigurationType::Dim == 2 )
          invariants.appendReference( chessImage );

        if ( argc > 6 ) {
          char fileName[1024];
          sprintf( fileName, "%s%d.pgm", argv[6], i + 1 );
          ConfigurationType::ArrayType image( fileName );
          image /= image.getMaxValue();
          RealType fac[3] = { sqrt( static_cast<RealType>( ConfigurationType::Dim ) ), sqrt( static_cast<RealType>( ConfigurationType::Dim ) ), 1 };
          for ( int k = 0; k < invariants[0].size(); k++ ) {
            for ( int j = 0; j < 3; j++ ) {
              invariants[j][k] *= image[k];
              invariants[j][k] += fac[j] * ( 1 - image[k] );
            }
            if ( image[k] < .5 ) {
              chessImage[k] += .5 * 255;
              chessImage[k] /= 2.;
            }
          }
        }

        // compute inversely deformed energy-image
        aol::MultiVector<RealType> invDefInvariants( invariants, aol::STRUCT_COPY );
        qc::InvDeformImage<ConfigurationType>( invariants, grid, invDefInvariants, d );
/*        // assign the energies the new values in the deformed configuration
        lengthInvariant  = invDefInvariants[0];
        surfaceInvariant = invDefInvariants[1];
        volumeInvariant  = invDefInvariants[2];
*/

        cerr<<lengthInvariant.getMaxValue()<<" "<<lengthInvariant.getMinValue()<<endl;
        cerr<<surfaceInvariant.getMaxValue()<<" "<<surfaceInvariant.getMinValue()<<endl;
        cerr<<volumeInvariant.getMaxValue()<<" "<<volumeInvariant.getMinValue()<<endl;


/*lengthInvariant.addToAll( -sqrt( static_cast<RealType>( ConfigurationType::Dim ) ) );
volumeInvariant.addToAll( -1. );*/


        //! save the hyperelastic energies and the chess pattern
        // save energies as .bz2-file
        char fileName[1024];
        sprintf( fileName, "%svolumeInvariant_%d.bz2", argv[3], i );
        volumeInvariant.save( fileName, qc::PGM_DOUBLE_BINARY );
        sprintf( fileName, "%slengthInvariant_%d.bz2", argv[3], i );
        lengthInvariant.save( fileName, qc::PGM_DOUBLE_BINARY );
        sprintf( fileName, "%ssurfaceInvariant_%d.bz2", argv[3], i );
        surfaceInvariant.save( fileName, qc::PGM_DOUBLE_BINARY );
        if ( ConfigurationType::Dim == 2 ) {
          // save energies as .png-file
          char fileName[1024];
          qc::MultiArray<unsigned char,2,3> colorImage( grid );
          qc::ScalarArray<RealType, qc::QC_2D> volumeInvariant2d( volumeInvariant, grid );
          convertToColor( volumeInvariant2d, colorImage, 0.95, 1.05 );
          sprintf( fileName, "%svolumeInvariant_%d.png", argv[3], i );
          colorImage.savePNG( fileName );
          qc::ScalarArray<RealType, qc::QC_2D> lengthInvariant2d( lengthInvariant, grid );
          convertToColor( lengthInvariant2d, colorImage, .95 * sqrt( static_cast<RealType>( ConfigurationType::Dim ) ), 1.05 * sqrt( static_cast<RealType>( ConfigurationType::Dim ) ) );
          sprintf( fileName, "%slengthInvariant_%d.png", argv[3], i );
          colorImage.savePNG( fileName );
          qc::ScalarArray<RealType, qc::QC_2D> surfaceInvariant2d( surfaceInvariant, grid );
          convertToColor( surfaceInvariant2d, colorImage, .6 * sqrt( static_cast<RealType>( ConfigurationType::Dim ) ), 1.4 * sqrt( static_cast<RealType>( ConfigurationType::Dim ) ) );
          sprintf( fileName, "%ssurfaceInvariant_%d.png", argv[3], i );
          colorImage.savePNG( fileName );
          sprintf( fileName, "%sdefChessPat_%d.png", argv[3], i );
          chessImage = invDefInvariants[3];
          chessImage.savePNG( fileName );
          // save a colour-bar (volumeInvariant2d here is just an auxiliary variable)
          qc::DataGenerator<ConfigurationType>( grid ).generateLineLevelset( volumeInvariant2d, 0, 0. );
          convertToColor( volumeInvariant2d, colorImage );
          sprintf( fileName, "%scolourBar.png", argv[3] );
          colorImage.savePNG( fileName );
        }
      }
    }
  }
  catch ( aol::Exception &el ) {
    // print out error message
    el.dump();
  }

  // pause before closing the program (to let the user first read the output)
  aol::callSystemPauseIfNecessaryOnPlatform();
  return EXIT_SUCCESS;
}
