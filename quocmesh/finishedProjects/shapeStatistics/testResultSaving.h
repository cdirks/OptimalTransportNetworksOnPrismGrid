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

#ifndef __TESTRESULTSAVING_H
#define __TESTRESULTSAVING_H

#include <aol.h>
#include <FEOpInterface.h>
#include <deformations.h>
#include <multiVector.h>

/**************************************************************
 * Auxiliary functions for testing.
 **************************************************************/

enum ImageSaveType {
  PGM,
  PNG
};

// only implemented for 2D
template<typename ConfiguratorType>
class StressRotator :
  public aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim,StressRotator<ConfiguratorType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const aol::DiscreteFunctionDefault<ConfiguratorType> _gradField;

public:
  StressRotator( const typename ConfiguratorType::InitType &Grid, const aol::Vector<RealType> &GradField ) :
    aol::FENonlinEvaluationOpInterface<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim,ConfiguratorType::Dim*ConfiguratorType::Dim,StressRotator<ConfiguratorType> >( Grid ),
    _gradField( Grid, GradField ) {}

  void getNonlinearity ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim*ConfiguratorType::Dim> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         typename aol::Vec<ConfiguratorType::Dim*ConfiguratorType::Dim,RealType> &NL ) const {
    NL.setZero();
    // compute normal \nu and stress \sigma
    aol::Vec<ConfiguratorType::Dim,RealType> normal;
    _gradField.evaluateGradientAtQuadPoint( El, QuadPoint, normal );
    aol::Vec<ConfiguratorType::Dim*ConfiguratorType::Dim,RealType> stress;
    DiscFunc.evaluateAtQuadPoint( El, QuadPoint, stress );
    // compute Q, which rotates the canonic first basis vector into \nu
    typename ConfiguratorType::MatType Q;
    Q.setCol( 0, normal );
    Q.set( 0, 1, -normal[1] );
    Q.set( 1, 1, normal[0] );
    // compute Q^T\sigma Q
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        for ( int k = 0; k < ConfiguratorType::Dim; k++ )
          for ( int l = 0; l < ConfiguratorType::Dim; l++ )
            NL[i*ConfiguratorType::Dim+j] += Q[k][i] * stress[k*ConfiguratorType::Dim+l] * Q[l][j];
  }
};

//! This function produces a 2d-image "ChessImage" with a chess pattern of "NumberOfStripes" horizontal or vertical stripes.
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

//! This function produces a 3d-image "ChessImage" with a chess pattern of "NumberOfStripes" horizontal or vertical stripes.
template <typename PixelValueType>
void generateChessPattern( const int NumberOfStripes, qc::ScalarArray<PixelValueType, qc::QC_3D> &ChessImage ){
  const int numX = ChessImage.getNumX();
  const int numY = ChessImage.getNumY();
  const int numZ = ChessImage.getNumZ();
  // compute stripe-width in pixels
  const int stripeWidth = static_cast<int>( numZ/NumberOfStripes );
  // make image black
  ChessImage.setZero();
  // for each pixel do...
  for ( int z = 0; z < numZ; z++ )
    for ( int y = 0; y < numY; y++ )
      for ( int x = 0; x < numX; x++ )
        // make white field
        if ( (y%(2*stripeWidth)<stripeWidth) == (z%(2*stripeWidth)<stripeWidth) )
          if ( x%(2*stripeWidth) < stripeWidth )
            ChessImage.set( x, y, 255 );
          else {}
        else if ( x%(2*stripeWidth) > stripeWidth )
            ChessImage.set( x, y, 255 );
}

//! This function converts the 2d-image "Image" (with ideally values in [0,1]) into a colored image "ColorImage", which can be saved as png-file.
template <typename RealType>
void convertToColor( const aol::Vector<RealType> &Image, qc::MultiArray<unsigned char,2,3> &ColorImage, const RealType Min = 0, const RealType Max = 1 ){
  // specify colors at certain values
  const int partitions = 4;
  RealType intervals[partitions+1] = {0, .25, .5, .75, 1};
  const RealType colorsRGB[partitions+1][3] = {{0,0,.5}, {0,.35,.15}, {.35,.35,.35}, {1,.35,.1}, {1,1,0}};

  // expand color intervals over the desired range
  for ( int i = 0; i <= partitions; i++ )
    intervals[i] = Min + intervals[i] * ( Max - Min );

  RealType rgb[3] = {0, 0, 0};

  // for each pixel do...
  for ( int i = 0; i < Image.size(); i++ ) {
    // get pixel value
    RealType b = Image[i];
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
    ColorImage[0][i] = static_cast< unsigned char> ( rgb[0] * 255 );
    ColorImage[1][i] = static_cast< unsigned char> ( rgb[1] * 255 );
    ColorImage[2][i] = static_cast< unsigned char> ( rgb[2] * 255 );
  }
}

template <typename ConfiguratorType>
inline void save( const aol::Vector<typename ConfiguratorType::RealType> &Image, const typename ConfiguratorType::InitType &Grid, const char* DestDir, const char* Name, const int Index = -1, const bool Mute = true, const bool Compress = ConfiguratorType::Dim == qc::QC_2D ? false : true ) {
  typename ConfiguratorType::ArrayType image( Image, Grid, aol::DEEP_COPY );
  char filename[1024];
  if ( Index >= 0 )
    sprintf( filename, "%s/%s%d.%s", DestDir, Name, Index, Compress ? "bz2" : "pgm" );
  else
    sprintf( filename, "%s/%s.%s", DestDir, Name, Compress ? "bz2" : "pgm" );
  typename ConfiguratorType::RealType min = Image.getMinValue(), max = Image.getMaxValue();
  if ( !Mute )
    cerr<<Name<<" "<<min<<" to "<<max<<endl;
  image.save( filename, qc::PGM_DOUBLE_BINARY );
}

template <typename ConfiguratorType, ImageSaveType saveType>
class ImageSaver {};

template <typename ConfiguratorType>
class ImageSaver<ConfiguratorType,PNG> {
public:
  inline static void saveImage( const typename ConfiguratorType::ArrayType &Image, const char* Name ) { Image.savePNG( Name ); }
};

template <typename ConfiguratorType>
class ImageSaver<ConfiguratorType,PGM> {
public:
  inline static void saveImage( const typename ConfiguratorType::ArrayType &Image, const char* Name ) { Image.save( Name, qc::PGM_UNSIGNED_CHAR_BINARY ); }
};

template <typename ConfiguratorType, ImageSaveType saveType>
inline void saveAsImage( const aol::Vector<typename ConfiguratorType::RealType> &Image, const typename ConfiguratorType::InitType &Grid, const char* DestDir, const char* Name, const int Index = -1, const typename ConfiguratorType::RealType LowerBound = 0, const typename ConfiguratorType::RealType UpperBound = 0, const bool Mute = true ) {
  typename ConfiguratorType::ArrayType image( Image, Grid, aol::DEEP_COPY );
  char filename[1024];
  if ( Index >= 0 )
    sprintf( filename, "%s/%s%d.%s", DestDir, Name, Index, saveType == PGM ? "pgm" : "png" );
  else
    sprintf( filename, "%s/%s.%s", DestDir, Name, saveType == PGM ? "pgm" : "png" );

  typename ConfiguratorType::RealType min = Image.getMinValue(), max = Image.getMaxValue();
  if ( !Mute )
    cerr<<Name<<" "<<min<<" to "<<max<<endl;
  if ( LowerBound < UpperBound ) {
    min = LowerBound;
    max = UpperBound;
  }
  image.addToAll( -min );
  image *= 255 / ( max - min );
  ImageSaver<ConfiguratorType,saveType>::saveImage( image, filename );
}

template <typename ConfiguratorType>
inline void saveAsColorImage( const aol::Vector<typename ConfiguratorType::RealType> &Image, const typename ConfiguratorType::InitType &Grid, const char* DestDir, const char* Name, const int Index = -1, const typename ConfiguratorType::RealType LowerBound = 0, const typename ConfiguratorType::RealType UpperBound = 0, const bool Mute = true ) {
  typename ConfiguratorType::ArrayType image( Image, Grid, aol::DEEP_COPY );
  char filename[1024];
  if ( Index >= 0 )
    sprintf( filename, "%s/%s%d.png", DestDir, Name, Index );
  else
    sprintf( filename, "%s/%s.png", DestDir, Name );

  typename ConfiguratorType::RealType min = Image.getMinValue(), max = Image.getMaxValue();
  if ( !Mute )
    cerr<<Name<<" "<<min<<" to "<<max<<endl;
  if ( LowerBound < UpperBound ) {
    min = LowerBound;
    max = UpperBound;
  }
  image.addToAll( -min );
  image *= 1 / ( max - min );

  qc::GridSize<ConfiguratorType::Dim> gridSize( Grid );
  qc::MultiArray<unsigned char,2,3> colorImage( gridSize );
  convertToColor( image, colorImage, 0., 1. );
  colorImage.savePNG( filename );
}

template <typename ConfiguratorType>
inline void saveDisplacement( const aol::MultiVector<typename ConfiguratorType::RealType> &Disp, const typename ConfiguratorType::InitType &Grid, const char* DestDir, const char* Name, const int Index = -1, const bool Mute = true ) {
  for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
    char filename[1024];
    if ( Index >= 0 )
      sprintf( filename, "%s%d", Name, Index );
    else
      sprintf( filename, "%s", Name );
    save<ConfiguratorType>( Disp[i], Grid, DestDir, filename, i, Mute, true );
  }
}

template <typename ConfiguratorType, ImageSaveType saveType>
inline void saveDisplacementAsChessPattern( const aol::MultiVector<typename ConfiguratorType::RealType> &Disp, const typename ConfiguratorType::InitType &Grid, const char* DestDir, const char* Name, const int Index = -1, const bool Mute = true ) {
  typename ConfiguratorType::ArrayType chessPattern( Grid ), defChessPattern( Grid );
  generateChessPattern<typename ConfiguratorType::RealType>( 10, chessPattern );
  qc::InvDeformImage<ConfiguratorType>( chessPattern, Grid, defChessPattern, Disp );
  saveAsImage<ConfiguratorType,saveType>( defChessPattern, Grid, DestDir, Name, Index, Mute );
}

template <typename ConfiguratorType>
inline void saveStressTensor( const aol::MultiVector<typename ConfiguratorType::RealType> &Stress, const typename ConfiguratorType::InitType &Grid, const char* DestDir, const char* Name, const int Index = -1, const bool Mute = true ) {
  for ( int i = 0; i < ConfiguratorType::Dim; i++ )
    for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
      char filename[1024];
      sprintf( filename, "%s_%d%d", Name, i, j );
      save<ConfiguratorType>( Stress[i*ConfiguratorType::Dim+j], Grid, DestDir, filename, Index, Mute );
    }
}

template <typename ConfiguratorType>
inline void saveStressTensor( const aol::MultiVector<typename ConfiguratorType::RealType> &Stress, const typename ConfiguratorType::InitType &Grid, const char* DestDir, typename ConfiguratorType::ArrayType &GradientField, const char* Name, const bool Mute = true ) {
  aol::MultiVector<typename ConfiguratorType::RealType> rotStress( Stress, aol::STRUCT_COPY );
  StressRotator<ConfiguratorType>( Grid, GradientField ).apply( Stress, rotStress );
  saveStressTensor<ConfiguratorType>( rotStress, Grid, DestDir, Name, -1, Mute );
}


#endif
