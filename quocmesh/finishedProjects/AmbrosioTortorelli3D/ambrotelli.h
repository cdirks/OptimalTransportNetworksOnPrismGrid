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

#ifndef __AMBROTELLI_H
#define __AMBROTELLI_H

#include <FEOpInterface.h>
#include <deformations.h>

#ifdef GUI_PRESENT
template <typename RealType, qc::Dimension dim>
class FXImageUpdater{
  FXApp *_appWindow;
  std::vector<FXImageFrame*> _imageFrame;
  std::vector<FXImage*> _image;
  std::vector<FXColor*> _imageData;
  const int _width;
public:
  FXImageUpdater( FXApp *AppWindow,
                  FXImageFrame *ImageFrame,
                  FXImage *Image,
                  FXColor *ImageData,
                  const int Width )
  :
    _appWindow(AppWindow),
    _imageFrame(1),
    _image(1),
    _imageData(1),
    _width(Width)
  {
    _imageFrame[0] = ImageFrame;
    _image[0] = Image;
    _imageData[0] = ImageData;
  }
  FXImageUpdater( FXApp *AppWindow,
                  std::vector<FXImageFrame*> &ImageFrame,
                  std::vector<FXImage*> &Image,
                  std::vector<FXColor*> &ImageData,
                  const int Width )
  :
    _appWindow(AppWindow),
    _imageFrame(ImageFrame),
    _image(Image),
    _imageData(ImageData),
    _width(Width)
  {
  }

  virtual ~FXImageUpdater( ) {
  }
  void fillFXColorWithScalarArray( FXColor *Color, const qc::GridDefinition &Grid, const aol::Vector<RealType> &Image ){
    if ( dim == qc::QC_2D ){
      FXuchar *pp = reinterpret_cast<FXuchar*>(Color);

      qc::ScalarArray<RealType, qc::QC_2D>* imageArray;
      if( Grid.getWidth() != _width ){
        qc::ScalarArray<RealType, qc::QC_2D> originalImageArray(Image, Grid);
        imageArray = new qc::ScalarArray<RealType, qc::QC_2D>(_width, _width);
        imageArray->resampleFrom(originalImageArray);
      }
      else{
        imageArray = new qc::ScalarArray<RealType, qc::QC_2D>(Image, Grid);
      }

      for( int i = 0; i < _width; i++){
        for( int j = 0; j < _width; j++,pp+=4){
          // Attention: x and y coordinates of a ScalarArray<QC_2D> and an FXImage are different!
          RealType value = aol::Min<RealType>( aol::Max<RealType>( imageArray->get(j,i), 0.), 1.);
          unsigned char temp = static_cast<unsigned char>( value*255.0);
          pp[0] = temp;
          pp[1] = temp;
          pp[2] = temp;
          pp[3] = 255;
        }
      }
      delete imageArray;
    }
    else
     cerr << "On screen rendering of 3D images is not supported.\n";
  }
  void fillFXColorWithTwoBlendedScalarArrays( FXColor *Color, const qc::GridDefinition &GridA, const aol::Vector<RealType> &ImageA, const qc::GridDefinition &GridB, const aol::Vector<RealType> &ImageB ){
    if ( dim == qc::QC_2D ){
      FXuchar *pp = reinterpret_cast<FXuchar*>(Color);

      qc::ScalarArray<RealType, qc::QC_2D>* imageArrayA;
      qc::ScalarArray<RealType, qc::QC_2D>* imageArrayB;
      if( GridA.getWidth() != _width ){
        qc::ScalarArray<RealType, qc::QC_2D> originalImageArrayA(ImageA, GridA);
        imageArrayA = new qc::ScalarArray<RealType, qc::QC_2D>(_width, _width);
        imageArrayA->resampleFrom(originalImageArrayA);
      }
      else{
        imageArrayA = new qc::ScalarArray<RealType, qc::QC_2D>(ImageA, GridA);
      }
      if( GridB.getWidth() != _width ){
        qc::ScalarArray<RealType, qc::QC_2D> originalImageArrayB(ImageB, GridB);
        imageArrayB = new qc::ScalarArray<RealType, qc::QC_2D>(_width, _width);
        imageArrayB->resampleFrom(originalImageArrayB);
      }
      else{
        imageArrayB = new qc::ScalarArray<RealType, qc::QC_2D>(ImageB, GridB);
      }

      for( int i = 0; i < _width; i++){
        for( int j = 0; j < _width; j++,pp+=4){
          const RealType valueA = aol::Min<RealType>( aol::Max<RealType>( imageArrayA->get(j,i), 0.), 1.);
          const RealType valueB = aol::Min<RealType>( aol::Max<RealType>( imageArrayB->get(j,i), 0.), 1.);
          // Attention: x and y coordinates of a ScalarArray<QC_2D> and an FXImage are different!
          const unsigned char tempA = static_cast<unsigned char>( valueA*255.0);
          const unsigned char tempB = static_cast<unsigned char>( valueB*255.0);
          pp[0] = tempA;
          pp[1] = tempB;
          pp[2] = tempA;
          pp[3] = 255;
        }
      }
      delete imageArrayA;
      delete imageArrayB;
    }
    else
     cerr << "On screen rendering of 3D images is not supported.\n";
  }
  void fillImageWithScalarArray( const int ImageNumber, const qc::GridDefinition &Grid, const aol::Vector<RealType> &Image ){
    fillFXColorWithScalarArray( _imageData[ImageNumber], Grid, Image );
  }
  void fillImageWithTwoBlendedScalarArrays( const int ImageNumber, const qc::GridDefinition &GridA, const aol::Vector<RealType> &ImageA, const qc::GridDefinition &GridB, const aol::Vector<RealType> &ImageB ){
    fillFXColorWithTwoBlendedScalarArrays( _imageData[ImageNumber], GridA, ImageA, GridB, ImageB );
  }
  void renderImage( const int ImageNumber ){
    _image[ImageNumber]->setData(_imageData[ImageNumber]);
    _image[ImageNumber]->render();
    _imageFrame[ImageNumber]->setImage(_image[ImageNumber]);
    _appWindow->repaint();
  }
  void updateImage( const qc::GridDefinition &Grid, const aol::Vector<RealType> &Image, const int ImageNumber = 0 ){
    fillFXColorWithScalarArray( _imageData[ImageNumber], Grid, Image );
    renderImage( ImageNumber );
  }
};

template <typename RealType, qc::Dimension dim>
void renderScalarArray( FXImage *GUIImage, FXColor *ImageData, const qc::GridDefinition &Grid, const aol::Vector<RealType> &Image, const int Max_depth ){
  if ( dim == qc::QC_2D ){
    FXuchar *pp = reinterpret_cast<FXuchar*>(ImageData);
    const int width = static_cast<int>(pow(2., Max_depth))+1;

    qc::ScalarArray<RealType, qc::QC_2D> imageArray(Image, Grid);
    qc::ScalarArray<RealType, qc::QC_2D> scaledImageArray(width, width);
    scaledImageArray.scale(imageArray);

    for( int i = 0; i < width; i++){
      for( int j = 0; j < width; j++,pp+=4){
        // Attention: x and y coordinates of a ScalarArray<QC_2D> and an FXImage are different!
        unsigned char temp = static_cast<unsigned char>( scaledImageArray.get(j,i)*255.0);
        pp[0] = temp;
        pp[1] = temp;
        pp[2] = temp;
        pp[3] = 255;
      }
    }
    GUIImage->setData(ImageData);
    GUIImage->render();
  }
  else
   cerr << "On screen rendering of 3D images is not supported.\n";
}
#endif

template <typename RealType, qc::Dimension dim>
void generateIdentity( const qc::GridDefinition &Grid, aol::MultiVector<RealType> &identity ){
  const RealType h= Grid.H();
  qc::GridDefinition::OldFullNodeIterator fnit;
  qc::FastILexMapper<dim> mapper(Grid);
  for ( fnit = Grid.begin(); fnit != Grid.end(); ++fnit ){
    for( int i = 0; i < dim; i++){
      identity[i][ mapper.getGlobalIndex(*fnit) ] = (*fnit)[i] * h;
    }
  }
}

#endif
