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

#ifndef __VTKVISUALIZE_H
#define __VTKVISUALIZE_H

#include <scalarArray.h>
#include <triangMesh.h>
#include <meshWithData.h>
#include <vtkIncludes.h>
#include "triangMeshToVTK.h"

void updateCamera ( const double* rotAngles, double Zoom, vtkRenderer *Renderer );

template<typename QuocDataType>
void convertQCArrayToVtkImageData ( qc::Array<QuocDataType> &Array, vtkImageData *ImageData ) {
  ImageData->SetScalarTypeToUnsignedShort();
  const int numX = Array.getNumX();
  const int numY = Array.getNumY();
  const int numZ = Array.getNumZ();
  ImageData->SetDimensions ( numX, numY, numZ );
  for ( int x = 0; x < numX; x++ ) {
    for ( int y = 0; y < numY; y++ ) {
      for ( int z = 0; z < numZ; z++ ) {
        * ( static_cast< unsigned short* > ( ImageData->GetScalarPointer ( x, numY - 1 - y, z ) ) ) = static_cast<unsigned short> ( 65535 * Array.get ( x, y, z ) );
      }
    }
  }
}

template<typename QuocDataType>
void convertQCArrayToVtkImageDataFloat ( qc::Array<QuocDataType> &Array, vtkImageData *ImageData ) {
  ImageData->SetScalarTypeToFloat();
  const int numX = Array.getNumX();
  const int numY = Array.getNumY();
  const int numZ = Array.getNumZ();
  ImageData->SetDimensions ( numX, numY, numZ );
  for ( int x = 0; x < numX; x++ ) {
    for ( int y = 0; y < numY; y++ ) {
      for ( int z = 0; z < numZ; z++ ) {
        * ( static_cast< float* > ( ImageData->GetScalarPointer ( numX - 1 - x, numY - 1 - y, z ) ) ) = static_cast<float> ( Array.get ( x, y, z ) );
      }
    }
  }
}

template <typename RealType>
void initColorTransferFunction ( vtkColorTransferFunction *ColorTF, const int ColorMap, const RealType MinValue, const RealType MaxValue, const RealType R = 1, const RealType G = 1, const RealType B = 1 ) {

  switch ( ColorMap ) {
  // Red to blue HSV
  case 0:
    ColorTF->AddRGBPoint ( MinValue + 0.00 * ( MaxValue - MinValue ),  0.0, 0.0, 1.0 ); // blue
    ColorTF->AddRGBPoint ( MinValue + 0.25 * ( MaxValue - MinValue ),  0.0, 1.0, 1.0 ); // cyan
    ColorTF->AddRGBPoint ( MinValue + 0.50 * ( MaxValue - MinValue ),  0.0, 1.0, 0.0 ); // green
    ColorTF->AddRGBPoint ( MinValue + 0.75 * ( MaxValue - MinValue ),  1.0, 1.0, 0.0 ); // yellow
    ColorTF->AddRGBPoint ( MinValue + 1.00 * ( MaxValue - MinValue ),  1.0, 0.0, 0.0 ); // red
    break;

  // Red to blue
  case 1:
    ColorTF->AddRGBPoint ( MinValue + 0.00 * ( MaxValue - MinValue ),  0.0, 0.0, 1.0 ); // blue
    ColorTF->AddRGBPoint ( MinValue + 1.00 * ( MaxValue - MinValue ),  1.0, 0.0, 0.0 ); // red
    break;

  // topographical
  case 2:
    ColorTF->AddRGBPoint ( MinValue +  0. / 11. * ( MaxValue - MinValue ), 0.0 , 0.0 , 0.5  );
    ColorTF->AddRGBPoint ( MinValue +  1. / 11. * ( MaxValue - MinValue ), 0.0 , 0.0 , 1.0  );
    ColorTF->AddRGBPoint ( MinValue +  2. / 11. * ( MaxValue - MinValue ), 0.4 , 0.4 , 1.0  );
    ColorTF->AddRGBPoint ( MinValue +  3. / 11. * ( MaxValue - MinValue ), 0.6 , 0.6 , 1.0  );
    ColorTF->AddRGBPoint ( MinValue +  4. / 11. * ( MaxValue - MinValue ), 0.3 , 1.0 , 0.3  );
    ColorTF->AddRGBPoint ( MinValue +  5. / 11. * ( MaxValue - MinValue ), 0.0 , 1.0 , 0.0  );
    ColorTF->AddRGBPoint ( MinValue +  6. / 11. * ( MaxValue - MinValue ), 0.0 , 0.5 , 0.0  );
    ColorTF->AddRGBPoint ( MinValue +  7. / 11. * ( MaxValue - MinValue ), 0.3 , 0.2 , 0.0  );
    ColorTF->AddRGBPoint ( MinValue +  8. / 11. * ( MaxValue - MinValue ), 0.5 , 0.3 , 0.0  );
    ColorTF->AddRGBPoint ( MinValue +  9. / 11. * ( MaxValue - MinValue ), 0.8 , 0.5 , 0.0  );
    ColorTF->AddRGBPoint ( MinValue + 10. / 11. * ( MaxValue - MinValue ), 1.0 , 0.7 , 0.5  );
    ColorTF->AddRGBPoint ( MinValue + 11. / 11. * ( MaxValue - MinValue ), 1.0 , 0.85, 0.75 );
    break;

  // blue to red HSV
  case 3:
    ColorTF->AddRGBPoint ( MinValue + 0.00 * ( MaxValue - MinValue ),  1.0, 0.0, 0.0 ); // red
    ColorTF->AddRGBPoint ( MinValue + 0.25 * ( MaxValue - MinValue ),  1.0, 1.0, 0.0 ); // yellow
    ColorTF->AddRGBPoint ( MinValue + 0.50 * ( MaxValue - MinValue ),  0.0, 1.0, 0.0 ); // green
    ColorTF->AddRGBPoint ( MinValue + 0.75 * ( MaxValue - MinValue ),  0.0, 1.0, 1.0 ); // cyan
    ColorTF->AddRGBPoint ( MinValue + 1.00 * ( MaxValue - MinValue ),  0.0, 0.0, 1.0 ); // blue
    break;

  case 4:
    ColorTF->AddRGBPoint ( MinValue + 0.00 * ( MaxValue - MinValue ),  0.0, 0.0, 0.0 ); // black
    ColorTF->AddRGBPoint ( MinValue + 1.00 * ( MaxValue - MinValue ),  1.0, 1.0, 1.0 ); // white
    break;

  case 5:
    ColorTF->AddRGBPoint ( MinValue + 0.00 * ( MaxValue - MinValue ),  1.0, 1.0, 1.0 ); // white
    ColorTF->AddRGBPoint ( MinValue + 1.00 * ( MaxValue - MinValue ),  0.0, 0.0, 0.0 ); // black
    break;

  case 6:
    ColorTF->AddRGBPoint ( MinValue + 0.00 * ( MaxValue - MinValue ),  1.0, 1.0, 1.0 ); // white
    ColorTF->AddRGBPoint ( MinValue + 1.00 * ( MaxValue - MinValue ),  0.0, 0.0, 1.0 ); // blue
    break;

  case 7:
    ColorTF->AddRGBPoint ( MinValue + 0.00 * ( MaxValue - MinValue ),  1.0, 1.0, 1.0 ); // white
    ColorTF->AddRGBPoint ( MinValue + 0.33 * ( MaxValue - MinValue ),  0.0, 1.0, 0.0 ); // green
    ColorTF->AddRGBPoint ( MinValue + 0.66 * ( MaxValue - MinValue ),  0.0, 0.0, 1.0 ); // blue
    ColorTF->AddRGBPoint ( MinValue + 1.00 * ( MaxValue - MinValue ),  1.0, 0.0, 0.0 ); // red
    break;

  case 8:
    ColorTF->AddRGBPoint ( MinValue + 0.00 * ( MaxValue - MinValue ),  1.0, 0.8, 0.0 ); // silver
    ColorTF->AddRGBPoint ( MinValue + 1.00 * ( MaxValue - MinValue ),  0.9, 0.9, 0.85 ); // gold
    break;

  case 9:
    ColorTF->AddRGBPoint ( MinValue + 0.00 * ( MaxValue - MinValue ),  R, G, B );
    ColorTF->AddRGBPoint ( MinValue + 1.00 * ( MaxValue - MinValue ),  R, G, B );
    break;

  default:
    cerr << "Invalid color transition" << endl;
  }
}

template <typename RealType>
class GraphToPolyData {
  const qc::ScalarArray<RealType, qc::QC_2D> _imageArray;
  qc::ScalarArray<RealType, qc::QC_2D> _imageDataArray;
  vtkPolyData *_polyData;
  vtkWarpScalar *_warp;
  vtkTexture *_texture;
  vtkTextureMapToPlane *_texturePlane;
public:
  const bool _loadGraphTexture, _loadGraphScalarData;
  //RealType _arrayMinValue, _arrayMaxValue;

  GraphToPolyData ( const char *GraphFilename, const char *GraphTextureFilename )
    : _imageArray ( GraphFilename ),
      _polyData ( NULL ),
      _warp ( vtkWarpScalar::New() ),
      _texture ( vtkTexture::New() ),
      _texturePlane ( vtkTextureMapToPlane::New() ),
      _loadGraphTexture ( GraphTextureFilename != NULL ),
      _loadGraphScalarData ( _loadGraphTexture && ( aol::fileNameEndsWith ( GraphTextureFilename, ".png" ) == false ) ) {

    vtkPlaneSource *plane = vtkPlaneSource::New();
    plane->SetResolution (_imageArray.getNumX(),_imageArray.getNumY());

    vtkTransform *transform = vtkTransform::New();
    //transform->Scale(10.0,10.0,1.0);
    transform->Translate(0.5,0.5,0.0);

    vtkTransformPolyDataFilter *transF = vtkTransformPolyDataFilter::New();
    transF->SetInputConnection(plane->GetOutputPort());
    transF->SetTransform(transform);
    transF->Update();

    // compute Bessel function and derivatives. This portion could be
    // encapsulated into source or filter object.
    //
    vtkPolyData *input = transF->GetOutput();
    int numPts = input->GetNumberOfPoints();
    vtkPoints *newPts = vtkPoints::New();
    newPts->SetNumberOfPoints(numPts);

    vtkFloatArray *derivs = vtkFloatArray::New();
    derivs->SetNumberOfTuples(numPts);

    vtkPolyData *bessel = vtkPolyData::New();
    bessel->CopyStructure(input);
    bessel->SetPoints(newPts);
    bessel->GetPointData()->SetScalars(derivs);

    if ( _loadGraphScalarData )
      _imageDataArray.load( GraphTextureFilename );

    double x[3];
    for (int i=0; i<numPts; i++){
      input->GetPoint(i,x);
      // UNUSED! const double r = sqrt(static_cast<double>(x[0]*x[0]) + x[1]*x[1]);
      x[2] = 0.5 * _imageArray.interpolate(x[0]*(_imageArray.getNumX()-1),(1-x[1])*(_imageArray.getNumY()-1));
      newPts->SetPoint(i,x);
      derivs->SetValue(i, _loadGraphScalarData ? _imageDataArray.interpolate(x[0]*(_imageDataArray.getNumX()-1),(1-x[1])*(_imageDataArray.getNumY()-1)) : x[2] );
    }

    // warp plane
    _warp->SetInput(bessel);
    _warp->XYPlaneOn();
    _warp->SetScaleFactor(1);
    if ( hasTexture() == false ) {
      _polyData = _warp->GetPolyDataOutput();
    }
    else
    {
      // Display a texture instead of scalar data.
      bessel->GetPointData()->SetScalars(NULL);
      // For some reason the z-scaling gets messed up if we remove the scalars and _warp->SetScaleFactor seems to have no effect.
      // Compensate for this manually.
      for (int i=0; i<numPts; i++){
        newPts->GetPoint(i,x);
        x[2] *= 2;
        newPts->SetPoint(i, x);
      }
      _texturePlane->AutomaticPlaneGenerationOff();
      _texturePlane->SetSRange ( 0, 1 );
      _texturePlane->SetTRange ( 0, 1 );
      _texturePlane->SetOrigin ( 0, 0, 0 );
      _texturePlane->SetPoint1 ( 1, 0, 0 );
      _texturePlane->SetPoint2 ( 0, 1, 0 );
      _texturePlane->SetInputConnection(_warp->GetOutputPort());
      vtkPNGReader *pngReader = vtkPNGReader::New();
      pngReader->SetFileName ( GraphTextureFilename );
      _texture->SetInputConnection(pngReader->GetOutputPort());

      _polyData = _texturePlane->GetPolyDataOutput();
      pngReader->Delete();
    }

     //reference counting - it's ok
    plane->Delete();
    transform->Delete();
    transF->Delete();
    newPts->Delete();
    derivs->Delete();
    bessel->Delete();
  }

  ~GraphToPolyData ( ) {
    _texturePlane->Delete();
    _texture->Delete();
    // If warp is deleted here, there is an error on exit.
    _warp->Delete();
  }

  vtkPolyData *getPolyDataPointer ( ) const {
    return _polyData;
  }

  vtkTexture *getTexturePointer ( ) const {
    return _texture;
  }

  bool hasTexture ( ) const {
    return ( _loadGraphScalarData == false ) && ( _loadGraphTexture == true );
  }

  bool hasScalarData ( ) const {
    return _loadGraphScalarData;
  }

  const qc::ScalarArray<RealType, qc::QC_2D> &getImageArrayReference ( ) const {
    return _imageArray;
  }

  const qc::ScalarArray<RealType, qc::QC_2D> &getImageDataArrayReference ( ) const {
    return _imageDataArray;
  }
};

/**
 * \author Berkels
 */
template <typename RealType>
class PLYToPolyData {
  vtkPolyData *_polyData;
  RealType _minValue, _maxValue;
public:

  PLYToPolyData ( const char *PLYFilename, const aol::PLY_FORMAT Format )
    : _polyData ( NULL ),
      _minValue ( aol::NumberTrait<RealType>::NaN ),
      _maxValue ( aol::NumberTrait<RealType>::NaN ) {
    aol::TriangMesh<RealType> mesh;
    mesh.loadFromStanfordOrUDPLY ( PLYFilename, Format );
    _polyData = aol::TriangMeshToVTK<RealType>::makePolyData ( mesh, _minValue, _maxValue );
  }

  ~PLYToPolyData ( ) {
    _polyData->Delete();
  }

  vtkPolyData *getPolyDataPointer ( ) const {
    return _polyData;
  }

  RealType getMinValue ( ) const {
    return _minValue;
  }

  RealType getMaxValue ( ) const {
    return _maxValue;
  }
};


/**
 * \author Heeren
 */
template <typename RealType>
class IsoSurfaceToPolyData {  
  vtkPolyData*      _polyData;  
  vtkImageData*     _levelData;
  vtkContourFilter* _contourFilter; 

public:

  IsoSurfaceToPolyData ( const char *FileName, RealType levelValueToDisplay )
    : _polyData ( NULL ),
      _levelData( vtkImageData::New() ),
      _contourFilter( vtkContourFilter::New() ){
	
    qc::ScalarArray<RealType, qc::QC_3D> volumeArray ( FileName );
    convertQCArrayToVtkImageDataFloat<RealType> ( volumeArray, _levelData );

    _contourFilter->SetInput ( _levelData );
    _contourFilter->SetValue ( 0, levelValueToDisplay );
  
    _polyData = _contourFilter->GetOutput();  
  }

  ~IsoSurfaceToPolyData ( ) { 
    _contourFilter->Delete(); 
    _levelData->Delete();                  
  }

  vtkPolyData *getPolyDataPointer ( ) const {
    return _polyData;
  }

};

/**
 * Based on VTK/Examples/Cxx/Utilities/OffScreenRendering.
 *
 * \todo Code partially overlaps with PolyDataQC. Separate the non-GUI code of PolyDataQC and use this here.
 *
 * \author Berkels
 */
class OffscreenPNGRenderer {
  vtkGraphicsFactory *_graphicsFactory;
  vtkImagingFactory *_imagingFactory;
  vtkActor *_actor;
  vtkPolyDataMapper *_mapper;
  vtkProperty *_property;
  vtkPolyDataNormals *_normals;
  vtkStripper *_stripper;
  aol::Vec2<int> _outputImageSize;
  aol::Vec3<double> _viewAngles;
public:
  OffscreenPNGRenderer ( )
    : _graphicsFactory ( vtkGraphicsFactory::New() ),
      _imagingFactory ( vtkImagingFactory::New() ),
      _actor ( vtkActor::New() ),
      _mapper ( vtkPolyDataMapper::New() ),
      _property ( vtkProperty::New() ),
      _normals ( vtkPolyDataNormals::New() ),
      _stripper ( vtkStripper::New() ),
      _outputImageSize ( 1000, 1000 ) {
    _graphicsFactory->SetOffScreenOnlyMode( 1);
    _graphicsFactory->SetUseMesaClasses( 1 );
    _imagingFactory->SetUseMesaClasses( 1 );
    _actor->SetMapper ( _mapper );
    _actor->SetProperty ( _property );
    _normals->SetFeatureAngle ( 60.0 );
}

  ~OffscreenPNGRenderer ( ) {
    _actor->Delete();
    _mapper->Delete();
    _property->Delete();
    _normals->Delete();
    _stripper->Delete();
    _imagingFactory->Delete();
    _graphicsFactory->Delete();
  }

  vtkPolyDataMapper *getMapperPointer ( ) const {
    return _mapper;
  }

  vtkActor *getActorPointer ( ) const {
    return _actor;
  }

  void usePolyData ( vtkPolyData *PolyData ) {
    _normals->SetInput ( PolyData );
    _stripper->SetInputConnection ( _normals->GetOutputPort() );
    _mapper->SetInputConnection ( _stripper->GetOutputPort() );
  }

  void setOutputImageSize ( const int SizeX, const int SizeY ) {
    _outputImageSize.set ( SizeX, SizeY );
  }

  void setViewAngles ( const aol::Vec3<double> &ViewAngles ) {
    _viewAngles = ViewAngles;
  }

  void renderViewToPNG ( const char *Filename ) const {
    vtkRenderer *renderer = vtkRenderer::New();
    vtkRenderWindow *renderWindow = vtkRenderWindow::New();
    renderWindow->SetSize ( _outputImageSize[0], _outputImageSize[1] );
    renderWindow->SetOffScreenRendering( 1 );
    renderWindow->AddRenderer(renderer);

    // Add the actor to the scene.
    renderer->AddActor(_actor);
    // Set background color to white.
    renderer->SetBackground(1,1,1);

    if ( _viewAngles.normSqr() != 0 )
      updateCamera ( _viewAngles.getData(), 1, renderer );

    renderWindow->Render();

    vtkWindowToImageFilter *windowToImageFilter = vtkWindowToImageFilter::New();
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->Update();

    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetFileName(Filename);
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();

    writer->Delete();
    windowToImageFilter->Delete();
    renderWindow->Delete();
    renderer->Delete();
  }

  template <typename RealType>
  void setColorTransferFunction ( const int ColorMap, const RealType MinValue, const RealType MaxValue, const RealType R = 1, const RealType G = 1, const RealType B = 1 ) {
    vtkColorTransferFunction *colorTF = vtkColorTransferFunction::New();
    initColorTransferFunction<RealType> ( colorTF, ColorMap, MinValue, MaxValue, R, G, B );
    getMapperPointer()->SetLookupTable ( colorTF );
    getMapperPointer()->SetColorModeToMapScalars();
    getMapperPointer()->InterpolateScalarsBeforeMappingOn();
    colorTF->Delete();
  }
};

#endif
