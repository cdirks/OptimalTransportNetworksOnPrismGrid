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

#ifndef __VTKFXVIEWER_H
#define __VTKFXVIEWER_H

#include <foxIncludes.h>
#include <vtkIncludes.h>

#include "vtkFXRenderWindowInteractor.h"
#include "FXVTKCanvas.h"
#include "vtkFXStereo.h"

// vtk equivalent of FXGLViewer -- some day, hopefully ...

#include "polyDataGroup.h"
#include "vtkFXUtils.h"

class PolyDataGroup;
class vtkFXGui;

class vtkFXViewer : public FXWindow {
  FXDECLARE ( vtkFXViewer )

public:
  enum { ID_BACKGROUND_COLOR = FXWindow::ID_LAST,
         ID_LOAD_VOLUME,
         ID_LOAD_UDPLY,
         ID_LOAD_STANFORDPLY,
         ID_LOAD_SRF,
         ID_LOAD_UDPLYSTRIPED,
         ID_LOAD_VTKPOLYDATA,
         ID_LOAD_TPCFE_PAR,
         ID_LOAD_LEVELSET,
         ID_LOAD_LEVELSET_TPCFE,
         ID_LOAD_LEVELSET_TEXTURED,
         ID_LOAD_LEVELSET_WITH_TEXTURECOORDS,
         ID_LOAD_GRAPH,
         ID_LOAD_LIDAR_POINT_CLOUD,
         ID_LOAD_OVERLAY_IMAGE_2D,
         ID_UPPER_OPACITY_VALUE,
         ID_LOWER_OPACITY_VALUE,
         ID_UPPER_TRANSFER_COLOR_VALUE,
         ID_LOWER_TRANSFER_COLOR_VALUE,
         ID_UPDATE_OVERLAY_IMAGE,
         ID_MAPPER_TEXTURE_2D,
         ID_MAPPER_TEXTURE_3D,
         ID_MAPPER_RAY_CAST_COMPOSITE,
         ID_MAPPER_RAY_CAST_MIP,
         ID_PNG_EXPORT_IMAGE,
         ID_PNG_HR_EXPORT_IMAGE,
         ID_PDF_EXPORT,
         ID_PLANE1_NORMALX,
         ID_PLANE1_NORMALY,
         ID_PLANE1_NORMALZ,
         ID_PLANE1_ORIGINX,
         ID_PLANE1_ORIGINY,
         ID_PLANE1_ORIGINZ,
         ID_TOGGLE_PLANE1,
         ID_MAKE_VIDEO,
         ID_LOAD_PLANE1,
         ID_SAVE_PLANE1,
         ID_LOAD_CAMERA,
         ID_SAVE_CAMERA,
         ID_RESET_CAMERA,
         ID_RESIZE_GLCANVAS,
         ID_REMOVE_ALL_PROPS_FROM_RENDERER,
         ID_REMOVE_MESH_FROM_RENDERER,
         ID_REMOVE_VOLUME_FROM_RENDERER,
         ID_STYLE_POINTS,
         ID_STYLE_WIREFRAME,
         ID_STYLE_SURFACE,
         ID_STYLE_FLAT,
         ID_STYLE_GOURAUD,
         ID_STYLE_PHONG,
         ID_STYLE_TOGGLE_SCALAR_VISIBILITY,
         ID_TOGGLE_GLOBAL_IMMEDIATE_MODE_RENDERING,
         ID_ISOLINE,
         ID_ISOSURFACE,
         ID_UPDATE_VISUAL,
         ID_UPDATE_CAMERA,
         ID_GET_CAMERA_ANGLES,
         ID_LIMIT_CAMERA_X_ANGLE,
         ID_LIGHT_SWITCH_UPDATE,
         ID_LIGHT_COLOR_UPDATE,
         ID_LIGHT_POSITION_UPDATE,
         ID_LAST
  };
  // TODO: sort these in some useful order

private:
  // think about how to sort these blocks in a useful way.
  // dialog stuff
  FXFileDialog _fileDialog;
  FXFileDialog _saveDialog;

  // color transfer for volume renderer
  FXfloat      _opacityFunctionUpperValue;
  FXfloat      _opacityFunctionLowerValue;
  FXfloat      _colorFunctionUpperValue;
  FXfloat      _colorFunctionLowerValue;
  vtkPiecewiseFunction  *_opacityFunction;
  vtkPiecewiseFunction  *_colorFunction;
  vtkVolume             *_volume;
  vtkVolumeMapper       *_volumeMapper;

  FXfloat _overlayImagePos[2];
  FXfloat _overlayImageOpacity;
  FXfloat _overlayImageSize;
public:
  FXDataTarget _overlayImageXPosTarget;
  FXDataTarget _overlayImageYPosTarget;
  FXDataTarget _overlayImageOpacityTarget;
  FXDataTarget _overlayImageSizeTarget;
private:
  vtkPoints *_overlayImagePoints;
  vtkPolyDataMapper2D *_overlayImageMapper;
  vtkTexturedActor2D *_overlayImageActor;

  // general color stuff
  FXColor      _backgroundColor;

  // clipping plane
  FXdouble     _plane1Normal[3];
  FXdouble     _plane1Origin[3];
  FXbool       _plane1Active;
  vtkPlane     *_plane1;

  FXdouble     _rotAngles[3];
  FXdouble     _zoom;
  vtkLight     *_light[8]; // openGL supports 8 light sources
  FXColor      _lightColor[8];
  FXdouble     _lightAngles[8][2];
  FXbool       _lightSwitch[8];

public:
  FXDataTarget _opacityFunctionUpperValueTarget;
  FXDataTarget _opacityFunctionLowerValueTarget;
  FXDataTarget _colorFunctionUpperValueTarget;
  FXDataTarget _colorFunctionLowerValueTarget;
  FXDataTarget _plane1NormalXTarget;
  FXDataTarget _plane1NormalYTarget;
  FXDataTarget _plane1NormalZTarget;
  FXDataTarget _plane1OriginXTarget;
  FXDataTarget _plane1OriginYTarget;
  FXDataTarget _plane1OriginZTarget;
  FXDataTarget _plane1ActiveTarget;
  FXDataTarget _backgroundColorTarget;
  FXDataTarget _rotXAngleTarget;
  FXDataTarget _rotYAngleTarget;
  FXDataTarget _rotZAngleTarget;
  FXDataTarget _zoomTarget;
  FXDataTarget *_lightColorTarget[8];
  FXDataTarget *_lightElevationTarget[8], *_lightAzimuthTarget[8];
  FXDataTarget *_lightSwitchTarget[8];

private:
  // general rendering
  vtkRenderWindow *_renWindow;
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  vtkStereoRenderer  *_renderer;
#else
  vtkRenderer        *_renderer;
#endif

  //std::auto_ptr<FXGLVisual> vis;
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  vtkRenderWindow    *_renWindowRight;
  FXVTKCanvasWithSecondCanvasUpdate   *_fxVTKCanvas;
  FXVTKCanvasWithSecondCanvasUpdate   *_fxVTKCanvasRight;
#else
  FXVTKCanvas                         *_fxVTKCanvas;
#endif
  FXGLVisual *_vis; // TODO: rename to _vis. think about whether we need an auto_ptr

  // data for poly renderer
  PolyDataGroup *_polyDataGroup;

public:
  // private std constructor
  // private copy constructor
  // private assignment operator

  vtkFXViewer ( FXComposite *Comp, vtkFXGui *MainWindow, FXApp *app );
  ~vtkFXViewer ( );

  void updateVTKCanvas();

  long onCmdUpdateVisual ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ );

  long onCmdUpperOpacityValue ( FXObject *sender, FXSelector sel, void *data );

  long onCmdUpdateCamera ( FXObject*, FXSelector, void* );
  long onCmdGetCameraAngles ( FXObject*, FXSelector, void* );
  long onCmdLimitCameraXAngle ( FXObject*, FXSelector, void* );

  long onCmdLightSwitchUpdate ( FXObject*, FXSelector, void* );
  long onCmdLightColorUpdate ( FXObject*, FXSelector, void* );
  long onCmdLightPositionUpdate ( FXObject*, FXSelector, void* );

  long onCmdChangeMapper ( FXObject* /*sender*/, FXSelector sel, void* /*data*/ );
  long onCmdBackgroundColor ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onCmdLoadPlane1 ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );
  long onCmdSavePlane1 ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onCmdLoadCamera ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );
  long onCmdSaveCamera ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );
  long onCmdResetCamera ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onCmdResizeGLCanvas ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onCmdRemoveAllPropsFromRenderer ( FXObject *sender, FXSelector sel, void *data );
  long onCmdRemoveObjectFromRenderer ( FXObject* /*sender*/, FXSelector sel, void* /*data*/ );

  long onCmdRenderStyle ( FXObject* /*sender*/, FXSelector sel, void* /*data*/ );
  long onUpdRenderStyle ( FXObject *sender, FXSelector sel, void* /*data*/ );

  void loadVolume ( const char * FileName, vtkImageData * & target, const bool LoadAsFloat = false );
  void loadTexture ( const char * Filename_x, const char * Filename_y, vtkImageData * & target );
  void loadAndShowVolume ( const char * FileName );
  void loadAndShowPLYStriped ( const char * FileName, const bool LeftView, const float StripeWidth );
  long onLoadVolumeClick ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long loadPlyClick ( const int Format );

  long onLoadUDPlyClick ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );
  long onLoadStanfordPlyClick ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onLoadSRFClick ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onLoadUDPlyStripedClick ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onLoadVTKPolyDataClick ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onLoadTPCFEParClick ( FXObject*, FXSelector, void* );

  long onLoadLevelsetTpcfeClick ( FXObject*, FXSelector, void* );

  long onLoadLevelsetClick ( FXObject*, FXSelector, void* );

  long onLoadLevelsetTexturedClick ( FXObject*, FXSelector, void* );

  long onLoadLevelsetWithTextureCoordsClick ( FXObject*, FXSelector, void* );

  long onLoadGraphClick ( FXObject*, FXSelector, void* );

  long onLoadLidarPointCloudClick ( FXObject*, FXSelector, void* );

  long onLoadOverlayImage2DClick ( FXObject*, FXSelector, void* ); 

  long onCmdUpdateOverlayImage ( FXObject*, FXSelector, void* ); 

  void exportCurrentViewToPng ( const char *FileName );

  void exportCurrentViewToHRPng ( const char *FileName );

  long onCmdIsoline  ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onCmdIsosurface ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ );

  long onExportPngClick ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );
  long onExportHRPngClick ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );
  long onExportPdfClick ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onCmdPlane1Normal ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onCmdTogglePlane1 ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  long onCmdMakeVideo ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ );

  void create_volume_pipeline( );

  void addMeshPolyDataToRenderer ( vtkPolyData* PolyData, const FXString & DataName );

  FXColor getLightColor ( const int i ) {
    return _lightColor[i];
  }

  PolyDataGroup *getPolyDataGroupPointer () {
    return _polyDataGroup;
  }

protected:
  vtkFXViewer() : FXWindow(), _fileDialog ( this, "Open" ), _saveDialog ( this, "Save Image" ),  _vis ( NULL )  { }

};


#endif
