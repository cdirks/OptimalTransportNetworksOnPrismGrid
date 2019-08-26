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

#include "vtkFXViewer.h"
#include <smallVec.h>
#include <scalarArray.h>
#include <triangMesh.h>
#include <levelsetToTriangMesh.h>
#include "triangMeshToVTK.h"
#include "tpcfeToVTK.h"
#include "vtkVisualize.h"


void vtkFXViewerKeypressCallbackFunction ( vtkObject* caller, long unsigned int /*eventId*/, void* clientData, void* /*callData*/ )
{
  vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
  if ( iren->GetKeyCode() == 'n' ) {

    FXVTKCanvas* canvas = static_cast<FXVTKCanvas*> ( clientData );

    // For some reason the mouse coordinates of the event are wrong, perhaps a bug in FXVTKCanvas.
    // To work around this, we shift the obtained coordinates based on the position of the canvas in the
    // main window (assuming that the windows returned by getParent() is a child of the main window).
    const int mouseX = iren->GetEventPosition()[0] - canvas->getX() + canvas->getParent()->getX();
    const int mouseY = iren->GetEventPosition()[1] + canvas->getY() + canvas->getParent()->getY();

    vtkPointPicker *picker = vtkPointPicker::New();
    cerr << "Picking at mouse position " << mouseX << " " << mouseY << endl;

    if ( picker->Pick ( mouseX, mouseY, 0, iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer() ) )
    {
      aol::Vec<3, double> pos ( picker->GetMapperPosition() );
      cerr << "Picked " << setprecision ( 20 ) << pos[0] << " " << pos[1] << " " << pos[2] << endl;
    }
    picker->Delete();
  }
}

FXDEFMAP ( vtkFXViewer ) vtkFXViewerMap[] =
  {
    FXMAPFUNCS ( SEL_COMMAND, vtkFXViewer::ID_UPPER_OPACITY_VALUE,              vtkFXViewer::ID_LOWER_TRANSFER_COLOR_VALUE, vtkFXViewer::onCmdUpperOpacityValue ),
    FXMAPFUNCS ( SEL_CHANGED, vtkFXViewer::ID_UPPER_OPACITY_VALUE,              vtkFXViewer::ID_LOWER_TRANSFER_COLOR_VALUE, vtkFXViewer::onCmdUpperOpacityValue ),
    FXMAPFUNCS ( SEL_COMMAND, vtkFXViewer::ID_MAPPER_TEXTURE_2D,                vtkFXViewer::ID_MAPPER_RAY_CAST_MIP, vtkFXViewer::onCmdChangeMapper ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_VOLUME,                      vtkFXViewer::onLoadVolumeClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_UDPLY,                       vtkFXViewer::onLoadUDPlyClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_STANFORDPLY,                 vtkFXViewer::onLoadStanfordPlyClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_SRF,                         vtkFXViewer::onLoadSRFClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_UDPLYSTRIPED,                vtkFXViewer::onLoadUDPlyStripedClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_VTKPOLYDATA,                 vtkFXViewer::onLoadVTKPolyDataClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_TPCFE_PAR,                   vtkFXViewer::onLoadTPCFEParClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_LEVELSET_TPCFE,              vtkFXViewer::onLoadLevelsetTpcfeClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_LEVELSET,                    vtkFXViewer::onLoadLevelsetClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_LEVELSET_TEXTURED,           vtkFXViewer::onLoadLevelsetTexturedClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_LEVELSET_WITH_TEXTURECOORDS, vtkFXViewer::onLoadLevelsetWithTextureCoordsClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_GRAPH,                       vtkFXViewer::onLoadGraphClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_LIDAR_POINT_CLOUD,           vtkFXViewer::onLoadLidarPointCloudClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_OVERLAY_IMAGE_2D,            vtkFXViewer::onLoadOverlayImage2DClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_UPDATE_OVERLAY_IMAGE,             vtkFXViewer::onCmdUpdateOverlayImage ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_PNG_EXPORT_IMAGE,                 vtkFXViewer::onExportPngClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_PNG_HR_EXPORT_IMAGE,              vtkFXViewer::onExportHRPngClick ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_PDF_EXPORT,                       vtkFXViewer::onExportPdfClick ),
    FXMAPFUNCS ( SEL_CHANGED, vtkFXViewer::ID_PLANE1_NORMALX,                   vtkFXViewer::ID_PLANE1_ORIGINZ, vtkFXViewer::onCmdPlane1Normal ),
    FXMAPFUNCS ( SEL_COMMAND, vtkFXViewer::ID_PLANE1_NORMALX,                   vtkFXViewer::ID_PLANE1_ORIGINZ, vtkFXViewer::onCmdPlane1Normal ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_TOGGLE_PLANE1,                    vtkFXViewer::onCmdTogglePlane1 ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_MAKE_VIDEO,                       vtkFXViewer::onCmdMakeVideo ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_PLANE1,                      vtkFXViewer::onCmdLoadPlane1 ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_SAVE_PLANE1,                      vtkFXViewer::onCmdSavePlane1 ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LOAD_CAMERA,                      vtkFXViewer::onCmdLoadCamera ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_SAVE_CAMERA,                      vtkFXViewer::onCmdSaveCamera ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_RESET_CAMERA,                     vtkFXViewer::onCmdResetCamera ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_BACKGROUND_COLOR,                 vtkFXViewer::onCmdBackgroundColor ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_RESIZE_GLCANVAS,                  vtkFXViewer::onCmdResizeGLCanvas ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_REMOVE_ALL_PROPS_FROM_RENDERER,   vtkFXViewer::onCmdRemoveAllPropsFromRenderer ),
    FXMAPFUNCS ( SEL_COMMAND, vtkFXViewer::ID_REMOVE_MESH_FROM_RENDERER,        vtkFXViewer::ID_REMOVE_VOLUME_FROM_RENDERER, vtkFXViewer::onCmdRemoveObjectFromRenderer ),
    FXMAPFUNCS ( SEL_COMMAND, vtkFXViewer::ID_STYLE_POINTS,                     vtkFXViewer::ID_TOGGLE_GLOBAL_IMMEDIATE_MODE_RENDERING, vtkFXViewer::onCmdRenderStyle ),
    FXMAPFUNCS ( SEL_UPDATE,  vtkFXViewer::ID_STYLE_POINTS,                     vtkFXViewer::ID_TOGGLE_GLOBAL_IMMEDIATE_MODE_RENDERING, vtkFXViewer::onUpdRenderStyle ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_ISOLINE,                          vtkFXViewer::onCmdIsoline ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_ISOSURFACE,                       vtkFXViewer::onCmdIsosurface ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_UPDATE_VISUAL,                    vtkFXViewer::onCmdUpdateVisual ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_UPDATE_CAMERA,                    vtkFXViewer::onCmdUpdateCamera ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_GET_CAMERA_ANGLES,                vtkFXViewer::onCmdGetCameraAngles ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LIMIT_CAMERA_X_ANGLE,             vtkFXViewer::onCmdLimitCameraXAngle ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LIGHT_SWITCH_UPDATE,              vtkFXViewer::onCmdLightSwitchUpdate ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LIGHT_COLOR_UPDATE,               vtkFXViewer::onCmdLightColorUpdate ),
    FXMAPFUNC  ( SEL_COMMAND, vtkFXViewer::ID_LIGHT_POSITION_UPDATE,            vtkFXViewer::onCmdLightPositionUpdate )
  };

FXIMPLEMENT ( vtkFXViewer, FXWindow, vtkFXViewerMap, ARRAYNUMBER ( vtkFXViewerMap ) )


vtkFXViewer::vtkFXViewer ( FXComposite *Comp, vtkFXGui *MainWindow, FXApp *app ) :
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  FXWindow ( MainWindow, DECOR_ALL, 0, 0, 3200, 1200 ),
#else
  FXWindow ( MainWindow, DECOR_ALL, 0, 0,  800,  600 ),
#endif
  _fileDialog ( this, "Open" ),
  _saveDialog ( this, "Save Image" ),
  _opacityFunctionUpperValue ( 1.0 ),
  _opacityFunctionLowerValue ( 0. ),
  _colorFunctionUpperValue ( 1.0 ),
  _colorFunctionLowerValue ( 0. ),
  _opacityFunction ( NULL ),
  _colorFunction ( NULL ),
  _volume ( vtkVolume::New() ),
  _volumeMapper ( vtkOpenGLVolumeTextureMapper2D::New() ),
  _overlayImageOpacity ( 0.5 ),
  _overlayImageSize ( 200 ),
  _overlayImageXPosTarget ( _overlayImagePos[0], this, ID_UPDATE_OVERLAY_IMAGE ),
  _overlayImageYPosTarget ( _overlayImagePos[1], this, ID_UPDATE_OVERLAY_IMAGE ),
  _overlayImageOpacityTarget ( _overlayImageOpacity, this, ID_UPDATE_OVERLAY_IMAGE ),
  _overlayImageSizeTarget ( _overlayImageSize, this, ID_UPDATE_OVERLAY_IMAGE ),
  _overlayImagePoints ( vtkPoints::New() ),
  _overlayImageMapper ( vtkPolyDataMapper2D::New() ),
  _overlayImageActor ( vtkTexturedActor2D::New() ),
  _backgroundColor ( FXRGB ( 0, 0, 0 ) ),
  _plane1Active ( false ),
  _plane1 ( vtkPlane::New() ),
  // FXDataTargets
  _opacityFunctionUpperValueTarget ( _opacityFunctionUpperValue, this, ID_UPPER_OPACITY_VALUE ),
  _opacityFunctionLowerValueTarget ( _opacityFunctionLowerValue, this, ID_LOWER_OPACITY_VALUE ),
  _colorFunctionUpperValueTarget ( _colorFunctionUpperValue, this, ID_UPPER_TRANSFER_COLOR_VALUE ),
  _colorFunctionLowerValueTarget ( _colorFunctionLowerValue, this, ID_LOWER_TRANSFER_COLOR_VALUE ),
  _plane1NormalXTarget ( _plane1Normal[0], this, ID_PLANE1_NORMALX ),
  _plane1NormalYTarget ( _plane1Normal[1], this, ID_PLANE1_NORMALY ),
  _plane1NormalZTarget ( _plane1Normal[2], this, ID_PLANE1_NORMALZ ),
  _plane1OriginXTarget ( _plane1Origin[0], this, ID_PLANE1_ORIGINX ),
  _plane1OriginYTarget ( _plane1Origin[1], this, ID_PLANE1_ORIGINY ),
  _plane1OriginZTarget ( _plane1Origin[2], this, ID_PLANE1_ORIGINZ ),
  _plane1ActiveTarget ( _plane1Active, this, ID_TOGGLE_PLANE1 ),
  _backgroundColorTarget ( _backgroundColor, this, ID_BACKGROUND_COLOR ),
  _rotXAngleTarget  ( _rotAngles[0], this, ID_UPDATE_CAMERA ),
  _rotYAngleTarget  ( _rotAngles[1], this, ID_UPDATE_CAMERA ),
  _rotZAngleTarget  ( _rotAngles[2], this, ID_UPDATE_CAMERA ),
  _zoomTarget ( _zoom, this, ID_UPDATE_CAMERA ),
  // end of data targets
  _renWindow ( NULL ),
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  _renderer ( new vtkStereoRenderer ),
  _renWindowRight ( NULL ),
#else
  _renderer ( vtkRenderer::New() ),
#endif
  _vis ( new FXGLVisual (app, VISUAL_DOUBLEBUFFER | VISUAL_STEREO ) ),
  _polyDataGroup ( new PolyDataGroup ( MainWindow, _renderer ) )
{
  _overlayImagePos[0] = 0.5;
  _overlayImagePos[1] = 0.5;

  _plane1Normal[0] = 1.;
  _plane1Normal[1] = 1.;
  _plane1Normal[2] = 0.;
  _plane1Origin[0] = 0.;
  _plane1Origin[1] = 0.;
  _plane1Origin[2] = 0.;
  _plane1->SetOrigin ( _plane1Origin );
  _plane1->SetNormal ( _plane1Normal );

  _rotAngles[0] =  0.0;
  _rotAngles[1] =  0.0;
  _rotAngles[2] =  0.0;

  _zoom = 1.0;

  _renderer->AutomaticLightCreationOff();

  for ( int i = 0; i < 8; ++i ) {
    _lightColor[i] = FXRGB ( 255, 255, 255 );
    _lightAngles[i][0] = 0.0;
    _lightAngles[i][1] = 0.0;
    _lightSwitch[i] = ( i == 0 );
    _lightColorTarget[i]     = new FXDataTarget ( _lightColor[i],     this, ID_LIGHT_COLOR_UPDATE );
    _lightElevationTarget[i] = new FXDataTarget ( _lightAngles[i][0], this, ID_LIGHT_POSITION_UPDATE );
    _lightAzimuthTarget[i]   = new FXDataTarget ( _lightAngles[i][1], this, ID_LIGHT_POSITION_UPDATE );
    _lightSwitchTarget[i]    = new FXDataTarget ( _lightSwitch[i],    this, ID_LIGHT_SWITCH_UPDATE );

    _light[i] = vtkLight::New(); //_renderer->MakeLight();
    _light[i]->SetLightTypeToCameraLight();
    _light[i]->SetDiffuseColor(  FXREDVAL ( _lightColor[i]  ) / 255.0, FXGREENVAL ( _lightColor[i]  ) / 255.0, FXBLUEVAL ( _lightColor[i]  ) / 255.0 );

    if ( i != 0 )
      _light[i]->SwitchOff();

    _renderer->AddLight( _light[i] );
  }

  FXHorizontalFrame *frame = new FXHorizontalFrame ( Comp, LAYOUT_FILL_X | LAYOUT_FILL_Y );
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  _fxVTKCanvas = new FXVTKCanvasWithSecondCanvasUpdate ( frame, _vis, NULL, 0, LAYOUT_FILL_X | LAYOUT_FILL_Y );
  _fxVTKCanvasRight = new FXVTKCanvasWithSecondCanvasUpdate ( frame, _vis, NULL, 0, LAYOUT_FILL_X | LAYOUT_FILL_Y );
  _fxVTKCanvas->setSecondVTKCanvasToUpdate ( _fxVTKCanvasRight );
  _fxVTKCanvasRight->setSecondVTKCanvasToUpdate ( _fxVTKCanvas );
#else
  _fxVTKCanvas = new FXVTKCanvas ( frame, _vis, NULL, 0, LAYOUT_FILL_X | LAYOUT_FILL_Y );
#endif

  vtkCallbackCommand *keypressCallback = vtkCallbackCommand::New();
  keypressCallback->SetCallback ( vtkFXViewerKeypressCallbackFunction );
  keypressCallback->SetClientData ( _fxVTKCanvas );
  _fxVTKCanvas->getInteractor ()->AddObserver ( vtkCommand::KeyPressEvent, keypressCallback );
  keypressCallback->Delete();

  create_volume_pipeline( );
  _renderer->SetBackground ( FXREDVAL ( _backgroundColor ) / 255.0, FXGREENVAL ( _backgroundColor ) / 255.0, FXBLUEVAL ( _backgroundColor ) / 255.0 );

}


vtkFXViewer::~vtkFXViewer() {
  // TODO: delete in correct order. are these all to be deleted?
  _opacityFunction->Delete();
  _colorFunction->Delete();
  _volume->Delete();
  _volumeMapper->Delete();
  _overlayImagePoints->Delete();
  _overlayImageMapper->Delete();
  _overlayImageActor->Delete();
  // _polyDataGroup depends on the renderer, so we have to delete it before deleting the renderer.
  delete _polyDataGroup;
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  delete _renderer;
#else
  _renderer->Delete();
#endif
  _renWindow->Delete();
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  _renWindowRight->Delete();
#endif
  _plane1->Delete();
  for (int i = 0; i < 8; ++i) {
    delete _lightColorTarget[i];
    delete _lightElevationTarget[i];
    delete _lightAzimuthTarget[i];
    delete _lightSwitchTarget[i];
    _light[i]->Delete();
  }
  delete _vis;
}


void vtkFXViewer::updateVTKCanvas() {
  _fxVTKCanvas->render();
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  _fxVTKCanvasRight->render();
#endif
}



long vtkFXViewer::onCmdUpdateVisual ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  updateVTKCanvas();
  return 1;
}

//! Small helper function for vtkFXViewer::onCmdUpperOpacityValue.
void setPiecewiseFunction ( vtkPiecewiseFunction *PiecewiseFunction, FXDataTarget &LowerValueTarget, FXDataTarget &UpperValueTarget ) {
  PiecewiseFunction->RemoveAllPoints();
  const FXfloat lowerVal = *( reinterpret_cast<FXfloat*>(LowerValueTarget.getData()) );
  const FXfloat upperVal = *( reinterpret_cast<FXfloat*>(UpperValueTarget.getData()) );
  PiecewiseFunction->AddPoint ( lowerVal * 65535, 0 );
  PiecewiseFunction->AddPoint ( upperVal * 65535, 1 );
}

long vtkFXViewer::onCmdUpperOpacityValue ( FXObject* /*sender*/, FXSelector sel, void* /*data*/ ) {

  switch ( FXSELID ( sel ) ) {
  case ID_UPPER_OPACITY_VALUE: {
    setPiecewiseFunction ( _opacityFunction, _opacityFunctionLowerValueTarget, _opacityFunctionUpperValueTarget );
    break;
  }
  case ID_LOWER_OPACITY_VALUE: {
    setPiecewiseFunction ( _opacityFunction, _opacityFunctionLowerValueTarget, _opacityFunctionUpperValueTarget );
    break;
  }
  case ID_UPPER_TRANSFER_COLOR_VALUE: {
    setPiecewiseFunction ( _colorFunction, _colorFunctionLowerValueTarget, _colorFunctionUpperValueTarget );
    break;
  }
  case ID_LOWER_TRANSFER_COLOR_VALUE: {
    setPiecewiseFunction ( _colorFunction, _colorFunctionLowerValueTarget, _colorFunctionUpperValueTarget );
    break;
  }
  default:
    throw aol::UnimplementedCodeException( "Selected FXSelector value not handled", __FILE__, __LINE__);
    break;
  }
  updateVTKCanvas();

  return 1;

}

long vtkFXViewer::onCmdUpdateCamera ( FXObject*, FXSelector, void* ) {
  // TODO: think about whether to treat each direction of rotation differently.

  cerr << _rotAngles[0] << " " << _rotAngles[1] << " " << _rotAngles[2] << endl;

#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  updateCamera ( _rotAngles, _zoom, _renderer->GetLeftRenderer() );
#else
  updateCamera ( _rotAngles, _zoom, _renderer );
#endif
  updateVTKCanvas(); // necessary??

  return 1;
}

long vtkFXViewer::onCmdGetCameraAngles ( FXObject*, FXSelector, void* ) {

  // OS is not sure whether this works correctly.

  vtkCamera *Camera = NULL;
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  Camera = _renderer->GetLeftRenderer()->GetActiveCamera();
#else
  Camera = _renderer->GetActiveCamera();
#endif

  FXdouble pos[3], foc[3], vup[3];
  Camera->GetPosition (  pos );
  //    cerr << "Pos: " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  Camera->GetFocalPoint ( foc );
  //    cerr << "Foc: " << foc[0] << " " << foc[1] << " " << foc[2] << endl;
  Camera->GetViewUp (    vup );
  //    cerr << "Vup: " << vup[0] << " " << vup[1] << " " << vup[2] << endl;

  FXdouble x = pos[0] - foc[0],  y = pos[1] - foc[1],  z = pos[2] - foc[2];

  cerr << "Co: " << x << " " << y << " " << z << endl;

  FXdouble phi = 0, theta = 0;

  if ( z > 0 )
    phi = atan ( x / z );
  else if ( z < 0 )
    phi = M_PI + atan ( x / z ); //  + M_PI;
  else if ( z == 0 ) {
    if ( x > 0 )
      phi = 0.5 * M_PI;
    else if ( x < 0 )
      phi = -0.5 * M_PI;
    else if ( x == 0 )
      phi = aol::NumberTrait<FXdouble>::NaN;
  }

  theta = 0.5 * M_PI - atan ( sqrt ( x * x + z * z ) / y );

  cerr << phi << " " << theta << " angles" << endl;

  aol::Vec3<FXdouble> yax ( 0.0, 1.0, 0.0 ), vup2 ( vup[0], vup[1], vup[2] );

  FXdouble elevation = theta;
  FXdouble azimuth   = phi;
  FXdouble roll      = acos ( ( vup2 * yax ) / ( vup2.norm() * yax.norm() ) );

  _rotAngles[0] = 180 * elevation / M_PI;
  _rotAngles[1] = 180 * azimuth   / M_PI;
  _rotAngles[2] = 180 * roll      / M_PI;

  // we still have the following problems: roll is not recovered correctly, theta may only be between 0 and 90 degrees

  return 1;
}

long vtkFXViewer::onCmdLimitCameraXAngle ( FXObject*, FXSelector, void* ) {
  if ( _rotAngles[0] < -90 ) {
    _rotAngles[0] = -180 - _rotAngles[0];
    _rotAngles[1] += 180;
  }
  else if ( _rotAngles[0] > 90 ) {
    _rotAngles[0] = 180 - _rotAngles[0];
    _rotAngles[1] += 180;
  }
  this->handle ( this, FXSEL ( SEL_COMMAND, ID_UPDATE_CAMERA ), NULL );
  return 1;
}

long vtkFXViewer::onCmdLightSwitchUpdate ( FXObject*, FXSelector, void* ) {

  for ( int i = 0; i < 8; ++i )
    if( _lightSwitch[i] )
      _light[i]->SwitchOn();
    else
      _light[i]->SwitchOff();

  updateVTKCanvas();

  return 1;
}

long vtkFXViewer::onCmdLightColorUpdate ( FXObject*, FXSelector, void* ) {

  for ( int i = 0; i < 8; ++i )
    _light[i]->SetDiffuseColor(  FXREDVAL ( _lightColor[i] ) / 255.0, FXGREENVAL ( _lightColor[i] ) / 255.0, FXBLUEVAL ( _lightColor[i] ) / 255.0 );

  updateVTKCanvas();

  return 1;
}

long vtkFXViewer::onCmdLightPositionUpdate ( FXObject*, FXSelector, void* ) {
  for ( int i = 0; i < 8; ++i )
    _light[i]->SetDirectionAngle( _lightAngles[i] );

  updateVTKCanvas();

  return 1;
}

long vtkFXViewer::onCmdChangeMapper ( FXObject* /*sender*/, FXSelector sel, void* /*data*/ ) {
  vtkVolumeMapper *newVolumeMapper = NULL;
  switch ( FXSELID ( sel ) ) {
  case ID_MAPPER_TEXTURE_2D: {
    newVolumeMapper = vtkOpenGLVolumeTextureMapper2D::New();
    break;
  }
  case ID_MAPPER_TEXTURE_3D: {
    newVolumeMapper = vtkOpenGLVolumeTextureMapper3D::New();
    break;
  }
  default:
    throw aol::UnimplementedCodeException( "Selected FXSelector value not handled", __FILE__, __LINE__);
    break;
  }
  if ( FXSELID ( sel ) == ID_MAPPER_RAY_CAST_COMPOSITE || FXSELID ( sel ) == ID_MAPPER_RAY_CAST_MIP ) {
    vtkVolumeRayCastMapper *mapper = vtkVolumeRayCastMapper::New();
    vtkVolumeRayCastFunction *rayCastFunction;
    if ( FXSELID ( sel ) == ID_MAPPER_RAY_CAST_COMPOSITE ) {
      rayCastFunction = vtkVolumeRayCastCompositeFunction::New();
    } else {
      rayCastFunction = vtkVolumeRayCastMIPFunction::New();
    }
    // Rendering artifacts of the vtkVolumeRayCastMapper can be reduced by reducing the sample distance.
    //mapper->SetSampleDistance ( 0.1 );
    mapper->SetVolumeRayCastFunction ( rayCastFunction );
    rayCastFunction->Delete();
    newVolumeMapper = mapper;
  }
  newVolumeMapper->SetInput ( _volumeMapper->GetInput() );
  if ( _plane1Active ) {
    newVolumeMapper->AddClippingPlane ( _plane1 );
  }
  // Now that we added the input of the old mapper to the new one, we don't need the old mapper anymore.
  _volumeMapper->Delete();
  _volumeMapper = newVolumeMapper;
  _volume->SetMapper ( _volumeMapper );
  updateVTKCanvas();
  return 1;
}

long vtkFXViewer::onCmdBackgroundColor ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ ) {
  _renderer->SetBackground ( FXREDVAL ( _backgroundColor ) / 255.0, FXGREENVAL ( _backgroundColor ) / 255.0, FXBLUEVAL ( _backgroundColor ) / 255.0 );
  return 1;
}

long vtkFXViewer::onCmdLoadPlane1 ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  ifstream infile ( "normal.dat", ofstream::binary );
  infile.read ( ( char* ) this->_plane1Normal, 3*sizeof ( FXdouble ) );
  infile.read ( ( char* ) this->_plane1Origin, 3*sizeof ( FXdouble ) );
  infile.close();
  this->handle ( this, FXSEL ( SEL_COMMAND, ID_PLANE1_NORMALX ), NULL );
  return 1;
}

long vtkFXViewer::onCmdSavePlane1 ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  ofstream outfile ( "normal.dat", ofstream::binary );
  outfile.write ( ( char* ) this->_plane1Normal, 3*sizeof ( FXdouble ) );
  outfile.write ( ( char* ) this->_plane1Origin, 3*sizeof ( FXdouble ) );
  outfile.close();
  return 1;
}

long vtkFXViewer::onCmdLoadCamera ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  loadCamera ( _renderer->GetLeftRenderer() );
#else
  loadCamera ( _renderer );
#endif
  updateVTKCanvas();

  return 1;
}

long vtkFXViewer::onCmdSaveCamera ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  saveCamera ( _renderer->GetLeftRenderer() );
#else
  saveCamera ( _renderer );
#endif
  return 1;
}

long vtkFXViewer::onCmdResetCamera ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  _renderer->ResetCamera();
  _renderer->GetActiveCamera()->SetViewUp ( 0., 1., 0. );

  updateVTKCanvas();
  return 1;
}

long vtkFXViewer::onCmdResizeGLCanvas ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  FXDialogBox resizeGLCanvas ( this, tr ( "Resize GLCanvas" ), DECOR_TITLE | DECOR_BORDER, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
  FXHorizontalFrame side ( &resizeGLCanvas, LAYOUT_SIDE_RIGHT | LAYOUT_FILL_X | LAYOUT_FILL_Y, 0, 0, 0, 0, 10, 10, 10, 10, 0, 0 );
  int widthInt = _fxVTKCanvas->getWidth();
  int heightInt = _fxVTKCanvas->getHeight();
  new FXLabel ( &side, "width", NULL, JUSTIFY_LEFT | ICON_BEFORE_TEXT | LAYOUT_FILL_X );
  FXDataTarget widthTarget ( widthInt );
  new FXTextField ( &side, 10, &widthTarget, FXDataTarget::ID_VALUE, TEXTFIELD_INTEGER | JUSTIFY_RIGHT | LAYOUT_CENTER_Y | LAYOUT_CENTER_X | FRAME_SUNKEN | FRAME_THICK | LAYOUT_FILL_ROW );
  new FXLabel ( &side, "height", NULL, JUSTIFY_LEFT | ICON_BEFORE_TEXT | LAYOUT_FILL_X );
  FXDataTarget heightTarget ( heightInt );
  new FXTextField ( &side, 10, &heightTarget, FXDataTarget::ID_VALUE, TEXTFIELD_INTEGER | JUSTIFY_RIGHT | LAYOUT_CENTER_Y | LAYOUT_CENTER_X | FRAME_SUNKEN | FRAME_THICK | LAYOUT_FILL_ROW );
  FXButton button ( &side, tr ( "&OK" ), NULL, &resizeGLCanvas, FXDialogBox::ID_ACCEPT, BUTTON_DEFAULT | FRAME_RAISED | FRAME_THICK | LAYOUT_RIGHT, 0, 0, 0, 0, 32, 32, 2, 2 );
  button.setFocus();
  resizeGLCanvas.execute ( PLACEMENT_OWNER );

  _fxVTKCanvas->resize ( widthInt, heightInt );
  this->layout();
  this->forceRefresh();
  this->repaint();
  return 1;
}

long vtkFXViewer::onCmdRemoveAllPropsFromRenderer ( FXObject *sender, FXSelector sel, void *data ) {
  _renderer->RemoveAllViewProps();
  _polyDataGroup->onCmdClear( sender, sel, data);
  updateVTKCanvas();
  return 1;
}

long vtkFXViewer::onCmdRemoveObjectFromRenderer ( FXObject* /*sender*/, FXSelector sel, void* /*data*/ ) {
  switch ( FXSELID ( sel ) ) {
  case ID_REMOVE_MESH_FROM_RENDERER: {
    FXMessageBox errorMsg ( this, "Sorry for this inconvenience", "Remove Mesh is not implemented yet!", NULL, MBOX_OK | DECOR_TITLE | DECOR_BORDER );
    errorMsg.execute();
    //if( _renderer->HasViewProp( _meshActor ) )
    //  _renderer->RemoveViewProp(_meshActor);
    break;
  }
  case ID_REMOVE_VOLUME_FROM_RENDERER: {
    if ( _renderer->HasViewProp ( _volume ) )
      _renderer->RemoveViewProp ( _volume );
    break;
  }
  default:
    throw aol::UnimplementedCodeException( "Selected FXSelector value not handled", __FILE__, __LINE__);
    break;
  }
  updateVTKCanvas();
  return 1;
}

long vtkFXViewer::onCmdRenderStyle ( FXObject* /*sender*/, FXSelector sel, void* /*data*/ ) {
  vtkProperty* meshProperty = _polyDataGroup->getCurrentProperty();
  if ( meshProperty ) {
    switch ( FXSELID ( sel ) ) {
    case ID_STYLE_POINTS: meshProperty->SetRepresentationToPoints(); break;
    case ID_STYLE_WIREFRAME: meshProperty->SetRepresentationToWireframe(); break;
    case ID_STYLE_SURFACE: meshProperty->SetRepresentationToSurface(); break;
    case ID_STYLE_FLAT: meshProperty->SetInterpolationToFlat(); break;
    case ID_STYLE_GOURAUD: meshProperty->SetInterpolationToGouraud(); break;
    case ID_STYLE_PHONG: meshProperty->SetInterpolationToPhong(); break;
    case ID_STYLE_TOGGLE_SCALAR_VISIBILITY: _polyDataGroup->getCurrentPolyDataMapper()->SetScalarVisibility ( 1 - _polyDataGroup->getCurrentPolyDataMapper()->GetScalarVisibility() ); break;
    case ID_TOGGLE_GLOBAL_IMMEDIATE_MODE_RENDERING: _polyDataGroup->getCurrentPolyDataMapper()->SetGlobalImmediateModeRendering ( 1 - _polyDataGroup->getCurrentPolyDataMapper()->GetGlobalImmediateModeRendering() ); break;
    default:
      throw aol::UnimplementedCodeException( "Selected FXSelector value not handled", __FILE__, __LINE__);
      break;
    }
  }
  updateVTKCanvas();
  return 1;
}

long vtkFXViewer::onUpdRenderStyle ( FXObject *sender, FXSelector sel, void* /*data*/ ) {
  vtkProperty* meshProperty = _polyDataGroup->getCurrentProperty();
  if ( meshProperty ) {
    FXSelector msg = FXSEL ( SEL_COMMAND, FXWindow::ID_UNCHECK );
    switch ( FXSELID ( sel ) ) {
    case ID_STYLE_POINTS: if ( meshProperty->GetRepresentation( ) == VTK_POINTS ) msg = FXSEL ( SEL_COMMAND, FXWindow::ID_CHECK ); break;
    case ID_STYLE_WIREFRAME: if ( meshProperty->GetRepresentation( ) == VTK_WIREFRAME ) msg = FXSEL ( SEL_COMMAND, FXWindow::ID_CHECK ); break;
    case ID_STYLE_SURFACE: if ( meshProperty->GetRepresentation( ) == VTK_SURFACE ) msg = FXSEL ( SEL_COMMAND, FXWindow::ID_CHECK ); break;
    case ID_STYLE_FLAT: if ( meshProperty->GetInterpolation( ) == VTK_FLAT ) msg = FXSEL ( SEL_COMMAND, FXWindow::ID_CHECK ); break;
    case ID_STYLE_GOURAUD: if ( meshProperty->GetInterpolation( ) == VTK_GOURAUD ) msg = FXSEL ( SEL_COMMAND, FXWindow::ID_CHECK ); break;
    case ID_STYLE_PHONG: if ( meshProperty->GetInterpolation( ) == VTK_PHONG ) msg = FXSEL ( SEL_COMMAND, FXWindow::ID_CHECK ); break;
    case ID_STYLE_TOGGLE_SCALAR_VISIBILITY: if ( _polyDataGroup->getCurrentPolyDataMapper()->GetScalarVisibility() ) msg = FXSEL ( SEL_COMMAND, FXWindow::ID_CHECK ); break;
    case ID_TOGGLE_GLOBAL_IMMEDIATE_MODE_RENDERING: if ( _polyDataGroup->getCurrentPolyDataMapper()->GetGlobalImmediateModeRendering() ) msg = FXSEL ( SEL_COMMAND, FXWindow::ID_CHECK ); break;
    default:
      throw aol::UnimplementedCodeException( "Selected FXSelector value not handled", __FILE__, __LINE__);
      break;
    }
    sender->handle ( this, msg, NULL );
    sender->handle ( this, FXSEL ( SEL_COMMAND, FXWindow::ID_ENABLE ), NULL );
  } else
    sender->handle ( this, FXSEL ( SEL_COMMAND, FXWindow::ID_DISABLE ), NULL );
  return 1;
}

void vtkFXViewer::loadVolume ( const char * FileName, vtkImageData * & target, const bool LoadAsFloat ) {
  cerr << "reading from " << FileName << endl;
  qc::ScalarArray<FXfloat, qc::QC_3D> volumeArray ( FileName );
  //  Upscale the input volume to enhance display quality
  //    qc::ScalarArray<FXfloat, qc::QC_3D> tempVolumeArray( FileName );
  //    qc::ScalarArray<FXfloat, qc::QC_3D> volumeArray( 512, 512, 512 );
  //    volumeArray.scale( tempVolumeArray );
  if ( !LoadAsFloat )
    volumeArray /= volumeArray.getMaxValue();
  if ( !LoadAsFloat )
    convertQCArrayToVtkImageData<FXfloat> ( volumeArray, target );
  else
    convertQCArrayToVtkImageDataFloat<FXfloat> ( volumeArray, target );
}

void vtkFXViewer::loadAndShowVolume ( const char * FileName ) {
  if ( aol::fileNameEndsWith ( FileName, ".par" ) ) {
    if ( !_volumeMapper->IsA ( "vtkOpenGLVolumeTextureMapper3D" ) )
      this->handle ( this, FXSEL ( SEL_COMMAND, ID_MAPPER_TEXTURE_3D ), NULL );

    aol::ParameterParser parser ( FileName );

    vtkImageData *imageColorR = vtkImageData::New();
    vtkImageData *imageColorG = vtkImageData::New();
    vtkImageData *imageColorB = vtkImageData::New();
    vtkImageData *imageOpacity = vtkImageData::New();

    loadVolume ( parser.getString( "volumeColorR" ).c_str(), imageColorR );
    loadVolume ( parser.getString( "volumeColorG" ).c_str(), imageColorG );
    loadVolume ( parser.getString( "volumeColorB" ).c_str(), imageColorB );
    loadVolume ( parser.getString( "volumeOpacity" ).c_str(), imageOpacity );

    vtkImageAppendComponents* imageAppendComponents = vtkImageAppendComponents::New();
    imageAppendComponents->AddInput ( imageColorR );
    imageAppendComponents->AddInput ( imageColorG );
    imageAppendComponents->AddInput ( imageColorB );
    imageAppendComponents->AddInput ( imageOpacity );

    imageColorR->Delete();
    imageColorG->Delete();
    imageColorB->Delete();
    imageOpacity->Delete();

    _volumeMapper->SetInputConnection ( imageAppendComponents->GetOutputPort() );
    _volume->GetProperty ()->IndependentComponentsOff ();
    imageAppendComponents->Delete();
  }
  else {
    vtkImageData *volumeData = vtkImageData::New();
    loadVolume ( FileName, volumeData );
    _volumeMapper->SetInput ( volumeData );
    volumeData->Delete();
  }
  _volume->SetMapper ( _volumeMapper );
  if ( !(_renderer->HasViewProp ( _volume )) )
    _renderer->AddViewProp ( _volume );
}

void vtkFXViewer::loadAndShowPLYStriped (const char *FileName, const bool LeftView, const float StripeWidth ) {
  aol::TriangMesh<FXfloat> surfMesh;
  surfMesh.loadFromUDPLY ( FileName );

  vtkPolyData* meshPolyData;
  FXfloat minC = aol::NumberTrait<FXfloat>::NaN, maxC = aol::NumberTrait<FXfloat>::NaN;

  meshPolyData = aol::TriangMeshToVTK<FXfloat>::makePolyData ( surfMesh, minC, maxC );

  vtkAppendPolyData *polyDataCollection = vtkAppendPolyData::New();
  buildStripedPolyData ( meshPolyData, polyDataCollection, LeftView, StripeWidth );
  meshPolyData->Delete();
  _polyDataGroup->setColorMapRange ( minC, maxC );
  FXString label = FileName;
  label += LeftView ? " left" : " right" ;
  addMeshPolyDataToRenderer ( polyDataCollection->GetOutput(), label );
  polyDataCollection->Delete();
}

long vtkFXViewer::onLoadVolumeClick ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  FXString outFile = _fileDialog.getOpenFilename ( this, "Open volume file", _fileDialog.getDirectory() + "/", "Quoc 3d volume file (*)" );
  _fileDialog.setDirectory ( getDirectoryFromString ( outFile ) );

  if ( outFile != FXString::null ) {
    loadAndShowVolume ( outFile.text() );
    _renderer->ResetCamera();
  }
  return 1;
}

long vtkFXViewer::loadPlyClick ( const int Format ) {
  FXString suffixField = "ply file (*.ply,*.ply.bz2";
  if ( Format == aol::UD_PLY )
    suffixField += ",*.udply";
  suffixField += ")";
  FXString* outFiles = _fileDialog.getOpenFilenames ( this, "Open ply files", _fileDialog.getDirectory() + "/", suffixField );

  if ( outFiles != NULL ) {
    _fileDialog.setDirectory ( getDirectoryFromString ( outFiles[0] ) );
    int i = 0;
    do {
      try {
        PLYToPolyData<FXfloat> plyToPolyData ( outFiles[i].text(), static_cast<aol::PLY_FORMAT> ( Format ) );
        _polyDataGroup->setColorMapRange ( plyToPolyData.getMinValue(), plyToPolyData.getMaxValue() );
        addMeshPolyDataToRenderer ( plyToPolyData.getPolyDataPointer(), outFiles[i] );
      }
      catch ( aol::Exception &el ) {
        el.dump();
        FXMessageBox::question ( this, MBOX_OK, "Warning", "%s", el.getMessage().c_str() );
        i++;
        continue;
      }
      i++;
    } while ( outFiles[i] != "" );
    _backgroundColor = FXRGB ( 255, 255, 255 );
    _backgroundColorTarget.handle ( this, FXSEL ( SEL_COMMAND, FXDataTarget::ID_VALUE ), NULL );
  }
  return 1;
}

long vtkFXViewer::onLoadUDPlyClick ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  return ( loadPlyClick ( aol::UD_PLY ) );
}

long vtkFXViewer::onLoadStanfordPlyClick ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  return ( loadPlyClick ( aol::STANFORD_PLY ) );
}

long vtkFXViewer::onLoadSRFClick ( FXObject* /*obj*/, FXSelector /*sel*/, void* /*data*/ ) {
  FXString suffixField = "srf file (*.srf)";
  FXString* outFiles = _fileDialog.getOpenFilenames ( this, "Open srf files", _fileDialog.getDirectory() + "/", suffixField );

  if ( outFiles != NULL ) {
    _fileDialog.setDirectory ( getDirectoryFromString ( outFiles[0] ) );
    int i = 0;
    do {
      aol::TriangMesh<FXfloat> surfMesh;

      try {
        surfMesh.loadFromSRF ( outFiles[i].text() );
      }
      catch ( aol::Exception &el ) {
        el.dump();
        FXMessageBox::question ( this, MBOX_OK, "Warning", "%s", el.getMessage().c_str() );
        i++;
        continue;
      }

      vtkPolyData* meshPolyData;
      FXfloat minC = aol::NumberTrait<FXfloat>::NaN, maxC = aol::NumberTrait<FXfloat>::NaN;

      meshPolyData = aol::TriangMeshToVTK<FXfloat>::makePolyData ( surfMesh, minC, maxC );

      _polyDataGroup->setColorMapRange ( minC, maxC );
      addMeshPolyDataToRenderer ( meshPolyData, outFiles[i] );
      meshPolyData->Delete();
      i++;
    } while ( outFiles[i] != "" );
    _backgroundColor = FXRGB ( 255, 255, 255 );
    _backgroundColorTarget.handle ( this, FXSEL ( SEL_COMMAND, FXDataTarget::ID_VALUE ), NULL );
  }
  return 1;
}

long vtkFXViewer::onLoadUDPlyStripedClick ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  FXString outFile = _fileDialog.getOpenFilename ( this, "Open ply-stripe-par file", _fileDialog.getDirectory() + "/", "Strip-par file (*.par)" );
  _fileDialog.setDirectory ( getDirectoryFromString ( outFile ) );

  if ( outFile != FXString::null ) {
    aol::ParameterParser parser ( outFile.text() );

    _fileDialog.setDirectory ( getDirectoryFromString ( outFile ) );

    if ( ( parser.hasVariable ( "setCurDirToParFileDir" ) == true ) && ( parser.getInt( "setCurDirToParFileDir" ) == 1 ) )
      aol::setCurrentDirectory ( getDirectoryFromString ( outFile ).text() );

    const float stripeWidth = parser.getReal<float> ( "stripewidth" );
    if ( parser.hasVariable ( "filename-left" ) )
      loadAndShowPLYStriped ( parser.getString ( "filename-left" ).c_str(), true, stripeWidth );
    if ( parser.hasVariable ( "filename-right" ) )
      loadAndShowPLYStriped ( parser.getString ( "filename-right" ).c_str(), false, stripeWidth );

    _backgroundColor = FXRGB ( 255, 255, 255 );
    _backgroundColorTarget.handle ( this, FXSEL ( SEL_COMMAND, FXDataTarget::ID_VALUE ), NULL );
  }
  return 1;
}

long vtkFXViewer::onLoadVTKPolyDataClick ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  FXString* outFiles = _fileDialog.getOpenFilenames ( this, "Open VTK polydata files", _fileDialog.getDirectory() + "/", "VTK polydata file (*.vtk)" );

  if ( outFiles != NULL ) {
    _fileDialog.setDirectory ( getDirectoryFromString ( outFiles[0] ) );
    int i = 0;
    do {
      vtkPolyDataReader *polyDataReader = vtkPolyDataReader::New();
      polyDataReader->SetFileName ( outFiles[i].text() );
      addMeshPolyDataToRenderer ( polyDataReader->GetOutput(), outFiles[i] );
      polyDataReader->Delete();
      i++;
    } while ( outFiles[i] != "" );
    _backgroundColor = FXRGB ( 255, 255, 255 );
    _backgroundColorTarget.handle ( this, FXSEL ( SEL_COMMAND, FXDataTarget::ID_VALUE ), NULL );
  }
  return 1;
}

long vtkFXViewer::onLoadTPCFEParClick ( FXObject*, FXSelector, void* ) {
  FXString outFile = _fileDialog.getOpenFilename ( this, "Open tpcfe-par file", _fileDialog.getDirectory() + "/", "TPCFE-par file (*.par)" );
  _fileDialog.setDirectory ( getDirectoryFromString ( outFile ) );

  if ( outFile != FXString::null ) {

    char parfilename[1024];
    sprintf ( parfilename, "%s", outFile.text() );

    tpcfe::TpcfeToVTK< tpcfe::CFEGrid< double, tpcfe::CFE_CD > > tpcfeToVTKConverter ( parfilename, _polyDataGroup );

    vtkPolyData *meshPolyData = NULL;
    try {
      meshPolyData = tpcfeToVTKConverter.makePolyData( );
    }
    catch ( aol::Exception &el ) {
      el.dump();
      FXMessageBox::question ( this, MBOX_OK, "Warning", "%s", el.getMessage().c_str() );
      return 1;
    }

    if ( meshPolyData != NULL ) {
      addMeshPolyDataToRenderer ( meshPolyData, outFile );
      meshPolyData->Delete();

      //      _polyDataGroup->getCurrentPolyDataMapper()->SetLookupTable( colorTF );
      _polyDataGroup->getCurrentPolyDataMapper()->SetColorModeToMapScalars();
    }
  }
  return 1;
}

long vtkFXViewer::onLoadLevelsetTpcfeClick ( FXObject*, FXSelector, void* ) {
  FXString* outFiles = _fileDialog.getOpenFilenames ( this, "Open levelset file", _fileDialog.getDirectory() + "/", "Quoc 3D data (*)" );

  if ( outFiles != NULL ) {
    _fileDialog.setDirectory ( getDirectoryFromString ( outFiles[0] ) );
    int i = 0;
    do {
      qc::ScalarArray<FXfloat, qc::QC_3D> levelset ( outFiles[i].text() );
      tpcfe::LevelsetToVTK<FXfloat> converter ( levelset );
      vtkPolyData *meshPolyData = converter.makePolyData( );
      addMeshPolyDataToRenderer ( meshPolyData, outFiles[i] );
      meshPolyData->Delete();
      i++;
    } while ( outFiles[i] != "" );

  }
  return 1;
}

long vtkFXViewer::onLoadLevelsetClick ( FXObject*, FXSelector, void* ) {
  FXString outFile = _fileDialog.getOpenFilename ( this, "Open levelset file", _fileDialog.getDirectory() + "/", "Quoc 3D Array (*)" );
  _fileDialog.setDirectory ( getDirectoryFromString ( outFile ) );

  if ( outFile != FXString::null ) {
    qc::ScalarArray<FXfloat, qc::QC_3D> levelset ( outFile.text() );

    aol::TriangMesh<FXfloat> surfMesh;
    cerr << "Converting zero levelset to TriangMesh ... ";
    qcsm::LevelsetToTriangMesh<FXfloat> levelsetToTriangMesh;
    levelsetToTriangMesh ( levelset, surfMesh );
    cerr << "done\n";

    FXfloat dummy1, dummy2;
    vtkPolyData *meshPolyData = aol::TriangMeshToVTK<FXfloat>::makePolyData ( surfMesh, dummy1, dummy2 );

    addMeshPolyDataToRenderer ( meshPolyData, outFile );
    meshPolyData->Delete();
  }
  return 1;
}

/// \author von Deylen
long vtkFXViewer::onLoadLevelsetTexturedClick ( FXObject*, FXSelector, void* ) {
  FXString parameterFilename = _fileDialog.getOpenFilename ( this, "Open parameter file", _fileDialog.getDirectory() + "/", "Quoc 3D Array (*)" );
  std::string paramFileDir = getDirectoryFromString ( parameterFilename ).text();
  _fileDialog.setDirectory ( paramFileDir.c_str() );

  // *** Parse parameter file ***
  if (parameterFilename == FXString::null)
    return 0;

  paramFileDir += "/";

  aol::ParameterParser parser( parameterFilename.text() );
  std::string levelsetFilename = paramFileDir + parser.getString("levelsetFilename");
  std::string textureFilename  = paramFileDir + parser.getString("textureFilename");
  double levelsetDisplayValue  = parser.getDouble("levelsetDisplayValue");

  // *** Load level and texture data ***
  qc::ScalarArray<FXfloat, qc::QC_3D> levelset ( levelsetFilename.c_str() );
  qc::ScalarArray<FXfloat, qc::QC_3D> texture  ( textureFilename .c_str() );

  // substract offset (given in parameter file) from given level values
  levelset.addToAll(-levelsetDisplayValue);

  // *** Extract TriangMesh from ScalarArray ***
  aol::TriangMesh<FXfloat> surface;
  qcsm::LevelsetToTriangMesh<FXfloat> converter;

  converter.apply(levelset, surface);

  surface.createVertexData();

  // *** Add texture data ***
  for ( int i = 0; i < surface.getNumVertices(); ++i ) {
    surface.getVertexData ()[i] = texture.interpolate_on01 ( surface.getVertex ( i ) );
  }

  // *** Convert TriangMesh into vtk surface ***
  vtkPolyData* meshPolyData;

  FXfloat minC = aol::NumberTrait<FXfloat>::NaN, maxC = aol::NumberTrait<FXfloat>::NaN;
  meshPolyData = aol::TriangMeshToVTK<FXfloat>::makePolyData ( surface, minC, maxC );

  _polyDataGroup->setColorMapRange ( minC, maxC );
  addMeshPolyDataToRenderer ( meshPolyData, parameterFilename );
  meshPolyData->Delete();
  return 1;
}

/// \author von Deylen
long vtkFXViewer::onLoadLevelsetWithTextureCoordsClick ( FXObject*, FXSelector, void* ) {
  FXString parameterFilename = _fileDialog.getOpenFilename ( this, "Open parameter file", _fileDialog.getDirectory() + "/", "*.par" );
  std::string paramFileDir = getDirectoryFromString ( parameterFilename ).text();
  _fileDialog.setDirectory ( paramFileDir.c_str() );

  // *** Parse parameter file ***
  if (parameterFilename == FXString::null)
    return 0;

  paramFileDir += "\\";

  aol::ParameterParser parser( parameterFilename.text() );
  string levelsetFilename = paramFileDir + parser.getString("levelsetFilename");
  string tCoordXFilename  = paramFileDir + parser.getString("textureCoordsXFilename");
  string tCoordYFilename  = paramFileDir + parser.getString("textureCoordsYFilename");
  string tBitmapFilename  = paramFileDir + parser.getString("textureBitmapFilename");
  double levelsetDisplayValue  = parser.getDouble("levelsetDisplayValue");

  _polyDataGroup->pushBackAndShow(new IsoSurfaceDataWithTextureMap ( _polyDataGroup, levelsetFilename.c_str(),
                                                                     levelsetFilename, tCoordXFilename, tCoordYFilename,
                                                                     tBitmapFilename, levelsetDisplayValue ) );
  _renderer->ResetCamera();
  return 1;
}

long vtkFXViewer::onLoadGraphClick ( FXObject*, FXSelector, void* ) {
  FXString outFile = _fileDialog.getOpenFilename ( this, "Open 2d quoc array", _fileDialog.getDirectory() + "/", "2d quoc array (*.pgm,*.bz2,*.dat,*.q2bz)" );
  _fileDialog.setDirectory ( getDirectoryFromString ( outFile ) );

  if ( outFile != FXString::null ) {
    const bool loadGraphTexture = ( MBOX_CLICKED_YES == FXMessageBox::question ( this, MBOX_YES_NO, "Load Graph Texture", "Load a texture for the graph?" ) );
    FXString textureFile = loadGraphTexture ? _fileDialog.getOpenFilename ( this, "Open 2d quoc array", _fileDialog.getDirectory() + "/", "2d quoc array (*.pgm,*.bz2,*.dat,*.q2bz,*.png)" ) : FXString::null;

    GraphToPolyData<FXfloat> graphToPolyData ( outFile.text(), ( textureFile != FXString::null ) ? textureFile.text() : NULL );

    _polyDataGroup->addPolyData ( graphToPolyData.getPolyDataPointer(), outFile );
    if ( graphToPolyData.hasTexture() == false ) {
      if ( graphToPolyData.hasScalarData() == false )
        _polyDataGroup->getCurrentPolyDataMapper()->SetScalarRange(graphToPolyData.getImageArrayReference().getMinValue(), graphToPolyData.getImageArrayReference().getMaxValue());
      else {
        _polyDataGroup->setColorMap ( 4 );
        _polyDataGroup->setColorMapRange ( graphToPolyData.getImageDataArrayReference().getMinValue(), graphToPolyData.getImageDataArrayReference().getMaxValue() );
      }
    }
    else
    {
      _polyDataGroup->getCurrentActor()->SetTexture ( graphToPolyData.getTexturePointer() );
    }
  }
  return 1;
}

long vtkFXViewer::onLoadLidarPointCloudClick ( FXObject*, FXSelector, void* ) {
  FXString outFile = _fileDialog.getOpenFilename ( this, "Open LIDAR file", _fileDialog.getDirectory() + "/", "LIDAR file (*.xyz)" );
  _fileDialog.setDirectory ( getDirectoryFromString ( outFile ) );

  if ( outFile != FXString::null ) {
    aol::Vector<double> verticeComponents;
    verticeComponents.read( outFile.text() );
    const int numPoints = verticeComponents.size() / 4;
    aol::MultiVector<double> verticesPlusIntensity( 4, numPoints );
    verticesPlusIntensity.copySplitTransposeFrom( verticeComponents );
    const aol::Vector<double> &intensitiesVec = verticesPlusIntensity[2];

    vtkPolyData *cloud = vtkPolyData::New();
    vtkPoints *points = vtkPoints::New();
    vtkFloatArray *intensities = vtkFloatArray::New();
    intensities ->SetNumberOfTuples(numPoints );

    for ( int i = 0; i < numPoints; ++i ) {
      points->InsertPoint(i,verticesPlusIntensity[0][i], verticesPlusIntensity[1][i], verticesPlusIntensity[2][i]);
      intensities->SetValue(i, intensitiesVec[i]);
    }

    cloud->SetPoints(points);
    points->Delete();
    cloud->GetPointData()->SetScalars(intensities);
    intensities->Delete();

    vtkVertexGlyphFilter* vertexGlyphFilter = vtkVertexGlyphFilter::New();
    vertexGlyphFilter->AddInput(cloud);
    vertexGlyphFilter->Update();
    cloud->Delete();

    _polyDataGroup->addPolyData ( vertexGlyphFilter->GetOutput(), outFile.text() );
    _polyDataGroup->setColorMapRange ( intensitiesVec.getMinValue(), intensitiesVec.getMaxValue() );
    vertexGlyphFilter->Delete();
  }
  return 1;
}

long vtkFXViewer::onLoadOverlayImage2DClick ( FXObject*, FXSelector, void* ) {
  FXString outFile = _fileDialog.getOpenFilename ( this, "Open PNG file", _fileDialog.getDirectory() + "/", "PNG file (*.png)" );
  _fileDialog.setDirectory ( getDirectoryFromString ( outFile ) );

  if ( outFile != FXString::null ) {
    vtkPNGReader *pngReader = vtkPNGReader::New();
    pngReader->SetFileName ( outFile.text() );
    vtkTexture *texture = vtkTexture::New();
    texture->SetInputConnection ( pngReader->GetOutputPort() );

    _overlayImageSize = aol::Min ( _renWindow->GetSize()[0], _renWindow->GetSize()[1] ) / 2;
    _overlayImagePoints->Reset();
    _overlayImagePoints->InsertPoint(0, 0., 0., 0.);
    _overlayImagePoints->InsertPoint(1, 1., 0., 0.);
    _overlayImagePoints->InsertPoint(2, 1., 1., 0.);
    _overlayImagePoints->InsertPoint(3, 0., 1., 0.);

    vtkCellArray *cells = vtkCellArray::New();
    cells->InsertNextCell(4);
    cells->InsertCellPoint(0);
    cells->InsertCellPoint(1);
    cells->InsertCellPoint(2);
    cells->InsertCellPoint(3);

    vtkFloatArray *tcoords = vtkFloatArray::New();
    tcoords->SetNumberOfComponents(2);
    tcoords->InsertNextTuple2(0.f, 0.f);
    tcoords->InsertNextTuple2(1.f, 0.f);
    tcoords->InsertNextTuple2(1.f, 1.f);
    tcoords->InsertNextTuple2(0.f, 1.f);

    vtkPolyData *textureCoords = vtkPolyData::New();
    textureCoords->SetPoints(_overlayImagePoints);
    textureCoords->SetPolys(cells);
    textureCoords->GetPointData()->SetTCoords(tcoords);

    _overlayImageMapper->SetInput( textureCoords );

    vtkCoordinate *coordinate = vtkCoordinate::New();
    coordinate->SetCoordinateSystemToViewport();

    _overlayImageMapper->SetTransformCoordinate ( coordinate );

    _overlayImageActor->SetTexture ( texture );
    _overlayImageActor->SetMapper ( _overlayImageMapper );
    _overlayImageActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    _renderer->AddActor ( _overlayImageActor );

    cells->Delete();
    tcoords->Delete();
    textureCoords->Delete();
    texture->Delete();
    pngReader->Delete();
    coordinate->Delete();

    onCmdUpdateOverlayImage ( NULL, 0, NULL );
  }
  return 1;
}

long vtkFXViewer::onCmdUpdateOverlayImage ( FXObject*, FXSelector, void* ) {
  _overlayImageActor->SetPosition ( _overlayImagePos[0], _overlayImagePos[1] );
  _overlayImageActor->GetProperty()->SetOpacity ( _overlayImageOpacity );
  _overlayImagePoints->SetPoint ( 1, _overlayImageSize, 0., 0. );
  _overlayImagePoints->SetPoint ( 2, _overlayImageSize, _overlayImageSize, 0. );
  _overlayImagePoints->SetPoint ( 3, 0., _overlayImageSize, 0. );
  updateVTKCanvas();
  return 1;
}

// Make sure the rendering window is set to offscreen rendering, otherwise the window will be drawn exactly as seen on the screen (e.g. partially hidden by other windows)
void vtkFXViewer::exportCurrentViewToPng ( const char *FileName ) {
  vtkWindowToImageFilter *windowToImageFilter = vtkWindowToImageFilter::New();
  windowToImageFilter->SetInputBufferTypeToRGBA();
  windowToImageFilter->SetInput ( _renWindow );
  vtkPNGWriter *pngWriter = vtkPNGWriter::New();
  pngWriter->SetInputConnection ( windowToImageFilter->GetOutputPort() );
  pngWriter->SetFileName ( FileName );
  pngWriter->Write();
  windowToImageFilter->Delete();
  pngWriter->Delete();
}


void vtkFXViewer::exportCurrentViewToHRPng ( const char *FileName ) {
  vtkRenderLargeImage *rendererLarge = vtkRenderLargeImage::New();

#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  rendererLarge->SetInput ( _renderer->GetLeftRenderer() );
#else
  rendererLarge->SetInput ( _renderer );
#endif
  rendererLarge->SetMagnification ( 4 );

  vtkPNGWriter *pngWriter = vtkPNGWriter::New();
  pngWriter->SetInputConnection ( rendererLarge->GetOutputPort() );
  pngWriter->SetFileName ( FileName );
  pngWriter->Write();

  rendererLarge->Delete();
  pngWriter->Delete();
}


long vtkFXViewer::onCmdIsoline ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  FXString outFile = _fileDialog.getOpenFilename ( this, "Open 2d quoc array", _fileDialog.getDirectory() + "/", "2d quoc array (*.pgm,*.bz2,*.dat,*.q2bz)" );

  if ( outFile != FXString::null ) {
    qc::ScalarArray<FXfloat, qc::QC_2D> imageArray ( outFile.text() );
    const FXfloat arrayMinValue = imageArray.getMinValue();
    const FXfloat arrayMaxValue = imageArray.getMaxValue();
    vtkImageData *imageData = vtkImageData::New();

    convertQCArrayToVtkImageDataFloat<FXfloat> ( imageArray, imageData );

    vtkMarchingSquares* isoXY = vtkMarchingSquares::New();
    isoXY->SetInput ( imageData );
    isoXY->GenerateValues ( 200, arrayMinValue, arrayMaxValue );

    // Transform to [0,1]^2
    vtkTransform *transform = vtkTransform::New();
    transform->Scale(1./imageArray.getNumX(),1./imageArray.getNumY(),1.0);

    vtkTransformPolyDataFilter *transF = vtkTransformPolyDataFilter::New();
    transF->SetInputConnection(isoXY->GetOutputPort());
    transF->SetTransform(transform);
    transF->Update();

    _polyDataGroup->addPolyData ( transF->GetOutput(), outFile );
    _polyDataGroup->setColorMapRange ( arrayMinValue, arrayMaxValue );

    imageData->Delete();
    isoXY->Delete();
    transform->Delete();
    transF->Delete();

    vtkScalarBarActor *scalarBarActor = vtkScalarBarActor::New();
    scalarBarActor->SetLookupTable ( _polyDataGroup->getCurrentPolyDataMapper()->GetLookupTable() );
    _renderer->AddViewProp ( scalarBarActor );
  }
  return 1;
}

long vtkFXViewer::onCmdIsosurface ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  FXString* outFiles = _fileDialog.getOpenFilenames ( this, "Open volume files", _fileDialog.getDirectory() + "/", "Quoc 3d volume file (*)" );

  // load
  if ( outFiles != NULL ) {
    _fileDialog.setDirectory ( getDirectoryFromString ( outFiles[0] ) );

    // get level value from user
    FXdouble levelValue = 0.5;
    FXInputDialog::getReal ( levelValue, this, "Select level value", "Which level surface shall be shown?" );

    int i = 0;
    do {
      // IsoSurfaceData takes ownership of the vtkImageData pointer, therefore we must not delete it here.
      vtkImageData * volumeData = vtkImageData::New();
      loadVolume ( outFiles[i].text(), volumeData, true );
      _polyDataGroup->pushBackAndShow(new IsoSurfaceData ( _polyDataGroup, outFiles[i], volumeData, levelValue ));
    } while ( outFiles[++i] != "" );
    // show
    _renderer->ResetCamera();
  }
  return 1;
}

long vtkFXViewer::onExportPngClick ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  FXString           filename;

  _saveDialog.setFilename ( _saveDialog.getDirectory() + "/" );
  if ( _saveDialog.execute() ) {
    if ( FXStat::exists ( _saveDialog.getFilename() ) ) {
      if ( MBOX_CLICKED_NO == FXMessageBox::question ( this, MBOX_YES_NO, "Overwrite Image", "Overwrite existing image?" ) ) return 1;
    }
    filename = _saveDialog.getFilename();
    _saveDialog.setDirectory ( getDirectoryFromString ( filename ) );

//       _renWindow->SetOffScreenRendering(1);
    exportCurrentViewToPng ( filename.text() );
//       _renWindow->SetOffScreenRendering(0);
  }
  return 1;
}

long vtkFXViewer::onExportHRPngClick ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  FXString           filename;

  _saveDialog.setFilename ( _saveDialog.getDirectory() + "/" );
  if ( _saveDialog.execute() ) {
    if ( FXStat::exists ( _saveDialog.getFilename() ) ) {
      if ( MBOX_CLICKED_NO == FXMessageBox::question ( this, MBOX_YES_NO, "Overwrite Image", "Overwrite existing image?" ) ) return 1;
    }
    filename = _saveDialog.getFilename();
    _saveDialog.setDirectory ( getDirectoryFromString ( filename ) );

    exportCurrentViewToHRPng ( filename.text() );
  }
  return 1;
}

long vtkFXViewer::onExportPdfClick ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  // GL2PS doesn't seem to like vtkStripper
  if ( _polyDataGroup->getCurrentPolyDataMapper() )
    _polyDataGroup->getCurrentPolyDataMapper()->SetInputConnection ( _polyDataGroup->getCurrentPolyDataNormals()->GetOutputPort() );

  vtkGL2PSExporter *gl2PSExporter = vtkGL2PSExporter::New();
  gl2PSExporter->DrawBackgroundOff();
  gl2PSExporter->SetFileFormatToPDF();
  gl2PSExporter->SetRenderWindow ( _renWindow );
  gl2PSExporter->SetFilePrefix ( "1" );
  gl2PSExporter->SetSortToBSP();
  //gl2PSExporter->Write3DPropsAsRasterImageOn();
  gl2PSExporter->Write();
  gl2PSExporter->Delete();

  if ( _polyDataGroup->getCurrentPolyDataMapper() )
    _polyDataGroup->getCurrentPolyDataMapper()->SetInputConnection ( _polyDataGroup->getCurrentStripper()->GetOutputPort() );
  return 1;
}

long vtkFXViewer::onCmdPlane1Normal ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  _plane1->SetOrigin ( _plane1Origin );
  FXdouble temp = sqrt ( pow ( _plane1Normal[0], 2 ) + pow ( _plane1Normal[1], 2 ) + pow ( _plane1Normal[2], 2 ) );
  for ( int i = 0; i < 3; i++ ) {
    _plane1Normal[i] /= temp;
  }
  _plane1->SetNormal ( _plane1Normal );
  updateVTKCanvas();
  return 1;
}

long vtkFXViewer::onCmdTogglePlane1 ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  _polyDataGroup->ClippingPlaneToAll ( _plane1, !!_plane1Active );

  if ( _plane1Active ) {
    _volumeMapper->AddClippingPlane ( _plane1 );
  } else {
    _volumeMapper->RemoveClippingPlane ( _plane1 );
  }

  //_volume->SetMapper( _volumeMapper );
  updateVTKCanvas();

  return 1;
}

long vtkFXViewer::onCmdMakeVideo ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  FXString* outFiles = _fileDialog.getOpenFilenames ( this, "Open volume files", _fileDialog.getDirectory() + "/", "Quoc 3d volume file (*)" );

  if ( outFiles != NULL ) {
    _renWindow->SetOffScreenRendering ( 1 );
    char fileName[1024];
    int numberOfFiles = 0;
    do {
      numberOfFiles++;
    } while ( outFiles[numberOfFiles] != "" );

    FXProgressDialog *progressdialog = new FXProgressDialog ( this, "Making video", "Rendering the frames of selected files\nusing off screen rendering.", DECOR_BORDER | DECOR_RESIZE | DECOR_TITLE | LAYOUT_EXPLICIT );
    progressdialog->setBarStyle ( PROGRESSBAR_PERCENTAGE | PROGRESSBAR_NORMAL );
    progressdialog->setTotal ( numberOfFiles );
    progressdialog->setProgress ( 0 );
    progressdialog->create();
    progressdialog->show ( PLACEMENT_SCREEN );

    int i = 0;
    do {
      /*
              this->raise();
              this->forceRefresh();
              this->update();
              this->repaint();
              progressdialog->layout();
              progressdialog->update();
      */
      progressdialog->repaint();
      loadAndShowVolume ( outFiles[i].text() );
      sprintf ( fileName, "%03d.png", i );
      exportCurrentViewToPng ( fileName );
      i++;
      progressdialog->setProgress ( i );
    } while ( outFiles[i] != "" );
    _renWindow->SetOffScreenRendering ( 0 );
    delete progressdialog;
  }
  return 1;
}


void vtkFXViewer::create_volume_pipeline( ) {

  _renWindow = vtkRenderWindow::New();
#ifndef SIDE_BY_SIDE_STEREO_RENDERING
  _renWindow->AddRenderer ( _renderer );
#else
  _renWindow -> StereoCapableWindowOn();
  _renWindow -> StereoRenderOn();
  _renWindow -> SetStereoTypeToLeft();
  _renWindow->AddRenderer ( _renderer->GetLeftRenderer() );

  _renWindowRight = vtkRenderWindow::New();
  _renWindowRight -> StereoCapableWindowOn();
  _renWindowRight -> StereoRenderOn();
  _renWindowRight -> SetStereoTypeToRight();
  _renWindowRight ->AddRenderer ( _renderer->GetRightRenderer() );
#endif

  //vtkInteractorStyle *style = vtkInteractorStyleTrackball::New();
  vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
  //vtkInteractorStyle *style = vtkInteractorStyleUnicam::New();

  vtkFXRenderWindowInteractor *rwi = _fxVTKCanvas->getInteractor();
  rwi->SetInteractorStyle ( style );

  rwi->SetRenderWindow ( _renWindow );
  rwi->Initialize();

#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  vtkFXRenderWindowInteractor *rwiRight = _fxVTKCanvasRight->getInteractor();
  vtkInteractorStyleTrackballCamera *styleRight = vtkInteractorStyleTrackballCamera::New();
  rwiRight->SetInteractorStyle ( styleRight );

  rwiRight->SetRenderWindow ( _renWindowRight );
  rwiRight->Initialize();
#endif

  _opacityFunction = vtkPiecewiseFunction::New();
  _opacityFunction->AddPoint ( 0, _opacityFunctionLowerValue );
  _opacityFunction->AddPoint ( 65535, _opacityFunctionUpperValue );
  _colorFunction = vtkPiecewiseFunction::New();
  _colorFunction->AddPoint ( 0, _colorFunctionLowerValue );
  _colorFunction->AddPoint ( 65535, _colorFunctionUpperValue );
  //vtkColorTransferFunction *negColor = vtkColorTransferFunction::New();
  //negColor->AddRGBPoint( 64, 1.0, 0.0, 0.0 );
  //negColor->AddRGBPoint( 128, 0.0, 0.0, 1.0 );
  //negColor->AddRGBPoint( 196, 0.0, 1.0, 0.0 );

  vtkVolumeProperty *negProperty = vtkVolumeProperty::New();
  negProperty->SetColor ( _colorFunction );
  negProperty->SetScalarOpacity ( _opacityFunction );
  negProperty->SetInterpolationTypeToLinear();
  //negProperty->ShadeOn();
  _volume->SetProperty ( negProperty );
  //_volume->SetScale( 1., 1., 1.5 );

  style->Delete();
  negProperty->Delete();

}

void vtkFXViewer::addMeshPolyDataToRenderer ( vtkPolyData* PolyData, const FXString & DataName ) {
  /*
      // Generate normals for the mesh
      _meshNormals->SetInput(_meshPolyData);
      _meshNormals->SetFeatureAngle(60.0);
      // Generate triangle strips
      _meshStripper->SetInputConnection(_meshNormals->GetOutputPort());

      _polyDataGroup->getCurrentPolyDataMapper()->SetInputConnection(_meshStripper->GetOutputPort());
  #ifdef SIDE_BY_SIDE_STEREO_RENDERING
      _polyDataGroup->getCurrentPolyDataMapper()->GlobalImmediateModeRenderingOn();
  #endif
      //_polyDataGroup->getCurrentPolyDataMapper()->SetInput(meshPolyData);
      _meshActor->SetMapper(_meshMapper);
      _meshActor->SetProperty( _meshProperty );

      if( !_renderer->HasViewProp( _meshActor ) )
        _renderer->AddViewProp(_meshActor);
  */
  _polyDataGroup->addPolyData ( PolyData, DataName );

  _renderer->ResetCamera();

}




