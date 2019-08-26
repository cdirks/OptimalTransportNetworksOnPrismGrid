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

#include "vtkFXGui.h"

vtkFXGui::vtkFXGui ( FXApp *app ) :
#ifdef SIDE_BY_SIDE_STEREO_RENDERING
  FXMainWindow ( app, "vtkFOX Viewer", NULL, NULL, DECOR_ALL, 0, 0, 3200, 1200 ),
#else
  FXMainWindow ( app, "vtkFOX Viewer", NULL, NULL, DECOR_ALL, 0, 0, 800, 600 ),
#endif
  _filemenu ( NULL ),
  _viewmenu ( NULL ),
  _viewer   ( NULL ) {


  FXDockSite *topdock = new FXDockSite ( this, DOCKSITE_NO_WRAP | LAYOUT_SIDE_TOP | LAYOUT_FILL_X ); // for menus
  FXDockSite *rightdock = new FXDockSite ( this, LAYOUT_SIDE_RIGHT | LAYOUT_FILL_Y );                // for control panel

  _viewer = new vtkFXViewer( this, this, app );

  FXMenuBar* menubar = new FXMenuBar ( topdock, LAYOUT_DOCK_SAME | LAYOUT_SIDE_TOP | LAYOUT_FILL_X | FRAME_RAISED );

  _filemenu = new FXMenuPane ( this );
  new FXMenuTitle ( menubar, "File", NULL, _filemenu );
  new FXMenuCommand ( _filemenu, "&Load Volume\t\tLoad 3d volume file.",                              NULL, _viewer, vtkFXViewer::ID_LOAD_VOLUME );
  new FXMenuCommand ( _filemenu, "Load UD-ply\tCtrl-u\t\tLoad UD-ply file.",                          NULL, _viewer, vtkFXViewer::ID_LOAD_UDPLY );
  new FXMenuCommand ( _filemenu, "Load Stanford-ply\t\t\tLoad Stanford-ply file.",              NULL, _viewer, vtkFXViewer::ID_LOAD_STANFORDPLY );
  new FXMenuCommand ( _filemenu, "Load BrainVoyager SRF\t\t\tLoad BrainVoyager SRF file.",              NULL, _viewer, vtkFXViewer::ID_LOAD_SRF );
  new FXMenuCommand ( _filemenu, "Load UD-ply(s) striped\tCtrl-u\t\tLoad striped view parameter file.", NULL, _viewer, vtkFXViewer::ID_LOAD_UDPLYSTRIPED );
  new FXMenuCommand ( _filemenu, "Load VTK polydata\tCtrl-u\t\tLoad VTK polydata file.",              NULL, _viewer, vtkFXViewer::ID_LOAD_VTKPOLYDATA );
  new FXMenuCommand ( _filemenu, "Load TPCFE-par\t\tLoad TPCFE-par file.",                            NULL, _viewer, vtkFXViewer::ID_LOAD_TPCFE_PAR );
  new FXMenuCommand ( _filemenu, "Load levelset via tpcfe\t\tLoad levelset from 3D data via tpcfe.",  NULL, _viewer, vtkFXViewer::ID_LOAD_LEVELSET_TPCFE );
  new FXMenuCommand ( _filemenu, "Load Levelset\t\tLoad Levelset from 3D Quoc Array.",                NULL, _viewer, vtkFXViewer::ID_LOAD_LEVELSET );
  new FXMenuCommand ( _filemenu, "Load textured Levelset via parameter file\t\tLoad textured Levelset and texture from 3D Quoc Array.", NULL, _viewer, vtkFXViewer::ID_LOAD_LEVELSET_TEXTURED );
  new FXMenuCommand ( _filemenu, "Load Levelset and texture coords via parameter file\t\tLoad Levelset and texture coords from 3D Quoc Array.", NULL, _viewer, vtkFXViewer::ID_LOAD_LEVELSET_WITH_TEXTURECOORDS );
  new FXMenuCommand ( _filemenu, "Load Graph\t\tLoad Graph from 2D Quoc Array.",                      NULL, _viewer, vtkFXViewer::ID_LOAD_GRAPH );
  new FXMenuCommand ( _filemenu, "Load Isosurface\t\tLoad Isosurface from 3D Quoc Array with VTK marching cubes.", NULL, _viewer, vtkFXViewer::ID_ISOSURFACE );
  new FXMenuCommand ( _filemenu, "Load LIDAR Point Cloud\t\tLoad LIDAR Point Cloud."       ,                      NULL, _viewer, vtkFXViewer::ID_LOAD_LIDAR_POINT_CLOUD );
  new FXMenuCommand ( _filemenu, "Load overlay image 2D\t\tLoad overlay image 2D."       ,                      NULL, _viewer, vtkFXViewer::ID_LOAD_OVERLAY_IMAGE_2D );
  new FXMenuCommand ( _filemenu, "Export PNG\t\tSave GLCanvas content as PNG image.",                 NULL, _viewer, vtkFXViewer::ID_PNG_EXPORT_IMAGE );
  new FXMenuCommand ( _filemenu, "Export high-res PNG\t\tSave GLCanvas content as PNG image in high resolution.", NULL, _viewer, vtkFXViewer::ID_PNG_HR_EXPORT_IMAGE );
  new FXMenuCommand ( _filemenu, "Export PDF\t\tSave GLCanvas content as PDF.",                       NULL, _viewer, vtkFXViewer::ID_PDF_EXPORT );
  new FXMenuCommand ( _filemenu, "Exit\tAlt-F4,Ctrl-q", NULL, getApp(), FXApp::ID_QUIT );

  _viewmenu = new FXMenuPane ( this );
  new FXMenuTitle ( menubar, "View", NULL, _viewmenu );
  new FXMenuCommand ( _viewmenu, "2D Textures\t\tChange volume renderer to 2D texture mode",           NULL, _viewer, vtkFXViewer::ID_MAPPER_TEXTURE_2D );
  new FXMenuCommand ( _viewmenu, "3D Textures\t\tChange volume renderer to 3D texture mode",           NULL, _viewer, vtkFXViewer::ID_MAPPER_TEXTURE_3D );
  new FXMenuCommand ( _viewmenu, "Ray Cast comp\t\tChange volume renderer to ray cast composite mode", NULL, _viewer, vtkFXViewer::ID_MAPPER_RAY_CAST_COMPOSITE );
  new FXMenuCommand ( _viewmenu, "Ray Cast MIP\t\tChange volume renderer to ray cast MIP mode",        NULL, _viewer, vtkFXViewer::ID_MAPPER_RAY_CAST_MIP );
  new FXMenuSeparator ( _viewmenu );
  new FXMenuRadio ( _viewmenu, "Points\t\tRender as points.",         _viewer, vtkFXViewer::ID_STYLE_POINTS, MENU_AUTOGRAY );
  new FXMenuRadio ( _viewmenu, "Wire Frame\t\tRender as wire frame.", _viewer, vtkFXViewer::ID_STYLE_WIREFRAME, MENU_AUTOGRAY );
  new FXMenuRadio ( _viewmenu, "Surface\t\tRender solid surface.",    _viewer, vtkFXViewer::ID_STYLE_SURFACE, MENU_AUTOGRAY );
  new FXMenuRadio ( _viewmenu, "Flat\t\t",                            _viewer, vtkFXViewer::ID_STYLE_FLAT, MENU_AUTOGRAY );
  new FXMenuRadio ( _viewmenu, "Gouraud\t\t",                         _viewer, vtkFXViewer::ID_STYLE_GOURAUD, MENU_AUTOGRAY );
  new FXMenuRadio ( _viewmenu, "Phong\t\t",                           _viewer, vtkFXViewer::ID_STYLE_PHONG, MENU_AUTOGRAY );
  new FXMenuRadio ( _viewmenu, "Colors\t\t",                          _viewer, vtkFXViewer::ID_STYLE_TOGGLE_SCALAR_VISIBILITY, MENU_AUTOGRAY );
  new FXMenuRadio ( _viewmenu, "Global Immediate Rendering\t\t",      _viewer, vtkFXViewer::ID_TOGGLE_GLOBAL_IMMEDIATE_MODE_RENDERING, MENU_AUTOGRAY );
  new FXMenuSeparator ( _viewmenu );

  FXDockBar *dockbar = new FXDockBar ( rightdock, LAYOUT_FILL_Y | LAYOUT_SIDE_RIGHT, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2 );
  new FXMenuCheck ( _viewmenu, "Control panel\tCtrl-e", dockbar, FXWindow::ID_TOGGLESHOWN );
  FXHorizontalFrame *dockframe = new FXHorizontalFrame ( dockbar, LAYOUT_SIDE_TOP | LAYOUT_FILL_X, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
  new FXDockTitle ( dockframe, "Controls", dockbar, FXToolBar::ID_TOOLBARGRIP, LAYOUT_FILL_X | FRAME_SUNKEN | JUSTIFY_CENTER_X );
  new FXMDIDeleteButton ( dockframe, dockbar, FXWindow::ID_HIDE, LAYOUT_FILL_Y );

  _panels = new FXTabBook ( dockbar, NULL, 0, LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0 );

  new FXTabItem ( _panels, "Miscellaneous\tFind better name for this tab.\tFind better name for this tab." );
  FXVerticalFrame *misc_frame = new FXVerticalFrame ( _panels, LAYOUT_FILL_X | LAYOUT_FILL_Y );

  FXMatrix *transfer = new FXMatrix ( misc_frame, 3, MATRIX_BY_COLUMNS | LAYOUT_FILL_X | LAYOUT_TOP | LAYOUT_LEFT, 0, 0, 0, 0, 10, 10, 10, 10 );
  new FXLabel ( transfer, "Opacity Upper:" );
  new FXTextField ( transfer, 10, &_viewer->_opacityFunctionUpperValueTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderOpacityUpper = new FXRealSlider ( transfer, &_viewer->_opacityFunctionUpperValueTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 100 );
  rsliderOpacityUpper->setRange ( -1.0, 1.0 );
  new FXLabel ( transfer, "Opacity Lower:" );
  new FXTextField ( transfer, 10, &_viewer->_opacityFunctionLowerValueTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderOpacityLower = new FXRealSlider ( transfer, &_viewer->_opacityFunctionLowerValueTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 5, 100 );
  rsliderOpacityLower->setRange ( -1.0, 1.0 );
  new FXLabel ( transfer, "Color Upper:" );
  new FXTextField ( transfer, 10, &_viewer->_colorFunctionUpperValueTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderColorUpper = new FXRealSlider ( transfer, &_viewer->_colorFunctionUpperValueTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 100 );
  rsliderColorUpper->setRange ( -1.0, 1.0 );
  new FXLabel ( transfer, "Color Lower:" );
  new FXTextField ( transfer, 10, &_viewer->_colorFunctionLowerValueTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderColorLower = new FXRealSlider ( transfer, &_viewer->_colorFunctionLowerValueTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 100 );
  rsliderColorLower->setRange ( -1.0, 1.0 );

  new FXButton ( transfer, "Make Video", NULL, _viewer, vtkFXViewer::ID_MAKE_VIDEO, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );

  new FXButton ( transfer, "Load camera", NULL, _viewer, vtkFXViewer::ID_LOAD_CAMERA, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );
  new FXButton ( transfer, "Save camera", NULL, _viewer, vtkFXViewer::ID_SAVE_CAMERA, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );
  new FXButton ( transfer, "Reset camera", NULL, _viewer, vtkFXViewer::ID_RESET_CAMERA, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );

  new FXButton ( transfer, "Resize GLCanvas", NULL, _viewer, vtkFXViewer::ID_RESIZE_GLCANVAS, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );

  new FXButton ( transfer, "Clear Objects", NULL, _viewer, vtkFXViewer::ID_REMOVE_ALL_PROPS_FROM_RENDERER, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );
  new FXButton ( transfer, "Remove Mesh", NULL, _viewer, vtkFXViewer::ID_REMOVE_MESH_FROM_RENDERER, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );
  new FXButton ( transfer, "Remove Volume", NULL, _viewer, vtkFXViewer::ID_REMOVE_VOLUME_FROM_RENDERER, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );

  new FXButton ( transfer, "Isoline", NULL, _viewer, vtkFXViewer::ID_ISOLINE, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );

  new FXLabel ( transfer, "Overlay X Pos:" );
  new FXTextField ( transfer, 10, &_viewer->_overlayImageXPosTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderOverlayImageXPos = new FXRealSlider ( transfer, &_viewer->_overlayImageXPosTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 100 );
  rsliderOverlayImageXPos->setRange ( 0.0, 1.0 );
  new FXLabel ( transfer, "Overlay Y Pos:" );
  new FXTextField ( transfer, 10, &_viewer->_overlayImageYPosTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderOverlayImageYPos = new FXRealSlider ( transfer, &_viewer->_overlayImageYPosTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 100 );
  rsliderOverlayImageYPos->setRange ( 0.0, 1.0 );
  new FXLabel ( transfer, "Overlay Size:" );
  new FXTextField ( transfer, 10, &_viewer->_overlayImageSizeTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderOverlayImageSize = new FXRealSlider ( transfer, &_viewer->_overlayImageSizeTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 100 );
  rsliderOverlayImageSize->setRange ( 0.0, 1000.0 );
  rsliderOverlayImageSize->setIncrement ( 1. );
  new FXLabel ( transfer, "Overlay Opacity:" );
  new FXTextField ( transfer, 10, &_viewer->_overlayImageOpacityTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderOverlayImageOpacity = new FXRealSlider ( transfer, &_viewer->_overlayImageOpacityTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 100 );
  rsliderOverlayImageOpacity->setRange ( 0.0, 1.0 );

  new FXTabItem ( _panels, "View\tClipping, Rotation etc\tClipping, Rotation etc" );
  FXVerticalFrame *vframe = new FXVerticalFrame ( _panels, LAYOUT_FILL_X | LAYOUT_FILL_Y );

  FXGroupBox *group = new FXGroupBox ( vframe, "Origin", GROUPBOX_TITLE_LEFT | FRAME_GROOVE | LAYOUT_FILL_X );
  FXMatrix *group_mat = new FXMatrix ( group, 3, MATRIX_BY_COLUMNS | LAYOUT_TOP | LAYOUT_LEFT, 0, 0, 0, 0, 2, 2, 2, 2 );
  new FXLabel ( group_mat, "X:" );
  new FXTextField ( group_mat, 10, &_viewer->_plane1OriginXTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderOriginX = new FXRealSlider ( group_mat, &_viewer->_plane1OriginXTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 200 );
  rsliderOriginX->setRange ( -250.0, 250.0 ); // think about whether we want 0..255 or 0..1 here!

  new FXLabel ( group_mat, "Y:" );
  new FXTextField ( group_mat, 10, &_viewer->_plane1OriginYTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderOriginY = new FXRealSlider ( group_mat, &_viewer->_plane1OriginYTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 200 );
  rsliderOriginY->setRange ( -250.0, 250.0 );
  new FXLabel ( group_mat, "Z:" );
  new FXTextField ( group_mat, 10, &_viewer->_plane1OriginZTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderOriginZ = new FXRealSlider ( group_mat, &_viewer->_plane1OriginZTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 200 );
  rsliderOriginZ->setRange ( -200.0, 200.0 );


  FXGroupBox *group_normal = new FXGroupBox ( vframe, "Normal", GROUPBOX_TITLE_LEFT | FRAME_GROOVE | LAYOUT_FILL_X );
  FXMatrix *group_mat_normal = new FXMatrix ( group_normal, 3, MATRIX_BY_COLUMNS | LAYOUT_TOP | LAYOUT_LEFT, 0, 0, 0, 0, 2, 2, 2, 2 );
  new FXLabel ( group_mat_normal, "X:" );
  new FXTextField ( group_mat_normal, 10, &_viewer->_plane1NormalXTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderNormalX = new FXRealSlider ( group_mat_normal, &_viewer->_plane1NormalXTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 200 );
  rsliderNormalX->setRange ( -2.0, 2.0 );
  new FXLabel ( group_mat_normal, "Y:" );
  new FXTextField ( group_mat_normal, 10, &_viewer->_plane1NormalYTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderNormalY = new FXRealSlider ( group_mat_normal, &_viewer->_plane1NormalYTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 200 );
  rsliderNormalY->setRange ( -2.0, 2.0 );
  new FXLabel ( group_mat_normal, "Z:" );
  new FXTextField ( group_mat_normal, 10, &_viewer->_plane1NormalZTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderNormalZ = new FXRealSlider ( group_mat_normal, &_viewer->_plane1NormalZTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 200 );
  rsliderNormalZ->setRange ( -2.0, 2.0 );


  new FXCheckButton ( vframe, "Plane1", &_viewer->_plane1ActiveTarget, FXDataTarget::ID_VALUE, ICON_BEFORE_TEXT );
  new FXButton ( vframe, "Load plane1", NULL, _viewer, vtkFXViewer::ID_LOAD_PLANE1, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );
  new FXButton ( vframe, "Save plane1", NULL, _viewer, vtkFXViewer::ID_SAVE_PLANE1, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );

  FXHorizontalFrame *rotZoomHFrame = new FXHorizontalFrame ( vframe, LAYOUT_FILL_X | FRAME_GROOVE );

  FXGroupBox *rotation = new FXGroupBox ( rotZoomHFrame, "Rotation", GROUPBOX_TITLE_LEFT | FRAME_GROOVE | LAYOUT_FILL_X );
  FXMatrix *rotation_mat = new FXMatrix ( rotation, 3, MATRIX_BY_COLUMNS | LAYOUT_TOP | LAYOUT_LEFT, 0, 0, 0, 0, 2, 2, 2, 2 );
  new FXLabel ( rotation_mat, "X:" );
  new FXTextField ( rotation_mat, 10, &_viewer->_rotXAngleTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXDial* x_dial = new FXDial ( rotation_mat, &_viewer->_rotXAngleTarget, FXDataTarget::ID_VALUE, FRAME_SUNKEN | FRAME_THICK | DIAL_CYCLIC | DIAL_HORIZONTAL | LAYOUT_FIX_WIDTH | LAYOUT_FIX_HEIGHT | LAYOUT_CENTER_Y, 0, 0, 160, 14, 0, 0, 0, 0 );
  x_dial->setTipText ( "Rotate about X" );
  x_dial->setRange ( -180, 180 );
  x_dial->setNotchOffset ( 900 );

  new FXLabel ( rotation_mat, "Y:" );
  new FXTextField ( rotation_mat, 10, &_viewer->_rotYAngleTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXDial* y_dial = new FXDial ( rotation_mat, &_viewer->_rotYAngleTarget, FXDataTarget::ID_VALUE, FRAME_SUNKEN | FRAME_THICK | DIAL_CYCLIC | DIAL_HORIZONTAL | LAYOUT_FIX_WIDTH | LAYOUT_FIX_HEIGHT | LAYOUT_CENTER_Y, 0, 0, 160, 14, 0, 0, 0, 0 );
  y_dial->setTipText ( "Rotate about Y" );
  y_dial->setRange ( -180, 180 );
  y_dial->setNotchOffset ( 900 );

  new FXLabel ( rotation_mat, "Z:" );
  new FXTextField ( rotation_mat, 10, &_viewer->_rotZAngleTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXDial* z_dial = new FXDial ( rotation_mat, &_viewer->_rotZAngleTarget, FXDataTarget::ID_VALUE, FRAME_SUNKEN | FRAME_THICK | DIAL_CYCLIC | DIAL_HORIZONTAL | LAYOUT_FIX_WIDTH | LAYOUT_FIX_HEIGHT | LAYOUT_CENTER_Y, 0, 0, 160, 14, 0, 0, 0, 0 );
  z_dial->setTipText ( "Rotate about Z" );
  z_dial->setRange ( -180, 180 );
  z_dial->setNotchOffset ( 900 );

  FXGroupBox *zoom = new FXGroupBox ( rotZoomHFrame, "Zoom", GROUPBOX_TITLE_LEFT | FRAME_GROOVE | LAYOUT_FILL_X );
  new FXTextField ( zoom, 10, &_viewer->_zoomTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderZoom = new FXRealSlider ( zoom, &_viewer->_zoomTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 50 );
  rsliderZoom->setRange ( 0.2, 5.0 );

  FXHorizontalFrame *camButonsHFrame = new FXHorizontalFrame ( rotation, LAYOUT_FILL_X );
  new FXButton ( camButonsHFrame, "GetCameraAngles", NULL, _viewer, vtkFXViewer::ID_GET_CAMERA_ANGLES,  LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );
  new FXButton ( camButonsHFrame, "LimitXAngle", NULL, _viewer, vtkFXViewer::ID_LIMIT_CAMERA_X_ANGLE,  LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );

  for ( int i = 0; i < 4; ++i ) { // finally use all 8
    char lightName[32];
    sprintf ( lightName, "Light %d:", i );

    FXGroupBox *light = new FXGroupBox ( vframe, lightName, GROUPBOX_TITLE_LEFT | FRAME_GROOVE | LAYOUT_FILL_X );

    FXMatrix *light_matS = new FXMatrix ( light, 3, MATRIX_BY_COLUMNS | LAYOUT_TOP | LAYOUT_LEFT, 0, 0, 0, 0, 2, 2, 2, 2 );
    new FXCheckButton ( light_matS, lightName, _viewer->_lightSwitchTarget[i], FXDataTarget::ID_VALUE, ICON_AFTER_TEXT );
    new FXLabel ( light_matS, "Color:" );
    new FXColorWell ( light_matS, _viewer->getLightColor(i), _viewer->_lightColorTarget[i], FXDataTarget::ID_VALUE, LAYOUT_FILL_ROW );

    FXMatrix *light_matP = new FXMatrix ( light, 3, MATRIX_BY_COLUMNS | LAYOUT_TOP | LAYOUT_LEFT, 0, 0, 0, 0, 2, 2, 2, 2 );
    new FXLabel ( light_matP, "Elevation:" );
    new FXTextField ( light_matP, 10, _viewer->_lightElevationTarget[i], FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
#ifndef SIDE_BY_SIDE_STEREO_RENDERING
    FXDial* le_dial = new FXDial ( light_matP, _viewer->_lightElevationTarget[i], FXDataTarget::ID_VALUE, FRAME_SUNKEN | FRAME_THICK | DIAL_CYCLIC | DIAL_HORIZONTAL | LAYOUT_FIX_WIDTH | LAYOUT_FIX_HEIGHT | LAYOUT_CENTER_Y, 0, 0, 160, 14, 0, 0, 0, 0 );
    le_dial->setTipText ( "Elevation" );
    le_dial->setRange ( -180, 180 );
    le_dial->setNotchOffset ( 900 );
#else
    new FXLabel ( light_matP, "wheel not available" ); // for some reason, this leads to transparent control panel and menu in stereo mode. TODO: find out why, fix this!
#endif
    new FXLabel ( light_matP, "Azimuth:" );
    new FXTextField ( light_matP, 10, _viewer->_lightAzimuthTarget[i], FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
#ifndef SIDE_BY_SIDE_STEREO_RENDERING
    FXDial* la_dial = new FXDial ( light_matP, _viewer->_lightAzimuthTarget[i], FXDataTarget::ID_VALUE, FRAME_SUNKEN | FRAME_THICK | DIAL_CYCLIC | DIAL_HORIZONTAL | LAYOUT_FIX_WIDTH | LAYOUT_FIX_HEIGHT | LAYOUT_CENTER_Y, 0, 0, 160, 14, 0, 0, 0, 0 );
    la_dial->setTipText ( "Azimuth" );
    la_dial->setRange ( -180, 180 );
    la_dial->setNotchOffset ( 900 );
#else
    new FXLabel ( light_matP, "wheel not available" );
#endif
  }

  new FXTabItem ( _panels, "Objects\t\tManaging objects." );
  FXVerticalFrame *objects = new FXVerticalFrame ( _panels, LAYOUT_SIDE_TOP | LAYOUT_FILL_X | LAYOUT_FILL_Y );
  FXHorizontalFrame *allObjectsControlFrame = new FXHorizontalFrame ( objects, LAYOUT_SIDE_TOP | LAYOUT_FILL_X );
  new FXButton ( allObjectsControlFrame, "&Clear all", NULL, _viewer->getPolyDataGroupPointer(), PolyDataGroup::ID_CLEAR, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );
  new FXButton ( allObjectsControlFrame, "Show all", NULL, _viewer->getPolyDataGroupPointer(), PolyDataGroup::ID_SHOW_ALL, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );
  FXHorizontalFrame *_singleActiveObjectControlFrame = new FXHorizontalFrame ( objects, LAYOUT_SIDE_TOP | LAYOUT_FILL_X );
  new FXLabel ( _singleActiveObjectControlFrame, "\"Single\" active object:" );
  FXSpinner *objectSpinner = new FXSpinner ( _singleActiveObjectControlFrame, 5, _viewer->getPolyDataGroupPointer(), PolyDataGroup::ID_SINGLE_ACTIVE_OBJECT_AND_HIDE_OTHERS, SPIN_CYCLIC | FRAME_SUNKEN | FRAME_THICK | LAYOUT_CENTER_Y | LAYOUT_FILL_ROW );
  objectSpinner->setTipText ( "Seclects the active object, makes it visible and hides all others." );
  new FXButton ( _singleActiveObjectControlFrame, "Save object to .VTK", NULL, _viewer->getPolyDataGroupPointer(), PolyDataGroup::ID_SAVE_ACTIVE_POLYDATA_TO_VTK_FILE, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );

  new FXLabel ( objects, "Currently active object:" );
  FXListBox* objectListBox = new FXListBox ( objects, _viewer->getPolyDataGroupPointer(), PolyDataGroup::ID_SINGLE_ACTIVE_OBJECT );
  objectListBox->setNumVisible ( 5 );

  new FXLabel ( objects, "Object visibility:" );
  FXScrollWindow *objectScrollArea = new FXScrollWindow ( objects, LAYOUT_FILL_X | LAYOUT_FILL_Y );
  FXVerticalFrame *objectsScrollFrame = new FXVerticalFrame ( objectScrollArea, LAYOUT_SIDE_TOP | LAYOUT_FILL_X | LAYOUT_FILL_Y );
  _viewer->getPolyDataGroupPointer()->setButtonList ( objectsScrollFrame );
  _viewer->getPolyDataGroupPointer()->setObjectDropDownList ( objectListBox );


  new FXTabItem ( _panels, "Colors\tColormapping\tColormapping" );
  FXVerticalFrame *_colorControlFrame = new FXVerticalFrame ( _panels, LAYOUT_SIDE_TOP | LAYOUT_FILL_X );
  FXMatrix *_ccfLM = new FXMatrix ( _colorControlFrame, 3, FRAME_THICK | FRAME_RAISED | MATRIX_BY_COLUMNS | LAYOUT_FILL_X | LAYOUT_TOP | LAYOUT_LEFT, 0, 0, 0, 0, 5, 5, 5, 5 );

  new FXLabel ( _ccfLM, "Background Color:" );
  new FXColorWell ( _ccfLM, FXRGB ( 0, 0, 0 ), &_viewer->_backgroundColorTarget, FXDataTarget::ID_VALUE );
  new FXButton ( _ccfLM, "Add zero data", NULL, _viewer->getPolyDataGroupPointer(), PolyDataGroup::ID_ADD_CONST_DATA_TO_ACTIVE_POLYDATA, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );

  new FXLabel ( _ccfLM, "Color Transition" );
  FXListBox* colormapList = new FXListBox ( _ccfLM, & ( _viewer->getPolyDataGroupPointer()->_colorMapTarget ), FXDataTarget::ID_VALUE );
  colormapList->appendItem ( "blue to red HSV" );
  colormapList->appendItem ( "blue to red" );
  colormapList->appendItem ( "topographical" );
  colormapList->appendItem ( "red to blue HSV" );
  colormapList->appendItem ( "black to white" );
  colormapList->appendItem ( "white to black" );
  colormapList->appendItem ( "white to blue" );
  colormapList->appendItem ( "white/green/blue/red" );
  colormapList->appendItem ( "silver to gold" );
  colormapList->appendItem ( "constant" );
  colormapList->setNumVisible ( 10 );

  new FXButton ( _ccfLM, "Maximize range", NULL, _viewer->getPolyDataGroupPointer(), PolyDataGroup::ID_MAXIMIZE_COLOR_RANGE_ACTIVE_POLYDATA, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );

  new FXLabel ( _ccfLM, "Lower value for colormap:" );
  new FXTextField ( _ccfLM, 10, & ( _viewer->getPolyDataGroupPointer()->_lowerValueForColormapTarget ), FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  new FXLabel ( _ccfLM, "" );

  new FXLabel ( _ccfLM, "Upper value for colormap:" );
  new FXTextField ( _ccfLM, 10, & ( _viewer->getPolyDataGroupPointer()->_upperValueForColormapTarget ), FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  new FXLabel ( _ccfLM, "" );

  new FXLabel ( _ccfLM, "Opacity:" );
  new FXTextField ( _ccfLM, 10, &_viewer->getPolyDataGroupPointer()->_opacityValueForColormapTarget, FXDataTarget::ID_VALUE, TEXTFIELD_REAL | JUSTIFY_RIGHT | FRAME_SUNKEN | FRAME_THICK );
  FXRealSlider *rsliderOpacity = new FXRealSlider ( _ccfLM, &_viewer->getPolyDataGroupPointer()->_opacityValueForColormapTarget, FXDataTarget::ID_VALUE, LAYOUT_CENTER_Y | LAYOUT_FILL_X | LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH, 0, 0, 100 );
  rsliderOpacity->setRange ( -0.0, 1.0 );

  new FXLabel ( _ccfLM, "Constant Object Color:" );
  new FXColorWell ( _ccfLM, FXRGB ( 0, 0, 0 ), &_viewer->getPolyDataGroupPointer()->_objectColorTarget, FXDataTarget::ID_VALUE );

  new FXButton ( _colorControlFrame, "Apply to all", NULL, _viewer->getPolyDataGroupPointer(), PolyDataGroup::ID_APPLY_FULL_COLORMAP_PROPERTIES_TO_ALL, LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );


  // shortcuts
  getAccelTable()->addAccel ( parseAccel ( "pgup" ), objectSpinner, FXSEL ( SEL_COMMAND, FXSpinner::ID_INCREMENT ) );
  getAccelTable()->addAccel ( parseAccel ( "pgdn" ), objectSpinner, FXSEL ( SEL_COMMAND, FXSpinner::ID_DECREMENT ) );
  getAccelTable()->addAccel ( parseAccel ( "space" ), objectSpinner, FXSEL ( SEL_COMMAND, FXSpinner::ID_INCREMENT ) );
  getAccelTable()->addAccel ( parseAccel ( "backspace" ), objectSpinner, FXSEL ( SEL_COMMAND, FXSpinner::ID_DECREMENT ) );
  getAccelTable()->addAccel ( parseAccel ( "Ctl-q" ), getApp(), FXSEL ( SEL_COMMAND, FXApp::ID_QUIT ) );

  // Make a tooltip
  // Necessary to show tool tips on controls set by setTipText.
  new FXToolTip(getApp());
}


vtkFXGui::~vtkFXGui ( ) {
  delete ( _filemenu );
  delete ( _viewmenu );
  // We need to delete the _viewer explicitly here (before the main window starts do destroy
  // all GUI elements, otherwise _viewer->_polyDataGroup->getObjectDropDownList()
  // returns an invalid pointer when _viewer is deleted automatically later causing
  // PolyDataQC::~PolyDataQC to crash when _polyDataGroup is destroyed.
  delete ( _viewer );
}


vtkFXGui::vtkFXGui( ) {}
