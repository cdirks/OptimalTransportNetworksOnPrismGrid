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

#include "polyDataQC.h"

#include <levelsetToTriangMesh.h>
#include "triangMeshToVTK.h"
#include "polyDataGroup.h"

// --------------------------------------------------------------------------

PolyDataQC::PolyDataQC ( PolyDataGroup *      owningGroup,
                         bool                 visible,
                         FXDataTarget *       visibleDataTarget,
                         FXFrame *            button,
                         FXString             description,
                         vtkPolyData *        polyData,
                         vtkActor *           actor,
                         vtkPolyDataMapper *  mapper,
                         vtkProperty *        properties,
                         vtkPolyDataNormals * normals,
                         vtkStripper *        stripper )

    : _owningGroup       ( owningGroup )
    , _visible           ( visible )
    , _visibleDataTarget ( visibleDataTarget )
    , _button            ( button )
    , _description       ( description )
    , _polyData          ( polyData )
    , _actor             ( actor )
    , _mapper            ( mapper )
    , _properties        ( properties )
    , _normals           ( normals )
    , _stripper          ( stripper )
    , _propertiesTab     ( NULL )
    , _propertiesFrame   ( NULL ) {}

// --------------------------------------------------------------------------

PolyDataQC::PolyDataQC ( PolyDataGroup *      owningGroup,
                         bool                 visible,
                         FXString             description,
                         vtkPolyData *        polyData,
                         vtkActor *           actor,
                         vtkPolyDataMapper *  mapper,
                         vtkProperty *        properties,
                         vtkPolyDataNormals * normals,
                         vtkStripper *        stripper )

    : _owningGroup       ( owningGroup )
    , _visible           ( visible )
    , _description       ( description )
    , _polyData          ( polyData )
    , _actor             ( actor )
    , _mapper            ( mapper )
    , _properties        ( properties )
    , _normals           ( normals )
    , _stripper          ( stripper )
    , _propertiesTab     ( NULL )
    , _propertiesFrame   ( NULL ) {

  // if necessary, create button and drow down list item
  FXComposite * buttonList = _owningGroup->getButtonList();
  FXButton * button = new FXButton ( buttonList, description, NULL, _owningGroup,
                           _owningGroup->getNewObjectVisibilityId(),
                           LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );
  button->create();
  button->recalc();
  _owningGroup->getObjectDropDownList()->appendItem ( description );

  // set member variables
  _visibleDataTarget = new FXDataTarget ( _visible );
  _button            = button;
}

// --------------------------------------------------------------------------

PolyDataQC::PolyDataQC ( PolyDataGroup * owningGroup,
                         FXString        description,
                         vtkPolyData *   polyData )

    : _owningGroup  ( owningGroup )
    , _visible      ( false )
    , _description  ( description )
    , _polyData     ( polyData )
    , _propertiesTab     ( NULL )
    , _propertiesFrame   ( NULL ) {

  // if necessary, create button and drow down list item
  FXComposite * buttonList = owningGroup->getButtonList();
  FXButton * button = new FXButton ( buttonList, description, NULL, _owningGroup,
                           _owningGroup->getNewObjectVisibilityId(),
                           LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X );
  button->create();
  button->recalc();
  owningGroup->getObjectDropDownList()->appendItem ( description );

  // set member variables
  _visibleDataTarget = new FXDataTarget ( _visible );
  _button            = button;
  _actor             = vtkActor::New();
  _mapper            = vtkPolyDataMapper::New();
  _properties        = vtkProperty::New();
  _normals           = vtkPolyDataNormals::New();
  _stripper          = vtkStripper::New();

  connectVTKPipeline();
  show();
}

// --------------------------------------------------------------------------

PolyDataQC::~PolyDataQC () {
  // find our entry in object drop down list
  FXListBox * list = _owningGroup->getObjectDropDownList();

  // list pointer is only NULL when the window is already being
  // destroyed. In this case, perform no destruction of
  // visual items at all.
  if (list) {
    hide();

    FXint n = list->getNumItems();
    for (FXint i = 0; i < n; ++i)
      if (list->getItem(i) == _description) {
        list->removeItem(i);
        break;
      }

    // Do not care for NULL pointers, as "delete NULL"
    // has no effect.
    delete _visibleDataTarget;
    delete _button;
    delete _propertiesTab;
    delete _propertiesFrame;
  }

  // delete VTK pipeline
  //! \note Do NOT delete _polyData. It will be deleted by VTK's reference counting.
  if (_actor)      _actor      ->Delete ();
  if (_mapper)     _mapper     ->Delete ();
  if (_properties) _properties ->Delete ();
  if (_normals)    _normals    ->Delete ();
  if (_stripper)   _stripper   ->Delete ();
}

// --------------------------------------------------------------------------

void PolyDataQC::connectVTKPipeline() {
  _normals->SetInput ( _polyData );
  _normals->SetFeatureAngle ( 60.0 );

  // Generate triangle strips
  _stripper->SetInputConnection ( _normals->GetOutputPort() );

  _mapper->SetInputConnection ( _stripper->GetOutputPort() );
  _mapper->GlobalImmediateModeRenderingOn();

  _actor->SetMapper ( _mapper );
  _actor->SetProperty ( _properties );
}

// --------------------------------------------------------------------------

void PolyDataQC::updateVTKCanvas() {
  _owningGroup->getMainWindow()->getViewer()->updateVTKCanvas();
}

// --------------------------------------------------------------------------

void PolyDataQC::show() {

  if (!_visible) {
    _owningGroup->getRenderer()->AddViewProp ( _actor );
    _visible = true;
  }
}

// --------------------------------------------------------------------------

void PolyDataQC::hide() {

  if (_visible) {
    _owningGroup->getRenderer()->RemoveViewProp ( _actor );
    _visible = false;
  }
}

// --------------------------------------------------------------------------

void PolyDataQC::activate() {
  FXString title = FXString ( "vtkFOX - File: " ).append ( _description );
  _owningGroup->getMainWindow()->setTitle ( title );

  if (_propertiesTab)
    _propertiesTab->show();
}

// --------------------------------------------------------------------------

void PolyDataQC::inactivate() {
  if (_propertiesTab)
    _propertiesTab->hide();
}

// ==========================================================================

IsoSurfaceData::IsoSurfaceData ( PolyDataGroup *      owningGroup,
                                 FXString             description,
                                 vtkImageData *       levelData,
                                 FXdouble             levelValueToDisplay )
    : PolyDataQC ( owningGroup, false, description )
    , _levelValueToDisplay ( levelValueToDisplay )
    , _levelData ( levelData )
    , _contourFilter ( vtkContourFilter::New() )
    , _levelValueTarget ( new FXDataTarget ( _levelValueToDisplay,
                                             owningGroup,
                                             PolyDataGroup::ID_POLYDATA_MESSAGE ) ) {

  // set up visualization pipeline
  _contourFilter->SetInput ( _levelData );
  _contourFilter->SetValue ( 0, _levelValueToDisplay );

  _normals->SetInput ( _contourFilter->GetOutput() );
  _normals->SetFeatureAngle ( 60.0 );

  _mapper->SetInput ( _normals->GetOutput() );
  _mapper->GlobalImmediateModeRenderingOn();

  _actor->SetMapper ( _mapper );
  _actor->SetProperty ( _properties );

  // *** create properties control frame ***
  FXTabBook * panels = _owningGroup->getControlPanel();
  _propertiesTab = new FXTabItem ( panels, "Level" );
  _propertiesFrame = new FXVerticalFrame ( panels,
                                           LAYOUT_SIDE_TOP |
                                           LAYOUT_FILL_X | LAYOUT_FILL_Y );
  // layout matrix 3 x 1
  FXMatrix * matrix
    = new FXMatrix ( _propertiesFrame,
                     3,
                     MATRIX_BY_COLUMNS | LAYOUT_FILL_X |
                     LAYOUT_TOP | LAYOUT_LEFT,
                     0, 0, 0, 0, 10, 10, 10, 10 );

  // description label
  new FXLabel ( matrix, "Level value:" );
  // text field for manual level value input
  new FXTextField ( matrix,
                    10,
                    _levelValueTarget,
                    FXDataTarget::ID_VALUE,
                    TEXTFIELD_REAL | JUSTIFY_RIGHT |
                    FRAME_SUNKEN | FRAME_THICK );

  // slider for level choice via mouse control
  FXRealSlider *rsliderLevelValue
    = new FXRealSlider ( matrix,
                         _levelValueTarget,
                         FXDataTarget::ID_VALUE,
                         LAYOUT_CENTER_Y | LAYOUT_FILL_X |
                         LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH,
                         0, 0, 100 );

  // set slider range according to existent level values
  double range[2];
  _levelData->GetScalarRange ( range );
  rsliderLevelValue->setRange ( range[0], range[1] );

  _propertiesTab->create();
  _propertiesFrame->create();
}

// --------------------------------------------------------------------------

IsoSurfaceData::~IsoSurfaceData() {
  if (_levelData)     _levelData->Delete();
  if (_contourFilter) _contourFilter->Delete();
}

// --------------------------------------------------------------------------

long IsoSurfaceData::
handle ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  _contourFilter->SetValue ( 0, _levelValueToDisplay );
  return 1;
}

// ==========================================================================

IsoSurfaceDataWithTextureMap::
IsoSurfaceDataWithTextureMap ( PolyDataGroup * owningGroup,
                               FXString        description,
                               std::string     levelsetFilename,
                               std::string     tCoordXFilename,
                               std::string     tCoordYFilename,
                               std::string     tBitmapFilename,
                               FXdouble        levelValueToDisplay )

    : PolyDataQC ( owningGroup, false, description, NULL )
    , _levelData ( levelsetFilename )
    , _tCoordX   ( tCoordXFilename )
    , _tCoordY   ( tCoordYFilename )
    , _levelValueToDisplay ( levelValueToDisplay )
    , _levelValueTarget ( new FXDataTarget ( _levelValueToDisplay,
                                             owningGroup,
                                             PolyDataGroup::ID_POLYDATA_MESSAGE ) ) {

  qc::ScalarArray<FXfloat, qc::QC_3D> levelDataForLevelValueZero (_levelData);
  levelDataForLevelValueZero.addToAll ( -_levelValueToDisplay );

  // create VTK pipeline
  _polyData = createSurfaceWithTCoordsFromZeroLevelSet
                  ( levelDataForLevelValueZero,
                    _tCoordX,
                    _tCoordY );

  _tBitmapReader = vtkBMPReader::New();
  _tBitmapReader->SetFileName ( tBitmapFilename.c_str() );

  _texture = vtkTexture::New();
  _texture->SetInputConnection ( _tBitmapReader->GetOutputPort() );
  _texture->InterpolateOn();

  connectVTKPipeline();
  _actor->SetTexture ( _texture );

  // create properties control frame
  FXTabBook * panels = _owningGroup->getControlPanel();
  _propertiesTab = new FXTabItem ( panels, "Textured Level Set" );

  _propertiesFrame = new FXVerticalFrame ( panels,
                                           LAYOUT_SIDE_TOP |
                                           LAYOUT_FILL_X | LAYOUT_FILL_Y );

  // layout matrix 3 x 1
  FXMatrix * matrix
    = new FXMatrix ( _propertiesFrame,
                     3,
                     MATRIX_BY_COLUMNS | LAYOUT_FILL_X |
                     LAYOUT_TOP | LAYOUT_LEFT,
                     0, 0, 0, 0, 10, 10, 10, 10 );

  // description label
  new FXLabel ( matrix, "Level value:" );
  // text field for manual level value input
  new FXTextField ( matrix,
                    10,
                    _levelValueTarget,
                    FXDataTarget::ID_VALUE,
                    TEXTFIELD_REAL | JUSTIFY_RIGHT |
                    FRAME_SUNKEN | FRAME_THICK );

  // slider for level choice via mouse control
  FXRealSlider *rsliderLevelValue
    = new FXRealSlider ( matrix,
                         _levelValueTarget,
                         FXDataTarget::ID_VALUE,
                         LAYOUT_CENTER_Y | LAYOUT_FILL_X |
                         LAYOUT_FILL_ROW | LAYOUT_FIX_WIDTH,
                         0, 0, 100 );

  // set slider range according to existent level values
  double range[2] = { _levelData.getMinValue(), _levelData.getMaxValue() };
  rsliderLevelValue->setRange ( range[0], range[1] );

  _propertiesTab->create();
  _propertiesFrame->create();
  show();
}
// --------------------------------------------------------------------------

vtkPolyData * IsoSurfaceDataWithTextureMap::
createSurfaceWithTCoordsFromZeroLevelSet ( const qc::ScalarArray<FXfloat, qc::QC_3D> & levelData,
                                           const qc::ScalarArray<FXfloat, qc::QC_3D> & tCoordX,
                                           const qc::ScalarArray<FXfloat, qc::QC_3D> & tCoordY ) {

  qcsm::LevelsetToTriangMesh<FXfloat> converter;

  // extract TriangMesh
  converter.apply(levelData, _levelSurface);

  // add texture data
  converter.addVertexData ( _levelSurface, tCoordX, "tCoordX" );
  converter.addVertexData ( _levelSurface, tCoordY, "tCoordY" );

  // convert TriangMesh into vtk surface, add texture coords
  FXfloat min, max;
  vtkPolyData * polyData
    = aol::TriangMeshToVTK<FXfloat>::makePolyData ( _levelSurface, min, max, false );
  return aol::TriangMeshToVTK<FXfloat>::setPolyDataTCoords ( polyData,
                                            _levelSurface.getVertexData(0),
                                            _levelSurface.getVertexData(1) );
}
// --------------------------------------------------------------------------

IsoSurfaceDataWithTextureMap::~IsoSurfaceDataWithTextureMap () {
  // destroy VTK objects (in fact, only decrease reference counts)
  _tBitmapReader->Delete();
  _texture->Delete();
}
// --------------------------------------------------------------------------

long IsoSurfaceDataWithTextureMap::
handle ( FXObject* sender, FXSelector /*sel*/, void* /*data*/ ) {
  if (sender == _levelValueTarget)
    updateLevelValueToDisplay();
  return 1;
}
// --------------------------------------------------------------------------

void IsoSurfaceDataWithTextureMap::updateLevelValueToDisplay() {
  qc::ScalarArray<FXfloat, qc::QC_3D> levelDataForLevelValueZero (_levelData);
  levelDataForLevelValueZero.addToAll ( -_levelValueToDisplay );
  _polyData = createSurfaceWithTCoordsFromZeroLevelSet
                    ( levelDataForLevelValueZero,
                      _tCoordX,
                      _tCoordY );
  connectVTKPipeline();
  updateVTKCanvas();
}
// --------------------------------------------------------------------------

