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

#ifndef __POLYDATAGROUP
#define __POLYDATAGROUP

#include <vector>

#include <foxIncludes.h>
#include <vtkIncludes.h>


#include "vtkFXRenderWindowInteractor.h"
#include "FXVTKCanvas.h"

#include "vtkFXViewer.h"
#include "vtkFXGui.h"
#include "polyDataQC.h"

// todo: get rid of this define! OS is not sure whether this really works: is the static member guaranteed to be initialized early enough for the enum?
// static const unsigned short VTKFOX_NUMBER_OF_MANAGEABLE_OBJECTS = 128;
#define VTKFOX_NUMBER_OF_MANAGEABLE_OBJECTS 128

// todo: think of a better solution.
using namespace std;

#ifdef SIDE_BY_SIDE_STEREO_RENDERING
#define RENDERER vtkStereoRenderer
#else
#define RENDERER vtkRenderer
#endif



// the polyDataGroup
// =================

//! manages all opened meshes
class PolyDataGroup: public FXObject {
  FXDECLARE ( PolyDataGroup )

public:
  // Messages for our class
  enum{
    ID_CLEAR,
    ID_SHOW_ALL,
    ID_APPLY_FULL_COLORMAP_PROPERTIES_TO_ALL,
    ID_COLORMAP_UPDATE,
    ID_SINGLE_ACTIVE_OBJECT_AND_HIDE_OTHERS,
    ID_SINGLE_ACTIVE_OBJECT,
    ID_POLYDATA_MESSAGE,
    ID_SAVE_ACTIVE_POLYDATA_TO_VTK_FILE,
    ID_MAXIMIZE_COLOR_RANGE_ACTIVE_POLYDATA,
    ID_ADD_CONST_DATA_TO_ACTIVE_POLYDATA,
    ID_SET_OBJECT_VISIBILITY,
    ID_LAST = ID_SET_OBJECT_VISIBILITY + VTKFOX_NUMBER_OF_MANAGEABLE_OBJECTS
  };

protected:
  vtkFXGui *     _mainWindow;
  RENDERER *_renderer;
  FXComposite* _buttonList;
  FXListBox* _objectDropDownList;

  FXint _singleActiveObject;
  FXint _colorMap;
  FXfloat _lowerValueForColormap, _upperValueForColormap, _opacityValueForColormap; // poor encapsulation TODO: write access functions.
  FXColor _objectColor;

  vector<PolyDataQC *> _polyData;

public:
  FXDataTarget _lowerValueForColormapTarget, _upperValueForColormapTarget, _opacityValueForColormapTarget;
  FXDataTarget _colorMapTarget;
  FXDataTarget _objectColorTarget;

protected:
  PolyDataGroup() {}

  // not implemented (should cause compile error if object is copied)
  PolyDataGroup ( const PolyDataGroup& FVC );
  PolyDataGroup& operator= ( const PolyDataGroup& FVC );

public:
  // Constructor
  PolyDataGroup ( vtkFXGui* MainWindow, RENDERER *Renderer )
    : _mainWindow ( MainWindow ),
      _renderer ( Renderer ),
      _buttonList ( NULL ),
      _singleActiveObject ( -1 ),
      _colorMap ( 0 ),
      _lowerValueForColormap ( 0.0 ),
      _upperValueForColormap ( 1.0 ),
      _opacityValueForColormap ( 1.0 ), // reasonable default??
      _objectColor ( FXRGB ( 255, 255, 255 ) ),
      _lowerValueForColormapTarget ( _lowerValueForColormap, this, ID_COLORMAP_UPDATE ),
      _upperValueForColormapTarget ( _upperValueForColormap, this, ID_COLORMAP_UPDATE ),
      _opacityValueForColormapTarget ( _opacityValueForColormap, this, ID_COLORMAP_UPDATE ),
      _colorMapTarget ( _colorMap, this, ID_COLORMAP_UPDATE ),
      _objectColorTarget ( _objectColor ) {}

  ~PolyDataGroup() {
    // invalidate drop down list pointer to prevent
    // PolyDataQC objects from accessing this already
    // destructed object.
    _objectDropDownList = NULL;
    deleteAllMeshData();
  }

  void deleteAllMeshData();

  void setSingleActiveObject ( const FXint number );

  //! set single object to be visible and active, other objects are switched off
  void setSingleActiveObjectAndHideOthers ( const FXint number );

  void setButtonList ( FXComposite *ButtonList ) {
    _buttonList = ButtonList;
  }

  void setObjectDropDownList ( FXListBox *ObjectDropDownList ) {
    _objectDropDownList = ObjectDropDownList;
  }

  FXComposite * getButtonList ()          {  return _buttonList;  }
  FXListBox *  getObjectDropDownList ()   {  return _objectDropDownList;  }
  vtkFXGui * getMainWindow ()             {  return _mainWindow;  }
  RENDERER * getRenderer ()               {  return _renderer;  }

  //! \note Memory management is handled by VTK's reference counting. If you don't need
  //! the vtkPolyData anywhere else, call Delete() on the pointer once it has been passed to PolyDataGroup.
  //!
  //! adds a new mesh which is set to be visible, all others are switched off, buttons are generated
  void addPolyData ( vtkPolyData* PolyData, const FXString & DataName );

  void setColorMapRange ( FXfloat lower, FXfloat upper ) {
    // The correspoding data target don't need to be explicitly notified of the changed values.
    // They have a reference to their watched value and always return the correct value when asked
    // with SEL_UPDATE. If the values are changed it does not then any messages though, so we need
    // to call doUpdateColormap() here.
    _lowerValueForColormap = lower;
    _upperValueForColormap = upper;
    this->doUpdateColormap();
  }

  void setColorMap ( FXint ColorMap ) {
    _colorMap = ColorMap;
    doUpdateColormap();
  }

  long onCmdClear                           ( FXObject*, FXSelector, void* );
  long onCmdShowAll                         ( FXObject*, FXSelector, void* );
  long onCmdApplyFullColormapPropertiesToAll( FXObject*, FXSelector, void* );
  long onCmdSingleActiveObjectAndHideOthers ( FXObject*, FXSelector, void* );
  long onUpdSingleActiveObject              ( FXObject*, FXSelector, void* );
  long onCmdSingleActiveObject              ( FXObject*, FXSelector, void* );
  long onUpdObjectVisibilityStatus          ( FXObject*, FXSelector, void* );
  long onCmdObjectVisibilityStatus          ( FXObject*, FXSelector, void* );
  long onCmdSaveActivePolydataToVTKFile     ( FXObject*, FXSelector, void* );
  long onCmdMaximizeColorRangeActivePolydata( FXObject*, FXSelector, void* );
  long onCmdAddConstDataToActivePolydata    ( FXObject*, FXSelector, void* );

  vtkPolyData* getCurrentPolyData() {
    if ( _singleActiveObject != -1 )
      return getPolyData(_singleActiveObject)->getPolyData();
    else
      return NULL;
  }

  vtkActor* getCurrentActor() {
    if ( _singleActiveObject != -1 )
      return getPolyData(_singleActiveObject)->getActor();
    else
      return NULL;
  }

  vtkPolyDataMapper* getCurrentPolyDataMapper() {
    if ( _singleActiveObject != -1 )
      return getPolyData(_singleActiveObject)->getMapper();
    else
      return NULL;
  }

  vtkProperty* getCurrentProperty() {
    if ( _singleActiveObject != -1 )
      return getPolyData(_singleActiveObject)->getProperties();
    else
      return NULL;
  }

  vtkPolyDataNormals* getCurrentPolyDataNormals() {
    if ( _singleActiveObject != -1 )
      return getPolyData(_singleActiveObject)->getNormals();
    else
      return NULL;
  }

  vtkStripper* getCurrentStripper() {
    if ( _singleActiveObject != -1 )
      return getPolyData(_singleActiveObject)->getStripper();
    else
      return NULL;
  }

  long onColormapUpdate ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*ptr*/ ) {
    return ( doUpdateColormap() );
  }

  long doUpdateColormap ( const int I = - 1 ) ;

  void ClippingPlaneToAll ( vtkPlane *_plane1, const bool active ) {
    for ( int i = 0; i < size()  ; ++i ) {
      if ( active ) {
        getPolyData(i)->getMapper()->AddClippingPlane( _plane1 );
      } else {
        getPolyData(i)->getMapper()->RemoveClippingPlane( _plane1 );
      }
    }
  }

  long polyDataMessage ( FXObject* sender, FXSelector sel, void* data ) {
    return getPolyData(_singleActiveObject)->handle ( sender, sel, data );
  }

  PolyDataQC * getPolyData(int i)             {  return (i < 0 ? NULL : _polyData[i]);  }
  int size() const                            {  return _polyData.size();  }
  void pushBackAndShow ( PolyDataQC * polyData )    {  _polyData.push_back ( polyData );   setSingleActiveObjectAndHideOthers ( size() - 1 ); }
  FXTabBook * getControlPanel()               {  return getMainWindow()->getControlPanel();  }

  int getNewObjectVisibilityId() const   {  return ID_SET_OBJECT_VISIBILITY + size();  }
};

#endif
