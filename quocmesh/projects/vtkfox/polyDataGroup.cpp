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

#include "polyDataGroup.h"
#include "vtkFXViewer.h"
#include "vtkVisualize.h"

// PolyDataGroup implementation
// ============================

FXDEFMAP ( PolyDataGroup ) PolyDataGroupMap[] = { //________Message_Type_____________________ID____________Message_Handler_______
                                                  FXMAPFUNC  ( SEL_COMMAND, PolyDataGroup::ID_CLEAR, PolyDataGroup::onCmdClear ),
                                                  FXMAPFUNC  ( SEL_COMMAND, PolyDataGroup::ID_SHOW_ALL, PolyDataGroup::onCmdShowAll ),
                                                  FXMAPFUNC  ( SEL_COMMAND, PolyDataGroup::ID_APPLY_FULL_COLORMAP_PROPERTIES_TO_ALL, PolyDataGroup::onCmdApplyFullColormapPropertiesToAll ),
                                                  FXMAPFUNC  ( SEL_COMMAND, PolyDataGroup::ID_COLORMAP_UPDATE, PolyDataGroup::onColormapUpdate ),
                                                  FXMAPFUNC  ( SEL_UPDATE,  PolyDataGroup::ID_SINGLE_ACTIVE_OBJECT_AND_HIDE_OTHERS, PolyDataGroup::onUpdSingleActiveObject ),
                                                  FXMAPFUNC  ( SEL_COMMAND, PolyDataGroup::ID_SINGLE_ACTIVE_OBJECT_AND_HIDE_OTHERS, PolyDataGroup::onCmdSingleActiveObjectAndHideOthers ),
                                                  FXMAPFUNC  ( SEL_UPDATE,  PolyDataGroup::ID_SINGLE_ACTIVE_OBJECT, PolyDataGroup::onUpdSingleActiveObject ),
                                                  FXMAPFUNC  ( SEL_CHANGED, PolyDataGroup::ID_SINGLE_ACTIVE_OBJECT, PolyDataGroup::onCmdSingleActiveObject ),
                                                  FXMAPFUNC  ( SEL_COMMAND, PolyDataGroup::ID_SINGLE_ACTIVE_OBJECT, PolyDataGroup::onCmdSingleActiveObject ),
                                                  FXMAPFUNCS ( SEL_UPDATE,  PolyDataGroup::ID_SET_OBJECT_VISIBILITY, PolyDataGroup::ID_SET_OBJECT_VISIBILITY + VTKFOX_NUMBER_OF_MANAGEABLE_OBJECTS - 1, PolyDataGroup::onUpdObjectVisibilityStatus ),
                                                  FXMAPFUNCS ( SEL_COMMAND, PolyDataGroup::ID_SET_OBJECT_VISIBILITY, PolyDataGroup::ID_SET_OBJECT_VISIBILITY + VTKFOX_NUMBER_OF_MANAGEABLE_OBJECTS - 1, PolyDataGroup::onCmdObjectVisibilityStatus ),
                                                  FXMAPFUNC  ( SEL_CHANGED, PolyDataGroup::ID_POLYDATA_MESSAGE,      PolyDataGroup::polyDataMessage ),
                                                  FXMAPFUNC  ( SEL_COMMAND, PolyDataGroup::ID_SAVE_ACTIVE_POLYDATA_TO_VTK_FILE, PolyDataGroup::onCmdSaveActivePolydataToVTKFile ),
                                                  FXMAPFUNC  ( SEL_COMMAND, PolyDataGroup::ID_MAXIMIZE_COLOR_RANGE_ACTIVE_POLYDATA, PolyDataGroup::onCmdMaximizeColorRangeActivePolydata ),
                                                  FXMAPFUNC  ( SEL_COMMAND, PolyDataGroup::ID_ADD_CONST_DATA_TO_ACTIVE_POLYDATA, PolyDataGroup::onCmdAddConstDataToActivePolydata ),
                                                };

FXIMPLEMENT ( PolyDataGroup, FXObject, PolyDataGroupMap, ARRAYNUMBER ( PolyDataGroupMap ) )

void PolyDataGroup::deleteAllMeshData() {
  for ( unsigned int i = 0; i < _polyData.size(); i++ )
    delete _polyData[i];
  _polyData.resize(0);
}

void PolyDataGroup::setSingleActiveObject ( const FXint number ) {
  if (number < 0 || number >= size())
    return;

  _singleActiveObject = number;

  for (int i = 0; i < size(); ++i) {
    if (i != number)
      getPolyData(i)->inactivate();
  }

  getPolyData(number)->activate();
  _mainWindow->handle ( this, FXSEL ( SEL_COMMAND, FXWindow::ID_UPDATE ), NULL );
}

void PolyDataGroup::setSingleActiveObjectAndHideOthers ( const FXint number ) {
  if ( number < 0 || number >= size() )
    return;

  for (int i = 0; i < size(); ++i) {
    if (i != number)
      getPolyData(i)->hide();
  }

  getPolyData(number)->show();
  setSingleActiveObject ( number );
}

void PolyDataGroup::addPolyData ( vtkPolyData* PolyData, const FXString & DataName ) {
  pushBackAndShow(new PolyDataQC(this, aol::getFileName ( DataName.text() ), PolyData));

  if ( PolyData->GetPointData()->GetScalars() )
    getCurrentPolyDataMapper()->SetScalarModeToUsePointData();

  if ( PolyData->GetCellData()->GetScalars() )
    getCurrentPolyDataMapper()->SetScalarModeToUseCellData();
  // note that this data is not displayed if a vtkStripper is used!

  if ( PolyData->GetPointData()->GetScalars() || PolyData->GetCellData()->GetScalars() )
    doUpdateColormap();

}

long PolyDataGroup::onCmdClear ( FXObject* /*sender*/, FXSelector, void* ) {
  _singleActiveObject = -1;
  deleteAllMeshData();
  _mainWindow->setTitle ( "vtkFOX - File: <no file loaded>" );
  _mainWindow->handle ( this, FXSEL ( SEL_COMMAND, FXWindow::ID_UPDATE ), NULL );
  return 1;
}

long PolyDataGroup::onCmdShowAll ( FXObject* /*sender*/, FXSelector, void* ) {
  for (int i = 0; i < size(); ++i)
    getPolyData(i)->show();
  _mainWindow->handle ( this, FXSEL ( SEL_COMMAND, FXWindow::ID_UPDATE ), NULL );
  return 1;
}

long PolyDataGroup::onCmdApplyFullColormapPropertiesToAll ( FXObject* /*sender*/, FXSelector, void* ) {
  for (int i = 0; i < size(); ++i)
    doUpdateColormap ( i );
  _mainWindow->handle ( this, FXSEL ( SEL_COMMAND, FXWindow::ID_UPDATE ), NULL );
  return 1;
}

long PolyDataGroup::onCmdSingleActiveObjectAndHideOthers ( FXObject* /*sender*/, FXSelector /*sel*/, void *ptr ) {
  int number = static_cast<FXint> ( reinterpret_cast<FXival> ( ptr ) );
  setSingleActiveObjectAndHideOthers ( number );
  return 1;
}

long PolyDataGroup::onUpdSingleActiveObject ( FXObject* sender, FXSelector, void* ) {
  FXint value = _singleActiveObject;
  sender->handle ( this, FXSEL ( SEL_COMMAND, FXTextField::ID_SETINTVALUE ), ( void* ) &value );
  return 1;
}

long PolyDataGroup::onCmdSingleActiveObject ( FXObject*, FXSelector, void* ptr ) {
  FXint number = static_cast<FXint> ( reinterpret_cast<FXival> ( ptr ) );
  setSingleActiveObject ( number );
  return 1;
}

long PolyDataGroup::onUpdObjectVisibilityStatus ( FXObject* sender, FXSelector sel, void* ) {
  FXint num = ( ( FXint ) FXSELID ( sel ) ) - ID_SET_OBJECT_VISIBILITY;
  if ( ( num >= size() ) || ( num < 0 ) )
    throw aol::OutOfBoundsException( "Out of bounds message revceived.", __FILE__, __LINE__);

  FXint value = getPolyData ( num )->isVisible();
  sender->handle ( this, FXSEL ( SEL_COMMAND, FXTextField::ID_SETINTVALUE ), ( void* ) &value );
  return 1;
}

long PolyDataGroup::onCmdObjectVisibilityStatus ( FXObject*, FXSelector sel, void* ) {
  FXint num = ( ( FXint ) FXSELID ( sel ) ) - ID_SET_OBJECT_VISIBILITY;
  if ( ( num >= size() ) || ( num < 0 ) )
    throw aol::OutOfBoundsException( "Out of bounds message revceived.", __FILE__, __LINE__);

  if ( getPolyData ( num )->isVisible() )
    getPolyData ( num )->hide();
  else
    getPolyData ( num )->show();
  _mainWindow->handle ( this, FXSEL ( SEL_COMMAND, FXWindow::ID_UPDATE ), NULL );
  return 1;
}

long PolyDataGroup::onCmdSaveActivePolydataToVTKFile ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  vtkPolyData *polyData = getCurrentPolyData ();
  if ( polyData )
    savePolyDataToVTKFile ( polyData, "out.vtk" );

  return 1;
}

long PolyDataGroup::onCmdMaximizeColorRangeActivePolydata ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  vtkPolyData *polyData = getCurrentPolyData ();
  if ( polyData && polyData->GetPointData() && polyData->GetPointData()->GetScalars() ) {
    double range[2];
    polyData->GetPointData()->GetScalars()->GetRange ( range );
    setColorMapRange ( range[0], range[1] );
  }

  return 1;
}

long PolyDataGroup::onCmdAddConstDataToActivePolydata ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {
  vtkPolyData *polyData = getCurrentPolyData ();
  if ( polyData && polyData->GetPointData() && ( polyData->GetPointData()->GetScalars() == NULL ) ) {
    vtkFloatArray *data = vtkFloatArray::New();
    const int numPoints = polyData->GetVerts()->GetSize();
    data->SetNumberOfValues ( numPoints );
    for ( int i = 0; i < numPoints; ++i )
      data->SetValue ( i, 0 );

    polyData->GetPointData()->SetScalars ( data );
    data->Delete();
  }

  return 1;
}

long PolyDataGroup::doUpdateColormap ( const int I ) {
  PolyDataQC *polyDataQC = getPolyData ( ( I < 0 ) ? _singleActiveObject : I );
  if ( polyDataQC && polyDataQC->getMapper() ) {
    vtkColorTransferFunction *colorTF = vtkColorTransferFunction::New();

    cerr << "Switching to colormap " << _colorMap << " in range " << _lowerValueForColormap << " ... " << _upperValueForColormap << endl;

    // These RGB values are just for the constant color case.
    const double redVal = FXREDVAL ( _objectColor ) / 255.0;
    const double greenVal = FXGREENVAL ( _objectColor ) / 255.0;
    const double blueVal = FXBLUEVAL ( _objectColor ) / 255.0;
    initColorTransferFunction<double> ( colorTF, _colorMap, _lowerValueForColormap, _upperValueForColormap, redVal, greenVal, blueVal );
    polyDataQC->getMapper()->SetLookupTable ( colorTF ); // TODO: need to delete old Lookup table??
    polyDataQC->getMapper()->SetColorModeToMapScalars();
    polyDataQC->getMapper()->InterpolateScalarsBeforeMappingOn();

    polyDataQC->getProperties()->SetOpacity ( _opacityValueForColormap );

    _mainWindow->handle ( this, FXSEL ( SEL_COMMAND, FXWindow::ID_UPDATE ), NULL );

    return ( 1 );
  }

  return ( 0 );

}

