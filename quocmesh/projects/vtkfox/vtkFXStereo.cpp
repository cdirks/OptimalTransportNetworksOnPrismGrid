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

// foxIncludes.h is needed to have <complex> included before aol.
#include <foxIncludes.h>

#include "vtkFXStereo.h"

#include <aol.h>
#include "vtkFXRenderWindowInteractor.h"

// FXVTKCanvasWithSecondCanvasUpdate implementation
// ================================================

FXDEFMAP ( FXVTKCanvasWithSecondCanvasUpdate ) FXVTKCanvasWithSecondCanvasUpdateMap[] =
  {
    FXMAPFUNC ( SEL_PAINT, 0, FXVTKCanvas::onPaint ),
    FXMAPFUNC ( SEL_CONFIGURE, 0, FXVTKCanvasWithSecondCanvasUpdate::onResize ),
    FXMAPFUNC ( SEL_LEFTBUTTONPRESS, 0, FXVTKCanvasWithSecondCanvasUpdate::onLeftButtonDown ),
    FXMAPFUNC ( SEL_LEFTBUTTONRELEASE, 0, FXVTKCanvasWithSecondCanvasUpdate::onLeftButtonUp ),
    FXMAPFUNC ( SEL_MIDDLEBUTTONPRESS, 0, FXVTKCanvasWithSecondCanvasUpdate::onMiddleButtonDown ),
    FXMAPFUNC ( SEL_MIDDLEBUTTONRELEASE, 0, FXVTKCanvasWithSecondCanvasUpdate::onMiddleButtonUp ),
    FXMAPFUNC ( SEL_RIGHTBUTTONPRESS, 0, FXVTKCanvasWithSecondCanvasUpdate::onRightButtonDown ),
    FXMAPFUNC ( SEL_RIGHTBUTTONRELEASE, 0, FXVTKCanvasWithSecondCanvasUpdate::onRightButtonUp ),
    FXMAPFUNC ( SEL_MOTION, 0, FXVTKCanvasWithSecondCanvasUpdate::onMotion ),
    FXMAPFUNC ( SEL_KEYPRESS, 0, FXVTKCanvas::onKeyboard )
  };

FXIMPLEMENT ( FXVTKCanvasWithSecondCanvasUpdate, FXVTKCanvas, FXVTKCanvasWithSecondCanvasUpdateMap, ARRAYNUMBER ( FXVTKCanvasWithSecondCanvasUpdateMap ) )


FXVTKCanvasWithSecondCanvasUpdate::FXVTKCanvasWithSecondCanvasUpdate ( FXComposite *p, FXGLVisual *vis, FXObject *tgt, FXSelector sel, FXuint opts, FXint x, FXint y, FXint w, FXint h )
    : FXVTKCanvas ( p, vis, tgt, sel, opts, x, y, w, h ),
  _secondVTKCanvasToUpdate ( NULL ) {}

void FXVTKCanvasWithSecondCanvasUpdate::create() {
  if ( _secondVTKCanvasToUpdate == NULL )
    throw aol::Exception ( "_secondVTKCanvasToUpdate not set", __FILE__, __LINE__ );
  FXVTKCanvas::create();
}

long FXVTKCanvasWithSecondCanvasUpdate::onResize ( FXObject*, FXSelector, void* ) {
  getInteractor()->UpdateSize ( getWidth(), getHeight() );
  _secondVTKCanvasToUpdate->getInteractor()->UpdateSize ( getWidth(), getHeight() );

  return 1;
}

long FXVTKCanvasWithSecondCanvasUpdate::onLeftButtonDown ( FXObject*, FXSelector, void *data ) {
  FXEvent *event = ( FXEvent* ) data;

  grab();

  getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  getInteractor()->InvokeEvent ( vtkCommand::LeftButtonPressEvent, NULL );

  _secondVTKCanvasToUpdate->getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  _secondVTKCanvasToUpdate->getInteractor()->InvokeEvent ( vtkCommand::LeftButtonPressEvent, NULL );

  setFocus();

  return 1;
}

long FXVTKCanvasWithSecondCanvasUpdate::onLeftButtonUp ( FXObject*, FXSelector, void *data ) {
  FXEvent *event = ( FXEvent* ) data;

  getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  getInteractor()->InvokeEvent ( vtkCommand::LeftButtonReleaseEvent, NULL );

  _secondVTKCanvasToUpdate->getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  _secondVTKCanvasToUpdate->getInteractor()->InvokeEvent ( vtkCommand::LeftButtonReleaseEvent, NULL );

  ungrab();

  return 1;
}

long FXVTKCanvasWithSecondCanvasUpdate::onMiddleButtonDown ( FXObject*, FXSelector, void *data ) {
  FXEvent *event = ( FXEvent* ) data;

  grab();

  getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  getInteractor()->InvokeEvent ( vtkCommand::MiddleButtonPressEvent, NULL );

  _secondVTKCanvasToUpdate->getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  _secondVTKCanvasToUpdate->getInteractor()->InvokeEvent ( vtkCommand::MiddleButtonPressEvent, NULL );

  setFocus();

  return 1;
}

long FXVTKCanvasWithSecondCanvasUpdate::onMiddleButtonUp ( FXObject*, FXSelector, void *data ) {
  FXEvent *event = ( FXEvent* ) data;

  getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  getInteractor()->InvokeEvent ( vtkCommand::MiddleButtonReleaseEvent, NULL );

  _secondVTKCanvasToUpdate->getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  _secondVTKCanvasToUpdate->getInteractor()->InvokeEvent ( vtkCommand::MiddleButtonReleaseEvent, NULL );

  ungrab();

  return 1;
}

long FXVTKCanvasWithSecondCanvasUpdate::onRightButtonDown ( FXObject*, FXSelector, void *data ) {
  FXEvent *event = ( FXEvent* ) data;

  grab();

  getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  getInteractor()->InvokeEvent ( vtkCommand::RightButtonPressEvent, NULL );

  _secondVTKCanvasToUpdate->getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  _secondVTKCanvasToUpdate->getInteractor()->InvokeEvent ( vtkCommand::RightButtonPressEvent, NULL );

  setFocus();

  return 1;
}

long FXVTKCanvasWithSecondCanvasUpdate::onRightButtonUp ( FXObject*, FXSelector, void *data ) {
  FXEvent *event = ( FXEvent* ) data;

  getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  getInteractor()->InvokeEvent ( vtkCommand::RightButtonReleaseEvent, NULL );

  _secondVTKCanvasToUpdate->getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  _secondVTKCanvasToUpdate->getInteractor()->InvokeEvent ( vtkCommand::RightButtonReleaseEvent, NULL );

  ungrab();

  return 1;
}

long FXVTKCanvasWithSecondCanvasUpdate::onMotion ( FXObject*, FXSelector, void *data ) {
  FXEvent *event = ( FXEvent* ) data;

  getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  getInteractor()->InvokeEvent ( vtkCommand::MouseMoveEvent, NULL );

  _secondVTKCanvasToUpdate->getInteractor()->SetEventInformationFlipY ( event->win_x, event->win_y );
  _secondVTKCanvasToUpdate->getInteractor()->InvokeEvent ( vtkCommand::MouseMoveEvent, NULL );
  return 1;
}
