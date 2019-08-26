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

#ifndef __VTKFXSTEREO_H
#define __VTKFXSTEREO_H

#include "FXVTKCanvas.h"
#include <vtkIncludes.h>

// the stereo renderer
// ===================

// only for stereo rendering
class vtkStereoRenderer {
  vtkRenderer *_rendererLeft;
  vtkRenderer *_rendererRight;

public:
  vtkStereoRenderer( )
      : _rendererLeft ( vtkRenderer::New() ),
      _rendererRight ( vtkRenderer::New() ) {
    _rendererRight->SetActiveCamera ( _rendererLeft->GetActiveCamera() );
  }

  ~vtkStereoRenderer( ) {
    _rendererLeft->Delete();
    _rendererRight->Delete();
  }

  void SetBackground ( double R, double G, double B ) {
    _rendererLeft->SetBackground ( R, G, B );
    _rendererRight->SetBackground ( R, G, B );
  }

  void ResetCamera() {
    _rendererLeft->ResetCamera();
    _rendererRight->ResetCamera();
  }

  vtkCamera *GetActiveCamera() {
    return _rendererLeft->GetActiveCamera();
  }

  void RemoveAllViewProps() {
    _rendererLeft->RemoveAllViewProps();
    _rendererRight->RemoveAllViewProps();
  }

  int HasViewProp ( vtkProp *P ) {
    return _rendererLeft->HasViewProp ( P );
  }

  void RemoveViewProp ( vtkProp *P ) {
    _rendererLeft->RemoveViewProp ( P );
    _rendererRight->RemoveViewProp ( P );
  }

  void AddViewProp ( vtkProp *P ) {
    _rendererLeft->AddViewProp ( P );
    _rendererRight->AddViewProp ( P );
  }

  void AutomaticLightCreationOff ( ) {
    _rendererLeft->AutomaticLightCreationOff();
    _rendererRight->AutomaticLightCreationOff();
  }

  void AddLight ( vtkLight *light ){
    _rendererLeft->AddLight( light );
    _rendererRight->AddLight( light );
  }

  vtkRenderer* GetLeftRenderer() {
    return _rendererLeft;
  }

  vtkRenderer* GetRightRenderer() {
    return _rendererRight;
  }
};

// FXVTKCanvasWithSecondCanvasUpdate class
// =======================================

/// only for stereo rendering
class FXVTKCanvasWithSecondCanvasUpdate : public FXVTKCanvas {
  FXDECLARE ( FXVTKCanvasWithSecondCanvasUpdate )

public:
  FXVTKCanvasWithSecondCanvasUpdate ( FXComposite *p, FXGLVisual *vis, FXObject *tgt = NULL, FXSelector sel = 0, FXuint opts = 0, FXint x = 0, FXint y = 0, FXint w = 0, FXint h = 0 );

  ~FXVTKCanvasWithSecondCanvasUpdate() {};

  void create();

  void resize ( FXint w, FXint h ) {
    FXWindow::resize ( w, h );
    _secondVTKCanvasToUpdate->FXWindow::resize ( w, h );
  }

  /*
    void render();
    long onPaint(FXObject *obj, FXSelector sel, void *data);
  */

  long onResize ( FXObject *obj, FXSelector sel, void *data );
  long onLeftButtonDown ( FXObject *obj, FXSelector sel, void *data );
  long onLeftButtonUp ( FXObject *obj, FXSelector sel, void *data );
  long onMiddleButtonDown ( FXObject *obj, FXSelector sel, void *data );
  long onMiddleButtonUp ( FXObject *obj, FXSelector sel, void *data );
  long onRightButtonDown ( FXObject *obj, FXSelector sel, void *data );
  long onRightButtonUp ( FXObject *obj, FXSelector sel, void *data );
  long onMotion ( FXObject *obj, FXSelector sel, void *data );
  //  long onKeyboard(FXObject *obj, FXSelector sel, void *data);

  void setSecondVTKCanvasToUpdate ( FXVTKCanvas *VTKCanvas ) {
    _secondVTKCanvasToUpdate = VTKCanvas;
  }

protected:
  FXVTKCanvasWithSecondCanvasUpdate() {}

  // not implemented (should cause compile error if object is copied)
  FXVTKCanvasWithSecondCanvasUpdate ( const FXVTKCanvasWithSecondCanvasUpdate& FVC );
  FXVTKCanvasWithSecondCanvasUpdate& operator= ( const FXVTKCanvasWithSecondCanvasUpdate& FVC );

private:
  FXVTKCanvas *_secondVTKCanvasToUpdate;

};

#endif
