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

#ifndef __VTKFXGUI_H
#define __VTKFXGUI_H

#include <foxIncludes.h>

class vtkFXViewer;

class vtkFXGui : public FXMainWindow {
public:
  vtkFXGui ( FXApp *app );
  ~vtkFXGui ( );

  FXTabBook *   getControlPanel ( ) {  return _panels;  }
  vtkFXViewer * getViewer ( )       {  return _viewer;  }

private:
  vtkFXGui();

  FXMenuPane*  _filemenu;
  FXMenuPane*  _viewmenu;
  vtkFXViewer* _viewer;
  FXTabBook*   _panels;
};

#endif
