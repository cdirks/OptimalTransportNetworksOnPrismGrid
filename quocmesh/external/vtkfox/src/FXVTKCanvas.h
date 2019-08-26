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

/*
* FXVTKCanvas Widget Definition
* 
* This class is intended to be a drop-in replacement for FXGLCanvas, with
* the exception that it allows for VTK rendering.  This class holds an
* instance of vtkFXRenderWindowInteractor, mapping all appropriate functions
* and events, to allow for proper VTK use.
* 
* Author: Doug Henry (brilligent@gmail.com)
*/


#ifndef FXVTKCanvas_DEF
#define FXVTKCanvas_DEF


#include <fx.h>
#include <FXGLCanvas.h>

class vtkFXRenderWindowInteractor;


class FXVTKCanvas : public FXGLCanvas
{
		FXDECLARE(FXVTKCanvas)
		
	public:
		FXVTKCanvas(FXComposite *p, FXGLVisual *vis, FXObject *tgt = NULL, FXSelector sel = 0, FXuint opts = 0, FXint x = 0, FXint y = 0, FXint w = 0, FXint h = 0);
		
		~FXVTKCanvas();
		
		vtkFXRenderWindowInteractor* getInteractor() const;
		void setInteractor(vtkFXRenderWindowInteractor *fxrwi);
		
		void create();
		void render();
		
		long onPaint(FXObject *obj, FXSelector sel, void *data);
		long onResize(FXObject *obj, FXSelector sel, void *data);
		
		long onLeftButtonDown(FXObject *obj, FXSelector sel, void *data);
		long onLeftButtonUp(FXObject *obj, FXSelector sel, void *data);
		long onMiddleButtonDown(FXObject *obj, FXSelector sel, void *data);
		long onMiddleButtonUp(FXObject *obj, FXSelector sel, void *data);
		long onRightButtonDown(FXObject *obj, FXSelector sel, void *data);
		long onRightButtonUp(FXObject *obj, FXSelector sel, void *data);
		long onMotion(FXObject *obj, FXSelector sel, void *data);
		long onKeyboard(FXObject *obj, FXSelector sel, void *data);
		
	protected:
		FXVTKCanvas()
		{}
		
		// not implemented (should cause compile error if object is copied)
		FXVTKCanvas(const FXVTKCanvas& FVC);
		FXVTKCanvas& operator=(const FXVTKCanvas& FVC);
		
	private:
		vtkFXRenderWindowInteractor *_fxrwi;
		void *_id;
		void *_display;
};


#endif
