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
* vtkFXRenderWindowInteractor Definition
* 
* This class is used to "close the loop" with VTK.  The FXVTKCanvas widget
* handles the fox side and translates events to VTK.  This class uses
* the canvas to update its own information.
* 
* Author: Doug Henry (brilligent@gmail.com)
*/


#ifndef vtkFXRenderWindowInteractor_DEF
#define vtkFXRenderWindowInteractor_DEF


#include <vtkRenderWindowInteractor.h>

class FXVTKCanvas;


/**
 * Interactor to allow VTK to be embedded in a FOX application.  This
 * class contains a reference to a FXVTKCanvas widget.
 */

class vtkFXRenderWindowInteractor : public vtkRenderWindowInteractor
{
	public:
		vtkFXRenderWindowInteractor();
		vtkFXRenderWindowInteractor(FXVTKCanvas *canvas);
		~vtkFXRenderWindowInteractor();
		
		void setCanvas(FXVTKCanvas *canvas);
		
		void Initialize();
		void Enable();
		void Disable();
		void Start();
		void SetRenderWindow(vtkRenderWindow *renderer);
		void UpdateSize(int w, int h);
		void TerminateApp();
		
		virtual char const* GetClassNameW() const { return "vtkFXRenderWindowInteractor"; }
		
	protected:
		// not implemented (should cause compile error)
		vtkFXRenderWindowInteractor(const vtkFXRenderWindowInteractor& RWI);
		vtkFXRenderWindowInteractor& operator=(const vtkFXRenderWindowInteractor& RWI);
		
	private:
		FXVTKCanvas *_canvas;
};


#endif
