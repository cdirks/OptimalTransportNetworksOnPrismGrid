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
* FXVTKCanvas Windows Specific Implementation
* 
* Author: Doug Henry (brilligent@gmail.com)
*/


#include "FXVTKCanvas.h"
#include "vtkFXRenderWindowInteractor.h"

#include <vtkRenderWindow.h>
#include <vtkCommand.h>


void FXVTKCanvas::render()
{
	makeCurrent();
	
	_fxrwi->GetRenderWindow()->SetWindowId(_id);
	//_fxrwi->GetRenderWindow()->SetDisplayId(_display);
	_fxrwi->Render();
	
	makeNonCurrent();
}
