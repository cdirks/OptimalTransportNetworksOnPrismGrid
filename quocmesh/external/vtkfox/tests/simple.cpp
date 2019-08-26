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

#include "vtkFXRenderWindowInteractor.h"
#include "FXVTKCanvas.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#define HAVE_GL_H
#include <fx3d.h>
#include <memory>


class Window : public FXMainWindow
{
		FXDECLARE(Window)
		
	public:
		Window(FXApp *app)
				: FXMainWindow(app, "vtkFOX Simple Test", NULL, NULL, DECOR_ALL, 50, 50, 640, 480)
				, vis(new FXGLVisual(app, VISUAL_DOUBLEBUFFER))
		{
			FXHorizontalFrame *frame = new FXHorizontalFrame(this, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			FXVTKCanvas *canvas = new FXVTKCanvas(frame, vis.get(), NULL, 0, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			
			//
			// do vtk stuff
			//
			
			vtkRenderer *ren = vtkRenderer::New();
			ren->SetBackground(0.5, 0.5, 0.5);
			
			vtkRenderWindow *renWindow = vtkRenderWindow::New();
			renWindow->AddRenderer(ren);
			
			canvas->getInteractor()->SetRenderWindow(renWindow);
			canvas->getInteractor()->Initialize();
		}
		
	protected:
		Window()
		{}
		
	private:
		std::auto_ptr<FXGLVisual> vis;
};


FXIMPLEMENT(Window, FXMainWindow, NULL, 0)


int main( int argc, char *argv[] )
{
	FXApp app("simple", "vtkFOX");
	app.init(argc, argv);
	
	Window *win = new Window(&app);
	
	app.create();
	win->show(PLACEMENT_SCREEN);
	
	return app.run();
}
