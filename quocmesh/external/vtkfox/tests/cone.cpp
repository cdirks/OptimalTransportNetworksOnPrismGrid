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
#include <vtkConeSource.h>
#include <vtkSuperquadricSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>

#define HAVE_GL_H
#include <fx3d.h>
#include <memory>


class Window : public FXMainWindow
{
		FXDECLARE(Window)
		
	public:
		enum
		{
		    ID_COLOR = FXMainWindow::ID_LAST
	};
	
		Window(FXApp *app)
				: FXMainWindow(app, "vtkFOX Cone Test", NULL, NULL, DECOR_ALL, 50, 50, 640, 480)
				, vis(new FXGLVisual(app, VISUAL_DOUBLEBUFFER))
		{
			FXHorizontalFrame *frame = new FXHorizontalFrame(this, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			FXPacker *canvasFrame = new FXPacker(frame, FRAME_THICK | FRAME_SUNKEN | LAYOUT_FILL_X | LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			FXVTKCanvas *canvas = new FXVTKCanvas(canvasFrame, vis.get(), NULL, 0, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			FXColorWell *colorWell = new FXColorWell(frame, FXRGB(50, 200, 50), this, ID_COLOR);
			
			create_cone_pipeline(canvas->getInteractor());
			
			FXColor color(colorWell->getRGBA());
			coneActor->GetProperty()->SetColor(FXREDVAL(color) / 255.0, FXGREENVAL(color) / 255.0, FXBLUEVAL(color) / 255.0);
		}
		
		long onColor(FXObject *obj, FXSelector sel, void *data)
		{
			FXColor color = (FXColor)(FXuval)data;
			coneActor->GetProperty()->SetColor(FXREDVAL(color) / 255.0, FXGREENVAL(color) / 255.0, FXBLUEVAL(color) / 255.0);
			
			return 1;
		}
		
		void create_cone_pipeline(vtkFXRenderWindowInteractor *rwi)
		{
			vtkRenderer *ren = vtkRenderer::New();
			ren->SetBackground(0.5, 0.5, 0.5);
			
			vtkRenderWindow *renWindow = vtkRenderWindow::New();
			renWindow->AddRenderer(ren);
			
			vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
			rwi->SetInteractorStyle(style);
			
			rwi->SetRenderWindow(renWindow);
			rwi->Initialize();
			
			vtkConeSource *cone = vtkConeSource::New();
			cone->SetResolution(50);
			vtkPolyDataMapper *coneMapper = vtkPolyDataMapper::New();
			coneMapper->SetInput(cone->GetOutput());
			coneActor = vtkActor::New();
			coneActor->SetMapper(coneMapper);
			
			ren->AddActor(coneActor);
			
			ren->Delete();
			renWindow->Delete();
			cone->Delete();
			coneMapper->Delete();
			coneActor->Delete();
		}
		
	protected:
		Window()
		{}
		
	private:
		std::auto_ptr<FXGLVisual> vis;
		vtkActor *coneActor;
};


FXDEFMAP(Window) WindowMap[] =
    {
        FXMAPFUNC(SEL_COMMAND, Window::ID_COLOR, Window::onColor)
    };
    
FXIMPLEMENT(Window, FXMainWindow, WindowMap, ARRAYNUMBER(WindowMap))


int main( int argc, char *argv[] )
{
	FXApp app("cone", "vtkFOX");
	app.init(argc, argv);
	
	Window *win = new Window(&app);
	
	app.create();
	win->show(PLACEMENT_SCREEN);
	
	return app.run();
}
