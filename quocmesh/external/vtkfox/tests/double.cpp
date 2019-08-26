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
#include <vtkQuadric.h>
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkOutlineFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkImageData.h>
#include <vtkInteractorStyleTrackballCamera.h>

#define HAVE_GL_H
#include <fx3d.h>
#include <memory>


class Window : public FXMainWindow
{
		FXDECLARE(Window)
		
	public:
		Window(FXApp *app)
				: FXMainWindow(app, "vtkFOX Double Test", NULL, NULL, DECOR_ALL, 50, 50, 640, 300)
				, vis(new FXGLVisual(app, VISUAL_DOUBLEBUFFER))
		{
			FXHorizontalFrame *frame = new FXHorizontalFrame(this, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			
			FXPacker *coneCanvasFrame = new FXPacker(frame, FRAME_THICK | FRAME_SUNKEN | LAYOUT_FILL_X | LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			FXVTKCanvas *coneCanvas = new FXVTKCanvas(coneCanvasFrame, vis.get(), NULL, 0, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			create_cone_pipeline(coneCanvas->getInteractor());
			
			FXPacker *quadCanvasFrame = new FXPacker(frame, FRAME_THICK | FRAME_SUNKEN | LAYOUT_FILL_X | LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			FXVTKCanvas *quadCanvas = new FXVTKCanvas(quadCanvasFrame, vis.get(), NULL, 0, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			create_quad_pipeline(quadCanvas->getInteractor());
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
			vtkActor *coneActor = vtkActor::New();
			coneActor->SetMapper(coneMapper);
			coneActor->GetProperty()->SetColor(0.33, 1.0, 0.33);
			ren->AddActor(coneActor);
			
			ren->Delete();
			renWindow->Delete();
			cone->Delete();
			coneMapper->Delete();
			coneActor->Delete();
		}
		
		void create_quad_pipeline(vtkFXRenderWindowInteractor *rwi)
		{
			vtkRenderer *ren = vtkRenderer::New();
			ren->SetBackground(0.5, 0.5, 0.5);
			
			vtkRenderWindow *renWindow = vtkRenderWindow::New();
			renWindow->AddRenderer(ren);
			
			vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
			rwi->SetInteractorStyle(style);
			
			rwi->SetRenderWindow(renWindow);
			rwi->Initialize();
			
			// -- create the quadric function object --
			
			// create the quadric function definition
			vtkQuadric *quadric = vtkQuadric::New();
			quadric->SetCoefficients(.5, 1, .2, 0, .1, 0, 0, .2, 0, 0);
			
			// sample the quadric function
			vtkSampleFunction *sample = vtkSampleFunction::New();
			sample->SetSampleDimensions(20, 20, 20);
			sample->SetImplicitFunction(quadric);
			
			// Create five surfaces F(x,y,z) = constant between range specified
			vtkContourFilter *contours = vtkContourFilter::New();
			contours->SetInput(sample->GetOutput());
			contours->GenerateValues(5, 0.0, 1.2);
			
			// map the contours to graphical primitives
			vtkPolyDataMapper *contMapper = vtkPolyDataMapper::New();
			contMapper->SetInput(contours->GetOutput());
			contMapper->SetScalarRange(0.0, 1.2);
			
			// create an actor for the contours
			vtkActor *contActor = vtkActor::New();
			contActor->SetMapper(contMapper);
			
			// -- create a box around the function to indicate the sampling volume --
			
			// create outline
			vtkOutlineFilter *outline = vtkOutlineFilter::New();
			outline->SetInput(sample->GetOutput());
			
			// map it to graphics primitives
			vtkPolyDataMapper *outlineMapper = vtkPolyDataMapper::New();
			outlineMapper->SetInput(outline->GetOutput());
			
			// create an actor for it
			vtkActor *outlineActor = vtkActor::New();
			outlineActor->SetMapper(outlineMapper);
			outlineActor->GetProperty()->SetColor(0, 0, 0);
			
			// -- render both of the objects --
			
			// add the actors to the scene
			ren->AddActor(contActor);
			ren->AddActor(outlineActor);
			ren->SetBackground(1, 1, 1); // Background color white
			
			quadric->Delete();
			sample->Delete();
			contours->Delete();
			contMapper->Delete();
			contActor->Delete();
			outline->Delete();
			outlineMapper->Delete();
			outlineActor->Delete();
			ren->Delete();
			renWindow->Delete();
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
	FXApp app("double", "vtkFOX");
	app.init(argc, argv);
	
	Window *win = new Window(&app);
	
	app.create();
	win->show(PLACEMENT_SCREEN);
	
	return app.run();
}
