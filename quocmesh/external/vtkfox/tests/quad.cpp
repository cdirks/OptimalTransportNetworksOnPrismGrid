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
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkQuadric.h>
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkOutlineFilter.h>
#include <vtkImageData.h>
#include <vtkInteractorStyleTrackballCamera.h>

#define HAVE_GL_H
#include <fx3d.h>
#include <memory>


class Window : public FXMainWindow
{
		FXDECLARE(Window)
		
	public:
		enum {
		    ID_NUM_CONTOUR = FXMainWindow::ID_LAST,
		    ID_RESOLUTION,
		    ID_LAST
	};
	
		Window(FXApp *app)
				: FXMainWindow(app, "vtkFOX Quadric Test", NULL, NULL, DECOR_ALL, 50, 50, 640, 480)
				, vis(new FXGLVisual(app, VISUAL_BEST|VISUAL_DOUBLEBUFFER))
		{
			FXHorizontalFrame *frame = new FXHorizontalFrame(this, LAYOUT_FILL_X | LAYOUT_FILL_Y, 0,0,0,0, 10,10,10,10, 10,10);
			FXPacker *canvasFrame = new FXPacker(frame, FRAME_THICK | FRAME_SUNKEN | LAYOUT_FILL_X | LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			canvas = new FXVTKCanvas(canvasFrame, vis.get(), NULL, 0, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			
			FXMatrix *toolframe = new FXMatrix(frame, 2, MATRIX_BY_COLUMNS | FRAME_GROOVE | LAYOUT_FILL_Y, 0,0,0,0, 5,5,5,5, 1,10);
			new FXLabel(toolframe, "Contours:", NULL, LABEL_NORMAL|LAYOUT_CENTER_Y);
			_numContours = new FXSpinner(toolframe, 3, this, ID_NUM_CONTOUR, SPIN_NORMAL|FRAME_SUNKEN|FRAME_THICK);
			new FXLabel(toolframe, "Resolution:", NULL, LABEL_NORMAL|LAYOUT_CENTER_Y);
			_resolution = new FXSpinner(toolframe, 3, this, ID_RESOLUTION, SPIN_NORMAL|FRAME_SUNKEN|FRAME_THICK);
			
			_numContours->setRange(3, 15);
			_numContours->setValue(5);
			
			_resolution->setRange(3, 50);
			_resolution->setValue(15);
			
			create_quad_pipeline(canvas->getInteractor());
		}
		
		long onNumContour(FXObject *obj, FXSelector sel, void *data)
		{
			contours->GenerateValues(_numContours->getValue(), 0.0, 1.2);
			
			canvas->render();
			
			return 1;
		}
		
		long onResolution(FXObject *obj, FXSelector sel, void *data)
		{
			sample->SetSampleDimensions(_resolution->getValue(), _resolution->getValue(), _resolution->getValue());
			
			canvas->render();
			
			return 1;
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
			sample = vtkSampleFunction::New();
			sample->SetSampleDimensions(_resolution->getValue(), _resolution->getValue(), _resolution->getValue());
			sample->SetImplicitFunction(quadric);
			
			// Create five surfaces F(x,y,z) = constant between range specified
			contours = vtkContourFilter::New();
			contours->SetInput(sample->GetOutput());
			contours->GenerateValues(_numContours->getValue(), 0.0, 1.2);
			
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
		FXVTKCanvas *canvas;
		FXSpinner *_numContours;
		FXSpinner *_resolution;
		
		vtkContourFilter *contours;
		vtkSampleFunction *sample;
};


FXDEFMAP(Window) WindowMap[] =
    {
        FXMAPFUNC(SEL_COMMAND, Window::ID_NUM_CONTOUR, Window::onNumContour),
        FXMAPFUNC(SEL_COMMAND, Window::ID_RESOLUTION, Window::onResolution)
    };
    
FXIMPLEMENT(Window, FXMainWindow, WindowMap, ARRAYNUMBER(WindowMap))


int main( int argc, char *argv[] )
{
	FXApp app("quad", "vtkFOX");
	app.init(argc, argv);
	
	Window *win = new Window(&app);
	
	app.create();
	win->show(PLACEMENT_SCREEN);
	
	return app.run();
}
