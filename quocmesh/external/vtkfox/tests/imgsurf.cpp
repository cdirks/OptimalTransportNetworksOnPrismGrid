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

#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkImageData.h>
#include <vtkImageActor.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkPolyData.h>
#include <vtkWarpScalar.h>
#include <vtkInteractorStyleTrackballCamera.h>

#define HAVE_GL_H
#include <fx3d.h>
#include <memory>


class Window : public FXMainWindow
{
		FXDECLARE(Window)
		
	public:
		enum {
		    ID_LAST = FXMainWindow::ID_LAST
	};
	
		Window(FXApp *app)
				: FXMainWindow(app, "vtkFOX Image Surface Test", NULL, NULL, DECOR_ALL, 50, 50, 640, 480)
				, vis(new FXGLVisual(app, VISUAL_BEST|VISUAL_DOUBLEBUFFER))
		{
			FXHorizontalFrame *frame = new FXHorizontalFrame(this, LAYOUT_FILL_X | LAYOUT_FILL_Y, 0,0,0,0, 10,10,10,10, 10,10);
			FXPacker *canvasFrame = new FXPacker(frame, FRAME_THICK | FRAME_SUNKEN | LAYOUT_FILL_X | LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			canvas = new FXVTKCanvas(canvasFrame, vis.get(), NULL, 0, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			
			create_data();
			create_pipeline(canvas->getInteractor());
		}
		
		~Window()
		{
			data->Delete();
		}
		
		void create_data()
		{
			_size = 100;
			float inc = 30 / float(_size);
			
			data = vtkImageData::New();
			data->SetDimensions(_size, _size, 1);
			data->SetScalarTypeToFloat();
			data->SetNumberOfScalarComponents(1);
			data->AllocateScalars();
			
			float *ptr = (float*)data->GetScalarPointer();
			float value;
			
			for (float i = -15; i <= 15; i += inc)
				for (float j = -15; j <= 15; j += inc)
				{
					value = (sin(i) / float(i)) + (sin(j) / float(j));
					
					_min = FXMIN(_min, value);
					_max = FXMAX(_max, value);
					
					*ptr++ = value;
				}
		}
		
		void create_pipeline(vtkFXRenderWindowInteractor *rwi)
		{
			vtkRenderer *ren = vtkRenderer::New();
			ren->SetBackground(0.5, 0.5, 0.5);
			
			vtkRenderWindow *renWindow = vtkRenderWindow::New();
			renWindow->AddRenderer(ren);
			
			vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
			rwi->SetInteractorStyle(style);
			
			rwi->SetRenderWindow(renWindow);
			rwi->Initialize();
			
			vtkImageDataGeometryFilter *geom = vtkImageDataGeometryFilter::New();
			geom->SetInput(data);
			
			vtkWarpScalar *warp = vtkWarpScalar::New();
			warp->SetInput(geom->GetOutput());
			warp->SetScaleFactor(_size / 10.0);
			
			vtkDataSetMapper *mapper = vtkDataSetMapper::New();
			mapper->SetInput(warp->GetPolyDataOutput());
			mapper->SetScalarRange(_min, _max);
			
			vtkActor *imageActor = vtkActor::New();
			imageActor->SetMapper(mapper);
			
			ren->AddActor(imageActor);
			
			mapper->Delete();
			geom->Delete();
			imageActor->Delete();
			ren->Delete();
			renWindow->Delete();
		}
		
	protected:
		Window()
		{}
		
	private:
		std::auto_ptr<FXGLVisual> vis;
		FXVTKCanvas *canvas;
		vtkImageData *data;
		int _size;
		float _max;
		float _min;
};


FXIMPLEMENT(Window, FXMainWindow, NULL, 0)


int main( int argc, char *argv[] )
{
	FXApp app("imgsurf", "vtkFOX");
	app.init(argc, argv);
	
	Window *win = new Window(&app);
	
	app.create();
	win->show(PLACEMENT_SCREEN);
	
	return app.run();
}
