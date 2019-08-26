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
#include <vtkXYPlotActor.h>
#include <vtkProperty2D.h>
#include <vtkFloatArray.h>
#include <vtkFieldData.h>
#include <vtkDataObject.h>
#include <vtkTextProperty.h>

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
				: FXMainWindow(app, "vtkFOX XY-Plot Test", NULL, NULL, DECOR_ALL, 50, 50, 640, 480)
				, vis(new FXGLVisual(app, VISUAL_DOUBLEBUFFER))
		{
			FXHorizontalFrame *frame = new FXHorizontalFrame(this, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			FXPacker *canvasFrame = new FXPacker(frame, FRAME_THICK | FRAME_SUNKEN | LAYOUT_FILL_X | LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			FXVTKCanvas *canvas = new FXVTKCanvas(canvasFrame, vis.get(), NULL, 0, LAYOUT_FILL_X | LAYOUT_FILL_Y);
			
			create_xyplot_pipeline(canvas->getInteractor());
		}
		
		void create_xyplot_pipeline(vtkFXRenderWindowInteractor *rwi)
		{
			vtkRenderer *ren = vtkRenderer::New();
			ren->SetBackground(1, 1, 1);
			
			vtkRenderWindow *renWindow = vtkRenderWindow::New();
			renWindow->AddRenderer(ren);
			
			//vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
			//rwi->SetInteractorStyle(style);
			
			rwi->SetRenderWindow(renWindow);
			rwi->Initialize();
			
			vtkFloatArray *xdata = vtkFloatArray::New();
			vtkFloatArray *ydata = vtkFloatArray::New();
			
			for (float x = -15.0; x <= 15.0; x += 0.3)
			{
				xdata->InsertNextValue(x);
				ydata->InsertNextValue(sin(x)/x);
			}
			
			vtkFieldData *data = vtkFieldData::New();
			data->AddArray(xdata);
			data->AddArray(ydata);
			
			vtkDataObject *dataobj = vtkDataObject::New();
			dataobj->SetFieldData(data);
			
			vtkXYPlotActor *plot = vtkXYPlotActor::New();
			plot->AddDataObjectInput(dataobj);
			
			//plot->SetTitle("Sample XY-Plot");
			plot->SetXTitle("x");
			plot->SetYTitle("sin(x)");
			plot->SetXValuesToValue();
			plot->SetDataObjectXComponent(0, 0);
			plot->SetDataObjectYComponent(0, 1);
			//plot->SetXRange(0, 2.0*M_PI);
			//plot->SetYRange( -1, 1);
			plot->PlotPointsOn();
			plot->PlotLinesOn();
			//plot->LegendOn();
			
			plot->GetProperty()->SetColor(0, 0, 0);
			plot->GetProperty()->SetPointSize(4);
			plot->GetAxisLabelTextProperty()->SetFontSize(2);
			plot->GetPositionCoordinate()->SetValue(0.03, 0.03, 0.0);
			plot->GetPosition2Coordinate()->SetValue(0.96, 0.96, 0.0);
			
			ren->AddActor2D(plot);
			
			xdata->Delete();
			ydata->Delete();
			data->Delete();
			dataobj->Delete();
			plot->Delete();
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
	FXApp app("xyplot", "vtkFOX");
	app.init(argc, argv);
	
	Window *win = new Window(&app);
	
	app.create();
	win->show(PLACEMENT_SCREEN);
	
	return app.run();
}
