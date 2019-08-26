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

#include <foxIncludes.h>

#include <aol.h>

typedef float RType;

// Main Window
class ScribbleWindow : public FXMainWindow {

  // Macro for class hierarchy declarations
  FXDECLARE(ScribbleWindow)

private:

  FXHorizontalFrame *contents;                // Content frame
  FXVerticalFrame   *canvasFrame;             // Canvas frame
  FXVerticalFrame   *buttonFrame;             // Button frame
  FXCanvas          *canvas;                  // Canvas to draw into
  int                mdflag;                  // Mouse button down?
  int                dirty;                   // Canvas has been painted?
  FXColor            drawColor;               // Color for the line
  FXDataTarget       color_target;
protected:
  ScribbleWindow(){}

public:

  // Message handlers
  long onPaint(FXObject*,FXSelector,void*);
  long onMouseDown(FXObject*,FXSelector,void*);
  long onMouseUp(FXObject*,FXSelector,void*);
  long onMouseMove(FXObject*,FXSelector,void*);
  long onCmdClear(FXObject*,FXSelector,void*);
  long onUpdClear(FXObject*,FXSelector,void*);
  long onCmdColor(FXObject*,FXSelector,void*);
  long onUpdColor(FXObject*,FXSelector,void*);

public:

  // Messages for our class
  enum{
    ID_CANVAS=FXMainWindow::ID_LAST,
    ID_CLEAR,
    ID_COLOR,
    ID_LAST
    };

public:

  // ScribbleWindow's constructor
  ScribbleWindow(FXApp* a);

  // Initialize
  virtual void create();

  virtual ~ScribbleWindow();
  };



// Message Map for the Scribble Window class
FXDEFMAP(ScribbleWindow) ScribbleWindowMap[]={

  //________Message_Type_____________________ID____________Message_Handler_______
  FXMAPFUNC(SEL_PAINT,             ScribbleWindow::ID_CANVAS, ScribbleWindow::onPaint),
  FXMAPFUNC(SEL_LEFTBUTTONPRESS,   ScribbleWindow::ID_CANVAS, ScribbleWindow::onMouseDown),
  FXMAPFUNC(SEL_LEFTBUTTONRELEASE, ScribbleWindow::ID_CANVAS, ScribbleWindow::onMouseUp),
  FXMAPFUNC(SEL_MOTION,            ScribbleWindow::ID_CANVAS, ScribbleWindow::onMouseMove),
  FXMAPFUNC(SEL_COMMAND,           ScribbleWindow::ID_CLEAR,  ScribbleWindow::onCmdClear),
  FXMAPFUNC(SEL_UPDATE,            ScribbleWindow::ID_CLEAR,  ScribbleWindow::onUpdClear),
  FXMAPFUNC(SEL_COMMAND,           ScribbleWindow::ID_COLOR,  ScribbleWindow::onCmdColor),
  FXMAPFUNC(SEL_UPDATE,            ScribbleWindow::ID_COLOR,  ScribbleWindow::onUpdColor),
  };



// Macro for the ScribbleApp class hierarchy implementation
FXIMPLEMENT(ScribbleWindow,FXMainWindow,ScribbleWindowMap,ARRAYNUMBER(ScribbleWindowMap))



// Construct a ScribbleWindow
ScribbleWindow::ScribbleWindow(FXApp *a):FXMainWindow(a,"Scribble Application",NULL,NULL,DECOR_ALL,0,0,800,600){

  contents=new FXHorizontalFrame(this,LAYOUT_SIDE_TOP|LAYOUT_FILL_X|LAYOUT_FILL_Y,0,0,0,0, 0,0,0,0);

  // LEFT pane to contain the canvas
  canvasFrame=new FXVerticalFrame(contents,FRAME_SUNKEN|LAYOUT_FILL_X|LAYOUT_FILL_Y|LAYOUT_TOP|LAYOUT_LEFT,0,0,0,0,10,10,10,10);

    // Label above the canvas
    new FXLabel(canvasFrame,"Canvas Frame",NULL,JUSTIFY_CENTER_X|LAYOUT_FILL_X);

    // Horizontal divider line
    new FXHorizontalSeparator(canvasFrame,SEPARATOR_GROOVE|LAYOUT_FILL_X);


    // Drawing canvas
    canvas=new FXCanvas(canvasFrame,this,ID_CANVAS,FRAME_SUNKEN|FRAME_THICK|LAYOUT_FILL_X|LAYOUT_FILL_Y|LAYOUT_FILL_ROW|LAYOUT_FILL_COLUMN);

  // RIGHT pane for the buttons
  buttonFrame=new FXVerticalFrame(contents,FRAME_SUNKEN|LAYOUT_FILL_Y|LAYOUT_TOP|LAYOUT_LEFT,0,0,0,0,10,10,10,10);

    // Label above the buttons
    new FXLabel(buttonFrame,"Button Frame",NULL,JUSTIFY_CENTER_X|LAYOUT_FILL_X);

    // Horizontal divider line
    new FXHorizontalSeparator(buttonFrame,SEPARATOR_RIDGE|LAYOUT_FILL_X);

    // Button to clear
    new FXButton(buttonFrame,"&Clear",NULL,this,ID_CLEAR,FRAME_THICK|FRAME_RAISED|LAYOUT_FILL_X|LAYOUT_TOP|LAYOUT_LEFT,0,0,0,0,10,10,5,5);

    // Exit button
    new FXButton(buttonFrame,"&Exit",NULL,getApp(),FXApp::ID_QUIT,FRAME_THICK|FRAME_RAISED|LAYOUT_FILL_X|LAYOUT_TOP|LAYOUT_LEFT,0,0,0,0,10,10,5,5);

  // Initialize private variables
  drawColor=FXRGB(255,0,0);
  mdflag=0;
  dirty=0;

    color_target.connect(drawColor);

    new FXColorWell(buttonFrame,0,&color_target,FXDataTarget::ID_VALUE);
    new FXColorWell(buttonFrame,0,&color_target,FXDataTarget::ID_VALUE);
    new FXColorWell(buttonFrame, drawColor, this, ID_COLOR);

  }


ScribbleWindow::~ScribbleWindow(){
  }


// Create and initialize
void ScribbleWindow::create(){

  // Create the windows
  FXMainWindow::create();

  // Make the main window appear
  show(PLACEMENT_SCREEN);
  }



// Mouse button was pressed somewhere
long ScribbleWindow::onMouseDown(FXObject*,FXSelector,void*){
  canvas->grab();

  // While the mouse is down, we'll draw lines
  mdflag=1;

  return 1;
  }



// The mouse has moved, draw a line
long ScribbleWindow::onMouseMove(FXObject*, FXSelector, void* ptr){
  FXEvent *ev=(FXEvent*)ptr;

  // Draw
  if(mdflag){

    // Get DC for the canvas
    FXDCWindow dc(canvas);

    // Set foreground color
    dc.setForeground(drawColor);

    // Draw line
    dc.drawLine(ev->last_x, ev->last_y, ev->win_x, ev->win_y);

    // We have drawn something, so now the canvas is dirty
    dirty=1;
    }
  return 1;
  }



// The mouse button was released again
long ScribbleWindow::onMouseUp(FXObject*,FXSelector,void* ptr){
  FXEvent *ev=(FXEvent*) ptr;
  canvas->ungrab();
  if(mdflag){
    FXDCWindow dc(canvas);

    dc.setForeground(drawColor);
    dc.drawLine(ev->last_x, ev->last_y, ev->win_x, ev->win_y);

    // We have drawn something, so now the canvas is dirty
    dirty=1;

    // Mouse no longer down
    mdflag=0;
    }
  return 1;
  }


// Paint the canvas
long ScribbleWindow::onPaint(FXObject*,FXSelector,void* ptr){
  FXEvent *ev=(FXEvent*)ptr;
  FXDCWindow dc(canvas,ev);
  dc.setForeground(canvas->getBackColor());
  dc.fillRectangle(ev->rect.x,ev->rect.y,ev->rect.w,ev->rect.h);
  return 1;
  }


// Handle the clear message
long ScribbleWindow::onCmdClear(FXObject*,FXSelector,void*){
  FXDCWindow dc(canvas);
  dc.setForeground(canvas->getBackColor());
  dc.fillRectangle(0,0,canvas->getWidth(),canvas->getHeight());
  dirty=0;
  return 1;
  }



// Update the clear button:- each gui element (widget) in FOX
// receives a message during idle processing asking it to be updated.
// For example, buttons can be sensitized or desensitized when the
// state of the application changes.
// In this case, we desensitize the sender (the clear button) when
// the canvas has already been cleared, and sensitize it when it has
// been painted (as indicated by the dirty flag).
long ScribbleWindow::onUpdClear(FXObject* sender,FXSelector,void*){

  if(dirty)
    sender->handle(this,FXSEL(SEL_COMMAND,FXWindow::ID_ENABLE),NULL);
  else
    sender->handle(this,FXSEL(SEL_COMMAND,FXWindow::ID_DISABLE),NULL);

  return 1;
  }

// Handle the color message
long ScribbleWindow::onCmdColor(FXObject*,FXSelector,void* data){
  drawColor = static_cast<FXColor>(reinterpret_cast<uintptr_t>(data));
  return 1;
  }

long ScribbleWindow::onUpdColor(FXObject* sender,FXSelector /*sel*/,void* /*data*/){
  sender->handle(this,FXSEL(SEL_COMMAND,FXWindow::ID_SETINTVALUE),(void*)&drawColor);
  return 1;
  }






int main( int argc, char **argv ) {
  try {

// Make application
  FXApp application("Scribble","FoxTest");

  // Start app
  application.init(argc,argv);

  // Scribble window
  new ScribbleWindow(&application);

  // Create the application's windows
  application.create();

  // Run the application
  return application.run();
  } catch ( aol::Exception &el ) {
    el.dump();
#ifdef __MINGW32_VERSION
    system("PAUSE");
#endif
  }
  return 0;
}
