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

#include <aol.h>

#include "vtkFXGui.h"

int main ( int argc , char **argv ) {

  try {

    FXApp app ( "NAME", "vtkFOX next generation" );
    app.init ( argc, argv );

    vtkFXGui *window = new vtkFXGui ( &app );

    app.create();

#if defined(SIDE_BY_SIDE_STEREO_RENDERING) || defined(WIN32)
    window->show ( PLACEMENT_SCREEN );
#else
    window->show ( PLACEMENT_MAXIMIZED );
#endif

    return app.run();

  } catch ( aol::Exception &ex ) {
    ex.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
    return ( EXIT_FAILURE );
  }

}
