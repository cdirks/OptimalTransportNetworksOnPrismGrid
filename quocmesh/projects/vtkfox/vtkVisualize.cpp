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

#include "vtkVisualize.h"

void updateCamera ( const double* rotAngles, double Zoom, vtkRenderer *Renderer ) {

  //   vtkCamera *Camera = Renderer->GetActiveCamera();
  //   cerr << "update camera" << endl;
  //   FXdouble pos[3], foc[3], vup[3], cart[3];
  //   Camera->GetPosition (  pos );
  //   //   cerr << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  //   Camera->GetFocalPoint( foc );
  //   //   cerr << foc[0] << " " << foc[1] << " " << foc[2] << endl;
  //   Camera->GetViewUp (    vup );
  //   //   cerr << vup[0] << " " << vup[1] << " " << vup[2] << endl;

  Renderer->GetActiveCamera()->SetPosition   ( 0.5, 0.5, 3.5 ); // approximately the initial position
  Renderer->GetActiveCamera()->SetFocalPoint ( 0.5, 0.5, 0.5 );
  Renderer->GetActiveCamera()->SetViewUp     ( 0.0, 1.0, 0.0 );
  Renderer->ResetCamera();
  Renderer->GetActiveCamera()->Roll (        rotAngles[2] );
  Renderer->GetActiveCamera()->Azimuth (     rotAngles[1] );
  Renderer->GetActiveCamera()->Elevation (   rotAngles[0] );
  Renderer->GetActiveCamera()->Zoom ( Zoom );
}
