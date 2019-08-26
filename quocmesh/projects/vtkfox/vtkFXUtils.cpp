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

#include "vtkFXUtils.h"

FXString getDirectoryFromString ( const FXString &String ) {
  FXString directory;
  if ( String.contains ( "/" ) > 0 )
    directory = String.left ( String.find_last_of ( "/" ) );
  if ( String.contains ( "\\" ) > 0 )
    directory = String.left ( String.find_last_of ( "\\" ) );
  return directory;
}

FXString getShortenedString ( const FXString &String ) {
  FXString shortVersion;
  if ( String.length() > 40 ) {
    shortVersion = "...";
    shortVersion.append ( String.right ( 40 ) );
  } else {
    shortVersion = String;
  }

  return ( shortVersion );
}


void saveCamera ( vtkRenderer *Renderer ) {
  FXdouble dA[3], dV;
  FXint iV;
  vtkCamera *Camera = Renderer->GetActiveCamera();
  ofstream outfile ( "camera.dat", ofstream::binary );

  Camera->GetPosition ( dA );
  outfile.write ( ( char* ) dA, 3*sizeof ( FXdouble ) );
  cerr << "Position = " << dA[0] << ", " << dA[1] << ", " << dA[2] << endl;
  Camera->GetFocalPoint ( dA );
  outfile.write ( ( char* ) dA, 3*sizeof ( FXdouble ) );
  cerr << "FocalPoint = " << dA[0] << ", " << dA[1] << ", " << dA[2] << endl;
  Camera->GetViewUp ( dA );
  outfile.write ( ( char* ) dA, 3*sizeof ( FXdouble ) );
  cerr << "ViewUp = " << dA[0] << ", " << dA[1] << ", " << dA[2] << endl;
  dV = Camera->GetViewAngle();
  outfile.write ( ( char* ) ( &dV ), sizeof ( FXdouble ) );
  cerr << "ViewAngle = " << dV << endl;
  iV = Camera->GetParallelProjection();
  outfile.write ( ( char* ) ( &iV ), sizeof ( FXint ) );
  cerr << "ParallelProjection = " << iV << endl;
  dV = Camera->GetParallelScale();
  outfile.write ( ( char* ) ( &dV ), sizeof ( FXdouble ) );
  cerr << "ParallelScale = " << dV << endl;

  outfile.close();
}

void loadCamera ( vtkRenderer *Renderer ) {
  FXdouble dA[3], dV;
  FXint iV;
  ifstream infile ( "camera.dat", ofstream::binary );
  vtkCamera *Camera = Renderer->GetActiveCamera();

  infile.read ( ( char* ) dA, 3*sizeof ( FXdouble ) );
  Camera->SetPosition ( dA );
  infile.read ( ( char* ) dA, 3*sizeof ( FXdouble ) );
  Camera->SetFocalPoint ( dA );
  infile.read ( ( char* ) dA, 3*sizeof ( FXdouble ) );
  Camera->SetViewUp ( dA );
  infile.read ( ( char* ) ( &dV ), sizeof ( FXdouble ) );
  Camera->SetViewAngle ( dV );
  infile.read ( ( char* ) ( &iV ), sizeof ( FXint ) );
  Camera->SetParallelProjection ( iV );
  infile.read ( ( char* ) ( &dV ), sizeof ( FXdouble ) );
  Camera->SetParallelScale ( dV );
  Camera->Modified();

  Renderer->ResetCameraClippingRange();

  vtkLightCollection* Lights = Renderer->GetLights();
  FXint i, n = Lights->GetNumberOfItems();
  vtkLight *Light;
  FXdouble pos[3], fp[3];

  Camera->GetPosition ( pos );
  Camera->GetFocalPoint ( fp );
  Lights->InitTraversal();
  for ( i = 0; i < n; i++ ) {
    Light = Lights->GetNextItem();
    Light->SetPosition ( pos );
    Light->SetFocalPoint ( fp );
    Light->Modified();
  }
  infile.close();
}

void buildStripedPolyData ( vtkPolyData *MeshPolyData, vtkAppendPolyData *PolyDataCollection, const bool LeftView, const float StripeWidth ) {
  vtkFloatArray *norms = vtkFloatArray::New();
  norms->SetNumberOfComponents ( 3 );
  norms->SetNumberOfTuples ( 2 );
  norms->SetTuple3 ( 0, -1, -1, 0 );
  norms->SetTuple3 ( 1, 1, 1, 0 );
  const float width = StripeWidth;
  const int offset = static_cast<int> ( LeftView );
  for ( int j = 0; (2*j+offset)*width < 1; ++j ) {
    vtkPoints *points = vtkPoints::New();
    points->SetNumberOfPoints ( 2 );
    points->SetPoint ( 0, (2*j+offset)*width, (2*j+offset)*width, 0 );
    points->SetPoint ( 1, (2*j+1+offset)*width, (2*j+1+offset)*width, 0 );

    vtkPlanes *planes = vtkPlanes::New();
    planes->SetPoints ( points );
    planes->SetNormals ( norms );

    vtkClipPolyData *clipper = vtkClipPolyData::New();
    clipper->SetInput ( MeshPolyData );
    clipper->SetClipFunction ( planes );
    clipper->InsideOutOn ();
    clipper->GenerateClippedOutputOn ( );

    PolyDataCollection->AddInput ( clipper->GetOutput() );
    planes->Delete();
    clipper->Delete();
  }
  norms->Delete();
}

void savePolyDataToVTKFile ( vtkPolyData *MeshPolyData, const char *FileName ) {
  vtkPolyDataWriter *polyDataWriter = vtkPolyDataWriter::New();
  polyDataWriter->SetInput ( MeshPolyData );
  polyDataWriter->SetFileName ( FileName );
  polyDataWriter->SetFileTypeToBinary();
  polyDataWriter->Write();
  polyDataWriter->Delete();
}
