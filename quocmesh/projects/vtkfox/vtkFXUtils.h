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

#ifndef __VTKFXUTILS_H
#define __VTKFXUTILS_H

#include <foxIncludes.h>
#include <vtkIncludes.h>
#include <scalarArray.h>

FXString getDirectoryFromString ( const FXString &String );

//! if the string is longer than 40 chars, then remove the leading part and replace it by "..."
FXString getShortenedString ( const FXString &String );

void saveCamera ( vtkRenderer *Renderer );
void loadCamera ( vtkRenderer *Renderer );
void buildStripedPolyData ( vtkPolyData *MeshPolyData, vtkAppendPolyData *PolyDataCollection, const bool LeftView, const float StripeWidth );
void savePolyDataToVTKFile ( vtkPolyData *MeshPolyData, const char *FileName );

#endif
