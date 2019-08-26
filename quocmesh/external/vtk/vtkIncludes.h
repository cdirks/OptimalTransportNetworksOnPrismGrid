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

#ifndef __VTKINCLUDES_H
#define __VTKINCLUDES_H

#ifdef __GNUC__
#pragma GCC system_header
#endif

// in alphabetical order:
#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkBMPReader.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkClipPolyData.h>
#include <vtkColorTransferFunction.h>
#include <vtkCommand.h>
#include <vtkConeSource.h>
#include <vtkContourFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkDecimatePro.h>
#include <vtkFloatArray.h>
#include <vtkGL2PSExporter.h>
#include <vtkGraphicsFactory.h>
#include <vtkImageAppendComponents.h>
#include <vtkImageData.h>
#include <vtkImagingFactory.h>
#include <vtkInteractorStyleTrackball.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleUnicam.h>
#include <vtkLight.h>
#include <vtkLightCollection.h>
#include <vtkMarchingCubes.h>
#include <vtkMarchingSquares.h>
#include <vtkMatrix4x4.h>
#include <vtkOpenGLVolumeTextureMapper2D.h>
#include <vtkOpenGLVolumeTextureMapper3D.h>
#include <vtkPNGReader.h>
#include <vtkPNGWriter.h>
#include <vtkPiecewiseFunction.h>
#include <vtkPlane.h>
#include <vtkPlanes.h>
#include <vtkPlaneSource.h>
#include <vtkPLYWriter.h>
#include <vtkPNMReader.h>
#include <vtkPointData.h>
#include <vtkPointPicker.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkRendererCollection.h>
#include <vtkRenderLargeImage.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkStripper.h>
#include <vtkSuperquadricSource.h>
#include <vtkTexture.h>
#include <vtkTexturedActor2D.h>
#include <vtkTextureMapToPlane.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkVolumeProperty.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeRayCastFunction.h>
#include <vtkVolumeRayCastMIPFunction.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkWarpScalar.h>
#include <vtkWindowToImageFilter.h>

#endif // __VTKINCLUDES_H
