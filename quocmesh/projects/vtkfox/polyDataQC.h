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

#ifndef __POLYDATAQC_H
#define __POLYDATAQC_H

#include <foxIncludes.h>
#include <vtkIncludes.h>

#include <scalarArray.h>
#include <triangMesh.h>

// forward declaration of the PolyDataGroup to avoid inclusion
// of "polyDataGroup.h"
class PolyDataGroup;

/**
 *  Base class for all displayable objects.
 *
 *  This class holds all properties that are stored inside
 *  vtkFOX to display a surface / a volume / etc. pp.
 *
 *  This class contains all non-GUI information associated
 *  to a displayable data set (surface, volume etc.).
 *  The pipeline goes as follows:
 *
 *  PolyData -> DataNormals -> Stripper -> PolyDataMapper
 *  -> Actor ( + Property) -> Renderer -> RenderWindow
 *
 *  For an explaination what these objects in the graphics
 *  pipeline are good for, look at the documentation of
 *  the data members below.
 *
 *  When one has properties that are
 *  only interesting to one special class of
 *  surfaces (e. g. the level value to be displayed), one
 *  can derive a subclass and endow it with additional
 *  handling routines. To display specialized properties
 *  or user abilities, place FOX widgets into
 *  _owningGroup.getIndivualPropertiesFrame(). Display
 *  and hide your widgets there everywhere activate() or
 *  inactivate() is called.
 *
 *  \author von Deylen
 */
class PolyDataQC {
public:
  //! "standard" constructor, getting only description and mesh data
  PolyDataQC ( PolyDataGroup *      owningGroup,
               FXString             description,
               vtkPolyData *        polyData );

  // Destructor
  virtual ~PolyDataQC ();

  //! makes this object visible, but does not activate.
  virtual void show();

  //! makes this object invisible, but does not inactivate.
  virtual void hide();

  //! performs special operations when this object is activated.
  //! Activating an object does not not imply that it is visible.
  virtual void activate();

  //! performs operations that are needed when this object no longer activated.
  virtual void inactivate();

  //! handle messages that have been sent to PolyDataGroup with Selector ID_POLYDATA_MESSAGE
  virtual long handle ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ ) {  return 1;  }

  // member getters
  bool                 isVisible () const             {  return ( _visible != 0 );  }
  FXDataTarget *       getVisibleDataTarget ()        {  return _visibleDataTarget;  }
  FXFrame *            getButton ()                   {  return _button;  }
  FXString             getDescription () const        {  return _description;  }
  vtkPolyData *        getPolyData ()                 {  return _polyData;  }
  vtkActor *           getActor ()                    {  return _actor;  }
  vtkPolyDataMapper *  getMapper ()                   {  return _mapper;  }
  vtkProperty *        getProperties ()               {  return _properties;  }
  vtkPolyDataNormals * getNormals ()                  {  return _normals;  }
  vtkStripper *        getStripper ()                 {  return _stripper;  }

protected:

  //! protected constructor for setting individual member values
  PolyDataQC ( PolyDataGroup *      owningPolyDataGroup,
               bool                 visible,
               FXDataTarget *       visibleDataTarget,
               FXFrame *            button,
               FXString             description,
               vtkPolyData *        polyData,
               vtkActor *           actor,
               vtkPolyDataMapper *  mapper,
               vtkProperty *        properties,
               vtkPolyDataNormals * normals,
               vtkStripper *        stripper );

  //! protected constructor for setting non-FOX member variables individually.
  //! VTK pipeline objects are created, but will not be
  //! connected to each other.
  //! To make up for this, call connectVTKPipeline(vtkPolyData *)
  //! anywhen afterwards.
  PolyDataQC ( PolyDataGroup *      owningGroup,
               bool                 visible,
               FXString             description,
               vtkPolyData *        polyData    = vtkPolyData::New(),
               vtkActor *           actor       = vtkActor::New(),
               vtkPolyDataMapper *  mapper      = vtkPolyDataMapper::New(),
               vtkProperty *        properties  = vtkProperty::New(),
               vtkPolyDataNormals * normals     = vtkPolyDataNormals::New(),
               vtkStripper *        stripper    = vtkStripper::New() );

  void connectVTKPipeline();
  void updateVTKCanvas();

  // ***** DATA MEMBERS *****

  //! The PolyDataGroup in which one's poly data list
  //! this object is an entry
  PolyDataGroup *      _owningGroup;

  //! "is visible" variable
  int                  _visible;
  // and corresponding data target
  FXDataTarget *       _visibleDataTarget;

  //! button corresponding to the object
  //! in the control panel "Objects".
  FXFrame *            _button;

  //! textual description (filename in most cases)
  //! displayed as button text and drop down list entry.
  FXString             _description;

  //! the given data, mostly output
  //! of an algorithm like "triangMeshToVTK" or
  //! "tpcfeToVTK".
  vtkPolyData *        _polyData;

  //! bundle for data and display properties,
  //! namely vtkPolyDataMapper and vtkProperty.
  //! (zoom, rotation etc).
  vtkActor *           _actor;

  //! mapping from maps polygonal data to graphics primitives
  vtkPolyDataMapper *  _mapper;

  //! set of display properties for the surface / volume
  vtkProperty *        _properties;

  //! filter that endows the polygonal surface with
  //! surface normals
  vtkPolyDataNormals * _normals;

  //! filter that collects individual triangles into
  //! strips, i. e. small lists of triangles,
  //! each consecutive pair of them sharing one
  //! edge.
  vtkStripper *        _stripper;

  //! pointer to vertical frame inside the "Controls" panel.
  //! Is not used in this base class, but may be filled with
  //! properties regulator in derived classes. If you want
  //! to use it, call "_owningGroup.getIndivualPropertiesPanel()"
  //! which creates a frame in the "Controls" panel and returns a
  //! pointer to it.
  FXTabItem *          _propertiesTab;
  FXVerticalFrame *    _propertiesFrame;

private:
  static vtkPolyDataNormals * createDataNormalsFilter ( vtkPolyData * polyData );
  void connectActorPropertiesStripper ();
};

// ==========================================================================

class IsoSurfaceData : public PolyDataQC {
public:
    IsoSurfaceData ( PolyDataGroup *      owningGroup,
                     FXString             description,
                     vtkImageData *       levelData,
                     FXdouble             levelValueToDisplay );
    ~IsoSurfaceData();

    long handle ( FXObject* /*sender*/, FXSelector /*sel*/, void* /*data*/ );

protected:
  double              _levelValueToDisplay;

  vtkImageData *      _levelData;
  vtkContourFilter *  _contourFilter;

  FXDataTarget *      _levelValueTarget;
};

// ==========================================================================

class IsoSurfaceDataWithTextureMap : public PolyDataQC {
public:
  IsoSurfaceDataWithTextureMap ( PolyDataGroup * owningGroup,
                                 FXString        description,
                                 std::string     levelsetFilename,
                                 std::string     tCoordXFilename,
                                 std::string     tCoordYFilename,
                                 std::string     tBitmapFilename,
                                 FXdouble        levelValueToDisplay );

  ~IsoSurfaceDataWithTextureMap();
  vtkPolyData * createSurfaceWithTCoordsFromZeroLevelSet ( const qc::ScalarArray<FXfloat, qc::QC_3D> & levelData,
                                                           const qc::ScalarArray<FXfloat, qc::QC_3D> & tCoordX,
                                                           const qc::ScalarArray<FXfloat, qc::QC_3D> & tCoordY );

  long handle ( FXObject* sender, FXSelector /*sel*/, void* /*data*/ );

  void updateLevelValueToDisplay();

protected:
  qc::ScalarArray<FXfloat, qc::QC_3D> _levelData;
  qc::ScalarArray<FXfloat, qc::QC_3D> _tCoordX;
  qc::ScalarArray<FXfloat, qc::QC_3D> _tCoordY;
  FXdouble                            _levelValueToDisplay;

  aol::TriangMesh<FXfloat>             _levelSurface;

  vtkImageReader *                     _tBitmapReader;
  vtkTexture *                         _texture;

  FXDataTarget *                       _levelValueTarget;
};

#endif
