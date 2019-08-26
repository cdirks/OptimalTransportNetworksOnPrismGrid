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

#include "tpcfeToVTK.h"

#include <tpCFEUtils.h>

namespace tpcfe {

template < class GridType >
class TriangulationToVTKConverter {
public:
  void setVTKPointsTriangles ( const typename CFEInterfaceTriangulationGenerator<GridType>::VertexVectorType &interfaceVertexVector, const typename CFEInterfaceTriangulationGenerator<GridType>::VertexVectorType &boundaryVertexVector,
                               const typename CFEInterfaceTriangulationGenerator<GridType>::TriangVectorType &interfaceTriangVector, const typename CFEInterfaceTriangulationGenerator<GridType>::TriangVectorType &boundaryTriangVector,
                               vtkPolyData *polyData ) {

    vtkPoints* points = vtkPoints::New();
    points->SetNumberOfPoints ( interfaceVertexVector.size() + boundaryVertexVector.size() );
    int counter = 0;

    for ( typename CFEInterfaceTriangulationGenerator<GridType>::VertexVectorType::const_iterator it = interfaceVertexVector.begin(); it != interfaceVertexVector.end(); ++it ) {
      points->SetPoint ( counter, ( *it ) [0], ( *it ) [1], ( *it ) [2] );
      ++counter;
    }

    for ( typename CFEInterfaceTriangulationGenerator<GridType>::VertexVectorType::const_iterator it = boundaryVertexVector.begin(); it != boundaryVertexVector.end(); ++it ) {
      points->SetPoint ( counter, ( *it ) [0], ( *it ) [1], ( *it ) [2] );
      ++counter;
    }

    polyData->SetPoints ( points );

    // need to delete points??
    points->Delete();


    polyData->Allocate ( interfaceTriangVector.size() + boundaryTriangVector.size() );

    for ( typename CFEInterfaceTriangulationGenerator<GridType>::TriangVectorType::const_iterator it = interfaceTriangVector.begin(); it != interfaceTriangVector.end(); ++it ) {
      vtkIdType triang[3] = { ( *it ) [0], ( *it ) [1], ( *it ) [2] };
      polyData->InsertNextCell ( VTK_TRIANGLE, 3, triang );
    }

    const int offset = interfaceVertexVector.size();

    for ( typename CFEInterfaceTriangulationGenerator<GridType>::TriangVectorType::const_iterator it = boundaryTriangVector.begin(); it != boundaryTriangVector.end(); ++it ) {
      vtkIdType triang[3] = { ( *it ) [0] + offset, ( *it ) [1] + offset, ( *it ) [2] + offset };
      polyData->InsertNextCell ( VTK_TRIANGLE, 3, triang );
    }

  } // end setVTKTriangles


  void setVTKVertexData ( const typename CFEInterfaceTriangulationGenerator<GridType>::VertexDataVectorType &interfaceVertexDataVector,
                          const typename CFEInterfaceTriangulationGenerator<GridType>::VertexDataVectorType &boundaryVertexDataVector,
                          vtkPolyData *polyData, PolyDataGroup *polyDataGroup ) {
    vtkFloatArray* f = vtkFloatArray::New();
    f->SetNumberOfValues ( interfaceVertexDataVector.size() + boundaryVertexDataVector.size() );
    int counter = 0;
    float min = 1e20, max = -1e20;

    for ( typename CFEInterfaceTriangulationGenerator<GridType>::VertexDataVectorType::const_iterator it = interfaceVertexDataVector.begin(); it != interfaceVertexDataVector.end(); ++it ) {
      float val = *it;
      if ( val > max )
        max = val;
      if ( val < min )
        min = val;
      f->SetValue ( counter, val );
      ++counter;
    }

    for ( typename CFEInterfaceTriangulationGenerator<GridType>::VertexDataVectorType::const_iterator it = boundaryVertexDataVector.begin(); it != boundaryVertexDataVector.end(); ++it ) {
      float val = *it;
      if ( val > max )
        max = val;
      if ( val < min )
        min = val;
      f->SetValue ( counter, val );
      ++counter;
    }

    polyData->GetPointData()->SetScalars ( f );

    // need to delete this??
    f->Delete();

    cerr << endl << endl << min << " " << max << endl << endl;

    polyDataGroup->setColorMapRange ( min, max );
  } // end setVTKVertexData

}
; // end class TriangulationToVTKConverter


template< typename GridType >
vtkPolyData* TpcfeToVTK<GridType>::makePolyData ( ) {

  char levelset_filename[1024];
  parser.getString ( "geometry_filename", levelset_filename );
  cerr << levelset_filename << endl;
  qc::ScalarArray<RealType, qc::QC_3D> levelset ( levelset_filename );

  GridType grid ( levelset.getSize() );
  grid.setDomainFrom ( levelset );
  grid.detectAndInitVirtualNodes();

  vtkPolyData*    polyData  = vtkPolyData::New();

  enum SHOW { DISPLAY_SCALAR_UNDEFORMED,
              DISPLAY_SCALAR_DEFORMED,
              DISPLAY_VONMISES_DEFORMED,
              DISPLAY_GEOMETRY };

  SHOW show = DISPLAY_GEOMETRY;
  if ( parser.hasVariable ( "scalar_data" ) && !( parser.hasVariable ( "deformations_fnmask" ) ) )
    show = DISPLAY_SCALAR_UNDEFORMED;
  else if ( parser.hasVariable ( "scalar_data" ) && parser.hasVariable ( "deformations_fnmask" ) )
    show = DISPLAY_SCALAR_DEFORMED;
  else if ( !( parser.hasVariable ( "scalar_data" ) ) && parser.hasVariable ( "deformations_fnmask" ) )
    show = DISPLAY_VONMISES_DEFORMED;
  else
    show = DISPLAY_GEOMETRY;

  switch ( show ) {
  case DISPLAY_SCALAR_UNDEFORMED:
    { cerr << "Loading and displaying scalar data ... " << endl;

      qc::ScalarArray<RealType, qc::QC_3D> scalarData ( grid );
      char scalarDataFN[1024];
      parser.getString ( "scalar_data", scalarDataFN );
      scalarData.load ( scalarDataFN );
      tpcfe::CFEInterfaceTriangulationWithScalarDataGenerator<GridType> itg ( grid, scalarData );
      itg.determineInterfaceAndBoundaryTriangulation();

      tpcfe::TriangulationToVTKConverter<GridType> ttv;
      ttv.setVTKPointsTriangles ( itg.getInterfaceVertexVectorRef(), itg.getBoundaryVertexVectorRef(), itg.getInterfaceTriangVectorRef(), itg.getBoundaryTriangVectorRef(), polyData );
      ttv.setVTKVertexData ( itg.getInterfaceVertexDataVectorRef(), itg.getBoundaryVertexDataVectorRef(), polyData , _polyDataGroup );
    }
    break;

  case DISPLAY_SCALAR_DEFORMED:
    { cerr << "Loading and displaying scalar data on deformed geometry ... " << endl;

      qc::ScalarArray<RealType, qc::QC_3D> scalarData ( grid );
      char scalarDataFN[1024];
      parser.getString ( "scalar_data", scalarDataFN );
      scalarData.load ( scalarDataFN );
      tpcfe::CFEInterfaceTriangulationWithScalarDataGenerator<GridType> itg ( grid, scalarData );
      itg.determineInterfaceAndBoundaryTriangulation();

      qc::MultiArray< RealType, 3 > deformation ( grid );
      char deformations_filenamemask[1024];
      parser.getString ( "deformations_fnmask", deformations_filenamemask );
      deformation.load ( deformations_filenamemask );

      RealType scale_def = static_cast<RealType> ( parser.getDouble ( "scaleto" ) / deformation.getMaxAbsValue() );
      itg.deformVertices ( deformation, scale_def );

      tpcfe::TriangulationToVTKConverter<GridType> ttv;
      ttv.setVTKPointsTriangles ( itg.getInterfaceVertexVectorRef(), itg.getBoundaryVertexVectorRef(), itg.getInterfaceTriangVectorRef(), itg.getBoundaryTriangVectorRef(), polyData );
      ttv.setVTKVertexData ( itg.getInterfaceVertexDataVectorRef(), itg.getBoundaryVertexDataVectorRef(), polyData , _polyDataGroup );
    }
    break;

  case DISPLAY_VONMISES_DEFORMED:
    { cerr << "Loading and displaying deformations" << endl;

      const RealType
        E      = static_cast<RealType> ( parser.getDouble ( "Emodulus" ) ),
        nu     = static_cast<RealType> ( parser.getDouble ( "nu" ) ),
        lambda = nu * E / ( ( 1.0 + nu ) * ( 1 - 2.0 * nu ) ),
        mu     = E / ( 2.0 * ( 1.0 + nu ) ),
        aleph  = static_cast<RealType> ( parser.getDouble ( "aleph" ) );

      qc::MultiArray< RealType, 3 > deformation ( grid );
      char deformations_filenamemask[1024];
      parser.getString ( "deformations_fnmask", deformations_filenamemask );
      deformation.load ( deformations_filenamemask );

      RealType scale_def = static_cast<RealType> ( parser.getDouble ( "scaleto" ) / deformation.getMaxAbsValue() );

      tpcfe::CFEInterfaceTriangulationWithVonMisesStressGenerator<GridType> itg ( grid, deformation, lambda, mu, aleph );
      itg.determineInterfaceAndBoundaryTriangulation();
      itg.deformVertices ( scale_def );

      tpcfe::TriangulationToVTKConverter<GridType> ttv;
      ttv.setVTKPointsTriangles ( itg.getInterfaceVertexVectorRef(), itg.getBoundaryVertexVectorRef(), itg.getInterfaceTriangVectorRef(), itg.getBoundaryTriangVectorRef(), polyData );

      ttv.setVTKVertexData ( itg.getInterfaceVertexDataVectorRef(), itg.getBoundaryVertexDataVectorRef(), polyData, _polyDataGroup );

    }
    break;

  case DISPLAY_GEOMETRY:
    { cerr << "Parameter file does not contain entry for data. Will just display geometry." << endl;

      LevelsetToVTK<RealType> ltov ( levelset );

      if ( parser.hasVariable ( "vtkfox_clip_min_x" ) ) { // assume other clipping parameters present as well
        polyData = ltov.makePolyData( aol::Vec3<int>( parser.getInt( "vtkfox_clip_min_x" ), parser.getInt( "vtkfox_clip_min_y" ), parser.getInt( "vtkfox_clip_min_z" ) ),
                                      aol::Vec3<int>( parser.getInt( "vtkfox_clip_max_x" ), parser.getInt( "vtkfox_clip_max_y" ), parser.getInt( "vtkfox_clip_max_z" ) ) );

      } else { // no clipping
        polyData = ltov.makePolyData();
      }

    }
    break;

  default:
    throw aol::UnimplementedCodeException( "Selected SHOW mode not implemented", __FILE__, __LINE__);
    break;
  }

  return ( polyData );
}



template< typename RealType >
vtkPolyData* LevelsetToVTK<RealType>::makePolyData ( ) {

  vtkPolyData*    polyData  = vtkPolyData::New();

  typedef tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > GridType;

  GridType grid ( _levelset.getSize() );
  grid.setDomainFrom ( _levelset );
  grid.detectVirtualNodes();

  tpcfe::CFEInterfaceTriangulationGenerator<GridType> itg ( grid );
  itg.determineInterfaceAndBoundaryTriangulation();

  tpcfe::TriangulationToVTKConverter<GridType> ttv;
  ttv.setVTKPointsTriangles ( itg.getInterfaceVertexVectorRef(), itg.getBoundaryVertexVectorRef(), itg.getInterfaceTriangVectorRef(), itg.getBoundaryTriangVectorRef(), polyData );

  return ( polyData );
}


template< typename RealType >
vtkPolyData* LevelsetToVTK<RealType>::makePolyData ( const aol::Vec3<int> &clip_lower, const aol::Vec3<int> &clip_upper ) {

  vtkPolyData*    polyData  = vtkPolyData::New();

  typedef tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > GridType;

  GridType grid ( _levelset.getSize() );
  grid.setDomainFrom ( _levelset );
  grid.detectVirtualNodes();

  tpcfe::CFEInterfaceTriangulationGenerator<GridType> itg ( grid );
  itg.determineInterfaceAndBoundaryTriangulation( clip_lower, clip_upper ); // note that this draws an "internal planar boundary" at the clipping position, but does NOT remove triangles from the interface triangulation!

  tpcfe::TriangulationToVTKConverter<GridType> ttv;
  ttv.setVTKPointsTriangles ( itg.getInterfaceVertexVectorRef(), itg.getBoundaryVertexVectorRef(), itg.getInterfaceTriangVectorRef(), itg.getBoundaryTriangVectorRef(), polyData );

  return ( polyData );
}

}

template class tpcfe::TpcfeToVTK< tpcfe::CFEGrid< double, tpcfe::CFE_CD > >;

template class tpcfe::LevelsetToVTK<float>;
// double not needed

