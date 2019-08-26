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

#ifndef __NARROW_H
#define __NARROW_H

#include <quoc.h>
#include <eikonalNA.h>
#include <scalarArray.h>
#include <levelSetDrawer.h>
#include <signedDistanceOp.h>
#include <FEOpInterface.h>
#include <baseFunctionSet.h>
#include <solver.h>
#include <mcm.h>
#include <levelSet.h>
#include <qmException.h>

#include <narrowBandGrid.h>
#include <narrowBandConfigurators.h>
#include <subGridSparseMatrix.h>

namespace nb {

template <typename RealType, typename NarrowBandGridType, qc::Dimension Dim>
class InitNarrowBand { };

/**
 * compute the signed distance function to the level lines of the function given by the argument
 * with specified isovalue and write the result in the destination vector.
 * 2d-only.
 * \author Droske
 * \todo Implement a better way to get a ConfiguratorType template for SignedDistanceBase. Just explicitly specifying QuocConfiguratorTraitMultiLin is not very clean.
 */
template <typename RealType, typename NarrowBandGridType>
class InitNarrowBand<RealType, NarrowBandGridType, qc::QC_2D> : public qc::SignedDistanceBase<qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> >, eik::EikonalNANarrowBand<RealType, NarrowBandGridType> > {
public:
  InitNarrowBand ( const qc::GridDefinition &Grid, NarrowBandGridType &NarrowBandGrid ) :
      qc::SignedDistanceBase<qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> >, eik::EikonalNANarrowBand<RealType, NarrowBandGridType> > ( Grid ) {
    this->_eik.setNarrowBandGrid ( NarrowBandGrid );
  }

  virtual ~InitNarrowBand( ) {}

  void operator() ( const qc::Array<RealType> &LevelSetFunction, RealType Epsilon ) const {
    this->_eik.reset();
    this->initializeBoundary ( LevelSetFunction );
    this->_eik.marchNeighborHood ( Epsilon );
  }
};

/**
 * compute the signed distance function to the level lines of the function given by the argument
 * with specified isovalue and write the result in the destination vector.
 * 3d-only.
 * \author Droske
 * \todo Implement a better way to get a ConfiguratorType template for SignedDistanceBase. Just explicitly specifying QuocConfiguratorTraitMultiLin is not very clean.
 */
template <typename RealType, typename NarrowBandGridType>
    class InitNarrowBand<RealType, NarrowBandGridType, qc::QC_3D> : public qc::SignedDistanceBase3D<qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType,qc::QC_3D, 3> >, eik::EikonalNANarrowBand3D<RealType, NarrowBandGridType> > {
public:
  InitNarrowBand ( const qc::GridDefinition &Grid, NarrowBandGridType &NarrowBandGrid ) :
    qc::SignedDistanceBase3D<qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType,qc::QC_3D, 3> >, eik::EikonalNANarrowBand3D<RealType, NarrowBandGridType> > ( Grid ) {
    this->_eik.setNarrowBandGrid ( NarrowBandGrid );
  }

  virtual ~InitNarrowBand( ) {}

  void operator() ( const qc::ScalarArray<RealType, qc::QC_3D> &LevelSetFunction, RealType Epsilon ) const {
    this->_eik.reset();
    this->initializeBoundary ( LevelSetFunction );
    this->_eik.marchNeighborHood ( Epsilon );
    for ( int i = 0; i < this->_eik.getTimeField().size(); i++ ) {
      if ( LevelSetFunction[i] < 0. ) this->_eik.getTimeField() [i] = - ( this->_eik.getTimeField() [i] - 1.0 );
      else this->_eik.getTimeField() [i] = this->_eik.getTimeField() [i] - 1.0;
    }
  }
};

/**
 * Initializes a narrow band, by inserting all elements where the value of
 * argument LevelsetFunction is lower or equal to argument Value at all
 * nodes of the element.
 *
 * \author Berkels
 */
template <typename NarrowConfiguratorType, typename FullConfiguratorType>
void initNarrowBandFromSubLevelset ( const aol::Vector<typename FullConfiguratorType::RealType> &LevelsetFunction,
                                     typename NarrowConfiguratorType::InitType &NarrowGrid,
                                     const typename FullConfiguratorType::RealType Value = 0 ) {
  // First create a mask of the sublevelset, so that we can use BitArray::elementTrue.
  qc::BitArray<FullConfiguratorType::Dim> mask ( qc::GridSize<FullConfiguratorType::Dim>::createFrom ( NarrowGrid ) );
  for ( int i = 0; i < mask.size(); ++i ) {
    if ( LevelsetFunction[i] <= Value )
      mask.set( i, true );
  }

  FullConfiguratorType config ( NarrowGrid.getFullGrid() );
  typedef typename FullConfiguratorType::ElementIteratorType IteratorType;
  const typename IteratorType::EndType end_it =  config.end();
  for ( IteratorType it = config.begin(); it != end_it; ++it ) {
    if ( mask.elementTrue ( *it ) )
      NarrowGrid.insert ( *it );
  }
}

template <qc::Dimension Dim>
class NarrowBandBoundaryMeshExport { };

template <>
class NarrowBandBoundaryMeshExport<qc::QC_3D> {
protected:
  NarrowBandMaskedGrid<double, qc::QC_3D> _narrow;
public:
  NarrowBandBoundaryMeshExport ( NarrowBandMaskedGrid<double, qc::QC_3D> &narrow )
      : _narrow ( narrow ) {}


  void write( ) {
    NarrowBandMaskedGrid<double, qc::QC_3D>::ITER_MODE itermode = _narrow.getIteratorMode( );

    _narrow.setIteratorMode ( NarrowBandMaskedGrid<double, qc::QC_3D>::EXT );

    write ( "outer_elements.pov", "neumann_faces.pov", "neumann_nodes.pov",
            _narrow.begin(), _narrow.end(), _narrow.ebegin(), _narrow.eend(), _narrow.bbegin(), _narrow.bend() );

    _narrow.setIteratorMode ( NarrowBandMaskedGrid<double, qc::QC_3D>::INT );

    write ( "inner_elements.pov", "dirichlet_faces.pov", "dirichlet_nodes.pov",
            _narrow.begin(), _narrow.end(), _narrow.eDbegin(), _narrow.eDend(), _narrow.bDbegin(), _narrow.bDend() );

    _narrow.setIteratorMode ( itermode );
  }

  void write ( const char *element_filename,
               const char *edge_filename,
               const char *node_filename,
               NarrowBandMaskedGrid<double, qc::QC_3D>::BeginIterType elementbegin,
               NarrowBandMaskedGrid<double, qc::QC_3D>::EndIterType elementend,
               NarrowBandMaskedGrid<double, qc::QC_3D>::eiterator edgebegin,
               NarrowBandMaskedGrid<double, qc::QC_3D>::eiterator edgeend,
               NarrowBandMaskedGrid<double, qc::QC_3D>::biterator nodesbegin,
               NarrowBandMaskedGrid<double, qc::QC_3D>::biterator nodesend ) {
    ofstream out_elements ( element_filename );
    ofstream out_edges ( edge_filename );
    ofstream out_nodes ( node_filename );


    vector<aol::Vec3<double> > vertices;
    vector<aol::Vec3<double> > normals;
    vector<aol::Vec3<int> > faceindices;

    const double h = _narrow.getFullGrid().H();
    qc::ScalarArray<int, qc::QC_3D> indices ( qc::GridSize<qc::QC_3D>::createFrom ( _narrow ) );

    indices.setAll ( -1 );
    int vindex = 0;

    int offsets[8][3] = { { 0, 0, 0 }
                          , { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 },
                          { 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 }
                        };


    out_elements << "union {\n";
    for ( NarrowBandMaskedGrid<double, qc::QC_3D>::iterator it = elementbegin; it != elementend; ++it ) {
      aol::Vec3<double> coord ( it->x(), it->y(), it->z() );
      coord[0] += 0.1 * h;
      coord[1] += 0.1 * h;
      coord[2] += 0.1 * h;
      coord *= h;
      out_elements << "box { <" << coord[0] << "," << coord[1] << "," << coord[2] << ">, ";
      coord[0] += 0.8 * h;
      coord[1] += 0.8 * h;
      coord[2] += 0.8 * h;
      out_elements << "<" << coord[0] << "," << coord[1] << "," << coord[2] << "> }\n";
    }
    out_elements << "}\n";

    for ( NarrowBandMaskedGrid<double, qc::QC_3D>::eiterator eit = edgebegin; eit != edgeend; ++eit ) {
      qc::CoordType nodes[4];
      for ( int i = 0; i < 4; i++ ) {
        nodes[i].set ( eit->_el.x() + offsets[ eit->_locNode[ i ] ][ 0 ],
                       eit->_el.y() + offsets[ eit->_locNode[ i ] ][ 1 ],
                       eit->_el.z() + offsets[ eit->_locNode[ i ] ][ 2 ] );

      }

      aol::Vec3<int> fi;

      for ( int i = 0; i < 3; i++ ) {
        if ( indices.get ( nodes[i] ) >= 0 ) {
          fi[i] = indices.get ( nodes[i] );
        } else {
          indices.set ( nodes[i], vindex );
          aol::Vec3<double> coord ( nodes[i][0], nodes[i][1], nodes[i][2] );
          coord *= h;
          vertices.push_back ( coord );
          fi[i] = vindex++;
        }
      }

      faceindices.push_back ( fi );

      if ( indices.get ( nodes[3] ) >= 0 ) {
        fi[0] = indices.get ( nodes[3] );
      } else {
        indices.set ( nodes[3], vindex );
        aol::Vec3<double> coord ( nodes[3][0], nodes[3][1], nodes[3][2] );
        coord *= h;
        vertices.push_back ( coord );
        fi[0] = vindex++;
      }

      faceindices.push_back ( fi );
    }
    out_edges << "mesh2 {\n  vertex_vectors {\n    " << vertices.size() << "," << endl;
    for ( int i = 0; i < static_cast< int > ( vertices.size() ) ; i++ ) {
      out_edges << "    <" << vertices[i][0] << ", " << vertices[i][1] << ", " << vertices[i][2] << ">";
      if ( i < static_cast< int > ( vertices.size() ) - 1 ) {
        out_edges << "," << endl;
      } else {
        out_edges << endl << "  }\n";
      }
    }
    out_edges << "  face_indices {\n    " << faceindices.size() << "," << endl;
    for ( int i = 0; i < static_cast<int> ( faceindices.size() ); i++ ) {
      out_edges << "    <" << faceindices[i][0] << ", " << faceindices[i][1] << ", " << faceindices[i][2] << ">";
      if ( i < static_cast< int > ( faceindices.size() ) - 1 ) {
        out_edges << "," << endl;
      } else {
        out_edges << endl << "  }\n";
      }
    }
    out_edges << "}";

    out_nodes << "union {\n";
    for ( NarrowBandMaskedGrid<double, qc::QC_3D>::biterator bit = nodesbegin; bit != nodesend; ++bit ) {
      aol::Vec3<double> coord ( bit->x(), bit->y(), bit->z() );
      coord *= h;
      out_nodes << "sphere { <" << coord[0] << "," << coord[1] << "," << coord[2] << ">, " << h * 0.15 << " } \n";
    }
    out_nodes << "}\n";
  }
};

} // end of namespace nb.

#endif
