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

/** \file
 *  \brief usage of Composite Finite Element iterators
 *
 *  Using an iterator over all elements and a nested iterator over all
 *  tetrahedra in an element, we print the list of all inner
 *  tetrahedra.
 *
 *  \author Schwen
 */

#include <tpCFEGrid.h>
#include <tpCFELevelsets.h>
#include <shapeLevelsetGenerator.h>

int main ( int, char** ) {
  try {

    typedef double RealType;
    typedef tpcfe::CFEGrid< RealType, tpcfe::CFE_CD> GridType;

    GridType grid ( 3 );

    qc::ScalarArray<RealType, qc::QC_3D> levelset( grid );
    qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );
    grid.setDomainFrom ( levelset );
    grid.detectAndInitVirtualNodes();

    std::ostringstream output;

    const qc::GridSize<qc::QC_3D> gridSize ( grid );
    for ( GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
      tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

      el.computeAssembleData ( grid ); // determine geometric data for the virtual subdivision of the element

      for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el, -1 /* inner tetrahedra */ ); tit.notAtEnd(); ++tit ) {
        const tpcfe::CFETetra<RealType> &tet = *tit;

        // print tet's vertices:
        for ( short v = 0; v < 4; ++v ) {
          aol::Vec3<RealType> vertex;
          tet.computeGlobalCoordinate( vertex, el, v );

          const int
            lIdx0 = ( *tit ) ( v, 0 ),
            lIdx1 = ( *tit ) ( v, 1 ),
            gIdx0 = el.globalIndex(lIdx0),
            gIdx1 = ( lIdx1 == tpcfe::NODE_NOT_VIRTUAL ? 0 : el.globalIndex(lIdx1) );
          const GridType::VNIndexType gIdx  = tpcfe::CFEVirtualNode< RealType, tpcfe::CFE_CD, RealType >::mapIndex( gIdx0, gIdx1 );

          output << ( lIdx1 == tpcfe::NODE_NOT_VIRTUAL ? gIdx0 : gIdx )  << " " << vertex << "   ";
        }

        output << endl;

      }
    }

    // cout << output.str();

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
