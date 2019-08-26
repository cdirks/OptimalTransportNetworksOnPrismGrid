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

// load a dataset, print all (regular and virtual) points and tetra to stdout (for piping into file)

#include <tpCFEGrid.h>

int main ( int, char** ) {
  try {

    typedef float RealType;
    typedef tpcfe::CFEGrid< RealType, tpcfe::CFE_CD> GridType;
    typedef aol::Vec< 4, uint64_t > TetraType;

    const int depth = 8;

    GridType grid ( depth );

    qc::ScalarArray<RealType, qc::QC_3D> levelset( grid );

    levelset.load( "../boneelast/datasets/961_T1Po1_5.dat.bz2" );

    grid.setDomainFrom ( levelset );
    grid.detectVirtualNodes();


    std::map< uint64_t, aol::Vec3<RealType> > vertex_map;
    std::vector< TetraType >   tetra_vector;

    aol::ProgressBar<> pb ( "Finding tetra" );
    pb.start ( grid.getNumberOfElements() );

    const qc::GridSize<qc::QC_3D> gridSize ( grid );
    for ( GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
      tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

      el.computeAssembleData ( grid );

      for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el, -1 ); tit.notAtEnd(); ++tit ) {

        TetraType tetra;
        // do something with the tetra:

        // put coordinates of vertices in vertex_map
        for ( int v = 0; v < 4; ++v ) {
          aol::Vec3<RealType> vertex;
          tit->computeGlobalCoordinate( vertex, el, v );

          const int
            lIdx0 = ( *tit ) ( v, 0 ),
            lIdx1 = ( *tit ) ( v, 1 ),
            gIdx0 = el.globalIndex(lIdx0),
            gIdx1 = ( lIdx1 == 11 ? 0 : el.globalIndex(lIdx1) );
          const uint64_t  gIdx  = tpcfe::CFEVirtualNode< RealType, tpcfe::CFE_CD, RealType >::mapIndex( gIdx0, gIdx1 );

          tetra[v] = gIdx;

          vertex_map[ gIdx ] = vertex;

        }

        tetra_vector.push_back( tetra );

      }
    }

    cerr << endl;

    std::map< uint64_t, int >  vertex_indexmap;
    std::vector< aol::Vec3<RealType> > vertex_vector;

    {
      // remap for printing
      int i = 0;
      for ( std::map< uint64_t, aol::Vec3<RealType> >::iterator it = vertex_map.begin(); it != vertex_map.end(); ++it ) {
        vertex_indexmap[ it->first ] = i;
        ++i;
        vertex_vector.push_back ( it->second );
      }
    }

    cout << "Number of vertices: " << vertex_vector.size() << endl;
    for ( std::vector< aol::Vec3<RealType> >::iterator it = vertex_vector.begin(); it != vertex_vector.end(); ++it ) {
      for( int i = 0; i < 3; ++i ) {
        cout << (*it)[i] << " ";
      }
      cout << endl;
    }

    cout << "Number of Tetrahedra: " << tetra_vector.size() << endl;
    for ( std::vector< TetraType >::iterator it = tetra_vector.begin(); it != tetra_vector.end(); ++it ) {
      for( int i = 0; i < 4; ++i ) {
        cout << vertex_indexmap[ (*it)[i] ] << " ";
      }
      cout << endl;
    }


  } catch ( aol::Exception &exc ) {

    exc.dump();
    return ( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );
}


