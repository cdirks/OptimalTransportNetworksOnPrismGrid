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

#include <narrow.h>

namespace nb {

template <typename _FullGridType>
void NarrowBandGrid<_FullGridType, qc::QC_2D>::extractEdges( ) {
  _edges.clear();
  // lookUp: the first two numbers per block are the (x,y) offsets of the left lower node of
  // the adjacent element. The last two are the nodenumbers of the corresponding edge.
  int lookUp[4][4] = {  { 1, 0, 1, 3 },     // right
                        { 0, 1, 3, 2 },     // above
                        { -1, 0, 2, 0 },    // left
                        { 0, -1, 0, 1 } };  // below
  // offsets contain the offset of the boundary node to the lower left node of the current element
  int offsets[4][2] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
  for ( OldFullElementIterator it = begin(); it != end(); ++it ) {
    for ( int i = 0; i < 4; i++ ) {
      ElementType el ( it->x() + lookUp[i][0], it->y() + lookUp[i][1], 0, it->level() );
      if ( !this->exists ( el ) ) {
        _edges.push_back ( Edge ( *it, lookUp[i][2], lookUp[i][3] ) );
        _boundaryNodes.push_back ( qc::CoordType ( el.x() + offsets[i][0], el.y() + offsets[i][1], 0 ) );
      }
    }
  }
}

template <typename _FullGridType>
void NarrowBandGrid<_FullGridType, qc::QC_3D>::extractEdges( ) {
  _edges.clear();
  int lookUp[6][7] = {  {  1,  0,  0, 1, 3, 5, 7 },
                        {  0,  1,  0, 2, 3, 6, 7 },
                        {  0,  0,  1, 4, 5, 6, 7 },
                        { -1,  0,  0, 0, 2, 4, 6 },
                        {  0, -1,  0, 0, 1, 4, 5 },
                        {  0,  0, -1, 0, 1, 2, 3 } };
  for ( OldFullElementIterator it = begin(); it != end(); ++it ) {
    for ( int i = 0; i < 6; i++ ) {
      ElementType el ( it->x() + lookUp[i][0],
                       it->y() + lookUp[i][1],
                       it->z() + lookUp[i][2], it->level() );
      if ( !this->exists ( el ) ) {

        // distinguish between the two edgesets here!!!

        _edges.push_back ( Face ( *it, lookUp[i][3], lookUp[i][4], lookUp[i][5], lookUp[i][6] ) );
      }
    }
  }
}

template <typename _RealType>
void NarrowBandMaskedGrid<_RealType, qc::QC_3D>::extractEdges( ) {
  qc::ScalarArray<int, qc::QC_3D> bNodeImage ( *this );
  bNodeImage.setZero();

  qc::BitArray<qc::QC_3D> nodeImage ( qc::GridSize<qc::QC_3D>::createFrom ( *this ) );

  _edges.clear();
  int lookUp[6][7] = {  {  1,  0,  0, 1, 3, 5, 7 },
                        {  0,  1,  0, 2, 3, 6, 7 },
                        {  0,  0,  1, 4, 5, 6, 7 },
                        { -1,  0,  0, 0, 2, 4, 6 },
                        {  0, -1,  0, 0, 1, 4, 5 },
                        {  0,  0, -1, 0, 1, 2, 3 } };
  int edgelookUp[6][4][3] = { { {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}},
                              { {0, 1, 0}, {0, 1, 1}, {1, 1, 0}, {1, 1, 1}},
                              { {0, 0, 1}, {0, 1, 1}, {1, 0 , 1}, {1, 1, 1}},
                              { {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}},
                              { {0, 0, 0}, {0, 0, 1}, {1, 0, 0}, {1, 0, 1}},
                              { {0, 0, 0}, {0, 1, 0}, {1, 0 , 0}, {1, 1, 0}}};

  setIteratorMode ( EXT );

  for ( iterator it = begin(); it != end(); ++it ) {
    for ( int i = 0; i < 6; i++ ) {
      qc::Element el ( it->x() + lookUp[i][0],
                       it->y() + lookUp[i][1],
                       it->z() + lookUp[i][2], it->level() );
      if ( !exists ( el ) ) {
        bool Dflag = true;
        for ( int j = 0; j < 4; j++ ) {
          qc::CoordType c ( it->x() + edgelookUp[i][j][0],
                            it->y() + edgelookUp[i][j][1],
                            it->z() + edgelookUp[i][j][2] );

          if ( _mask_img.get ( c ) >= 0.0 ) {
            Dflag = false;
            if ( !bNodeImage.get ( c ) ) {
              bNodeImage.set ( c, 1 );
              _boundaryNodes.push_back ( c );
            }
          } else {
            if ( !bNodeImage.get ( c ) ) {
              bNodeImage.set ( c, 1 );
              _boundaryNodesD.push_back ( c );
            }
          }
        }
        if ( Dflag ) {
          _edgesD.push_back ( Face ( *it, lookUp[i][3], lookUp[i][4], lookUp[i][5], lookUp[i][6] ) );
        } else {
          _edges.push_back ( Face ( *it, lookUp[i][3], lookUp[i][4], lookUp[i][5], lookUp[i][6] ) );
        }
      }
    }
  }


  _hashInterior = _hash;
  _hash.erase ( _hash.begin(), _hash.end() );

  setIteratorMode ( INT );

  const int w = 1;

  for ( iterator it = begin(); it != end(); ++it ) {
    for ( int i = -w; i <= w; i++ ) {
      for ( int j = -w; j <= w; j++ ) {
        for ( int k = -w; k <= w; k++ ) {
          qc::Element el ( it->x() + i, it->y() + j, it->z() + k, it->level() );
          if ( el[0] >= 0 && el[0] <= this->getWidth() - 2 &&
               el[1] >= 0 && el[1] <= this->getWidth() - 2 &&
               el[2] >= 0 && el[2] <= this->getWidth() - 2 ) {
            if ( !exists ( el ) ) {
              _hash.insert ( el );
            }
          }
        }
      }
    }
  }


  nodeImage.setAll ( false );

  int numInnerNodes = 0;
  int numOuterNodes = 0;

  setIteratorMode ( INT );
  for ( iterator it = begin(); it != end(); ++it ) {
    for ( qc::ElementNodeIterator<qc::QC_3D> nit = this->elementNodeBegin<qc::QC_3D> ( *it ); nit != this->elementNodeEnd<qc::QC_3D> ( *it ); ++nit )  {
      if ( !nodeImage.get ( *nit ) ) { nodeImage.set ( *nit, true ); numInnerNodes++; }
    }
  }


  nodeImage.setAll ( false );
  setIteratorMode ( EXT );
  for ( iterator it = begin(); it != end(); ++it ) {
    for ( qc::ElementNodeIterator<qc::QC_3D> nit = this->elementNodeBegin<qc::QC_3D> ( *it ); nit != this->elementNodeEnd<qc::QC_3D> ( *it ); ++nit )  {
      if ( !nodeImage.get ( *nit ) ) { nodeImage.set ( *nit, true ); numOuterNodes++; }
    }
  }

  cerr << "\n\nnumInnerNode = " << numInnerNodes << "   numOuterNodes = " << numOuterNodes << endl << endl;


  setIteratorMode ( EXT );

  bNodeImage.setZero();
  for ( iterator it = begin(); it != end(); ++it ) {
    for ( int i = 0; i < 6; i++ ) {
      qc::Element el ( it->x() + lookUp[i][0],
                       it->y() + lookUp[i][1],
                       it->z() + lookUp[i][2],
                       it->level() );

      if ( !exists ( el ) ) {
        for ( int j = 0; j < 4; j++ ) {
          qc::CoordType c ( it->x() + edgelookUp[i][j][0],
                            it->y() + edgelookUp[i][j][1],
                            it->z() + edgelookUp[i][j][2] );
          if ( !bNodeImage.get ( c ) ) {
            bNodeImage.set ( c, 1 );
            _outerBoundaryNodes.push_back ( c );
          }
        }
      }
    }
  }
}

template class NarrowBandGrid<qc::GridDefinition, qc::QC_2D>;
template class NarrowBandGrid<qc::GridDefinition, qc::QC_3D>;
template class NarrowBandGrid<qc::RectangularGrid<qc::QC_2D>, qc::QC_2D>;
template class NarrowBandMaskedGrid<double, qc::QC_3D>;

}   // end namespace
