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

#include <tpCFELevelsets.h>

#include <iterators.h>
#include <randomGenerator.h>

namespace tpcfe {

template < class RealType >
void PoempelRemoverBase<RealType>::distanceSearch ( qc::ScalarArray<distType, qc::QC_3D> &distances ) {
  qc::BitArray<qc::QC_3D> insidePoints ( distances.getNumX(), distances.getNumY(), distances.getNumZ() );
  insidePoints.setAll ( true );
  this->distanceSearch ( distances, insidePoints );
}


template < class RealType >
void PoempelRemoverBase<RealType>::distanceSearch ( qc::ScalarArray<distType, qc::QC_3D> &distances, const qc::BitArray<qc::QC_3D> &insidePoints ) {
  bool change_made = true;
  distType step = 0;

  const int
  nx = distances.getNumX(),
  ny = distances.getNumY(),
  nz = distances.getNumZ();

#ifdef VERBOSE
  cerr << "PoempelRemover: Distance search ";
#endif

  do {
    change_made = false;

    for ( qc::RectangularIterator<qc::QC_3D> bit ( distances ); bit.notAtEnd(); ++bit ) {
      const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];

      if ( distances.get ( *bit ) == step ) {

        // need to check all 15 AFE directions, not only along cube edges!
        for ( int dir = 0; dir < NUM_NEIGHBORS; ++dir ) {
          qc::CoordType nbp ( i + CFELookup<RealType>::_tetNbOffsets[dir][0], j + CFELookup<RealType>::_tetNbOffsets[dir][1], k + CFELookup<RealType>::_tetNbOffsets[dir][2] ); // indices of neighboring point

          if ( ( nbp[0] >= 0 && nbp[0] < nx && nbp[1] >= 0 && nbp[1] < ny && nbp[2] >= 0 && nbp[2] < nz ) && insidePoints.get ( nbp ) && ( distances.get ( nbp ) == DISTANCE_UNSET ) ) {
            distances.set ( nbp, step + 1 );
            change_made |= true;
          }
        }

      }
    }

    ++step;
#ifdef VERBOSE
    cerr << ".";
#endif
  } while ( change_made );

  for ( int i = 0; i < distances.size(); ++i ) {
    if ( distances[i] == DISTANCE_UNSET ) {
      distances[i] = DISTANCE_FAR ;
    }
  }

#ifdef VERBOSE
  cerr << " distance search done." << endl;
#endif
}


template< typename RealType >
void TopBottomPoempelRemover<RealType>::applySingle ( qc::ScalarArray<RealType, qc::QC_3D> &levelset ) {
  const int
  nx = levelset.getNumX(),
  ny = levelset.getNumY(),
  nz = levelset.getNumZ();

  qc::BitArray<qc::QC_3D> point_interesting ( nx, ny, nz );
  for ( int i = 0; i < levelset.size(); ++i ) {
    point_interesting.set ( i, ( levelset[i] < 0.0 ) );
  }

  qc::ScalarArray<distType, qc::QC_3D> connectivity ( nx, ny, nz );
  connectivity.setAll ( this->DISTANCE_UNSET );

  for ( int i = 0; i < nx; ++i ) {
    for ( int j = 0; j < ny; ++j ) {
      if ( point_interesting.get ( i, j, 0 ) ) {
        connectivity.set ( i, j, 0,   0 );
      }
      if ( point_interesting.get ( i, j, nz - 1 ) ) {
        connectivity.set ( i, j, nz - 1,   0 );
      }
    }
  }

  this->distanceSearch ( connectivity, point_interesting );

  // remove components not connected to top / bottom plate

  aol::ProgressBar<> pb ( "Removing poempels" );
  pb.start ( nx * ny * nz );

  for ( qc::RectangularIterator<qc::QC_3D> it ( levelset ); !it.atEnd(); ++it ) {

    // this does not shift the interface because nodes modified here are not close to the interface (else they'd be connected)

    if ( ( connectivity.get ( *it ) == this->DISTANCE_FAR ) && ( levelset.get ( *it ) < 0 ) ) {
      levelset.set ( *it,  this->levelsetOutside );
#ifdef VERBOSE
      cerr << *it  << " detected as part of poempel, removed." << endl;
#endif
    }

  }

  pb.finish();
}


template< typename RealType >
void AllFacePoempelRemover<RealType>::applySingle ( qc::ScalarArray<RealType, qc::QC_3D> &levelset ) {
  const aol::Vec3<int> size = levelset.getSize();

  for ( short face = 0; face < 3; ++face ) {
    for ( short boundaryPos = 0; boundaryPos < size[face]; boundaryPos += ( size[face] - 1 ) ) {

      qc::BitArray<qc::QC_3D> point_interesting ( size );
      for ( int i = 0; i < levelset.size(); ++i ) {
        point_interesting.set ( i, ( levelset[i] < 0.0 ) );
      }

      qc::ScalarArray<distType, qc::QC_3D> connectivity ( size );
      connectivity.setAll ( this->DISTANCE_UNSET );

      for ( int i = 0; i < size[ ( face+1 ) %3]; ++i ) {
        for ( int j = 0; j < size[ ( face+2 ) %3]; ++j ) {
          qc::CoordType position;
          position [ ( face + 0 ) % 3 ] = boundaryPos;
          position [ ( face + 1 ) % 3 ] = i;
          position [ ( face + 2 ) % 3 ] = j;

          if ( point_interesting.get ( position ) ) {
            connectivity.set ( position,   0 );
          }
        }
      }

      this->distanceSearch ( connectivity, point_interesting );

      // remove components not connected to top / bottom plate

      aol::ProgressBar<> pb ( "Removing face-poempels" );
      pb.start ( size[0]*size[1]*size[2] );

      for ( qc::RectangularIterator<qc::QC_3D> it ( levelset ); it.notAtEnd(); ++it ) {
        pb++;
        // this does not shift the interface because nodes modified here are not close to the interface (else they'd be connected)

        if ( ( connectivity.get ( *it ) == this->DISTANCE_FAR ) && ( levelset.get ( *it ) < 0 ) ) {
          levelset.set ( *it,  this->levelsetOutside );
#ifdef VERBOSE
          cerr << *it  << " detected as part of poempel, removed." << endl;
#endif
        }

      }

      pb.finish();

    }
  }
}


template< typename RealType >
void SmallComponentsPoempelRemover<RealType>::applySingle ( qc::ScalarArray<RealType, qc::QC_3D> &levelset ) {
  applySingleWithMinSize ( levelset, levelset.size() );
}


template< typename RealType >
void SmallComponentsPoempelRemover<RealType>::applySingleWithMinSize ( qc::ScalarArray<RealType, qc::QC_3D> &levelset, const int MinSize ) {
  const int nx = levelset.getNumX(), ny = levelset.getNumY(), nz = levelset.getNumZ();

  qc::BitArray<qc::QC_3D> point_interesting ( nx, ny, nz );
  for ( int i = 0; i < levelset.size(); ++i ) {
    point_interesting.set ( i, ( levelset[i] < 0.0 ) );
  }

  qc::BitArray<qc::QC_3D> pointStillToBeTreated ( nx, ny, nz );
  pointStillToBeTreated.setAll ( true );

  qc::ScalarArray< distType, qc::QC_3D > largestConnectivityComponent ( nx, ny, nz );
  int largestConnectivityComponentSize = -1;

  int pointNo = 0;

  while ( pointNo < nx * ny * nz ) {
    if ( point_interesting.get ( pointNo ) && pointStillToBeTreated.get ( pointNo ) ) {

      qc::ScalarArray<distType, qc::QC_3D> conn ( nx, ny, nz );
      conn.setAll ( this->DISTANCE_UNSET );
      conn.set ( pointNo, 0 );

      this->distanceSearch ( conn, point_interesting );

      int pointCounter = 0;
      for ( int i = 0; i < conn.size(); ++i ) {
        pointStillToBeTreated.set ( i, pointStillToBeTreated.get ( i ) && conn[i] == this->DISTANCE_FAR ); // untreated points to which we have not yet found a finite distance still to be treated
        if ( ( conn[i] >= 0 ) && ( conn[i] != this->DISTANCE_FAR ) ) {
          ++pointCounter;
        }
      }
#ifdef VERBOSE
      cerr << aol::intFormat ( pointCounter ) << " ";
#endif
      if ( pointCounter > largestConnectivityComponentSize ) {
#ifdef VERBOSE
        cerr << "is largest so far";
#endif
        largestConnectivityComponent = conn;
        largestConnectivityComponentSize = pointCounter;
      }
#ifdef VERBOSE
      cerr << endl;
#endif
    }

    ++pointNo;

    if ( largestConnectivityComponentSize >= MinSize )
      break;
  }

#ifdef VERBOSE
  aol::ProgressBar<> pb ( "Removing poempels" );
  pb.start ( nx * ny * nz );
#endif
  for ( qc::RectangularIterator<qc::QC_3D> it ( levelset ); !it.atEnd(); ++it ) {

    // this does not shift the interface because nodes modified here are not close to the interface (else they'd be connected)

    if ( ( largestConnectivityComponent.get ( *it ) == this->DISTANCE_FAR ) && ( levelset.get ( *it ) < 0 ) ) {
      levelset.set ( *it,  this->levelsetOutside );
#ifdef VERBOSE
      cerr << *it  << " detected as part of poempel, removed." << endl;
#endif
    }

  }
#ifdef VERBOSE
  pb.finish();
#endif
}



template class PoempelRemoverBase<float>;
template class PoempelRemoverBase<double>;
template class PoempelRemoverBase<long double>;

template class TopBottomPoempelRemover<float>;
template class TopBottomPoempelRemover<double>;
template class TopBottomPoempelRemover<long double>;

template class AllFacePoempelRemover<float>;
template class AllFacePoempelRemover<double>;
template class AllFacePoempelRemover<long double>;

template class SmallComponentsPoempelRemover<float>;
template class SmallComponentsPoempelRemover<double>;
template class SmallComponentsPoempelRemover<long double>;


// end namespace
}
