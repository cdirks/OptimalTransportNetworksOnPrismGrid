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

#ifndef __TPCFELEVELSETS_H
#define __TPCFELEVELSETS_H

#include <tpCFEGrid.h>

#include <scalarArray.h>
#include <multiArray.h>

// maybe move to separate file (deformations?)
namespace tpcfe {

//! Set an shift in direction dir as deformation
template< typename GridType, typename DataType >
void setDirShift ( const GridType &grid, qc::MultiArray<DataType, 3> &Dirichlet_BC, qc::BitArray<qc::QC_3D> &Dirichlet_Mask, const int dir, const DataType value ) {
  Dirichlet_BC.setZero();
  Dirichlet_Mask.setAll ( false );

  for ( int i = 0; i < grid.getNumX(); ++i ) {
    for ( int j = 0; j < grid.getNumY(); ++j ) {
      Dirichlet_Mask.set ( i, j, 0, true ); // set all nodes in top and bottom plate to Dirichet to make sure they are Dirichlet on coarsened grids.
      if ( grid.isDomainNode ( aol::Vec3<short> ( i, j, 0 ) ) ) {
        Dirichlet_BC.comp ( dir ).set ( i, j, 0, 0.0 );
      }

      Dirichlet_Mask.set ( i, j, grid.getNumZ() - 1, true );
      if ( grid.isDomainNode ( aol::Vec3<short> ( i, j, grid.getNumZ() - 1 ) ) ) {
        Dirichlet_BC.comp ( dir ).set ( i, j, grid.getNumZ() - 1, value );
      }

    }
  }

}


//! Set an x-shift as deformation (shearing)
template< typename GridType, typename DataType >
void setXShift ( const GridType &grid, qc::MultiArray<DataType, 3> &Dirichlet_BC, qc::BitArray<qc::QC_3D> &Dirichlet_Mask, const DataType value ) {
  setDirShift ( grid, Dirichlet_BC, Dirichlet_Mask, 0, value );
}


//! Set an z-shift as deformation (compression)
template< typename GridType, typename DataType >
void setZShift ( const GridType &grid, qc::MultiArray<DataType, 3> &Dirichlet_BC, qc::BitArray<qc::QC_3D> &Dirichlet_Mask, const DataType value ) {
  setDirShift ( grid, Dirichlet_BC, Dirichlet_Mask, 2, value );
}


template< typename GridType, typename DataType >
void setCompression ( const GridType &grid, qc::MultiArray<DataType, 3> &Dirichlet_BC, qc::BitArray<qc::QC_3D> &Dirichlet_Mask, const char dir, const DataType value ) {
  set_general_shrearing ( grid, Dirichlet_BC, Dirichlet_Mask, dir, dir, value );
}


template< typename GridType, typename DataType >
void setGeneralShearing ( const GridType &grid, qc::MultiArray<DataType, 3> &Dirichlet_BC, qc::BitArray<qc::QC_3D> &Dirichlet_Mask, const char fix_dir, const char shift_dir, const DataType value ) {
  Dirichlet_BC.setZero();
  Dirichlet_Mask.setAll ( false );

  qc::GridSize<qc::QC_3D> size = grid.getSize();

  for ( int i = 0; i < size[ ( fix_dir + 1 ) % 3 ]; ++i ) {
    for ( int j = 0; j < size[ ( fix_dir + 2 ) % 3 ]; ++j ) {
      aol::Vec3<short> pos1, pos2;

      switch ( fix_dir ) {
        case 0:
          pos1.set ( 0,           i,           j           );
          pos2.set ( size[0] - 1, i,           j           );
          break;
        case 1:
          pos1.set ( j,           0,           i           );
          pos2.set ( j,           size[1] - 1, i           );
          break;
        case 2:
          pos1.set ( i,           j,           0           );
          pos2.set ( i,           j,           size[2] - 1 );
          break;
        default:
          throw aol::Exception ( "tpcfe::setGeneralShearing: Illegal compression/shearing direction specified", __FILE__, __LINE__ );
      }

      Dirichlet_Mask.set ( pos1, true );
      Dirichlet_Mask.set ( pos2, true );

      if ( grid.isDomainNode ( pos2 ) ) {
        Dirichlet_BC.comp ( shift_dir ).set ( pos2, value );
      }

    }
  }

}


template< typename GridType, typename DataType >
void setTorsion ( const GridType &grid, qc::MultiArray<DataType, 3> &Dirichlet_BC, qc::BitArray<qc::QC_3D> &Dirichlet_Mask, const DataType angle ) {
  Dirichlet_BC.setZero();
  Dirichlet_Mask.setAll ( false );

  const DataType
  co = cos ( angle ),
  si = sin ( angle );

  const int mZ = grid.getNumZ() - 1;

  for ( int i = 0; i < grid.getNumX(); ++i ) { // may also work for non-cube domains
    for ( int j = 0; j < grid.getNumY(); ++j ) {

      Dirichlet_Mask.set ( i, j, 0, true );
      Dirichlet_Mask.set ( i, j, mZ, true );

      if ( grid.isDomainNode ( aol::Vec3<short> ( i, j, mZ ) ) ) {

        const DataType
        x = i / ( grid.getNumX() - 1.0 ) - 0.5,
        y = j / ( grid.getNumY() - 1.0 ) - 0.5;

        Dirichlet_BC.comp ( 0 ).set ( i, j, mZ, - (  co * x + si * y - x ) );
        Dirichlet_BC.comp ( 1 ).set ( i, j, mZ, - ( -si * x + co * y - y ) );
      }
    }
  }
}

}


namespace tpcfe {

/** A abstract basis class to remove certain components from a volume
 *  (0-sublevelset) that, being disconnected from boundary conditions,
 *  make the problem underdetermined.
 *  Implement apply in derived classes.
 *  \author Schwen
 */

template < class RealType >
class PoempelRemoverBase {
protected:
  typedef signed short distType;

  enum {
    DISTANCE_FAR         = 1 << 14, // should be big enough ...
    DISTANCE_IRRELEVANT  = -2,
    DISTANCE_UNSET       = -10
  };

  RealType levelsetOutside;

public:
  PoempelRemoverBase ( RealType _levelsetOutside = 1.0 ) : levelsetOutside ( _levelsetOutside ) {
    // nothing else to be done ...
  }

  virtual ~PoempelRemoverBase ( ) {
  }

public:
  virtual void applySingle ( qc::ScalarArray<RealType, qc::QC_3D> &levelset ) = 0;

protected:
  /** Do breadth-first search of interesting points
   *  \param distances must contain 0 for point to start from (or multiple points to start from), DISTANCE_UNSET at all other points
   *  \param insidePoints must contain true for all points that are considered "inside", false for all other nodes. Both arrays must be of same size
   */
  void distanceSearch ( qc::ScalarArray<distType, qc::QC_3D> &distances, const qc::BitArray<qc::QC_3D> &insidePoints );

  //! distanceSearch for full connectivity, i. e. all points inside
  void distanceSearch ( qc::ScalarArray<distType, qc::QC_3D> &distances );
};


/** Remove all parts of the structure that are not connected to the top or bottom plate. */
template < class RealType >
class TopBottomPoempelRemover : public PoempelRemoverBase<RealType> {
  typedef typename PoempelRemoverBase<RealType>::distType distType;

public:
  TopBottomPoempelRemover ( RealType _levelsetOutside = 1.0 ) : PoempelRemoverBase<RealType> ( _levelsetOutside ) {}

  void applySingle ( qc::ScalarArray<RealType, qc::QC_3D> &levelset );
};


/** Remove all parts of the structure that are not connected to all six plates */
template < class RealType >
class AllFacePoempelRemover : public PoempelRemoverBase<RealType> {
  typedef typename PoempelRemoverBase<RealType>::distType distType;

public:
  AllFacePoempelRemover ( RealType _levelsetOutside = 1.0 ) : PoempelRemoverBase<RealType> ( _levelsetOutside ) {}

  void applySingle ( qc::ScalarArray<RealType, qc::QC_3D> &levelset );
};

/** Removes all components that are not one biggest connected component */
template < class RealType >
class SmallComponentsPoempelRemover : public PoempelRemoverBase<RealType> {
  typedef typename PoempelRemoverBase<RealType>::distType distType;

public:
  SmallComponentsPoempelRemover ( RealType _levelsetOutside = 1.0 ) : PoempelRemoverBase<RealType> ( _levelsetOutside ) {}

  void applySingle ( qc::ScalarArray<RealType, qc::QC_3D> &levelset );

  void applySingleWithMinSize ( qc::ScalarArray<RealType, qc::QC_3D> &levelset, const int MinSize );
};

} // end namespace

#endif
