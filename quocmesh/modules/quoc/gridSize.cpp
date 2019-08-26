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

#include <gridSize.h>

namespace qc {

// --------------------------------------------------------------------------
template <>
GridSize<QC_2D>::GridSize ( const GridSize<QC_2D> &Size ) : CoordType() {
  for (int i = 0; i < 2; ++i)
    (*this)[i] = Size[i];
  CoordType::operator [] (2) = 1;
}

// --------------------------------------------------------------------------
template <>
GridSize<QC_3D>::GridSize ( const GridSize<QC_3D> &Size ) : CoordType() {
  for (int i = 0; i < 3; ++i)
    (*this)[i] = Size[i];
}

// ----------------------------------------------------------------------------------------------
template <>
GridSize<QC_2D>::GridSize ( GridSize<QC_2D>::IndexType NumX,
                            GridSize<QC_2D>::IndexType NumY ) {
  (*this)[0] = NumX;
  (*this)[1] = NumY;
  for (int i = 2; i < dim; ++i)
    CoordType::operator [] (i) = 1;
}
// ----------------------------------------------------------------------------------------------
template <>
GridSize<QC_3D>::GridSize ( GridSize<QC_3D>::IndexType NumX,
                            GridSize<QC_3D>::IndexType NumY,
                            GridSize<QC_3D>::IndexType NumZ ) {
  (*this)[0] = NumX;
  (*this)[1] = NumY;
  (*this)[2] = NumZ;
  for (int i = 3; i < dim; ++i)
    CoordType::operator [] (i) = 1;
}
// ----------------------------------------------------------------------------------------------
template <>
int GridSize<QC_1D>::getNumberOfNodes () const {
  return getNumX();
}
// ----------------------------------------------------------------------------------------------
template <>
int GridSize<QC_2D>::getNumberOfNodes () const {
  return getNumX() * getNumY();
}
// ----------------------------------------------------------------------------------------------
template <>
int GridSize<QC_3D>::getNumberOfNodes () const {
  return getNumX() * getNumY() * getNumZ();
}
// ----------------------------------------------------------------------------------------------
template <Dimension _Dim>
int GridSize<_Dim>::getNumberOfDofs () const {
  return GridSize<_Dim>::getNumberOfNodes ();
}
// ----------------------------------------------------------------------------------------------
template <Dimension _Dim>
typename GridSize<_Dim>::IndexType GridSize<_Dim>::getNumX() const {
  return (*this)[0];
}
// ----------------------------------------------------------------------------------------------
template <Dimension _Dim>
typename GridSize<_Dim>::IndexType GridSize<_Dim>::getNumY() const {
  if ( Dim < QC_2D )
    throw aol::DimensionMismatchException ( "GridSize<d>::getNumY() called "
                      "for Dimension d < 2.", __FILE__, __LINE__ );
  return (*this)[1];
}
// ----------------------------------------------------------------------------------------------
template <Dimension _Dim>
typename GridSize<_Dim>::IndexType GridSize<_Dim>::getNumZ() const {
  if ( Dim < QC_3D )
    throw aol::DimensionMismatchException ( "GridSize<d>::getNumZ() called "
                      "for Dimension d < 3.", __FILE__, __LINE__ );
  return (*this)[2];
}
// ----------------------------------------------------------------------------------------------
template <Dimension _Dim>
typename GridSize<_Dim>::IndexType GridSize<_Dim>::operator [] (int i) const {
#ifdef BOUNDS_CHECK
  if (i < 0 || i >= _Dim)
    throw aol::DimensionMismatchException ( "GridSize::operator[] called "
                      "with argument >= Dim.", __FILE__, __LINE__ );
#endif
  return CoordType::operator [] (i);
}
// ----------------------------------------------------------------------------------------------
template <Dimension _Dim>
typename GridSize<_Dim>::IndexType & GridSize<_Dim>::operator [] (int i) {
#ifdef BOUNDS_CHECK
  if (i < 0 || i >= _Dim)
    throw aol::DimensionMismatchException ( "GridSize::operator[] called "
                      "with argument >= Dim.", __FILE__, __LINE__ );
#endif
  return CoordType::operator [] (i);
}
// ==============================================================================================

// explicit instantiations

template class GridSize<QC_1D>;
template class GridSize<QC_2D>;
template class GridSize<QC_3D>;

// ----------------------------------------------------------------------------------------------

} // end namespace qc
