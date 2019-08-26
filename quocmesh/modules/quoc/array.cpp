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

#include <array.h>

//-----------------------------------------------------------------------------------------------------
template<typename _DataType>
void qc::Array<_DataType>::putSlice ( const qc::Comp Comp, int Index, const qc::Array<_DataType> &Src ) {
  int X, Y, Z;
  switch ( Comp ) {
    case QC_X:
      if ( this->getNumY() != Src.getNumX() || this->getNumZ() != Src.getNumY() ) {
        if ( !this->quietMode )  cerr << "getNumY() = " << this->getNumY()
          << "getNumZ() = " << this->getNumZ()
          << "Src.getNumX() = " << Src.getNumX()
          << "Src.getNumY() = " << Src.getNumY() << endl;
        throw aol::Exception ( "dimensions don't match", "qc::Array<DataType, qc::QC_3D>::getSlice" );
      }
      for ( Y = 0; Y < this->getNumY(); ++Y ) {
        for ( Z = 0; Z < this->getNumZ(); ++Z ) {
          set ( Index, Y, Z, Src.get ( Y, Z ) );
        }
      }
      break;
    case QC_Y:
      if ( this->getNumX() != Src.getNumX() || this->getNumZ() != Src.getNumY() ) {
        if ( !this->quietMode ) cerr << "getNumX() = " << this->getNumX()
          << "getNumZ() = " << this->getNumZ()
          << "Src.getNumX() = " << Src.getNumX()
          << "Src.getNumY() = " << Src.getNumY() << endl;
        throw aol::Exception ( "dimensions don't match", "qc::Array<DataType, qc::QC_3D>::getSlice" );
      }
      for ( X = 0; X < this->getNumX(); ++X ) {
        for ( Z = 0; Z < this->getNumZ(); ++Z ) {
          set ( X, Index, Z, Src.get ( X, Z ) );
        }
      }
      break;
    case qc::QC_Z:
      if ( this->getNumX() != Src.getNumX() || this->getNumY() != Src.getNumY() ) {
        if ( !this->quietMode ) cerr << "getNumX() = " << this->getNumX()
          << "getNumY() = " << this->getNumY()
          << "Src.getNumX() = " << Src.getNumX()
          << "Src.getNumY() = " << Src.getNumY() << endl;
        throw aol::Exception ( "dimensions don't match", "qc::Array<DataType, qc::QC_3D>::getSlice" );
      }
      for ( X = 0; X < this->getNumX(); ++X ) {
        for ( Y = 0; Y < this->getNumY(); ++Y ) {
          set ( X, Y, Index, Src.get ( X, Y ) );
        }
      }
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unknown qc::Comp", __FILE__, __LINE__ );
  }
}

//-----------------------------------------------------------------------------------------------------
template<typename _DataType>
void qc::Array<_DataType>::getSlice ( const qc::Comp Comp, int Index, qc::Array<_DataType> &Dest ) const {
  int X, Y, Z;
  switch ( Comp ) {
    case QC_X:
      if ( this->getNumY() != Dest.getNumX() || this->getNumZ() != Dest.getNumY() ) {
        if ( !this->quietMode ) cerr << "getNumY() = " << this->getNumY()
          << "getNumZ() = " << this->getNumZ()
          << "Dest.getNumX() = " << Dest.getNumX()
          << "Dest.getNumY() = " << Dest.getNumY() << endl;
        throw aol::Exception ( "dimensions don't match", "qc::Array<DataType, qc::QC_3D>::getSlice" );
      }
      for ( Y = 0; Y < this->getNumY(); ++Y ) {
        for ( Z = 0; Z < this->getNumZ(); ++Z ) {
          Dest.set ( Y, Z, this->get ( Index, Y, Z ) );
        }
      }
      break;
    case QC_Y:
      if ( this->getNumX() != Dest.getNumX() || this->getNumZ() != Dest.getNumY() ) {
        if ( !this->quietMode ) cerr << "getNumX() = " << this->getNumX()
          << "getNumZ() = " << this->getNumZ()
          << "Dest.getNumX() = " << Dest.getNumX()
          << "Dest.getNumY() = " << Dest.getNumY() << endl;
        throw aol::Exception ( "dimensions don't match", "qc::Array<DataType, qc::QC_3D>::getSlice" );
      }
      for ( X = 0; X < this->getNumX(); ++X ) {
        for ( Z = 0; Z < this->getNumZ(); ++Z ) {
          Dest.set ( X, Z, this->get ( X, Index, Z ) );
        }
      }
      break;
    case qc::QC_Z:
      if ( this->getNumX() != Dest.getNumX() || this->getNumY() != Dest.getNumY() ) {
        if ( !this->quietMode ) cerr << "getNumX() = " << this->getNumX()
          << "getNumY() = " << this->getNumY()
          << "Dest.getNumX() = " << Dest.getNumX()
          << "Dest.getNumY() = " << Dest.getNumY() << endl;
        throw aol::Exception ( "dimensions don't match", "qc::Array<DataType, qc::QC_3D>::getSlice" );
      }
      for ( X = 0; X < this->getNumX(); ++X ) {
        for ( Y = 0; Y < this->getNumY(); ++Y ) {
          Dest.set ( X, Y, this->get ( X, Y, Index ) );
        }
      }
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unknown qc::Comp", __FILE__, __LINE__ );
  }
}


template class qc::Array <unsigned char>;
template class qc::Array <signed char>;
template class qc::Array <short>;
template class qc::Array <unsigned short>;
template class qc::Array <int>;
template class qc::Array <unsigned int>;
template class qc::Array <float>;
template class qc::Array <double>;
template class qc::Array <long double>;

namespace qc {

template<> const float qc::Array<unsigned char >::DIFF_STD_SIGMA = 0.2;
template<> const float qc::Array<signed char   >::DIFF_STD_SIGMA = 0.2;
template<> const float qc::Array<short         >::DIFF_STD_SIGMA = 0.2;
template<> const float qc::Array<unsigned short>::DIFF_STD_SIGMA = 0.2;
template<> const float qc::Array<int           >::DIFF_STD_SIGMA = 0.2;
template<> const float qc::Array<unsigned int  >::DIFF_STD_SIGMA = 0.2;
template<> const float qc::Array<float         >::DIFF_STD_SIGMA = 0.2;
template<> const float qc::Array<double        >::DIFF_STD_SIGMA = 0.2;
template<> const float qc::Array<long double   >::DIFF_STD_SIGMA = 0.2;

}
