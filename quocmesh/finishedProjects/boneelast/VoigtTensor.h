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

#ifndef __VOIGTTENSOR_H
#define __VOIGTTENSOR_H

#include <vec.h>
#include <matrix.h>
#include <scalarArray.h>
#include <arrayExtensions.h>

// "really-Voigt numbering: xx, yy, zz, yz, xz, xy
static const int Voigtmap[3][3] = { { 0, 5, 4 },
                                    { 5, 1, 3},
                                    { 4, 3, 2 } };


template< typename RealType >
class VoigtElastOp : public aol::FullMatrix<RealType> {
public:
  VoigtElastOp ( ) : aol::FullMatrix<RealType> ( 6, 6 ) {
  }

  // Constructor in the isotropic case
  VoigtElastOp ( const RealType E, const RealType nu, const bool /*avoid calling constructor with two integers*/ ) : aol::FullMatrix<RealType> ( 6, 6 ) {
    const RealType
      denom = ( 1 + nu ) * ( 1 - 2 * nu ),
      upperDiag = E * ( 1 - nu ) / denom,
      upperOffD = E * nu / denom,
      lowerDiag = E / ( 2 * ( 1 + nu ) ); //  think about the strange factor sqrt(2)!

    this->set ( 0, 0, upperDiag );    this->set ( 0, 1, upperOffD );    this->set ( 0, 2, upperOffD );
    this->set ( 1, 0, upperOffD );    this->set ( 1, 1, upperDiag );    this->set ( 1, 2, upperOffD );
    this->set ( 2, 0, upperOffD );    this->set ( 2, 1, upperOffD );    this->set ( 2, 2, upperDiag );
    this->set ( 3, 3, lowerDiag );    this->set ( 4, 4, lowerDiag );    this->set ( 5, 5, lowerDiag );
  }

  VoigtElastOp ( const char* const filename ) : aol::FullMatrix<RealType> ( 6, 6 ) {
    ifstream inpu ( filename );

    if ( !inpu.good() )
      throw aol::Exception ( "VoigtElastOp: could not open file", __FILE__, __LINE__ );

    for ( int i = 0; i < 6; ++i ) {
      for ( int j = 0; j < 6; ++j ) {
        RealType value;
        inpu >> value;
        this->set ( i, j, value );
      }
    }
  }

  void get4thOrderTensor ( aol::Matrix33<RealType> fourthOrderTensor[3][3] ) { // yes, this is call-by-reference
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        for ( int k = 0; k < 3; ++k )
          for ( int l = 0; l < 3; ++l )
            fourthOrderTensor[i][j][k][l] = this->get( Voigtmap[i][j] , Voigtmap[k][l] );

  }

  void printTensor ( ostream &out = cout ) const {
    for ( int i = 0; i < 6; ++i ) {
      for ( int j = 0; j < 6; ++j ) {
        out << aol::longScientificFormat ( this->get ( i, j ) ) << " ";
      }
      out << endl;
    }
  }

  void printLatexTensor ( ostream &out = cout ) const {
    for ( int i = 0; i < 6; ++i ) {
      for ( int j = 0; j < 6; ++j ) {
        out << aol::mixedFormat ( this->get ( i, j ) ) << " " << ( j != 5 ?  "&" : "\\\\" );
      }
      out << endl;
    }
  }

  RealType FrobeniusNormOrthoEntries ( ) const {
    return ( sqrt ( aol::Sqr ( this->get( 0, 0 ) ) + aol::Sqr ( this->get( 0, 1 ) ) + aol::Sqr ( this->get( 0, 2 ) ) +
                    aol::Sqr ( this->get( 1, 0 ) ) + aol::Sqr ( this->get( 1, 1 ) ) + aol::Sqr ( this->get( 1, 2 ) ) +
                    aol::Sqr ( this->get( 2, 0 ) ) + aol::Sqr ( this->get( 2, 1 ) ) + aol::Sqr ( this->get( 2, 2 ) ) +
                    4 * aol::Sqr ( this->get( 3, 3 ) ) + 4 * aol::Sqr ( this->get( 4, 4 ) ) + 4 * aol::Sqr ( this->get( 5, 5 ) ) ) );
  }

  void restrictToOrtho ( ) {
    for ( short i = 0; i < 3; ++i ) {
      for ( short j = 3; j < 6; ++j ) {
        this->set ( i, j, 0 );
        this->set ( j, i, 0 );
      }
    }
    for ( short i = 3; i < 6; ++i ) {
      for ( short j = 3; j < 6; ++j ) {
        if ( i != j ) {
          this->set ( i, j, 0 );
        }
      }
    }
  }

};

// forward rotation by three angles in yz, xz and xy plane (around x, y and z axis); angles in rad
template< typename RealType >
void getRotationMatrix ( const aol::Vec3<RealType> &angles, aol::Matrix33<RealType> &rotMat ) {
  aol::Matrix33<RealType> rotXY, rotXZ, rotYZ;
  rotYZ.setRotationAboutX ( angles[0] );
  rotXZ.setRotationAboutY ( angles[1] );
  rotXY.setRotationAboutZ ( angles[2] );

  rotMat.setIdentity();
  rotMat *= rotYZ;
  rotMat *= rotXZ;
  rotMat *= rotXY;
}

// backward rotation in opposite order; angles in rad
template< typename RealType >
void getBackRotationMatrix ( const aol::Vec3<RealType> &angles, aol::Matrix33<RealType> &rotMat ) {
  aol::Matrix33<RealType> rotXY, rotXZ, rotYZ;
  rotYZ.setRotationAboutX ( -angles[0] );
  rotXZ.setRotationAboutY ( -angles[1] );
  rotXY.setRotationAboutZ ( -angles[2] );

  rotMat.setIdentity();
  rotMat *= rotXY;
  rotMat *= rotXZ;
  rotMat *= rotYZ;
}


template< typename RealType >
void applyRotation ( const aol::Matrix33<RealType> &rotation, const VoigtElastOp<RealType> &arg, VoigtElastOp<RealType> &dest ) {

 aol::Matrix33<RealType> argNonV[3][3], destNonV[3][3];
  VoigtElastOp<RealType> denom;

  // Voigt expand
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int k = 0; k < 3; ++k )
        for ( int l = 0; l < 3; ++l )
          argNonV[i][j][k][l] = arg.get( Voigtmap[i][j] , Voigtmap[k][l] );

  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int k = 0; k < 3; ++k )
        for ( int l = 0; l < 3; ++l )
          for ( int m = 0; m < 3; ++m )
            for ( int n = 0; n < 3; ++n )
              for ( int o = 0; o < 3; ++o )
                for ( int p = 0; p < 3; ++p )
                  destNonV[i][j][k][l] += rotation[i][m] * rotation[j][n] * rotation[k][o] * rotation[l][p] * argNonV[m][n][o][p];

  // Voigt reduce
  for ( int i = 0; i < 3; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      for ( int k = 0; k < 3; ++k ) {
        for ( int l = 0; l < 3; ++l ) {
          dest.add ( Voigtmap[i][j], Voigtmap[k][l],  destNonV[i][j][k][l] );
          denom.add ( Voigtmap[i][j], Voigtmap[k][l], 1.0 );
        }
      }
    }
  }

  // averageing:
  for ( int i = 0; i < 6; ++i )
    for ( int j = 0; j < 6; ++j )
      dest.set ( i, j, dest.get ( i, j ) / denom.get( i, j ) );

}


// rotate 6x6 Voigt matrix arg to dest by angles; angles in rad
template< typename RealType >
void applyRotation ( const aol::Vec3<RealType> &angles, const VoigtElastOp<RealType> &arg, VoigtElastOp<RealType> &dest, const bool inverseRotation ) {
  dest.setZero();

  aol::Matrix33<RealType> rotation;
  if ( !inverseRotation ) {
    getRotationMatrix ( angles, rotation );
  } else {
    getBackRotationMatrix ( angles, rotation );
  }
  applyRotation ( rotation, arg, dest );

}

// apply forward rotation of tensor by angles; angles in rad
template< typename RealType >
void applyForwardRotation ( const aol::Vec3<RealType> &angles, const VoigtElastOp<RealType> &arg, VoigtElastOp<RealType> &dest ) {
  applyRotation ( angles, arg, dest, false );
}

// apply backward rotation of tensor by angles; angles in rad
template< typename RealType >
void applyBackwardRotation ( const aol::Vec3<RealType> &angles, const VoigtElastOp<RealType> &arg, VoigtElastOp<RealType> &dest ) {
  applyRotation ( angles, arg, dest, true );
}

#endif
