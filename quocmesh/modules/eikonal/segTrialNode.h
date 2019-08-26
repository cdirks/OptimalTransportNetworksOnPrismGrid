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

#ifndef __SEGTRIALNODE_H
#define __SEGTRIALNODE_H

#include <aol.h>
#include <segDefs.h>

namespace eik {

class Eikonal3d;
class MetaEstimator;

enum QC_SEED_TYPE {
  QC_POS = 1,
  QC_NEG = 2
};

inline void offsetCoords ( int Dir, int &X, int &Y, int &Z, int Step  ) {
  switch ( Dir & TN_DIR_ALL ) {
  case TN_DIR_BACK:
    if ( Dir & TN_DIR_OPP ) Z -= Step;
    else Z += Step;
    break;
  case TN_DIR_UP:
    if ( Dir & TN_DIR_OPP ) Y -= Step;
    else Y += Step;
    break;
  case TN_DIR_RIGHT:
    if ( Dir & TN_DIR_OPP ) X -= Step;
    else X += Step;
    break;
  default:
    throw aol::UnimplementedCodeException ( "Unhandled switch case!", __FILE__, __LINE__ );
  }
}

class Eikonal3d;

template <class T>
class TrialNode {
public:
  TrialNode ( int X, int Y, int Z, QC_SEED_TYPE Seed_Type = QC_POS )
      : _x ( X ), _y ( Y ), _z ( Z ),
      val ( -1.0 ),
      seed_type ( Seed_Type ),
      hangNodeLevel ( TN_HN_NONE_SEG ) {
    for ( int i = 0; i <= TN_DIR_MAX_NUM ; i++ ) {
      dirLevel[ i ] = TN_HN_NONE_SEG;
    }
    faceHangNode = TN_HN_NONE_SEG;
    edgeHangNode = TN_HN_NONE_SEG;

    // TEST
    if ( X == 84 && Y == 110 && Z == 168 ) {
      void debughalt( );
    }
  }

  TrialNode( )
      : _x ( 0 ), _y ( 0 ), _z ( 0 ),
      val ( -1.0 ),
      seed_type ( QC_POS ),
      hangNodeLevel ( TN_HN_NONE_SEG ) {
    for ( int i = 0; i <= TN_DIR_MAX_NUM; i++ ) {
      dirLevel[ i ] = TN_HN_NONE_SEG;
    }
    faceHangNode = TN_HN_NONE_SEG;
    edgeHangNode = TN_HN_NONE_SEG;

  }

  void debughalt( );

  int consistencyCheck ( int MaxLevel, MetaEstimator &estimator, Eikonal3d &segmentor );

  inline void addHangNode ( int HNValue, int HNLevel );

  int minDirLevel( ) const {
    int minLevel = dirLevel[ 0 ];
    for ( int i = 1; i <= TN_DIR_MAX_NUM; i++ ) {
      if ( dirLevel[ i ] < minLevel )
        minLevel = dirLevel[ i ];
    }
    return minLevel;
  }


  int maxDirLevel( ) const {
    int maxLevel = dirLevel[ 0 ];

    for ( int i = 1; i <= TN_DIR_MAX_NUM; i++ ) {
      if ( dirLevel[ i ] > maxLevel )
        maxLevel = dirLevel[ i ];
    }
    return maxLevel;
  }

  int maxDirLevel ( int Dir ) const {
    int maxLevel = -1;
    int thisdir = ( Dir & TN_DIR_ALL );
    int dira = 0, dirb = 0;


    switch ( thisdir ) {
    case TN_DIR_BACK:
      dira = TN_DIR_RIGHT; dirb = TN_DIR_UP;
      break;
    case TN_DIR_RIGHT:
      dira = TN_DIR_BACK; dirb = TN_DIR_UP;
      break;
    case TN_DIR_UP:
      dira = TN_DIR_BACK; dirb = TN_DIR_RIGHT;
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unhandled switch case!", __FILE__, __LINE__ );
    }
    for ( int dir = 0; dir < 4; dir++ ) {
      int dirc;
      if ( ! ( Dir & TN_DIR_OPP ) )
        dirc = thisdir;
      else
        dirc = 0;


      if ( dir & 1 )
        dirc |= dira;
      if ( dir & 2 )
        dirc |= dirb;

      if ( dirLevel[ dirc ] > maxLevel )
        maxLevel = dirLevel[ dirc ];
    }
    return maxLevel;
  }

  int minDirLevel ( int Dir ) const {
    int minLevel = -1;
    int thisdir = ( Dir & TN_DIR_ALL );
    int dira = 0, dirb = 0;

    switch ( thisdir ) {
    case TN_DIR_BACK:
      dira = TN_DIR_RIGHT; dirb = TN_DIR_UP;
      break;
    case TN_DIR_RIGHT:
      dira = TN_DIR_BACK; dirb = TN_DIR_UP;
      break;
    case TN_DIR_UP:
      dira = TN_DIR_BACK; dirb = TN_DIR_RIGHT;
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unhandled switch case!", __FILE__, __LINE__ );
    }
    for ( int dir = 0; dir < 4; dir++ ) {
      int dirc;
      if ( ! ( Dir & TN_DIR_OPP ) )
        dirc = thisdir;
      else
        dirc = 0;


      if ( dir & 1 )
        dirc |= dira;

      if ( dir & 2 )
        dirc |= dirb;

      if ( dirLevel[ dirc ] > minLevel )
        minLevel = dirLevel[ dirc ];
    }

    return minLevel;
  }

bool operator== ( const TrialNode<T> &Other ) const { return ( val == Other.val ); }

  bool operator< ( const TrialNode<T> &Other ) const { return ( val < Other.val ); }

  bool operator> ( const TrialNode<T> &Other ) const { return ( val > Other.val ); }


  int x() const { return _x; }
  int y() const { return _y; }
  int z() const { return _z; }

  int _x, _y, _z;
  T val;
  //  T curv;
  //int curvSet;
  //float dx, dy, dz;

  QC_SEED_TYPE seed_type;

  signed char faceHangNode;
  signed char edgeHangNode;

  short hangNodeLevel;
  short dirLevel[ 8 ];

  void dump( ) const {

    // cout << " - Node: (" << x << ", " << y << ", " << z << ") val = " << val << endl;
    for ( int dir = 0; dir <= TN_DIR_MAX_NUM; dir++ ) {
      cout << "   dirLevel[ ";
      parsedir ( dir );
      cout << " ] = " <<  static_cast<int> ( dirLevel[ dir ] ) << endl;
    }
    cout << "   faceHangNode = " <<  static_cast<int> ( faceHangNode )
    << "   edgeHangNode = " <<  static_cast<int> ( edgeHangNode ) << endl;

    cerr << "hangNodeLevel = " << hangNodeLevel << endl;
  }

  static void parsedir ( int dir ) {
    if ( dir & TN_DIR_RIGHT )
      cout << "RIGHT ";
    else
      cout << "LEFT  ";

    if ( dir & TN_DIR_UP )
      cout << "UP    ";
    else
      cout << "LOW   ";

    if ( dir & TN_DIR_BACK )
      cout << "BACK  ";
    else
      cout << "FRONT ";
  }
};

template <class T>
void TrialNode<T>::addHangNode ( int HNValue, int HNLevel ) {
  if ( hangNodeLevel != TN_HN_NONE_SEG && hangNodeLevel != HNLevel ) {
    // cerr << "WARNING: setting different hangnodelevel\n";
  }
  hangNodeLevel = HNLevel;

  if ( HNValue < 12 ) {
    edgeHangNode = HNValue;
  } else {
    if ( faceHangNode != TN_HN_NONE_SEG && faceHangNode != HNValue ) {
      cerr << "WARNING: setting different hangNode value!\n";
    }
    faceHangNode = HNValue;
  }
}

}   // end namespace eik

#endif
