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

#include <eikonal3d.h>

namespace eik {

template <class T>
int TrialNode<T>::consistencyCheck ( int MaxLevel, MetaEstimator &estimator, Eikonal3d &segmentor ) {
  int fail = 0;
  int dir;

  if ( faceHangNode != TN_HN_NONE_SEG ) {
    int lev = MaxLevel - hangNodeLevel;
    int coarseStep = 1 << ( lev );
    int HalfStep = coarseStep >> 1;
    int cX = x() - segmentor.tn_hn_offs[ faceHangNode ][ 0 ] * HalfStep;
    int cY = y() - segmentor.tn_hn_offs[ faceHangNode ][ 1 ] * HalfStep;
    int cZ = z() - segmentor.tn_hn_offs[ faceHangNode ][ 2 ] * HalfStep;
    if ( !estimator.checkElement ( cX, cY, cZ, hangNodeLevel ) ) {
      cerr << " - consistency check failed: false hangNode detected! \n";
      fail = 4;
    }
  }
  for ( dir = 0; dir <= TN_DIR_MAX_NUM; dir++ ) {
    if ( dirLevel[ dir ] != TN_HN_NONE_SEG ) {
      int step = 1 << ( MaxLevel - dirLevel[ dir ] );
      int coarseStep = step << 1;

      // check if Node in grid of stepsize
      if ( x() % step != 0
           || y() % step != 0
           || z() % step != 0 ) {
        fail = 1;
        break;
      }
      // check Estimator
      if ( dirLevel[ dir ] < MaxLevel
           && !estimator.checkElement ( dir, x(), y(), z(), dirLevel[ dir ] ) ) {
        cerr << "   dir = ";
        parsedir ( dir );
        cerr << " dirLevel = " << dirLevel[ dir ] << " estimator false for element!\n";
        fail = 2;
        break;
      }
      if ( dirLevel[ dir ] > 0
           && ( x() % coarseStep == 0 )
           && ( y() % coarseStep == 0 )
           && ( z() % coarseStep == 0 )
           && estimator.checkElement ( dir, x(), y(), z(), dirLevel[ dir ] - 1 ) ) {
        cerr << "   dir = ";
        parsedir ( dir );
        cerr << " dirLevel = " << dirLevel[ dir ] << " estimator true for parent element!\n";
        fail = 3;
        break;
      }
    }
  }
  if ( !fail ) {
    return !fail;
  }
  cerr << "\n - consistency check failed. " << fail << "\n";
  dump( );

  if ( fail == 1 ) {
    for ( dir = 0; dir <= TN_DIR_MAX_NUM; dir++ ) {
      if ( dirLevel[ dir ] != TN_HN_NONE_SEG ) {
        int step = 1 << ( MaxLevel - dirLevel[ dir ] );
        if ( x() % step != 0
             || y() % step != 0
             || z() % step != 0 ) {
          cerr << "   ";
          parsedir ( dir );
          cerr << "  level = " << dirLevel[ dir ] << " step = "
          << ( 1 << (MaxLevel - dirLevel[ dir ]) ) << endl;
        }
      }
    }
  }
  return !fail;
}


template <class T>
void TrialNode<T>::debughalt( ) {
  int i = 0;
  i++;
  cerr << i << endl;
}

template class TrialNode<float>;

}   // end namespace eik
