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

#ifndef __LEVELSETS_H
#define __LEVELSETS_H

// this file provides test levelsets for tpcfe

// TODO: think about where to move this code

// TODO: use brickIterator instead of three nested loops! maybe allow Vectors if nested in class containing grid.

#include <tpCFEStandardOp.h>
#include <tpCFELevelsets.h>

#include <scalarArray.h>
#include <progressBar.h>
#include <bitVector.h>
#include <randomGenerator.h>

#include "tpcfe_utils.h"

namespace tpcfe {

template< typename DataType >
void generate_slot_rhs ( qc::ScalarArray<DataType, qc::QC_3D> &rhs, const DataType width = 0.1, const DataType depth = 0.8 ) {

  const DataType de = depth * rhs.getNumX() - width * rhs.getNumZ() ;


  for ( int k = 0; k < rhs.getNumZ(); ++k ) {
    for ( int j = 0; j < rhs.getNumY(); ++j ) {
      for ( int i = 0; i < rhs.getNumX(); ++i ) {

        DataType value = 0.0;

        if ( i > de ) {
          // "base"
          value = 0.0;
        } else {
          if ( k < 0.5 ) {
            // "upper part of horseshoe"
            value = 1.0;
          } else {
            // "lower part of horseshoe"
            value = -1.0;
          }
        }

        rhs.set ( i, j, k,  value );
      }
    }
  }

}

} // end namespace tpcfe

#endif
