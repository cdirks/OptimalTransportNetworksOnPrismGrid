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

#ifndef __SIGNEDDISTANCESWEEPING_H
#define __SIGNEDDISTANCESWEEPING_H

#include <op.h>
#include <quoc.h>
#include <scalarArray.h>
#include <bitArray.h>


namespace qc {

/** Operator to compute signed distance functions for zero level sets, using the sweeping method.
 *  This operator is not restricted to cubic domains.
 */
template< typename RealType >
class SignedDistanceSweepingOp3D : public aol::Op< qc::ScalarArray<RealType, qc::QC_3D> > {
protected:
  const RealType _delta;
public:
  SignedDistanceSweepingOp3D ( ) : _delta ( 1.0e-12f ) {
    // do nothing
  }

  explicit SignedDistanceSweepingOp3D ( const RealType Delta ) : _delta ( Delta ) {
    // do nothing
  }

  void applySingle ( qc::ScalarArray<RealType, qc::QC_3D> &Arg ) const {
    qc::ScalarArray<RealType, qc::QC_3D> tmp ( Arg, aol::DEEP_COPY );
    this->apply ( tmp, Arg );
  }

  //! compute distance information only on subset of relevant entries
  void applySingle ( qc::ScalarArray<RealType, qc::QC_3D> &Arg, const qc::BitArray<qc::QC_3D> &pIsRelevant ) const {
    qc::ScalarArray<RealType, qc::QC_3D> tmp ( Arg, aol::STRUCT_COPY );
    applyAdd ( Arg, tmp, &pIsRelevant );
    Arg = tmp;
  }

protected:
  //! applyAdd does not really make sense for signed distance functions
  void applyAdd ( const qc::ScalarArray<RealType, qc::QC_3D> &Arg, qc::ScalarArray<RealType, qc::QC_3D> &Dest ) const {
    applyAdd ( Arg, Dest, NULL );
  }

  void applyAdd ( const qc::ScalarArray<RealType, qc::QC_3D> &Arg, qc::ScalarArray<RealType, qc::QC_3D> &Dest, const qc::BitArray<qc::QC_3D>* const pIsRelevant ) const ;

};

}

#endif
