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

#include <configurators.h>
/* ******************************************************************************************
 * testAnisoDerivative.cpp
 * This small program tests the derivative of an anisotropie. The anisotropie itself has
 * to be called via a function gammaNorm(...), the derivative via gammaFirstDerivative(...).
 * Author of this small program: Oliver Nemitz,
 * aol::testFirstDerivativeVector: Benjamin Berkels, validator: Marc Droske
 * Last date: 09.08.2006
****************************************************************************************** */

#include <aol.h>
#include <quoc.h>
#include <op.h>
#include <anisotropies.h>
#include <gradientflow.h>
#include <gridBase.h>
#include <scalarArray.h>
#include <mcm.h>
#include <Willmore.h>

typedef double RealType;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D, 3> > ConfType;

// typedef qcIsotropy<RealType> AnisoType;
typedef qc::Pedestal2d<RealType> AnisoType;

using namespace aol::color;
using namespace aol;


// The scalar energy op, this operator just uses the qc::GammaIntegrationOp,
// which integrates over the anisotropie. The qc::GammaIntegrationOp can't be
// used, because it is derived from FENonlinIntegrationVectorGeneralInterface,
// which uses a multivector as DomType.
template <typename AnisoType, typename DomType, typename RangeType=DomType>
class AnisoEnergyOp : public Op<DomType,RangeType> {

protected:
  const AnisoType &_anisotropy;
  const qc::GridDefinition &_grid;

public:
  AnisoEnergyOp( const qc::GridDefinition &grid, const AnisoType &Anisotropy ) :
    _anisotropy( Anisotropy ), _grid ( grid ) { }

  void applyAdd( const DomType &Arg, RangeType &Dest ) const {
    qc::GammaIntegrationOp< ConfType, AnisoType > E( _grid, _anisotropy );

    MultiVector<RealType> tmp(0,0);
    tmp.appendReference( Arg );
    E.applyAdd( tmp, Dest );
  }
};


// -------------------------------------------------------------------------------
// ------------------------ the main program -------------------------------------
// -------------------------------------------------------------------------------
int main( int /*argc*/, char **/*argv*/)
{
  try
    {
      // the direction of the derivative
      qc::ScalarArray<double, qc::QC_2D> position( "../Data2D/Circle/Sphere2Podest_Rad01_Abs005_100.bz2" );
      int N = position.getNumX ();
      int d = qc::logBaseTwo (N);
      qc::GridDefinition grid( d, qc::QC_2D );

      // the anisotropy
      AnisoType aniso( 1.047, 0.05, 0.05 );

      // declare the energy-op and it's first variation (qc::AnisotropyIntegrationOp)
      AnisoEnergyOp<AnisoType, aol::Vector<RealType>, aol::Scalar<RealType> > E( grid, aniso );
      qc::AnisotropyIntegrationOp< ConfType, AnisoType > DE( grid, aniso );


      // --------------------------- now the test of the derivative -----------------------
      char baseFileName[1024];
      sprintf( baseFileName, "EnergiePlot" );
      cerr << "testing the derivative ...";
//       aol::testFirstDerivativeVector<ConfType>( grid, position, E, DE, baseFileName, 0.001 );


      cerr << blue << "done! Thanx for using testAnisoDerivative!\n" << reset;
    }
  catch ( aol::Exception &ex )
    {
      ex.dump();
    }
}
