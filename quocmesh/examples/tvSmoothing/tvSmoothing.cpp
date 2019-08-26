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

#include <anisoStiffOps.h>
#include <configurators.h>
#include <gradientDescent.h>

using namespace aol;
using namespace qc;

typedef double                                                      RType;
typedef GaussQuadrature<RType, QC_2D, 3>                            QuadType;
typedef QuocConfiguratorTraitMultiLin<RType, QC_2D, QuadType>       ConfType;
//typedef RegularSimplexConfigurator<RType, QC_2D, QuadType, GridDefinition> ConfType;
typedef ScalarArray<RType, QC_2D>                                   Image;
typedef Vector<RType>                                               VType;

// ------------------------------------------------------------------------------------------------

class TVNorm : public FENonlinIntegrationScalarInterface<ConfType, TVNorm> {

public:
  TVNorm ( const ConfType::InitType & grid,
           RType regularizationEpsilon )
    : FENonlinIntegrationScalarInterface<ConfType, TVNorm> ( grid )
    , _regularizationEpsilon ( regularizationEpsilon )
  {}

  RType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfType> &DiscFuncs,
                               const ConfType::ElementType &El,
                               int QuadPoint, const ConfType::DomVecType & /*RefCoord*/ ) const {

    Vec<ConfType::DomDim, RType> grad;
    DiscFuncs.evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
    return sqrt ( grad.normSqr() + _regularizationEpsilon );
  }

protected:
  RType _regularizationEpsilon;
};
// ------------------------------------------------------------------------------------------------

class TVEnergyFirstVariation : public FELinScalarWeightedStiffInterface<ConfType, TVEnergyFirstVariation> {

public:
  TVEnergyFirstVariation ( const ConfType::InitType & grid, const VType & phi,
                           RType regularizationEpsilon )
    : FELinScalarWeightedStiffInterface<ConfType, TVEnergyFirstVariation> ( grid )
    , _phi ( grid, phi )
    , _regularizationEpsilon ( regularizationEpsilon )
  {}

  RealType getCoeff ( const ConfType::ElementType &El, int QuadPoint, const DomVecType & /* RefCoord */ ) const {
    Vec<ConfType::DomDim, RealType> grad;
    _phi.evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
    return 1. / sqrt ( grad.normSqr() + _regularizationEpsilon );
  }

protected:
  DiscreteFunctionDefault<ConfType> _phi;
  RType _regularizationEpsilon;
};
// ------------------------------------------------------------------------------------------------

class TVSmoothingEnergy : public Op<VType, Scalar<RType> > {
public:
  TVSmoothingEnergy ( const ConfType::InitType & grid,
                      const Image & phi_0,
                      RType mu )
    : _grid ( grid )
    , _phi_0 ( phi_0 )
    , _massOp ( grid, ASSEMBLED )
    , _mu ( mu )
  {}

  virtual void applyAdd ( const VType & phi, Scalar<RType> & dest ) const {

    Scalar<RType> ret;

    // TV norm
    TVNorm tvNorm ( _grid, /* norm regularization epsilon = */ 1E-5 );
    tvNorm.applyAdd ( phi, ret );
    ret *= _mu;

    // fidelity term
    VType phi_minus_phi_0 ( phi ),
          M_diff ( phi, STRUCT_COPY );
    phi_minus_phi_0 -= _phi_0;
    _massOp.apply ( phi_minus_phi_0, M_diff );
    ret += 0.5 * (M_diff * phi_minus_phi_0);

    dest += ret;
  }

protected:
  const ConfType::InitType & _grid;
  const Image & _phi_0;
  MassOp<ConfType> _massOp;
  RType _mu;
};
// ------------------------------------------------------------------------------------------------

class TVSmoothingEnergyFirstVariation : public Op<VType, VType> {
public:
  TVSmoothingEnergyFirstVariation ( const ConfType::InitType & grid,
                                    const Image & phi_0,
                                    RType mu )
    : _grid ( grid )
    , _phi_0 ( phi_0 )
    , _massOp ( grid, ASSEMBLED )
    , _mu ( mu )
  {}

  virtual void applyAdd ( const VType & phi, VType & dest ) const {

    VType ret ( dest, STRUCT_COPY );

    // TVNorm
    TVEnergyFirstVariation tvStiffOp ( _grid, phi, /* norm regularization epsilon = */ 1E-5 );
    tvStiffOp.apply ( phi, ret );
    ret *= _mu;

    // fidelity term
    VType phi_minus_phi_0 ( phi );
    phi_minus_phi_0 -= _phi_0;
    _massOp.applyAdd ( phi_minus_phi_0, ret );

    dest += ret;
  }

protected:
  const ConfType::InitType & _grid;
  const Image & _phi_0;
  MassOp<ConfType> _massOp;
  RType _mu;
};
// ------------------------------------------------------------------------------------------------

class GradientDescentWithPictureSaving : public QuocGradientDescent<ConfType, VType> {

public:
  GradientDescentWithPictureSaving ( const ConfType::InitType & grid,
                                     const Op<VType, Scalar<RType> > & E,
                                     const Op<VType, VType> & DE,
                                     int maxIter,
                                     string filename )
    : QuocGradientDescent<ConfType, VType> ( grid, E, DE, maxIter )
    , _grid ( grid )
    , _filename ( filename ) {

    setWriteTimeSteps ( true );
    setWriteAllTimeSteps ( true );
  }

protected:
  virtual void writeTimeStep ( const VType & dest, const int numIter ) const {

    Image phi ( dest, _grid, DEEP_COPY );
    cerr << "iterate range: " << phi.getMinValue() << " .. " << phi.getMaxValue() << endl;
    phi *= 255. / phi.getMaxAbsValue();
    phi.save ( strprintf(_filename.c_str(), numIter).c_str(),  PGM_UNSIGNED_CHAR_ASCII );
  }

  const ConfType::InitType & _grid;
  string _filename;
};
// ------------------------------------------------------------------------------------------------

int main ( int argc, char ** argv ) {

  if ( argc != 4 ) {
    cerr << "Usage: " << argv[0] << " <file.pgm> <mu> <maxIterations>" << endl
         << "Try "    << argv[0] << " ../testdata/image_257_noisy.pgm 0.1 10" << endl;
    exit(-1);
  }

  try {
    Image phi_0 ( argv[1] ),
          phi ( phi_0 );
    RType mu = atof ( argv[2] );
    int maxIter = atoi ( argv[3] );

    phi_0 /= phi_0.getMaxAbsValue();

    ConfType::InitType grid ( phi_0.getSize() );
    TVSmoothingEnergy E ( grid, phi_0, mu );
    TVSmoothingEnergyFirstVariation DE ( grid, phi_0, mu );

    GradientDescentWithPictureSaving gradDescent ( grid, E, DE, maxIter, "step%03i.pgm" );
    gradDescent.apply ( phi_0, phi );
  }
  catch (Exception & exc) {
    exc.dump();
    return -1;
  }

  return 0;
}
// ------------------------------------------------------------------------------------------------
