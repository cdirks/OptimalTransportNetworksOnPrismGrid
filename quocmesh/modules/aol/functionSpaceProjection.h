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

#ifndef __FUNCTIONSPACEPROJECTION_H
#define __FUNCTIONSPACEPROJECTION_H

#include <aol.h>
#include <FEOpInterface.h>
#include <discreteFunction.h>
#include <op.h>
#include <solver.h>

namespace aol {

/**
  * class for projecting discrete functions onto different FE-spaces.
  * \author Droske
  */
template <typename ConfigDomainType, typename ConfigRangeType>
class FE_L2Projection : public Op<aol::Vector<typename ConfigDomainType::RealType>, aol::Vector<typename ConfigRangeType::RealType> > {
protected:
  const typename ConfigDomainType::InitType &_initializer;
public:
  FE_L2Projection ( const typename ConfigDomainType::InitType &initializer ) :
      _initializer ( initializer ) {}

  virtual ~FE_L2Projection() {}

  void applyAdd ( const aol::Vector<typename ConfigDomainType::RealType> &arg,
                  aol::Vector<typename ConfigRangeType::RealType> &dest ) const {
    typedef typename ConfigDomainType::QuadType QType;
    typedef typename ConfigRangeType::QuadType QTypeRange;

    typedef typename ConfigDomainType::RealType RealType;

    DiscreteFunctionDefault<ConfigDomainType> discFuncArg ( _initializer, arg );
    DiscreteFunctionDefault<ConfigRangeType> discFuncDest ( _initializer, dest );

    const ConfigRangeType &config = discFuncDest.getConfigurator();

    if ( typeid ( QType ) != typeid ( QTypeRange ) ) {
      throw aol::Exception ( "this version of l2-projection only works when the same quadrature rules are used for both function spaces.\n" );
    }

    const int numDofsRange = config.getNumGlobalDofs();

    aol::Vector<RealType> rhs ( numDofsRange );
    rhs.setZero();

    for ( typename ConfigDomainType::ElementIteratorType elit = config.begin(); elit !=  config.end(); ++elit ) {
      // make quadrature
      const int numLocalDofs = config.getNumLocalDofs ( *elit );
      for ( int q = 0; q < QType::numQuadPoints; q++ ) {
        RealType argAtQuadPoint = discFuncArg.evaluateAtQuadPoint ( *elit, q );
        const RealType weight = config.getBaseFunctionSet ( *elit ).getWeight ( q );

        RealType volume = ( ConfigDomainType::Dim == qc::QC_2D ) ? aol::Sqr ( config.H ( *elit ) ) : aol::Cub ( config.H ( *elit ) );

        for ( int i = 0; i < numLocalDofs; i++ ) {
          rhs[config.localToGlobal ( *elit, i ) ] +=  volume * weight * argAtQuadPoint * config.getBaseFunctionSet ( *elit ).evaluate ( i, q );
        }
      }
    }

    aol::MassOp<ConfigRangeType> massOp ( config.getInitializer(), ONTHEFLY );
    aol::CGInverse<aol::Vector<RealType> > cg ( massOp, 1e-16, 1000 );
    cg.applyAdd ( rhs, dest );
  }

};

} // end namespace aol

#endif
