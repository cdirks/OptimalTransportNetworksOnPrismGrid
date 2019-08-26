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

#ifndef __JUMPCOEFFMG_H
#define __JUMPCOEFFMG_H

#ifdef USE_EXTERNAL_HYPRE

#include <tpCFEMultigrid.h>
#include <boomerAMG.h>

namespace tpcfe {

template< class CFEOpType, aol::OperatorType RePrOpType, mg::CoarseningMode CMode >
class CFEHybridMultigrid : public CFEMultigrid<CFEOpType, RePrOpType, CMode> {
  typedef typename CFEOpType::RealType      RealType;
  typedef typename CFEOpType::VectorType    VectorType;
  typedef typename CFEOpType::MatrixType    MatrixType;
  typedef typename CFEOpType::GridType      GridType;

public:
  CFEHybridMultigrid ( const GridType &Grid, MatrixType &fineMatrix,
                       int ExplicitLevel = 2,
                       int PreSmooth = 3, int PostSmooth = 3,
                       aol::GaussSeidelSweepingMode GSmode = aol::GAUSS_SEIDEL_FORWARD,
                       RealType Relax = 1.0,
                       mg::MGCycleMode CycleMode = mg::MULTIGRID_V_CYCLE,
                       RealType CycleEps = 1.0e-16,
                       int MaxCycle = 100,
                       aol::StoppingMode stop = aol::STOPPING_UNSET,
                       const RealType coarsenDropSmallEntries = 0 ) : CFEMultigrid<CFEOpType, RePrOpType, CMode> ( Grid, fineMatrix, ExplicitLevel, PreSmooth, PostSmooth, GSmode, Relax, CycleMode, CycleEps, MaxCycle, stop, coarsenDropSmallEntries ) {
  }

public:
  virtual void solveCoarseProblem ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    if ( Level != this->_explicitLevel )
      throw aol::Exception ("solveCoarseProblem must not be called on other than explicitLevel", __FILE__, __LINE__ );
    mg::BoomerAMGSolver< MatrixType, RealType > AMGSolver ( this->getOperatorRef ( this->_explicitLevel ), this->coarseSolverEps, this->coarseSolverSteps );
    AMGSolver.apply ( Rhs, X );
  }
};

} // end namespace

#endif

#endif
