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

#ifndef __MPIMSSEG_H
#define __MPIMSSEG_H

#ifdef MPI__

#include <mpiIncludes.h>
#include "ghostExchange.h"
#include "mpiUtils.h"

/**
 * \author Springer, Berkels
 */
template <typename ConfiguratorType>
class FirstOrderPrimalTwoPhaseMSSegmentorEngineMPI : public qc::FirstOrderPrimalTwoPhaseMSSegmentorEngine<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

public:
  FirstOrderPrimalTwoPhaseMSSegmentorEngineMPI ( const typename ConfiguratorType::InitType &Initializer,
                                                 const RealType Gamma,
                                                 const ArrayType &Indicator1Plus2,
                                                 const ArrayType &Indicator2 )
    : qc::FirstOrderPrimalTwoPhaseMSSegmentorEngine<ConfiguratorType> ( Initializer, Gamma, Indicator1Plus2, Indicator2 ) {}

  void minimizeMPI ( ArrayType &U, const int ( &procneigh ) [3][2], qc::MultiArray<RealType, ConfiguratorType::Dim> *PDual = NULL ) const {
    qc::MultiArray<RealType, ConfiguratorType::Dim> p ( this->_grid );
    if ( PDual == NULL )
      this->initializeDualVariable ( p );
    else
      p = *PDual;
    qc::MultiArray<RealType, ConfiguratorType::Dim> divTemp ( this->_grid );
    ArrayType uOld ( U, aol::FLAT_COPY );
    ArrayType uNew ( this->_grid );
    int num_bins[3] = {uNew.getNumX(), uNew.getNumY(), uNew.getNumZ() };
    RealType tau = aol::ZOTrait<RealType>::one / 8;
    RealType sigma = tau;
    RealType theta;
    //    const RealType gammaHDependent = this->_gamma / this->_grid.H();
    const RealType gammaHDependent = this->_gamma;
    const RealType pockGamma =  0.7 / gammaHDependent;

    //    aol::ProgressBar<> progressBar ( "Minimizing" );
    //    progressBar.start ( this->_maxIterations );
    //    progressBar.display ( cerr );

    for ( int iterations = 0; iterations < this->_maxIterations; ++iterations ) {
      qc::ghost_exchange ( uOld.getData(), num_bins, procneigh );
      //      qc::ghost_exchange(uOld[1].getData(), num_bins, procneigh);
      //      qc::ghost_exchange(uOld[2].getData(), num_bins, procneigh);

      // Update p
      qc::calculateCentralFDGradientMPI_PBC<RealType, ConfiguratorType::Dim> ( uOld, divTemp );


      p.addMultiple ( divTemp, sigma );


      this->applyResolventOfAdjointRegTermSingle ( sigma, p );

      qc::ghost_exchange ( p[0].getData(), num_bins, procneigh );
      qc::ghost_exchange ( p[1].getData(), num_bins, procneigh );
      qc::ghost_exchange ( p[2].getData(), num_bins, procneigh );

      // Calc uNew
      qc::calculateCentralFDDivergenceMPI_PBC<RealType, ConfiguratorType::Dim> ( p, uNew );
      uNew.scaleAndAdd ( tau, uOld );

      this->applyResolventOfDataTermSingle ( tau / gammaHDependent, uNew );

      // Update parameters
      theta = 1 / sqrt ( 1 + 2 * pockGamma * tau );
      tau *= theta;
      sigma /= theta;

      uNew.scaleAndAddMultiple ( ( 1 + theta ), uOld, -theta );

      double max_change = 0.0;
      //iterate over inner points
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int k = 1; k < num_bins[2] - 1; ++k ) {
        for ( int j = 1; j < num_bins[1] - 1; ++j ) {
          for ( int i = 1; i < num_bins[0] - 1; ++i ) {
            double tmp = abs ( uOld.get ( i, j, k ) - uNew.get ( i, j, k ) ) ;
            if ( tmp > max_change )
#ifdef _OPENMP
#pragma omp critical
#endif
            {
              if ( tmp > max_change )
                max_change = tmp;
            }
          }
        }
      }

      double final_change;
      MPI_Allreduce ( ( void* ) &max_change, ( void* ) &final_change, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

      uOld = uNew;

      // If the change is small enough, we consider the algorithm to be converged
      if ( final_change < this->_stopEpsilon ) {
        break;
      }

    }

    if ( PDual != NULL )
      *PDual = p;
  }
};

/**
 * \author Springer, Berkels
 */
template <typename ConfiguratorType, int ImageDimension>
class MPIPiecewiseConstantTwoPhaseMSSegmentor : public qc::PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, ImageDimension> {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  const int ( &_procneigh ) [3][2];
public:
  MPIPiecewiseConstantTwoPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                                            const RealType Gamma,
                                            aol::MultiVector<RealType> &ImageMVec,
                                            const int ( &Procneigh ) [3][2] )
    : qc::PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, ImageDimension> ( Initializer, Gamma, ImageMVec ),
      _procneigh ( Procneigh ) {}

private:

  void doSegment ( ArrayType &Segmentation, qc::MultiArray<RealType, ConfiguratorType::Dim> *PDual ) const {
    ArrayType indicator2 ( this->_grid );
    ArrayType indicator1Plus2 ( this->_grid );

    this->generateIndicatorFunction ( 1, indicator2 );
    this->generateIndicatorFunction ( 0, indicator1Plus2 );
    indicator1Plus2 += indicator2;
    FirstOrderPrimalTwoPhaseMSSegmentorEngineMPI<ConfiguratorType> engine ( this->_grid, this->_gamma, indicator1Plus2, indicator2 );
    engine.setMaxIterations ( this->getMaxIterations ( ) );
    engine.setStopEpsilon ( this->getStopEpsilon( ) );
    if ( this->getStepSaverPointer ( ) )
      engine.setStepSaverReference ( * ( this->getStepSaverPointer ( ) ) );
    engine.minimizeMPI ( Segmentation, _procneigh, PDual );
  }
};
#endif

#endif // __MPIMSSEG_H
