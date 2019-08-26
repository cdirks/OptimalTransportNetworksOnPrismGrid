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

/**
 * \file
 *
 * \brief GPU implementation of examples/MumfordShahSegmentation/DualizedBinaryMSSeg.cpp
 *
 * \author Berkels
 */

#include <vec.h>
#include <scalarArray.h>
#include <finiteDifferences.h>
#include "cudaVector.h"
#include "cudaMatrix.h"

extern void MSSegGPU(float* indicator1, float* indicator2, int nx, int ny, float lambda, float tol, float timeStep, int maxIter);

/**
 * \author Berkels
 */
template <typename ConfiguratorType, int ImageDimension>
class GPUPiecewiseConstantTwoPhaseMSSegmentor : public qc::PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, ImageDimension> {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
public:
  GPUPiecewiseConstantTwoPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                                            const RealType Gamma,
                                            aol::MultiVector<RealType> &ImageMVec )
    : qc::PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, ImageDimension> ( Initializer, Gamma, ImageMVec ) {}

private:
  void doSegment ( ArrayType &Segmentation, qc::MultiArray<RealType, ConfiguratorType::Dim> * /*PDual*/ ) const {
    ArrayType indicator2 ( this->_grid );
    ArrayType indicator1 ( this->_grid );
    this->generateIndicatorFunction ( 1, indicator2 );
    this->generateIndicatorFunction ( 0, indicator1 );
    MSSegGPU( indicator1.getData(), indicator2.getData(), indicator1.getNumX(), indicator1.getNumY(), this->_gamma, this->_stopEpsilon, this->_tau, this->getMaxIterations ( ) );
    Segmentation = indicator1;
  }

};

typedef float RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType;

int main( int, char** ) {

  try {
    qc::ScalarArray<float, qc::QC_2D> input ( "../../../input/2D/cameraman.pgm" );
    input.scaleValuesTo01();

    ConfType::InitType grid ( input.getSize() );
    aol::MultiVector<RType> inputMVec ( 0, 0 );
    inputMVec.appendReference ( input );

    GPUPiecewiseConstantTwoPhaseMSSegmentor<ConfType, 1>  segmentor ( grid, 0.005 ,inputMVec );
    segmentor.setStopEpsilon ( 0.000001 );
    segmentor.setMaxIterations ( 10000 );

    qc::ScalarArray<RType, ConfType::Dim> temp ( grid );

    aol::StopWatch watch;
    watch.start();

    segmentor.segmentAndAdjustGrayValues ( temp );
    temp.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
    temp.savePNG ( "u.png" );
    temp.threshold ( 0.5, 0, 1 );
    temp.savePNG ( "segmented.png" );

    watch.stop();
    watch.printReport( cerr );
  }
  catch(std::exception &ex){
    cerr << aol::color::error << endl;
    cerr << "\n\nstd::exception caught:\n";
    cerr << ex.what () << endl;
  }
  catch(aol::Exception &ex){
    cerr << aol::color::error << endl;
    cerr << "\n\naol::Exception caught:\n";
    ex.dump ();
  }
  catch (...){
    cerr << aol::color::error << endl;
    cerr << "\n\nUnknown exception caught.\n";
  }
  cerr << aol::color::reset;
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
