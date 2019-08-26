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
 * \brief Segments an image into two piecewise constant regions with the Mumford Shah model.
 *
 * Segments an image into two piecewise constant regions with the Mumford Shah model. A global MS
 * minimizer is obtained by using a quadratic Esedoglu model and a primal / dual algorithm propsed
 * by Antonin Chambolle and Thomas Pock in
 * {A first-order primal-dual algorithm for convex problems with applications to imaging}.
 *
 * \author Berkels
 */

#include <aol.h>
#include <parameterParser.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <ChanVese.h>
#include <generator.h>
#include <chanVeseDescent.h>
#include <finiteDifferences.h>
#include <firstOrderTVAlgos.h>

template <typename RealType>
RealType naiveNearestNeighborSearchDistance ( const std::set<aol::Vec<2, short> > &PointSet, const aol::Vec2<short> &Point ) {
  RealType distance = aol::NumberTrait<RealType>::Inf;
  for ( std::set<aol::Vec<2, short> >::const_iterator it = PointSet.begin(); it != PointSet.end(); ++it ) {
    // short is too small to store the distance.
    const aol::Vec2<RealType> distanceVec ( (*it)[0] - Point[0], (*it)[1] - Point[1] );
    distance = aol::Min ( distance, distanceVec.normSqr() );
  }
  return sqrt ( distance );
}

void cleanMask ( qc::BitArray<qc::QC_2D> &Mask, const int DropComponentsSmallerThan, const bool DropBoundaryComponents ) {
  qc::GridStructure grid ( qc::GridSize2d::createFrom ( Mask ) );
  qc::ScalarArray<int, qc::QC_2D> labelArray ( grid );
  const int numberOfLabels = qc::ConnectedComponentsLabeler::doLabel ( Mask, labelArray );

  std::set<int> boundaryComponents;
  if ( DropBoundaryComponents ) {
    for ( qc::GridStructure::AllBoundaryNodeIterator it = grid; it.notAtEnd(); ++it ) {
      if ( Mask.get ( *it ) )
        boundaryComponents.insert ( labelArray.get ( *it ) );
    }
  }

  for ( int i = 0; i < numberOfLabels; ++i ) {
    const int componentIndex = i + 1;
    const int componentSize = labelArray.numOccurence ( componentIndex );
    if ( ( componentSize < DropComponentsSmallerThan ) || ( boundaryComponents.find ( componentIndex ) != boundaryComponents.end() ) ) {
      for ( int j = 0; j < Mask.size(); ++j ) {
        if ( labelArray[j] == ( i + 1 ) )
          Mask.set ( j, false );
      }
    }
  }
}

void convertMaskToSet ( const qc::BitArray<qc::QC_2D> &Mask, std::set<aol::Vec<2, short> > &MaskSet ) {
  MaskSet.clear();
  for ( qc::RectangularIterator<qc::QC_2D, aol::Vec<2, short> > it ( Mask ); it.notAtEnd(); ++it ) {
    if ( Mask.get ( *it ) )
      MaskSet.insert ( *it );
  }
}

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {

  try {
    aol::ParameterParser parser( argc, argv, "MSSeg.par" );

    aol::StopWatch watch;
    watch.start();

    typedef RType RealType;
    typedef ConfType ConfiguratorType;

    const RealType gamma = parser.getDouble( "gamma" );

    qc::DefaultArraySaver<RealType, ConfiguratorType::Dim> saver;
    saver.initFromParser( parser, true );

    const int imageDimension = 1;

    qc::MultiArray<RealType, ConfiguratorType::Dim, imageDimension> u0 ( parser.getString( "input-image" ).c_str() );
    if ( parser.hasVariable( "thresholdInputAt" ) )
      u0.threshold ( parser.getDouble ( "thresholdInputAt" ), 0, 1 );
    else
      u0 /= u0.getMaxValue();
    u0.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
    u0.savePNG( saver.createSaveName( "", ".png", -1, "input" ).c_str() );
    ConfiguratorType::InitType grid ( qc::GridSize<ConfType::Dim> ( u0.getSize() ) );

    qc::PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, imageDimension, qc::FirstOrderPrimalTwoPhaseMSSegmentor<ConfiguratorType> > segmentor ( grid, gamma, u0 );
    qc::ScalarArray<RealType, ConfiguratorType::Dim> temp ( grid );

    segmentor.setMaxIterations ( parser.getInt( "numSteps") );
    segmentor.setStopEpsilon ( parser.getDouble ( "epsilon" ) );
    if ( parser.hasVariable( "outerIterations"  ) )
      segmentor.setOuterIterations( parser.getInt ( "outerIterations" ) );

    if ( parser.hasVariable( "initialGrayValues" ) ) {
      aol::MultiVector<RealType> initialGrayVaues;
      parser.getRealMultiVec( "initialGrayValues", initialGrayVaues );
      segmentor.getMeanValuesReference() = initialGrayVaues;
    }

    qc::MultiArray<RealType, ConfiguratorType::Dim> pDual ( grid );
    segmentor.segmentAndAdjustGrayValues ( temp, &pDual );

    saver.saveStep ( temp, -1, "u" );
    pDual.save ( ( saver.createSaveName ( "", "_%d", -1, "p" ) + qc::getDefaultArraySuffix ( ConfiguratorType::Dim ) ).c_str(), qc::PGM_DOUBLE_BINARY );

    if ( parser.hasVariable( "dropComponentsSizeThreshold" ) ) {
      const int dropComponentsSmallerThan = grid.getNumberOfNodes() * parser.getDouble( "dropComponentsSizeThreshold" );
      qc::BitArray<qc::QC_2D> mask ( grid );
      mask.thresholdFrom ( temp, 0.5 );
      mask.invert();

      qc::BitArray<qc::QC_2D> maskMajor ( mask );
      cleanMask ( maskMajor, grid.getNumberOfNodes() * parser.getDouble( "dropDistantComponentsMajorSizeThreshold" ), true );

      cleanMask ( mask, dropComponentsSmallerThan, parser.checkAndGetBool ( "dropBoundaryComponents" ) );

      if ( parser.hasVariable( "dropDistantComponentsDistance" ) ) {
        std::set<aol::Vec<2, short> > maskMajorSet;
        convertMaskToSet ( maskMajor, maskMajorSet );

        const RType dropDistantComponentsDistance = parser.getDouble( "dropDistantComponentsDistance" ) / grid.H();
        qc::ScalarArray<int, qc::QC_2D> labelArray ( grid );
        const int numberOfLabels = qc::ConnectedComponentsLabeler::doLabel ( mask, labelArray );
        aol::Vector<RType> labelDistances ( numberOfLabels );
        labelDistances.setAll ( aol::NumberTrait<RType>::Inf );

        for ( qc::RectangularIterator<qc::QC_2D, aol::Vec<2, short> > it ( labelArray ); it.notAtEnd(); ++it ) {
          const int labelNum = labelArray.get( *it ) - 1;
          if ( ( labelNum >= 0 ) && ( labelDistances [ labelNum ] > 0 ) )
            labelDistances [ labelNum ] = aol::Min ( labelDistances [ labelNum ], naiveNearestNeighborSearchDistance<RType> ( maskMajorSet, *it ) );
        }

        for ( int i = 0; i < numberOfLabels; ++i ) {
          const int componentIndex = i + 1;
          if ( labelDistances [ componentIndex - 1 ] > dropDistantComponentsDistance ) {
            for ( int j = 0; j < mask.size(); ++j ) {
              if ( labelArray[j] == ( componentIndex ) )
                mask.set ( j, false );
            }
          }
        }
      }

      mask.invert();
      mask.save ( saver.createSaveName ( "", ".pgm", -1, "maskCleaned" ).c_str() );
      maskMajor.invert();
      maskMajor.save ( saver.createSaveName ( "", ".pgm", -1, "maskMajor" ).c_str() );
    }

    if ( ConfiguratorType::Dim == qc::QC_2D ) {
      qc::QuocTimestepSaver<ConfiguratorType> timestepSaver( parser.getInt( "saveOffset"), parser.getInt( "numSaveFirst"), "segmentation", true );
      timestepSaver.setSaveDirectory ( saver.getSaveDirectory() );
      aol::MultiVector<RealType> levelsetFunctions( 0, 0 );
      levelsetFunctions.appendReference( temp );
      timestepSaver.savePiecewiseConstantImageSegmented ( segmentor.getMeanValuesReference(), levelsetFunctions, grid, 0.5 );
    }

    watch.stop();
    watch.printReport( cerr );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
