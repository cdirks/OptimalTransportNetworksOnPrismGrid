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

#include <aol.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <configurators.h>
#include <jointHistogram.h>
#include <mutualInformation.h>
#include "utilities.h"

template <typename ConfiguratorType>
class MutualInformationAnalyzer{
public:
  typedef typename ConfiguratorType::RealType RealType;
  enum ANALYZE_TYPE {
    TRANSLATION,
    ROTATION,
    GAUSSIAN_NOISE,
    EQUALLY_DISTRIBUTED_NOISE
  };
protected:
  ANALYZE_TYPE _analyzeType;
  const RealType _beta;
  const int _intensityLevel;
  const int _gridLevel;
  const qc::GridDefinition _grid;
  char _outputDir[1024];
public:
  MutualInformationAnalyzer( ANALYZE_TYPE AnalyzeType,
                             const RealType Beta,
                             const int IntensityLevel,
                             const char* OutputDir,
                             const int GridLevel = 8 ) :
    _analyzeType( AnalyzeType ),
    _beta( Beta ),
    _intensityLevel( IntensityLevel ),
    _gridLevel( GridLevel ),
    _grid( _gridLevel, qc::QC_2D ){
    sprintf( _outputDir, "%s", OutputDir );
  }

  void analyze( ) {
    char fn[1024];
    int numberOfSamples = 1;
    switch( _analyzeType ) {
    case TRANSLATION:
      {
        numberOfSamples = static_cast<int>(pow(2., _gridLevel));
      }
      break;
    case ROTATION:
      {
        numberOfSamples = 8*8+4*8+1;
      }
      break;
    case GAUSSIAN_NOISE:
      {
        numberOfSamples = static_cast<int>(pow(2., _gridLevel));
      }
      break;
    case EQUALLY_DISTRIBUTED_NOISE:
      {
        numberOfSamples = static_cast<int>(pow(2., _gridLevel));
      }
      break;
    }

    char namemienergydat[1024];
    sprintf( namemienergydat, "%smienergy.dat", _outputDir );
    ofstream outMIEnergy( namemienergydat, ios::out );
    char namementropydat[1024];
    sprintf( namementropydat, "%sentropy.dat", _outputDir );
    ofstream outJointEntropy( namementropydat, ios::out );

    typename ConfiguratorType::ArrayType r0( _grid );
    typename ConfiguratorType::ArrayType t0( _grid );

    // Generating images
    fillWithBox<RealType>( r0, 128, 0.);
    sprintf( fn, "%sreference", _outputDir );
    qc::writeImage<RealType>( _grid, r0, fn);
    r0 /= r0.getMaxValue();

    for( int i = 0; i < numberOfSamples; i++ ){

      switch( _analyzeType ) {
      case TRANSLATION:
        {
          fillWithBox<RealType>( t0, i, 0);
        }
        break;
      case ROTATION:
        {
          fillWithBox<RealType>( t0, 128, i*aol::NumberTrait<RealType>::pi/(8.*8.));
        }
        break;
      case GAUSSIAN_NOISE:
        {
          fillWithBox<RealType>( t0, 128, 0);
          aol::NoiseOperator<RealType> noiseOP(aol::NoiseOperator<RealType>::GAUSSIAN, 0., i*0.05);
          noiseOP.applyAdd(t0, t0);
        }
        break;
      case EQUALLY_DISTRIBUTED_NOISE:
        {
          fillWithBox<RealType>( t0, 128, 0);
          aol::NoiseOperator<RealType> noiseOP(aol::NoiseOperator<RealType>::EQUALLY_DISTRIBUTED, 0., i*0.05);
          noiseOP.applyAdd(t0, t0);
        }
        break;
      }

      sprintf( fn, "%stemplate%03d", _outputDir, i);
      qc::writeImage<RealType>( _grid, t0, fn);
      t0 /= t0.getMaxValue();
      qc::MIRegistrationEnergyWithRegardToPhi<ConfiguratorType> E(_grid, r0, t0, _intensityLevel, _beta);
      aol::MultiVector<RealType> zeroTransformation( qc::QC_2D, _grid.getNumberOfNodes() );
      aol::Scalar<RealType> energy;
      E.apply( zeroTransformation, energy );
      cerr << "Mutual Information energy : " << energy[0] << endl;
      outMIEnergy << i << " " << energy[0] << endl;

      qc::JointHistogram<RealType> Joint(r0, t0, _intensityLevel, _beta);
      Joint.computeJointHistogram();
      Joint.computeMutualInformation();
      const RealType jointEntropy = Joint.computeEntropyOfJointHistogram();
      cerr << "Joint Entropy : " << jointEntropy << endl;
      outJointEntropy << i << " " << jointEntropy << endl;

      qc::ScalarArray<RealType, qc::QC_2D> jointHistogram(Joint.getJointHistogram());
      jointHistogram *=255.0/jointHistogram.getMaxValue();
      clipAndLogaritmicScaleTo01<RealType> (jointHistogram, 0., 1., 1000.);

      typename ConfiguratorType::ArrayType combined;
      joinThreeArrays2D<RealType>( r0, t0, jointHistogram, combined, 10 );
      combined.setQuietMode( true );
      combined.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
      sprintf( fn, "%scombined%03d.pgm", _outputDir, i);
      combined.save( fn, qc::PGM_UNSIGNED_CHAR_BINARY );

      jointHistogram.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
      sprintf( fn, "%sHistogram%03d.pgm", _outputDir, i);
      jointHistogram.save( fn, qc::PGM_UNSIGNED_CHAR_BINARY );

      qc::ScalarArray<RealType, qc::QC_2D> mutualInformation(Joint.getMutualInformation());

      mutualInformation *=1.0/mutualInformation.getMaxValue();
      clipAndLogaritmicScaleTo01<RealType> (mutualInformation, 0., 1., 1000.);

      mutualInformation.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
      sprintf( fn, "%sMI%03d.pgm", _outputDir, i);
      mutualInformation.save( fn, qc::PGM_UNSIGNED_CHAR_BINARY );
    }
    outMIEnergy.close();
    outJointEntropy.close();
    aol::Plotter<RealType> miEnergyPlotter ( namemienergydat );
    miEnergyPlotter.genPNG();
    aol::Plotter<RealType> jointEntropyPlotter ( namementropydat );
    jointEntropyPlotter.genPNG();
  }
};

typedef double RType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfTypeFeatureDomain;

int main( int argc, char **argv ) {

  try {
    if ( argc < 2 ) {
      cerr << "USAGE: " << argv[0] << " <Bin> <beta>\n";
      return EXIT_FAILURE;
    }

    // get the parameters
    int intensityLevel = atoi( argv[1] );

    ConfType::RealType Beta = atof( argv[2] );
/*
    MutualInformationAnalyzer<ConfType> miAnalyzer( MutualInformationAnalyzer<ConfType>::ROTATION, Beta, intensityLevel, "test/" );
    miAnalyzer.analyze();
*/
    char fn[1024];
    char outputDir[1024];
    sprintf( outputDir, "test/" );
    int cur_level = 8;
    qc::GridDefinition grid( cur_level, ConfType::Dim );

    ConfType::ArrayType r0( grid );
    ConfType::ArrayType t0( grid );

    // loading images
    r0.load( "../../../input/2D/rdata.pgm" );
    //r0.load( "input/Jasmin_r" );
    r0 /= r0.getMaxValue();

    t0.load( "../../../input/2D/tdata.pgm" );
    //t0.load( "input/Jasmin_t" );
    t0 /= t0.getMaxValue();

    qc::MIRegistrationConfigurator<ConfType> miRegisConfig ( intensityLevel, Beta );
    qc::MIRegistrationEnergyWithRegardToPhi<ConfType> E(grid, r0, t0, miRegisConfig);
    aol::MultiVector<RType> zeroTransformation( ConfType::Dim, grid.getNumberOfNodes() );
    aol::Scalar<RType> energy;
    E.apply( zeroTransformation, energy );
    cerr << "Mutual Information energy : " << energy[0] << endl;

    qc::JointHistogram<RType> Joint(r0,t0,intensityLevel,Beta);
    Joint.computeJointHistogram();
    Joint.computeMutualInformation();
    cerr << "Joint Entropy : " << Joint.computeEntropyOfJointHistogram() << endl;

    // get the 2D array of table
    qc::ScalarArray<RType, qc::QC_2D> jointHistogram(Joint.getJointHistogram());
    const RType jointHistogramSum = jointHistogram.sum();
    jointHistogram *=255.0/jointHistogram.getMaxValue();
    clipAndLogaritmicScaleTo01<RType> (jointHistogram, 0., 1.);

    if( cur_level!= intensityLevel )
      cerr << "Can not combine images, because cur_level!= intensityLevel\n";
    else{
      ConfType::ArrayType combined;
      joinThreeArrays2D<RType>( r0, t0, jointHistogram, combined, 10 );
      combined.setQuietMode( true );
      combined.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
      sprintf( fn, "%scombined.pgm", outputDir );
      combined.save( fn, qc::PGM_UNSIGNED_CHAR_BINARY );
    }

    // output the table
    jointHistogram.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
    sprintf( fn, "%sHistogram.pgm", outputDir );
    jointHistogram.save( fn, qc::PGM_UNSIGNED_CHAR_BINARY );

    // marginalization
    const int numberOfIntensityValues = static_cast<int>(pow((2.0),intensityLevel)+1);
    aol::Vector<RType> Pr(numberOfIntensityValues);
    Joint.marginalizeVector(Pr, 0);
    aol::Vector<RType> Pt(numberOfIntensityValues);
    Joint.marginalizeVector(Pt, 1);
    qc::ScalarArray<RType, qc::QC_2D> marTable_r (numberOfIntensityValues);
    Joint.marginalize(marTable_r, 0);
    qc::ScalarArray<RType, qc::QC_2D> marTable_t (numberOfIntensityValues);
    Joint.marginalize(marTable_t, 1);

    cerr << "sum of Pr = "<<Pr.sum()<<endl;
    cerr << "sum of Pt = "<<Pt.sum()<<endl;
    cerr << "sum of Prt= "<<jointHistogramSum<<endl;

    // output both maginal histograms
    sprintf( fn, "%sHisto_r.dat", outputDir );
    qc::writeHistogram(fn, marTable_r,0);
    sprintf( fn, "%sHisto_t.dat", outputDir );
    qc::writeHistogram(fn, marTable_t,1);

    qc::ScalarArray<RType, qc::QC_2D> histogram_r(marTable_r);
    histogram_r *=255.0/histogram_r.getMaxValue();
    histogram_r.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
    sprintf( fn, "%sHisto_r.pgm", outputDir);
    histogram_r.save( fn, qc::PGM_UNSIGNED_CHAR_BINARY );

    qc::ScalarArray<RType, qc::QC_2D> histogram_t(marTable_t);
    histogram_t *=255.0/histogram_t.getMaxValue();
    histogram_t.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
    sprintf( fn, "%sHisto_t.pgm", outputDir);
    histogram_t.save( fn, qc::PGM_UNSIGNED_CHAR_BINARY );

    // get the 2D array of table
    qc::ScalarArray<RType, qc::QC_2D> mutualInformation(Joint.getMutualInformation());

    // output the table
    mutualInformation *=1.0/mutualInformation.getMaxValue();
    clipAndLogaritmicScaleTo01<RType> (mutualInformation, 0., 1.);

    mutualInformation.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
    sprintf( fn, "%sMI.pgm", outputDir );
    mutualInformation.save( fn, qc::PGM_UNSIGNED_CHAR_BINARY );
  }//try
  catch ( aol::Exception &el ){
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
