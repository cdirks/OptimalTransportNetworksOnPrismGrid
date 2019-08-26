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

#include <parameterParser.h>
#include <scalarArray.h>
#include <tpCFEGrid.h>
#include <quocMultigrid.h>
#include <bitVector.h>
#include <multiArray.h>
#include <parameterParser.h>

#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMultigrid.h>
#include <tpCFEUtils.h>

#include "computeForce.h"
#include "thresholding.h"


typedef double RealType;
typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD, tpcfe::IsotropicElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> > ConfiguratorType;

int main ( int, char ** ) {
  
  try {

#if 1
    {
      const std::string
        scan_filenamemask = "recon_1007_L4_PMMA/full_crop2/1007_L4_PMMA_rec%04d.pgm",
        roiFilename = "out/1007_L4_PMMA_90mu.dat.bz2";
      
      const int
        scan_size_x      = 1824,
        scan_size_y      = 1461,
        scan_first_slice = 200,
        scan_last_slice  = 898;
      
      qc::ScalarArray< unsigned char, qc::QC_3D > loadImage ( scan_size_x, scan_size_y, scan_last_slice - scan_first_slice + 1 );
      cerr << loadImage.size() << endl;
      loadImage.setQuietMode ( true );
      cerr << "loading raw data ...";
      loadImage.loadSlices ( scan_filenamemask.c_str(), qc::QC_Z, scan_first_slice, scan_last_slice );
      cerr << " done." << endl;
      
      const unsigned int resampleFactor = 3;
      
      qc::ScalarArray< RealType, qc::QC_3D > ROI ( 1824/resampleFactor, 1461/resampleFactor, 699/resampleFactor );
      
      for ( qc::RectangularIterator< qc::QC_3D > rit ( ROI ); rit.notAtEnd(); ++rit ) {
        RealType averageVal = 0;
        for ( qc::RectangularIterator<qc::QC_3D> bit ( qc::CoordType( 0, 0, 0 ), qc::CoordType( resampleFactor, resampleFactor, resampleFactor ) ); bit.notAtEnd(); ++bit ) {
          averageVal += loadImage.get ( static_cast<short>( resampleFactor ) * (*rit) + (*bit) );
        }
        ROI.set ( *rit, averageVal / aol::Cub ( resampleFactor ) );
      }
      
      ROI.saveSlices ( "out/ROI_%03d.pgm", qc::QC_Z, qc::PGM_UNSIGNED_CHAR_BINARY, NULL, aol::CLIP_THEN_SCALE, 0, 255 );
      
      cerr << "Minimal and maximal values are: " << ROI.getMinValue() << " " << ROI.getMaxValue() << endl;
      
      cerr << "Save dataset ...";
      ROI.save ( roiFilename.c_str(), qc::SaveTypeTrait<RealType>::BinarySaveType );
      cerr << " done." << endl;
    }
#endif

#if 1
    {
      const std::string
        roiFilename = "out/1007_L4_PMMA_90mu.dat.bz2",
        segFilename = "out/1007_L4_PMMA_90mu-seg.dat.bz2";
      
      qc::ScalarArray < RealType, qc::QC_3D > ROI ( roiFilename.c_str() ), levelset ( ROI, aol::STRUCT_COPY );
      
      const RealType
        resolution = 90.0e-6, 
        aleph = resolution * ROI.getNumX();
      
      RealType threshold = RiCaThreshold<RealType> ( ROI, 0, 255, aol::Cub ( aleph ) );
      
      thresholdAndDepoempelize< RealType > ( ROI, threshold, aol::Cub ( aleph ), levelset, true );
      
      qc::ScalarArray<RealType, qc::QC_3D> savePattern ( levelset, aol::STRUCT_COPY );
      for ( int i = 0; i < levelset.size(); ++i ) {
        savePattern[i] = ( levelset[i] < 0 ? 255. : 0. );
      }
      savePattern.saveSlices ( "out/interior_%03d.pgm", qc::QC_Z, qc::PGM_UNSIGNED_CHAR_BINARY, NULL, aol::CLIP_THEN_SCALE, 0, 255 );
      
      levelset.save ( segFilename.c_str(), qc::SaveTypeTrait<RealType>::BinarySaveType );    
    }
#endif

#if 1
    {
      const std::string
        segFilename = "out/1007_L4_PMMA_90mu-seg.dat.bz2";
      
      qc::ScalarArray < RealType, qc::QC_3D > levelset ( segFilename );
      
      const RealType
        resolution = 90.0e-6, 
        aleph = resolution * levelset.getNumX(),

        E           = 13.0e9,                          // in Pa = N / m^2
        nu          = 0.32,                            // in 1
        lambda      = tpcfe::computeLambda ( E, nu ),
        mu          = tpcfe::computeMu ( E, nu );

      
      GridType grid ( qc::GridSize<qc::QC_3D>::createFrom ( levelset ) );
      
      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes();
      
      qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
      
      qc::MultiArray<RealType, 3> dirichletBCs ( grid ), u ( grid ), rhs ( grid );
      
      DirichletMask.setAll ( false );
      
      setZShift  ( grid, dirichletBCs, DirichletMask, - 1.0e-4 ); // compression
      
      grid.setDirichletMask ( DirichletMask );
      grid.setDOFMaskFromDirichletAndDomainNodeMask();
      
      cerr << "Using 3 x " << grid.getNumDOFs() << " = " << 3 * grid.getNumDOFs() << " degrees of freedom, total number of nodes in cube: " << grid.getNumberOfNodes() << endl;
      cerr << "Total inner volume = " << aol::detailedFormat ( grid.getTotalInnerVolume() * aol::Cub ( aleph ) ) << " m^3." << endl;
      
      cerr << "Setting up elasticity operator (based on mixed derivative ops)" << endl;
      tpcfe::CFEElastOp< ConfiguratorType > elastop ( grid, lambda, mu );
      
      u -= dirichletBCs;
      
      elastop.restrictNonDomainEntries();
      
      grid.restrictToDomain ( u );
      
      elastop.apply ( u, rhs );
      
      elastop.restrictDirichletEntries();
      grid.restrictDirichletNodes ( rhs );
      
      elastop.getBlockMatrixRef().getReference(0,0).printStatistics();
      qc::MultiArray<RealType, 3> soln ( grid );
      aol::StopWatch timerSolve;
      timerSolve.start();
      
      aol::SSORPreconditioner < aol::MultiVector<RealType>, ConfiguratorType::MatrixType > prec ( elastop.getBlockMatrixRef() );
      aol::PCGInverse< aol::MultiVector<RealType> > pcgsolver ( elastop.getBlockMatrixRef(), prec, 1.0e-12, 50000 );
      pcgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_RIGHT_HAND_SIDE );
    
      cerr << "memusage = " << ( aol::memusage() >> 20 ) << " MiB" << endl;
    
      pcgsolver.apply( rhs, soln );
    
      timerSolve.stop();
      cerr << "PCG(non-block SSOR) solver took " << timerSolve.elapsedWallClockTime() << " s (wall clock time)" << endl;
    
      soln += dirichletBCs;
    
      soln.save ( "out/1007_L4_PMMA_90mu_def%d.dat.bz2", qc::PGM_DOUBLE_BINARY );

      {
        const unsigned int surfaceDirection = 2;
        const int surfacePosition = -1;
        const RealType geomTolerance = 1.0e-8f;
        aol::Vec3<RealType> surfaceOuterNormal; surfaceOuterNormal[2] = 1.0f;
      
        tpcfe::IsotropicElasticityCoefficient<RealType> ENuMinus ( E, nu ), ENuPlus ( aol::NumberTrait<RealType>::NaN, aol::NumberTrait<RealType>::NaN );
        qc::AArray< tpcfe::IsotropicElasticityCoefficient<RealType>, qc::QC_3D > isoCoeff ( grid );
        tpcfe::setCoeffForLevelset ( isoCoeff, levelset, ENuMinus, ENuPlus );
      
        aol::Vector<double> dummy;
        tpcfe::computeForce ( grid, soln, surfaceDirection, surfaceOuterNormal, geomTolerance, surfacePosition, aleph, isoCoeff, cerr, dummy );
      }

    }
#endif

  } catch ( aol::Exception e ) {
    e.dump();
  }

  return ( EXIT_SUCCESS );
}
