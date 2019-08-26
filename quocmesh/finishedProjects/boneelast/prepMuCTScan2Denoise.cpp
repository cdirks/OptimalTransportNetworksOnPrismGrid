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

// 2nd program to load muCT datasets obtained by UW:
// denoise dataset

// #define SAVE_PGMS 1

#include <parameterParser.h>
#include <anisoStiffOps.h>
#include <solver.h>
#include <preconditioner.h>
#include <FEOpInterface.h>
#include <configurators.h>
#include <quocMatrices.h>
#include <tpCFEGrid.h>
#include <tpCFELevelsets.h>

using namespace std;

typedef double RealType;

int main ( int argc, char **argv ) {
  try{

    aol::ParameterParser params ( argc == 2 ? argv[1] : "par/prepareMuctScan.par" );

    char roiFilename[1024], denFilename[1024];
    params.getString ( "roiFilename", roiFilename );
    params.getString ( "denFilename", denFilename );

    qc::ScalarArray<RealType, qc::QC_3D> image ( roiFilename );
    image.setQuietMode ( true );

    tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > grid ( qc::GridSize<qc::QC_3D>::createFrom ( image ) );
    cerr << "grid.H() = " << grid.H() << endl;
    grid.setDomainFrom ( image );
    grid.detectAndInitVirtualNodes();
    cerr << "Volume before denoising: " << grid.getTotalInnerVolume() << endl;

#ifdef SAVE_PGMS
    image.saveSlices ( "out/Noisy_%03d.pgm", qc::QC_Z, qc::PGM_UNSIGNED_CHAR_ASCII, NULL, aol::CLIP_THEN_SCALE, 0, 255 );
#endif

    typedef qc::RectangularGridConfigurator< RealType, qc::QC_3D, aol::GaussQuadrature< RealType, qc::QC_3D, 3 > > ConfigType;

    if ( params.hasVariable ( "numAnisoSteps" ) && params.getInt ( "numAnisoSteps" ) > 0 ) {
      qc::AnisotropicDiffusion3DLevelSetStiffOp< ConfigType, qc::GeneralLinearSmoothOp< ConfigType > > L ( grid, aol::ONTHEFLY );
      aol::MassOp< ConfigType > M ( grid, aol::ONTHEFLY );

      RealType
        tau           = static_cast<RealType> ( params.getDouble ( "tau_factor" ) ) * grid.H(),
        lambda        = static_cast<RealType> ( params.getDouble ( "lambda" ) ),
        rho_factor    = static_cast<RealType> ( params.getDouble ( "rho_factor" ) ),
        sigma_factor  = static_cast<RealType> ( params.getDouble ( "sigma_factor" ) );

      L.setLambda ( lambda );
      L.setStructSmooth ( sigma_factor, rho_factor );

      qc::MultilinFEBandMatrix< RealType, qc::QC_3D > mat ( grid );

      cerr << "computing " <<  params.getInt ( "numAnisoSteps" ) << " steps of anisotropic diffusion" << endl;
      for ( int iter = 0 ; iter < params.getInt ( "numAnisoSteps" ) ; ++iter ) {
        qc::ScalarArray<RealType, qc::QC_3D> dwarf ( grid );

        mat.setZero();

        L.setImageReference ( image );
        L.computeStructureTensor ( image );
        cerr << "assembling L ";
        L.assembleAddMatrix ( mat );
        cerr << "scaling with tau ";
        mat *= ( tau );
        cerr << "assembling M ";
        M.assembleAddMatrix ( mat );

        cerr << "applying M ";
        M.apply ( image, dwarf );

        cerr << "setting up solver ";
        aol::DiagonalPreconditioner< aol::Vector< RealType> > precond ( mat );
        aol::PCGInverse< aol::Vector<RealType> > solver ( mat, precond, 1.0e-16, 10000 ); // read epsilon and max_num_iter from parameter file?
        solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
        cerr << endl;

        solver.apply ( dwarf, image );

#ifdef SAVE_PGMS
        if ( iter % 1 == 0 ) {
          char fnmask[1024];
          sprintf ( fnmask, "out/Denoised%02d%%03d.pgm", iter );
          cerr << fnmask << endl;
          image.saveSlices ( fnmask, qc::QC_Z, qc::PGM_UNSIGNED_CHAR_ASCII, NULL, aol::CLIP_THEN_SCALE, 0, 255 );
        }
#endif
      }
    }

    if ( params.hasVariable ( "numLinsmoothSteps" ) && params.getInt ( "numLinsmoothSteps" ) > 0 ) {

      qc::GeneralLinearSmoothOp< ConfigType > lso ( grid, params.getDouble ( "linsmoothSigmaFactor" ) * grid.H() );
      cerr << "computing " <<  params.getInt ( "numLinsmoothSteps" ) << " steps of isotropic diffusion" << endl;
      for ( int iter = 0 ; iter < params.getInt ( "numLinsmoothSteps" ) ; ++iter ) {
        lso.applySingle ( image );

#ifdef SAVE_PGMS
        if ( iter % 1 == 0 ) {
          char fnmask[1024];
          sprintf ( fnmask, "out/IsoDenoised%02d%%03d.pgm", iter );
          cerr << fnmask << endl;
          image.saveSlices ( fnmask, qc::QC_Z, qc::PGM_UNSIGNED_CHAR_ASCII, NULL, aol::CLIP_THEN_SCALE, 0, 255 );
        }
#endif
      }
    }

#ifdef SAVE_PGMS
    image.saveSlices ( "out/Denoised_%03d.pgm", qc::QC_Z, qc::PGM_UNSIGNED_CHAR_ASCII, NULL, aol::CLIP_THEN_SCALE, 0, 255 );
#endif

    {
      tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > grid2 ( qc::GridSize<qc::QC_3D>::createFrom ( image ) );
      grid2.setDomainFrom ( image );
      grid2.detectAndInitVirtualNodes();
      cerr << "Volume after denoising: " << grid2.getTotalInnerVolume() << endl;
    }

    if ( params.hasVariable ( "removePoempels" ) && params.getInt ( "removePoempels" ) != 0 ) {
      // tpcfe::AllFacePoempelRemover<RealType> poempelRemover;
      tpcfe::SmallComponentsPoempelRemover<RealType> poempelRemover;
      poempelRemover.applySingle ( image );

      tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > grid2 ( qc::GridSize<qc::QC_3D>::createFrom ( image ) );
      grid2.setDomainFrom ( image );
      grid2.detectAndInitVirtualNodes();
      cerr << "Volume after removing Poempels: " << grid2.getTotalInnerVolume() << endl;
    }


    cerr << "Save dataset ...";
    image.save ( denFilename, qc::PGM_FLOAT_BINARY );
    cerr << " done." << endl;

  } catch ( aol::Exception &e ) {
    e.dump();
  }

  return ( EXIT_SUCCESS );
}

