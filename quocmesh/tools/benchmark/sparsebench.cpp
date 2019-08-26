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

/** \file
 *  \brief benchmark for sparse matrix classes
 *
 *  Performs matrix operations on different matrix types (for 2D grids
 *  and 3D grids; default depths are 9 and 6) and measures runtime.
 *  Runtime is printed in seconds and relative to aol::SparseMatrix.
 *  With parameters bench file <filename>, uses default parameters and
 *  appends relative result to benchmark table in given file Usage:
 *  sparsebench <depth2d> <depth3d> or sparsebench bench file
 *  <filename>
 *
 *  \author Schwen
 */


#include <quoc.h>
#include <scalarArray.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <preconditioner.h>
#include <sparseMatrices.h>
#include <bandMatrix.h>
#include <mcm.h>
#include <CSRMatrix.h>
#include <UGBMatrix.h>
#include <configurators.h>

const int times_assemble = 10;
const int times_apply    = 100;
const int times_delRC    = 100;

static double ref_inst, ref_assm, ref_appl, ref_reop, ref_drwc; // global variables are not nice, but useful here.

template <typename MatType, typename OpType>
double profileAssemble ( MatType &mat, const OpType &massOp ) {

  mat.setZero();

  aol::StopWatch watch;
  watch.start();
  for ( int i = 0; i < times_assemble; ++i ) {
    massOp.assembleAddMatrix ( mat );
  };
  watch.stop();

  return ( watch.elapsedCpuTime() / times_assemble );
}


template <typename MatType>
double profileApply ( const MatType &mat ) {

  aol::Vector<double> v ( mat.getNumRows() ); v.setZero();
  aol::Vector<double> w ( mat.getNumRows() ); w.setZero();

  for ( int i = 0; i < mat.getNumRows(); ++i ) {
    v.set ( i, sin ( static_cast<double>(i) ) );
  };

  aol::StopWatch watch;
  watch.start();
  for ( int i = 0; i < times_apply; ++i ) {
    mat.apply ( v, w );
  }
  watch.stop();

  return ( watch.elapsedCpuTime() / times_apply );
}


template <typename MatType>
double profileRowEntryOp ( const MatType & mat ) {
  aol::RowEntryOp< double, MatType > REOp ( mat, mat.getNumRows() );

  aol::Vector<double> v ( mat.getNumRows() ); v.setZero();
  aol::Vector<double> w ( mat.getNumRows() ); w.setZero();

  for ( int i = 0; i < mat.getNumRows(); ++i ) {
    v.set ( i, sin ( static_cast<double>(i) ) );
  };

  aol::StopWatch watch;
  watch.start();
  for ( int i = 0; i < times_apply; ++i ) {
    REOp.apply ( v, w );
  }
  watch.stop();

  return ( watch.elapsedCpuTime() / times_apply );

}


template <typename MatType>
double profileSetRowColToZero ( MatType &mat ) {

  aol::StopWatch watch;
  watch.start();
  for ( int i = 0; i < times_delRC / 2; ++i ) {
    mat.setRowColToZero (       mat.getNumRows()   / 3 );  // twice at somewhat central location
    mat.setRowColToZero ( ( 2 * mat.getNumRows() ) / 3 );  // to test rowwise and colwise access
  }
  watch.stop();

  return ( watch.elapsedCpuTime() / times_delRC );
}


template< qc::Dimension Dim, typename MatrixType, typename GridType, typename OpType >
void benchmarkMatrix ( const char* name, const GridType &grid, const OpType &op, bool isReference = false ) {
  static aol::MixedFormat sformat ( 6, 8 ), rformat ( 4, 2 );

  aol::StopWatch timer; timer.start();
  MatrixType mat ( qc::GridSize<Dim>::createFrom ( grid ) );
  timer.stop();
  const double inst = timer.elapsedCpuTime();
  const double assm = profileAssemble (   mat, op );
  const double appl = profileApply (      mat );
  const double reop = profileRowEntryOp ( mat );
  const double drwc = profileSetRowColToZero (  mat );

  if ( isReference ) {
    ref_inst = inst; ref_assm = assm; ref_appl = appl; ref_reop = reop; ref_drwc = drwc;
  }

  cerr << "Profiling " << name << " create    : " << sformat ( inst ) << " seconds: " << rformat ( 10 * ref_inst / inst ) << " rupsi" << endl
       << "Profiling " << name << " assemble  : " << sformat ( assm ) << " seconds: " << rformat ( 10 * ref_assm / assm ) << " rupsi" << endl
       << "Profiling " << name << " apply     : " << sformat ( appl ) << " seconds: " << rformat ( 10 * ref_appl / appl ) << " rupsi" << endl
       << "Profiling " << name << " RowEntryOp: " << sformat ( reop ) << " seconds: " << rformat ( 10 * ref_reop / reop ) << " rupsi" << endl
       << "Profiling " << name << " delR/C    : " << sformat ( drwc ) << " seconds: " << rformat ( 10 * ref_drwc / drwc ) << " rupsi" << endl
       << endl;
}

template< qc::Dimension Dim, typename MatrixType, typename GridType, typename OpType >
void benchmarkMatrixR ( const char* name, const GridType &grid, const OpType &op ) {
  static aol::MixedFormat sformat ( 6, 8 ), rformat ( 4, 2 );

  aol::StopWatch timer; timer.start();
  MatrixType mat ( qc::GridSize<Dim>::createFrom ( grid ) );
  timer.stop();
  const double inst = timer.elapsedCpuTime();
  const double assm = profileAssemble (   mat, op );
  const double appl = profileApply (      mat );
  const double reop = profileRowEntryOp ( mat );


  cerr << "Profiling " << name << " create    : " << sformat ( inst ) << " seconds: " << rformat ( 10 * ref_inst / inst ) << " rupsi" << endl
       << "Profiling " << name << " assemble  : " << sformat ( assm ) << " seconds: " << rformat ( 10 * ref_assm / assm ) << " rupsi" << endl
       << "Profiling " << name << " apply     : " << sformat ( appl ) << " seconds: " << rformat ( 10 * ref_appl / appl ) << " rupsi" << endl
       << "Profiling " << name << " RowEntryOp: " << sformat ( reop ) << " seconds: " << rformat ( 10 * ref_reop / reop ) << " rupsi" << endl
       << "Profiling " << name << " delR/C    : not available " << endl
       << endl;
}


int main ( int argc, char **argv ) {
  const int MY_CACHE_SIZE = 1024 * 1024; // 1024 KB, to be divided by number of bands (rounded up to nearest power of 2) times sizeof(double), for use with UniGridCSR_Matrix

  bool bench = false;
  string resultfilename = "/dev/null";

  if ( aol::checkForBenchmarkArguments ( argc, argv, resultfilename ) )
    bench = true;

  if ( argc == 1 ) bench = true;

  int depth_2d, depth_3d;
  if ( argc == 3 ) {
    depth_2d = atoi ( argv[1] );
    depth_3d = atoi ( argv[2] );
  } else {
    depth_2d = 9;
    depth_3d = 6;
  };

  try {
    cerr << "Profiling different types of matrices:" << endl
         << "- creating one instance of each matrix" << endl
         << "- assembling massOp into matrix " << times_assemble << " times" << endl
         << "- calling apply (matrix-vector multiplication) " << times_apply << " times" << endl
         << "- applying via RowEntryOp (tests makeRowEntries) " << times_apply << " times" << endl
         << "- deleting row and col " << times_delRC << " times" << endl << endl ;

    aol::StopWatch total_time;
    total_time.start();

    {

      cerr << "Performing tests in 2D, grid level " << depth_2d << endl
           << "=====================================" << endl << endl;

      qc::GridDefinition grid2d ( depth_2d, qc::QC_2D );

      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature< double, qc::QC_2D, 3 > > ConfigType2d;

      aol::MassOp< ConfigType2d >  massOp_2d ( grid2d, aol::ONTHEFLY  );


      benchmarkMatrix < qc::QC_2D, aol::SparseMatrix<double>,                                                     qc::GridDefinition, aol::MassOp<ConfigType2d> > ( "SparseMatrix 2D                  ", grid2d, massOp_2d, true );
      benchmarkMatrix < qc::QC_2D, qc::UniformGridSparseMatrix<double>,                                           qc::GridDefinition, aol::MassOp<ConfigType2d> > ( "UniformGridSparseMatrix 2D       ", grid2d, massOp_2d );
      benchmarkMatrix < qc::QC_2D, qc::MultilinFEBandMatrix<double, qc::QC_2D>,                                   qc::GridDefinition, aol::MassOp<ConfigType2d> > ( "MultilinFEBandMatrix             ", grid2d, massOp_2d );
      benchmarkMatrix < qc::QC_2D, qc::FastUniformGridMatrix<double, qc::QC_2D>,                                  qc::GridDefinition, aol::MassOp<ConfigType2d> > ( "FastUniformGridMatrix 2D         ", grid2d, massOp_2d );
      benchmarkMatrixR< qc::QC_2D, qc::FastAssembleUniformGridMatrix<double>,                                     qc::GridDefinition, aol::MassOp<ConfigType2d> > ( "FastAssembleUniformGridMatrix 2D ", grid2d, massOp_2d );
      benchmarkMatrix < qc::QC_2D, qc::UGBMatrix < 9, 4, MY_CACHE_SIZE / ( 16*8 ), double, qc::GridDefinition > , qc::GridDefinition, aol::MassOp<ConfigType2d> > ( "UGBMatrix 2D                     ", grid2d, massOp_2d );
      benchmarkMatrixR< qc::QC_2D, qc::UniGridCSR_Matrix<qc::QC_2D>,                                              qc::GridDefinition, aol::MassOp<ConfigType2d> > ( "UniGridCSRMatrix 2D              ", grid2d, massOp_2d );

    }

    {

      cerr << "Performing tests in 3D, grid level " << depth_3d << endl
           << "=====================================" << endl << endl;

      qc::GridDefinition grid3d ( depth_3d, qc::QC_3D );

      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > ConfigType3d;

      aol::MassOp< ConfigType3d >   massOp_3d ( grid3d, aol::ONTHEFLY );

      benchmarkMatrix < qc::QC_3D, aol::SparseMatrix<double>,                                                      qc::GridDefinition, aol::MassOp<ConfigType3d> > ( "SparseMatrix 3D            ", grid3d, massOp_3d, true );
      benchmarkMatrix < qc::QC_3D, qc::UniformGridSparseMatrix<double>,                                            qc::GridDefinition, aol::MassOp<ConfigType3d> > ( "UniformGridSparseMatrix 3D ", grid3d, massOp_3d );
      benchmarkMatrix < qc::QC_3D, qc::MultilinFEBandMatrix<double, qc::QC_3D>,                                    qc::GridDefinition, aol::MassOp<ConfigType3d> > ( "MultilinFEBandMatrix       ", grid3d, massOp_3d );
      benchmarkMatrix < qc::QC_3D, qc::FastUniformGridMatrix<double, qc::QC_3D>,                                   qc::GridDefinition, aol::MassOp<ConfigType3d> > ( "FastUniformGridMatrix 3D   ", grid3d, massOp_3d );
      benchmarkMatrix < qc::QC_3D, qc::UGBMatrix < 27, 8, MY_CACHE_SIZE / ( 32*8 ), double, qc::GridDefinition > , qc::GridDefinition, aol::MassOp<ConfigType3d> > ( "UGBMatrix 3D               ", grid3d, massOp_3d );
      benchmarkMatrixR< qc::QC_3D, qc::UniGridCSR_Matrix<qc::QC_3D>,                                               qc::GridDefinition, aol::MassOp<ConfigType3d> > ( "UniGridCSRMatrix 3D        ", grid3d, massOp_3d );

    }

    total_time.stop();

    total_time.printReport ( cerr );

    if ( bench ) {

      double nupsi = 100 * 220.257 / total_time.elapsedCpuTime(), wupsi = 100 * 220.257 / total_time.elapsedWallClockTime();

      // log benchmark results
      aol::logBenchmarkResult ( "sparsebench", nupsi, wupsi, resultfilename );
      // Number of Uniform grid heat conduction timesteps Per Second In appropriate units: 100 nupsi = brahma, single job.
      // in terms of cputime and wall clock time
    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}
