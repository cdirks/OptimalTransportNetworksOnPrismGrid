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
 * \brief Calculates one heat equation step on the CPU and on the GPU on a 2D quadratic
 *        grid using structured sparse matrices.
 *
 * \author Berkels
 */

#include <generator.h>
#include <configurators.h>
#include "cudaMatrix.h"

/**
 * Rudimentary example code that shows how to use CSR matrixes on the GPU.
 *
 * \note Requires CUDA 3.2 or newer
 *
 * \author Berkels
 */
/* Since the final release of CUDA 3.2 was not ready (only a RC was available) when
 * the following code was written, it is commented out.
template <typename RealType>
void testGPUCSRMatrix ( const int numRows, const int numCols ) {
  // Generate a sparse matrix
  aol::SparseMatrix<RealType> mat ( numRows, numCols );
  mat.set(0, 0, 4.0);
  mat.set(0, 1, 1.0);
  for(int i = 1; i < (numRows - 1); i++){
    mat.set(i, i-1, 1.0);
    mat.set(i, i  , 40.0);
    mat.set(i, i+1, 1.0);
  }
  mat.set(numRows-1, numRows-2, 1.0);
  mat.set(numRows-1, numRows-1, 4.0);

  // Manually generate the same matrix in CSR format.
  aol::Vector<RealType> values ( 3 *numRows - 2 );
  aol::Vector<int> rowStartingIndices ( numRows + 1 );
  aol::Vector<int> columnIndices ( 3 * numRows - 2 );
  rowStartingIndices[0] = 0;
  values[0] = 4;
  values[1] = 1;
  columnIndices[0] = 0;
  columnIndices[1] = 1;
  int index=2;
  for(int i = 1; i < numRows-1; i++){
    rowStartingIndices[i] = index;
    values[index] = 1;
    values[index+1] = 40;
    values[index+2] = 1;
    columnIndices[index] = i-1;
    columnIndices[index+1] = i;
    columnIndices[index+2] = i+1;
    index += 3;
  }
  rowStartingIndices[numRows-1] = index;
  values[index] = 1;
  values[index+1] = 4;
  columnIndices[index] = numRows - 2;
  columnIndices[index+1] = numRows - 1;
  rowStartingIndices[numRows] = index+2;

  // Generate a right hand side.
  aol::Vector<RealType> vec ( numCols );
  for ( int j = 0; j < numCols; ++j )
    vec[j] = static_cast<RealType> ( j );

  // Vector to store the solution.
  aol::Vector<RealType> destVec ( numRows );

  // Invert the matrix.
  aol::DiagonalPreconditioner<aol::Vector<RealType> > prec ( mat );
  aol::PCGInverse<aol::Vector<RealType> > solver ( mat, prec, 1e-32 );
  aol::StopWatch watch;
  watch.start();
  solver.apply ( vec, destVec );
  watch.stop();
  watch.printReport( cerr );

  // Allocate memory on the CPU and copy all data to the GPU.
  cuda::CUDASparseHandle sparseHandle;
  cuda::GPUCSRMatrix<RealType> gpuCSRMat ( numRows, numCols, values.size(), values.getData(), rowStartingIndices.getData(), columnIndices.getData(), sparseHandle );
  cuda::GPUDiagonalMatrix<RealType> gpuPrecMat ( numRows );
  gpuPrecMat.copyFromHost ( prec.getPreconditionerMatrixReference().getDiagVectorReference().getData() );
  cuda::GPUVector<RealType> gpuDestVec ( numRows );
  gpuDestVec.setZero ( );
  cuda::GPUVector<RealType> gpuVec ( numCols );
  gpuVec.copyFromHost ( vec.getData() );

  // Invert the matrix on the GPU.
  aol::PCGInverse<cuda::GPUVector<RealType>, cuda::GPUCSRMatrix<RealType>, cuda::GPUDiagonalMatrix<RealType> > gpuSolver ( gpuCSRMat, gpuPrecMat, 1e-32 );
  watch.start();
  gpuSolver.apply ( gpuVec, gpuDestVec );
  watch.stop();
  watch.printReport( cerr );

  aol::Vector<RealType> destVecTest ( numRows );
  gpuDestVec.copyToHost ( destVecTest.getData() );
  destVecTest -= destVec;
  cerr << "Relative difference of CPU and CPU solution: " << ( destVecTest.normSqr() / destVec.normSqr() );
}
*/

typedef float RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType;

int main( int, char** ) {

  try {
    aol::StopWatch watch;

    cuda::initLibCublas();
    // The GPU solver performs best if the grid can be divided exactly in blocks, i.e. 256 = 16*16.
    const qc::RectangularGrid<qc::QC_2D> grid( qc::GridSize<qc::QC_2D> ( 256, 256 ) );
    ConfType::ArrayType uOld( grid );
    uOld.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
    ConfType::ArrayType uNew( grid );
    uNew.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
    ConfType::ArrayType rhs( grid );
    qc::DataGenerator<ConfType> generator ( grid );
    generator.generateRectangular ( uOld );
    uOld.savePNG ( "test-input.png" );

    // Assemble the system matrix.
    RType tau = 0.01;
    aol::StiffOp<ConfType> stiff( grid );
    aol::MassOp<ConfType> mass( grid );
    qc::FastUniformGridMatrix<RType, qc::QC_2D> mat( grid.getSize() );
    stiff.assembleAddMatrix( mat );
    mat *= ( tau );
    mass.assembleAddMatrix( mat );
    // Assemble the right hand side.
    mass.apply( uOld, rhs );

    // Diagonal-preconditioning doesn't help when inverting (M+0.01*L), but is used here
    // to show that the GPU solver supports it.
    aol::DiagonalPreconditioner<aol::Vector<RType> > prec ( mat );
    aol::PCGInverse<aol::Vector<RType> > inv ( mat, prec, 1e-16, 1000 );
    inv.setStopping( aol::STOPPING_ABSOLUTE );
    watch.start();
    inv.apply( rhs, uNew );
    watch.stop();
    cerr << "\nCGInverse CPU\n";
    watch.printReport( cerr );
    uNew.savePNG ( "test-solution-CPU.png" );

    const int numRows = grid.getNumberOfNodes();
    // Copy rhs to the GPU.
    cuda::GPUVector<RType> gpuRhs ( numRows );
    gpuRhs.copyFromHost ( rhs.getData() );
    // Initialize the solution vector in the GPU.
    cuda::GPUVector<RType> gpuUNew ( numRows );
    gpuUNew.setZero();

    // To copy the structured sparse system matrix to the GPU, we first copy the data
    // to a single memory block on the CPU and then copy the whole block to the GPU.
    aol::Vector<RType> rowVector ( 9 * numRows );
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < 9; j++ )
        rowVector[i*9+j] = mat.getInternalDataReference( j / 3 )[i][j%3];

    cuda::GPUDiagonalMatrix<RType> gpuPrecMat ( numRows );
    gpuPrecMat.copyFromHost ( prec.getPreconditionerMatrixReference().getDiagVectorReference().getData() );
    cuda::GPU2DQuadraticGridMatrix<RType> gpuMat ( numRows );
    gpuMat.copyFromHost ( rowVector.getData() );

    aol::PCGInverse<cuda::GPUVector<RType>, cuda::GPU2DQuadraticGridMatrix<RType>, cuda::GPUDiagonalMatrix<RType> > gpuSolver ( gpuMat, gpuPrecMat, 1e-16, 1000 );
    watch.start();
    gpuSolver.apply ( gpuRhs, gpuUNew );
    watch.stop();
    cerr << "\nCGInverse GPU\n";
    watch.printReport( cerr );

    gpuUNew.copyToHost ( uNew.getData() );
    uNew.savePNG ( "test-solution-GPU.png" );

    return 0;
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
