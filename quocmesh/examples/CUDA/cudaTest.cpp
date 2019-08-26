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
 * \brief GPU implementation of a Jacobi Solver for full matrices.
 *
 * \author Berkels, Boerdgen
 */

#include <vec.h>
#include <matrix.h>
#include <Newton.h>
#include <scalarArray.h>

typedef float RType;

extern void cudaJacobiSolver ( float *pMatrix, float *pRhs, float *pSolution, const int Length );

int main( int, char** ) {

  try {
    const int size = 1024;

    const RType diag = 100;

    aol::Vector<RType> rhs ( size );
    rhs.setAll ( 1 );
    aol::Vector<RType> solution ( size );

    aol::FullMatrix<RType> fullMatrix ( size, size );
    fullMatrix.set (0, 0, diag );
    fullMatrix.set (0, 1, -1 );
    fullMatrix.set (size -1, size -1, diag );
    fullMatrix.set (size -1, size -2, -1 );
    for ( int i = 1; i < size -1 ; ++i ) {
      fullMatrix.set (i, i, diag );
      fullMatrix.set (i, i - 1, -1 );
      fullMatrix.set (i, i + 1, -1 );
    }

    aol::JacobiInverse< aol::Vector<RType>, aol::FullMatrix<RType> > jacobiSolver ( fullMatrix, 1e-12, 1000 );
    jacobiSolver.apply ( rhs, solution );

    aol::Vector<RType> cudaSolution ( size );

    cudaJacobiSolver ( fullMatrix.getDataVectorReference().getData(), rhs.getData(), cudaSolution.getData(), size );

    cudaSolution -= solution;
    cerr << "solution difference norm aol::JacobiInverse / cudaJacobiSolver: " << cudaSolution.norm() << endl << endl;
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
