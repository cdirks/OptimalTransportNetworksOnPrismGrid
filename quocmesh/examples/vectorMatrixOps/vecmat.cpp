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
 *  \brief usage of aol::Matrices and aol::Vectors
 *  \author Schwen
 */

#include <vec.h>
#include <matrix.h>
#include <sparseMatrices.h>

#include <gridBase.h>
#include <quocMatrices.h>

int main ( int, char**) {
  try {

    int size = 16;

    aol::Vector<double> vec ( size );   // create Vector of given size

    vec[2] = 42.0;                      // set and access its entries
    vec[5] = vec[2];

    // here are some methods a Vector provides:
    cout << "Norm = " << vec.norm() << ", maximal value = " << vec.getMaxValue() << endl;

    vec.resize ( 20 );                  // change size of vector, preserving old contents

    vec.setZero();                      // set all entries to zero

    vec.reallocate ( size );                 // change size, deleting old contents (allocating new memory)


    aol::FullMatrix<double> full_Matrix ( size, size );  // FullMatrix is a matrix class for dense matrices, all entries are stored.

    full_Matrix.setZero();                               // set all entries to zero
    full_Matrix.setIdentity();                           // obvious
    full_Matrix *= 2;                                    // multiply by two, standard operator *= can be used here, too.
    full_Matrix.set ( 2, 5, 42.0 );                      // accessing matrix entries
    full_Matrix.set ( 2, 6, full_Matrix.get ( 2, 5 ) );
    full_Matrix.add ( 2, 6, 23.0 );

    aol::Vector<double> prod ( size );
    full_Matrix.apply ( vec, prod );                     // matrix-vector-multiplication: prod = M * vec


    aol::SparseMatrix<double> sparse_Matrix ( size, size ); // this is a class for general unstructured sparse matrices, only nonzero entries are stored.
    sparse_Matrix.setIdentity();
    // filling it with the same values as above works exactly the same way.

    aol::DiagonalMatrix<double> diagonal_Matrix ( size );        // ... and many more Matrix classes are availabe. We're not giving a complete list here.
    diagonal_Matrix.setIdentity();


    qc::GridDefinition grid ( 4, qc::QC_2D );
    qc::UniformGridSparseMatrix<double> ugs_Matrix ( grid );  // This is one specially structured matrix that can create itself for a given grid.
    ugs_Matrix.set ( 3, 1, 1.0 );

    qc::ScalarArray<double, qc::QC_2D> array ( grid ), arrayprod ( grid );  // A ScalarArray<QC_2D> is a subclass of aol::Vector and can be used in the same way.
    array.set ( 2, 3, 20.0 );
    ugs_Matrix.apply ( array, arrayprod );

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}

