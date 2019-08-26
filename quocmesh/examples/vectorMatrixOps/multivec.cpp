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
 *  \brief usage of aol::MultiVectors and appropriate block operators
 *  \author Schwen
 */


#include <vec.h>
#include <matrix.h>
#include <sparseMatrices.h>

int main ( int, char**) {
  try {

    { // first, we play around with MultiVectors:

      aol::MultiVector<double> multivec ( 3, 5 ); // create MultiVector with 3 components of size 5.

      multivec[0].setZero();                      // this is how to access one of the MultiVector's components, a Vector
      multivec[1][2] = 2.0;                       // this is how to access an element

      aol::Vector<double> vec ( 16 );
      aol::Vector<double> anothervec ( 32 );

      multivec.appendReference ( vec );           // now multivec has one more component, a reference to vec (note that its size is different from the other components)
      multivec.appendReference ( anothervec );

      bool samelength = multivec.allDimsEqual();       // can you guess the result? :-)
      cout << "Same length? " << samelength << endl;

    }

    { // here's how BlockMatrix-MultiVector multiplication works:
      aol::MultiVector<double> multivec ( 3, 5 ), multiprod ( 2, 4 );
      aol::BlockMatrix< aol::SparseMatrix<double> > blockmat ( 2, 3, 4, 5 );
      blockmat.apply ( multivec, multiprod );
    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}

