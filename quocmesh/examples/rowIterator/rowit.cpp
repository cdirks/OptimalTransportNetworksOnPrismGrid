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
 *  \brief undocumented
 *
 *  Not very useful and might be moved to vectormatrixops example.
 *
 *  \author Droske, Lenz, von Deylen
 */

#include <quoc.h>
#include <sparseMatrices.h>


int main() {

  aol::SparseMatrix<double> m( 100, 100 );

  m.set( 12, 15, 2. );

  m.set( 12, 18, 5. );

  m.set( 12, 12, 1. );



  vector<aol::Row<double>::RowEntry> vec;
  m.getRow(12).makeRowEntries( vec, 12 );

  for ( vector<aol::Row<double>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
    cerr << it->col << " " << it->value << endl;
  }

  return 0;
}

