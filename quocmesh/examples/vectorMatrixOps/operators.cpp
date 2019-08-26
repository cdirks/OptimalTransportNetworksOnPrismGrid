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
 *  \brief usage of some basic aol operators
 *  \author Schwen
 */

#include <vec.h>
#include <op.h>
#include <sparseMatrices.h>

int main ( int, char**) {
  try {

    aol::Vector<double> ArgumentVector ( 10 ), DestinationVector ( 10 ), SomeVector ( 10 );
    ArgumentVector[2] = 1.0;
    SomeVector[7] = aol::NumberTrait<double>::pi;

    aol::IdentityOp< aol::Vector<double> > Identity_Operator;                  // behaves as the identity
    Identity_Operator.apply ( ArgumentVector, DestinationVector );

    aol::NullOp< aol::Vector<double> > Null_Operator;                          // always maps to zero
    Null_Operator.apply ( ArgumentVector, DestinationVector );

    aol::ConstOp< aol::Vector<double> > Const_Operator ( SomeVector) ;         // always maps to same Vector
    Const_Operator.apply ( ArgumentVector, DestinationVector );

    aol::NoiseOperator<double> Noise_Operator;                                 // creates noise
    Noise_Operator.apply ( ArgumentVector, DestinationVector );                // in general, apply means: Destination = Operator ( Argument )
    Noise_Operator.applyAdd ( ArgumentVector, DestinationVector );             // whereas applyAdd means: Destination += Operator ( Argument )

    // The NoiseOperator is a so-called BiOp which means that the same instance can be applied to Vectors and MultiVectors.
    // Moreover, it can be applied to vectors of arbitrary size.
    aol::MultiVector<double> mv1 ( 5, 5 ), mv2 ( 5, 5 );
    Noise_Operator.apply( mv1, mv2 );


    // It is possible to compare whether operators are identical up to some tolerance:
    aol::SparseMatrix<double> Identity_Matrix ( 10, 10 );
    Identity_Matrix.setIdentity();

    cerr << "Identity = Identity? " << aol::compareOps ( Identity_Operator, Identity_Matrix, 10, 10, 1e-20 ) << endl; // this should be the case
    // cerr << "Noise = Noise?" << endl; cerr << aol::compareOps ( Noise_Operator, Noise_Operator, 10, 10, 1e-20 ) << endl; // these are not identical, and we do not want to see the error output


    // the CompositeOp successively applies multiple operators in the same space dimension:
    aol::CompositeOp< aol::Vector<double> > Composite_Op;
    Composite_Op.appendReference ( Const_Operator ); // applied first
    Composite_Op.appendReference ( Noise_Operator ); // applied second
    Composite_Op.apply ( ArgumentVector, DestinationVector );


    // the LinCompOp computes the linear combination of the application of different ops:
    aol::LinCombOp< aol::Vector<double> > LinComb_Op;
    LinComb_Op.appendReference ( Identity_Operator );      // implicitely with factor 1
    LinComb_Op.appendReference ( Noise_Operator, 0.01 );   // Now LinComp_Op = Identity_Operator + 0.01 * Noise_Operator
    LinComb_Op.apply( ArgumentVector, DestinationVector );


    // the BlockOp can be used in combination with MultiVectors. All unset entries are treated as zero blocks.
    aol::BlockOp<double> Block_Op( 2, 2 );
    Block_Op.setReference ( 0, 0, Identity_Operator );
    Block_Op.setReference ( 1, 1, Composite_Op );

    aol::MultiVector<double> ArgumentMultiVector ( 2, 10 ), DestinationMultiVector ( 2, 10 );
    Block_Op.apply( ArgumentMultiVector, DestinationMultiVector );

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
