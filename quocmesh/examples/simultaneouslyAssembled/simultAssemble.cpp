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
 * \brief a small example program for the .... (drum roll) SimultaneouslyAssembledLinCombOp!
 *
 * a small example program for the .... (drum roll) SimultaneouslyAssembledLinCombOp!
 * This operator is similar to the standard LinCombOp, but the whole assembly is done within
 * one traverse of the grid. I.e. in each element the contribution of each operator is computed
 * and added to the matrix. This is good if the traverse of the grid is expensive (like in the
 * rle-structure), but even for the simple quocmesh data structure it is faster if more than one
 * operator contributes to be assembly.
 *
 * \author Oliver Nemitz
 */

#include <configurators.h>
#include <gridBase.h>
#include <fastUniformGridMatrix.h>
#include <scalarArray.h>
#include <implicitSurfaceFEOps.h>
#include <op.h>
#include <simultaneouslyAssembled.h>
#include <mcm.h>

typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D,3> > ConfigType;
typedef aol::GaussQuadrature< double, qc::QC_1D, 7 > BoundaryQuadratureType;

using namespace aol::color;


// now the main program
int main( int, char** ) {

  try {
    cerr << "Level? ";
    int level;
    cin >> level;

    // this is required for the constructor of the single operators
    qc::GridDefinition grid( level, qc::QC_2D );
    ConfigType config( grid );

    // just some examples
    double epsilon = grid.H();
    double tau = aol::Sqr( grid.H() );

    // define a vector as weight for the mcm-ops and the result-vector
    qc::ScalarArray< double, qc::QC_2D > Phi( grid ), result( grid );
    Phi.setAll( 1. );


    // now we define some operators, which will be simultaneously assembled.
    aol::MassOp<ConfigType> massOp( grid );
    qc::MCMMassOp<ConfigType> mcmMassOp( grid, aol::ONTHEFLY, epsilon );
    mcmMassOp.setImageReference( Phi );

    aol::StiffOp<ConfigType> stiffOp( grid );
    qc::MCMStiffOp<ConfigType> mcmStiffOp( grid, aol::ONTHEFLY, epsilon );
    mcmStiffOp.setImageReference( Phi );

    qc::LaplaceBeltramiOp<ConfigType> projOp( grid, aol::ASSEMBLED, epsilon );
    projOp.setImageReference( Phi );

    // declare the main operator
    aol::SimultaneouslyAssembledLinCombOp<ConfigType> SAOp( grid );

    //! append the single operators: Therefore you have to know from which interface the operator
    //! is derived and which method it overloads (e.g. getCoeff, getCoeffMatrix, getNonLinearitiy...)
    //! You can find out this by looking into the documentation.
    //! Example: The MCMMassOp is derived from the aol::FELinScalarWeightedMassInterface-class and has a method
    //! getCoeff. Hence we call the "appendGetCoeffOpReference" of the simultaneouslyAssembled...Op
    //! and give "aol::SCALED_MASS_INT" as assemble-type.
    SAOp.appendGetCoeffOpReference( massOp, aol::SCALED_MASS_INT );
    SAOp.appendGetCoeffOpReference( mcmMassOp, aol::SCALED_MASS_INT, tau );
    SAOp.appendGetCoeffOpReference( stiffOp, aol::SCALED_STIFF_INT, tau );
    SAOp.appendGetCoeffOpReference( mcmStiffOp, aol::SCALED_STIFF_INT, tau );
    SAOp.appendGetCoeffMatrixOpReference( projOp, aol::MATRIX_STIFF_INT, tau );

    // declare a matrix in the usual way
    qc::FastUniformGridMatrix< double, qc::QC_2D > matrix( grid );

    // ... and assemble the operator
    cerr << "Assembling...";
    matrix.setZero();
    SAOp.assembleAddMatrix( matrix );          // new operator

    // now you can apply this operator to some data
    cerr << "done. Applying...";
    SAOp.apply( Phi, result );
    cerr << "done. I hope this program was helpful...\n";

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}

