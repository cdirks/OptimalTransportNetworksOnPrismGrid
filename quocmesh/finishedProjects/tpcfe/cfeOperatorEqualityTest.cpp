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

// this program tests tpCFE Mass and StiffOps for whether their on-the-fly and assembled variants are equal.

#include <FEOpInterface.h>
#include <prolongation.h>
#include <restriction.h>
#include <tpCFEStandardOp.h>

#include "affine.h"


typedef double RealType;

static const RealType DIFF_COEFF_PLUS  =  1.0;
static const RealType DIFF_COEFF_MINUS = 10.0;
static const RealType TANG_A           =  0.0;
static const RealType TANG_B           =  0.0;

static const RealType tol = 1e-15;

static const tpcfe::ConstraintType AT = tpcfe::CFE_LIEHR;

static const int DEPTH = 2;

typedef tpcfe::CFEConfigurator < tpcfe::CFEGrid< RealType, AT> , aol::SparseMatrix<double> >   ConfiguratorType;
typedef ConfiguratorType::MatrixType   MatrixType;
typedef tpcfe::CFEMassOp < ConfiguratorType >                            MassOpType;
typedef tpcfe::CFEStiffOp < ConfiguratorType >                           StiffOpType;
typedef tpcfe::CFEMassOpWI < ConfiguratorType >                          WMassOpType;
typedef tpcfe::CFEStiffOpWI < ConfiguratorType >                         WStiffOpType;

int main ( int, char** ) {
  try {

    for ( int TYPE = 0; TYPE < 7; ++TYPE ) {

      tpcfe::CFEGrid < RealType, AT >  grid ( DEPTH );

      qc::ScalarArray<RealType, qc::QC_3D>
        my_interface ( grid ),        // spherical interface
        boundary_condition ( grid ),  // boundary values
        orig ( grid );                // theoretical solution U
      qc::AArray<RealType, qc::QC_3D> diff_coeff ( grid );          // diffusion coefficient at grid points


      tpcfe::AffineFunction<RealType>        affineFct;

      affineFct.setCoeffs ( DIFF_COEFF_PLUS, DIFF_COEFF_MINUS, TANG_A, TANG_B );

      affineFct.setTypeNumber ( TYPE );

      affineFct.setWidth ( ( 1 << DEPTH ) + 1 );
      std::cerr << affineFct << std::endl;

      // even if we're passing references, the structure becomes clearer if all data is assembled in the order in which it is needed.

      affineFct.createInterface ( &my_interface, &orig, &diff_coeff );


      grid.addStructureFrom ( my_interface );
      grid.detectAndInitVirtualNodes ( diff_coeff );

      // Initialize the grid
      cerr << "Comparing cfe mass / stiff operators ..." << endl;
      const int nn = grid.getNumberOfNodes();
      ConfiguratorType _dummy_config ( grid );

      {
        MatrixType matrix ( grid );

        MassOpType                    massOp_ass ( grid, aol::ASSEMBLED );
        MassOpType                    massOp_otf ( grid, aol::ONTHEFLY );
        massOp_ass._quietMode = true;        massOp_otf._quietMode = true;

        const bool mass_otf_ass = compareOps ( massOp_otf, massOp_ass, nn, nn, tol );
        matrix.setZero();
        massOp_otf.assembleAddMatrix ( matrix );
        const bool mass_otf_mat = compareOps ( massOp_otf, matrix, nn, nn, tol );
        matrix.setZero();
        massOp_ass.assembleAddMatrix ( matrix );
        const bool mass_ass_mat = compareOps ( massOp_ass, matrix, nn, nn, tol );

        cerr << "MassOp on-the-fly vs assembled:       " << mass_otf_ass << endl;
        cerr << "MassOp on-the-fly vs matrix variant:  " << mass_otf_mat << endl;
        cerr << "MassOp assembled vs matrix variant:   " << mass_ass_mat << endl;
      }

      {
        MatrixType matrix ( grid );

        StiffOpType                    stiffOp_ass ( grid, aol::ASSEMBLED );
        StiffOpType                    stiffOp_otf ( grid, aol::ONTHEFLY );
        stiffOp_ass._quietMode = true;        stiffOp_otf._quietMode = true;

        const bool stiff_otf_ass = compareOps ( stiffOp_otf, stiffOp_ass, nn, nn, tol );
        matrix.setZero();
        stiffOp_otf.assembleAddMatrix ( matrix );
        const bool stiff_otf_mat = compareOps ( stiffOp_otf, matrix, nn, nn, tol );
        matrix.setZero();
        stiffOp_ass.assembleAddMatrix ( matrix );
        const bool stiff_ass_mat = compareOps ( stiffOp_ass, matrix, nn, nn, tol );

        cerr << "StiffOp on-the-fly vs assembled:       " << stiff_otf_ass << endl;
        cerr << "StiffOp on-the-fly vs matrix variant:  " << stiff_otf_mat << endl;
        cerr << "StiffOp assembled vs matrix variant:   " << stiff_ass_mat << endl;
      }

      {
        MatrixType matrix ( grid );

        WMassOpType                    WmassOp_ass ( diff_coeff, grid, aol::ASSEMBLED );
        WMassOpType                    WmassOp_otf ( diff_coeff, grid, aol::ONTHEFLY );
        WmassOp_ass._quietMode = true;        WmassOp_otf._quietMode = true;

        const bool Wmass_otf_ass = compareOps ( WmassOp_otf, WmassOp_ass, nn, nn, tol );
        matrix.setZero();
        WmassOp_otf.assembleAddMatrix ( matrix );
        const bool Wmass_otf_mat = compareOps ( WmassOp_otf, matrix, nn, nn, tol );
        matrix.setZero();
        WmassOp_ass.assembleAddMatrix ( matrix );
        const bool Wmass_ass_mat = compareOps ( WmassOp_ass, matrix, nn, nn, tol );

        cerr << "WMassOp on-the-fly vs assembled:       " << Wmass_otf_ass << endl;
        cerr << "WMassOp on-the-fly vs matrix variant:  " << Wmass_otf_mat << endl;
        cerr << "WMassOp assembled vs matrix variant:   " << Wmass_ass_mat << endl;
      }

      {
        MatrixType matrix ( grid );

        WStiffOpType                    WstiffOp_ass ( diff_coeff, grid, aol::ASSEMBLED );
        WStiffOpType                    WstiffOp_otf ( diff_coeff, grid, aol::ONTHEFLY );
        WstiffOp_ass._quietMode = true;        WstiffOp_otf._quietMode = true;

        const bool Wstiff_otf_ass = compareOps ( WstiffOp_otf, WstiffOp_ass, nn, nn, tol );
        matrix.setZero();
        WstiffOp_otf.assembleAddMatrix ( matrix );
        const bool Wstiff_otf_mat = compareOps ( WstiffOp_otf, matrix, nn, nn, tol );
        matrix.setZero();
        WstiffOp_ass.assembleAddMatrix ( matrix );
        const bool Wstiff_ass_mat = compareOps ( WstiffOp_ass, matrix, nn, nn, tol );

        cerr << "WStiffOp on-the-fly vs assembled:       " << Wstiff_otf_ass << endl;
        cerr << "WStiffOp on-the-fly vs matrix variant:  " << Wstiff_otf_mat << endl;
        cerr << "WStiffOp assembled vs matrix variant:   " << Wstiff_ass_mat << endl;
      }

    }

  } catch ( aol::Exception &exc ) {

    exc.dump();
    return ( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );
}
