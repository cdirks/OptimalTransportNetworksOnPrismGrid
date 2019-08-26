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

/* *******************************************************************
 * 1D ELLAM matrices for the case 0 <= rho <= 1
 * function to set up block operator for 1D ELLAM on vessel trees
 * ******************************************************************* */

#ifndef __ELLAM_MATRICES_H
#define __ELLAM_MATRICES_H

#include <quocMatrices.h>

#include "vesseltree.h"

// TODO: see which parameters for the functions are really used. remove others.
// use template instead of double.

namespace {

inline double I11 ( const double h, const double rho ) {
  return ( h * ( - rho*rho*rho / 6.0 + rho*rho / 2.0 - rho / 2.0 + 1.0 / 6.0 ) );
}

inline double I21 ( const double h, const double rho ) {
  return ( h * ( rho*rho*rho / 6.0 - rho / 2.0 + 1.0 / 3.0 ) );
}

inline double I22 ( const double h, const double rho ) {
  return ( h * ( rho*rho*rho / 6.0 - rho*rho + rho ) );
}

inline double I23 ( const double h, const double rho ) {
  return ( h * ( rho*rho*rho / 6.0 - rho / 2.0 + 1.0 / 3.0 ) );
}

inline double I32 ( const double h, const double rho ) {
  return ( h * ( - rho*rho*rho / 6.0 + rho*rho / 2.0 ) );
}

inline double I33 ( const double h, const double rho ) {
  return ( h * ( - rho*rho*rho / 6.0 - rho*rho / 2.0 + rho / 2.0 + 1.0 / 6.0 ) );
}

inline double I34 ( const double h, const double rho ) {
  return ( h * ( - rho*rho*rho / 6.0 + rho*rho / 2.0 ) );
}

inline double I44 ( const double h, const double rho ) {
  return ( h * ( rho*rho*rho / 6.0 ) );
}

// note: here we use upwind/downwind, that is upflow/downflow

inline double J22 ( const double h_upw, const double rho_upw, const double rho_downw ) {
  return ( h_upw * ( rho_upw - 1.0 / 2.0 * rho_upw * rho_downw - 1.0 / 2.0 * rho_upw * rho_upw + 1.0 / 6.0 * rho_downw * rho_upw * rho_upw ) );
}

inline double J32 ( const double h_upw, const double rho_upw, const double rho_downw ) {
  return ( rho_downw / rho_upw * I32 ( h_upw, rho_upw ) );
}

inline double J34 ( const double h_upw, const double rho_upw, const double rho_downw ) {
  return ( h_upw * ( 1.0 / 2.0 * rho_upw * rho_upw - 1.0 / 6.0 * rho_downw * rho_upw * rho_upw ) );
}

inline double J44 ( const double h_upw, const double rho_upw, const double rho_downw ) {
  return ( rho_downw / rho_upw * I44 ( h_upw, rho_upw ) );
}

}
// up/down refers to structure in the tree, *not* the direction of flow!

// LHS setup

void set_M_fl_matrix ( const double h, const double A, aol::Matrix<double> &M ) {
  // has first and last unknown
  if ( ( M.getNumRows() == M.getNumCols() ) && ( M.getNumRows() >= 3 ) ) {
    const int n = M.getNumCols();

    M.set ( 0, 0, h / 3.0 );
    M.set ( 0, 1, h / 6.0 );

    for ( int i = 1; i < ( n - 1 ); i++ ) {
      M.set ( i, i - 1, h / 6.0 );
      M.set ( i,   i, h * 2.0 / 3.0 );
      M.set ( i, i + 1, h / 6.0 );
    }

    M.set ( n - 1, n - 2, h / 6.0 );
    M.set ( n - 1, n - 1, h / 3.0 );

    M *= A;

  } else {
    cerr << "Called set_M_inflow_matrix with a non-square or insufficiently large matrix. Doing nothing. " << endl;
  }
}

void set_M_ol_matrix ( const double h, const double A, aol::Matrix<double> &M ) {
  // has only last unknown, not first one
  if ( ( M.getNumRows() == M.getNumCols() ) && ( M.getNumRows() >= 3 ) ) {
    const int n = M.getNumCols();

    M.set ( 0, 0, h * 2.0 / 3.0 );
    M.set ( 0, 1, h / 6.0 );

    for ( int i = 1; i < ( n - 1 ); i++ ) {
      M.set ( i, i - 1, h / 6.0 );
      M.set ( i,   i, h * 2.0 / 3.0 );
      M.set ( i, i + 1, h / 6.0 );
    }

    M.set ( n - 1, n - 2, h / 6.0 );
    M.set ( n - 1, n - 1, h / 3.0 );

    M *= A;

  } else {
    cerr << "Called set_M_matrix with a non-square or insufficiently large matrix. Doing nothing. " << endl;
  }
}


void set_M_of_matrix ( const double h, const double A, aol::Matrix<double> &M ) {
  // has only first unknown, not last one
  if ( ( M.getNumRows() == M.getNumCols() ) && ( M.getNumRows() >= 3 ) ) {
    const int n = M.getNumCols();

    M.set ( 0, 0, h / 3.0 );
    M.set ( 0, 1, h / 6.0 );

    for ( int i = 1; i < ( n - 1 ); i++ ) {
      M.set ( i, i - 1, h / 6.0 );
      M.set ( i,   i, h * 2.0  / 3.0 );
      M.set ( i, i + 1, h / 6.0 );
    }

    M.set ( n - 1, n - 2, h / 6.0 );
    M.set ( n - 1, n - 1, 2.0 * h / 3.0 );

    M *= A;

  } else {
    cerr << "Called set_M_matrix with a non-square or insufficiently large matrix. Doing nothing. " << endl;
  }
}


// RHS setup

void set_Me_fl_matrix ( const double h, const double rho, const double A, aol::Matrix<double> &Me ) {
  // has first and last unknown

  if ( ( Me.getNumRows() == Me.getNumCols() ) && ( Me.getNumRows() > 3 ) ) {
    const int n = Me.getNumCols();

    Me.set ( 0, 0, I23 ( h, rho ) ); // non-standard
    Me.set ( 0, 1, I11 ( h, rho ) );

    Me.set ( 1, 0, I33 ( h, rho ) + I34 ( h, rho ) ); // non-standard
    Me.set ( 1, 1, I21 ( h, rho ) + I22 ( h, rho ) + I23 ( h, rho ) );
    Me.set ( 1, 2, I11 ( h, rho ) );

    for ( int i = 2; i < ( n - 1 ); i++ ) {
      Me.set ( i, i - 2, I44 ( h, rho ) );
      Me.set ( i, i - 1, I32 ( h, rho ) + I33 ( h, rho ) + I34 ( h, rho ) );
      Me.set ( i, i  , I21 ( h, rho ) + I22 ( h, rho ) + I23 ( h, rho ) );
      Me.set ( i, i + 1, I11 ( h, rho ) );
    }

    Me.set ( n - 1, n - 3, I44 ( h, rho ) );
    Me.set ( n - 1, n - 2, I32 ( h, rho ) + I33 ( h, rho ) );    // I34 lives on \Gamma_out
    Me.set ( n - 1, n - 1, I21 ( h, rho ) );          // I22 lives on \Gamma_out

    Me *= A;

  } else {
    cerr << "Called set_Me_inflow_matrix with a non-square matrix. Doing nothing. " << endl;
  }
}

void set_Me_ol_matrix ( const double h, const double rho, const double A, aol::Matrix<double> &Me ) {
  // has only last unknown, not first one

  if ( ( Me.getNumRows() == Me.getNumCols() ) && ( Me.getNumRows() > 3 ) ) {
    const int n = Me.getNumCols();

    Me.set ( 0, 0, I21 ( h, rho ) + I22 ( h, rho ) + I23 ( h, rho ) );
    Me.set ( 0, 1, I11 ( h, rho ) );

    Me.set ( 1, 0, I32 ( h, rho ) + I33 ( h, rho ) + I34 ( h, rho ) );
    Me.set ( 1, 1, I21 ( h, rho ) + I22 ( h, rho ) + I23 ( h, rho ) );
    Me.set ( 1, 2, I11 ( h, rho ) );

    for ( int i = 2; i < ( n - 1 ); i++ ) {
      Me.set ( i, i - 2, I44 ( h, rho ) );
      Me.set ( i, i - 1, I32 ( h, rho ) + I33 ( h, rho ) + I34 ( h, rho ) );
      Me.set ( i, i  , I21 ( h, rho ) + I22 ( h, rho ) + I23 ( h, rho ) );
      Me.set ( i, i + 1, I11 ( h, rho ) );
    }

    Me.set ( n - 1, n - 3, I44 ( h, rho ) );
    Me.set ( n - 1, n - 2, I32 ( h, rho ) + I33 ( h, rho ) );    // I34 lives on \Gamma_out
    Me.set ( n - 1, n - 1, I21 ( h, rho ) );          // I22 lives on \Gamma_out

    Me *= A;

  } else {
    cerr << "Called set_M_matrix with a non-square matrix. Doing nothing. " << endl;
  }
}


void set_Me_of_matrix ( const double h, const double rho, const double A, aol::Matrix<double> &Me ) {
  // has only first unknown, not last one
  if ( ( Me.getNumRows() == Me.getNumCols() ) && ( Me.getNumRows() >= 3 ) ) {
    const int n = Me.getNumCols();

    Me.set ( 0, 0, I23 ( h, rho ) );
    Me.set ( 0, 1, I11 ( h, rho ) );

    Me.set ( 1, 0, I33 ( h, rho ) + I34 ( h, rho ) );
    Me.set ( 1, 1, I21 ( h, rho ) + I22 ( h, rho ) + I23 ( h, rho ) );
    Me.set ( 1, 2, I11 ( h, rho ) );

    for ( int i = 2; i < ( n - 1 ); i++ ) {
      Me.set ( i, i - 2, I44 ( h, rho ) );
      Me.set ( i, i - 1, I32 ( h, rho ) + I33 ( h, rho ) + I34 ( h, rho ) );
      Me.set ( i, i  , I21 ( h, rho ) + I22 ( h, rho ) + I23 ( h, rho ) );
      Me.set ( i, i + 1, I11 ( h, rho ) );
    }

    Me.set ( n - 1, n - 3, I44 ( h, rho ) );
    Me.set ( n - 1, n - 2, I32 ( h, rho ) + I33 ( h, rho ) + I34 ( h, rho ) );
    Me.set ( n - 1, n - 1, I21 ( h, rho ) + I22 ( h, rho ) + I23 ( h, rho ) );

    Me *= A;

  } else {
    cerr << "Called set_Me_inflow_matrix with a non-square matrix. Doing nothing. " << endl;
  }
}


// LHS split terms

void add_M_split_bif_terms ( const double h_down, const double A_down, aol::Matrix<double> &M ) {
  const int
  m = M.getNumRows(),
      n = M.getNumCols();

  M.add ( m - 1, n - 1,  A_down * h_down / 3.0 );
  // brown
}

void set_C_split_u_matrix ( const double h_down, const double A_down, aol::SparseMatrix<double> &C ) {
  const int m = C.getNumRows();
  C.set ( m - 1, 0,  A_down * h_down / 6.0 );
  // magenta
}

void set_C_split_l_matrix ( const double h_down, const double A_down, aol::SparseMatrix<double> &C ) {
  const int n = C.getNumCols();
  C.set ( 0, n - 1,  A_down * h_down / 6.0 );
  // orange
}


// RHS split terms

void add_Me_split_bif_terms ( const double h_up, const double rho_up, const double A_up, const double h_down, const double rho_down, const double A_down, const double theta_down, aol::Matrix<double> &Me ) {
  const int
  m = Me.getNumRows(),
      n = Me.getNumCols();

  Me.add ( m - 1, n - 2,   A_up * theta_down * J34 ( h_up, rho_up, rho_down ) );
  // violet

  Me.add ( m - 1, n - 1, ( A_up * theta_down * J22 ( h_up, rho_up, rho_down ) +
                           A_down * I23 ( h_down, rho_down ) )                    );
  // blue
}

void set_Ce_split_u_matrix ( const double h_down, const double rho_down,  const double A_down, aol::SparseMatrix<double> &C ) {
  const int m = C.getNumRows();

  C.set ( m - 1, 0,  A_down * I11 ( h_down, rho_down ) );
  // green

}

void set_Ce_split_l_matrix ( const double h_up, const double rho_up, const double A_up, const double h_down, const double rho_down, const double A_down, const double theta_down, aol::SparseMatrix<double> &C ) {
  const int n = C.getNumCols();

  C.set ( 0, n - 2,   A_up * theta_down * J44 ( h_up, rho_up, rho_down ) );
  // light green

  C.set ( 0, n - 1, ( A_up * theta_down * J32 ( h_up, rho_up, rho_down ) +
                      A_down * I33 ( h_down, rho_down ) +
                      A_down * I34 ( h_down, rho_down )   ) );
  // red

  C.set ( 1, n - 1,   A_down * I44 ( h_down, rho_down ) );
  // black
}


// LHS combine terms

void add_M_combn_bif_terms ( const double h_down, const double A_down, aol::Matrix<double> &M ) {

  M.add ( 0, 0,      A_down * h_down / 3.0 );
  // brown
}

void set_C_combn_u_matrix ( const double h_down, const double A_down, aol::SparseMatrix<double> &C ) {
  const int n = C.getNumCols();
  C.set ( 0, n - 1,     A_down * h_down / 6.0 );
  // magenta
}

void set_C_combn_l_matrix ( const double h_down, const double A_down, aol::SparseMatrix<double> &C ) {
  const int m = C.getNumRows();
  C.set ( m - 1, 0,     A_down * h_down / 6.0 );
  // orange
}


// RHS combine terms

void add_Me_combn_bif_terms ( const double rho_up, const double h_down, const double rho_down, const double A_down,  aol::Matrix<double> &Me ) {
  Me.add ( 0, 0,   ( A_down * I21 ( h_down, rho_down ) +
                     A_down * J22 ( h_down, rho_down, rho_up ) ) );
  // blue

  Me.add ( 1, 0,     A_down * J32 ( h_down, rho_down, rho_up ) );
  // red
}

void set_Ce_combn_u_matrix ( const double rho_up,  const double h_down, const double rho_down, const double A_down, aol::SparseMatrix<double> &C ) {
  const int n = C.getNumCols();

  C.set ( 0, n - 2,   A_down * I44 ( h_down, rho_down ) );
  // green

  C.set ( 0, n - 1, ( A_down * I32 ( h_down, rho_down ) +
                      A_down * I33 ( h_down, rho_down ) +
                      A_down * J34 ( h_down, rho_down, rho_up ) ) );
  // violet

  C.set ( 1, n - 1,    A_down * J44 ( h_down, rho_down, rho_up ) );
  // light green

}

void set_Ce_combn_l_matrix ( const double h_down, const double rho_down, const double A_down, aol::SparseMatrix<double> &C ) {
  const int m = C.getNumRows();

  C.set ( m - 1, 0,    A_down * I11 ( h_down, rho_down ) );
  // black
}

// other stuff:
// set up block operators for ellam advection on tree.

#ifdef VERBOSE
void check_rho ( double rho ) {
  if ( rho < 1e-5 || rho > 0.9999 )
    cerr << "rho = " << rho << " very small / big " << endl;
}
#else
void check_rho ( double /*rho*/ ) {
}
#endif


void build_blockops_recurse ( geomtree<double> &tree, aol::Vector<int> &lenv, int arc, double tau, aol::BlockOp< double, aol::Matrix<double> > &M_block, aol::BlockOp< double, aol::Matrix<double> > &Me_block, bool splitting ) {
  // treat this arc
  double h = -1e20;
  if ( arc == 0 ) {
    h = tree.length_of_arc ( arc ) / ( lenv[arc] - 1.0 );
  } else {
    h = tree.length_of_arc ( arc ) / ( lenv[arc] - 0.0 );
  }

  const double
  speed = tree.arcs[arc].get_speed(),
          Ap = tree.arcs[arc].get_cross_sect_area();

  double rho = speed * tau / h;
  check_rho ( rho );

  {
    aol::TriBandMatrix<double>* M = new aol::TriBandMatrix<double> ( lenv[arc] );
    aol::LQuadBandMatrix<double>* Me = new aol::LQuadBandMatrix<double> ( lenv[arc] );

    if ( arc == 0 ) { // on root segment, both cases

      set_M_fl_matrix ( h, Ap, ( *M ) );
      set_Me_fl_matrix ( h, rho, Ap, ( *Me ) );

    } else { // not on root segment.
      if ( splitting ) {

        set_M_ol_matrix ( h, Ap, ( *M ) );
        set_Me_ol_matrix ( h, rho, Ap, ( *Me ) );

      } else { // combining

        set_M_of_matrix ( h, Ap, ( *M ) );
        set_Me_of_matrix ( h, rho, Ap, ( *Me ) );

      }
    }

    // special treatment of terminal arcs -- outflow!

    M_block.setReference ( arc, arc,  *M );
    Me_block.setReference ( arc, arc, *Me );
  }

  //  cerr << "Setting " << arc << endl;

  const int
  dau_1 = tree.arcs[arc].dau_1,
          dau_2 = tree.arcs[arc].dau_2;

  if ( ( dau_1 != -1 ) && ( dau_2 != -1 ) ) {

#ifdef VERBOSE
    cerr << arc << " splits into " << dau_1 << " and " << dau_2 << endl;
#endif

    // arc to bifurcation node
    double // check if all these are necessary.
      hp = h,
      rhop = rho,
      speed_d = tree.arcs[dau_1].get_speed(),
      speed_e = tree.arcs[dau_2].get_speed(),
      hd = tree.length_of_arc ( dau_1 ) / ( lenv[dau_1] - 0.0 ),
      he = tree.length_of_arc ( dau_2 ) / ( lenv[dau_2] - 0.0 ),
      rhod = speed_d * tau / hd,
      rhoe = speed_e * tau / he,
      Ad = tree.arcs[dau_1].get_cross_sect_area(),
      Ae = tree.arcs[dau_2].get_cross_sect_area(),
      theta_d = ( Ad * speed_d ) / ( Ap * speed ),
      theta_e = ( Ae * speed_e ) / ( Ap * speed );

    check_rho ( rhod );
    check_rho ( rhoe );

    aol::SparseMatrix<double>* M_d_below  = new aol::SparseMatrix<double> ( lenv[dau_1], lenv[arc] );
    aol::SparseMatrix<double>* M_d_above  = new aol::SparseMatrix<double> ( lenv[arc], lenv[dau_1] );
    aol::SparseMatrix<double>* Me_d_below = new aol::SparseMatrix<double> ( lenv[dau_1], lenv[arc] );
    aol::SparseMatrix<double>* Me_d_above = new aol::SparseMatrix<double> ( lenv[arc], lenv[dau_1] );
    aol::SparseMatrix<double>* M_e_below  = new aol::SparseMatrix<double> ( lenv[dau_2], lenv[arc] );
    aol::SparseMatrix<double>* M_e_above  = new aol::SparseMatrix<double> ( lenv[arc], lenv[dau_2] );
    aol::SparseMatrix<double>* Me_e_below = new aol::SparseMatrix<double> ( lenv[dau_2], lenv[arc] );
    aol::SparseMatrix<double>* Me_e_above = new aol::SparseMatrix<double> ( lenv[arc], lenv[dau_2] );

    if ( splitting ) {

      add_M_split_bif_terms ( hd, Ad, ( *M_block.getPointer ( arc, arc ) ) );
      add_M_split_bif_terms ( he, Ae, ( *M_block.getPointer ( arc, arc ) ) );
      add_Me_split_bif_terms ( hp, rhop, Ap,  hd, rhod, Ad, theta_d,  ( *Me_block.getPointer ( arc, arc ) ) );
      add_Me_split_bif_terms ( hp, rhop, Ap,  he, rhoe, Ae, theta_e,  ( *Me_block.getPointer ( arc, arc ) ) );

      set_C_split_u_matrix  ( hd, Ad, *M_d_above );
      set_C_split_u_matrix  ( he, Ae, *M_e_above );
      set_C_split_l_matrix  ( hd, Ad, *M_d_below );
      set_C_split_l_matrix  ( he, Ae, *M_e_below );

      set_Ce_split_u_matrix ( hd, rhod, Ad, *Me_d_above );
      set_Ce_split_u_matrix ( he, rhoe, Ae, *Me_e_above );
      set_Ce_split_l_matrix ( hp, rhop, Ap,  hd, rhod, Ad, theta_d,  *Me_d_below );
      set_Ce_split_l_matrix ( hp, rhop, Ap,  he, rhoe, Ae, theta_e,  *Me_e_below );

    } else { // cominbation

      add_M_combn_bif_terms ( hd, Ad, ( *M_block.getPointer ( arc, arc ) ) );
      add_M_combn_bif_terms ( he, Ae, ( *M_block.getPointer ( arc, arc ) ) );
      add_Me_combn_bif_terms ( rhop,  hd, rhod, Ad,  ( *Me_block.getPointer ( arc, arc ) ) );
      add_Me_combn_bif_terms ( rhop,  he, rhoe, Ae,  ( *Me_block.getPointer ( arc, arc ) ) );

      set_C_combn_u_matrix  ( hd, Ad, *M_d_above );
      set_C_combn_u_matrix  ( he, Ae, *M_e_above );
      set_C_combn_l_matrix  ( hd, Ad, *M_d_below );
      set_C_combn_l_matrix  ( he, Ae, *M_e_below );

      set_Ce_combn_u_matrix ( rhop,  hd, rhod, Ad, *Me_d_above );
      set_Ce_combn_u_matrix ( rhop,  he, rhoe, Ae, *Me_e_above );
      set_Ce_combn_l_matrix ( hd, rhod, Ad, *Me_d_below );
      set_Ce_combn_l_matrix ( he, rhoe, Ae, *Me_e_below );

    }

    M_block.setReference ( dau_1, arc, *M_d_below );
    M_block.setReference ( arc, dau_1, *M_d_above );
    Me_block.setReference ( dau_1, arc, *Me_d_below );
    Me_block.setReference ( arc, dau_1, *Me_d_above );

    M_block.setReference ( dau_2, arc, *M_e_below );
    M_block.setReference ( arc, dau_2, *M_e_above );
    Me_block.setReference ( dau_2, arc, *Me_e_below );
    Me_block.setReference ( arc, dau_2, *Me_e_above );

    build_blockops_recurse ( tree, lenv, dau_1, tau, M_block, Me_block, splitting );
    build_blockops_recurse ( tree, lenv, dau_2, tau, M_block, Me_block, splitting );

  } else if ( ( dau_1 != -1 ) && ( dau_2 == -1 ) ) {
    // arc to intermediate node
    cerr << " ATTENTION: the case of intermediate nodes is not implemented correctly yet" << endl;

  } else if ( ( dau_1 == -1 ) && ( dau_2 == -1 ) ) {
    // arc to terminal node
    // do nothing else.
  } else {
    cerr << "usevesseltree_ellamadvect_blockop<real_type>::build_blockops_recurse: Houston, we have a problem ... " << endl;
  }
}


void delete_blockop_entries ( aol::BlockOp< double, aol::Matrix<double> > &delop ) {
  for ( int i = 0; i < delop.getNumRows(); i++ )
    for ( int j = 0; j < delop.getNumCols(); j++ )
      if ( delop.getPointer ( i, j ) ) {
        delete ( delop.getPointer ( i, j ) );
        delop.unset ( i, j );
      }
}


// for testing ellam:

void plot_ellam ( int t, aol::Vector<double> &vec_1, aol::Vector<double> &vec_2, aol::Vector<double> &vec_3 ) {
  char fn[32];
  char fnp[32];
  sprintf ( fn, "output/output%04d.dat", t );
  sprintf ( fnp, "output/output%04d.plot", t );
  FILE * outfile;
  FILE * outfilep;
  outfile = fopen ( fn, "w" );
  outfilep = fopen ( fnp, "w" );
#if 0
  fprintf ( outfilep, "set terminal png\nset output \"output%04d.png\"\nplot [1:%d][-0.5:1.5] \"output%04d.dat\" w l\n",
            t, vec_1.size() + vec_2.size() + vec_3.size() + 60, t );
#else
  fprintf ( outfilep, "set terminal postscript eps color solid\nset output \"output%04d.eps\"\nplot [1:%d][-0.5:1.5] \"output%04d.dat\" notitle w l \n",
            t, vec_1.size() + vec_2.size() + vec_3.size() + 60, t );
#endif
  fclose ( outfilep );

#if 0
  for ( int i = 0; i < vec_1.size(); i++ ) {
    fprintf ( outfile, "%f\n", ( vec_1[i] ) );
  }

  for ( int i = 0; i < 20; i++ ) {
    fprintf ( outfile, "-0.1\n" );
  }

  for ( int i = 0; i < vec_2.size(); i++ ) {
    fprintf ( outfile, "%f\n", ( vec_2[i] ) );
  }

  for ( int i = 0; i < 20; i++ ) {
    fprintf ( outfile, "-0.1\n" );
  }

  for ( int i = 0; i < vec_3.size(); i++ ) {
    fprintf ( outfile, "%f\n", ( vec_3[i] ) );
  }

  for ( int i = 0; i < 20; i++ ) {
    fprintf ( outfile, "-0.1\n" );
  }
#else
  for ( int i = 0; i < vec_1.size(); i++ ) {
    fprintf ( outfile, "%d %f\n", i, vec_1[i] );
  }

  fprintf ( outfile, "\n" );

  for ( int i = 0; i < vec_2.size(); i++ ) {
    fprintf ( outfile, "%d %f\n", 20 + vec_1.size() + i, vec_2[i] );
  }

  fprintf ( outfile, "\n" );

  for ( int i = 0; i < vec_3.size(); i++ ) {
    fprintf ( outfile, "%d %f\n", 40 + vec_1.size() + vec_2.size() + i, vec_3[i] );
  }
#endif
  fclose ( outfile );
}


#endif

