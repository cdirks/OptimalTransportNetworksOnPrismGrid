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
 * Convergence experiment for 1D ELLAM on a single bifurcation
 * ******************************************************************* */

#include "sparseMatrices.h"
#include "solver.h"
#include "ellam_matrices.h"

#define test_3

template class aol::GenBandMatrix<double>;

double uzero ( double x ) {
  return ( exp ( - ( x - 0.5 ) * ( x - 0.5 ) / 0.125 ) );
}


int main() {
  try {
    int i = 7;
    //    int j = 0;
    // for(int i = 4; i < 11; i++){
    //    for(int j = 0; j < 4; j++){
    const int
      np = ( 1 << i ) + 1,
      nd = np - 1,
      ne = np - 1;

#ifdef test_1a
    const double
      lp = 1.0,
      ld = 1.0,
      le = 1.0,
      hp = lp / ( static_cast<double> ( np ) - 1 ),
      hd = ld / ( static_cast<double> ( nd ) ),
      he = le / ( static_cast<double> ( ne ) ),
      Ap = 1.0,
      Ad = 0.5,
      Ae = 0.5,
      theta_d = 0.5,
      theta_e = 1.0 - theta_d,
      vp = 0.8,
      vd = theta_d * vp * Ap / Ad,
      ve = theta_e * vp * Ap / Ae,
      tau = hp,
      rhop = vp * tau / hp,
      rhod = vd * tau / hd,
      rhoe = ve * tau / he;

    const int timesteps = ( ( 1 << i ) + ( 1 << ( i - 2 ) ) );
#endif
#ifdef test_1b
    const double
      lp = 1.0,
      ld = 1.0,
      le = 1.0,
      hp = lp / ( static_cast<double> ( np ) - 1 ),
      hd = ld / ( static_cast<double> ( nd ) ),
      he = le / ( static_cast<double> ( ne ) ),
      Ap = 1.0,
      Ad = 0.5,
      Ae = 0.5,
      theta_d = 0.5,
      theta_e = 1.0 - theta_d,
      vp = 0.8,
      vd = theta_d * vp * Ap / Ad,
      ve = theta_e * vp * Ap / Ae,
      tau = hp, // (1 << j), //  * hp * 32, // hp,
      rhop = vp * tau / hp,
      rhod = vd * tau / hd,
      rhoe = ve * tau / he;

    const int timesteps = 2 * ( ( 1 << i ) + ( 1 << ( i - 2 ) ) );
#endif
#ifdef test_2
    const double
      lp = 1.0,
      ld = 1.0,
      le = 1.0,
      hp = lp / ( static_cast<double> ( np ) - 1 ),
      hd = ld / ( static_cast<double> ( nd ) ),
      he = le / ( static_cast<double> ( ne ) ),
      Ap = 1.0,
      Ad = 1.0,
      Ae = 1.0,
      theta_d = 0.75,
      theta_e = 1.0 - theta_d,
      vp = 0.8,
      vd = theta_d * vp * Ap / Ad,
      ve = theta_e * vp * Ap / Ae,
      tau = hp, // (1 << j), //  * hp * 32, // hp,
      rhop = vp * tau / hp,
      rhod = vd * tau / hd,
      rhoe = ve * tau / he;

    const int timesteps = ( ( 1 << i ) + ( 1 << ( i - 2 ) ) );
#endif
#ifdef test_3
    const double
      lp = 1.15,
      ld = 1.34,
      le = 0.85,
      hp = lp / ( static_cast<double> ( np ) - 1 ),
      hd = ld / ( static_cast<double> ( nd ) ),
      he = le / ( static_cast<double> ( ne ) ),
      Ap = 1.23,
      Ad = 1.47,
      Ae = 0.82,
      theta_d = 0.63,
      theta_e = 1.0 - theta_d,
      vp = 0.8,
      vd = theta_d * vp * Ap / Ad,
      ve = theta_e * vp * Ap / Ae,
      tau = hp, // (1 << j), //  * hp * 32, // hp,
      rhop = vp * tau / hp,
      rhod = vd * tau / hd,
      rhoe = ve * tau / he;

    const int timesteps = ( 1 << i );
#endif


#if 0
    cerr << "Ap = " << Ap << ", Ad = " << Ad << ", Ae = " << Ae << endl;
    cerr << "lp = " << lp << ", ld = " << ld << ", le = " << le << endl;
    cerr << "hp = " << hp << ", hd = " << hd << ", he = " << he << endl;
    cerr << "np = " << np << ", nd = " << nd << ", ne = " << ne << endl;
    cerr << "vp = " << vp << ", vd = " << vd << ", ve = " << ve << endl;
    cerr << "rhop = " << rhop << ", rhod = " << rhod << ", rhoe = " << rhoe << endl;
    cerr << "theta_d = " << theta_d << ", theta_e = " << theta_e << endl;
#endif

    aol::Vector<double> errors_l2 ( timesteps ), errors_atend ( timesteps );

    aol::TriBandMatrix<double>                              // LHS matrices
    M_p ( np ),
    M_d ( nd ),
    M_e ( ne );

    aol::LQuadBandMatrix<double>                            // RHS matrices (ELLAM matrices)
    Me_p ( np ),
    Me_d ( nd ),
    Me_e ( ne );

    aol::Vector<double>
    up_old ( np ), ud_old ( nd ), ue_old ( ne ),
    up_new ( np ), ud_new ( nd ), ue_new ( ne ),
    bprhs ( np ),  bdrhs ( nd ),  berhs ( ne ),
    dwarfp ( np ), dwarfd ( nd ), dwarfe ( ne ),
    sol_p ( np ),  sol_d ( nd ),  sol_e ( ne );

    aol::Vector<double> Rp_in ( timesteps + 1 );

    up_old.setZero();
    up_new.setZero();
    ud_old.setZero();
    ud_new.setZero();
    ue_old.setZero();
    ue_new.setZero();

    bprhs.setZero();
    bdrhs.setZero();
    berhs.setZero();

    set_M_fl_matrix ( hp, Ap, M_p );
    set_M_ol_matrix ( hd, Ad, M_d );
    set_M_ol_matrix ( he, Ae, M_e );

    set_Me_fl_matrix ( hp, rhop, Ap, Me_p );
    set_Me_ol_matrix ( hd, rhod, Ad, Me_d );
    set_Me_ol_matrix ( he, rhoe, Ae, Me_e );

    aol::SparseMatrix<double>
    M_u_d ( np, nd ),
    M_l_d ( nd, np ),
    M_u_e ( np, ne ),
    M_l_e ( ne, np ),

    Me_u_d ( np, nd ),
    Me_l_d ( nd, np ),
    Me_u_e ( np, ne ),
    Me_l_e ( ne, np );

    add_M_split_bif_terms ( hd, Ad, M_p );
    add_M_split_bif_terms ( he, Ae, M_p );
    add_Me_split_bif_terms ( hp, rhop, Ap,  hd, rhod, Ad, theta_d, Me_p );
    add_Me_split_bif_terms ( hp, rhop, Ap,  he, rhoe, Ae, theta_e, Me_p );

    set_C_split_u_matrix  ( hd, Ad, M_u_d );
    set_C_split_u_matrix  ( he, Ae, M_u_e );
    set_C_split_l_matrix  ( hd, Ad, M_l_d );
    set_C_split_l_matrix  ( he, Ae, M_l_e );

    set_Ce_split_u_matrix ( hd, rhod, Ad, Me_u_d );
    set_Ce_split_u_matrix ( he, rhoe, Ae, Me_u_e );
    set_Ce_split_l_matrix ( hp, rhop, Ap,  hd, rhod, Ad, theta_d, Me_l_d );
    set_Ce_split_l_matrix ( hp, rhop, Ap,  he, rhoe, Ae, theta_e, Me_l_e );

    aol::MultiVector<double> u_old ( 0, 0 ), u_new ( 0, 0 ), b_rhs ( 0, 0 ), dwarf ( 0, 0 );
    u_old.appendReference ( up_old );
    u_old.appendReference ( ud_old );
    u_old.appendReference ( ue_old );

    u_new.appendReference ( up_new );
    u_new.appendReference ( ud_new );
    u_new.appendReference ( ue_new );

    b_rhs.appendReference ( bprhs  );
    b_rhs.appendReference ( bdrhs  );
    b_rhs.appendReference ( berhs  );

    dwarf.appendReference ( dwarfp );
    dwarf.appendReference ( dwarfd );
    dwarf.appendReference ( dwarfe );

    aol::BlockOp< double > M_block ( 3, 3 );
    M_block.setReference ( 0, 0, M_p );
    M_block.setReference ( 1, 1, M_d );
    M_block.setReference ( 2, 2, M_e );
    M_block.setReference ( 0, 1, M_u_d );
    M_block.setReference ( 0, 2, M_u_e );
    M_block.setReference ( 1, 0, M_l_d );
    M_block.setReference ( 2, 0, M_l_e );

    aol::CGInverse< aol::MultiVector< double > > solver ( M_block, 1e-18 );
    solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    solver.setQuietMode ( true );

    aol::BlockOp< double > Me_block ( 3, 3 );
    Me_block.setReference ( 0, 0, Me_p );
    Me_block.setReference ( 1, 1, Me_d );
    Me_block.setReference ( 2, 2, Me_e );
    Me_block.setReference ( 0, 1, Me_u_d );
    Me_block.setReference ( 0, 2, Me_u_e );
    Me_block.setReference ( 1, 0, Me_l_d );
    Me_block.setReference ( 2, 0, Me_l_e );

    // initial values:
    for ( int ii = 0; ii < np; ii++ )   up_new[ii] = uzero ( ii * hp );
    for ( int ii = 0; ii < nd; ii++ )   ud_new[ii] = uzero ( lp + ( ii + 1 ) * hd  );
    for ( int ii = 0; ii < ne; ii++ )   ue_new[ii] = uzero ( lp + ( ii + 1 ) * he  );

    for ( int t = 0; t < timesteps + 1; t++ ) {
      Rp_in[t] =  uzero ( -vp * t * tau );
    }

    plot_ellam ( np, u_new[0], u_new[1], u_new[2] );

    aol::StopWatch timer;
    timer.start();
    for ( int t = 1; t < timesteps + 1; t++ ) {
      u_old = u_new;

      plot_ellam ( t - 1, u_new[0], u_new[1], u_new[2] );

      bprhs[0] = Ap * hp * ( Rp_in[t-1] * ( rhop / 2.0 - rhop * rhop / 3.0 ) + Rp_in[t] * ( rhop / 2.0 - rhop * rhop / 6.0 ) );
      bprhs[1] = Ap * hp * ( Rp_in[t-1] *   rhop * rhop / 3.0                + Rp_in[t] *   rhop * rhop / 6.0                );

      Me_block.apply ( u_old, dwarf );
      dwarf += b_rhs;
      solver.apply ( dwarf, u_new );

      // compute solution

      for ( int ii = 0; ii < np; ii++ )
        sol_p[ii] = uzero ( ii * hp - vp * t * tau );
      for ( int ii = 0; ii < nd; ii++ )
        sol_d[ii] = uzero ( lp + ( ( ( ii + 1 ) * hd - vd * t * tau > 0 ) ? ( ii + 1 ) * hd - vd * t * tau : ( ii + 1 ) * hd * vp / vd - vp * t * tau ) );
      for ( int ii = 0; ii < ne; ii++ )
        sol_e[ii] = uzero ( lp + ( ( ( ii + 1 ) * he - ve * t * tau > 0 ) ? ( ii + 1 ) * he - ve * t * tau : ( ii + 1 ) * he * vp / ve - vp * t * tau ) );

      //        plot_ellam(t, sol_p, sol_d, sol_e);

      sol_p -= u_new[0];
      sol_d -= u_new[1];
      sol_e -= u_new[2];

      //        plot_ellam(t, sol_p, sol_d, sol_e);

      errors_l2[t-1] = sol_p.norm() / ( np - 1.0 ) + sol_d.norm() / ( nd - 1.0 ) + sol_e.norm() / ( ne - 1.0 ); // these are positive!
      errors_atend[t-1] = theta_e * sol_d[nd-1] + theta_d * sol_e[ne-1];
    }
    timer.stop();

    plot_ellam ( np + timesteps + 1, u_new[0], u_new[1], u_new[2] );

    double error1 = errors_l2.getMaxValue(),
    error2 = errors_atend.norm() / errors_atend.size(),
                             error = error1 + error2;

    cerr << hp << " " << tau << " " << error << " " << timer.elapsedCpuTime() <<  endl;

    cout.precision ( 10 );
    //      cout << log(hp * hp + tau) << " " << log( error ) << endl;
    cout << log ( hp ) << " " << log ( error ) << endl;
    //    }
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return ( 0 );
}


