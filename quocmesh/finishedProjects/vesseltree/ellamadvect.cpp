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
 * 1D ELLAM on a single bifurcation
 * ******************************************************************* */

#include "sparseMatrices.h"
#include "solver.h"
#include "ellam_matrices.h"
#include "parameterParser.h"

// this is advection using ellam on a tree with a single bifurcation.

// #define VERBOSE 1
// #define INIT_PROFILE_SPLIT 1
// #define INIT_PROFILE_COMBN 1

template class aol::GenBandMatrix<double>;

inline double u_surr_e ( int /*timestep*/, int /*gridpoint*/ ) {
  return ( 1.0 );
}

int main() {
  try {

    // ---------- read parameters ----------

    aol::ParameterParser parser ( "par/ellam.par" );

    const bool
      splitting = static_cast< bool > ( parser.getInt ( "splitting" ) );   // splitting or combining of flow

    const int
      np = parser.getInt ( "np" ),                          // number of unknowns
      nd = parser.getInt ( "nd" ),
      ne = parser.getInt ( "ne" ),

      timesteps = parser.getInt ( "timesteps" ),            // number of time steps to compute
      save_every = parser.getInt ( "save_every" );

    const double
      lp = parser.getDouble ( "lp" ),                       // lengths
      ld = parser.getDouble ( "ld" ),
      le = parser.getDouble ( "le" ),

      hp = lp / ( np - 1.0 ),                               // grid spacing
      hd = ld / ( nd * 1.0 ),
      he = le / ( ne * 1.0 ),

      Ap = parser.getDouble ( "Ap" ),                       // cross section area
      Ad = parser.getDouble ( "Ad" ),
      Ae = parser.getDouble ( "Ae" ),

      theta_d = parser.getDouble ( "theta" ),               // flow splitting ratio
      theta_e = 1.0 - theta_d,

      vp = parser.getDouble ( "vp" ),                       // velocities
      vd = theta_d * vp * Ap / Ad,
      ve = theta_e * vp * Ap / Ae,

      tau  = 1.0 / parser.getDouble ( "oneovertau" ),       // time step

      rhop = vp * tau / hp,                                 // velocity in grid cells / time step
      rhod = vd * tau / hd,
      rhoe = ve * tau / he;

#ifdef VERBOSE
    cerr << "Ap = " << Ap << ", Ad = " << Ad << ", Ae = " << Ae << endl;
    cerr << "lp = " << lp << ", ld = " << ld << ", le = " << le << endl;
    cerr << "hp = " << hp << ", hd = " << hd << ", he = " << he << endl;
    cerr << "np = " << np << ", nd = " << nd << ", ne = " << ne << endl;
    cerr << "vp = " << vp << ", vd = " << vd << ", ve = " << ve << endl;
    cerr << "rhop = " << rhop << ", rhod = " << rhod << ", rhoe = " << rhoe << endl;
    cerr << "theta_d = " << theta_d << ", theta_e = " << theta_e << endl;
#endif

    // ---------- set up matrices and vectors ----------

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
    effp ( np ), effd ( nd ), effe ( ne );

    aol::Vector<double> Rp_in ( timesteps + 1 );

    up_old.setZero();  up_new.setZero();  ud_old.setZero();
    ud_new.setZero();  ue_old.setZero();  ue_new.setZero();
    bprhs.setZero();   bdrhs.setZero();   berhs.setZero();

    cerr << "Assembling LHS matrices ... ";
    if ( splitting ) {
      set_M_fl_matrix ( hp, Ap, M_p );
      set_M_ol_matrix ( hd, Ad, M_d );
      set_M_ol_matrix ( he, Ae, M_e );
    } else {
      set_M_fl_matrix ( hp, Ap, M_p );
      set_M_of_matrix ( hd, Ad, M_d );
      set_M_of_matrix ( he, Ae, M_e );
    }

    cerr << " assembling RHS matrices ... ";
    if ( splitting ) {
      set_Me_fl_matrix ( hp, rhop, Ap, Me_p );
      set_Me_ol_matrix ( hd, rhod, Ad, Me_d );
      set_Me_ol_matrix ( he, rhoe, Ae, Me_e );
    } else {
      set_Me_fl_matrix ( hp, rhop, Ap, Me_p );
      set_Me_of_matrix ( hd, rhod, Ad, Me_d );
      set_Me_of_matrix ( he, rhoe, Ae, Me_e );
    }

    cerr << "assembling coupling matrices ... " ;
    aol::SparseMatrix<double>
    M_u_d ( np, nd ),
    M_l_d ( nd, np ),
    M_u_e ( np, ne ),
    M_l_e ( ne, np ),

    Me_u_d ( np, nd ),
    Me_l_d ( nd, np ),
    Me_u_e ( np, ne ),
    Me_l_e ( ne, np );

    if ( splitting ) {

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

    } else {

      add_M_combn_bif_terms ( hd, Ad, M_p );
      add_M_combn_bif_terms ( he, Ae, M_p );
      add_Me_combn_bif_terms ( rhop,  hd, rhod, Ad, Me_p );
      add_Me_combn_bif_terms ( rhop,  he, rhoe, Ae, Me_p );

      set_C_combn_u_matrix  ( hd, Ad, M_u_d );
      set_C_combn_u_matrix  ( he, Ae, M_u_e );
      set_C_combn_l_matrix  ( hd, Ad, M_l_d );
      set_C_combn_l_matrix  ( he, Ae, M_l_e );

      set_Ce_combn_u_matrix ( rhop, hd, rhod, Ad, Me_u_d );
      set_Ce_combn_u_matrix ( rhop, he, rhoe, Ae, Me_u_e );
      set_Ce_combn_l_matrix ( hd, rhod, Ad, Me_l_d );
      set_Ce_combn_l_matrix ( he, rhoe, Ae, Me_l_e );

    }

    cerr << " done." << endl;

    cerr << "Setting up MultiVectors ... ";
    aol::MultiVector<double> u_old ( 0, 0 ), u_new ( 0, 0 ), b_rhs ( 0, 0 ), dwarf ( 0, 0 ), eff ( 0, 0 );
    u_old.appendReference ( up_old );  u_old.appendReference ( ud_old );  u_old.appendReference ( ue_old );
    u_new.appendReference ( up_new );  u_new.appendReference ( ud_new );  u_new.appendReference ( ue_new );
    b_rhs.appendReference ( bprhs  );  b_rhs.appendReference ( bdrhs  );  b_rhs.appendReference ( berhs  );
    dwarf.appendReference ( dwarfp );  dwarf.appendReference ( dwarfd );  dwarf.appendReference ( dwarfe );
    eff.appendReference (   effp   );  eff.appendReference (   effd   );  eff.appendReference (   effe   );

    cerr << " setting up BlockOps ... ";

    aol::BlockOp< double > M_block ( 3, 3 );
    M_block.setReference ( 0, 0, M_p );        M_block.setReference ( 0, 1, M_u_d );       M_block.setReference ( 0, 2, M_u_e );
    M_block.setReference ( 1, 0, M_l_d );      M_block.setReference ( 1, 1, M_d );
    M_block.setReference ( 2, 0, M_l_e );                                               M_block.setReference ( 2, 2, M_e );

    aol::CGInverse< aol::MultiVector< double > > M_inverse ( M_block );
    M_inverse.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    M_inverse.setQuietMode ( true );

    aol::BlockOp< double > Me_block ( 3, 3 );
    Me_block.setReference ( 0, 0, Me_p );      Me_block.setReference ( 0, 1, Me_u_d );     Me_block.setReference ( 0, 2, Me_u_e );
    Me_block.setReference ( 1, 0, Me_l_d );    Me_block.setReference ( 1, 1, Me_d );
    Me_block.setReference ( 2, 0, Me_l_e );                                             Me_block.setReference ( 2, 2, Me_e );

#ifdef VERBOSE
    if ( np < 10 ) {
      M_block.getReference ( 0, 0 ).print ( cout );
      M_block.getReference ( 0, 1 ).print ( cout );
      M_block.getReference ( 0, 2 ).print ( cout );
      M_block.getReference ( 1, 0 ).print ( cout );
      M_block.getReference ( 1, 1 ).print ( cout );
      M_block.getReference ( 2, 0 ).print ( cout );
      M_block.getReference ( 2, 2 ).print ( cout );

      Me_block.getReference ( 0, 0 ).print ( cout );
      Me_block.getReference ( 0, 1 ).print ( cout );
      Me_block.getReference ( 0, 2 ).print ( cout );
      Me_block.getReference ( 1, 0 ).print ( cout );
      Me_block.getReference ( 1, 1 ).print ( cout );
      Me_block.getReference ( 2, 0 ).print ( cout );
      Me_block.getReference ( 2, 2 ).print ( cout );
    }
#endif

    cerr << " done. " << endl;


    // ---------- initial and boundary (inflow) data  ----------

    // initial values:
    cerr << "Setting initial values ... ";

    if ( splitting ) {

      for ( int i = 0; i < np; i++ ) up_new[i] = 1.0;
      for ( int i = 0; i < nd; i++ ) ud_new[i] = 1.0;
      for ( int i = 0; i < ne; i++ ) ue_new[i] = 1.0 - ( i + 1 ) / ( ne * 1.0 );

#ifdef INIT_PROFILE_SPLIT
      if ( np > 100 ) {
        for ( int i = 20; i < 40; i++ ) {
          up_new[i] = 1.0;
        }
        for ( int i = 0; i < 20; i++ ) {
          up_new[80 - i] = 1.0 - i / 20.0;
          up_new[80 + i] = 1.0 - i / 20.0;
        }
      }
#endif
    } else {

      for ( int i = 0; i < np; i++ ) up_new[i] = 0.0;
      for ( int i = 0; i < nd; i++ ) ud_new[i] = 0.0; //1.0 - i / ( nd - 1.0 );
      for ( int i = 0; i < ne; i++ ) ue_new[i] = 0.0; // (1.0 * i) / (ne - 1.0);

#ifdef INIT_PROFILE_COMBN
      if ( nd > 100 ) {
        for ( int i = 20; i < 40; i++ ) {
          ud_new[i] = 1.0;
        }
        for ( int i = 0; i < 20; i++ ) {
          ud_new[80 - i] = 1.0 - i / 20.0;
          ud_new[80 + i] = 1.0 - i / 20.0;
        }
      }
      if ( ne > 100 ) {
        for ( int i = 20; i < 40; i++ ) {
          ue_new[i] = 1.0;
        }
        for ( int i = 0; i < 20; i++ ) {
          ue_new[80 - i] = 1.0 - i / 20.0;
          ue_new[80 + i] = 1.0 - i / 20.0;
        }
      }
#endif
    }

    // inflow data:
    for ( int t = 0; t < timesteps + 1; t++ ) {
      Rp_in[t] =  1.0 - 0.5 * ( 1.0 - cos ( ( 1.0 * t ) / 20.0 ) );
    }
    cerr << " done." << endl;

    aol::StopWatch timer;
    timer.start();


    // ---------- time stepping ----------

    for ( int t = 0; t < timesteps; t++ ) {
#ifdef VERBOSE
      cerr << "time step " << t;
#endif
      // determine temperature on outflow segment.
      //       // TEST OUTPUT
      //       for(int i = 0 ; i < ne ; ++i){
      //                cerr << u_new[2][i] << " " << (ne - i - 1.0) << " " << u_new[2][i] * ( ne * 1.0 ) / ( ne - 1.0 - i ) << endl;
      //       }
      //       cerr << endl;

      u_old = u_new;

      b_rhs.setZero();
      dwarf.setZero();
      eff.setZero();

      // compute rhs source term:
      // 1. inflow terms
      if ( splitting ) {

        bprhs[0] += Ap * hp * ( Rp_in[t] * ( rhop / 2.0 - rhop * rhop / 3.0 ) + Rp_in[t+1] * ( rhop / 2.0 - rhop * rhop / 6.0 ) );
        bprhs[1] += Ap * hp * ( Rp_in[t] *   rhop * rhop / 3.0                + Rp_in[t+1] *   rhop * rhop / 6.0                );

      } else {

//         bdrhs[0] += Ad * hd * ( 1.0 * ( rhod / 2.0 - rhod*rhod / 3.0 ) + 1.0 * ( rhod / 2.0 - rhod*rhod / 6.0 ) );
//         bdrhs[1] += Ad * hd * ( 1.0 *   rhod*rhod / 3.0                + 1.0 *   rhod*rhod / 6.0                );

        bdrhs[0] += Ad * hd * ( Rp_in[t] * ( rhod / 2.0 - rhod * rhod / 3.0 ) + Rp_in[t+1] * ( rhod / 2.0 - rhod * rhod / 6.0 ) );
        bdrhs[1] += Ad * hd * ( Rp_in[t] *   rhod * rhod / 3.0                + Rp_in[t+1] *   rhod * rhod / 6.0                );

        berhs[0] += Ae * he * ( Rp_in[t] * ( rhoe / 2.0 - rhoe * rhoe / 3.0 ) + Rp_in[t+1] * ( rhoe / 2.0 - rhoe * rhoe / 6.0 ) );
        berhs[1] += Ae * he * ( Rp_in[t] *   rhoe * rhoe / 3.0                + Rp_in[t+1] *   rhoe * rhoe / 6.0                );

      }

      // 2. pseudo-source for outflow/inflow segment
      if ( splitting ) {
        aol::CompositeOp< aol::Vector<double> > M_sc, Me_sc;
        aol::DiagonalMatrix<double> diag_m ( ne ), diag_me ( ne );

        for ( int i = 0; i < ne; i++ ) {
          diag_m.set ( i, i, ( i == ( ne - 1 ) ? 0.0 : ne / ( ne - i - 1.0 ) ) );
          diag_me.set ( i, i, ne / ( ne - i - 1.0 + rhoe ) );
        }

        M_sc.appendReference (  ( *M_block.getPointer ( 2, 2 ) ) );
        M_sc.appendReference (  diag_m );

        Me_sc.appendReference ( ( *Me_block.getPointer ( 2, 2 ) ) );
        Me_sc.appendReference ( diag_me );

        M_sc.applyAdd ( u_old[2], dwarfe );
        Me_sc.applyAdd ( u_old[2], dwarfe );
        dwarfe *= ve * tau / ( - 2.0 * le );

        berhs += dwarfe;
        dwarfe.setZero();

      } else {

//           for(int i = 0; i < ne; i++){
//             effe[i] = ( ne + 1.0 ) * ve / ( ne * le) * u_surr_e(t, i) ;
//           }

//          // use constant during time step or 2 pt temporal quadrature.
//          M_block.get(2,2) -> applyAdd(effe, dwarfe);
//          Me_block.get(2,2)-> applyAdd(effe, dwarfe);
//          // dwarfe *= ( tau  );
//          dwarfe *= ( tau / 2 );

//          berhs += dwarfe;
//          effe.setZero();
//          dwarfe.setZero();

      }


      // 3. source term
      if ( splitting ) {

        if ( nd > 60 ) {
          for ( int i = 10; i < 50; i++ ) {
            effd[i] = + 0.3;
          }
        }

        M_block.applyAdd ( eff, dwarf );
        Me_block.applyAdd ( eff, dwarf );
        dwarf *= ( tau / 2 );

        b_rhs += dwarf;
        dwarf.setZero();

      } else {

        //          if(nd > 60){
        //            for(int i = 10; i < 50; i++){
        //              effd[i] = -0.3;
        //            }
        //          }

        //         M_block.applyAdd(eff, dwarf);
        //         Me_block.applyAdd(eff, dwarf);
        //         dwarf *= (tau / 2);

        //         b_rhs += dwarf;
        //         dwarf.setZero();

      }


      // only save every save_every timesteps
      if ( ! ( t % save_every ) ) {
        plot_ellam ( t, u_new[0], u_new[1], u_new[2] );
      }

#ifdef VERBOSE
      cerr << ": "; // printed after "time step t" -- now preparations are complete.
#endif

      // compute time step
      Me_block.apply ( u_old, dwarf );
      dwarf += b_rhs;
      M_inverse.apply ( dwarf, u_new );
#ifdef VERBOSE
      cerr << endl;
#endif

    }

    timer.stop();
    cerr << "took " << timer.elapsedCpuTime() << " seconds." << endl;

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
  return ( 0 );
}

