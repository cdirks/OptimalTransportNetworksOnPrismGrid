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
 * tissue_vessel_geometry is a class storing data and providing
 * functions for which both tissue and vessels are needed.
 * It also provides functions for computing mutual source terms
 * for the coupled advection-diffusion problem.
 * ******************************************************************* */

#ifndef __TISSUE_VESSEL_GEOMETRY_H
#define __TISSUE_VESSEL_GEOMETRY_H

#include "vesseltree.h"
#include "gridBase.h"
#include "epswriter.h"
#include "quoc.h"


template <typename real_type>
class tissue_vessel_geometry {
public:
  // data:
  vesseltree<real_type>               vessels;
  qc::GridDefinition                  tissue;

  //  move this into geomtree
  unsigned int                        A_narcs;
  unsigned int                        V_narcs;
  int                                 vessellevel;

  aol::Vector<int>                    A_lenv;
  aol::Vector<int>                    V_lenv;

  aol::Vector<int>                    A_term;
  aol::BitVector                       A_is_term;
  aol::Vector<int>                    V_term;
  aol::BitVector                       V_is_term;

  aol::MultiVector<real_type>         atree_data;
  aol::MultiVector<real_type>         vtree_data;
  qc::ScalarArray<real_type, qc::QC_2D>  tissue_data;
  real_type                           tissue_h;
  real_type                           smallest_vessel_h;

  static const int ACOL, VCOL, TOPVIEW_RD;

  tissue_vessel_geometry() :
      tissue ( 0 ),
      vessellevel ( 0 ),
      A_term ( 0 ),
      V_term ( 0 ),
      V_is_term ( 0 ),
      atree_data ( ),
      vtree_data ( ),
      tissue_data ( tissue ),
      tissue_h ( tissue.H() ),
      smallest_vessel_h ( 0.0 ) {
    // workaround!
  }

  tissue_vessel_geometry ( char* atree_fn, char* vtree_fn, int _tissuelevel, int _vessellevel ) :
      tissue ( _tissuelevel, qc::QC_2D ),
      vessellevel ( _vessellevel ),
      A_term ( 0 ),
      A_is_term ( 0 ),
      V_term ( 0 ),
      V_is_term ( 0 ),
      atree_data ( ),
      vtree_data ( ),
      tissue_data ( tissue ),
      tissue_h ( tissue.H() ),
      smallest_vessel_h ( 2.0 ) {
    // set up discretization for tissue in initialization list

    vessels.read_trees ( atree_fn, vtree_fn );

    discretize_vessels();

    vessels.update_mutual_grid_segment_information ( tissue );

  }

  tissue_vessel_geometry ( char* atree_fn, char* vtree_fn, int _tissuelevel ) :
      tissue ( _tissuelevel, qc::QC_2D ),
      vessellevel ( _tissuelevel ),
      A_term ( 0 ),
      A_is_term ( 0 ),
      V_term ( 0 ),
      V_is_term ( 0 ),
      atree_data ( ),
      vtree_data ( ),
      tissue_data ( tissue ),
      tissue_h ( tissue.H() ),
      smallest_vessel_h ( 2.0 ) {
    // set up discretization for tissue in initialization list

    vessels.read_trees ( atree_fn, vtree_fn );

    discretize_vessels();

    vessels.update_mutual_grid_segment_information ( tissue );

  }


  void save();
  void plot_timestep();
  void save_timestep();

  real_type evaluate_aseg_at_lambda ( int seg, real_type lambda ) {
    return ( evaluate_mvec_on_atree_at_lambda ( atree_data, seg, lambda ) );
  }

  real_type evaluate_vseg_at_lambda ( int seg, real_type lambda ) {
    return ( evaluate_mvec_on_vtree_at_lambda ( vtree_data, seg, lambda ) );
  }

  real_type evaluate_mvec_on_atree_at_lambda ( aol::MultiVector<real_type> &mvec, int seg, real_type lambda );
  real_type evaluate_mvec_on_vtree_at_lambda ( aol::MultiVector<real_type> &mvec, int seg, real_type lambda );

  void compute_data_phi_products ( qc::ScalarArray<real_type, qc::QC_2D> &rhs ) {
    rhs.setZero();
    add_mvec_phi_products_onetree ( atree_data, rhs, true  );
    add_mvec_phi_products_onetree ( vtree_data, rhs, false );
  }

  void add_mvec_phi_products ( aol::MultiVector<real_type> &A_mvec, aol::MultiVector<real_type> &V_mvec, qc::ScalarArray<real_type, qc::QC_2D> &rhs ) {
    add_mvec_phi_products_onetree ( A_mvec, rhs, true  );
    add_mvec_phi_products_onetree ( V_mvec, rhs, false );
  }

  void add_mvec_phi_products_onetree ( aol::MultiVector<real_type> &mvec, qc::ScalarArray<real_type, qc::QC_2D> &rhs, bool on_atree );


  void add_atree_outflow_rhs ( aol::MultiVector<real_type> &rhs, aol::MultiVector<real_type> &mutsrc, const aol::BlockOp< real_type, aol::Matrix<real_type> > &M_block, const aol::BlockOp< real_type, aol::Matrix<real_type> > &Me_block, const real_type tau, const real_type C_tis, const real_type C_ves );

  void add_vtree_inflow_rhs ( aol::MultiVector<real_type> &rhs, aol::MultiVector<real_type> &mutsrc, const aol::BlockOp< real_type, aol::Matrix<real_type> > &M_block, const aol::BlockOp< real_type, aol::Matrix<real_type> > &Me_block, const real_type tau, const real_type C_tis, const real_type C_ves );

  void save_advect_step ( char* ofn, aol::colorTrans colTrans, real_type zzoom, real_type color_scale, aol::FullMatrix<real_type> &rot, aol::Vector<real_type> &shift, bool temperature = false, real_type min_speed = 0.0, bool do_save_atree = true, bool do_save_vtree = true );

  void save_advect_step_topview_recursive ( char* ofn, aol::colorTrans colTrans, const bool temperatures = false, const real_type min_speed = 0.0, const real_type ctr = 0.5, const real_type scl = 2.0 );

  void add_ves_tis_transfer_rhs_onetree ( aol::MultiVector<real_type> &brhsmv, aol::MultiVector<real_type> &mutsrc, const aol::BlockOp< real_type, aol::Matrix<real_type> > &M_block, const aol::BlockOp< real_type, aol::Matrix<real_type> > &Me_block, const real_type tau, const real_type C_tis, const real_type C_ves, const real_type perc, const bool on_atree );

  void add_ves_tis_transfer_rhs ( aol::MultiVector<real_type> &A_brhsmv, aol::MultiVector<real_type> &A_mutsrc, const aol::BlockOp< real_type, aol::Matrix<real_type> > &A_M_block, const aol::BlockOp< real_type, aol::Matrix<real_type> > &A_Me_block, aol::MultiVector<real_type> &V_brhsmv, aol::MultiVector<real_type> &V_mutsrc, const aol::BlockOp< real_type, aol::Matrix<real_type> > &V_M_block, const aol::BlockOp< real_type, aol::Matrix<real_type> > &V_Me_block, const real_type tau, const real_type C_tis, const real_type C_ves, const real_type perc ) {

    add_ves_tis_transfer_rhs_onetree ( A_brhsmv, A_mutsrc, A_M_block, A_Me_block, tau, C_tis, C_ves, perc, true  );
    add_ves_tis_transfer_rhs_onetree ( V_brhsmv, V_mutsrc, V_M_block, V_Me_block, tau, C_tis, C_ves, perc, false );
  }

  void add_corresponding_tissue_sources ( aol::MultiVector<real_type> &A_mutsrc, aol::MultiVector<real_type> &V_mutsrc, qc::ScalarArray<real_type, qc::QC_2D> &tis_rhs, const real_type C_tis, const real_type C_ves );

  void save_all_data ( char* ofn, const int prec = 6 );
  void read_all_data ( char* ifn );

  void adapt_speeds_up_from_termsegm() {
    // workaround! assumptions: see geomtree.cpp!
    vessels.atree.adapt_speeds_up_recurse ( 0 );
    vessels.vtree.adapt_speeds_up_recurse ( 0 );
  }

protected:

  void save_advect_step_atree ( aol::EpsWriter &saver, real_type zzoom, real_type color_scale, aol::FullMatrix<real_type> &rot, aol::Vector<real_type> &shift, bool temperatures, real_type min_speed );
  void save_advect_step_vtree ( aol::EpsWriter &saver, real_type zzoom, real_type color_scale, aol::FullMatrix<real_type> &rot, aol::Vector<real_type> &shift, bool temperatures, real_type min_speed );

  /*
    void save_advect_step_topview_atree( aol::EpsWriter &saver, real_type color_scale, bool temperatures, real_type min_speed );
    void save_advect_step_topview_vtree( aol::EpsWriter &saver, real_type color_scale, bool temperatures, real_type min_speed );
  */

  void save_advect_step_topview_atree_recurse ( aol::EpsWriter &saver, aol::EpsWriter &alpha_saver, int cur_arc, const bool temperatures, const real_type min_speed, const real_type ctr, const real_type scl );
  void save_advect_step_topview_vtree_recurse ( aol::EpsWriter &saver, aol::EpsWriter &alpha_saver, int cur_arc, const bool temperatures, const real_type min_speed, const real_type ctr, const real_type scl );

  void _sas_trapezoid ( aol::Vector<real_type> &g1, aol::Vector<real_type> &g2, aol::Vector<real_type> &g3, aol::Vector<real_type> &g4, aol::Vector<real_type> &shift, aol::FullMatrix<real_type> &rot, real_type v_avg, real_type color_scale, aol::EpsWriter &saver );
  void _sas_segment ( aol::Vector<real_type> &g1, aol::Vector<real_type> &g2, aol::Vector<real_type> &shift, aol::FullMatrix<real_type> &rot, int rd, int acol, aol::EpsWriter &saver );

  void discretize_vessels();

  real_type temp_surr_vtree ( const int arc, const int gridpt );

private:
  tissue_vessel_geometry ( tissue_vessel_geometry &/*other*/ ) :
    tissue ( 0 ), atree_data ( ), vtree_data ( ) {} // no copy constructor

};

#endif


