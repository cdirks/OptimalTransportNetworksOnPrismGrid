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
 * nodes and arcs used in geomtree
 * auxiliary classes for interaction between trees and grids
 * geomtree: binary geometric trees
 * vesseltree: pair of geomtrees, functions to be called on both trees
 * other general stuff
 * ******************************************************************* */

#ifndef __VESSELTREE_H
#define __VESSELTREE_H

// stl classes
#include <vector>
#include <fstream>
#include <iostream>

// quoc classes
#include <gridBase.h>
#include <array.h>
#include <vec.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <scalarArray.h>

// vesseltree classes
#include "epswriter.h"


// #define VERBOSE

/* Some definitions
   ================
*/

// #define VERBOSE_GROWING

template <typename real_type>
inline const real_type sqr ( real_type a ) {
  return ( a*a );
}

template <typename real_type>
inline const real_type nrnd ( real_type lower, real_type upper ) {
  // return normally distributed random in [lower, upper]
  return ( lower + ( upper - lower ) * ( ( static_cast<real_type> ( rand() ) / RAND_MAX ) ) );
}


/* Auxiliary Classes: Nodes and arcs
   =================================
*/

template <typename real_type>
class node {
public:
  real_type x;
  real_type y;
  signed char level; // I hope this is sufficient - else use signed short or signed int ...

  node() {
    x = 0.0; y = 0.0; level = 0;
  }

  node ( real_type _x, real_type _y, signed char _level ) {
    replaceby ( _x, _y, _level );
  }

  node ( const node& othernode ) {
    replaceby ( othernode.x, othernode.y, othernode.level );
  }

  const node& operator= ( const node &othernode ) {
    replaceby ( othernode.x, othernode.y, othernode.level );
    return ( *this );
  }

  void replaceby ( real_type _x, real_type _y ) {
    x = _x;
    y = _y;
  }

  void replaceby ( real_type _x, real_type _y, signed char _level ) {
    x = _x;
    y = _y;
    level = _level;
  }

  void set_level ( signed char _level ) {
    level = _level;
  }

  bool is_different_from ( node &othernode ) const {
    // attention: this only compares the geometric position, not the level
    return ( ( x != othernode.x ) || ( y != othernode.y ) );
  }
};

template <typename real_type>
class arc {
public:
  int                     init;        // index of initial node
  int                     term;        // index of terminal node
  int                     parent;      // index of parent arc (-1: no parent arc (root) )
  int                     dau_1;       // index of one daughter arc (-1: no daughter arc(terminal segment) )
  int                     dau_2;       // index of the other daughter arc (or intermediate node if only 2nd.)
  int                     sibling;     // index of sister arc (-1: no sister arc)
  real_type               rad;    // radius
  real_type               speed;  // speed of blood flowing through
  real_type               theta;  // flow splitting ratio to parent arc

  arc() {
    init = 0; term = 0; parent = -1; dau_1 = -1; dau_2 = -1; sibling = -1; rad = 0; speed = 0; theta = 0; value = 0;
  }

  arc ( int _init, int _term, real_type _rad, real_type _speed = 0.0, real_type _theta = 0.0, real_type _value = 0.0 ) {
    replaceby ( _init, _term, -1, -1, -1, -1, _rad, _speed, _theta, _value ); // fixme
  }

  arc ( const arc& otherarc ) {
    replaceby ( otherarc.init, otherarc.term, otherarc.parent, otherarc.dau_1, otherarc.dau_2, otherarc.sibling, otherarc.rad, otherarc.speed, otherarc.theta, otherarc.value );
  }

  const arc& operator= ( const arc &otherarc ) {
    replaceby ( otherarc.init, otherarc.term, otherarc.parent, otherarc.dau_1, otherarc.dau_2, otherarc.sibling, otherarc.rad, otherarc.speed, otherarc.theta, otherarc.value );
    return ( *this );
  }

  void replaceby ( int _init, int _term, real_type _rad, real_type _speed = 0.0, real_type _theta = 0.0, real_type _value = 0.0 ) {
    replaceby ( _init, _term, -1, -1, -1, -1, _rad, _speed, _theta, _value );
  }

  void replaceby ( int _init, int _term, int _parent, int _dau_1, int _dau_2, int _sibling, real_type _rad, real_type _speed, real_type _theta, real_type _value ) {
    init    = _init ;
    term    = _term ;
    parent  = _parent;
    dau_1   = _dau_1;
    dau_2   = _dau_2;
    sibling = _sibling;
    rad     = _rad  ;
    speed   = _speed;
    theta   = _theta;
    value   = _value;
  }

  void set_rad ( real_type _rad ) {
    rad = _rad;
  }

  void set_speed ( real_type _speed ) {
    speed = _speed;
  }

  void set_theta ( real_type _theta ) {
    theta = _theta;
  }

  void set_value ( real_type _value ) {
    value = _value;
  }

  real_type get_value() {
    return ( value );
  }

  real_type get_speed() {
    return ( speed );
  }

  real_type get_cross_sect_area() {
    return ( 2.0 * rad );
  }

  real_type evaluate_at ( real_type /*xco*/, real_type /*yco*/, aol::Vector<real_type> &/*values*/ ) {
    // project point (xco, yco) onto segment, compute local coordinate on segment, do pcw lin interpolation of values at that point.
    return ( 0.0 );
  }

private:
  real_type value; // used for pressure / heat conduction as a single flow/temperature value on the segment.

};



/* Auxiliary Classes: Store information about segments intersecting grid cells
   ===========================================================================
*/


template <typename real_type>
class segment_in_cell {
public:
  int arcno;
  real_type
  x_in, // local coordinate of grid cell where segment enters etc.
  y_in,
  lambda_in, // local coordinate on segment where grid cell entered
  x_out,
  y_out,
  lambda_out;

  segment_in_cell() :
      arcno ( -1 ),
      x_in ( -1.0 ),
      y_in ( -1.0 ),
      lambda_in ( -1.0 ),
      x_out ( 2.0 ),
      y_out ( 2.0 ),
      lambda_out ( 2.0 ) {}

  segment_in_cell ( int _arcno, real_type _x_in, real_type _y_in, real_type _x_out, real_type _y_out ) {
    replaceby ( _arcno, _x_in, _y_in, _x_out, _y_out );
  }

  segment_in_cell ( int _arcno, real_type _x_in, real_type _y_in, real_type _lambda_in, real_type _x_out, real_type _y_out, real_type _lambda_out ) {
    replaceby ( _arcno, _x_in, _y_in, _lambda_in, _x_out, _y_out, _lambda_out );
  }

  segment_in_cell ( const segment_in_cell& other_sic ) {
    replaceby ( other_sic.arcno, other_sic.x_in, other_sic.y_in, other_sic.lambda_in, other_sic.x_out, other_sic.y_out, other_sic.lambda_out );
  }

  const segment_in_cell& operator= ( const segment_in_cell &other_sic ) {
    replaceby ( other_sic.arcno, other_sic.x_in, other_sic.y_in, other_sic.lambda_in, other_sic.x_out, other_sic.y_out, other_sic.lambda_out );
    return ( *this );
  }

  void replaceby ( int _arcno, real_type _x_in, real_type _y_in, real_type _x_out, real_type _y_out ) {
    replaceby ( _arcno, _x_in, _y_in, -1.0, _x_out, _y_out, 2.0 );
  }

  void replaceby ( int _arcno, real_type _x_in, real_type _y_in, real_type _lambda_in, real_type _x_out, real_type _y_out, real_type _lambda_out ) {
    arcno = _arcno;
    x_in = _x_in;
    y_in = _y_in;
    lambda_in = _lambda_in;
    x_out = _x_out;
    y_out = _y_out;
    lambda_out = _lambda_out;
  }

  void set_lambdas ( real_type _lambda_in, real_type _lambda_out ) {
    lambda_in = _lambda_in;
    lambda_out = _lambda_out;
  }

};

template <typename real_type>
class segment_in_cell_storage {
private:
  int m;
  int n;
  int s;
  std::vector< std::vector< segment_in_cell<real_type> > > data;

public:
  segment_in_cell_storage() {
    m = 0;
    n = 0;
    s = 0;
  }

  // maybe write constructor that takes a grid.
  segment_in_cell_storage ( int _m, int _n ) {
    resize ( _m, _n );
  }

  void resize ( int _m, int _n ) {
    m = _m;
    n = _n;
    s = m * n;

    data.clear();
    std::vector< segment_in_cell<real_type> > dummy;
    for ( int i = 0; i < s; i++ ) {
      data.push_back ( dummy );
    }
  }


  void append ( const int i, const int j, segment_in_cell<real_type> sic ) {
#ifdef VERBOSE
    cerr << "segment_in_cell_storage::append: "
    << i << " " << j << ": " << sic.arcno << " " << sic.x_in << " " << sic.y_in << " " << sic.lambda_in << " "
    << sic.x_out << " " << sic.y_out << " " << sic.lambda_out ;
#endif
    // check if nonzero length inside grid cell
    if (  ( abs ( sic.x_in - sic.x_out ) + abs ( sic.y_in - sic.y_out ) ) > 1e-8 ) {

      data[m*i+j].push_back ( sic );

#ifdef VERBOSE
      cerr << endl;
    } else {
      cerr << " zero length. not appended." << endl;
#endif
    }
  }

  typename std::vector< segment_in_cell<real_type> >::iterator at_begin ( const int i, const int j ) {
    return ( data[m*i+j].begin() );
  }

  typename std::vector< segment_in_cell<real_type> >::iterator at_end ( const int i, const int j ) {
    return ( data[m*i+j].end() );
  }
};

template <typename real_type>
class cell_under_segment {
public:
  int
  ix,
  jy;
  segment_in_cell<real_type> sic;

  cell_under_segment() :
      ix ( 0 ),
      jy ( 0 ),
      sic ( NULL ) {}

  cell_under_segment ( const int _ix, const int _jy, segment_in_cell<real_type> &_sic ) :
      ix (  _ix  ),
      jy (  _jy  ),
      sic ( _sic ) {}

  cell_under_segment ( const cell_under_segment &other_cus ) {
    ix  = other_cus.ix;
    jy  = other_cus.jy;
    sic = other_cus.sic;
  }

  const cell_under_segment& operator= ( const cell_under_segment &other_cus ) {
    ix  = other_cus.ix;
    jy  = other_cus.jy;
    sic = other_cus.sic;
    return ( *this );
  }
};

template <typename real_type>
class cell_under_segment_storage {
private:
  int s;
  std::vector< std::vector< cell_under_segment<real_type> > > data;

public:
  cell_under_segment_storage() :
      s ( 0 ) {}

  void resize ( int _s ) {
    s = _s;
    data.clear();
    for ( int i = 0; i < s; i++ ) {
      std::vector< cell_under_segment<real_type> > dummy;
      data.push_back ( dummy );
    }
  }

  void append ( int ix, int jy, segment_in_cell< real_type> &sic ) {
    cell_under_segment<real_type> cus ( ix, jy, sic );
    data[ sic.arcno ].push_back ( cus );
  }

  typename std::vector< cell_under_segment<real_type> >::iterator at_begin ( const int S ) {
    return ( data[S].begin() );
  }

  typename std::vector< cell_under_segment<real_type> >::iterator at_end ( const int S ) {
    return ( data[S].end() );
  }
};


/* A geometric binary tree with arc data
   =====================================
*/

template <typename real_type>
class geomtree {
public:
  //  enum treetypes { atree, vtree, unknown };

  std::vector< node<real_type> > nodes; // nodes[0] is always the root node
  std::vector<  arc<real_type> >  arcs; // arcs[0] is always the root segment
  segment_in_cell_storage<real_type> segments_in_grid;
  cell_under_segment_storage<real_type> grid_under_segments;

  //  treetypes type;

  geomtree(); // standard constructor

  ~geomtree();

  void clear();

  // I/O
  void read_tree ( char* filename );
  void save_tree ( char fname[] );


  // service for grids
  void update_mutual_grid_segment_information ( const qc::GridDefinition &grid, bool is_atree ) {
    update_segments_in_grid ( grid );
    update_grid_in_segments ( grid, is_atree );
  }

  void update_segments_in_grid ( const qc::GridDefinition &grid );
protected:
  void update_grid_in_segments ( const qc::GridDefinition &grid, bool is_atree ); // may not be called alone. update_segments_in_grid must be called before
public:


  // service for velocity determinaton
  void get_terminal_arc_indices ( aol::Vector<int> &term_arc_indices, aol::BitVector &is_terminal_arc );
  void set_terminal_speeds ( aol::Vector<int> &term, aol::Vector<real_type> &s, int offset );
  void set_speeds_from_term();


  // service for heat conduction problem
  void compute_rhs_fphi ( qc::GridDefinition &grid, qc::ScalarArray<real_type, qc::QC_2D> &rhs );

  // other service
  void remove_slow_arcs ( geomtree<real_type> &dest_tree, real_type threshold, int arc_index, int cp_init_index );
  void set_all_speeds ( real_type value );
  //  void plot_tree_speeds( aol::EpsWriter *saver, const real_type v_min = (1.0/1024.0), const real_type v_max = 1.0 );
  void plot_tree_speeds ( aol::EpsWriter &saver, const real_type v_min = ( 1.0 / 1024.0 ), const real_type v_max = 1.0 );

  double adapt_speeds_up_recurse ( int current_arc );

  // general stuff
  int get_max_level();

  inline const real_type getXCO ( const int i ) {
    return ( nodes[i].x );
  }
  inline const real_type getYCO ( const int i ) {
    return ( nodes[i].y );
  }

  inline const real_type length_of_arc ( const int i ) {
    return ( sqrt ( sqr ( nodes[ arcs[i].term ].x - nodes[ arcs[i].init ].x ) + sqr ( nodes[ arcs[i].term ].y - nodes[ arcs[i].init ].y ) ) );
  }

  void set_par_dau_sib();

  inline node<real_type> get_midpoint_of_arc ( const int i ) {
    node<real_type> dummy ( 0.5 * ( nodes[ arcs[i].init].x + nodes[ arcs[i].term ].x ),
                            0.5 * ( nodes[ arcs[i].init].y + nodes[ arcs[i].term ].y ), nodes[ arcs[i].init].level );
    return ( dummy );
  }

  typename std::vector< node<real_type> >::iterator nodes_at_begin() {
    return ( nodes.begin() );
  }

  typename std::vector< node<real_type> >::iterator nodes_at_end() {
    return ( nodes.end() );
  }

  typename std::vector< arc<real_type> >::iterator arcs_at_begin() {
    return ( arcs.begin() );
  }

  typename std::vector< arc<real_type> >::iterator arcs_at_end() {
    return ( arcs.end() );
  }

  // internal functions
protected:
  void compute_rhs_fphi_single_term ( int ix, int jy, real_type xa, real_type ya, real_type xb, real_type yb, real_type value, qc::GridDefinition &grid, qc::ScalarArray<real_type, qc::QC_2D> &rhs );

  void update_segments_in_grid_go_E ( int current_index, real_type h, int ix_init, int jy_init, int ix_term, int jy_term, real_type xg_init, real_type yg_init, real_type xg_term, real_type yg_term );
  void update_segments_in_grid_go_W ( int current_index, real_type h, int ix_init, int jy_init, int ix_term, int jy_term, real_type xg_init, real_type yg_init, real_type xg_term, real_type yg_term );
  void update_segments_in_grid_go_N ( int current_index, real_type h, int ix_init, int jy_init, int ix_term, int jy_term, real_type xg_init, real_type yg_init, real_type xg_term, real_type yg_term );
  void update_segments_in_grid_go_S ( int current_index, real_type h, int ix_init, int jy_init, int ix_term, int jy_term, real_type xg_init, real_type yg_init, real_type xg_term, real_type yg_term );
  void update_segments_in_grid_stay ( int current_index, real_type h, int ix_init, int jy_init, int ix_term, int jy_term, real_type xg_init, real_type yg_init, real_type xg_term, real_type yg_term );

  real_type set_speeds_from_term_recurse ( int arc_index );

  void get_max_level_recurse ( int cur_arc, int cur_level, int &max_level_so_far );

  inline real_type speedfunc ( double arg ) {
    return ( abs ( log ( arg ) ) );
  }

  inline real_type speedfct (  int arcn ) {
    return ( speedfunc ( arcs[arcn].get_speed() ) );
  }


private:
  geomtree ( geomtree &other_tree ); // copy constructor

};

/* A pair of geomtrees with functions to be used on both trees
   ===========================================================
*/

template <typename real_type>
class vesseltree {

public:
  vesseltree() {
  }

  void compute_rhs_fphi ( qc::GridDefinition &grid, qc::ScalarArray<real_type, qc::QC_2D> &rhs ) {
    rhs.setZero();
    atree.compute_rhs_fphi ( grid, rhs );
    vtree.compute_rhs_fphi ( grid, rhs );
  }

  void read_trees ( char* atree_filename, char* vtree_filename ) {
    cerr << "Reading atree from file " << atree_filename << " ... ";
    atree.read_tree ( atree_filename );
    cerr << "done. Reading vtree from file " << vtree_filename << " ... ";
    vtree.read_tree ( vtree_filename );
    cerr << "done." << endl;

    set_par_dau_sib();
  }

  void save_trees ( int seed ) {
    char aname[1024], vname[1024];
    sprintf ( aname, "output/a_%d_%d.tree", seed, atree.nodes.size() );
    sprintf ( vname, "output/v_%d_%d.tree", seed, vtree.nodes.size() );

    save_trees ( aname, vname );
  }

  void save_trees ( char* atree_filename, char* vtree_filename ) {
    atree.save_tree ( atree_filename );
    vtree.save_tree ( vtree_filename );
  }

  void clear() {
    atree.clear();
    vtree.clear();
  }

  void set_terminal_speeds ( aol::Vector<int> &a_term, aol::Vector<int> &v_term, aol::Vector<real_type> &s ) {
    atree.set_terminal_speeds ( a_term, s, 0 );
    vtree.set_terminal_speeds ( v_term, s, a_term.size() ); // problem here.

    atree.set_speeds_from_term();
    vtree.set_speeds_from_term();
  }

  void set_par_dau_sib() {
    atree.set_par_dau_sib();
    vtree.set_par_dau_sib();
  }

  void update_segments_in_grid ( qc::GridDefinition &grid ) {
    atree.update_segments_in_grid ( grid );
    vtree.update_segments_in_grid ( grid );
  }

  void update_mutual_grid_segment_information ( const qc::GridDefinition &grid ) {
    atree.update_mutual_grid_segment_information ( grid, true );
    vtree.update_mutual_grid_segment_information ( grid, false );
  }

  void remove_slow_arcs ( vesseltree<real_type> &dest_vesseltree, real_type threshold ) {
    atree.set_par_dau_sib();
    vtree.set_par_dau_sib();
    dest_vesseltree.clear();

    dest_vesseltree.atree.nodes.push_back ( atree.nodes[0] );
    dest_vesseltree.atree.nodes.push_back ( atree.nodes[1] );
    dest_vesseltree.atree.arcs.push_back (  atree.arcs[0]  );
    atree.remove_slow_arcs ( dest_vesseltree.atree, threshold, 0, 1 );

    dest_vesseltree.vtree.nodes.push_back ( vtree.nodes[0] );
    dest_vesseltree.vtree.nodes.push_back ( vtree.nodes[1] );
    dest_vesseltree.vtree.arcs.push_back (  vtree.arcs[0]  );
    vtree.remove_slow_arcs ( dest_vesseltree.vtree, threshold, 0, 1 );
  }

  void plot_tree_speeds ( char* ofn, const real_type v_min = ( 1.0 / 1024.0 ), const real_type v_max = 1.0 ) {
    aol::colorTrans colTrans = aol::HSV_BLUE_TO_RED;
    aol::EpsWriter::writeScalePPM ( colTrans, "output/speedscale.ppm" );

    //    EpsWriter* saver = new EpsWriter( ofn, colTrans );
    aol::EpsWriter saver ( ofn, colTrans );

    atree.plot_tree_speeds ( saver, v_min, v_max );
    vtree.plot_tree_speeds ( saver, v_min, v_max );

    //    saver->close();
  }


  geomtree<real_type> atree; // arteries
  geomtree<real_type> vtree; // veins

private:
  vesseltree ( vesseltree &othertree ); // copy constructor

};

/* Other general stuff
   ===================
*/

// TODO: redesign classes so that this need not be inline. currently,
// inlining is necessary beause otherwise these functions are multiply
// defined.


void normalize_speed ( geomtree<double> &tree, double maxi = 1.0 );

// save ppm image: values from ctr-scl to ctr+scl
void img_scaled_save_ctrscl ( const qc::ScalarArray<double, qc::QC_2D> &img, char* filename, char* filename_color, const double ctr, double scl, aol::colorTrans colTrans );

void img_scaled_save ( const qc::ScalarArray<double, qc::QC_2D> &img, char* filename, char* filename_color, aol::colorTrans colTrans );

void img_scaled_save_scl ( const qc::ScalarArray<double, qc::QC_2D> &img, char* filename, char* filename_color, const double scl, aol::colorTrans colTrans );

void img_scaled_save_ctr ( const qc::ScalarArray<double, qc::QC_2D> &img, char* filename, char* filename_color, const double ctr, aol::colorTrans colTrans );

void project_ ( aol::Vector<double> &Vec );


namespace aol {
template < typename RealType, typename VecType = aol::Vector<RealType>, typename OpType = Op<VecType> >
class LagrangeZeroProjectSolver : public aol::InverseOp<VecType> {
public:
  LagrangeZeroProjectSolver ( const OpType &_Op, int _n, RealType _epsilon = 1e-16, int _maxIter = 500, bool _solverquiet = true ) :
      n ( _n ), epsilon ( _epsilon ), maxIter ( _maxIter ), solverquiet ( _solverquiet ), oop ( &_Op ) {
    onerow = new  aol::FullMatrix<RealType> ( 1, n );
    onecol = new  aol::FullMatrix<RealType> ( n, 1 );
    zero11 = new  aol::FullMatrix<RealType> ( 1, 1 );
    AE = new aol::BlockOp<RealType> ( 2, 2 );

    // ATTENTION: should this be just all-1 vector or mass matrix times such vector??
    onerow->setZero();
    onecol->setZero();
    for ( int i = 0; i < n; i++ ) {
      onerow->set ( 0, i, 1.0 );
    }
    for ( int i = 0; i < n; i++ ) {
      onecol->set ( i, 0, 1.0 );
    }
    zero11->setZero();
    AE->setReference ( 0, 0, *oop );
    AE->setReference ( 1, 0, *onerow );
    AE->setReference ( 0, 1, *onecol );
    AE->setReference ( 1, 1, *zero11 );

    AEin = new aol::CGInverse< aol::MultiVector <RealType> > ( *AE, epsilon, maxIter );
    AEin->setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    AEin->setQuietMode ( solverquiet );

  }

  ~LagrangeZeroProjectSolver() {
    delete ( onerow );
    delete ( onecol );
    delete ( zero11 );
    delete (   AE   );
    delete (  AEin  );
  }

  virtual void apply ( const VecType &Arg, VecType &Dest ) const {
    Dest.setZero();
    applyAdd ( Arg, Dest );
  }

  virtual void applyAdd ( const VecType &Arg, VecType &Dest ) const {
    // this is badly implemented. use pointers as member variables!
    // aol::GMRESInverse< aol::MultiVector <RealType> > AEin(AE, epsilon, 50, maxIter); // untested
    // aol::BiCGInverse< aol::MultiVector <RealType> > AEin(AE, AE, epsilon, maxIter);

    aol::MultiVector<RealType> DestE ( 1, Dest.size() ), ArgE ( 1, Arg.size() );
    aol::Vector<RealType> ArgA ( 1 ), DestA ( 1 );

    ArgE[0] = Arg;
    DestE[0] = Dest;

    ArgA.setZero(); DestA.setZero();
    ArgE.appendReference ( ArgA );
    DestE.appendReference ( DestA );

    AEin->apply ( ArgE, DestE ); // sic.
    Dest += DestE[0];
  }

protected:
  int           n;
  RealType      epsilon;
  int           maxIter;
  bool          solverquiet;
  OpType const  *oop;

  aol::FullMatrix<RealType>*  onerow;
  aol::FullMatrix<RealType>*  onecol;
  aol::FullMatrix<RealType>*  zero11;
  aol::BlockOp<RealType>*     AE;
  aol::CGInverse< aol::MultiVector <RealType> > * AEin;

private:
  LagrangeZeroProjectSolver(); // do not allow use of standard constructor

};
} // end namespace aol

#endif

