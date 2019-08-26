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

@file

@brief The example shows how to solve a mixed boundary value problem in linearized elasticity using the
2D collocation boundary element implementation in module bem.

The example considers a square with dirichlet boundary conditions on top and bottom specifying compression, and natural boundary
conditions on the sides. The Solution uses the class BlockSolver that contains the discretization of the direct
boundary integral equation and constraints specifying the boundary conditions. The unknown vectors has length @f$4n@f$ and has the form
@f{eqnarray*}
  \Big( u_1 (\xi_1), u_1 (\xi_2), \ldots u_1 (\xi_n), & & u_2 (\xi_1), u_2 (\xi_2), \ldots u_2 (\xi_n), \\
        t_1 (\xi_1), t_1 (\xi_2), \ldots t_1 (\xi_n), & & t_2 (\xi_1), t_2 (\xi_2), \ldots t_2 (\xi_n) \Big)\,,
@f}
where @f$u@f$ is the displacement and @f$t = C \epsilon \nu@f$ the traction on the boundary and @f$\xi_i@f$ are the collocation nodes.

Finally, the deformed square is plotted as a postscript image (using gnuplot).

@author Lenz

*/

#include "elasticOps.h"
#include "boundary.h"
#include "parametric.h"
#include "blocksolver.h"

int main ( int /* argc */, char* /* argv */ [] ) {
  try {

    int side = 30, n = 4 * side;                                                                              // points per side of square, total points
    bm::Boundary<bm::ParaParticle<double> > boundary;                                                         // boundary container
    bm::ParaParticle<double> square ( aol::Vec2<double> ( 0.5, 0.5 ), aol::Vec2<double> ( 0.4, 0.4 ), side ); // square boundary component, starts counter-clockwise in bottom left corner
    boundary.push_back ( square );                                                                            // insert boundary component into container

    aol::ElasticTensor<double> elast ( 3.0, 1.0, 1.1 );                                                        // parameters of elasticity tensor with cubic symmetry (c11, c12, c44 in Voigt notation)
    bm::ElasticGreen<double> green ( elast );                                                                 // compute coefficients of fundamental solution

    aol::BlockSolver<double> system;                                                                          // solver for special type of block systems
    aol::Vector<double> sol ( 4*n );                                                                          // contains four sub-vectors ( x-component of boundary stress,
    int tracx = 0, tracy = n, dispx = 2*n, dispy = 3*n;                                                       //                             y-component of boundary stress,
                                                                                                              //                             x-component of displacement,
                                                                                                              //                             y-component of displacement),
                                                                                                              // each then indexed by element numbers

    bm::ElasticOperator<bm::ParaParticle<double>,bm::ElasticSolution> singop ( boundary, green, true );       // single layer operator
    bm::ElasticOperator<bm::ParaParticle<double>,bm::ElasticSolutionDiff> doubop ( boundary, green, true );   // double layer operator
    aol::Vector<double> rhs ( 2*n ); rhs.setZero ();                                                          // right hand side is zero
    system.appendBlockRef ( singop, doubop, rhs );                                                            // direct boundary integral equation

    system.appendConstBlock ( 0, dispx, dispx + side, 0 );                                                    // (0,0) Dirichlet values on bottom
    system.appendConstBlock ( 0, dispy, dispy + side, 0 );
    system.appendConstBlock ( 0, dispx + 2 * side, dispx + 3 * side, 0 );                                     // (0,1) Dirichlet values on top
    system.appendConstBlock ( 0, dispy + 2 * side, dispy + 3 * side, 0.2 );
    system.appendConstBlock ( 0, tracx + side, tracx + 2 * side, 0 );                                         // zero Neumann values on right side
    system.appendConstBlock ( 0, tracy + side, tracy + 2 * side, 0 );
    system.appendConstBlock ( 0, tracx + 3 * side, tracx + 4 * side, 0 );                                     // zero Neumann values on left side
    system.appendConstBlock ( 0, tracy + 3 * side, tracy + 4 * side, 0 );

    system.apply_fast ( system.getRHS_fast (), sol );                                                         // solve equation system with boundary values given above

/*                                       // text output
    aol::Vector<double> tmp ( n );                                                                            // temporary values
    sol.getBlock ( dispx, tmp ); cerr << "X displacement: " << tmp << endl;                                   // print x displacements
    sol.getBlock ( dispy, tmp ); cerr << "Y displacement: " << tmp << endl;                                   // print y displacements
    sol.getBlock ( tracx, tmp ); cerr << "X traction:     " << tmp << endl;                                   // print x tractions
    sol.getBlock ( tracy, tmp ); cerr << "Y traction:     " << tmp << endl;                                   // print y tractions
*/

    boundary.setPlotPoints(true);                        // also mark discretisation points
    boundary.setFormat( bm::Boundary<bm::ParaParticle<double> >::POSTSCRIPT );              // choose postscript as output format
    {
  aol::opfstream os ( "| gnuplot > org.eps" );                       // plot original square
        os << boundary;
    }

    aol::MultiVector<double> v1 ( 2, n );                      // get displacement
    sol.getBlock ( dispx, v1[0] );
    sol.getBlock ( dispy, v1[1] );
    boundary.move( v1, true );                          // move square

    aol::opfstream os ( "| gnuplot > out.eps" );                           // plot displaced square
    os << boundary;

  } catch ( aol::Exception& ) {
    return 23;
  }

  return 0;
}
