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

// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: hs071_nlp.cpp 1324 2008-09-16 14:19:26Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16
//
// Slightly modified to adapt to QuocMesh

#include <aol.h>

#include "ipoptExampleNLP.hpp"

using namespace Ipopt;

// constructor
HS071_NLP::HS071_NLP()
{}

//destructor
HS071_NLP::~HS071_NLP()
{}

// returns the size of the problem
bool HS071_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in HS071_NLP.hpp has 4 variables, x[0] through x[3]
  n = 4;

  // one equality constraint and one inequality constraint
  m = 2;

  // in this example the jacobian is dense and contains 8 nonzeros
  nnz_jac_g = 8;

  // the hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = 10;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool HS071_NLP::get_bounds_info(Index HIDE_PARAMETER_IN_OPT_MODE(n), Number* x_l, Number* x_u,
                                Index HIDE_PARAMETER_IN_OPT_MODE(m), Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  QUOC_ASSERT(n == 4);
  QUOC_ASSERT(m == 2);

  // the variables have lower bounds of 1
  for (Index i=0; i<4; i++) {
    x_l[i] = 1.0;
  }

  // the variables have upper bounds of 5
  for (Index i=0; i<4; i++) {
    x_u[i] = 5.0;
  }

  // the first constraint g1 has a lower bound of 25
  g_l[0] = 25;
  // the first constraint g1 has NO upper bound, here we set it to 2e19.
  // Ipopt interprets any number greater than nlp_upper_bound_inf as
  // infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
  // is 1e19 and can be changed through ipopt options.
  g_u[0] = 2e19;

  // the second constraint g2 is an equality constraint, so we set the
  // upper and lower bound to the same value
  g_l[1] = g_u[1] = 40.0;

  return true;
}

// returns the initial point for the problem
bool HS071_NLP::get_starting_point(Index /*n*/, bool HIDE_PARAMETER_IN_OPT_MODE(init_x), Number* x,
                                   bool HIDE_PARAMETER_IN_OPT_MODE(init_z), Number* /*z_L*/, Number* /*z_U*/,
                                   Index /*m*/, bool HIDE_PARAMETER_IN_OPT_MODE(init_lambda),
                                   Number* /*lambda*/)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  QUOC_ASSERT(init_x == true);
  QUOC_ASSERT(init_z == false);
  QUOC_ASSERT(init_lambda == false);

  // initialize to the given starting point
  x[0] = 1.0;
  x[1] = 5.0;
  x[2] = 5.0;
  x[3] = 1.0;

  return true;
}

// returns the value of the objective function
bool HS071_NLP::eval_f(Index HIDE_PARAMETER_IN_OPT_MODE(n), const Number* x, bool /*new_x*/, Number& obj_value)
{
  QUOC_ASSERT(n == 4);

  obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool HS071_NLP::eval_grad_f(Index HIDE_PARAMETER_IN_OPT_MODE(n), const Number* x, bool /*new_x*/, Number* grad_f)
{
  QUOC_ASSERT(n == 4);

  grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
  grad_f[1] = x[0] * x[3];
  grad_f[2] = x[0] * x[3] + 1;
  grad_f[3] = x[0] * (x[0] + x[1] + x[2]);

  return true;
}

// return the value of the constraints: g(x)
bool HS071_NLP::eval_g(Index HIDE_PARAMETER_IN_OPT_MODE(n), const Number* x, bool /*new_x*/, Index HIDE_PARAMETER_IN_OPT_MODE(m), Number* g)
{
  QUOC_ASSERT(n == 4);
  QUOC_ASSERT(m == 2);

  g[0] = x[0] * x[1] * x[2] * x[3];
  g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];

  return true;
}

// return the structure or values of the jacobian
bool HS071_NLP::eval_jac_g(Index /*n*/, const Number* x, bool /*new_x*/,
                           Index /*m*/, Index /*nele_jac*/, Index* iRow, Index *jCol,
                           Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    // this particular jacobian is dense
    iRow[0] = 0;
    jCol[0] = 0;
    iRow[1] = 0;
    jCol[1] = 1;
    iRow[2] = 0;
    jCol[2] = 2;
    iRow[3] = 0;
    jCol[3] = 3;
    iRow[4] = 1;
    jCol[4] = 0;
    iRow[5] = 1;
    jCol[5] = 1;
    iRow[6] = 1;
    jCol[6] = 2;
    iRow[7] = 1;
    jCol[7] = 3;
  }
  else {
    // return the values of the jacobian of the constraints

    values[0] = x[1]*x[2]*x[3]; // 0,0
    values[1] = x[0]*x[2]*x[3]; // 0,1
    values[2] = x[0]*x[1]*x[3]; // 0,2
    values[3] = x[0]*x[1]*x[2]; // 0,3

    values[4] = 2*x[0]; // 1,0
    values[5] = 2*x[1]; // 1,1
    values[6] = 2*x[2]; // 1,2
    values[7] = 2*x[3]; // 1,3
  }

  return true;
}

//return the structure or values of the hessian
bool HS071_NLP::eval_h(Index /*n*/, const Number* x, bool /*new_x*/,
                       Number obj_factor, Index /*m*/, const Number* lambda,
                       bool /*new_lambda*/, Index HIDE_PARAMETER_IN_OPT_MODE(nele_hess), Index* iRow,
                       Index* jCol, Number* values)
{
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense
    Index idx=0;
    for (Index row = 0; row < 4; row++) {
      for (Index col = 0; col <= row; col++) {
        iRow[idx] = row;
        jCol[idx] = col;
        idx++;
      }
    }

    QUOC_ASSERT(idx == nele_hess);
  }
  else {
    // return the values. This is a symmetric matrix, fill the lower left
    // triangle only

    // fill the objective portion
    values[0] = obj_factor * (2*x[3]); // 0,0

    values[1] = obj_factor * (x[3]);   // 1,0
    values[2] = 0.;                    // 1,1

    values[3] = obj_factor * (x[3]);   // 2,0
    values[4] = 0.;                    // 2,1
    values[5] = 0.;                    // 2,2

    values[6] = obj_factor * (2*x[0] + x[1] + x[2]); // 3,0
    values[7] = obj_factor * (x[0]);                 // 3,1
    values[8] = obj_factor * (x[0]);                 // 3,2
    values[9] = 0.;                                  // 3,3


    // add the portion for the first constraint
    values[1] += lambda[0] * (x[2] * x[3]); // 1,0

    values[3] += lambda[0] * (x[1] * x[3]); // 2,0
    values[4] += lambda[0] * (x[0] * x[3]); // 2,1

    values[6] += lambda[0] * (x[1] * x[2]); // 3,0
    values[7] += lambda[0] * (x[0] * x[2]); // 3,1
    values[8] += lambda[0] * (x[0] * x[1]); // 3,2

    // add the portion for the second constraint
    values[0] += lambda[1] * 2; // 0,0

    values[2] += lambda[1] * 2; // 1,1

    values[5] += lambda[1] * 2; // 2,2

    values[9] += lambda[1] * 2; // 3,3
  }

  return true;
}

bool HS071_NLP::intermediate_callback( AlgorithmMode /*mode*/, Index /*iter*/, Number /*obj_value*/, Number /*inf_pr*/, Number /*inf_du*/,
                                      Number /*mu*/, Number /*d_norm*/, Number /*regularization_size*/, Number /*alpha_du*/, Number /*alpha_pr*/,
                                      Index /*ls_trials*/, const IpoptData* /*ip_data*/, IpoptCalculatedQuantities* /*ip_cq*/ )
{
  // print intermediate results etc.
  // return false to request termination in case something went terribly wrong

  return true;
}

void HS071_NLP::finalize_solution(SolverReturn /*status*/,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* /*lambda*/,
                                  Number obj_value,
                                  const IpoptData* /*ip_data*/,
                                  IpoptCalculatedQuantities* /*ip_cq*/)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  std::cout << "\n\nSolution of the primal variables, x\n";
  for (Index i=0; i<n; i++) {
    std::cout << "x[" << i << "] = " << x[i] << "\n";
  }

  std::cout << "\n\nSolution of the bound multipliers, z_L and z_U\n";
  for (Index i=0; i<n; i++) {
    std::cout << "z_L[" << i << "] = " << z_L[i] << "\n";
  }
  for (Index i=0; i<n; i++) {
    std::cout << "z_U[" << i << "] = " << z_U[i] << "\n";
  }

  std::cout << "\n\nObjective value\n";
  std::cout << "f(x*) = " << obj_value << "\n";

  std::cout << "\nFinal value of the constraints:\n";
  for (Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << "\n";
  }
}
