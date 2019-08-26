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

// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: hs071_main.cpp 1428 2009-04-18 00:41:05Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10

#include "ipoptExampleNLP.hpp"

using namespace Ipopt;

int main( void )
{
  // Create an instance of your nlp...
  // SmartPtr is a pointer implementation of Ipopt, comparable to STL's auto_ptr
  SmartPtr<TNLP> mynlp = new HS071_NLP();

  // Create an instance of the IpoptApplication
  SmartPtr<IpoptApplication> app = new IpoptApplication();

  // Set some options, see documentation, uncomment following line to see a list of all options
//   app->Options()->SetStringValue( "print_options_documentation", "yes" );
  app->Options()->SetStringValue(  "derivative_test",              "second-order" ); // or: first-order
  app->Options()->SetStringValue(  "derivative_test_print_all",    "no" );
  app->Options()->SetNumericValue( "derivative_test_tol",          1e-6 );
//   app->Options()->SetStringValue( "hessian_approximation", "limited-memory" );
//   app->Options()->SetNumericValue( "tol", 1e-4 );
//   app->Options()->SetNumericValue( "dual_inf_tol", 2.0 );
//   app->Options()->SetStringValue( "linear_solver",  ma57 );

  // Initialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    std::cerr << "\n\n*** Error during initialization!\n";
    return static_cast<int>( status );
  }

  // Optimize!
  status = app->OptimizeTNLP(mynlp);

  if (status == Solve_Succeeded) {
    // Retrieve some statistics about the solve
    Index iter_count = app->Statistics()->IterationCount();
    std::cout << "\n\n*** The problem solved in " << iter_count << " iterations!";

    Number final_obj = app->Statistics()->FinalObjective();
    std::cout << "\n\n*** The final value of the objective function is " << final_obj << ".\n\n";
  }

  return static_cast<int>( status );
}
