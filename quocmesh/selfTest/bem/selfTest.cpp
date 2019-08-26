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

#include "ahmed.h"
#include "bemesh.h"
#include "bemOps.h"
#include "blocksolver.h"
#include "boundary.h"
#include "curvatureOps.h"
#include "elasticFunSol.h"
#include "elasticOps.h"
#include "elasticTensor.h"
#include "laplaceFunSol.h"
#include "parametric.h"
#include "particle.h"
#include "rectangle.h"
#include "segment.h"
#include "typedParser.h"

int main ()
{
  aol::printSelfTestSuccessMessage ( "--                        BEM Self Test Successful                            --" );
  // aol::printSelfTestFailureMessage ( "!!                        BEM Self Test FAILED                                !!" );

  return 0;
}
