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

#ifndef __MPIINCLUDES_H
#define __MPIINCLUDES_H

#ifdef __GNUC__
#pragma GCC system_header
#endif

// The mpi headers use old style casts in defines, so flagging them as system headers is not enough.
WARNING_OFF(old-style-cast)
#include <mpi.h>

#endif // __MPIINCLUDES_H
