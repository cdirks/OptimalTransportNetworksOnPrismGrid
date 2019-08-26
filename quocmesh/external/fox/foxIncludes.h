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

#ifndef __FOXINCLUDES_H
#define __FOXINCLUDES_H

#ifdef __GNUC__
// This warning has to be turned off not only inside the fox headers, but also inside
// code that uses fox. Otherwise the defines like FXDECLARE that need to be used to
// create new classes derived from fox classes cause warnings.
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC system_header
#endif

// [BB] complex needs to be included here, otherwise including of complex in aol.h
// fails under VC++.
#include <complex>
#include <fx.h>
#include <fx3d.h>

#ifdef PI
#undef PI
#endif
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

#ifdef __GNUC__
/* #pragma GCC diagnostic warning "-Wold-style-cast" */
#endif

// small difference between 1.6 (default version) and 1.7
#if ( ( FOX_MAJOR == 1 ) && ( FOX_MINOR == 7 ) )
#define VISUAL_DOUBLEBUFFER VISUAL_DOUBLE_BUFFER 
#endif

#endif
