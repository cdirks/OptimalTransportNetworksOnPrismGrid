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

#ifndef __CIMGINCLUDES_H
#define __CIMGINCLUDES_H

#include <platformDependent.h>

#ifdef __GNUC__
WARNING_OFF(old-style-cast)
#pragma GCC system_header
#endif

#ifdef USE_LIB_PNG
// Enable LibPNG support
#define cimg_use_png 1
#endif
#include <CImg.h>

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
/* WARNING_ON(old-style-cast) */
#endif

#endif // __CIMGINCLUDES_H
