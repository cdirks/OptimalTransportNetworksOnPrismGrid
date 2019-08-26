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

#ifndef __GRAPE2_H
#define __GRAPE2_H

#ifdef USE_EXTERNAL_GRAPE

#ifdef __GNUC__
#pragma GCC system_header
#endif

#define G_NO_CONST 1

#define class GRAPE_class
#define CLASS GRAPE_CLASS
#define private GRAPE_private
#define explicit GRAPE_explicit
#define EXPLICIT GRAPE_EXPLICIT
#define friend GRAPE_friend
#define GRAPE_METHOD_RETURN ___GRAPE_METHOD_RETURN
#define GRAPE_METHOD_CODE ___GRAPE_METHOD_CODE
#define G_CPP
extern "C" {

  #include <grape.h>

}
#undef GRAPE_METHOD_RETURN
#undef GRAPE_METHOD_CODE
#undef class
#undef CLASS
#undef private
#undef explicit
#undef EXPLICIT
#undef friend

#undef GRAPE
#ifdef G_CODE_CHECK
#define GRAPE(obj,meth) (*grape((obj),const_cast<char*>(meth),G_HEADER_CODE,#meth ZEROSTR2[0]))
#else
#define GRAPE(obj,meth) (*grape((obj),const_cast<char*>(meth),#meth ZEROSTR2[0]))
#endif

#else /* USE_EXTERNAL_GRAPE */

#error GRAPE external is needed to compile module grape

#endif /* USE_EXTERNAL_GRAPE */

#endif
