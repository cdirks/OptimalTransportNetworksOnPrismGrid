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

#ifndef __SEGDEFS_H
#define __SEGDEFS_H

#include <aol.h>
#include <ocTree.h>

// definitions used for segmentation in Eikonal...
// and in qc::OcTree

namespace eik {

const int TN_HN_NONE_SEG         = -1;

const int TN_HN_EDGE_FRONT_UP    = 0;
const int TN_HN_EDGE_FRONT_LOW   = 1;
const int TN_HN_EDGE_FRONT_LEFT  = 2;
const int TN_HN_EDGE_FRONT_RIGHT = 3;

const int TN_HN_EDGE_BACK_UP     = 4;
const int TN_HN_EDGE_BACK_LOW    = 5;
const int TN_HN_EDGE_BACK_LEFT   = 6;
const int TN_HN_EDGE_BACK_RIGHT  = 7;

const int TN_HN_EDGE_RIGHT_UP    = 8;
const int TN_HN_EDGE_RIGHT_LOW   = 9;
const int TN_HN_EDGE_LEFT_UP     = 10;
const int TN_HN_EDGE_LEFT_LOW    = 11;

const int TN_HN_MAX_EDGE_INDEX   = 11;

const int TN_HN_FACE_UP          = 12;
const int TN_HN_FACE_LOW         = 13;
const int TN_HN_FACE_LEFT        = 14;
const int TN_HN_FACE_RIGHT       = 15;
const int TN_HN_FACE_FRONT       = 16;
const int TN_HN_FACE_BACK        = 17;

const int TN_HN_MAX_NUM          = 17;

const int TN_DIR_X     = qc::TN_DIR_RIGHT;
const int TN_DIR_Y     = qc::TN_DIR_UP;
const int TN_DIR_Z     = qc::TN_DIR_BACK;

const int TN_DIR_RIGHT = TN_DIR_X;
const int TN_DIR_UP    = TN_DIR_Y;
const int TN_DIR_BACK  = TN_DIR_Z;
const int TN_DIR_ALL   = TN_DIR_UP | TN_DIR_BACK | TN_DIR_RIGHT;
const int TN_DIR_OPP   = 8;   // opposite direction, e.g. TN_DIR_RIGHT | TN_DIR_OPP means left

const int TN_DIR_FRONT_LOW_LEFT   = 0;
const int TN_DIR_FRONT_LOW_RIGHT  = TN_DIR_RIGHT;
const int TN_DIR_FRONT_UP_LEFT    = TN_DIR_UP;
const int TN_DIR_FRONT_UP_RIGHT   = TN_DIR_UP | TN_DIR_RIGHT;
const int TN_DIR_BACK_LOW_LEFT    = TN_DIR_BACK;
const int TN_DIR_BACK_LOW_RIGHT   = TN_DIR_BACK | TN_DIR_RIGHT;
const int TN_DIR_BACK_UP_LEFT     = TN_DIR_UP | TN_DIR_BACK;
const int TN_DIR_BACK_UP_RIGHT    = TN_DIR_UP | TN_DIR_BACK | TN_DIR_RIGHT;
const int TN_DIR_MAX_NUM          = TN_DIR_UP | TN_DIR_BACK | TN_DIR_RIGHT;

typedef ostream outstream;
typedef istream instream;
typedef ofstream outfilestream;

}

#endif
