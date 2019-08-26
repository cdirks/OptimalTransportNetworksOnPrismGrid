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

#ifndef __GHOSTEXCHANGE_H
#define __GHOSTEXCHANGE_H

namespace qc {

enum Face {GH_LEFT, GH_RIGHT, GH_TOP, GH_BUTTOM, GH_FRONT, GH_BACK};

/**
 * \param buffer must be of size num_bins_local[] * num_bins_local[]
 */
void pack_buffer ( const double*bins_local, const int num_bins_local[3], double* buffer, Face face );

void unpack_buffer ( double*bins_local, const int num_bins_local[3], const double* buffer, Face face );

void exchange_faces ( double* bins_local, int num_bins_local[3], int count, int nbrL, int nbrR, Face faceL, Face faceR );

void ghost_exchange ( double*bins_local, int num_bins_local[3], const int ( &procneigh ) [3][2] );

}

#endif // __GHOSTEXCHANGE_H
