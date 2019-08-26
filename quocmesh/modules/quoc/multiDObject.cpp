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

#include <multiDObject.h>

namespace qc {

template<>
void qc::MultiDStorageObject<qc::QC_3D>::reallocate ( const GridStructure &Grid ) {
  QUOC_ASSERT ( Grid.getDimOfWorld() == QC_3D );
  this->reallocate ( Grid.getNumX(), Grid.getNumY(), Grid.getNumZ() );
}


template<>
void qc::MultiDStorageObject<qc::QC_3D>::resize ( const GridStructure &Grid ) {
  QUOC_ASSERT ( Grid.getDimOfWorld() == QC_3D );
  resize ( Grid.getNumX(), Grid.getNumY(), Grid.getNumZ() );
}

}
