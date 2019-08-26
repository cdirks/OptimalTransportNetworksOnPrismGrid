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

#ifndef __MEMORYMANAGER_H
#define __MEMORYMANAGER_H

#include <aol.h>

namespace aol {

/** Class for allocating and deallocating memory (16-byte aligned in
 *  case SSE is used).  A certain amount of memory is kept here rather
 *  than returned to the system to improve performance (unless
 *  explicitly disabled by definig DO_NOT_USE_MEMORYMANAGER).  Newest
 *  blocks of memory deallocated are recycled first, oldest blocks are
 *  dropped first.
 *
 *  \todo consistent use of size_t and int data types
 *
 *  \author Schwen (MEVIS), based on older code
 */
class MemoryManager {
public:
  struct MMStore {
    int size;
    void* pBlock;
    MMStore ( const int Size, void* PBlock ) : size ( Size ), pBlock ( PBlock ) { };
  };

private:
  static list< MMStore > _storage;

  static int
    _numStored,       //!< number of memory blocks being stored
    _memusage,        //!< memory being used
    _maxRetain,       //!< Maximum number of memory blocks to retain in manager
    _memusageLimit;   //!< maximum memory usage

public:
  //! Return how much memory is used by the MemoryManager
  static int memoryManagerMemoryUsage ();

  //! Return a pointer to memory of at least size Length * PointeeSize, write actual length to Length
  //! \warning The user has to store the actual length returned and pass this value to deallocate later, otherwise memory will get lost.
  static void* allocateAtLeast ( int& Length, const size_t PointeeSize );

  //! Vector manager, store for reuse
  static void deallocate ( void *Ptr, const int Length, const size_t PointeeSize );

  //! Clean reuse storage
  //! \todo rename method
  static void deleteUnlocked ( const int NumToRetain = 0 );

  //! set maximum number of vectors to be kept, does not affect current size
  static void setMaxRetain ( const int MaxRetain );

  //! set maximum amount of memory to be kept, does not affect current size
  static void setMemusageLimit ( const int MaxMemusage );

  //! to free used memory on exit to prevent memory leak, should not be called directly except by atexit
  static void clearOnExit ( ) {
    deleteUnlocked();
  }
};

}

#endif
