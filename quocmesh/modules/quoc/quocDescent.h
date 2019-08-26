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

#ifndef __QUOCDESCENT_H
#define __QUOCDESCENT_H

#include <quoc.h>
#include <gridBase.h>

namespace qc {

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class MultilevelDescentInterface {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  typedef typename ConfiguratorType::InitType InitType;
  const int _maxDepth;
  const InitType _grid;
  const InitType *_curGrid;
  int _curLevel;
  bool _logLevelStart;
  vector<const InitType*> _grids;
public:
  MultilevelDescentInterface ( const int MaxDepth )
   : _maxDepth ( MaxDepth ),
     _grid ( MaxDepth, ConfiguratorType::Dim ),
     _curGrid ( NULL ),
     _curLevel ( MaxDepth ),
     _logLevelStart ( false ) {

    for ( int level = 0; level <= _grid.getGridDepth(); level++ ) {
      _grids.push_back( new InitType( level, ConfiguratorType::Dim ) );
    }

    _curGrid = _grids[ _curLevel ];
  }

  virtual ~MultilevelDescentInterface( ) {
    for ( int level = 0; level <= _grid.getGridDepth(); level++ ) {
     delete _grids[ level ];
   }
  }

  RealType H() const {
    return _curGrid->H();
  }

  int getLevel( ) const {
    return _curLevel;
  }

  int getMaxGridDepth( ) const {
    return _grid.getGridDepth();
  }

  const InitType& getInitializerRef( ) const {
    return _grid;
  }

  const InitType& getCurrentGrid () const {
    return *this->_curGrid;
  }

  virtual void setLevel( const int Level ) {
    QUOC_ASSERT ( (Level <= _grid.getGridDepth()) && (Level >= 0) );
    _curLevel = Level;
    _curGrid = _grids[ _curLevel ];
  }

  void setLogLevelStart ( const bool LogLevelStart ) {
    _logLevelStart = LogLevelStart;
  }

  virtual void prolongate( ) = 0;

  virtual void descentOnCurrentGrid( ) = 0;

  void solve( const int StartLevel, const int StopLevel, const int LevelIncrement = 1 ) {
    setLevel( StartLevel );
    for ( int level = StartLevel; level <= StopLevel; ) {
      if ( _logLevelStart ) {
        cerr << "\n--------------------------------------------------------------------------------\n";
        cerr << "Descent on level " << getLevel() << " started";
        cerr << "\n";
        cerr << "--------------------------------------------------------------------------------\n\n";
      }

      descentOnCurrentGrid();
      for ( int i = 0; i < LevelIncrement; i++ ) {
        if ( level < StopLevel )
          prolongate( );
        ++level;
      }
    }
  }
};

} // end namespace qc

#endif // __QUOCDESCENT_H
