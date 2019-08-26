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

#ifndef __QUADTREE_H
#define __QUADTREE_H

#include <bitVector.h>

namespace qc {

class QuadTree;

// this class can only be used by the QuadTree;

class QuadGridField {
  friend class QuadTree;

protected:
  QuadGridField ( int Step, int NumX, int NumY );

  ~QuadGridField();

  void setElement ( int X, int Y );

  void clearElement ( int X, int Y );

  bool checkElement ( int X, int Y );

  void clearAll();

  void setAll() {
    bitField->setAll();
  }

private:

  int step, numX, numY;

  aol::BitVector *bitField;

};

typedef QuadGridField* PQuadGridField;

class QuadTree {
public:
  QuadTree ( int GridDepth );

  ~QuadTree();

  void setElement ( int X, int Y, int Level );

  void setLowLeftElement ( int X, int Y, int Level ) {
    setElement ( X - getStep ( Level ), Y - getStep ( Level ), Level );
  }

  void setLowRightElement ( int X, int Y, int Level ) {
    setElement ( X, Y - getStep ( Level ), Level );
  }

  void setUpLeftElement ( int X, int Y, int Level ) {
    setElement ( X - getStep ( Level ), Y, Level );
  }

  void setUpRightElement ( int X, int Y, int Level ) {
    setElement ( X, Y, Level );
  }


  void setElementAndUnite ( int X, int Y, int Level );

  void setElementAndClearChildren ( int X, int Y, int Level );

  bool checkForUnite ( int X, int Y, int Level );

  void clearElement ( int X, int Y );

  void clearElement ( int X, int Y, int Level );

  bool getElement ( int X, int Y, int Level );

  bool upNeighbour ( int X, int Y, int Level );

  bool lowNeighbour ( int X, int Y, int Level );

  bool leftNeighbour ( int X, int Y, int Level );

  bool rightNeighbour ( int X, int Y, int Level );

  void   uniteElements ( int X, int Y, int Level );

  void clearAll();

  int getStep ( int Level ) {
    if ( Level > maxLevel ) {
      cerr << "ERROR in QuadTree::getStep! Level = " << Level << " out of range..\n";
    }
    return ( 1 << ( maxLevel - Level ) );
  }

  // this function checks if the node is member of the grid at a certain level
  int nodeInGrid ( int X, int Y, int Level ) {
    if ( Level == maxLevel ) return 1;
    if ( Level < 0 || Level > maxLevel ) return 0;

    int mask = 1;
    for ( int level = maxLevel; level > Level; level-- ) {
      if ( X & mask || Y & mask ) return 0;
      mask <<= 1;
    }
    return 1;
  }

  // returns the coarsest(lowest) number of level in which the node (X,Y) appears
  int getCoarsestLevel ( int X, int Y );

  int getMaxLevel() {
    return maxLevel;
  }

  int getWidth() {
    return width;
  }

  int getHeight() {
    return height;
  }

  void statistics ( void );

  void setAllFine() {
    gridFields[maxLevel]->setAll();
  }

protected:

  QuadGridField **gridFields;

  int width, height, maxLevel;
};

}

#endif
