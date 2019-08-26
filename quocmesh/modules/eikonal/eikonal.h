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

#ifndef __EIKONAL_H
#define __EIKONAL_H

#include <quoc.h>
#include <quadTree.h>
#include <qmHeap.h>
#include <estimator2d.h>

namespace eik {

const int TN_HN_NONE  = 0;
const int TN_HN_UP    = 1;
const int TN_HN_LOW   = 2;
const int TN_HN_LEFT  = 3;
const int TN_HN_RIGHT = 4;
const int TN_HN_INTERPOLATED = 8;

struct MetaEstimatorStruct {
  int grid_depth;

  int checkElement ( int /*X*/, int /*Y*/, int Level ) {
    return ( Level == grid_depth );
  }

  int checkLowLeftElement ( int X, int Y, int Level ) {
    int step = 1 << ( grid_depth - Level );
    return checkElement ( X - step, Y - step, Level );
  }

  int checkLowRightElement ( int X, int Y, int Level ) {
    int step = 1 << ( grid_depth - Level );
    return checkElement ( X, Y - step, Level );
  }

  int checkUpLeftElement ( int X, int Y, int Level ) {
    int step = 1 << ( grid_depth - Level );
    return checkElement ( X - step, Y, Level );
  }

  int checkUpRightElement ( int X, int Y, int Level ) {
    return checkElement ( X, Y, Level );
  }
};

template <typename RealType = float>
class Eikonal {
public:

class TrialNode : public aol::Vec3<short> {
  public:
    TrialNode ( aol::Vec3<short> Coord )
        : aol::Vec3<short> ( Coord ),
        upLeftLevel ( -1 ),
        upRightLevel ( -1 ),
        lowLeftLevel ( -1 ),
        lowRightLevel ( -1 ),
        hangNodeLevel ( -1 ),
        hangNode ( TN_HN_NONE ),
        val ( -1.0 )  {}

    TrialNode ( short X, short Y )
        : aol::Vec3<short> ( X, Y, 0 ),
        upLeftLevel ( -1 ),
        upRightLevel ( -1 ),
        lowLeftLevel ( -1 ),
        lowRightLevel ( -1 ),
        hangNodeLevel ( -1 ),
        hangNode ( TN_HN_NONE ),
        val ( -1.0 )  {}

    TrialNode( )
        : upLeftLevel ( -1 ),
        upRightLevel ( -1 ),
        lowLeftLevel ( -1 ),
        lowRightLevel ( -1 ),
        hangNodeLevel ( -1 ),
        hangNode ( TN_HN_NONE ),
        val ( -1.0 ) {}

    bool operator== ( const TrialNode &Other ) const {
      return ( val == Other.val );
    }

    bool operator< ( const TrialNode &Other ) const {
      return ( val < Other.val );
    }

    bool operator> ( const TrialNode &Other ) const {
      return ( val > Other.val );
    }

    signed char upLeftLevel;
    signed char upRightLevel;
    signed char lowLeftLevel;
    signed char lowRightLevel;
    signed char hangNodeLevel;
    signed char hangNode;
    RealType val;
  };

  Eikonal ( const qc::GridDefinition &Grid )
      : grid ( Grid ) ,
      timeField ( NULL ),
      scaleTime ( 1.0f ),
      order ( 1 ), quadTree ( NULL ), image ( NULL ) {
    elementEstimator.grid_depth = Grid.getGridDepth();
    quadTree = new qc::QuadTree ( Grid.getGridDepth() );
    timeField = new qc::Array<RealType> ( Grid );
    timeField->setAll ( -1. );
  }

  ~Eikonal( ) {
    if ( quadTree != NULL )
      delete quadTree;
    if ( timeField != NULL )
      delete timeField;
  }

  void reset( ) {
    trialNodeHeap.erase ( trialNodeHeap.begin( ), trialNodeHeap.end( ) );
    quadTree->clearAll( );
    timeField->setAll ( -1.0f );
  }

  void setOrder ( int Order ) { order = Order; }

  void setSeedPoint ( const aol::Vec3<short> Point ) {
    short NewX = 0, NewY = 0;

    for ( int Level = 0; Level <= grid.getGridDepth(); Level++ ) {
      int Step = grid.h[ Level ];
      NewX = Point.x() - ( Point.x() % Step );
      NewY = Point.y() - ( Point.y() % Step );
      if ( elementEstimator.checkElement ( NewX, NewY, Level ) ) {
        cerr << "Level = " << Level << endl;
        break;
      }
    }

    TrialNode Node ( NewX, NewY );
    Node.val = 0.0f;

    trialNodeHeap.push ( Node );
  }

  void makeTrialNodeActive ( TrialNode &Node );

  RealType solveQuadratic ( RealType v1, RealType v2, RealType h1, RealType h2, RealType F = 1.0f );

  RealType solveQuadratic ( RealType a, RealType b, RealType c );

  void createElementUpRight ( TrialNode &Node ) {
    if ( Node.upLeftLevel != -1 && Node.lowRightLevel != -1 ) {
      if ( aol::Abs ( Node.upLeftLevel - Node.lowRightLevel ) == 2 ) {
        Node.upRightLevel = ( Node.upLeftLevel + Node.lowRightLevel ) >> 1;
        quadTree->setElement ( Node.x(), Node.y(), Node.upRightLevel );
        updateElementToTrialNodes ( Node.x(), Node.y(), Node.upRightLevel );
        return;
      }
    }

    int RangeUp = 1000, RangeLow = -1;
    if ( Node.upLeftLevel != -1 ) {
      if ( RangeUp > Node.upLeftLevel + 1 )
        RangeUp = Node.upLeftLevel + 1;
      if ( RangeLow < Node.upLeftLevel - 1 )
        RangeLow = Node.upLeftLevel - 1;
    }
    if ( Node.lowRightLevel != -1 ) {
      if ( RangeUp > Node.lowRightLevel + 1 )
        RangeUp = Node.lowRightLevel + 1;
      if ( RangeLow < Node.lowRightLevel - 1 )
        RangeLow = Node.lowRightLevel - 1;
    }

    if ( RangeLow < 0 ) RangeLow = 0;
    if ( RangeUp > grid.getGridDepth() ) RangeUp = grid.getGridDepth();

    for ( int level = RangeLow; level < RangeUp; level++ ) {
      if ( quadTree->nodeInGrid ( Node.x(), Node.y(), level )
           && elementEstimator.checkUpRightElement ( Node.x(), Node.y(), level ) ) {
        quadTree->setElement ( Node.x(), Node.y(), level );
        updateElementToTrialNodes ( Node.x(), Node.y(), level );
        Node.upRightLevel = level;
        return;
      }
    }

    quadTree->setElement ( Node.x(), Node.y(), RangeUp );
    updateElementToTrialNodes ( Node.x(), Node.y(), RangeUp );
    Node.upRightLevel = RangeUp;
  }

  void createElementUpLeft ( TrialNode &Node ) {
    if ( Node.lowLeftLevel != -1 && Node.upRightLevel != -1 ) {
      if ( aol::Abs ( Node.lowLeftLevel - Node.upRightLevel ) == 2 ) {
        Node.upLeftLevel = ( Node.lowLeftLevel + Node.upRightLevel ) >> 1;
        quadTree->setUpLeftElement ( Node.x(), Node.y(), Node.upLeftLevel );
        int step = 1 << ( grid.getGridDepth()  - Node.upLeftLevel );
        updateElementToTrialNodes ( Node.x() - step, Node.y(), Node.upLeftLevel );
        return;
      }
    }

    int RangeUp = 1000, RangeLow = -1;
    if ( Node.lowLeftLevel != -1 ) {
      if ( RangeUp > Node.lowLeftLevel + 1 )
        RangeUp = Node.lowLeftLevel + 1;
      if ( RangeLow < Node.lowLeftLevel - 1 )
        RangeLow = Node.lowLeftLevel - 1;
    }
    if ( Node.upRightLevel != -1 ) {
      if ( RangeUp > Node.upRightLevel + 1 )
        RangeUp = Node.upRightLevel + 1;
      if ( RangeLow < Node.upRightLevel - 1 )
        RangeLow = Node.upRightLevel - 1;
    }

    if ( RangeLow < 0 ) RangeLow = 0;
    if ( RangeUp > grid.getGridDepth() ) RangeUp = grid.getGridDepth();

    for ( int level = RangeLow; level < RangeUp; level++ ) {
      int step = 1 << ( grid.getGridDepth() - level );

      if ( quadTree->nodeInGrid ( Node.x() - step, Node.y(), level )
           && elementEstimator.checkUpLeftElement ( Node.x(), Node.y(), level ) ) {
        quadTree->setUpLeftElement ( Node.x(), Node.y(), level );
        updateElementToTrialNodes ( Node.x() - step, Node.y(), level );
        Node.upLeftLevel = level;
        return;
      }
    }

    int step = 1 << ( grid.getGridDepth() - RangeUp );
    quadTree->setUpLeftElement ( Node.x(), Node.y(), RangeUp );
    updateElementToTrialNodes ( Node.x() - step, Node.y(), RangeUp );
    Node.upLeftLevel = RangeUp;
  }


  void createElementLowLeft ( TrialNode &Node ) {
    if ( Node.lowRightLevel != -1 && Node.upLeftLevel != -1 ) {
      if ( aol::Abs ( Node.lowRightLevel - Node.upLeftLevel ) == 2 ) {
        Node.lowLeftLevel = ( Node.lowRightLevel + Node.upLeftLevel ) >> 1;
        quadTree->setLowLeftElement ( Node.x(), Node.y(), Node.lowLeftLevel );
        int step = 1 << ( grid.getGridDepth() - Node.lowLeftLevel );
        updateElementToTrialNodes ( Node.x() - step, Node.y() - step, Node.lowLeftLevel );
        return;
      }
    }

    int RangeUp = 1000, RangeLow = -1;
    if ( Node.lowRightLevel != -1 ) {
      if ( RangeUp > Node.lowRightLevel + 1 )
        RangeUp = Node.lowRightLevel + 1;
      if ( RangeLow < Node.lowRightLevel - 1 )
        RangeLow = Node.lowRightLevel - 1;
    }
    if ( Node.upLeftLevel != -1 ) {
      if ( RangeUp > Node.upLeftLevel + 1 )
        RangeUp = Node.upLeftLevel + 1;
      if ( RangeLow < Node.upLeftLevel - 1 )
        RangeLow = Node.upLeftLevel - 1;
    }

    if ( RangeLow < 0 ) RangeLow = 0;
    if ( RangeUp > grid.getGridDepth() ) RangeUp = grid.getGridDepth() ;

    for ( int level = RangeLow; level < RangeUp; level++ ) {
      int step = 1 << ( grid.getGridDepth() - level );

      if ( quadTree->nodeInGrid ( Node.x() - step, Node.y() - step, level )
           && elementEstimator.checkLowLeftElement ( Node.x(), Node.y(), level ) ) {
        quadTree->setLowLeftElement ( Node.x(), Node.y(), level );
        updateElementToTrialNodes ( Node.x() - step, Node.y() - step, level );
        Node.lowLeftLevel = level;
        return;
      }
    }

    int step = 1 << ( grid.getGridDepth() - RangeUp );
    quadTree->setLowLeftElement ( Node.x(), Node.y(), RangeUp );
    updateElementToTrialNodes ( Node.x() - step, Node.y() - step, RangeUp );
    Node.lowLeftLevel = RangeUp;
  }

  void createElementLowRight ( TrialNode &Node ) {
    if ( Node.lowLeftLevel != -1 && Node.upRightLevel != -1 ) {
      if ( aol::Abs ( Node.lowLeftLevel - Node.upRightLevel ) == 2 ) {
        Node.lowRightLevel = ( Node.lowLeftLevel + Node.upRightLevel ) >> 1;
        quadTree->setLowRightElement ( Node.x(), Node.y(), Node.lowRightLevel );
        int step = 1 << ( grid.getGridDepth() - Node.lowRightLevel );
        updateElementToTrialNodes ( Node.x(), Node.y() - step, Node.lowRightLevel );
        return;
      }
    }

    int RangeUp = 1000, RangeLow = -1;
    if ( Node.lowLeftLevel != -1 ) {
      if ( RangeUp > Node.lowLeftLevel + 1 )
        RangeUp = Node.lowLeftLevel + 1;
      if ( RangeLow < Node.lowLeftLevel - 1 )
        RangeLow = Node.lowLeftLevel - 1;
    }
    if ( Node.upRightLevel != -1 ) {
      if ( RangeUp > Node.upRightLevel + 1 )
        RangeUp = Node.upRightLevel + 1;
      if ( RangeLow < Node.upRightLevel - 1 )
        RangeLow = Node.upRightLevel - 1;
    }

    if ( RangeLow < 0 ) RangeLow = 0;
    if ( RangeUp > grid.getGridDepth() ) RangeUp = grid.getGridDepth();

    for ( int level = RangeLow; level < RangeUp; level++ ) {
      int step = 1 << ( grid.getGridDepth() - level );

      if ( quadTree->nodeInGrid ( Node.x(), Node.y() - step, level )
           && elementEstimator.checkLowRightElement ( Node.x(), Node.y(), level ) ) {
        quadTree->setLowRightElement ( Node.x(), Node.y(), level );
        updateElementToTrialNodes ( Node.x(), Node.y() - step, level );
        Node.lowRightLevel = level;
        return;
      }
    }

    int step = 1 << ( grid.getGridDepth() - RangeUp );
    quadTree->setLowRightElement ( Node.x(), Node.y(), RangeUp );
    updateElementToTrialNodes ( Node.x(), Node.y() - step, RangeUp );
    Node.lowRightLevel = RangeUp;
  }

  void updateElementToTrialNodes ( int X, int Y, int Level ) {
    typename vector<TrialNode>::iterator it;

    int FullStep = grid.h[ Level ];
    int HalfStep = FullStep >> 1;

    // for ( it = trialNodeHeap.data.begin( ); it != trialNodeHeap.data.end( ); ++it ) {
    for ( it = trialNodeHeap.begin( ); it != trialNodeHeap.end( ); ++it ) {
      TrialNode &Node = *it;

      if ( Node.x() == X && Node.y() == Y ) {
        Node.upRightLevel = Level;
      } else if ( Node.x() == X + FullStep && Node.y() == Y ) {
        Node.upLeftLevel = Level;
      } else if ( Node.x() == X + FullStep && Node.y() == Y + FullStep ) {
        Node.lowLeftLevel = Level;
      } else if ( Node.x() == X && Node.y() == Y + FullStep ) {
        Node.lowRightLevel = Level;
      }

      // now check hanging Nodes
      else if ( Node.x() == X + HalfStep && Node.y() == Y ) {
        Node.hangNode = TN_HN_UP;
      } else if ( Node.x() == X + HalfStep && Node.y() == Y + FullStep ) {
        Node.hangNode = TN_HN_LOW;
      } else if ( Node.x() == X && Node.y() == Y + HalfStep ) {
        Node.hangNode = TN_HN_RIGHT;
      } else if ( Node.x() == X + FullStep && Node.y() == Y + HalfStep ) {
        Node.hangNode = TN_HN_LEFT;
      }

    }
  }

  int searchOrAppendTrialNode ( TrialNode *&Node, int X, int Y ) {
    // for ( int i = 0; i < trialNodeHeap.data.size( ); i++ ) {
    for ( int i = 0; i < trialNodeHeap.size( ); i++ ) {
      if ( trialNodeHeap[ i ].x() == X && trialNodeHeap[ i ].y() == Y ) {
        Node = &trialNodeHeap[ i ];
        return i;
      }
    }
    TrialNode newNode ( X, Y );
    //trialNodeHeap.push_back( newNode );
    trialNodeHeap.append ( newNode );
    int LastIndex = trialNodeHeap.size( ) - 1;
    Node = & ( trialNodeHeap[ LastIndex ] );
    return LastIndex;
  }

  void marchStep( ) {
    if ( trialNodeHeap.size( ) == 0 ) return;
    TrialNode popNode;
    trialNodeHeap.pop ( popNode );

    makeTrialNodeActive ( popNode );
  }



  void trialNodesExport( ) {
    for ( int i = 0; i < trialNodeHeap.size( ); i++ ) {
      timeField->set ( trialNodeHeap[ i ].x(),
                       trialNodeHeap[ i ].y(),
                       trialNodeHeap[ i ].val );
    }
  }

  void dumpTrialNodes( ) {
    typename vector<TrialNode>::const_iterator it;
    for ( it = trialNodeHeap.begin( ); it != trialNodeHeap.end( ); ++it ) {
      const TrialNode &node = *it;
      cerr << "trialNode x = " << node.x() << "  y = " << node.y() << "  val = " << node.val << endl;
      cerr << "   UL = " << static_cast<int> ( node.upLeftLevel )
      << "  UR = " << static_cast<int> ( node.upRightLevel )
      << "  LL = " << static_cast<int> ( node.lowLeftLevel )
      << "  LR = " << static_cast<int> ( node.lowRightLevel )
      << "  hangNode = ";
      if ( node.hangNode == TN_HN_NONE ) {
        cerr << " NONE";
      } else if ( node.hangNode == TN_HN_UP ) {
        cerr << "UP   " << static_cast<int> ( node.hangNodeLevel );
      } else if ( node.hangNode == TN_HN_LOW ) {
        cerr << "LOW   " << static_cast<int> ( node.hangNodeLevel );
      } else if ( node.hangNode == TN_HN_LEFT ) {
        cerr << "LEFT   " << static_cast<int> ( node.hangNodeLevel );
      } else if ( node.hangNode == TN_HN_RIGHT ) {
        cerr << "RIGHT   " << static_cast<int> ( node.hangNodeLevel );
      }
      cerr << endl;
    }
  }

  RealType getSpeed ( int, int ) {
    return 1.0;
  }

  RealType getTime ( const aol::Vec3<short>& Node ) {
    return timeField->get ( Node );
  }

  void downwindLeftRight ( int X, int Y, int Step, RealType UpwindVal, RealType &Value );

  void downwindUpDown ( int X, int Y, int Step, RealType UpwindVal, RealType &Value );

  void downwindLeftRight ( int X, int Y, int Step, int UpwindY, RealType &Value );

  void downwindUpDown ( int X, int Y, int Step, int UpwindX, RealType &Value );

  void downwindLeftRight2ndOrder ( int X, int Y, int Step, int UpwindY, RealType &Value );

  void downwindUpDown2ndOrder ( int X, int Y, int Step, int UpwindX, RealType &Value );

  int getGridDepth( ) { return grid.getGridDepth(); }

  qc::QuadTree& getTree( ) { return *quadTree; }

  qc::Heap<TrialNode>& getTrialNodeHeap( ) { return trialNodeHeap; }

  qc::Array<RealType>& getTimeField( ) { return *timeField; }

protected:
  const qc::GridDefinition &grid;
  qc::Array<RealType> *timeField;
  RealType scaleTime;
  int order;

  qc::QuadTree *quadTree;
  qc::Array<RealType> *image;
  MetaEstimatorStruct   elementEstimator;
  qc::Heap<TrialNode> trialNodeHeap;
};

}   // end namespace eik

#endif

