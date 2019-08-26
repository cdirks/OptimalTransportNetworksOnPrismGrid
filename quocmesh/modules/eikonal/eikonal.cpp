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

#include <eikonal.h>

namespace eik {

template <typename RealType>
void Eikonal<RealType>::makeTrialNodeActive ( TrialNode &Node ) {
  if ( Node.xref() == 0 || Node.xref() == ( grid.getWidth() - 1 ) || Node.y() == 0 || Node.yref() == ( grid.getWidth() - 1 ) ) {
    RealType v = timeField->get ( Node.xref(), Node.yref() );
    if ( ( v > 0.0 && v < Node.val ) || v < 0.0 ) {
      timeField->set ( Node.xref(), Node.yref(), Node.val );
      return;
    }
  }
  if ( Node.xref() < 0 || Node.xref() > ( grid.getWidth() - 1 ) || Node.yref() < 0 || Node.yref() > ( grid.getWidth() - 1 ) )
    return;


  if ( Node.upLeftLevel == -1 &&
       Node.upRightLevel == -1 &&
       Node.lowLeftLevel == -1 &&
       Node.lowRightLevel == -1 ) {
    // create all surrounding elements
    // only needed for startup
    int mask = 0, step = 1;
    for ( int level = grid.getGridDepth(); level >= 0; level-- ) {

      if ( ! ( Node.xref() & mask ) && ! ( Node.yref() & mask )
           && elementEstimator.checkElement ( Node.xref(), Node.yref(), level ) ) {
        Node.upRightLevel = level;
        break;
      }
      mask |= step;
      step <<= 1;
    }

    quadTree->setUpRightElement ( Node.xref(), Node.yref(), Node.upRightLevel );

    // now UpLeft
    if ( quadTree->nodeInGrid ( Node.xref(), Node.yref(), Node.upRightLevel - 1 )
         && elementEstimator.checkUpLeftElement ( Node.xref(), Node.yref(), Node.upRightLevel - 1 ) ) {
      Node.upLeftLevel = Node.upRightLevel - 1;
    } else if ( quadTree->nodeInGrid ( Node.xref(), Node.yref(), Node.upRightLevel )
                && elementEstimator.checkUpLeftElement ( Node.xref(), Node.yref(), Node.upRightLevel ) ) {
      Node.upLeftLevel = Node.upRightLevel;
    } else if ( quadTree->nodeInGrid ( Node.xref(), Node.yref(), Node.upRightLevel + 1 )
                && elementEstimator.checkUpLeftElement ( Node.xref(), Node.yref(), Node.upRightLevel + 1 ) ) {
      Node.upLeftLevel = Node.upRightLevel + 1;
    } else {
      cerr << "Something strange has happened here! UPLEFT\n";
    }

    quadTree->setUpLeftElement ( Node.xref(), Node.yref(), Node.upLeftLevel );

    // no LowRight
    if ( quadTree->nodeInGrid ( Node.xref(), Node.yref(), Node.upRightLevel - 1 )
         && elementEstimator.checkLowRightElement ( Node.xref(), Node.yref(), Node.upRightLevel - 1 ) ) {
      Node.lowRightLevel = Node.upRightLevel - 1;
    } else if ( quadTree->nodeInGrid ( Node.xref(), Node.yref(), Node.upRightLevel )
                && elementEstimator.checkLowRightElement ( Node.xref(), Node.yref(), Node.upRightLevel ) ) {
      Node.lowRightLevel = Node.upRightLevel;
    } else if ( quadTree->nodeInGrid ( Node.xref(), Node.yref(), Node.upRightLevel + 1 )
                && elementEstimator.checkLowRightElement ( Node.xref(), Node.yref(), Node.upRightLevel + 1 ) ) {
      Node.lowRightLevel = Node.upRightLevel + 1;
    } else {
      cerr << "Something strange has happened here! LowRight\n";
    }

    quadTree->setLowRightElement ( Node.xref(), Node.yref(), Node.lowRightLevel );

    // no LowLeft
    if ( quadTree->nodeInGrid ( Node.xref(), Node.yref(), Node.upLeftLevel - 1 )
         && elementEstimator.checkLowLeftElement ( Node.xref(), Node.yref(), Node.upLeftLevel - 1 ) ) {
      Node.lowLeftLevel = Node.upLeftLevel - 1;
    } else if ( quadTree->nodeInGrid ( Node.xref(), Node.yref(), Node.upLeftLevel )
                && elementEstimator.checkLowLeftElement ( Node.xref(), Node.yref(), Node.upLeftLevel ) ) {
      Node.lowLeftLevel = Node.upLeftLevel;
    } else if ( quadTree->nodeInGrid ( Node.xref(), Node.yref(), Node.upLeftLevel + 1 )
                && elementEstimator.checkLowLeftElement ( Node.xref(), Node.yref(), Node.upLeftLevel + 1 ) ) {
      Node.lowLeftLevel = Node.upLeftLevel + 1;
    } else {
      cerr << "Something strange has happened here! LowLeft\n";
    }

    quadTree->setLowLeftElement ( Node.xref(), Node.yref(), Node.lowLeftLevel );
  }

  if ( Node.hangNode == TN_HN_NONE ) {
    // check for HangNode in left direction
    if ( Node.lowLeftLevel == -1 && Node.upLeftLevel == -1
         && ( Node.upRightLevel == -1 || Node.lowRightLevel == -1
              || Node.upRightLevel == Node.lowRightLevel ) ) {

      int coarseLevel = ( Node.upRightLevel > Node.lowRightLevel ) ?
                        ( Node.upRightLevel - 1 ) : ( Node.lowRightLevel - 1 );

      int HalfStep = 1 << ( grid.getGridDepth() - coarseLevel - 1 );

      if ( quadTree->nodeInGrid ( Node.xref(), Node.yref() - HalfStep, coarseLevel )
           && elementEstimator.checkUpLeftElement ( Node.xref(), Node.yref() - HalfStep, coarseLevel ) ) {
        if ( !quadTree->getElement ( Node.xref() - ( HalfStep << 1 ), Node.yref() - HalfStep, coarseLevel ) ) {
          quadTree->setElement ( Node.xref() - ( HalfStep << 1 ), Node.yref() - HalfStep, coarseLevel );
          updateElementToTrialNodes ( Node.xref() - ( HalfStep << 1 ), Node.yref() - HalfStep, coarseLevel );
        }
        Node.hangNode = TN_HN_LEFT;
        Node.hangNodeLevel = coarseLevel;
      }
    }
    // check for HangNode in right direction
    else if ( Node.lowRightLevel == -1 && Node.upRightLevel == -1
              && ( Node.upLeftLevel == -1 || Node.lowLeftLevel == -1
                   || Node.upLeftLevel == Node.lowLeftLevel ) ) {

      int coarseLevel = ( Node.upLeftLevel > Node.lowLeftLevel ) ?
                        ( Node.upLeftLevel - 1 ) : ( Node.lowLeftLevel - 1 );

      int HalfStep = 1 << ( grid.getGridDepth() - coarseLevel - 1 );

      if ( quadTree->nodeInGrid ( Node.xref(), Node.yref() - HalfStep, coarseLevel )
           && elementEstimator.checkUpRightElement ( Node.xref(), Node.yref() - HalfStep, coarseLevel ) ) {
        if ( !quadTree->getElement ( Node.xref(), Node.yref() - HalfStep, coarseLevel ) ) {
          quadTree->setElement ( Node.xref(), Node.yref() - HalfStep, coarseLevel );
          updateElementToTrialNodes ( Node.xref(), Node.yref() - HalfStep, coarseLevel );
        }
        Node.hangNode = TN_HN_RIGHT;
        Node.hangNodeLevel = coarseLevel;
      }
    }
    // check for hangNode in down direction
    else if ( Node.lowLeftLevel == -1 && Node.lowRightLevel == -1
              && ( Node.upLeftLevel == -1 || Node.upRightLevel == -1
                   || Node.upLeftLevel == Node.upRightLevel ) ) {

      int coarseLevel = ( Node.upLeftLevel > Node.upRightLevel ) ?
                        ( Node.upLeftLevel - 1 ) : ( Node.upRightLevel - 1 );

      int HalfStep = 1 << ( grid.getGridDepth()  - coarseLevel - 1 );

      if ( quadTree->nodeInGrid ( Node.xref() - HalfStep, Node.yref(), coarseLevel )
           && elementEstimator.checkLowRightElement ( Node.xref() - HalfStep, Node.yref(), coarseLevel ) ) {
        if ( !quadTree->getElement ( Node.xref() - HalfStep, Node.yref() - ( HalfStep << 1 ), coarseLevel ) ) {
          quadTree->setElement ( Node.xref() - HalfStep, Node.yref() - ( HalfStep << 1 ), coarseLevel );
          updateElementToTrialNodes ( Node.xref() - HalfStep, Node.yref() - ( HalfStep << 1 ), coarseLevel );
        }
        Node.hangNode = TN_HN_LOW;
        Node.hangNodeLevel = coarseLevel;
      }
    }
    // check for hangNode in up direction
    else if ( Node.upLeftLevel == -1 && Node.upRightLevel == -1
              && ( Node.lowLeftLevel == -1 || Node.lowRightLevel == -1
                   || Node.lowLeftLevel == Node.lowRightLevel ) ) {

      int coarseLevel = ( Node.lowLeftLevel > Node.lowRightLevel ) ?
                        ( Node.lowLeftLevel - 1 ) : ( Node.lowRightLevel - 1 );

      int HalfStep = 1 << ( grid.getGridDepth() - coarseLevel - 1 );

      if ( quadTree->nodeInGrid ( Node.xref() - HalfStep, Node.yref(), coarseLevel )
           && elementEstimator.checkUpRightElement ( Node.xref() - HalfStep, Node.yref(), coarseLevel ) ) {
        if ( !quadTree->getElement ( Node.xref() - HalfStep, Node.yref(), coarseLevel ) ) {
          quadTree->setElement ( Node.xref() - HalfStep, Node.yref(), coarseLevel );
          updateElementToTrialNodes ( Node.xref() - HalfStep, Node.yref(), coarseLevel );
        }
        Node.hangNode = TN_HN_UP;
        Node.hangNodeLevel = coarseLevel;
      }
    }
  } // end if ( Node.hangNode == TN_HN_NONE )

  // if this is a hanging node, then fill up the quadtree around this node.
  if ( Node.hangNode != TN_HN_NONE ) {
    cerr << "There's a hanging node\n";

    if ( Node.hangNode == TN_HN_LEFT ) {
      if ( Node.lowRightLevel == -1 ) {
        Node.lowRightLevel = Node.upRightLevel;
        int Step = quadTree->getStep ( Node.lowRightLevel );

        quadTree->setLowRightElement ( Node.xref(), Node.yref(), Node.lowRightLevel );
        updateElementToTrialNodes ( Node.xref(), Node.yref() - Step, Node.lowRightLevel );
      } else if ( Node.upRightLevel == -1 ) {
        Node.upRightLevel = Node.lowRightLevel;

        quadTree->setUpRightElement ( Node.xref(), Node.yref(), Node.upRightLevel );
        updateElementToTrialNodes ( Node.xref(), Node.yref(), Node.upRightLevel );
      }
    } else if ( Node.hangNode == TN_HN_RIGHT ) {
      if ( Node.lowLeftLevel == -1 ) {
        Node.lowLeftLevel = Node.upLeftLevel;
        int Step = quadTree->getStep ( Node.lowLeftLevel );

        quadTree->setLowLeftElement ( Node.xref(), Node.yref(), Node.lowLeftLevel );
        updateElementToTrialNodes ( Node.xref() - Step, Node.yref() - Step, Node.lowLeftLevel );
      } else if ( Node.upLeftLevel == -1 ) {
        Node.upLeftLevel = Node.lowLeftLevel;
        int Step = quadTree->getStep ( Node.upLeftLevel );

        quadTree->setUpLeftElement ( Node.xref(), Node.yref(), Node.upLeftLevel );
        updateElementToTrialNodes ( Node.xref() - Step, Node.yref(), Node.upLeftLevel );
      }
    } else if ( Node.hangNode == TN_HN_LOW ) {
      if ( Node.upLeftLevel == -1 ) {
        Node.upLeftLevel = Node.upRightLevel;
        int Step = quadTree->getStep ( Node.upLeftLevel );

        quadTree->setUpLeftElement ( Node.xref(), Node.yref(), Node.upLeftLevel );
        updateElementToTrialNodes ( Node.xref() - Step, Node.yref(), Node.upLeftLevel );
      } else if ( Node.upRightLevel == -1 ) {
        Node.upRightLevel = Node.upLeftLevel;

        quadTree->setUpRightElement ( Node.xref(), Node.yref(), Node.upRightLevel );
        updateElementToTrialNodes ( Node.xref(), Node.yref(), Node.upRightLevel );
      }
    } else if ( Node.hangNode == TN_HN_UP ) {
      if ( Node.lowLeftLevel == -1 ) {
        Node.lowLeftLevel = Node.lowRightLevel;
        int Step = quadTree->getStep ( Node.lowLeftLevel );

        quadTree->setLowLeftElement ( Node.xref(), Node.yref(), Node.lowLeftLevel );
        updateElementToTrialNodes ( Node.xref() - Step, Node.yref() - Step, Node.lowLeftLevel );
      } else if ( Node.lowRightLevel == -1 ) {
        Node.lowRightLevel = Node.lowLeftLevel;
        int Step = quadTree->getStep ( Node.lowRightLevel );

        quadTree->setLowRightElement ( Node.xref(), Node.yref(), Node.lowRightLevel );
        updateElementToTrialNodes ( Node.xref(), Node.yref() - Step, Node.lowRightLevel );
      }
    }
  } // endif ( Node.hangNode != TN_HN_NONE )

  int hasBeenSet = 0;

  if ( Node.hangNode == TN_HN_NONE ) {
    while ( Node.upLeftLevel == -1
            || Node.upRightLevel == -1
            || Node.lowLeftLevel == -1
            || Node.lowRightLevel == -1 ) {

      if ( Node.upRightLevel != -1
           || Node.lowLeftLevel != -1 ) {
        if ( Node.upLeftLevel == - 1 ) createElementUpLeft ( Node );
        if ( Node.lowRightLevel == - 1 ) createElementLowRight ( Node );
      }

      if ( Node.upLeftLevel != -1
           || Node.lowRightLevel != - 1 ) {
        if ( Node.upRightLevel == - 1 ) createElementUpRight ( Node );
        if ( Node.lowLeftLevel == - 1 ) createElementLowLeft ( Node );
      }
    }
    RealType v = timeField->get ( Node.xref(), Node.yref() );
    if ( v < 0.0f || Node.val < v ) {
      timeField->set ( Node.xref(), Node.yref(), Node.val );
      hasBeenSet = 1;
    }
  } else {
    return;
  }

  if ( hasBeenSet && Node.upLeftLevel != Node.upRightLevel ) {
    int Level = aol::Min ( Node.upLeftLevel, Node.upRightLevel );
    int Step = 1 << ( grid.getGridDepth() - Level );

    RealType v;
    if ( ( v = timeField->get ( Node.xref(), Node.yref() + Step ) ) >= 0.0f ) {

      timeField->set ( Node.xref(), Node.yref() + ( Step >> 1 ), 0.5f * ( Node.val + v ) );

      if ( Node.upLeftLevel < Node.upRightLevel ) {
        // go to the right
        TrialNode *newNode;
        int Index = searchOrAppendTrialNode ( newNode,
                                              Node.xref() + ( Step >> 1 ),
                                              Node.yref() + ( Step >> 1 ) );
        newNode->lowLeftLevel = newNode->upLeftLevel = Level + 1;
        downwindUpDown ( Node.xref() + ( Step >> 1 ),
                         Node.yref() + ( Step >> 1 ),
                         Step >> 1, 0.5f * ( Node.val + v ), newNode->val );
        trialNodeHeap.upHeap ( Index );
      } else {
        TrialNode *newNode;
        int Index = searchOrAppendTrialNode ( newNode,
                                              Node.xref() - ( Step >> 1 ),
                                              Node.yref() + ( Step >> 1 ) );
        newNode->lowRightLevel = newNode->upRightLevel = Level + 1;
        downwindUpDown ( Node.xref() - ( Step >> 1 ),
                         Node.yref() + ( Step >> 1 ),
                         Step >> 1, 0.5f * ( Node.val + v ), newNode->val );
        trialNodeHeap.upHeap ( Index );
      }
    }
  }

  if ( hasBeenSet && Node.lowLeftLevel != Node.lowRightLevel ) {
    int Level = aol::Min ( Node.lowLeftLevel, Node.lowRightLevel );
    int Step = 1 << ( grid.getGridDepth() - Level );

    RealType v;
    if ( ( v = timeField->get ( Node.xref(), Node.yref() - Step ) ) >= 0.0f ) {

      //RealType d = timeField->get( Node.xref(), Node.yref() - ( Step >> 1 )) ;

      timeField->set ( Node.xref(), Node.yref() - ( Step >> 1 ), 0.5f * ( Node.val + v ) );

      //cerr << "Interpolated value = " << 0.5f * (Node.val + v) << " old = " << d << endl;

      if ( Node.lowLeftLevel < Node.lowRightLevel ) {
        // go to the right
        TrialNode *newNode;
        int Index = searchOrAppendTrialNode ( newNode,
                                              Node.xref() + ( Step >> 1 ),
                                              Node.yref() - ( Step >> 1 ) );
        newNode->lowLeftLevel = newNode->upLeftLevel = Level + 1;
        downwindUpDown ( Node.xref() + ( Step >> 1 ),
                         Node.yref() - ( Step >> 1 ),
                         Step >> 1, 0.5f * ( Node.val + v ), newNode->val );
        trialNodeHeap.upHeap ( Index );
      } else {
        TrialNode *newNode;
        int Index = searchOrAppendTrialNode ( newNode,
                                              Node.xref() - ( Step >> 1 ),
                                              Node.yref() - ( Step >> 1 ) );
        newNode->lowRightLevel = newNode->upRightLevel = Level + 1;
        downwindUpDown ( Node.xref() - ( Step >> 1 ),
                         Node.yref() - ( Step >> 1 ),
                         Step >> 1, 0.5f * ( Node.val + v ), newNode->val );
        trialNodeHeap.upHeap ( Index );
      }
    }
  }

  if ( hasBeenSet && Node.upLeftLevel != Node.lowLeftLevel ) {
    int Level = aol::Min ( Node.upLeftLevel, Node.lowLeftLevel );
    int Step = 1 << ( grid.getGridDepth() - Level );

    RealType v;
    if ( ( v = timeField->get ( Node.xref() - Step, Node.yref() ) ) >= 0.0f ) {

      timeField->set ( Node.xref() - ( Step >> 1 ), Node.yref(), 0.5f * ( Node.val + v ) );

      // cerr << "Interpolated value = " << 0.5f * (Node.val + v) << endl;

      if ( Node.upLeftLevel < Node.lowLeftLevel ) {
        // go down
        TrialNode *newNode;
        int Index = searchOrAppendTrialNode ( newNode, Node.xref() - ( Step >> 1 ), Node.yref() - ( Step >> 1 ) );
        newNode->upLeftLevel = newNode->upRightLevel = Level + 1;
        downwindLeftRight ( Node.xref() - ( Step >> 1 ),
                            Node.yref() - ( Step >> 1 ),
                            Step >> 1, 0.5f * ( Node.val + v ), newNode->val );
        trialNodeHeap.upHeap ( Index );

      } else {
        TrialNode *newNode;
        int Index = searchOrAppendTrialNode ( newNode,
                                              Node.xref() - ( Step >> 1 ),
                                              Node.yref() + ( Step >> 1 ) );
        newNode->lowLeftLevel = newNode->lowRightLevel = Level + 1;
        downwindLeftRight ( Node.xref() - ( Step >> 1 ),
                            Node.yref() + ( Step >> 1 ),
                            Step >> 1, 0.5f * ( Node.val + v ), newNode->val );
        trialNodeHeap.upHeap ( Index );
      }
    }
  }

  if ( hasBeenSet && Node.upRightLevel != Node.lowRightLevel ) {
    int Level = aol::Min ( Node.upRightLevel, Node.lowRightLevel );
    int Step = 1 << ( grid.getGridDepth() - Level );

    RealType v;
    if ( ( v = timeField->get ( Node.xref() + Step, Node.yref() ) ) >= 0.0f ) {

      timeField->set ( Node.xref() + ( Step >> 1 ), Node.yref(), 0.5f * ( Node.val + v ) );

      // cerr << "Interpolated value = " << 0.5f * (Node.val + v) << endl;

      if ( Node.upRightLevel < Node.lowRightLevel ) {
        // go down
        TrialNode *newNode;
        int Index = searchOrAppendTrialNode ( newNode, Node.xref() + ( Step >> 1 ), Node.yref() - ( Step >> 1 ) );
        newNode->upLeftLevel = newNode->upRightLevel = Level + 1;
        downwindLeftRight ( Node.xref() + ( Step >> 1 ), Node.yref() - ( Step >> 1 ), Step >> 1, 0.5f * ( Node.val + v ), newNode->val );
        trialNodeHeap.upHeap ( Index );
      } else {
        TrialNode *newNode;
        int Index = searchOrAppendTrialNode ( newNode, Node.xref() + ( Step >> 1 ), Node.yref() + ( Step >> 1 ) );
        newNode->lowLeftLevel = newNode->lowRightLevel = Level + 1;
        downwindLeftRight ( Node.xref() + ( Step >> 1 ), Node.yref() + ( Step >> 1 ), Step >> 1, 0.5f * ( Node.val + v ), newNode->val );
        trialNodeHeap.upHeap ( Index );
      }
    }
  }

  // jetzt wird's erst einigermassen spannend:
  // kreiere neue TrialNodes und schiebe sie in den heap,
  // beziehungsweise update der Zeitwerte.


  if ( Node.hangNode != TN_HN_UP ) {
    //int level = imax( Node.upLeftLevel, Node.upRightLevel );
    int level = aol::Min ( Node.upLeftLevel, Node.upRightLevel );

    int step = grid.h[ level ];

    if ( Node.y() + step <= grid.getWidth() - 1 && timeField->get ( Node.x(), Node.y() + step ) < 0.0f ) {

      TrialNode *newNode;
      int Index = searchOrAppendTrialNode ( newNode, Node.xref(), Node.yref() + step );

      // first case: emanating from a regular node
      newNode->lowLeftLevel = Node.upLeftLevel;
      newNode->lowRightLevel = Node.upRightLevel;
      if ( newNode->lowLeftLevel > newNode->lowRightLevel ) {
        quadTree->setLowLeftElement ( newNode->x(), newNode->y(), newNode->lowLeftLevel );
        updateElementToTrialNodes ( newNode->x() - ( 1 << ( grid.getGridDepth() - newNode->lowLeftLevel ) ),
                                    newNode->y() - ( 1 << ( grid.getGridDepth() - newNode->lowLeftLevel ) ),
                                    newNode->lowLeftLevel );
      } else if ( newNode->lowRightLevel > newNode->lowLeftLevel ) {
        quadTree->setLowRightElement ( newNode->x(), newNode->y(), newNode->lowRightLevel );
        updateElementToTrialNodes ( newNode->x(),
                                    newNode->y() - ( 1 << ( grid.getGridDepth() - newNode->lowRightLevel ) ),
                                    newNode->lowRightLevel );
      }

      switch ( order ) {
      case 1:
        downwindLeftRight2ndOrder ( newNode->x(), newNode->y(), step, Node.yref(), newNode->val );
        break;
      case 2:
        downwindLeftRight ( newNode->x(), newNode->y(), step, Node.yref(), newNode->val );
        break;
      default:
        downwindLeftRight ( newNode->x(), newNode->y(), step, Node.val, newNode->val );
      }
      trialNodeHeap.upHeap ( Index );
    }
  }

  if ( Node.hangNode != TN_HN_LOW ) {
    //int level = imax( Node.upLeftLevel, Node.upRightLevel );
    int level = aol::Min ( Node.lowLeftLevel, Node.lowRightLevel );

    int step = 1 << ( grid.getGridDepth() - level );

    if ( Node.yref() - step >= 0 && timeField->get ( Node.xref(), Node.yref() - step ) < 0.0f ) {

      TrialNode *newNode;
      int Index = searchOrAppendTrialNode ( newNode, Node.xref(), Node.yref() - step );

      // first case: emanating from a regular node
      newNode->upLeftLevel = Node.lowLeftLevel;
      newNode->upRightLevel = Node.lowRightLevel;
      if ( newNode->upLeftLevel > newNode->upRightLevel ) {
        quadTree->setUpLeftElement ( newNode->x(), newNode->y(), newNode->upLeftLevel );
        updateElementToTrialNodes ( newNode->x() - ( 1 << ( grid.getGridDepth() - newNode->upLeftLevel ) ),
                                    newNode->y(),
                                    newNode->upLeftLevel );
      } else if ( newNode->upRightLevel > newNode->upLeftLevel ) {
        quadTree->setUpRightElement ( newNode->x(), newNode->y(), newNode->upRightLevel );
        updateElementToTrialNodes ( newNode->x(),
                                    newNode->y(),
                                    newNode->upRightLevel );
      }

      switch ( order ) {
      case 1:
        downwindLeftRight2ndOrder ( newNode->x(), newNode->y(), step, Node.yref(), newNode->val );
        break;
      case 2:
        downwindLeftRight ( newNode->x(), newNode->y(), step, Node.yref(), newNode->val );
        break;
      default:
        downwindLeftRight ( newNode->x(), newNode->y(), step, Node.val, newNode->val );
      }
      trialNodeHeap.upHeap ( Index );
    }
  }

  if ( Node.hangNode != TN_HN_LEFT ) {
    //int level = imax( Node.upLeftLevel, Node.upRightLevel );
    int level = aol::Min ( Node.lowLeftLevel, Node.upLeftLevel );

    int step = 1 << ( grid.getGridDepth() - level );

    if ( Node.xref() - step >= 0 && timeField->get ( Node.xref() - step, Node.yref() ) < 0.0f ) {

      TrialNode *newNode;
      int Index = searchOrAppendTrialNode ( newNode, Node.xref() - step, Node.yref() );

      // first case: emanating from a regular node
      newNode->upRightLevel = Node.upLeftLevel;
      newNode->lowRightLevel = Node.lowLeftLevel;
      if ( newNode->lowRightLevel > newNode->upRightLevel ) {
        quadTree->setLowRightElement ( newNode->x(), newNode->y(), newNode->lowRightLevel );
        updateElementToTrialNodes ( newNode->x(),
                                    newNode->y() - ( 1 << ( grid.getGridDepth() - newNode->lowRightLevel ) ),
                                    newNode->lowRightLevel );
      } else if ( newNode->upRightLevel > newNode->lowRightLevel ) {
        quadTree->setUpRightElement ( newNode->x(), newNode->y(), newNode->upRightLevel );
        updateElementToTrialNodes ( newNode->x(),
                                    newNode->y(),
                                    newNode->upRightLevel );
      }


      switch ( order ) {
      case 1:
        downwindUpDown2ndOrder ( newNode->x(), newNode->y(), step, Node.xref(), newNode->val );
        break;
      case 2:
        downwindUpDown ( newNode->x(), newNode->y(), step, Node.xref(), newNode->val );
        break;
      default:
        downwindUpDown ( newNode->x(), newNode->y(), step, Node.val, newNode->val );
      }
      trialNodeHeap.upHeap ( Index );
    }
  }

  if ( Node.hangNode != TN_HN_RIGHT ) {
    //int level = imax( Node.upLeftLevel, Node.upRightLevel );
    int level = aol::Min ( Node.lowRightLevel, Node.upRightLevel );

    int step = 1 << ( grid.getGridDepth() - level );

    if ( Node.xref() + step <= grid.getWidth() - 1 && timeField->get ( Node.xref() + step, Node.yref() ) < 0.0f ) {

      TrialNode *newNode;
      int Index = searchOrAppendTrialNode ( newNode, Node.xref() + step, Node.yref() );

      // first case: emanating from a regular node
      newNode->upLeftLevel = Node.upRightLevel;
      newNode->lowLeftLevel = Node.lowRightLevel;
      if ( newNode->lowLeftLevel > newNode->upLeftLevel ) {
        quadTree->setLowLeftElement ( newNode->x(), newNode->y(), newNode->lowLeftLevel );
        updateElementToTrialNodes ( newNode->x() - ( 1 << ( grid.getGridDepth() - newNode->lowLeftLevel ) ),
                                    newNode->y() - ( 1 << ( grid.getGridDepth() - newNode->lowLeftLevel ) ),
                                    newNode->lowLeftLevel );
      } else if ( newNode->upLeftLevel > newNode->lowLeftLevel ) {
        quadTree->setUpLeftElement ( newNode->x(), newNode->y(), newNode->upLeftLevel );
        updateElementToTrialNodes ( newNode->x() - ( 1 << ( grid.getGridDepth() - newNode->upLeftLevel ) ),
                                    newNode->y(),
                                    newNode->upLeftLevel );
      }


      switch ( order ) {
      case 1:
        downwindUpDown2ndOrder ( newNode->x(), newNode->y(), step, Node.xref(), newNode->val );
        break;
      case 2:
        downwindUpDown ( newNode->x(), newNode->y(), step, Node.xref(), newNode->val );
        break;
      default:
        downwindUpDown ( newNode->x(), newNode->y(), step, Node.val, newNode->val );
      }
      trialNodeHeap.upHeap ( Index );
    }
  }
}

template <typename RealType>
RealType Eikonal<RealType>::solveQuadratic ( RealType v1, RealType v2, RealType h1, RealType h2, RealType F ) {
  /* solve ( T - v1 )^2 / h1^2 + ( T - v2 )^2/h2^2 = 1.0/F^2
     <=> (T^2 - 2*T*v1 + v1^2) * h2^2 + ( T^2 - 2*T*v2 + v2^2) * h1^2 = h1^2 * h2 ^2 / F^2
     <=> (h1^2 + h2^2) * T^2 - 2 * (v1 * h2^2 + v2 * h1^2 ) * T+ h2^2 * v1^2 + h1^2 * v2^2 - h1^2 * h2 ^2 / F^2 = 0 */


  RealType d1 = h1 * h1;
  RealType d2 = h2 * h2;
  RealType a = d1 + d2;
  RealType b = -2 * ( v1 * d2 + v2 * d1 );
  RealType c =  d2 * v1 * v1 + d1 * v2 * v2 - d1 * d2 / ( F * F );

  RealType q = b * b - 4 * a * c;

  if ( q < 0.0f || a == 0.0f ) {
    /*
    cerr << "q = " << q << " v1 = "
    << v1 << "  v2 = " << v2
    << "  h1 = " << h1 << "  h2 = "
    << h2 << "  F = " << F << endl;
    */
    return -1.0f;
  }

  RealType r = sqrt ( q );

  if ( r > b ) {
    return ( r - b ) / ( 2.0f * a );
  } else if ( r + b < 0.0f ) {
    return ( - r - b ) / ( 2.0f * a );
  } else {
    return 0.0f;
  }
  /*s = ( r - b ) / ( 2.0f * a );
  if ( s < 0.0f ) {
    s = ( - r - b ) / ( 2.0f * a );
  }*/
  /*cerr << "SolveQuadratic: v1 = " << v1 << "  v2 = " << v2
    << "  h1 = " << h1 << "  h2 = " << h2 << "  s = " << s << endl;

    return s;
  */


}

template <typename RealType>
void Eikonal<RealType>::downwindUpDown ( int X, int Y, int Step, RealType UpwindVal, RealType &Value ) {
  RealType F = getSpeed ( X, Y );
  RealType H = RealType ( Step ) / RealType ( ( grid.getWidth() - 1 ) );

  RealType v = timeField->get ( X, Y + Step );
  if ( v > 0.0f ) {
    RealType newH = RealType ( Step ) / RealType ( ( grid.getWidth() - 1 ) );
    RealType newV = solveQuadratic ( UpwindVal, v, H, newH, F );
    if ( ( Value < 0.0f || Value > newV ) && newV >= 0.0f ) {
      Value = newV;
    }
  }

  v = timeField->get ( X, Y - Step );
  if ( v > 0.0f ) {
    RealType newH = RealType ( Step ) / RealType ( ( grid.getWidth() - 1 ) );
    RealType newV = solveQuadratic ( UpwindVal, v, H, newH, F );
    if ( ( Value < 0.0f || Value > newV ) && newV >= 0.0f ) {
      Value = newV;
    }
  }

  if ( Value <= 0.0f ) {
    Value = UpwindVal + H / F;
  }
}

template <typename RealType>
void Eikonal<RealType>::downwindLeftRight ( int X, int Y, int Step, RealType UpwindVal, RealType &Value ) {
  RealType F = getSpeed ( X, Y );
  RealType H = RealType ( Step ) / RealType ( ( grid.getWidth() - 1 ) );

  RealType v = timeField->get ( X - Step, Y );
  if ( v > 0.0f ) {
    RealType newH = RealType ( Step ) / RealType ( ( grid.getWidth() - 1 ) );
    RealType newV = solveQuadratic ( UpwindVal, v, H, newH, F );
    if ( ( Value < 0.0f || Value > newV ) && newV >= 0.0f ) {
      Value = newV;
    }
  }

  v = timeField->get ( X + Step, Y );
  if ( v > 0.0f ) {
    RealType newH = RealType ( Step ) / RealType ( ( grid.getWidth() - 1 ) );
    RealType newV = solveQuadratic ( UpwindVal, v, H, newH, F );
    if ( ( Value < 0.0f || Value > newV ) && newV >= 0.0f ) {
      Value = newV;
    }
  }

  if ( Value <= 0.0f ) {
    Value = UpwindVal + H / F;
  }
}

template <typename RealType>
void Eikonal<RealType>::downwindLeftRight ( int X, int Y, int Step,
                                            int UpwindY, RealType &Value ) {

  RealType h = RealType ( Step ) / RealType ( ( grid.getWidth() - 1 ) );
  RealType a, b, c, t;

  RealType UpwindVal = timeField->get ( X, UpwindY );
  RealType F = getSpeed ( X, Y );

  if ( UpwindVal < 0.0 ) cerr << "\n\n\n>>>> ERROR Upwindval < 0\n\n\n";

  RealType v = timeField->get ( X + Step, Y );
  if ( v >= 0.0f ) {
    RealType OrgVal = timeField->get ( X + Step, UpwindY );

    a = 1.0f;
    b = -2.0f * OrgVal;

    t = v - UpwindVal;
    c = 0.5f * ( ( t - OrgVal ) * ( t - OrgVal ) + ( t + OrgVal ) * ( t + OrgVal ) ) - 2.0f * h * h / ( F * F );

    RealType newV = solveQuadratic ( a, b, c );
    if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
      Value = newV;
    }
  }

  v = timeField->get ( X - Step, Y );
  if ( v >= 0.0f ) {
    RealType OrgVal = timeField->get ( X - Step, UpwindY );

    a = 1.0f;
    b = -2.0f * OrgVal;

    t = v - UpwindVal;
    c = 0.5f * ( ( t - OrgVal ) * ( t - OrgVal ) + ( t + OrgVal ) * ( t + OrgVal ) ) - 2.0f * h * h / ( F * F );

    RealType newV = solveQuadratic ( a, b, c );
    if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
      Value = newV;
    }
  }
  if ( Value < 0.0f ) {
    Value = UpwindVal + h / F;
  }
}

template <typename RealType>
void Eikonal<RealType>::downwindUpDown ( int X, int Y, int Step,
                                         int UpwindX, RealType &Value ) {

  RealType h = RealType ( Step ) / RealType ( ( grid.getWidth() - 1 ) );
  RealType a, b, c, t;

  RealType F = getSpeed ( X, Y );
  RealType UpwindVal = timeField->get ( UpwindX, Y );

  if ( UpwindVal < 0.0 ) cerr << "\n\n\n>>>> ERROR Upwindval < 0\n\n\n";

  RealType v = timeField->get ( X, Y + Step );
  if ( v >= 0.0 ) {
    RealType OrgVal = timeField->get ( UpwindX, Y + Step );

    a = 1.0f;
    b = -2 * OrgVal;

    t = v - UpwindVal;
    c = 0.5f * ( ( t - OrgVal ) * ( t - OrgVal ) + ( t + OrgVal ) * ( t + OrgVal ) ) -  2.0f * h * h / ( F * F );

    RealType newV = solveQuadratic ( a, b, c );
    if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
      Value = newV;
    }
  }

  v = timeField->get ( X, Y - Step );
  if ( v >= 0.0f ) {
    RealType OrgVal = timeField->get ( UpwindX, Y - Step );

    a = 1.0f;
    b = -2.0f * OrgVal;

    t = v - UpwindVal;
    c = 0.5f * ( ( t - OrgVal ) * ( t - OrgVal ) + ( t + OrgVal ) * ( t + OrgVal ) ) - 2.0f * h * h / ( F * F );

    RealType newV = solveQuadratic ( a, b, c );
    if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
      Value = newV;
    }
  }
  if ( Value < 0.0f ) {
    Value = UpwindVal + h / F;
  }
}

template <typename RealType>
RealType Eikonal<RealType>::solveQuadratic ( RealType a, RealType b, RealType c ) {
  RealType d = b * b - 4 * a * c;

  if ( d < 0.0 ) {
    cerr << "ERROR in solveQuadratic: not solvable.\n";
    return -1.0;
  }

  d = sqrt ( d );

  RealType r1, r2;
  r1 = 0.5f * ( -b + d ) / a;
  r2 = 0.5f * ( -b - d ) / a;

  if ( r1 < 0.0  ) return r2;
  else if ( r2 < 0.0 ) return r1;
  else return aol::Max ( r1, r2 );
}

template <typename RealType>
void Eikonal<RealType>::downwindLeftRight2ndOrder ( int X, int Y, int Step,
                                                    int UpwindY, RealType &Value ) {

  RealType h = RealType ( Step ) / RealType ( ( grid.getWidth() - 1 ) );
  RealType a, b, c, ax, ay, bx, by;

  RealType UpwindY1 = timeField->get ( X, UpwindY );
  RealType UpwindY2 = timeField->get ( X, UpwindY + UpwindY - Y );
  RealType F = getSpeed ( X, Y );

  if ( UpwindY1 < 0.0 ) cerr << "\n\n\n>>>> ERROR Upwindval < 0\n\n\n";

  if ( UpwindY2 < UpwindY1 && UpwindY2 > 0.0 ) {
    ay = 1.5;
    by = 2 * UpwindY1 - 0.5f * UpwindY2;
  } else {
    ay = 1.0;
    by = UpwindY1;
  }

  RealType x1 = timeField->get ( X + Step, Y );
  if ( x1 >= 0.0 ) {
    RealType x2 = timeField->get ( X + Step + Step, Y );
    if ( x2 < x1 && x2 > 0.0 ) {
      ax = 1.5;
      bx = 2 * x1 - 0.5f * x2;
    } else {
      ax = 1.0;
      bx = x1;
    }
    a = ax * ax + ay * ay;
    b = -2 * ( ax * bx + ay * by );
    c = bx * bx + by * by - h * h / ( F * F );

    RealType newV = solveQuadratic ( a, b, c );
    if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
      Value = newV;
    }
  }

  x1 = timeField->get ( X - Step, Y );
  if ( x1 >= 0.0f ) {
    RealType x2 = timeField->get ( X - Step - Step, Y );
    if ( x2 < x1 && x2 > 0.0 ) {
      ax = 1.5;
      bx = 2 * x1 - 0.5f * x2;
    } else {
      ax = 1.0;
      bx = x1;
    }
    a = ax * ax + ay * ay;
    b = -2 * ( ax * bx + ay * by );
    c = bx * bx + by * by - h * h / ( F * F );

    RealType newV = solveQuadratic ( a, b, c );
    if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
      Value = newV;
    }
  }

  if ( Value < 0.0f ) {
    if ( UpwindY2 < UpwindY1 && UpwindY2 > 0.0 ) {
      Value = ( 2 * h / F + 4 * UpwindY1  -  UpwindY2 ) / 3.0f;
    } else {
      Value = UpwindY1 + h / F;
    }
  }
}

template <typename RealType>
void Eikonal<RealType>::downwindUpDown2ndOrder ( int X, int Y, int Step,
                                                 int UpwindX, RealType &Value ) {

  RealType h = RealType ( Step ) / RealType ( ( grid.getWidth() - 1 ) );
  RealType a, b, c, ax, bx, ay, by;

  RealType UpwindX1 = timeField->get ( UpwindX, Y );
  RealType UpwindX2 = timeField->get ( UpwindX + UpwindX - X, Y );
  RealType F = getSpeed ( X, Y );

  if ( UpwindX1 < 0.0 ) cerr << "\n\n\n>>>> ERROR Upwindval < 0\n\n\n";

  if ( UpwindX2 < UpwindX1 && UpwindX2 > 0.0 ) {
    ax = 1.5;
    bx = 2 * UpwindX1 - 0.5f * UpwindX2;
  } else {
    ax = 1.0;
    bx = UpwindX1;
  }

  RealType y1 = timeField->get ( X, Y + Step );
  if ( y1 >= 0.0f ) {
    RealType y2 = timeField->get ( X, Y + Step + Step );
    if ( y2 < y1 && y2 > 0.0 ) {
      ay = 1.5;
      by = 2 * y1 - 0.5f * y2;
    } else {
      ay = 1.0;
      by = y1;
    }
    a = ax * ax + ay * ay;
    b = -2 * ( ax * bx + ay * by );
    c = bx * bx + by * by - h * h / ( F * F );

    RealType newV = solveQuadratic ( a, b, c );
    if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
      Value = newV;
    }
  }

  y1 = timeField->get ( X, Y - Step );
  if ( y1 >= 0.0f ) {
    RealType y2 = timeField->get ( X, Y - Step - Step );
    if ( y2 < y1 && y2 > 0.0 ) {
      ay = 1.5;
      by = 2 * y1 - 0.5f * y2;
    } else {
      ay = 1.0;
      by = y1;
    }
    a = ax * ax + ay * ay;
    b = -2 * ( ax * bx + ay * by );
    c = bx * bx + by * by - h * h / ( F * F );

    RealType newV = solveQuadratic ( a, b, c );
    if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
      Value = newV;
    }
  }

  if ( Value < 0.0f ) {
    if ( UpwindX2 < UpwindX1 && UpwindX2 > 0.0 ) {
      Value = ( 2 * h / F + 4 * UpwindX1  -  UpwindX2 ) / 3.0f;
    } else {
      Value = UpwindX1 + h / F;
    }
  }
}

template class Eikonal<float>;

}   // end namespace eik
