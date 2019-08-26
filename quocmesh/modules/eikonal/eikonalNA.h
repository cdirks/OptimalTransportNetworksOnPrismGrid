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

#ifndef __EIKONALNA_H
#define __EIKONALNA_H

#include <quoc.h>
#include <qmHeap.h>
#include <gridBase.h>
#include <scalarArray.h>
#include <anisotropies.h>

namespace eik {

template <typename RealType, class EikonalImp, qc::Dimension, typename GridType = qc::GridDefinition>
class EikonalNAInt { };


template <typename RealType, class EikonalImp, typename GridType>
class EikonalNAInt< RealType, EikonalImp , qc::QC_2D, GridType > {
  /*** Barton Nackman ***/
  EikonalImp &asImp( ) { return static_cast<EikonalImp&> ( *this ); }
  const EikonalImp &asImp( ) const { return static_cast<const EikonalImp&> ( *this ); }

public:

class TrialNode : public aol::Vec3<short> {
  public:
    TrialNode ( aol::Vec3<short> Coord )
        : aol::Vec3<short> ( Coord ),
        val ( -1.0 )  {}

    TrialNode ( short X, short Y )
        : aol::Vec3<short> ( X, Y, 0 ),
        val ( -1.0 )  {}

    TrialNode( )
        : val ( -1.0 ) {}

    bool operator== ( const TrialNode &Other ) const {
      return ( val == Other.val );
    }

    bool operator< ( const TrialNode &Other ) const {
      return ( val < Other.val );
    }

    bool operator> ( const TrialNode &Other ) const {
      return ( val > Other.val );
    }

    RealType val;

  };

  //! constructor with grid definition
  EikonalNAInt ( const GridType &Grid )
      : grid ( Grid ) ,
      timeField ( NULL ),
      indexField ( Grid ),
      order ( 1 ) {
    timeField = new qc::ScalarArray<RealType, qc::QC_2D> ( Grid );
    trialNodeHeap.setIndexField ( &indexField );

    reset();
  }

  ~EikonalNAInt( ) {
    if ( timeField != NULL )
      delete timeField;
  }

  void reset( ) {
    trialNodeHeap.erase ( trialNodeHeap.begin( ), trialNodeHeap.end( ) );
    timeField->setAll ( -1.0f );
    indexField.setAll ( -1 );
  }

  void setOrder ( int Order ) { order = Order; }

  void setSeedPoint ( const aol::Vec3<short> Point, RealType Value = 0. ) {
    TrialNode *Node;
    int index = searchOrAppendTrialNode ( Node, Point.x(), Point.y() );
    if ( Node->val > 0. ) {
      Node->val = aol::Min ( Node->val, Value );
    } else {
      Node->val = Value;
    }

    trialNodeHeap.upHeap ( index );
  }

  void march( ) {
    while ( trialNodeHeap.size() > 0 ) {
      marchStep();
    }
  }


  void trialNodesExport( ) {
    timeField->setAll ( 1. );
    for ( int i = 0; i < trialNodeHeap.size( ); i++ ) {
      timeField->set ( trialNodeHeap[ i ].x(),
                       trialNodeHeap[ i ].y(),
                       trialNodeHeap[ i ].val );
    }
  }

  void checkTrialNodes( ) {
    for ( int i = 0; i < trialNodeHeap.size( ); i++ ) {
      cerr << trialNodeHeap[i].x() << " " << trialNodeHeap[i].x() << " " << trialNodeHeap[i].val << endl;
      if ( trialNodeHeap[ i ].val <= 0. ) {
        cerr << "trialNode smaller than zero!\n";
        trialNodeHeap[ i ].val = - trialNodeHeap[ i ].val;
      }
    }
  }

protected:

  void updateExtensionVelocity ( const aol::Vec2<short> &origNode,
                                 const aol::Vec2<short> &newNode,
                                 const aol::Vec2<short> &thirdNode,
                                 RealType newValue ) {
    // barton-nackman
    asImp().updateExtensionVelocity ( origNode, newNode, thirdNode, newValue );
  }

  void updateExtensionVelocity ( const aol::Vec2<short> &newNode, const aol::Vec2<short> &oldNode ) {
    asImp().updateExtensionVelocity ( newNode, oldNode );
  }

  void newNodeInserted ( const aol::Vec2<short> &newNode ) {
    asImp().newNodeInserted ( newNode );
  }

  void makeTrialNodeActive ( TrialNode &Node );

  RealType solveQuadratic ( RealType v1, RealType v2, RealType h1, RealType h2, RealType F = 1.0f );

  inline RealType solveQuadratic ( RealType a, RealType b, RealType c );

  int searchOrAppendTrialNode ( TrialNode *&Node, int X, int Y ) {
#if 0
    for ( int i = 0; i < trialNodeHeap.size( ); i++ ) {
      if ( trialNodeHeap[ i ].x() == X && trialNodeHeap[ i ].y() == Y ) {
        Node = &trialNodeHeap[ i ];
        return i;
      }
    }
#endif
    int i;
    if ( ( i = indexField.get ( X, Y ) ) != -1 ) {
      Node = &trialNodeHeap[ i ];
      return i;
    }
    TrialNode newNode ( X, Y );
    trialNodeHeap.append ( newNode );
#if 0
    int LastIndex = trialNodeHeap.size( ) - 1;
    Node = & ( trialNodeHeap[ LastIndex ] );
#endif
    int LastIndex = indexField.get ( X, Y );
    Node = & ( trialNodeHeap[ LastIndex ] );
    return LastIndex;
  }

  void marchStep( ) {
    if ( trialNodeHeap.size( ) == 0 ) return;
    TrialNode popNode;
    trialNodeHeap.pop ( popNode );

    makeTrialNodeActive ( popNode );
  }

  RealType getSpeed ( int X, int Y ) {
    return asImp().getSpeed ( X, Y );
  }

  RealType getTime ( const aol::Vec3<short>& Node ) {
    return timeField->get ( Node );
  }

  void downwindFirstOrder ( const aol::Vec2<short> &origNode,
                            const aol::Vec2<short> &newNode,
                            const aol::Vec2<short> &thirdNode,
                            RealType &Value );

  void downwindLeftRight ( int X, int Y, int Step, RealType UpwindVal, RealType &Value );

  void downwindUpDown ( int X, int Y, int Step, RealType UpwindVal, RealType &Value );

  void downwindLeftRight ( int X, int Y, int Step, int UpwindY, RealType &Value );

  void downwindUpDown ( int X, int Y, int Step, int UpwindX, RealType &Value );

  void downwindLeftRight2ndOrder ( int X, int Y, int Step, int UpwindY, RealType &Value );

  void downwindUpDown2ndOrder ( int X, int Y, int Step, int UpwindX, RealType &Value );

public:
  int getGridDepth( ) { return grid.getGridDepth(); }

  qc::IndexedHeap<TrialNode>& getTrialNodeHeap( ) { return trialNodeHeap; }

  qc::ScalarArray<RealType, qc::QC_2D>& getTimeField( ) { return *timeField; }

protected:
  const GridType &grid;
  qc::ScalarArray<RealType, qc::QC_2D> *timeField;
  qc::ScalarArray<int, qc::QC_2D> indexField;
  int order;
  qc::IndexedHeap<TrialNode> trialNodeHeap;

private:
  template< typename DataType >
  inline bool checkB ( const DataType s, const int width ) {
    return ( ( s >= aol::ZTrait<DataType>::zero ) && ( s <= static_cast<DataType> ( width ) ) );
  }

};


template <typename RealType, class EikonalImp, typename GridType>
RealType EikonalNAInt<RealType, EikonalImp, qc::QC_2D, GridType>::solveQuadratic ( RealType a, RealType b, RealType c ) {
  RealType d = b * b - 4.0 * a * c;

  if ( d < 0.0 ) {
    cerr << "ERROR in solveQuadratic: not solvable.\n";
    return -1.0;
  }

  d = sqrt ( d );

  RealType r1, r2;
  r1 = 0.5 * ( -b + d ) / a;
  r2 = 0.5 * ( -b - d ) / a;

  if ( r1 < 0.0  ) return r2;
  else if ( r2 < 0.0 ) return r1;
  else return aol::Max ( r1, r2 );
}

template <typename RealType, class EikonalImp, typename GridType>
void EikonalNAInt<RealType, EikonalImp, qc::QC_2D, GridType>::makeTrialNodeActive ( TrialNode &Node ) {
  const int widthX = grid.getNumX() - 1;
  const int widthY = grid.getNumY() - 1;

  RealType v = timeField->get ( Node.x(), Node.y() );
  if ( v < 0.0f || Node.val < v ) {
    timeField->set ( Node.x(), Node.y(), Node.val );
  }

  // jetzt wird's erst einigermassen spannend:
  // kreiere neue TrialNodes und schiebe sie in den heap,
  // beziehungsweise update der Zeitwerte.
  for ( int dir = 0; dir < 4; dir++ ) {
    TrialNode *newNode = NULL;
    int index = -1;
    aol::Vec2<short> oldNode ( Node.x(), Node.y() );
    aol::Vec2<short> trial1, trial2;

    if ( dir == 0 && checkB ( Node.x() + 1, widthX ) && timeField->get ( Node.x() + 1, Node.y() ) < 0. ) {
      index = searchOrAppendTrialNode ( newNode, Node.x() + 1, Node.y() );

    }
    else if ( dir == 1 && checkB ( Node.x() - 1, widthX ) && timeField->get ( Node.x() - 1, Node.y() ) < 0. ) {
      index = searchOrAppendTrialNode ( newNode, Node.x() - 1, Node.y() );

    }
    else if ( dir == 2 && checkB ( Node.y() + 1, widthY ) && timeField->get ( Node.x(), Node.y() + 1 ) < 0. ) {
      index = searchOrAppendTrialNode ( newNode, Node.x(), Node.y() + 1 );

    }
    else if ( dir == 3 && checkB ( Node.y() - 1, widthY ) && timeField->get ( Node.x(), Node.y() - 1 ) < 0. ) {
      index = searchOrAppendTrialNode ( newNode, Node.x(), Node.y() - 1 );
    }

    if ( index >= 0 ) {
      aol::Vec2<short> newNode2 ( newNode->x(), newNode->y() );
      if ( dir < 2 ) {
        trial1.set ( newNode->x(), Node.y() + 1 );
        trial2.set ( newNode->x(), Node.y() - 1 );
      } else {
        trial1.set ( Node.x() + 1, newNode->y() );
        trial2.set ( Node.x() - 1, newNode->y() );
      }
      if ( checkB ( trial1.x(), widthX ) && checkB ( trial1.y(), widthY ) ) {
        downwindFirstOrder ( oldNode, newNode2, trial1, newNode->val );
      }
      if ( checkB ( trial2.x(), widthX ) && checkB ( trial2.y(), widthY ) ) {
        downwindFirstOrder ( oldNode, newNode2, trial2, newNode->val );
      }
      if ( newNode->val < 0. ) {
        newNode->val = timeField->get ( oldNode ) + grid.H() / getSpeed ( newNode->x(), newNode->y() );

#if 0
        RealType p = 4. / 3.;
        qcLpAnisotropy<RealType> aniso ( 2, p, grid.H() );
        aol::Vec2<RealType> a ( grid.H(), 0. );
        newNode->val = aniso.gammaNorm ( a ) + timeField->get ( oldNode );
#endif

        updateExtensionVelocity ( newNode2, oldNode );
      }
      trialNodeHeap.upHeap ( index );
    }
  }
  return;

  if ( Node.y() + 1 <= grid.getNumY() - 1 && timeField->get ( Node.x(), Node.y() + 1 ) < 0.0f ) {

    TrialNode *newNode;
    int Index = searchOrAppendTrialNode ( newNode, Node.x(), Node.y() + 1 );

    switch ( order ) {
    case 1:
      downwindLeftRight2ndOrder ( newNode->x(), newNode->y(), 1, Node.y(), newNode->val );
      break;
    case 2:
      downwindLeftRight ( newNode->x(), newNode->y(), 1, Node.y(), newNode->val );
      break;
    default:
      downwindLeftRight ( newNode->x(), newNode->y(), 1, Node.val, newNode->val );
    }
    trialNodeHeap.upHeap ( Index );
  }

  if ( Node.y() - 1 >= 0 && timeField->get ( Node.x(), Node.y() - 1 ) < 0.0f ) {
    TrialNode *newNode;
    int Index = searchOrAppendTrialNode ( newNode, Node.x(), Node.y() - 1 );

    switch ( order ) {
    case 1:
      downwindLeftRight2ndOrder ( newNode->x(), newNode->y(), 1, Node.y(), newNode->val );
      break;
    case 2:
      downwindLeftRight ( newNode->x(), newNode->y(), 1, Node.y(), newNode->val );
      break;
    default:
      downwindLeftRight ( newNode->x(), newNode->y(), 1, Node.val, newNode->val );
    }
    trialNodeHeap.upHeap ( Index );
  }

  if ( Node.x() - 1 >= 0 && timeField->get ( Node.x() - 1, Node.y() ) < 0.0f ) {

    TrialNode *newNode;
    int Index = searchOrAppendTrialNode ( newNode, Node.x() - 1, Node.y() );

    switch ( order ) {
    case 1:
      downwindUpDown2ndOrder ( newNode->x(), newNode->y(), 1, Node.x(), newNode->val );
      break;
    case 2:
      downwindUpDown ( newNode->x(), newNode->y(), 1, Node.x(), newNode->val );
      break;
    default:
      downwindUpDown ( newNode->x(), newNode->y(), 1, Node.val, newNode->val );
    }
    trialNodeHeap.upHeap ( Index );
  }

  if ( Node.x() + 1 <= grid.getNumX() - 1 && timeField->get ( Node.x() + 1, Node.y() ) < 0.0f ) {

    TrialNode *newNode;
    int Index = searchOrAppendTrialNode ( newNode, Node.x() + 1, Node.y() );

    switch ( order ) {
    case 1:
      downwindUpDown2ndOrder ( newNode->x(), newNode->y(), 1, Node.x(), newNode->val );
      break;
    case 2:
      downwindUpDown ( newNode->x(), newNode->y(), 1, Node.x(), newNode->val );
      break;
    default:
      downwindUpDown ( newNode->x(), newNode->y(), 1, Node.val, newNode->val );
    }
    trialNodeHeap.upHeap ( Index );
  }
}

template <typename RealType, class EikonalImp, typename GridType>
RealType EikonalNAInt<RealType, EikonalImp, qc::QC_2D, GridType>::solveQuadratic ( RealType v1, RealType v2, RealType h1, RealType h2, RealType F ) {
  /* solve ( T - v1 )^2 / h1^2 + ( T - v2 )^2/h2^2 = 1.0/F^2
     <=> (T^2 - 2*T*v1 + v1^2) * h2^2 + ( T^2 - 2*T*v2 + v2^2) * h1^2 = h1^2 * h2 ^2 / F^2
     <=> (h1^2 + h2^2) * T^2 - 2 * (v1 * h2^2 + v2 * h1^2 ) * T+ h2^2 * v1^2 + h1^2 * v2^2 - h1^2 * h2 ^2 / F^2 = 0 */


  RealType d1 = h1 * h1;
  RealType d2 = h2 * h2;
  RealType a = d1 + d2;
  RealType b = -2.0 * ( v1 * d2 + v2 * d1 );
  RealType c =  d2 * v1 * v1 + d1 * v2 * v2 - d1 * d2 / ( F * F );

  RealType q = b * b - 4.0 * a * c;

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


template <typename RealType, class EikonalImp, typename GridType>
void EikonalNAInt<RealType, EikonalImp, qc::QC_2D, GridType>::downwindFirstOrder ( const aol::Vec2<short> &origNode,
                                                                         const aol::Vec2<short> &newNode,
                                                                         const aol::Vec2<short> &thirdNode,
                                                                         RealType &Value ) {
  RealType F = getSpeed ( newNode.x(), newNode.y() );
  const RealType H = grid.H();
  RealType UpwindVal = timeField->get ( origNode.x(), origNode.y() );

  RealType v = timeField->get ( thirdNode.x(), thirdNode.y() );
  if ( v > 0. ) {
    RealType newV = solveQuadratic ( UpwindVal, v, H, H, F );

#if 0  // only a test
    RealType p = 4. / 3.;
    qcLpAnisotropy<RealType> aniso ( 2, p, H );
    // Newton
    for ( int i = 0; i < 20; i++ ) {
      aol::Vec2<RealType> a ( newV - UpwindVal, newV - v ), b;
      aniso.gammaFirstDerivative ( a, b );
      RealType g = b[0] + b[1];
      newV += ( aniso.gammaNorm ( a ) - H ) / g ;
    }
#endif
    if ( ( Value < 0. || Value > newV ) && newV >= 0. ) {
      Value = newV;
      updateExtensionVelocity ( origNode, newNode, thirdNode, newV );
    }
  }
}

template <typename RealType, class EikonalImp, typename GridType>
void EikonalNAInt<RealType, EikonalImp, qc::QC_2D, GridType>::downwindUpDown ( int X, int Y, int Step, RealType UpwindVal, RealType &Value ) {
  RealType F = getSpeed ( X, Y );
  RealType H = RealType ( Step ) * grid.H();;

  RealType v = timeField->get ( X, Y + Step );
  if ( v > 0.0f ) {
    RealType newH = RealType ( Step ) * grid.H();;
    RealType newV = solveQuadratic ( UpwindVal, v, H, newH, F );
    if ( ( Value < 0.0f || Value > newV ) && newV >= 0.0f ) {
      Value = newV;
    }
  }

  v = timeField->get ( X, Y - Step );
  if ( v > 0.0f ) {
    RealType newH = RealType ( Step ) * grid.H();;
    RealType newV = solveQuadratic ( UpwindVal, v, H, newH, F );
    if ( ( Value < 0.0f || Value > newV ) && newV >= 0.0f ) {
      Value = newV;
    }
  }

  if ( Value <= 0.0f ) {
    Value = UpwindVal + H / F;
  }
}

template <typename RealType, class EikonalImp, typename GridType>
void EikonalNAInt<RealType, EikonalImp, qc::QC_2D, GridType>::downwindLeftRight ( int X, int Y, int Step, RealType UpwindVal, RealType &Value ) {
  RealType F = getSpeed ( X, Y );
  RealType H = RealType ( Step ) * grid.H();;

  RealType v = timeField->get ( X - Step, Y );
  if ( v > 0.0f ) {
    RealType newH = RealType ( Step ) * grid.H();;
    RealType newV = solveQuadratic ( UpwindVal, v, H, newH, F );
    if ( ( Value < 0.0f || Value > newV ) && newV >= 0.0f ) {
      Value = newV;
    }
  }

  v = timeField->get ( X + Step, Y );
  if ( v > 0.0f ) {
    RealType newH = RealType ( Step ) * grid.H();;
    RealType newV = solveQuadratic ( UpwindVal, v, H, newH, F );
    if ( ( Value < 0.0f || Value > newV ) && newV >= 0.0f ) {
      Value = newV;
    }
  }

  if ( Value <= 0.0f ) {
    Value = UpwindVal + H / F;
  }
}

template <typename RealType, class EikonalImp, typename GridType>
void EikonalNAInt<RealType, EikonalImp, qc::QC_2D, GridType>::downwindLeftRight ( int X, int Y, int Step,
                                                                        int UpwindY, RealType &Value ) {

  RealType h = RealType ( Step ) * grid.H();;
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

template <typename RealType, class EikonalImp, typename GridType>
void EikonalNAInt<RealType, EikonalImp, qc::QC_2D, GridType>::downwindUpDown ( int X, int Y, int Step,
                                                                     int UpwindX, RealType &Value ) {

  RealType h = RealType ( Step ) * grid.H();;
  RealType a, b, c, t;

  RealType F = getSpeed ( X, Y );
  RealType UpwindVal = timeField->get ( UpwindX, Y );

  if ( UpwindVal < 0.0 ) cerr << "\n\n\n>>>> ERROR Upwindval < 0\n\n\n";

  RealType v = timeField->get ( X, Y + Step );
  if ( v >= 0.0 ) {
    RealType OrgVal = timeField->get ( UpwindX, Y + Step );

    a = 1.0f;
    b = -2.0 * OrgVal;

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





template <typename RealType, class EikonalImp, typename GridType>
void EikonalNAInt<RealType, EikonalImp, qc::QC_2D, GridType>::downwindLeftRight2ndOrder ( int X, int Y, int Step,
                                                                                int UpwindY, RealType &Value ) {

  RealType h = RealType ( Step ) * grid.H();;
  RealType a, b, c, ax, ay, bx, by;

  int UpwindYY = UpwindY + UpwindY - Y;
  const int widthX = grid.getNumX() - 1;
  const int widthY = grid.getNumY() - 1;

  RealType UpwindY1 = ( checkB ( UpwindY, widthY ) ) ? timeField->get ( X, UpwindY ) : -1.;
  RealType UpwindY2 = ( checkB ( UpwindYY, widthY ) ) ? timeField->get ( X, UpwindYY ) : -1.;
  RealType F = getSpeed ( X, Y );

  if ( UpwindY1 < 0.0 ) cerr << "\n\n\n>>>> ERROR Upwindval < 0\n\n\n";

  if ( UpwindY2 < UpwindY1 && UpwindY2 > 0.0 ) {
    ay = 1.5;
    by = 2.0 * UpwindY1 - 0.5 * UpwindY2;
  } else {
    ay = 1.0;
    by = UpwindY1;
  }

  RealType x1 = ( checkB ( X + Step, widthX ) ) ? timeField->get ( X + Step, Y ) : -1.;
  if ( x1 >= 0.0 ) {
    RealType x2 = ( checkB ( X + Step + Step, widthX ) ) ? timeField->get ( X + Step + Step, Y ) : -1.;
    if (  x2 > 0.0 ) {
      if ( x2 < x1 ) {
        ax = 1.5;
        bx = 2.0 * x1 - 0.5 * x2;
      } else {
        ax = 1.0;
        bx = x1;
      }
      a = ax * ax + ay * ay;
      b = -2.0 * ( ax * bx + ay * by );
      c = bx * bx + by * by - h * h / ( F * F );

      RealType newV = solveQuadratic ( a, b, c );
      if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
        Value = newV;
      }
    }
  }

  x1 = ( checkB ( X - Step, widthX ) ) ? timeField->get ( X - Step, Y ) : -1.;
  if ( x1 >= 0.0f ) {
    RealType x2 = ( checkB ( X - Step - Step, widthX ) ) ? timeField->get ( X - Step - Step, Y ) : -1.;
    if ( x2 > 0.0 ) {
      if ( x2 < x1 ) {
        ax = 1.5;
        bx = 2.0 * x1 - 0.5 * x2;
      } else {
        ax = 1.0;
        bx = x1;
      }
      a = ax * ax + ay * ay;
      b = -2.0 * ( ax * bx + ay * by );
      c = bx * bx + by * by - h * h / ( F * F );

      RealType newV = solveQuadratic ( a, b, c );
      if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
        Value = newV;
      }
    }
  }

  if ( Value < 0.0f ) {
    if ( UpwindY2 < UpwindY1 && UpwindY2 > 0.0 ) {
      Value = ( 2.0 * h / F + 4.0 * UpwindY1  -  UpwindY2 ) / 3.0;
    } else {
      Value = UpwindY1 + h / F;
    }
  }
}

template <typename RealType, class EikonalImp, typename GridType>
void EikonalNAInt<RealType, EikonalImp, qc::QC_2D, GridType>::downwindUpDown2ndOrder ( int X, int Y, int Step,
                                                                             int UpwindX, RealType &Value ) {

  RealType h = RealType ( Step ) * grid.H();;
  RealType a, b, c, ax, bx, ay, by;

  int UpwindXX = UpwindX + UpwindX - X;
  const int widthX = grid.getNumX() - 1;
  const int widthY = grid.getNumY() - 1;

  RealType UpwindX1 = ( checkB ( UpwindX, widthX ) ) ? timeField->get ( UpwindX, Y ) : -1.;
  RealType UpwindX2 = ( checkB ( UpwindXX, widthX ) ) ? timeField->get ( UpwindXX, Y ) : -1.;
  RealType F = getSpeed ( X, Y );

  if ( UpwindX1 < 0.0 ) cerr << "\n\n\n>>>> ERROR Upwindval < 0\n\n\n";

  if ( UpwindX2 < UpwindX1 && UpwindX2 > 0.0 ) {
    ax = 1.5;
    bx = 2.0 * UpwindX1 - 0.5 * UpwindX2;
  } else {
    ax = 1.0;
    bx = UpwindX1;
  }

  RealType y1 = ( checkB ( Y + Step, widthY ) ) ? timeField->get ( X, Y + Step ) : -1.;
  if ( y1 >= 0.0f ) {
    RealType y2 = ( checkB ( Y + Step + Step, widthY ) ) ? timeField->get ( X, Y + Step + Step ) : -1;
    if ( y2 > 0.0 ) {
      if ( y2 < y1 ) {
        ay = 1.5;
        by = 2.0 * y1 - 0.5 * y2;
      } else {
        ay = 1.0;
        by = y1;
      }
      a = ax * ax + ay * ay;
      b = -2.0 * ( ax * bx + ay * by );
      c = bx * bx + by * by - h * h / ( F * F );

      RealType newV = solveQuadratic ( a, b, c );
      if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
        Value = newV;
      }
    }
  }

  y1 = ( checkB ( Y - Step, widthY ) ) ? timeField->get ( X, Y - Step ) : -1.;
  if ( y1 >= 0.0f ) {
    RealType y2 = ( checkB ( Y - Step - Step, widthY ) ) ? timeField->get ( X, Y - Step - Step ) : -1.;
    if ( y2 >= 0. ) {
      if ( y2 < y1 ) {
        ay = 1.5;
        by = 2.0 * y1 - 0.5 * y2;
      } else {
        ay = 1.0;
        by = y1;
      }
      a = ax * ax + ay * ay;
      b = -2.0 * ( ax * bx + ay * by );
      c = bx * bx + by * by - h * h / ( F * F );

      RealType newV = solveQuadratic ( a, b, c );
      if ( newV > 0.0 && ( Value < 0.0f || Value > newV ) ) {
        Value = newV;
      }
    }
  }

  if ( Value < 0.0f ) {
    if ( UpwindX2 < UpwindX1 && UpwindX2 > 0.0 ) {
      Value = ( 2.0 * h / F + 4.0 * UpwindX1  -  UpwindX2 ) / 3.0;
    } else {
      Value = UpwindX1 + h / F;
    }
  }
}


template <typename RealType, typename GridType = qc::GridDefinition>
class EikonalNA : public EikonalNAInt<RealType, EikonalNA<RealType, GridType> , qc::QC_2D, GridType> {
public:
  EikonalNA ( const GridType &Grid )
      : EikonalNAInt<RealType, EikonalNA<RealType, GridType> , qc::QC_2D, GridType> ( Grid ) {}

  RealType getSpeed ( int , int ) {
    return 1.;
  }

  void updateExtensionVelocity ( const aol::Vec2<short> &, const aol::Vec2<short> &,  const aol::Vec2<short> &, RealType ) { /* do nothing, compiler will optimize it away */ }

  void updateExtensionVelocity ( const aol::Vec2<short> &, const aol::Vec2<short> & ) { }

  void newNodeInserted ( const aol::Vec2<short> & ) { }
};

/*template <typename RealType, qc::Dimension Dim>
  class EikonalNAWithExtensionVelocity { };*/

template <typename RealType>
//class EikonalNAWithExtensionVelocity<RealType, qc::QC_2D> : public EikonalNAInt<RealType, EikonalNAWithExtensionVelocity<RealType, qc::QC_2D> , qc::QC_2D> {
class EikonalNAWithExtensionVelocity : public EikonalNAInt<RealType, EikonalNAWithExtensionVelocity<RealType> , qc::QC_2D > {
public:
  EikonalNAWithExtensionVelocity ( const qc::GridDefinition &Grid )
      : EikonalNAInt<RealType, EikonalNAWithExtensionVelocity<RealType> , qc::QC_2D > ( Grid ), _extVelocity ( NULL ) {}

  void setExtVelocityField ( qc::ScalarArray<RealType, qc::QC_2D> &extVelocity ) { _extVelocity = &extVelocity; }

  RealType getSpeed ( int /*X*/, int /*Y*/ ) {
    return 1.;
  }

  void updateExtensionVelocity ( const aol::Vec2<short> &oldNode, const aol::Vec2<short> &newNode,  const aol::Vec2<short> &thirdNode, RealType newT ) {
    if ( _extVelocity == NULL ) {
      throw aol::Exception ( "set velocity field first.", __FILE__, __LINE__ );
    }
    RealType ta = this->timeField->get ( oldNode );
    RealType tb = this->timeField->get ( thirdNode );
    RealType fa = _extVelocity->get ( oldNode );
    RealType fb = _extVelocity->get ( thirdNode );
    _extVelocity->set ( newNode.x(), newNode.y(), ( fa* ( newT - ta ) + fb* ( newT - tb ) ) / ( 2.*newT - ta - tb ) );
  }

  void updateExtensionVelocity ( const aol::Vec2<short> &newNode, const aol::Vec2<short> &oldNode ) {
    if ( _extVelocity == NULL ) {
      throw aol::Exception ( "set velocity field first.", __FILE__, __LINE__ );
    }
    _extVelocity->set ( newNode.x(), newNode.y() , _extVelocity->get ( oldNode.x(), oldNode.y() ) );
  }

  qc::ScalarArray<RealType, qc::QC_2D> *_extVelocity;

  void newNodeInserted ( const aol::Vec2<short> & ) { }
};

template <typename RealType, typename NarrowBandGridType>
class EikonalNANarrowBand : public EikonalNAInt<RealType, EikonalNANarrowBand<RealType, NarrowBandGridType>, qc::QC_2D > {
protected:
  NarrowBandGridType *_narrowBand;
public:
  EikonalNANarrowBand ( const qc::GridDefinition &grid )
      : EikonalNAInt<RealType, EikonalNANarrowBand<RealType, NarrowBandGridType>, qc::QC_2D > ( grid ),
      _narrowBand ( NULL ) {}

  void setNarrowBandGrid ( NarrowBandGridType& narrowBandGrid ) {
    _narrowBand = &narrowBandGrid;
  }

  RealType getSpeed ( int /*X*/, int /*Y*/ ) {
    return 1.;
  }

  void updateExtensionVelocity ( const aol::Vec2<short> &, const aol::Vec2<short> &,  const aol::Vec2<short> &, RealType ) { /* do nothing, compiler will optimize it away */ }

  void updateExtensionVelocity ( const aol::Vec2<short> &, const aol::Vec2<short> & ) { }

  void newNodeInserted ( const aol::Vec2<short> &newNode ) {
    qc::GridDefinition::FullNodeNeighborIterator<qc::QC_2D> eit;
    qc::CoordType c ( newNode[0], newNode[1] );
    for ( eit.start ( this->grid, c ); eit != this->grid.end(); ++eit ) {
      if ( checkElement ( *eit ) ) {
        markElement ( *eit );
      }
    }
  }

  void marchNeighborHood ( RealType Epsilon ) {
    while ( this->trialNodeHeap.size( ) ) {
      typename EikonalNAInt<RealType, EikonalNANarrowBand<RealType, NarrowBandGridType>, qc::QC_2D >::TrialNode popNode;
      this->trialNodeHeap.pop ( popNode );
      this->makeTrialNodeActive ( popNode );

      if ( popNode.val > Epsilon ) {
        return;
      } else {
        newNodeInserted ( aol::Vec2<short> ( popNode[0], popNode[1] ) );
      }
    }
  }

protected:
  bool checkElement ( const qc::Element &el ) {
    for ( qc::ElementNodeIterator<qc::QC_2D> nit = qc::GridDefinition::elementNodeBegin<qc::QC_2D> ( el ); nit != qc::GridDefinition::elementNodeEnd<qc::QC_2D> ( el ); ++nit ) {
      if ( this->timeField->get ( *nit ) < 0 ) return false;
    }
    return true;
  }

  void markElement ( const qc::Element &el ) {
    _narrowBand->insert ( el );
  }
};

}   // end namespace eik

#endif
