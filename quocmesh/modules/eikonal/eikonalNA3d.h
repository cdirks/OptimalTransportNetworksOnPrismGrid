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

#ifndef __EIKONALNA3D_H
#define __EIKONALNA3D_H

#include <eikonalNA.h>

namespace eik {

template <typename RealType, class EikonalImp, typename GridType>
class EikonalNAInt<RealType, EikonalImp , qc::QC_3D, GridType> {
  /*** Barton Nackman ***/
  EikonalImp &asImp( ) { return static_cast<EikonalImp&> ( *this ); }
  const EikonalImp &asImp( ) const { return static_cast<const EikonalImp&> ( *this ); }

  static const qc::CoordType directionLookUp[6];
  static const qc::CoordType trialLookUp[3][4];

public:

class TrialNode : public aol::Vec3<short> {
  public:

    TrialNode ( aol::Vec3<short> Coord )
        : aol::Vec3<short> ( Coord ), val ( aol::NumberTrait<RealType>::NaN )  { // val( -1.0 )  {
    }

    TrialNode ( short X, short Y , short Z )
        : aol::Vec3<short> ( X, Y, Z ), val ( aol::NumberTrait<RealType>::NaN )  {}

    TrialNode( )
        : val ( aol::NumberTrait<RealType>::NaN ) {}

    bool operator== ( const TrialNode &Other ) const {
      return ( testop() == Other.testop() );
    }

    bool operator< ( const TrialNode &Other ) const {
      return ( testop() < Other.testop() );
    }

    bool operator> ( const TrialNode &Other ) const {
      return ( testop() > Other.testop() );
    }

    RealType val;
    RealType testop() const {
      if ( aol::isNaN ( val ) ) return -1.0;
      return aol::Abs ( val );
    }
  };

  //! constructor with grid definition
  EikonalNAInt ( const GridType &Grid ) : grid ( Grid ), timeField ( NULL ), indexField ( Grid ) {
    timeField = new qc::ScalarArray<RealType, qc::QC_3D> ( Grid );
    trialNodeHeap.setIndexField ( &indexField );

    reset();
  }

  ~EikonalNAInt( ) {
    if ( timeField != NULL ) {
      delete timeField;
      timeField = NULL;
    }
  }

  void reset( ) {
    trialNodeHeap.erase ( trialNodeHeap.begin( ), trialNodeHeap.end( ) );
    timeField->setAll ( aol::NumberTrait<RealType>::NaN );
    //timeField->setAll( -1.0f );
    indexField.setAll ( -1 );
  }

  void setOrder ( int Order ) { this->order = Order; }

  //! kann ein punkt so wohl innerhalb als auch auserhalb sein ? (Node->val = copysign( aol::Min( aol::Abs(Node->val) , aol::Abs(Value) ) , Value ); )
  void setSeedPoint ( const aol::Vec3<short> Point, RealType Value = 0. ) {
    TrialNode *Node;
    int index = searchOrAppendTrialNode ( Node, Point.x(), Point.y() , Point.z() );

    if ( !aol::isNaN ( Node->val ) ) {
      Node->val = static_cast<RealType>(copysign ( aol::Min ( aol::Abs ( Node->val ) , aol::Abs ( Value ) ) , Value ));
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
                       trialNodeHeap[ i ].z(),
                       trialNodeHeap[ i ].val );
    }
  }

  void checkTrialNodes( ) {
    for ( int i = 0; i < trialNodeHeap.size( ); i++ ) {
      cerr << trialNodeHeap[i].x() << " " << trialNodeHeap[i].y() << " " << trialNodeHeap[i].z() << " " << trialNodeHeap[i].val << endl;
      /*if ( trialNodeHeap[ i ].val <= 0. ) {
      cerr << "trialNode smaller than zero!\n";
      trialNodeHeap[ i ].val = - trialNodeHeap[ i ].val;
      }*/
    }
  }

  //protected:

  void updateExtensionVelocity ( const aol::Vec3<short> &origNode,
                                 const aol::Vec3<short> &newNode,
                                 const aol::Vec3<short> &thirdNode,
                                 const aol::Vec3<short> &fourthNode,
                                 RealType newValue ) {
    // barton-nackman
    asImp().updateExtensionVelocity ( origNode, newNode, thirdNode, fourthNode, newValue );
  }

  void updateExtensionVelocity ( const aol::Vec3<short> &origNode,
                                 const aol::Vec3<short> &newNode,
                                 const aol::Vec3<short> &thirdNode,
                                 RealType newValue ) {
    // barton-nackman
    asImp().updateExtensionVelocity ( origNode, newNode, thirdNode, newValue );
  }

  void updateExtensionVelocity ( const aol::Vec3<short> &oldNode, const aol::Vec3<short> &newNode , RealType newV ) {
    asImp().updateExtensionVelocity ( oldNode, newNode , newV );
  }

  void makeTrialNodeActive ( TrialNode &Node );

  RealType solveQuadratic ( RealType v1, RealType v2, RealType h1, RealType h2, RealType F = 1.0f );

  RealType solveQuadratic ( RealType v1, RealType v2, RealType v3, RealType h1, RealType h2, RealType h3, RealType F = 1.0f );

  int searchOrAppendTrialNode ( TrialNode *&Node, int X, int Y , int Z ) {
    int i;
    if ( ( i = indexField.get ( X, Y, Z ) ) != -1 ) {
      Node = &trialNodeHeap[ i ];
      return i;
    }
    TrialNode newNode ( X, Y, Z );
    trialNodeHeap.append ( newNode );

    int LastIndex = indexField.get ( X, Y, Z );
    Node = & ( trialNodeHeap[ LastIndex ] );
    return LastIndex;
  }

  void marchStep( ) {
    if ( trialNodeHeap.size( ) == 0 ) return;
    TrialNode popNode;
    trialNodeHeap.pop ( popNode );

    makeTrialNodeActive ( popNode );
  }

  RealType getSpeed ( int X, int Y, int Z )  const {
    return asImp().getSpeed ( X, Y, Z );
  }

  RealType getTime ( const aol::Vec3<short>& Node ) {
    return timeField->get ( Node );
  }

  void downwindFirstOrder ( const aol::Vec3<short> &origNode,
                            const aol::Vec3<short> &newNode,
                            const aol::Vec3<short> &thirdNode,
                            const aol::Vec3<short> &fourthNode,
                            RealType &Value );

public:
  int getGridDepth( ) { return grid.getGridDepth(); }

  qc::IndexedHeap<TrialNode>& getTrialNodeHeap( ) { return trialNodeHeap; }

  qc::ScalarArray<RealType, qc::QC_3D>& getTimeField( ) { return *timeField; }

protected:
  const GridType &grid;
  qc::ScalarArray<RealType, qc::QC_3D> *timeField;
  qc::ScalarArray<int, qc::QC_3D> indexField;
  qc::IndexedHeap<TrialNode> trialNodeHeap;

private:
  template< typename DataType >
  inline bool checkB ( const DataType s, const int width ) {
    return ( ( s >= aol::ZTrait<DataType>::zero ) && ( s <= static_cast<DataType> ( width ) ) );
  }

};

template <typename RealType, class EikonalImp, typename GridType>
const qc::CoordType EikonalNAInt<RealType, EikonalImp , qc::QC_3D, GridType>::directionLookUp[6] = {qc::CoordType ( 1, 0, 0 ), qc::CoordType ( 0, 1, 0 ), qc::CoordType ( 0, 0, 1 ), qc::CoordType ( -1, 0, 0 ), qc::CoordType ( 0, -1, 0 ), qc::CoordType ( 0, 0, -1 ) };

template <typename RealType, class EikonalImp, typename GridType>
const qc::CoordType EikonalNAInt<RealType, EikonalImp , qc::QC_3D, GridType>::trialLookUp[3][4] = { {
      qc::CoordType ( 0, 1, 0 ),
      qc::CoordType ( 0, 0, 1 ),
      qc::CoordType ( 0, -1, 0 ),
      qc::CoordType ( 0, 0, -1 )
    },
    { qc::CoordType ( 1, 0, 0 ),
      qc::CoordType ( 0, 0, 1 ),
      qc::CoordType ( -1, 0, 0 ),
      qc::CoordType ( 0, 0, -1 ) },
    { qc::CoordType ( 1, 0, 0 ),
      qc::CoordType ( 0, 1, 0 ),
      qc::CoordType ( -1, 0, 0 ),
      qc::CoordType ( 0, -1, 0 ) }
                                                                                  };


template <typename RealType, class EikonalImp, typename GridType>
RealType EikonalNAInt<RealType, EikonalImp , qc::QC_3D, GridType>::solveQuadratic ( RealType v1, RealType v2, RealType v3,
                                                                          RealType h1, RealType h2, RealType h3, RealType F ) {
  /* solve ( T - v1 )^2 / h1^2 + ( T - v2 )^2/h2^2+ ( T - v3 )^2/h3^2 = 1.0/F^2 */
  const RealType hh1 = h1 * h1;
  const RealType hh2 = h2 * h2;
  const RealType hh3 = h3 * h3;

  const RealType div = 2.f * F * ( hh1 * hh2 + hh1 * hh3 + hh2 * hh3 );
  const RealType b = 2.f * F * ( hh1 * hh2 * v3 + hh1 * hh3 * v2 + hh2 * hh3 * v1 );

  const RealType a_sqrt = ( 2.f * F * F * hh1 * hh1 * hh2 * v3 * hh3 * v2 ) +
                          ( 2.f * F * F * hh1 * hh2 * hh2 * v3 * hh3 * v1 ) +
                          ( 2.f * F * F * hh1 * hh3 * hh3 * v2 * hh2 * v1 ) +
                          ( hh1 * hh1 * hh2 * hh2 * hh3 ) -
                          ( hh2 * hh3 * hh3 * F * F * hh1 * v2 * v2 ) +
                          ( hh2 * hh2 * hh3 * F * F * hh1 * v3 * v3 ) +
                          ( hh1 * hh1 * hh3 * hh3 * hh2 ) -
                          ( hh1 * hh3 * hh3 * F * F * hh2 * v1 * v1 ) -
                          ( hh1 * hh1 * hh3 * F * F * v3 * v3 ) +
                          ( hh1 * hh1 * hh2 * hh2 * hh3 ) -
                          ( hh1 * hh2 * hh2 * F * F * hh3 * v1 * v1 ) -
                          ( hh1 * hh1 * hh2 * F * F * hh3 * v2 * v2 );

  if ( a_sqrt < 0 ) return aol::NumberTrait<RealType>::NaN; //throw aol::Exception( "negative arg of sqrt.", __FILE__, __LINE__ );

  const RealType r = 2.f * sqrt ( a_sqrt );

  const RealType T1 = ( b + r ) / div;
  const RealType T2 = ( b - r ) / div;
  //    T=( b (+-) 2.*sqrt( a_sqrt ) ) / div;
  if ( T1 < 0 && T2 < 0 )  throw aol::Exception ( "T is neg", __FILE__, __LINE__ );
  if ( T1 < 0 ) return T2;
  else if ( T2 < 0 ) return T1;
  else return aol::Min ( T1, T2 );
}


template <typename RealType, class EikonalImp, typename GridType>
RealType EikonalNAInt<RealType, EikonalImp, qc::QC_3D, GridType>::solveQuadratic ( RealType v1, RealType v2, RealType h1, RealType h2, RealType F ) {
  /* solve ( T - v1 )^2 / h1^2 + ( T - v2 )^2/h2^2 = 1.0/F^2
     <=> (T^2 - 2*T*v1 + v1^2) * h2^2 + ( T^2 - 2*T*v2 + v2^2) * h1^2 = h1^2 * h2 ^2 / F^2
     <=> (h1^2 + h2^2) * T^2 - 2 * (v1 * h2^2 + v2 * h1^2 ) * T+ h2^2 * v1^2 + h1^2 * v2^2 - h1^2 * h2 ^2 / F^2 = 0 */


  RealType d1 = h1 * h1;
  RealType d2 = h2 * h2;
  RealType a = d1 + d2;
  RealType b = -2.0f * ( v1 * d2 + v2 * d1 );
  RealType c =  d2 * v1 * v1 + d1 * v2 * v2 - d1 * d2 / ( F * F );

  RealType q = b * b - 4.0f * a * c;

  if ( q < 0.0f || a == 0.0f ) {
    return aol::NumberTrait<RealType>::NaN;//./0.;//-1.0f;
  }

  RealType r = sqrt ( q );

  if ( r > b ) {
    return ( r - b ) / ( 2.0f * a );
  } else if ( r + b < 0.0f ) {
    return ( - r - b ) / ( 2.0f * a );
  } else {
    return 0.0f;
  }
}


template <typename RealType, class EikonalImp, typename GridType>
void EikonalNAInt<RealType, EikonalImp, qc::QC_3D, GridType>::makeTrialNodeActive ( TrialNode &Node ) {
  const int widthX = grid.getNumX() - 1;
  const int widthY = grid.getNumY() - 1;
  const int widthZ = grid.getNumZ() - 1;

  RealType v = timeField->get ( Node.x(), Node.y(), Node.z() );
  //if ( v < 0.0f || Node.val < v ) {
  if ( aol::isNaN ( v ) || aol::Abs ( Node.val ) < aol::Abs ( v ) ) {
    timeField->set ( Node.x(), Node.y(), Node.z(), Node.val );
  }

  // jetzt wird's erst einigermassen spannend:
  // kreiere neue TrialNodes und schiebe sie in den heap,
  // beziehungsweise update der Zeitwerte.

  for ( int dir = 0; dir < 6; dir++ ) {

    TrialNode *newNode = NULL;
    int index = -1;
    aol::Vec3<short> oldNode ( Node.x(), Node.y(), Node.z() );
    aol::Vec3<short> trial1, trial2;

    const qc::CoordType dirCoord ( Node.x() + directionLookUp[dir].x(), Node.y() + directionLookUp[dir].y(), Node.z() + directionLookUp[dir].z() );

    if ( checkB ( dirCoord.x(), widthX ) && checkB ( dirCoord.y(), widthY ) && checkB ( dirCoord.z(), widthZ ) &&
         aol::isNaN ( timeField->get ( dirCoord.x(), dirCoord.y(), dirCoord.z() ) ) ) {
      index = searchOrAppendTrialNode ( newNode, dirCoord.x(), dirCoord.y(), dirCoord.z() );
    }

    if ( index >= 0 ) {
      aol::Vec3<short> newNode2 ( newNode->x(), newNode->y(), newNode->z() );
      for ( int trial = 0 ; trial != 4 ; trial++ ) {

        trial1.set ( ( newNode->x() + trialLookUp[dir%3][trial].x() ),
                     ( newNode->y() + trialLookUp[dir%3][trial].y() ),
                     ( newNode->z() + trialLookUp[dir%3][trial].z() ) );

        trial2.set ( ( newNode->x() + trialLookUp[dir%3][ ( trial+1 ) %4].x() ),
                     ( newNode->y() + trialLookUp[dir%3][ ( trial+1 ) %4].y() ),
                     ( newNode->z() + trialLookUp[dir%3][ ( trial+1 ) %4].z() ) );


        if ( checkB ( trial1.x(), widthX ) && checkB ( trial1.y(), widthY ) && checkB ( trial1.z(), widthZ ) &&
             checkB ( trial2.x(), widthX ) && checkB ( trial2.y(), widthY ) && checkB ( trial2.z(), widthZ ) ) {
          downwindFirstOrder ( oldNode, newNode2, trial1, trial2, newNode->val );
        }
      }
      //if ( newNode->val < 0. ) {
      if ( aol::isNaN ( newNode->val ) ) {
        newNode->val = static_cast<RealType>(aol::Abs ( timeField->get ( oldNode ) ) + grid.H() / getSpeed ( newNode->x(), newNode->y() , newNode->z() ));
        newNode->val = static_cast<RealType>(copysign ( newNode->val , timeField->get ( oldNode ) ));
        updateExtensionVelocity ( oldNode, newNode2 , newNode->val );
      }
      trialNodeHeap.upHeap ( index );
    }
  }
}

template <typename RealType, class EikonalImp, typename GridType>
void EikonalNAInt<RealType, EikonalImp, qc::QC_3D, GridType>::downwindFirstOrder ( const aol::Vec3<short> &origNode,
                                                                         const aol::Vec3<short> &newNode,
                                                                         const aol::Vec3<short> &thirdNode,
                                                                         const aol::Vec3<short> &fourthNode,
                                                                         RealType &Value ) {
  RealType F = getSpeed ( newNode.x(), newNode.y() , newNode.z() );
  const RealType H = static_cast<RealType>(grid.H());
  RealType UpwindVal = aol::Abs ( timeField->get ( origNode.x(), origNode.y(), origNode.z() ) );
  //RealType UpwindVal = ( timeField->get( origNode.x(), origNode.y(), origNode.z() ) );

  RealType v1 = aol::Abs ( timeField->get ( thirdNode.x(), thirdNode.y(), thirdNode.z() ) );
  //RealType v1 = timeField->get( thirdNode.x(), thirdNode.y(), thirdNode.z() ) ;
  RealType v2 = aol::Abs ( timeField->get ( fourthNode.x(), fourthNode.y(), fourthNode.z() ) );

  //if( v1 > 0. && v2 > 0. ){
  if ( !aol::isNaN ( v1 ) && !aol::isNaN ( v2 ) && v1 != 0 && v2 != 0 ) {
    RealType newV = solveQuadratic ( UpwindVal, v1, v2, H, H, H, F );
    //if (( Value < 0. || Value > newV ) && newV >= 0. ){
    if ( ( aol::isNaN ( Value ) || aol::Abs ( Value ) > newV ) && !aol::isNaN ( newV ) ) {
      Value = static_cast<RealType>(copysign ( newV ,  Value ));
      updateExtensionVelocity ( origNode, newNode, thirdNode, fourthNode, static_cast<RealType>(copysign ( newV ,  Value )) );
    }
  }
  else {
    //if( v1 > 0.){
    if ( !aol::isNaN ( v1 ) && v1 != 0 ) {
      RealType newV = solveQuadratic ( UpwindVal, v1, H, H, F );
      //if (( Value < 0. || Value > newV ) && newV >= 0. ){
      if ( ( aol::isNaN ( Value ) || aol::Abs ( Value ) > newV ) && !aol::isNaN ( newV ) ) {
        Value = static_cast<RealType>(copysign ( newV ,  Value ));
        updateExtensionVelocity ( origNode, newNode, thirdNode, static_cast<RealType>(copysign ( newV ,  Value )) );
      }
      //}else if( v2 > 0.){
    }
    else if ( !aol::isNaN ( v2 ) && v2 != 0 ) {
      RealType newV = solveQuadratic ( UpwindVal, v2, H, H, F );
      //if (( Value < 0. || Value > newV ) && newV >= 0. ){
      if ( ( aol::isNaN ( Value ) || aol::Abs ( Value ) > newV ) && !aol::isNaN ( newV ) ) {
        Value = static_cast<RealType>(copysign ( newV ,  Value ));
        updateExtensionVelocity ( origNode, newNode, fourthNode,  static_cast<RealType>(copysign ( newV ,  Value )) );
      }
    }
  }
}


template <typename RealType, typename GridType = qc::GridDefinition>
class EikonalNA3D : public EikonalNAInt<RealType, EikonalNA3D<RealType, GridType> , qc::QC_3D, GridType> {
public:
 // ! constructor with grid definition
  EikonalNA3D ( const GridType &Grid )
      : EikonalNAInt<RealType, EikonalNA3D<RealType, GridType> , qc::QC_3D, GridType> ( Grid ) {}

  RealType getSpeed ( int /*X*/ , int /*Y*/, int /*Z*/ ) const {
    return 1.;
  }

  void updateExtensionVelocity ( const aol::Vec3<short> &, const aol::Vec3<short> &,  const aol::Vec3<short> &, const aol::Vec3<short> &, RealType ) { /* do nothing, compiler will optimize it away */ }

  void updateExtensionVelocity ( const aol::Vec3<short> &, const aol::Vec3<short> &,  const aol::Vec3<short> &, RealType ) { /* do nothing, compiler will optimize it away */ }

  void updateExtensionVelocity ( const aol::Vec3<short> &, const aol::Vec3<short> & , RealType ) { }
};



template <typename RealType>
class EikonalNA3DWithExtVelocity : public EikonalNAInt<RealType, EikonalNA3DWithExtVelocity<RealType> , qc::QC_3D > {
public:
 // ! constructor with grid definition
  EikonalNA3DWithExtVelocity ( const qc::GridDefinition &Grid )
      : EikonalNAInt<RealType, EikonalNA3DWithExtVelocity<RealType> , qc::QC_3D > ( Grid ) {}

   void setExtVelocityField( qc::ScalarArray<RealType, qc::QC_3D> &extVelocity ) { _extVelocity = &extVelocity; }


  void updateExtensionVelocity ( const aol::Vec3<short> &origNode, const aol::Vec3<short> &newNode, const aol::Vec3<short> &thirdNode, const aol::Vec3<short> &fourthNode, RealType newV){

    if (_extVelocity == NULL ) {
      throw aol::Exception( "set velocity field first. \n ", __FILE__, __LINE__);
    }

  RealType ta = this->timeField->get ( origNode );
  RealType tb = this->timeField->get ( thirdNode );
  RealType tc = this->timeField->get ( fourthNode );

  RealType fa = _extVelocity->get ( origNode );
  RealType fb = _extVelocity->get ( thirdNode );
  RealType fc = _extVelocity->get ( fourthNode );

  _extVelocity->set( newNode.x(), newNode.y(), newNode.z(),  (fa* (ta- newV) + fb* (tb- newV) + fc* (tc- newV)) / (ta + tb + tc - (3. * newV))  );

  }

 void updateExtensionVelocity ( const aol::Vec3<short> &origNode, const aol::Vec3<short> &newNode, const aol::Vec3<short> &thirdNode, RealType newV){

    if (_extVelocity == NULL ) {
      throw aol::Exception( "set velocity field first. \n ", __FILE__, __LINE__);
    }

  RealType ta = this->timeField->get ( origNode );
  RealType tb = this->timeField->get ( thirdNode );
  //RealType tc = timeField->get ( fourthNode );

  RealType fa = _extVelocity->get ( origNode );
  RealType fb = _extVelocity->get ( thirdNode );
  //RealType fc = _extVelocity->get ( fourthNode );

  _extVelocity->set( newNode.x(), newNode.y(), newNode.z(),  (fa* (ta- newV) + fb* (tb- newV) ) / (ta + tb  -(2. * newV))  );

  }

 void updateExtensionVelocity ( const aol::Vec3<short> &origNode, const aol::Vec3<short> &newNode, RealType /*newV*/){

    if (_extVelocity == NULL ) {
      throw aol::Exception( "set velocity field first. \n ", __FILE__, __LINE__);
    }

   //RealType ta = timeField->get ( origNode );
  //RealType tb = timeField->get ( thirdNode );
  //RealType tc = timeField->get ( fourthNode );

  RealType fa = _extVelocity->get ( origNode );
  //RealType fb = _extVelocity->get ( thirdNode );
  //RealType fc = _extVelocity->get ( fourthNode );

  _extVelocity->set( newNode.x(), newNode.y(), newNode.z(), fa );

  }


  qc::ScalarArray<RealType, qc::QC_3D> *_extVelocity;


};

/*

template<typename RealType>
class EikonalNA3DWithExtVelocity : public EikonalNAInt<RealType, EikonalNA3DWithExtVelocity< RealType>, qc::QC_3D> {
public:


  EikonalNA3DWithExtVelocity ( const qc::GridDefinition &Grid )
 //   : EikonalNAInt<RealType, EikonalNA3DWithExtVelocity<RealType>, qc::QC_3D> (Grid){}

  void setExtVelocityField( qc::ScalarArray<RealType, qc::QC_3D> &extVelocity ) { _extVelocity = &extVelocity; }


  void updateExtensionVelocity ( const aol::Vec3<short> &origNode, const aol::Vec3<short> &newNode, const aol::Vec3<short> &thirdNode, const aol::Vec3<short> &fourthNode, RealType newV){

    if (_extVelocity == NULL ) {
      throw aol::Exception( "set velocity field first. \n ", __FILE__, __LINE__);
    }

  RealType ta = timeField->get ( origNode );
  RealType tb = timeField->get ( thirdNode );
  RealType tc = timeField->get ( fourthNode );

  RealType fa = _extVelocity->get ( origNode );
  RealType fb = _extVelocity->get ( thirdNode );
  RealType fc = _extVelocity->get ( fourthNode );

  _extVelocity->set( newNode.x(), newNode.y(), newNode.z(),  (fa* (ta- newV) + fb* (tb- newV) + fc* (tc- newV)) / (ta + tb + tc - (3. * newV))  );

  }


 void updateExtensionVelocity ( const aol::Vec3<short> &origNode, const aol::Vec3<short> &newNode, const aol::Vec3<short> &thirdNode, RealType newV){

    if (_extVelocity == NULL ) {
      throw aol::Exception( "set velocity field first. \n ", __FILE__, __LINE__);
    }

  RealType ta = timeField->get ( origNode );
  RealType tb = timeField->get ( thirdNode );
  //RealType tc = timeField->get ( fourthNode );

  RealType fa = _extVelocity->get ( origNode );
  RealType fb = _extVelocity->get ( thirdNode );
  //RealType fc = _extVelocity->get ( fourthNode );

  _extVelocity->set( newNode.x(), newNode.y(), newNode.z(),  (fa* (ta- newV) + fb* (tb- newV) ) / (ta + tb  -(2. * newV))  );

  }

 void updateExtensionVelocity ( const aol::Vec3<short> &origNode, const aol::Vec3<short> &newNode, RealType newV){

    if (_extVelocity == NULL ) {
      throw aol::Exception( "set velocity field first. \n ", __FILE__, __LINE__);
    }

   //RealType ta = timeField->get ( origNode );
  //RealType tb = timeField->get ( thirdNode );
  //RealType tc = timeField->get ( fourthNode );

  RealType fa = _extVelocity->get ( origNode );
  //RealType fb = _extVelocity->get ( thirdNode );
  //RealType fc = _extVelocity->get ( fourthNode );

  _extVelocity->set( newNode.x(), newNode.y(), newNode.z(), fa );

  }

  qc::ScalarArray<RealType, qc::QC_3D> *_extVelocity;
};


*/



template <typename RealType, typename NarrowBandGridType>
class EikonalNANarrowBand3D : public EikonalNAInt<RealType, EikonalNANarrowBand3D<RealType, NarrowBandGridType>, qc::QC_3D > {
protected:
  NarrowBandGridType *_narrowBand;
public:
  EikonalNANarrowBand3D ( const qc::GridDefinition &grid )
      : EikonalNAInt<RealType, EikonalNANarrowBand3D<RealType, NarrowBandGridType>, qc::QC_3D > ( grid ),
      _narrowBand ( NULL ) {}

  void setNarrowBandGrid ( NarrowBandGridType& narrowBandGrid ) {
    _narrowBand = &narrowBandGrid;
  }

  RealType getSpeed ( int /*X*/, int /*Y*/, int /*Z*/ ) const {
    return 1.0;
  }

  void updateExtensionVelocity ( const aol::Vec3<short> &, const aol::Vec3<short> &,  const aol::Vec3<short> &, const aol::Vec3<short> &, RealType ) {}

  void updateExtensionVelocity ( const aol::Vec3<short> &, const aol::Vec3<short> &,  const aol::Vec3<short> &, RealType ) {}

  void updateExtensionVelocity ( const aol::Vec3<short> &, const aol::Vec3<short> &, RealType ) {}

  void newNodeInserted ( const aol::Vec3<short> &newNode ) {
    // cerr << "inserted node: " << newNode << endl;
    qc::GridDefinition::FullNodeNeighborIterator<qc::QC_3D> eit;
    qc::CoordType c ( newNode[0], newNode[1], newNode[2] );
    for ( eit.start ( this->grid, c ); eit != this->grid.end(); ++eit ) {
      if ( checkElement ( *eit ) ) {
        markElement ( *eit );
      }

    }
  }

  void marchNeighborHood ( RealType Epsilon ) {
    while ( this->trialNodeHeap.size( ) ) {
      typename EikonalNAInt<RealType, EikonalNANarrowBand3D<RealType, NarrowBandGridType>, qc::QC_3D >::TrialNode popNode;
      //TrialNode popNode;
      this->trialNodeHeap.pop ( popNode );
      this->makeTrialNodeActive ( popNode );

      if ( popNode.val - 1.0 > Epsilon ) {
        return;
      } else {
        newNodeInserted ( aol::Vec3<short> ( popNode[0], popNode[1], popNode[2] ) );
      }
    }
  }

protected:
  bool checkElement ( const qc::Element &el ) {
    qc::ElementNodeIterator<qc::QC_3D> nit =  qc::GridDefinition::elementNodeBegin<qc::QC_3D> ( el );
    for ( ; nit != qc::GridDefinition::elementNodeEnd<qc::QC_3D> ( el ); ++nit ) {
      if ( this->timeField->get ( *nit ) < 0 ) return false;
    }
    return true;
  }

  void markElement ( const qc::Element &el ) {
    // cerr << " inserting " << el;
    if ( el.x() >= 0 && el.y() >= 0 && el.z() >= 0 ) {
      _narrowBand->insert ( el );
    }
  }
};

#if 0
#include <narrowBandGrid.h>

//! eikonal solver for NarrowBand \author Notthoff
template <typename RealType>
class EikonalNB : public EikonalNAInt<RealType, EikonalNB<RealType> , qc::QC_3D > {
  const RealType &eps;
  narr::qcNarrowBandGrid<RealType> &_narGrid;
public:
  EikonalNB ( narr::qcNarrowBandGrid<RealType> &Grid , const RealType &Eps )
      : EikonalNAInt<RealType, EikonalNB<RealType> , qc::QC_3D > ( Grid ) , eps ( Eps ), _narGrid ( Grid ) {}


  RealType getSpeed ( int X, int Y, int Z ) const {
    return 1.0;
  }

  /* used to update the hash_set of the narrowband and the hash_map which store the distance funktion in the narrow band */
  void updateExtensionVelocity ( const aol::Vec3<short> &origNode, const aol::Vec3<short> &newNode,  const aol::Vec3<short> &thirdNode, const aol::Vec3<short> &fourthNode, RealType newV ) {
    this->updateExtensionVelocity ( origNode, newNode, newV );
  }

  void updateExtensionVelocity ( const aol::Vec3<short> &origNode, const aol::Vec3<short> &newNode,  const aol::Vec3<short> &thirdNode, RealType newV ) {
    this->updateExtensionVelocity ( origNode, newNode, newV );

  }

  void updateExtensionVelocity ( const aol::Vec3<short> &origNode, const aol::Vec3<short> &newNode , RealType newV ) {
    if ( aol::Abs ( timeField->get ( origNode ) ) < eps ) {
      qc::Element tmp_el ( origNode.x(), origNode.y(), origNode.z(), _narGrid.getGridDepth(), 9 );
      _narGrid.getNarrowSet()->insert ( tmp_el );
      ( *_narGrid.getNarrowMap() ) [origNode] = timeField->get ( origNode );
      if ( aol::Abs ( newV ) > eps  && !aol::isNaN ( newV ) ) {
        ( *_narGrid.getNarrowBoundMap() ) [newNode] = newV;
        ( *_narGrid.getNarrowBoundMap() ).erase ( origNode );
      }
    }
  }
};
#endif

}   // end namespace eik

#endif

