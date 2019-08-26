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

#ifndef __EIKONAL3D_H
#define __EIKONAL3D_H

#include <aol.h>

#include <ocTree.h>
#include <estimator3d.h>
#include <scalarArray.h>
#include <bitVector.h>
#include <qmHeap.h>

#include <segDefs.h>
#include <segTrialNode.h>

// originates from file SegEstimator3d.h:

namespace eik {

const int EST_GREY_UP   = 1;
const int EST_GREY_DOWN = 2;
const int EST_GRAD      = 3;

class SegEstimator3d : public qc::Estimator3d<unsigned char> {
public:
  SegEstimator3d ( int Size, int ID ) :
      qc::Estimator3d<unsigned char> ( Size ),
      id ( ID ) {}

  void         clear( );

  void         setThreshold ( unsigned char Thres ) { Threshold = Thres; }

  unsigned char         getThreshold( ) { return Threshold; }

  int          checkElement ( int X, int Y, int Z, int Level, unsigned char Threshold ) {
    return qc::Estimator3d<unsigned char>::checkElement ( X, Y, Z, Level, Threshold );
  }

  int          checkElement ( int X, int Y, int Z, int Level ) {
    return qc::Estimator3d<unsigned char>::checkElement ( X, Y, Z, Level, Threshold );
  }

  int          checkElement ( const qc::Element &El ) {
    return qc::Estimator3d<unsigned char>::checkElement ( El, Threshold );
  }

  void         setImagePointer ( qc::ScalarArray<float, qc::QC_3D> *OrigImage ) {
    origImage = OrigImage;
  }

  void         retrieveOriginalImage( );

  void makeGradients( );

  int  getID( ) const { return id; }

protected:
  const int id;
  unsigned char  Threshold;
  qc::ScalarArray<float, qc::QC_3D> *origImage; // initialization??
};


// originates from file MetaEstimator.h:
const int  MAX_NUM_ESTIMATORS = 10;

class Traversable {
  virtual ~Traversable() {}
public:
  virtual void actionOnElement ( int X, int Y, int Z, int Level ) = 0;
};

class MetaEstimator {

protected:
  struct SegEstItem {
    SegEstimator3d *Estimator;
    unsigned char Threshold;
    int id;
  };

public:
  MetaEstimator ( int Size ) :
      size ( Size ),
      numEstimators ( 0 ) {
    maxLevel = -1;
  }

  ~MetaEstimator() {
    cerr << "MetaEstimator::~MetaEstimator (destructor) called, doing nothing" << endl;
  }

  void addEstimator ( SegEstimator3d *Estimator ) {
    estimatorList[ numEstimators   ].Estimator = Estimator;
    estimatorList[ numEstimators++ ].id = Estimator->getID( );

    if ( maxLevel != -1 && maxLevel != Estimator->getMaxLevel( ) ) {
      cerr << "WARNING: Putting incompatible Estimators into MetaEstimator.\n";
    } else {
      maxLevel = Estimator->getMaxLevel( );
      cerr << "- adding an estimator to metaestimator. " << endl;
    }
  }

  void setMaxLevel ( int MaxLevel ) { maxLevel = MaxLevel; }

  void traverseAllElements ( Traversable& Trav ) {
    recursiveTraverse ( 0, 0, 0, 0, Trav );
  }

  void setThreshold ( int ID, unsigned char Threshold ) {
    for ( int i = 0; i < numEstimators; i++ ) {
      if ( estimatorList[ i ].id == ID ) {
        estimatorList[ i ].Threshold = Threshold;
        return;
      }
    }
    /*cerr << "- WARNING in setThreshold: Estimator with ID " << ID << " not found.\n";*/
  }

  int  checkElementSure ( int Dir, int X, int Y, int Z, int Level ) {
    if ( Level == maxLevel ) return 1;
    else return checkElement ( Dir, X, Y, Z, Level );
  }

  int  checkElementSure ( int X, int Y, int Z, int Level ) {
    if ( Level == maxLevel ) return 1;
    else
      return checkElement ( X, Y, Z, Level );
  }

  int  checkElement ( int Dir, int X, int Y, int Z, int Level ) {

#ifdef NOT_ADAPTIVE
    return 0;
#endif

    int step = 1 << ( maxLevel - Level );
    int checkX = X,
                 checkY = Y,
                          checkZ = Z;

    if ( ! ( Dir & TN_DIR_RIGHT ) )
      checkX -= step;

    if ( ! ( Dir & TN_DIR_UP ) )
      checkY -= step;

    if ( ! ( Dir & TN_DIR_BACK ) )
      checkZ -= step;

    return checkElement ( checkX, checkY, checkZ, Level );
  }

  int  checkElement ( int X, int Y, int Z, int Level ) {

#ifdef NOT_ADAPTIVE
    return 0;
#endif

    if ( numEstimators == 0 ) {
      return ( Level == maxLevel );
    }

    for ( int i = 0; i < numEstimators; i++ ) {
      if ( estimatorList[ i ].Estimator->checkElement ( X, Y, Z, Level,
                                                        estimatorList[ i ].Threshold ) == 0 ) {
        return 0;
      }
    }
    return 1;
  }

  int getSize( ) { return size; }

protected:

  void recursiveTraverse ( int X, int Y, int Z, int Level, Traversable &Trav ) {
    if ( Level == maxLevel || checkElement ( X, Y, Z, Level ) ) {
      Trav.actionOnElement ( X, Y, Z, Level );
    } else {
      int hStep = 1 << ( maxLevel - ( Level + 1 ) );
      recursiveTraverse ( X        , Y        , Z        , Level + 1, Trav );
      recursiveTraverse ( X + hStep, Y        , Z        , Level + 1, Trav );
      recursiveTraverse ( X        , Y + hStep, Z        , Level + 1, Trav );
      recursiveTraverse ( X + hStep, Y + hStep, Z        , Level + 1, Trav );

      recursiveTraverse ( X        , Y        , Z + hStep, Level + 1, Trav );
      recursiveTraverse ( X + hStep, Y        , Z + hStep, Level + 1, Trav );
      recursiveTraverse ( X        , Y + hStep, Z + hStep, Level + 1, Trav );
      recursiveTraverse ( X + hStep, Y + hStep, Z + hStep, Level + 1, Trav );
    }
  }

  SegEstItem estimatorList[ MAX_NUM_ESTIMATORS ];
  int size;
  int numEstimators;
  int maxLevel;
};


const int NOF     = 0;
const int UPF     = 1;
const int DOWNF   = 2;
const int LEFTF   = 4;
const int RIGHTF  = 8;
const int FRONTF  = 16;
const int BACKF   = 32;
const int ALLF    = 63;

// originates from file Segment3d.h:

#if 0
// this seems not to be used anyway.
inline void swap_bytes_4 ( void *v ) {
  char a, b, c, d;
  a = static_cast<char*> ( v ) [0]; // does this work or should it be a reinterpret cast?
  b = static_cast<char*> ( v ) [1];
  c = static_cast<char*> ( v ) [2];
  d = static_cast<char*> ( v ) [3];

  static_cast<char*> ( v ) [0] = d;
  static_cast<char*> ( v ) [1] = c;
  static_cast<char*> ( v ) [2] = b;
  static_cast<char*> ( v ) [3] = a;
}
#endif


const bool FMM_WARNINGS = false;

typedef TrialNode<float> TrialNodeFloat;
typedef TrialNode<double> TrialNodeDouble;

class InitInteriorTimeFieldTraversable;
class InitBoundaryTimeFieldTraversable;

inline int iabs ( int a ) { return ( a > 0 ) ? a : -a; }
inline int imax ( int a, int b ) { return ( a > b ) ? a : b; }
inline int imin ( int a, int b ) { return ( a < b ) ? a : b; }

class Eikonal3d {
  friend class InitInteriorTimeFieldTraversable;
  friend class InitBoundaryTimeFieldTraversable;

public:

  // lookuptables:
  int tn_nb_dirs[ 8 ][ 3 ];   // this array stores for each direction the 3 "neighbourdirections"
  int tn_hn_offs[ 18 ][ 3 ];   // this array stores for a hanging node direction the offset from the lowleftfront position of the element
  int tn_hn_nb_cell_dirs[ 18 ];
  int tn_dir_ordering[ 8 ];   // this array gives an ordering of the 8 directions such that they are "neighboured"

  Eikonal3d( );

  virtual ~Eikonal3d();

  void setMetaEstimator ( MetaEstimator *Estimator ) {
    estimator = Estimator;
  }

  int  getGridDepth( ) { return GRID_DEPTH; }

  void init ( int NumX, int NumY, int NumZ );

  void refine ( int X, int Y, int Z, int level );

  void recMarkTouchedCell ( int X, int Y, int Z, int level );

  void reset( );

  qc::ScalarArray<unsigned char, qc::QC_3D>* getComplementField( ) { return complementField; }

  void initComplementField( );

  void addSeedPoint ( int X, int Y, int Z );

  void initSeedElement ( TrialNode<float> &Node );

  void makeTrialNodeActive ( TrialNode<float> &Node );

  void localHangNodeGridUpdate ( TrialNode<float> &Node );

  static inline float solveQuadratic ( float a, float b, float c );

  static inline float solveQuadratic ( float v1, float v2, float h1, float h2, float F = 1.0 );

  static inline float solveQuadratic ( float v1, float v2, float v3,
                                       float h1, float h2, float h3, float F = 1.0 );

  inline void updateElementToTrialNodes ( int Dir, int X, int Y, int Z, int Level );

  void updateElementToTrialNodes ( int X, int Y, int Z, int Level );

  inline int searchOrAppendTrialNode ( TrialNode<float> *&Node, int X, int Y, int Z );

  void segmentStep( ) {
    if ( trialNodeHeap.size( ) == 0 ) return;
    TrialNode<float> popNode;
    trialNodeHeap.pop ( popNode );
    makeTrialNodeActive ( popNode );
  }

  inline void createElement ( int Direction, TrialNode<float> &Node );

  void dumpTrialNodes( );

  void dumpTrialNode ( int X, int Y, int Z );

  void dumpMinTrialNode( );

  virtual float getSpeed ( int X, int Y, int Z , int Filter = 0 ) = 0;
  // a suitable getSpeed function should be implemented in classes derived from Eikonal3d.
  //  return 1.0f;
  // }

void setScaleTime ( float Scale ) { scaleTime = Scale; }

  float getTime ( int X, int Y, int Z ) { return timeField->get ( X, Y, Z ); }

  void  setTime ( int X, int Y, int Z, float T ) { timeField->set ( X, Y, Z, T ); }

  float peekMinTime( ) {
    if ( trialNodeHeap.size( ) == 0 ) return 0.0f;

    TrialNodeFloat node;
    trialNodeHeap.peekMin ( node );
    return node.val;
  }

  void trialNodesConsistencyCheck( );

  inline void downwind ( int Dir, int X, int Y, int Z, int Step, float &Value, float MeanCurv = 0.0f );

  void exportBitMask ( qc::ScalarArray<unsigned char, qc::QC_3D> &Array, float Time );

  void exportToArray3d ( qc::ScalarArray<float, qc::QC_3D> *Estimator );

  void exportToArray3d ( qc::ScalarArray<unsigned char, qc::QC_3D> *Estimator );

  void  generateTimeSliceX ( qc::ScalarArray<float, qc::QC_2D> *TimeSlice, int X );

  void  generateTimeSliceY ( qc::ScalarArray<float, qc::QC_2D> *TimeSlice, int Y );

  void  generateTimeSliceZ ( qc::ScalarArray<float, qc::QC_2D> *TimeSlice, int Z );

  inline int checkElementDataSet ( int X, int Y, int Z, int Level );

  inline int checkElementSomeDataSet ( int X, int Y, int Z, int Level );

  float getRelTimeArea ( int X, int Y, int Z, int Level );

qc::ScalarArray<float, qc::QC_3D>* getTimeField( ) { return timeField; }

  void loadTimeField ( char *FileName ) {
    timeField->load ( FileName );
  }

  void setTimeField ( float *Data ) {
    cerr << "max = " << timeField->getMaxValue( );
    cerr << "min = " << timeField->getMinValue( );
    timeField->readFromBuffer ( Data );

    cerr << "max = " << timeField->getMaxValue( );
    cerr << "min = " << timeField->getMinValue( );

    ocTree->clearAll( );
    recMarkTouchedCell ( 0, 0, 0, 0 );
  }

  void cont_dist( );

  void accept ( float T );

  void erase ( float T );

  qc::ScalarArray<unsigned char, qc::QC_3D> *getTagField( ) { return tagField; }

  qc::ScalarArray<unsigned char, qc::QC_3D> *getSeedTypeField( ) { return seedTypeField; }



protected:

  void fillTimeField( );
  void createLookUpTables( );

  /* MEMBER VARIABLES */
  float scaleTime;
  float thres;
  int numX, numY, numZ;
  int GRID_DEPTH;

  qc::OcTree *ocTree;
  MetaEstimator  *estimator;
  qc::ScalarArray<float, qc::QC_3D> *timeField;
  qc::ScalarArray<float, qc::QC_3D> *finalTimeField;
  qc::ScalarArray<float, qc::QC_3D> *tmpTimeField;
  qc::ScalarArray<unsigned char, qc::QC_3D>  *tagField;
  qc::ScalarArray<int, qc::QC_3D>   *indexField;
  qc::ScalarArray<unsigned char, qc::QC_3D>  *complementField;
  qc::ScalarArray<unsigned char, qc::QC_3D>  *seedTypeField;
  qc::IndexedHeap<TrialNodeFloat> trialNodeHeap;
  vector<TrialNodeFloat> hangNodes;
  vector<TrialNodeFloat> interpolatedHangNodes;

  // changed slist to list. does it still work?
  // slist<TrialNodeFloat> undoList;
  list<TrialNodeFloat> undoList;

  deque<TrialNodeFloat> actList;

  //  aol::Vec3<short> firstSeedPoint; // necessary??
};


void Eikonal3d::updateElementToTrialNodes ( int Dir, int X, int Y, int Z, int Level ) {
  int Step = 1 << ( GRID_DEPTH - 1 - Level );
  int NewX = X, NewY = Y, NewZ = Z;
  if ( ! ( Dir & TN_DIR_RIGHT ) ) {
    NewX -= Step;
  }
  if ( ! ( Dir & TN_DIR_UP ) ) {
    NewY -= Step;
  }
  if ( ! ( Dir & TN_DIR_BACK ) ) {
    NewZ -= Step;
  }
  updateElementToTrialNodes ( NewX, NewY, NewZ, Level );
}


int Eikonal3d::searchOrAppendTrialNode ( TrialNode<float> *&Node,
                                         int X, int Y, int Z ) {
  int i = indexField->get ( X, Y, Z );
  if ( i != -1 ) {
    Node = &trialNodeHeap[ i ];
    return i;
  }

  TrialNodeFloat newNode ( X, Y, Z );
  trialNodeHeap.append ( newNode );
  int LastIndex = trialNodeHeap.size( ) - 1;
  Node = & ( trialNodeHeap[ LastIndex ] );
  return LastIndex;
}


float Eikonal3d::solveQuadratic ( float a, float b, float c ) {
  float d = b * b - 4 * a * c;

  if ( d < 0.0f ) {
    cerr << "ERROR in solveQuadratic: not solvable\n";
    return -1.0;
  }
  return ( -b + sqrt ( d ) ) / ( 2.0f * a );
}

float Eikonal3d::solveQuadratic ( float v1, float v2, float h1, float h2, float F ) {
  /* solve ( T - v1 )^2 / h1^2 + ( T - v2 )^2/h2^2 = 1.0/F^2
     <=> (T^2 - 2*T*v1 + v1^2) * h2^2 + ( T^2 - 2*T*v2 + v2^2) * h1^2 = h1^2 * h2 ^2 / F^2
     <=> (h1^2 + h2^2) * T^2 - 2 * (v1 * h2^2 + v2 * h1^2 ) * T+ h2^2 * v1^2 + h1^2 * v2^2 - h1^2 * h2 ^2 / F^2 = 0 */


  float d1 = h1 * h1;
  float d2 = h2 * h2;
  float a = d1 + d2;
  float b = -2 * ( v1 * d2 + v2 * d1 );
  float c =  d2 * v1 * v1 + d1 * v2 * v2 - d1 * d2 / ( F * F );

  float q = b * b - 4 * a * c;

  if ( q < 0.0f || a == 0.0 ) {
    return -1.0;
  }

  float r = sqrt ( q );

  if ( r > b ) {
    return ( r - b ) / ( 2 * a );
  } else if ( r + b < 0.0 ) {
    return ( - r - b ) / ( 2 * a );
  } else {
    return -1.0;
  }
}

float Eikonal3d::solveQuadratic ( float v1, float v2, float v3,
                                  float h1, float h2, float h3,
                                  float F ) {
  float d1 = 1.0f / ( h1 * h1 );
  float d2 = 1.0f / ( h2 * h2 );
  float d3 = 1.0f / ( h3 * h3 );
  float a = d1 + d2 + d3;
  float b = -2 * ( v1 * d1 + v2 * d2 + v3 * d3 );
  float c =  d1 * v1 * v1 + d2 * v2 * v2 + d3 * v3 * v3 - 1.0f / ( F * F );

  float q = b * b - 4 * a * c;

  if ( q < 0.0f || a == 0.0 ) {
    return -1.0;
  }

  float r = sqrt ( q );

  if ( r > b ) {
    return ( r - b ) / ( 2 * a );
  } else if ( r + b < 0.0 ) {
    return ( - r - b ) / ( 2 * a );
  } else {
    return -1.0;
  }
}

void Eikonal3d::downwind ( int Dir, int X, int Y, int Z, int Step, float &Value, float /*MeanCurv*/ ) {
  // Direction should be one of TN_DIR_...
  if ( Dir < 0 || Dir > TN_DIR_MAX_NUM ) {
    cerr << "ERROR in Eikonal3d::downwind: Direction out of range!\n";
  }

  float F = getSpeed ( X, Y, Z, 0 );

  float H = ( Step ) / float ( numX );

  int NewX = X, NewY = Y, NewZ = Z;

  if ( Dir & TN_DIR_RIGHT ) {
    NewX += Step;
  } else {
    NewX -= Step;
  }
  if ( Dir & TN_DIR_UP ) {
    NewY += Step;
  } else {
    NewY -= Step;
  }
  if ( Dir & TN_DIR_BACK ) {
    NewZ += Step;
  } else {
    NewZ -= Step;
  }

  //TrialNode<float>::parsedir( Dir );
  //cerr << " ==> ( " << NewX << ", " << NewY << ", " << NewZ << ") ";

  /*
    float T010 = timeField->get ( X   , NewY, Z    );
    float T001 = timeField->get ( X   ,    Y, NewZ );
    float T011 = timeField->get ( X   , NewY, NewZ );

    float T100 = timeField->get ( NewX,    Y, Z    );
    float T110 = timeField->get ( NewX, NewY, Z    );
    float T101 = timeField->get ( NewX,    Y, NewZ );
    float T111 = timeField->get ( NewX, NewY, NewZ );

    T111 = 0.0f;
    T100 = T010 = T001 = sqrt ( 2.0f );
    T011 = T101 = T110 = 1.0f;
  */

  float vX = timeField->get ( NewX, Y, Z );
  float vY = timeField->get ( X, NewY, Z );
  float vZ = timeField->get ( X, Y, NewZ );

  if ( vX < 0.0 && vY < 0.0 && vZ < 0.0 ) {
    if ( FMM_WARNINGS )
      cerr << "all vX vY vZ < 0.0\n";
    return;
  }

  if ( vX >= 0.0f && vY >= 0.0 && vZ >= 0.0 ) {
    float newV = solveQuadratic ( vX, vY, vZ, H, H, H, F );
    if ( ( Value < 0.0 || Value > newV ) && newV >= 0.0f ) {
      Value = newV;
    }
  } else if ( vX < 0.0 && vY >= 0.0 && vZ >= 0.0 ) {
    float newV = solveQuadratic ( vY, vZ, H, H, F );
    if ( ( Value < 0.0 || Value > newV ) && newV >= 0.0 ) {
      Value = newV;
    }
  } else if ( vX >= 0.0 && vY < 0.0 && vZ >= 0.0 ) {
    float newV = solveQuadratic ( vX, vZ, H, H, F );
    if ( ( Value < 0.0 || Value > newV ) && newV >= 0.0 ) {
      Value = newV;
    }
  } else if ( vX >= 0.0 && vY >= 0.0 && vZ < 0.0 ) {
    float newV = solveQuadratic ( vX, vY, H, H, F );
    if ( ( Value < 0.0 || Value > newV ) && newV >= 0.0 ) {
      Value = newV;
    }
  } else if ( vX >= 0.0 && vY < 0.0 && vZ < 0.0 ) {
    float newV = vX + H / F;
    if ( ( Value < 0.0 || Value > newV ) && newV >= 0.0 )
      Value = newV;
  } else if ( vX < 0.0 && vY >= 0.0 && vZ < 0.0 ) {
    float newV = vY + H / F;
    if ( ( Value < 0.0 || Value > newV ) && newV >= 0.0 )
      Value = newV;
  } else if ( vX < 0.0 && vY < 0.0 && vZ >= 0.0 ) {
    float newV = vZ + H / F;
    if ( ( Value < 0.0 || Value > newV ) && newV >= 0.0 )
      Value = newV;
  }
}

int Eikonal3d::checkElementDataSet ( int X, int Y, int Z, int Level ) {
  int Step = 1 << ( GRID_DEPTH - Level );

  return ( timeField->get ( X, Y, Z ) >= 0.0f
           && timeField->get ( X + Step, Y, Z ) >= 0.0f
           && timeField->get ( X + Step, Y + Step, Z ) >= 0.0f
           && timeField->get ( X, Y + Step, Z ) >= 0.0f
           && timeField->get ( X, Y, Z + Step ) >= 0.0f
           && timeField->get ( X + Step, Y + Step, Z + Step ) >= 0.0f
           && timeField->get ( X + Step, Y, Z + Step ) >= 0.0f
           && timeField->get ( X, Y + Step, Z + Step ) >= 0.0f );
}

int Eikonal3d::checkElementSomeDataSet ( int X, int Y, int Z, int Level ) {
  int Step = 1 << ( GRID_DEPTH - Level );

  return ( timeField->get ( X, Y, Z ) >= 0.0
           || timeField->get ( X + Step, Y, Z ) >= 0.0
           || timeField->get ( X + Step, Y + Step, Z ) >= 0.0
           || timeField->get ( X, Y + Step, Z ) >= 0.0
           || timeField->get ( X, Y, Z + Step ) >= 0.0
           || timeField->get ( X + Step, Y + Step, Z + Step ) >= 0.0
           || timeField->get ( X + Step, Y, Z + Step ) >= 0.0
           || timeField->get ( X, Y + Step, Z + Step ) >= 0.0 );
}

void Eikonal3d::createElement ( int Direction, TrialNode<float> &Node ) {
  if ( Node.dirLevel[ Direction ] != TN_HN_NONE_SEG )
    return;

  int d0 = Node.dirLevel[ tn_nb_dirs[ Direction ][ 0 ] ],
           d1 = Node.dirLevel[ tn_nb_dirs[ Direction ][ 1 ] ],
                d2 = Node.dirLevel[ tn_nb_dirs[ Direction ][ 2 ] ],
                     d = -1;

  if ( d0 == TN_HN_NONE_SEG && d1 == TN_HN_NONE_SEG && d2 == TN_HN_NONE_SEG ) {
    cerr << "ERROR in Eikonal3d::createElement: all neighbouring dirLevels equal TN_HN_NONE_SEG!\n";
    Node.dump( );
    cerr << "Direction = ";
    TrialNode<float>::parsedir ( Direction );
    cerr << endl;
    abort( );
  }

  if ( d0 != TN_HN_NONE_SEG && d1 != TN_HN_NONE_SEG && iabs ( d0 - d1 ) == 2 ) {
    d = ( d0 + d1 ) >> 1;
  } else if ( d0 != TN_HN_NONE_SEG && d2 != TN_HN_NONE_SEG && iabs ( d0 - d2 ) == 2 ) {
    d = ( d0 + d2 ) >> 1;
  } else if ( d1 != TN_HN_NONE_SEG && d2 != TN_HN_NONE_SEG && iabs ( d1 - d2 ) == 2 ) {
    d = ( d1 + d2 ) >> 1;
  }

  if ( d != TN_HN_NONE_SEG ) {
    Node.dirLevel[ Direction ] = d;
    ocTree->setElement ( Direction, Node.x(), Node.y(), Node.z(), d );
    updateElementToTrialNodes ( Direction, Node.x(), Node.y(), Node.z(), d );
    return;
  }

  int RangeUp = 1000, RangeLow = -1;
  if ( d0 != TN_HN_NONE_SEG ) {
    if ( RangeUp > d0 + 1 )
      RangeUp = d0 + 1;
    if ( RangeLow < d0 - 1 )
      RangeLow = d0 - 1;
  }
  if ( d1 != TN_HN_NONE_SEG ) {
    if ( RangeUp > d1 + 1 )
      RangeUp = d1 + 1;
    if ( RangeLow < d1 - 1 )
      RangeLow = d1 - 1;
  }
  if ( d2 != TN_HN_NONE_SEG ) {
    if ( RangeUp > d2 + 1 )
      RangeUp = d2 + 1;
    if ( RangeLow < d2 - 1 )
      RangeLow = d2 - 1;
  }

  if ( RangeLow < 0 ) RangeLow = 0;
  if ( RangeUp > GRID_DEPTH ) RangeUp = GRID_DEPTH;

  for ( int level = RangeLow; level < RangeUp; level++ ) {
    if ( ocTree->nodeInGrid ( Node.x(), Node.y(), Node.z(), level )
         && estimator->checkElementSure ( Direction, Node.x(), Node.y(), Node.z(), level ) ) {

      ocTree->setElement ( Direction, Node.x(), Node.y(), Node.z(), level );
      updateElementToTrialNodes ( Direction, Node.x(), Node.y(), Node.z(), level );
      Node.dirLevel[ Direction ] = level;
      return;
    }
  }

  // else rangeUp on finest level
  ocTree->setElement ( Direction, Node.x(), Node.y(), Node.z(), RangeUp );
  updateElementToTrialNodes ( Direction, Node.x(), Node.y(), Node.z(), RangeUp );
  Node.dirLevel[ Direction ] = RangeUp;
}

}   // end namespace eik

#endif
