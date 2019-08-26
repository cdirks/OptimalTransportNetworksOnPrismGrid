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

#include <eikonal3d.h>

#include <quoc.h>
#include <aol.h>

namespace eik {

static QC_SEED_TYPE seed_type = QC_POS;

Eikonal3d::Eikonal3d( )
    : scaleTime ( 1.0f ),
    ocTree ( NULL ),
    timeField ( NULL ),
    finalTimeField ( NULL ),
    tmpTimeField ( NULL ),
    tagField ( NULL ),
    indexField ( NULL ),
    complementField ( NULL ),
    seedTypeField ( NULL )
    {
  createLookUpTables( );
}

Eikonal3d::~Eikonal3d( ) {
  if ( ocTree != NULL )
    delete ocTree;
  if ( timeField != NULL )
    delete timeField;
  if ( finalTimeField != NULL )
    delete finalTimeField;
  if ( tmpTimeField != NULL )
    delete tmpTimeField;
  if ( tagField != NULL )
    delete tagField;
  if ( indexField != NULL )
    delete indexField;

  if ( seedTypeField != NULL )
    delete seedTypeField;
}


void Eikonal3d::createLookUpTables( ) {
  cout << "- creating lookuptable \"tn_dir_order\"..." << endl;
  tn_dir_ordering[ 0 ] = TN_DIR_BACK_UP_RIGHT;
  tn_dir_ordering[ 1 ] = TN_DIR_BACK_UP_LEFT;
  tn_dir_ordering[ 2 ] = TN_DIR_BACK_LOW_LEFT;
  tn_dir_ordering[ 3 ] = TN_DIR_BACK_LOW_RIGHT;
  tn_dir_ordering[ 4 ] = TN_DIR_FRONT_LOW_RIGHT;
  tn_dir_ordering[ 5 ] = TN_DIR_FRONT_LOW_LEFT;
  tn_dir_ordering[ 6 ] = TN_DIR_FRONT_UP_LEFT;
  tn_dir_ordering[ 7 ] = TN_DIR_FRONT_UP_RIGHT;

  cout << "- creating lookuptable \"tn_nb_dirs\"..." << endl;
  for ( int i = 0; i < 8; i++ ) {
    for ( int j = 0; j < 3; j++ ) {
      int d = 1 << j;
      if ( i & d ) {
        tn_nb_dirs[ i ][ j ] = i & ~d;
      } else {
        tn_nb_dirs[ i ][ j ] = i | d;
      }
    }
  }

  cout << "- creating lookuptable \"tn_hn_nb_cell_dirs\"..." << endl;

  // first the hanging nodes on edges
  tn_hn_nb_cell_dirs[ TN_HN_EDGE_FRONT_UP ] = 0xFF &
                                              ( ~ ( ( 1 << TN_DIR_BACK_LOW_LEFT )
                                                    | ( 1 << TN_DIR_BACK_LOW_RIGHT )
                                                    | ( 1 << TN_DIR_FRONT_UP_LEFT )
                                                    | ( 1 << TN_DIR_FRONT_UP_RIGHT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_FRONT_LOW ] = 0xFF &
                                               ( ~ ( ( 1 << TN_DIR_BACK_UP_LEFT )
                                                     | ( 1 << TN_DIR_BACK_UP_RIGHT )
                                                     | ( 1 << TN_DIR_FRONT_LOW_LEFT )
                                                     | ( 1 << TN_DIR_FRONT_LOW_RIGHT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_FRONT_LEFT ] = 0xFF &
                                                ( ~ ( ( 1 << TN_DIR_BACK_UP_RIGHT )
                                                      | ( 1 << TN_DIR_BACK_LOW_RIGHT )
                                                      | ( 1 << TN_DIR_FRONT_LOW_LEFT )
                                                      | ( 1 << TN_DIR_FRONT_UP_LEFT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_FRONT_RIGHT ] = 0xFF &
                                                 ( ~ ( ( 1 << TN_DIR_BACK_UP_LEFT )
                                                       | ( 1 << TN_DIR_BACK_LOW_LEFT )
                                                       | ( 1 << TN_DIR_FRONT_LOW_RIGHT )
                                                       | ( 1 << TN_DIR_FRONT_UP_RIGHT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_BACK_UP ] = 0xFF &
                                             ( ~ ( ( 1 << TN_DIR_FRONT_LOW_LEFT )
                                                   | ( 1 << TN_DIR_FRONT_LOW_RIGHT )
                                                   | ( 1 << TN_DIR_BACK_UP_RIGHT )
                                                   | ( 1 << TN_DIR_BACK_UP_LEFT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_BACK_LOW ] = 0xFF &
                                              ( ~ ( ( 1 << TN_DIR_FRONT_UP_LEFT )
                                                    | ( 1 << TN_DIR_FRONT_UP_RIGHT )
                                                    | ( 1 << TN_DIR_BACK_LOW_LEFT )
                                                    | ( 1 << TN_DIR_BACK_LOW_RIGHT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_BACK_LEFT ] = 0xFF &
                                               ( ~ ( ( 1 << TN_DIR_FRONT_UP_RIGHT )
                                                     | ( 1 << TN_DIR_FRONT_LOW_RIGHT )
                                                     | ( 1 << TN_DIR_BACK_UP_LEFT )
                                                     | ( 1 << TN_DIR_BACK_LOW_LEFT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_BACK_RIGHT ] = 0xFF &
                                                ( ~ ( ( 1 << TN_DIR_FRONT_UP_LEFT )
                                                      | ( 1 << TN_DIR_FRONT_LOW_LEFT )
                                                      | ( 1 << TN_DIR_BACK_UP_RIGHT )
                                                      | ( 1 << TN_DIR_BACK_LOW_RIGHT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_RIGHT_UP ] = 0xFF &
                                              ( ~ ( ( 1 << TN_DIR_FRONT_LOW_LEFT )
                                                    | ( 1 << TN_DIR_BACK_LOW_LEFT )
                                                    | ( 1 << TN_DIR_FRONT_UP_RIGHT )
                                                    | ( 1 << TN_DIR_BACK_UP_RIGHT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_RIGHT_LOW ] = 0xFF &
                                               ( ~ ( ( 1 << TN_DIR_FRONT_UP_LEFT )
                                                     | ( 1 << TN_DIR_BACK_UP_LEFT )
                                                     | ( 1 << TN_DIR_BACK_LOW_RIGHT )
                                                     | ( 1 << TN_DIR_FRONT_LOW_RIGHT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_LEFT_UP ] = 0xFF &
                                             ( ~ ( ( 1 << TN_DIR_FRONT_LOW_RIGHT )
                                                   | ( 1 << TN_DIR_BACK_LOW_RIGHT )
                                                   | ( 1 << TN_DIR_FRONT_UP_LEFT )
                                                   | ( 1 << TN_DIR_BACK_UP_LEFT ) ) );

  tn_hn_nb_cell_dirs[ TN_HN_EDGE_LEFT_LOW ] = 0xFF &
                                              ( ~ ( ( 1 << TN_DIR_FRONT_UP_RIGHT )
                                                    | ( 1 << TN_DIR_BACK_UP_RIGHT )
                                                    | ( 1 << TN_DIR_BACK_LOW_LEFT )
                                                    | ( 1 << TN_DIR_FRONT_LOW_LEFT ) ) );

  // hanging nodes on faces:
  tn_hn_nb_cell_dirs[ TN_HN_FACE_UP ] = ( 1 << TN_DIR_BACK_UP_LEFT ) |
                                        ( 1 << TN_DIR_BACK_UP_RIGHT ) |
                                        ( 1 << TN_DIR_FRONT_UP_RIGHT ) |
                                        ( 1 << TN_DIR_FRONT_UP_LEFT );

  tn_hn_nb_cell_dirs[ TN_HN_FACE_LOW ] = ( 1 << TN_DIR_BACK_LOW_LEFT ) |
                                         ( 1 << TN_DIR_BACK_LOW_RIGHT ) |
                                         ( 1 << TN_DIR_FRONT_LOW_RIGHT ) |
                                         ( 1 << TN_DIR_FRONT_LOW_LEFT );

  tn_hn_nb_cell_dirs[ TN_HN_FACE_RIGHT ] = ( 1 << TN_DIR_BACK_LOW_RIGHT ) |
                                           ( 1 << TN_DIR_BACK_UP_RIGHT ) |
                                           ( 1 << TN_DIR_FRONT_LOW_RIGHT ) |
                                           ( 1 << TN_DIR_FRONT_UP_RIGHT );

  tn_hn_nb_cell_dirs[ TN_HN_FACE_LEFT ] = ( 1 << TN_DIR_BACK_LOW_LEFT ) |
                                          ( 1 << TN_DIR_BACK_UP_LEFT ) |
                                          ( 1 << TN_DIR_FRONT_LOW_LEFT ) |
                                          ( 1 << TN_DIR_FRONT_UP_LEFT );

  tn_hn_nb_cell_dirs[ TN_HN_FACE_BACK ] = ( 1 << TN_DIR_BACK_LOW_LEFT ) |
                                          ( 1 << TN_DIR_BACK_LOW_RIGHT ) |
                                          ( 1 << TN_DIR_BACK_UP_LEFT ) |
                                          ( 1 << TN_DIR_BACK_UP_RIGHT );

  tn_hn_nb_cell_dirs[ TN_HN_FACE_FRONT ] = ( 1 << TN_DIR_FRONT_LOW_LEFT ) |
                                           ( 1 << TN_DIR_FRONT_LOW_RIGHT ) |
                                           ( 1 << TN_DIR_FRONT_UP_LEFT ) |
                                           ( 1 << TN_DIR_FRONT_UP_RIGHT );

  cout << "- creating lookuptable \"tn_hn_offs\"..." << endl;
  tn_hn_offs[ TN_HN_EDGE_FRONT_UP ][ 0 ] = 1;
  tn_hn_offs[ TN_HN_EDGE_FRONT_UP ][ 1 ] = 2;
  tn_hn_offs[ TN_HN_EDGE_FRONT_UP ][ 2 ] = 0;

  tn_hn_offs[ TN_HN_EDGE_FRONT_LOW ][ 0 ] = 1;
  tn_hn_offs[ TN_HN_EDGE_FRONT_LOW ][ 1 ] = 0;
  tn_hn_offs[ TN_HN_EDGE_FRONT_LOW ][ 2 ] = 0;

  tn_hn_offs[ TN_HN_EDGE_FRONT_LEFT ][ 0 ] = 0;
  tn_hn_offs[ TN_HN_EDGE_FRONT_LEFT ][ 1 ] = 1;
  tn_hn_offs[ TN_HN_EDGE_FRONT_LEFT ][ 2 ] = 0;

  tn_hn_offs[ TN_HN_EDGE_FRONT_RIGHT ][ 0 ] = 2;
  tn_hn_offs[ TN_HN_EDGE_FRONT_RIGHT ][ 1 ] = 1;
  tn_hn_offs[ TN_HN_EDGE_FRONT_RIGHT ][ 2 ] = 0;

  tn_hn_offs[ TN_HN_EDGE_BACK_UP ][ 0 ] = 1;
  tn_hn_offs[ TN_HN_EDGE_BACK_UP ][ 1 ] = 2;
  tn_hn_offs[ TN_HN_EDGE_BACK_UP ][ 2 ] = 2;

  tn_hn_offs[ TN_HN_EDGE_BACK_LOW ][ 0 ] = 1;
  tn_hn_offs[ TN_HN_EDGE_BACK_LOW ][ 1 ] = 0;
  tn_hn_offs[ TN_HN_EDGE_BACK_LOW ][ 2 ] = 2;

  tn_hn_offs[ TN_HN_EDGE_BACK_LEFT ][ 0 ] = 0;
  tn_hn_offs[ TN_HN_EDGE_BACK_LEFT ][ 1 ] = 1;
  tn_hn_offs[ TN_HN_EDGE_BACK_LEFT ][ 2 ] = 2;

  tn_hn_offs[ TN_HN_EDGE_BACK_RIGHT ][ 0 ] = 2;
  tn_hn_offs[ TN_HN_EDGE_BACK_RIGHT ][ 1 ] = 1;
  tn_hn_offs[ TN_HN_EDGE_BACK_RIGHT ][ 2 ] = 2;

  tn_hn_offs[ TN_HN_EDGE_RIGHT_UP ][ 0 ] = 2;
  tn_hn_offs[ TN_HN_EDGE_RIGHT_UP ][ 1 ] = 2;
  tn_hn_offs[ TN_HN_EDGE_RIGHT_UP ][ 2 ] = 1;

  tn_hn_offs[ TN_HN_EDGE_RIGHT_LOW ][ 0 ] = 2;
  tn_hn_offs[ TN_HN_EDGE_RIGHT_LOW ][ 1 ] = 0;
  tn_hn_offs[ TN_HN_EDGE_RIGHT_LOW ][ 2 ] = 1;

  tn_hn_offs[ TN_HN_EDGE_LEFT_UP ][ 0 ] = 0;
  tn_hn_offs[ TN_HN_EDGE_LEFT_UP ][ 1 ] = 2;
  tn_hn_offs[ TN_HN_EDGE_LEFT_UP ][ 2 ] = 1;

  tn_hn_offs[ TN_HN_EDGE_LEFT_LOW ][ 0 ] = 0;
  tn_hn_offs[ TN_HN_EDGE_LEFT_LOW ][ 1 ] = 0;
  tn_hn_offs[ TN_HN_EDGE_LEFT_LOW ][ 2 ] = 1;

  tn_hn_offs[ TN_HN_FACE_UP ][ 0 ] = 1;
  tn_hn_offs[ TN_HN_FACE_UP ][ 1 ] = 2;
  tn_hn_offs[ TN_HN_FACE_UP ][ 2 ] = 1;

  tn_hn_offs[ TN_HN_FACE_LOW ][ 0 ] = 1;
  tn_hn_offs[ TN_HN_FACE_LOW ][ 1 ] = 0;
  tn_hn_offs[ TN_HN_FACE_LOW ][ 2 ] = 1;

  tn_hn_offs[ TN_HN_FACE_RIGHT ][ 0 ] = 2;
  tn_hn_offs[ TN_HN_FACE_RIGHT ][ 1 ] = 1;
  tn_hn_offs[ TN_HN_FACE_RIGHT ][ 2 ] = 1;

  tn_hn_offs[ TN_HN_FACE_LEFT ][ 0 ] = 0;
  tn_hn_offs[ TN_HN_FACE_LEFT ][ 1 ] = 1;
  tn_hn_offs[ TN_HN_FACE_LEFT ][ 2 ] = 1;

  tn_hn_offs[ TN_HN_FACE_FRONT ][ 0 ] = 1;
  tn_hn_offs[ TN_HN_FACE_FRONT ][ 1 ] = 1;
  tn_hn_offs[ TN_HN_FACE_FRONT ][ 2 ] = 0;

  tn_hn_offs[ TN_HN_FACE_BACK ][ 0 ] = 1;
  tn_hn_offs[ TN_HN_FACE_BACK ][ 1 ] = 1;
  tn_hn_offs[ TN_HN_FACE_BACK ][ 2 ] = 2;

}

void Eikonal3d::initComplementField( ) {
  if ( complementField == NULL ) {
    complementField = new qc::ScalarArray<unsigned char, qc::QC_3D> ( numX + 1, numY + 1, numZ + 1 );
  }
  complementField->setAll ( 0 );
}

void Eikonal3d::recMarkTouchedCell ( int X, int Y, int Z, int level ) {
  int hs = 1 << ( GRID_DEPTH - level - 1 );

  if ( ( estimator->checkElementSure ( X, Y, Z, level ) &&
         checkElementSomeDataSet ( X, Y, Z, level ) ) ||
       level == GRID_DEPTH ) {
    ocTree->setElement ( X, Y, Z, level );
  } else {
    recMarkTouchedCell ( X, Y, Z, level + 1 );
    recMarkTouchedCell ( X + hs, Y, Z, level + 1 );
    recMarkTouchedCell ( X, Y + hs, Z, level + 1 );
    recMarkTouchedCell ( X + hs, Y + hs, Z, level + 1 );
    recMarkTouchedCell ( X, Y, Z + hs, level + 1 );
    recMarkTouchedCell ( X + hs, Y, Z + hs, level + 1 );
    recMarkTouchedCell ( X, Y + hs, Z + hs, level + 1 );
    recMarkTouchedCell ( X + hs, Y + hs, Z + hs, level + 1 );
  }
}

void Eikonal3d::refine ( int X, int Y, int Z, int level ) {
  int step = 1 << ( GRID_DEPTH - level );
  for ( int x = 0; x <= step; x += step ) {
    for ( int y = 0; y <= step; y += step ) {
      for ( int z = 0; z <= step; z += step ) {
        if ( estimator->checkElementSure ( X + x, Y + y, Z + z, level ) || level == GRID_DEPTH ) {
          if ( level < GRID_DEPTH )
            ocTree->setElement ( X + x, Y + y, Z + z, level );
        } else {
          refine ( X + x, Y + y, Z + z, level + 1 );
        }
      }
    }
  }
}

void Eikonal3d::reset( ) {
  trialNodeHeap.erase ( trialNodeHeap.begin( ), trialNodeHeap.end( ) );
  hangNodes.erase ( hangNodes.begin( ), hangNodes.end( ) );
  interpolatedHangNodes.erase ( interpolatedHangNodes.begin( ), interpolatedHangNodes.end( ) );
  ocTree->clearAll( );
  timeField->setAll ( -1.0f );
  indexField->setAll ( -1 );
  seedTypeField->setAll ( 0 );

  seed_type = QC_POS;

}

void Eikonal3d::addSeedPoint ( int X, int Y, int Z ) {

  cerr << "Adding SeedPoint/DrySeedPoint at position " << X << " " << Y << " " << Z << endl;

  int NewX = X,
             NewY = Y,
                    NewZ = Z,
                           Step = 1;

  for ( int Level = GRID_DEPTH; Level > 0; Level -- ) {
    if ( estimator->checkElementSure ( X - ( X % Step ),
                                       Y - ( Y % Step ),
                                       Z - ( Z % Step ),
                                       Level ) ) {
      NewX = X - ( X % Step );
      NewY = Y - ( Y % Step );
      NewZ = Z - ( Z % Step );
      Step <<= 1;
    } else {
      break;
    }
  }

  TrialNode<float> Node ( NewX, NewY, NewZ );

  Node.seed_type = seed_type;

  cerr << "seed_type = " << seed_type << endl;

  if ( seed_type == QC_POS ) {
    seed_type = QC_NEG;
  } else {
    seed_type = QC_POS;
  }

  Node.val = peekMinTime( );
  initSeedElement ( Node );
  trialNodeHeap.push ( Node );
}

void Eikonal3d::init ( int NumX, int NumY, int NumZ ) {
  numX = NumX;
  numY = NumY;
  numZ = NumZ;

  GRID_DEPTH = qc::logBaseTwo ( NumX );

  cerr << "- GRID_DEPTH = " << GRID_DEPTH << endl;

  if ( ocTree != NULL )
    delete ocTree;
  ocTree = new qc::OcTree ( numX, numY, numZ, GRID_DEPTH );
  ocTree->clearAll( );
  if ( timeField != NULL )
    delete timeField;

  timeField = new qc::ScalarArray<float, qc::QC_3D> ( numX + 1, numY + 1, numZ + 1 );
  // timeField->initL2Projection( 3 ); // according to Marc this is not necessary
  timeField->setAll ( -1.0 );
  if ( seedTypeField != NULL )
    delete seedTypeField;

  seedTypeField = new qc::ScalarArray<unsigned char, qc::QC_3D> ( numX + 1, numY + 1, numZ + 1 );

  seedTypeField->setAll ( 0 );
  if ( finalTimeField != NULL )
    delete finalTimeField;

#ifdef SUCCESSIVE_SEGMENT
  finalTimeField = new qc::ScalarArray<float, qc::QC_3D> ( numX + 1, numY + 1, numZ + 1 );
  finalTimeField->setAll ( 100.0f );
#endif
  if ( tmpTimeField != NULL )
    delete tmpTimeField;

#ifdef SUCCESSIVE_SEGMENT
  tmpTimeField = new qc::ScalarArray<float, qc::QC_3D> ( numX + 1, numY + 1, numZ + 1 );
#endif
  if ( tagField != NULL )
    delete tagField;

  tagField = new qc::ScalarArray<unsigned char, qc::QC_3D> ( numX + 1, numY + 1, numZ + 1 );
  tagField->setAll ( 0 );
  if ( indexField != NULL )
    delete indexField;

  indexField = new qc::ScalarArray<int, qc::QC_3D> ( numX + 1, numY + 1, numZ + 1 );
  indexField->setAll ( -1 );
  trialNodeHeap.setIndexField ( indexField );
  cerr << "- MarchSegment3d has been successfully initialized!\n";
}

void Eikonal3d::dumpTrialNodes( ) {
  vector<TrialNodeFloat>::const_iterator it;
  for ( it = trialNodeHeap.begin( ); it != trialNodeHeap.end( ); ++it ) {
    const TrialNode<float> &node = *it;
    node.dump( );
  }
}

void Eikonal3d::dumpTrialNode ( int X, int Y, int Z ) {
  vector<TrialNodeFloat>::const_iterator it;
  for ( it = trialNodeHeap.begin( ); it != trialNodeHeap.end( ); ++it ) {
    const TrialNode<float> &node = *it;
    if ( node.x() == X && node.y() == Y && node.z() == Z ) {
      node.dump( );
      break;
    }
  }
  cerr << "TrialNode not found." << endl;
}

void Eikonal3d::dumpMinTrialNode( ) {
  TrialNode<float> Node;
  trialNodeHeap.peekMin ( Node );
  Node.dump( );
}


void Eikonal3d::initSeedElement ( TrialNode<float> &Node ) {
  for ( int dir = 0; dir <= TN_DIR_MAX_NUM; dir++ ) {
    if ( Node.dirLevel[ dir ] != -1 ) {
      cerr << "WARNING in initSeedElement: Node already has neighbours!\n";
      return;
    }
  }
  for ( int level = GRID_DEPTH; level >= 0; level-- ) {
    if ( ocTree->nodeInGrid ( Node.x(), Node.y(), Node.z(), level )
         && estimator->checkElementSure ( tn_dir_ordering[ 0 ], Node.x(), Node.y(), Node.z(), level ) ) {
      Node.dirLevel[ tn_dir_ordering[ 0 ] ] = level;
    }
  }

  ocTree->setElement ( tn_dir_ordering[ 0 ],
                       Node.x(), Node.y(), Node.z(),
                       Node.dirLevel[ tn_dir_ordering[ 0 ] ] );

  // now the other 7 directions
  for ( int dind = 1; dind <= TN_DIR_MAX_NUM; dind++ ) {
    int prevLevel = Node.dirLevel[ tn_dir_ordering[ dind - 1 ]];
    if ( prevLevel > 0 &&
         ocTree->nodeInGrid ( Node.x(), Node.y(), Node.z(), prevLevel - 1 )
         && estimator->checkElementSure ( tn_dir_ordering[ dind ],
                                          Node.x(), Node.y(), Node.z(), prevLevel - 1 ) ) {
      Node.dirLevel[ tn_dir_ordering[ dind ] ] = prevLevel - 1;
    } else if ( ( ocTree->nodeInGrid ( Node.x(), Node.y(), Node.z(), prevLevel )
      && estimator->checkElementSure ( tn_dir_ordering[ dind ],
               Node.x(), Node.y(), Node.z(), prevLevel ) )
                || ( prevLevel == GRID_DEPTH ) ) {
      // ( A && B || C )  is ( (A&&B) || C ) but OS has no idea whether this is correct here ...
      Node.dirLevel[ tn_dir_ordering[ dind ] ] = prevLevel;
    } else {
      Node.dirLevel[ tn_dir_ordering[ dind ] ] = prevLevel + 1;
    }
    ocTree->setElement ( tn_dir_ordering[ dind ],
                         Node.x(), Node.y(), Node.z(),
                         Node.dirLevel[ tn_dir_ordering[ dind ] ] );
  }
}

// void Eikonal3d::makeTrialNodeActive( TrialNode<float> &Node, int SweepMode=0 )
void Eikonal3d::makeTrialNodeActive ( TrialNode<float> &Node ) {
  // how to detect hanging nodes on edges?
  // for each node find element offset of coarse element and test estimator.
  // how to determine level of coarse element?
  // one of Node.dirLevel must be different from -1 (6 directions for edges, 4 for faces)
  // take maximum value and subtract one

  // here detect hanging nodes;

#ifdef DEBUG
  cerr << " - Eikonal3d::makeTrialNodeActive: \n";
  Node.dump( );
#endif

  if ( Node.val < 0.0 )
    return;

  if ( Node.x() == 0 || Node.x() == numX ||
       Node.y() == 0 || Node.y() == numY ||
       Node.z() == 0 || Node.z() == numZ  ) {
    timeField->set ( Node.x(), Node.y(), Node.z(), Node.val );
    return;
  }

  // detected all hanging nodes here
  // choose maximum(finest) level of all dirlevels
  int max = -1;
  int dir;
  for ( dir = 0; dir <= TN_DIR_MAX_NUM; dir++ ) {
    if ( Node.dirLevel[ dir ] != -1 && Node.dirLevel[ dir ] > max ) {
      max = Node.dirLevel[ dir ];
    }
  }
  if ( max < 0 ) {
    cerr << "ERROR in Eikonal3d::makeTrialNodeActive: all dirLevels < 0 \n";
    abort( );
  }
  for ( int hn = 0; hn <= TN_HN_MAX_NUM; hn++ ) {
    if ( max > 0 ) {
      int coarseLevel = max - 1;
      int HalfStep = 1 << ( GRID_DEPTH - max );
      int ccX = Node.x() - tn_hn_offs[ hn ][ 0 ] * HalfStep,
                ccY = Node.y() - tn_hn_offs[ hn ][ 1 ] * HalfStep,
                      ccZ = Node.z() - tn_hn_offs[ hn ][ 2 ] * HalfStep;
      if ( ocTree->nodeInGrid ( ccX, ccY, ccZ, coarseLevel )
           && estimator->checkElementSure ( ccX, ccY, ccZ, coarseLevel ) ) {
        if ( !ocTree->getElement ( ccX, ccY, ccZ, coarseLevel ) ) {
          ocTree->setElement ( ccX, ccY, ccZ, coarseLevel );
          updateElementToTrialNodes ( ccX, ccY, ccZ, coarseLevel );
        }
        Node.addHangNode ( hn, coarseLevel );
        hangNodes.push_back ( Node );
        // cerr << "HangNode Detected : " << hn << endl;
      }
    }
  }

#ifdef DEBUG
  if ( !Node.consistencyCheck ( GRID_DEPTH, *estimator, *this ) ) {
    abort( );
  }
#endif

  localHangNodeGridUpdate ( Node );

  int hasBeenSet = 0;
  // fill up Tree for normal nodes
  if ( Node.faceHangNode == TN_HN_NONE_SEG
       && Node.edgeHangNode == TN_HN_NONE_SEG ) {
    int firstNonEmptyDirection = -1;
    for ( dir = 0; dir <= TN_DIR_MAX_NUM; dir++ ) {
      if ( Node.dirLevel[ tn_dir_ordering[ dir ] ] != -1 ) {
        firstNonEmptyDirection = dir;
        break;
      }
    }
    if ( firstNonEmptyDirection == -1 ) {
      cerr << "ERROR in makeTrialNodeActive: all dirLevels of Node are -1\n";
      abort();
    }

    for ( dir = firstNonEmptyDirection; dir <= TN_DIR_MAX_NUM; dir++ ) {
      createElement ( tn_dir_ordering[ dir ], Node );
    }
    for ( dir = 0; dir < firstNonEmptyDirection; dir++ ) {
      createElement ( tn_dir_ordering[ dir ], Node );
    }

#ifdef DEBUG
    for ( dir = 0; dir <= TN_DIR_MAX_NUM; dir++ ) {
      if ( Node.dirLevel[ dir ] == TN_HN_NONE_SEG ) {
        cerr << "\n\nCould not set level for ";
        TrialNode<float>::parsedir ( dir );
        cerr << ".\n";
        Node.dump( );
        abort( );
      }
    }
#endif
    float v = timeField->get ( Node.x(), Node.y(), Node.z() );
    if ( ( v < 0.0 || Node.val < v ) && Node.val >= 0.0 ) {
      timeField->set ( Node.x(), Node.y(), Node.z(), Node.val );
      seedTypeField->set ( Node.x(), Node.y(), Node.z(), Node.seed_type );
      hasBeenSet = 1;
    }
  } else {
    /*if ( timeField->get( Node.x(), Node.y(), Node.z() ) < 0.0 ) {
      addNotInterpolatedHangNode( Node );
    }*/
    return;
  }

  //cerr << " - Node after createElements: \n";
  //Node.dump( );
#ifdef DEBUG
  if ( !Node.consistencyCheck ( GRID_DEPTH, *estimator, *this ) ) {
    abort( );
  }
#endif
  // - interpolation of hanging nodes

  // check for hanging nodes on the 8 surrounding edges
  if ( hasBeenSet ) {
    for ( int axis = 0; axis < 3; axis++ ) {

      int thisdir = 1 << axis;

      int dira = 1 << ( axis + 1 ) % 3;
      int dirb = 1 << ( axis + 2 ) % 3;

      for ( int ori = 0; ori < 2; ori++ ) {
        int minLevel = GRID_DEPTH + 1;
        int maxLevel = -1;
        int initDir;
        if ( ori == 0 )
          initDir = 0;
        else
          initDir = thisdir;

        for ( int i = 0; i < 4; i++ ) {

          int dir = initDir;
          if ( i & 1 )
            dir |= dira;
          if ( i & 2 )
            dir |= dirb;

          int dl = Node.dirLevel[ dir ];
          if ( dl != TN_HN_NONE_SEG && dl < minLevel ) {
            minLevel = dl;
          }
          if ( dl != TN_HN_NONE_SEG && dl > maxLevel ) {
            maxLevel = dl;
          }
          if ( minLevel != GRID_DEPTH + 1 && minLevel != maxLevel ) {
            int checkX = Node.x(),
                         checkY = Node.y(),
                                  checkZ = Node.z();
            int step = 1 << ( GRID_DEPTH - minLevel );
            if ( thisdir == TN_DIR_RIGHT ) {
              if ( ori == 0 )
                checkX -= step;
              else
                checkX += step;
            } else if ( thisdir == TN_DIR_UP ) {
              if ( ori == 0 )
                checkY -= step;
              else
                checkY += step;
            } else if ( thisdir == TN_DIR_BACK ) {
              if ( ori == 0 )
                checkZ -= step;
              else
                checkZ += step;
            }

            float v;
            if ( ( v = timeField->get ( checkX, checkY, checkZ ) ) >= 0.0 ) {
              int intX, intY, intZ;
              intX = ( checkX + Node.x() ) >> 1;
              intY = ( checkY + Node.y() ) >> 1;
              intZ = ( checkZ + Node.z() ) >> 1;
              timeField->set ( intX, intY, intZ, ( v + Node.val ) *0.5f );
              interpolatedHangNodes.push_back ( TrialNode<float> ( intX, intY, intZ ) );

              if ( minLevel < GRID_DEPTH - 1 ) {
                int stepX = ( checkX - Node.x() ) >> 2;
                int stepY = ( checkY - Node.y() ) >> 2;
                int stepZ = ( checkZ - Node.z() ) >> 2;

                intX = Node.x() + stepX;
                intY = Node.y() + stepY;
                intZ = Node.z() + stepZ;
                timeField->set ( intX, intY, intZ, 0.75f * Node.val + 0.25f * v );
                interpolatedHangNodes.push_back ( TrialNode<float> ( intX, intY, intZ ) );

                intX = Node.x() + 3 * stepX;
                intY = Node.y() + 3 * stepY;
                intZ = Node.z() + 3 * stepZ;
                timeField->set ( intX, intY, intZ, 0.25f * Node.val + 0.75f * v );
                interpolatedHangNodes.push_back ( TrialNode<float> ( intX, intY, intZ ) );
              }
            }
          }
        }
      }
    }
  }
#ifdef DEBUG
  if ( !Node.consistencyCheck ( GRID_DEPTH, *estimator, *this ) ) {
    abort( );
  }
#endif

  // check for hanging nodes on the 12 surrounding faces

  if ( hasBeenSet ) {
    for ( int i = 0; i < 3; i++ ) {
      int dira = ( 1 << i );
      int dirb = ( 1 << ( ( i + 1 ) % 3 ) );
      int dirc = ( 1 << ( ( i + 2 ) % 3 ) );

      for ( int j = 0; j < 4; j++ ) {
        int dir = 0;
        if ( j & 1 ) dir |= dira;
        if ( j & 2 ) dir |= dirb;

        int tdira = dir;
        int tdirb = ( dir | dirc );

        // fi there is a hanging node on the face, levels must be different
        if ( Node.dirLevel[ tdira ] != Node.dirLevel[ tdirb ]
             && Node.dirLevel[ tdira ] != TN_HN_NONE_SEG
             && Node.dirLevel[ tdirb ] != TN_HN_NONE_SEG ) {
          int coarseLevel = imin ( Node.dirLevel[ tdira ], Node.dirLevel[ tdirb ] );
          int Step = 1 << ( GRID_DEPTH - coarseLevel );
          float v[ 3 ];
          int checkX, checkY, checkZ;
          int interpolate = 1;

          for ( int k = 0; k < 3; k++ ) {
            // make linear combinations dira, dirb, dira + dirb
            int kaux = k + 1;
            checkX = Node.x();
            checkY = Node.y();
            checkZ = Node.z();

            if ( kaux & 1 ) {
              if ( dira == TN_DIR_RIGHT ) {
                if ( j & 1 ) {
                  checkX += Step;
                } else {
                  checkX -= Step;
                }
              } else if ( dira == TN_DIR_UP ) {
                if ( j & 1 ) {
                  checkY += Step;
                } else {
                  checkY -= Step;
                }
              } else if ( dira == TN_DIR_BACK ) {
                if ( j & 1 ) {
                  checkZ += Step;
                } else {
                  checkZ -= Step;
                }
              }
            }
            if ( kaux & 2 ) {
              if ( dirb == TN_DIR_RIGHT ) {
                if ( j & 2 ) {
                  checkX += Step;
                } else {
                  checkX -= Step;
                }
              } else if ( dirb == TN_DIR_UP ) {
                if ( j & 2 ) {
                  checkY += Step;
                } else {
                  checkY -= Step;
                }
              } else if ( dirb == TN_DIR_BACK ) {
                if ( j & 2 ) {
                  checkZ += Step;
                } else {
                  checkZ -= Step;
                }
              }
            }
            v[ k ] = timeField->get ( checkX, checkY, checkZ );
            if ( v[ k ] < 0.0 ) {
              interpolate = 0;
              break;
            }
          }
          if ( interpolate ) {
            float newVal = 0.25f * ( v[0] + v[1] + v[2] + Node.val );

            // if interpolate is 1, checkX...checkY contain coordinates
            // of node across the face, thus:
            int intX = ( Node.x() + checkX ) >> 1,
                       intY = ( Node.y() + checkY ) >> 1,
                              intZ = ( Node.z() + checkZ ) >> 1;

            timeField->set ( intX, intY, intZ, newVal );

            int sdir = 0;
            int newNodeX = intX,
                           newNodeY = intY,
                                      newNodeZ = intZ;
            if ( Node.dirLevel[ tdira ] < Node.dirLevel[ tdirb ] ) {
              // downwind in direction of dirc
              switch ( dirc ) {
              case TN_DIR_RIGHT:
                newNodeX += ( Step >> 1 );
                break;
              case TN_DIR_UP:
                newNodeY += ( Step >> 1 );
                break;
              case TN_DIR_BACK:
                newNodeZ += ( Step >> 1 );
                break;
              default:
                cerr << "ERROR: dirc none of { TN_DIR_RIGHT, TN_DIR_UP, TN_DIR_BACK }\n";
                abort( );
              }
              sdir = 0;
            } else {
              // downwind in direction of -dirc
              switch ( dirc ) {
              case TN_DIR_RIGHT:
                newNodeX -= ( Step >> 1 );
                break;
              case TN_DIR_UP:
                newNodeY -= ( Step >> 1 );
                break;
              case TN_DIR_BACK:
                newNodeZ -= ( Step >> 1 );
                break;
              default:
                cerr << "ERROR: dirc none of { TN_DIR_RIGHT, TN_DIR_UP, TN_DIR_BACK }\n";
                abort( );
              }
              sdir = dirc;
            }
            TrialNode<float> *newNode = NULL;

            int index;
            index = searchOrAppendTrialNode ( newNode, newNodeX, newNodeY, newNodeZ );
            newNode->seed_type = Node.seed_type;

            for ( int k = 0; k < 4; k++ ) {
              int newDir = sdir;
              if ( k & 1 ) {
                newDir |= dira;
              }
              if ( k & 2 ) {
                newDir |= dirb;
              }
              ocTree->setElement ( newDir, newNodeX, newNodeY, newNodeZ, coarseLevel + 1 );
              if ( !estimator->checkElementSure ( newDir, newNodeX, newNodeY, newNodeZ, coarseLevel + 1 ) ) {
                cerr << "ERROR in makeTrialNodeActive:: estimator false for element #1\n";
              }
              updateElementToTrialNodes ( newDir, newNodeX, newNodeY, newNodeZ, coarseLevel + 1 );
              newNode->dirLevel[ newDir ] = coarseLevel + 1;

              downwind ( newDir, newNodeX, newNodeY, newNodeZ, ( Step >> 1 ), newNode->val, 0.0f );
            }
            trialNodeHeap.upHeap ( index );
          }
        }
      }
    }
  }

  // now look in every direction and add new Trial Nodes or update the old values.

  for ( int i = 0; i < 8; i++ ) {
    if ( Node.dirLevel[ i ] > GRID_DEPTH ) {
      cerr << "ERROR: Node.dirLevel[ " << i << " ]= " << Node.dirLevel[ i ] << endl;
      Node.dump( );
    }
  }

  for ( int axis = 0; axis < 3; axis++ ) {

    int thisdir = 1 << axis;

    int dira = 1 << ( ( axis + 1 ) % 3 );
    int dirb = 1 << ( ( axis + 2 ) % 3 );

    for ( int ori = 0; ori < 2; ori++ ) {
      int initDir;
      if ( ori == 0 )
        initDir = 0;
      else
        initDir = thisdir;

      int minLevelDirs[ 4 ];
      int numMinLevelDirs = 0;

      int minLevel = GRID_DEPTH + 1;
      for ( int i = 0; i < 4; i++ ) {

        int dir = initDir;
        if ( i & 1 )
          dir |= dira;
        if ( i & 2 )
          dir |= dirb;

        if ( Node.dirLevel[ dir ] < minLevel
             && Node.dirLevel[ dir ] >= 0 ) {
          minLevelDirs[ 0 ] = dir;
          numMinLevelDirs = 1;
          minLevel = Node.dirLevel[ dir ];
        } else if ( Node.dirLevel[ dir ] == minLevel ) {
          minLevelDirs[ numMinLevelDirs++ ] = dir;
        }
      }
      if ( minLevel == GRID_DEPTH + 1 ) {
        cerr << "ERROR minLevel == " << minLevel << endl;
      }
      int step = 1 << ( GRID_DEPTH - minLevel );
      int newNodeX = Node.x(),
                     newNodeY = Node.y(),
                                newNodeZ = Node.z();

      if ( thisdir == TN_DIR_RIGHT ) {
        if ( ori == 0 )
          newNodeX -= step;
        else
          newNodeX += step;
      } else if ( thisdir == TN_DIR_UP ) {
        if ( ori == 0 )
          newNodeY -= step;
        else
          newNodeY += step;
      } else if ( thisdir == TN_DIR_BACK ) {
        if ( ori == 0 )
          newNodeZ -= step;
        else
          newNodeZ += step;
      }

      if ( newNodeX >= 0 && newNodeX <= numX &&
           newNodeY >= 0 && newNodeY <= numY &&
           newNodeZ >= 0 && newNodeZ <= numZ  ) {

        if ( timeField->get ( newNodeX, newNodeY, newNodeZ ) < 0.0 ) {

          TrialNode<float> *newNode = NULL;
          int index;
          index = searchOrAppendTrialNode ( newNode, newNodeX, newNodeY, newNodeZ );

          // at least one of dirlevel from newNode must be initialized

          newNode->seed_type = Node.seed_type;

#ifdef DEBUG
          cerr << "step = " << step << " numMinLevelDirs = " << numMinLevelDirs;
          cerr << " minLevel = " << minLevel << endl;
          cerr << ">>> thisDir = " << thisdir << endl;
#endif
          if ( newNode->faceHangNode == TN_HN_NONE_SEG
               && newNode->edgeHangNode == TN_HN_NONE_SEG ) {
            for ( int j = 0; j < numMinLevelDirs; j++ ) {
              int oDir = 0;
              if ( ori == 0 ) {
                oDir = minLevelDirs[ j ] | thisdir;
              } else {
                oDir = minLevelDirs[ j ] & ( ~thisdir );
              }
              newNode->dirLevel[ oDir ] = minLevel;
              ocTree->setElement ( oDir, newNodeX, newNodeY, newNodeZ, minLevel );
#ifdef DEBUG
              cerr << "Setting element at ( " << newNodeX << " ,"
              << newNodeY << ", " << newNodeZ << ") in dir ";
              TrialNode<float>::parsedir ( oDir );
              cerr << endl;
#endif
              updateElementToTrialNodes ( oDir, newNodeX, newNodeY, newNodeZ, minLevel );

              downwind ( oDir, newNodeX, newNodeY, newNodeZ, step, newNode->val, 0. );
            }

            trialNodeHeap.upHeap ( index );
          }
#ifdef DEBUG
          if ( !newNode->consistencyCheck ( GRID_DEPTH, *estimator, *this ) ) {
            cerr << "Failure #3";
          }
#endif
        }
      }
    }
  }
#ifdef DEBUG
  if ( !Node.consistencyCheck ( GRID_DEPTH, *estimator, *this ) ) {
    cerr << "Failure #4\n";
  }
#endif
}

void Eikonal3d::localHangNodeGridUpdate ( TrialNode<float> &Node ) {
  if ( Node.faceHangNode == TN_HN_NONE_SEG ) {
#ifdef DEBUG
    cerr << "Eikonal3d::localHangNodeGridUpdate: nothing to do... \n";
#endif
    return;
  }

  unsigned char mask = tn_hn_nb_cell_dirs[ Node.faceHangNode ];

  int level = Node.hangNodeLevel + 1;

  for ( int dir = 0; dir <= TN_DIR_MAX_NUM; dir++ ) {
    if ( mask & ( 1 << dir ) ) {
      ocTree->setElement ( dir, Node.x(), Node.y(), Node.z(), level );
      // TEST
      if ( level < GRID_DEPTH
           && !estimator->checkElementSure ( dir, Node.x(), Node.y(), Node.z(), level ) ) {
        cerr << "ERROR in Eikonal3d::localHangNodeGridUpdate: Estimator false for element.\n";
        return;
      }
      updateElementToTrialNodes ( dir, Node.x(), Node.y(), Node.z(), level );
      Node.dirLevel[ dir ] = level;
    }
  }
}

void Eikonal3d::updateElementToTrialNodes ( int X, int Y, int Z, int Level ) {
  return;
  int FullStep = 1 << ( GRID_DEPTH - Level );
  int HalfStep = FullStep >> 1;
  int dir;

  if ( Level > 0 ) {
    int coarseStep = FullStep << 1;
    int cX = X - ( X % coarseStep );
    int cY = Y - ( Y % coarseStep );
    int cZ = Z - ( Z % coarseStep );

    if ( estimator->checkElementSure ( cX, cY, cZ, Level - 1 ) ) {
      cerr << "ERROR in Eikonal3d::updateElementToTrialNodes: Estimator true for parent element!\n";
    }
  }

  if ( Level < GRID_DEPTH
       && !estimator->checkElementSure ( X, Y, Z, Level ) ) {
    cerr << ">>>>>>>>> ERROR in Eikonal3d::updateElementToTrialNodes: Estimator false for element.\n";
    return;
  }

  for ( dir = 0; dir <= TN_DIR_MAX_NUM; dir++ ) {
    int NewX = X, NewY = Y, NewZ = Z;
    if ( ! ( dir & TN_DIR_BACK ) ) {
      NewZ += FullStep;
    }
    if ( ! ( dir & TN_DIR_UP ) ) {
      NewY += FullStep;
    }
    if ( ! ( dir & TN_DIR_RIGHT ) ) {
      NewX += FullStep;
    }

    int i = indexField->get ( NewX, NewY, NewZ );
    if ( i > 0 ) {
      TrialNodeFloat &Node = trialNodeHeap[ i ];
      if ( Node.x() != NewX
           || Node.y() != NewY
           || Node.z() != NewZ ) {
        cerr << "\nERROR indexField incorrect! i = " << i << "\n";
        Node.dump();
      }
#ifdef DEBUG
      if ( !Node.consistencyCheck ( GRID_DEPTH, *estimator, *this ) ) {
        cerr << "Fail in updateElementToTrialNodes: before #1\n";
        cerr << "qc::Element = " << X << ", " << Y << ", " << Z << " Level = " << Level << endl;
        cerr << "NewX = " << NewX << " NewY = " << NewY << " NewZ = " << NewZ << endl;
        abort( );
        Node.debughalt( );
      }
#endif
      Node.dirLevel[ dir ] = Level;
#ifdef DEBUG
      if ( !Node.consistencyCheck ( GRID_DEPTH, *estimator, *this ) ) {
        cerr << "Fail in updateElementToTrialNodes: #1\n";
        cerr << "qc::Element = " << X << ", " << Y << ", " << Z << " Level = " << Level << endl;
        cerr << "NewX = " << NewX << " NewY = " << NewY << " NewZ = " << NewZ << endl;
        abort( );
      }
#endif
    }
  }

  if ( Level < GRID_DEPTH ) {
    for ( dir = 0; dir <= TN_HN_MAX_NUM; dir++ ) {
      int NodeX = X + tn_hn_offs[ dir ][ 0 ] * HalfStep;
      int NodeY = Y + tn_hn_offs[ dir ][ 1 ] * HalfStep;
      int NodeZ = Z + tn_hn_offs[ dir ][ 2 ] * HalfStep;

      int i = indexField->get ( NodeX, NodeY, NodeZ );
      if ( i > 0 ) {
        TrialNodeFloat &Node = trialNodeHeap[ i ];
#ifdef DEBUG
        if ( !Node.consistencyCheck ( GRID_DEPTH, *estimator, *this ) ) {
          cerr << "Fail in updateElementToTrialNodes: before #2\n";
        }
#endif
        Node.addHangNode ( dir, Level );
        hangNodes.push_back ( Node );
#ifdef DEBUG
        if ( !Node.consistencyCheck ( GRID_DEPTH, *estimator, *this ) ) {
          cerr << "Fail in updateElementToTrialNodes: #2\n";
        }
#endif
      }
    }
  }
}

void Eikonal3d::exportBitMask ( qc::ScalarArray<unsigned char, qc::QC_3D> &Array, float Time ) {

  Array.setAll ( 0 );

  for ( unsigned int level = 0; level <= ocTree->getMaxLevel( ); level++ ) {
    int Step = 1 << ( GRID_DEPTH - level );

    for ( int X = 0; X < numX; X += Step ) {
      for ( int Y = 0; Y < numY; Y += Step ) {
        for ( int Z = 0; Z < numZ; Z += Step ) {
          if ( ocTree->getElement ( X, Y, Z, level ) &&
               checkElementDataSet ( X, Y, Z, level ) ) {

            float fll, flr, ful, fur, bll, blr, bul, bur;

            fll = timeField->get ( X       , Y       , Z        );
            flr = timeField->get ( X + Step, Y       , Z        );
            ful = timeField->get ( X       , Y + Step, Z        );
            fur = timeField->get ( X + Step, Y + Step, Z        );
            bll = timeField->get ( X       , Y       , Z + Step );
            blr = timeField->get ( X + Step, Y       , Z + Step );
            bul = timeField->get ( X       , Y + Step, Z + Step );
            bur = timeField->get ( X + Step, Y + Step, Z + Step );

            for ( int lX = 0; lX < Step; lX++ ) {
              for ( int lY = 0; lY < Step; lY++ ) {
                for ( int lZ = 0; lZ < Step; lZ++ ) {
                  float alphaX = float ( lX ) / float ( Step ),
                                 alphaY = float ( lY ) / float ( Step ),
                                          alphaZ = float ( lZ ) / float ( Step );

                  float ll = ( 1.0f - alphaZ ) * fll + alphaZ * bll;
                  float lr = ( 1.0f - alphaZ ) * flr + alphaZ * blr;
                  float ul = ( 1.0f - alphaZ ) * ful + alphaZ * bul;
                  float ur = ( 1.0f - alphaZ ) * fur + alphaZ * bur;

                  float l = ( 1.0f - alphaY ) * ll + alphaY * ul;
                  float r = ( 1.0f - alphaY ) * lr + alphaY * ur;

                  if ( X + lX < numX && Y + lY < numY && Z + lZ < numZ &&
                       ( ( 1.0 - alphaX ) * l + alphaX * r <= Time ) )
                    Array.set ( X + lX, Y + lY, Z + lZ, 255 );
                }
              }
            }
          }
        }
      }
    }
  }
}


void Eikonal3d::exportToArray3d ( qc::ScalarArray<unsigned char, qc::QC_3D> *Array ) {

  Array->setAll ( 255 );
  float max_rec = 1.0f / timeField->getMaxValue( );

  for ( unsigned int level = 0; level <= ocTree->getMaxLevel( ); level++ ) {
    int Step = 1 << ( GRID_DEPTH - level );

    for ( int X = 0; X < numX; X += Step ) {
      for ( int Y = 0; Y < numY; Y += Step ) {
        for ( int Z = 0; Z < numZ; Z += Step ) {
          if ( ocTree->getElement ( X, Y, Z, level ) &&
               checkElementDataSet ( X, Y, Z, level ) ) {

            float fll, flr, ful, fur, bll, blr, bul, bur;

            fll = timeField->get ( X, Y, Z );
            flr = timeField->get ( X + Step, Y, Z );
            ful = timeField->get ( X, Y + Step, Z );
            fur = timeField->get ( X + Step, Y + Step, Z );
            bll = timeField->get ( X, Y, Z + Step );
            blr = timeField->get ( X + Step, Y, Z + Step );
            bul = timeField->get ( X, Y + Step, Z + Step );
            bur = timeField->get ( X + Step, Y + Step, Z + Step );

            for ( int lX = 0; lX < Step; lX++ ) {
              for ( int lY = 0; lY < Step; lY++ ) {
                for ( int lZ = 0; lZ < Step; lZ++ ) {
                  float alphaX = float ( lX ) / float ( Step ),
                                 alphaY = float ( lY ) / float ( Step ),
                                          alphaZ = float ( lZ ) / float ( Step );

                  float ll = ( 1.0f - alphaZ ) * fll + alphaZ * bll;
                  float lr = ( 1.0f - alphaZ ) * flr + alphaZ * blr;
                  float ul = ( 1.0f - alphaZ ) * ful + alphaZ * bul;
                  float ur = ( 1.0f - alphaZ ) * fur + alphaZ * bur;

                  float l = ( 1.0f - alphaY ) * ll + alphaY * ul;
                  float r = ( 1.0f - alphaY ) * lr + alphaY * ur;

                  if ( X + lX < numX && Y + lY < numY && Z + lZ < numZ )
                    Array->set ( X + lX, Y + lY, Z + lZ, static_cast<unsigned char> ( 255.0f * ( ( 1.0 - alphaX ) * l + alphaX * r ) * max_rec ) );
                }
              }
            }
          }
        }
      }
    }
  }
  vector<TrialNodeFloat>::iterator it;
  for ( it = trialNodeHeap.begin( ); it != trialNodeHeap.end( ); ++it ) {
    TrialNodeFloat &cur = ( *it );
    float val = 255.0f * cur.val * max_rec;
    if ( val > 255.0f ) val = 255.0f;
    Array->set ( cur.x(), cur.y(), cur.z(), static_cast<unsigned char> ( val ) );
  }
}

void Eikonal3d::exportToArray3d ( qc::ScalarArray<float, qc::QC_3D> *Array ) {

  float max = timeField->getMaxValue( );
  Array->setAll ( max );

  for ( unsigned int level = 0; level <= ocTree->getMaxLevel( ); level++ ) {
    int Step = 1 << ( GRID_DEPTH - level );

    for ( int X = 0; X < numX; X += Step ) {
      for ( int Y = 0; Y < numY; Y += Step ) {
        for ( int Z = 0; Z < numZ; Z += Step ) {
          if ( ocTree->getElement ( X, Y, Z, level ) &&
               checkElementDataSet ( X, Y, Z, level ) ) {

            float fll, flr, ful, fur, bll, blr, bul, bur;

            fll = timeField->get ( X, Y, Z );
            flr = timeField->get ( X + Step, Y, Z );
            ful = timeField->get ( X, Y + Step, Z );
            fur = timeField->get ( X + Step, Y + Step, Z );
            bll = timeField->get ( X, Y, Z + Step );
            blr = timeField->get ( X + Step, Y, Z + Step );
            bul = timeField->get ( X, Y + Step, Z + Step );
            bur = timeField->get ( X + Step, Y + Step, Z + Step );

            for ( int lX = 0; lX < Step; lX++ ) {
              for ( int lY = 0; lY < Step; lY++ ) {
                for ( int lZ = 0; lZ < Step; lZ++ ) {
                  float alphaX = float ( lX ) / float ( Step ),
                                 alphaY = float ( lY ) / float ( Step ),
                                          alphaZ = float ( lZ ) / float ( Step );

                  float ll = ( 1.0f - alphaZ ) * fll + alphaZ * bll;
                  float lr = ( 1.0f - alphaZ ) * flr + alphaZ * blr;
                  float ul = ( 1.0f - alphaZ ) * ful + alphaZ * bul;
                  float ur = ( 1.0f - alphaZ ) * fur + alphaZ * bur;

                  float l = ( 1.0f - alphaY ) * ll + alphaY * ul;
                  float r = ( 1.0f - alphaY ) * lr + alphaY * ur;

                  if ( X + lX < numX && Y + lY < numY && Z + lZ < numZ )
                    Array->set ( X + lX, Y + lY, Z + lZ, ( ( 1.0f - alphaX ) * l + alphaX * r ) );
                }
              }
            }
          }
        }
      }
    }
  }
  vector<TrialNodeFloat>::iterator it;
  for ( it = trialNodeHeap.begin( ); it != trialNodeHeap.end( ); ++it ) {
    Array->set ( ( *it ).x(), ( *it ).y(), ( *it ).z(), ( *it ).val );
  }
}

void Eikonal3d::fillTimeField( ) {
  for ( unsigned int level = 0; level <= ocTree->getMaxLevel( ); level++ ) {
    int Step = 1 << ( GRID_DEPTH - level );

    for ( int X = 0; X < numX; X += Step ) {
      for ( int Y = 0; Y < numY; Y += Step ) {
        for ( int Z = 0; Z < numZ; Z += Step ) {
          if ( ocTree->getElement ( X, Y, Z, level ) &&
               checkElementDataSet ( X, Y, Z, level ) ) {

            float fll, flr, ful, fur, bll, blr, bul, bur;

            fll = timeField->get ( X, Y, Z );

            flr = timeField->get ( X + Step, Y, Z );

            ful = timeField->get ( X, Y + Step, Z );

            fur = timeField->get ( X + Step, Y + Step, Z );

            bll = timeField->get ( X, Y, Z + Step );

            blr = timeField->get ( X + Step, Y, Z + Step );

            bul = timeField->get ( X, Y + Step, Z + Step );

            bur = timeField->get ( X + Step, Y + Step, Z + Step );

            for ( int lX = 0; lX < Step; lX++ ) {
              for ( int lY = 0; lY < Step; lY++ ) {
                for ( int lZ = 0; lZ < Step; lZ++ ) {
                  float alphaX = float ( lX ) / float ( Step ),
                                 alphaY = float ( lY ) / float ( Step ),
                                          alphaZ = float ( lZ ) / float ( Step );

                  float ll = ( 1.0f - alphaZ ) * fll + alphaZ * bll;
                  float lr = ( 1.0f - alphaZ ) * flr + alphaZ * blr;
                  float ul = ( 1.0f - alphaZ ) * ful + alphaZ * bul;
                  float ur = ( 1.0f - alphaZ ) * fur + alphaZ * bur;

                  float l = ( 1.0f - alphaY ) * ll + alphaY * ul;
                  float r = ( 1.0f - alphaY ) * lr + alphaY * ur;

                  timeField->set ( X + lX, Y + lY, Z + lZ, ( 1.0f - alphaX ) * l + alphaX * r );
                }
              }
            }
          }
        }
      }
    }
  }
}

void Eikonal3d::generateTimeSliceX ( qc::ScalarArray<float, qc::QC_2D> *TimeSlice, int X ) {
  TimeSlice->setAll ( -1.0f );
  int cX, cY, cZ;
  for ( int Level = 0; Level <= GRID_DEPTH; Level++ ) {
    int Step = 1 << ( GRID_DEPTH - Level );
    cX = X - ( X % Step );
    float alphaX = 1.0f - ( static_cast<float> ( X - cX ) ) / static_cast<float> ( Step );
    for ( int Y = 0; Y < numY; Y += Step ) {
      cY = Y;
      for ( int Z = 0; Z < numZ; Z += Step ) {
        cZ = Z;
        if ( ocTree->getElement ( cX, cY, cZ, Level ) &&
             checkElementDataSet ( cX, cY, cZ, Level ) ) {

          float t11 = timeField->get ( cX, cY, cZ ) * alphaX +
                      timeField->get ( cX + Step, cY, cZ ) * ( 1.0f - alphaX );

          float t12 = timeField->get ( cX, cY, cZ + Step ) * alphaX
                      + timeField->get ( cX + Step, cY, cZ + Step ) * ( 1.0f - alphaX );

          float t21 = timeField->get ( cX, cY + Step, cZ ) * alphaX
                      + timeField->get ( cX + Step, cY + Step, cZ ) * ( 1.0f - alphaX );

          float t22 = timeField->get ( cX, cY + Step, cZ + Step ) * alphaX
                      + timeField->get ( cX + Step, cY + Step, cZ + Step ) * ( 1.0f - alphaX );

          for ( int lY = 0; lY <= Step; lY++ ) {
            float alphaY = 1.0f - ( static_cast<float> ( lY ) ) / ( static_cast<float> ( Step ) );
            for ( int lZ = 0; lZ <= Step; lZ++ ) {
              float alphaZ = 1.0f - ( static_cast<float> ( lZ ) ) / ( static_cast<float> ( Step ) );
              float v = ( alphaY * alphaZ * t11 + alphaY * ( 1.0f - alphaZ ) * t12
                          + ( 1.0f - alphaY ) * alphaZ * t21 + ( 1.0f - alphaY ) * ( 1.0f - alphaZ ) * t22 );
              TimeSlice->set ( Y + lY, Z + lZ, v );
            }
          }
        }
      }
    }
  }
}

void Eikonal3d::generateTimeSliceY ( qc::ScalarArray<float, qc::QC_2D> *TimeSlice, int Y ) {
  TimeSlice->setAll ( -1.0f );
  int cX, cY, cZ;
  for ( int Level = 0; Level <= GRID_DEPTH; Level++ ) {
    int Step = 1 << ( GRID_DEPTH - Level );
    cY = Y - ( Y % Step );
    float alphaY = 1.0f - ( static_cast<float> ( Y - cY ) ) / static_cast<float> ( Step );
    for ( int X = 0; X < numX; X += Step ) {
      cX = X - ( X % Step );
      for ( int Z = 0; Z < numZ; Z += Step ) {
        cZ = Z - ( Z % Step );
        if ( ocTree->getElement ( cX, cY, cZ, Level ) &&
             checkElementDataSet ( cX, cY, cZ, Level ) ) {

          float t11 = timeField->get ( cX, cY, cZ ) * alphaY +
                      timeField->get ( cX, cY + Step, cZ ) * ( 1.0f - alphaY );

          float t12 = timeField->get ( cX, cY, cZ + Step ) * alphaY +
                      timeField->get ( cX, cY + Step, cZ + Step ) * ( 1.0f - alphaY );

          float t21 = timeField->get ( cX + Step, cY, cZ ) * alphaY +
                      timeField->get ( cX + Step, cY + Step, cZ ) * ( 1.0f - alphaY );

          float t22 = timeField->get ( cX + Step, cY, cZ + Step ) * alphaY +
                      timeField->get ( cX + Step, cY + Step, cZ + Step ) * ( 1.0f - alphaY );

          for ( int lX = 0; lX <= Step; lX++ ) {
            float alphaX = 1.0f - static_cast<float> ( lX ) / static_cast<float> ( Step );
            for ( int lZ = 0; lZ <= Step; lZ++ ) {
              float alphaZ = 1.0f - static_cast<float> ( lZ ) / static_cast<float> ( Step );
              float v = ( alphaX * alphaZ * t11 + alphaX * ( 1.0f - alphaZ ) * t12
                          + ( 1.0f - alphaX ) * alphaZ * t21 + ( 1.0f - alphaX ) * ( 1.0f - alphaZ ) * t22 );
              TimeSlice->set ( X + lX, Z + lZ, v );
            }
          }
        }
      }
    }
  }
}
void Eikonal3d::generateTimeSliceZ ( qc::ScalarArray<float, qc::QC_2D> *TimeSlice, int Z ) {
  TimeSlice->setAll ( -1.0f );
  int cX, cY, cZ;
  for ( int Level = 0; Level <= GRID_DEPTH; Level++ ) {
    int Step = 1 << ( GRID_DEPTH - Level );
    cZ = Z - ( Z % Step );
    float alphaZ = 1.0f - static_cast<float> ( Z - cZ ) / static_cast<float> ( Step );
    for ( int X = 0; X < numX; X += Step ) {
      cX = X - ( X % Step );
      for ( int Y = 0; Y < numZ; Y += Step ) {
        cY = Y - ( Y % Step );
        if ( ocTree->getElement ( cX, cY, cZ, Level ) &&
             checkElementDataSet ( cX, cY, cZ, Level ) ) {

          float t11 = timeField->get ( cX, cY, cZ ) * alphaZ +
                      timeField->get ( cX, cY, cZ + Step ) * ( 1.0f - alphaZ );

          float t12 = timeField->get ( cX, cY + Step, cZ ) * alphaZ +
                      timeField->get ( cX, cY + Step, cZ + Step ) * ( 1.0f - alphaZ );

          float t21 = timeField->get ( cX + Step, cY, cZ ) * alphaZ +
                      timeField->get ( cX + Step, cY, cZ + Step ) * ( 1.0f - alphaZ );

          float t22 = timeField->get ( cX + Step, cY + Step, cZ ) * alphaZ +
                      timeField->get ( cX + Step, cY + Step, cZ + Step ) * ( 1.0f - alphaZ );

          for ( int lX = 0; lX <= Step; lX++ ) {
            float alphaX = 1.0f - static_cast<float> ( lX ) / static_cast<float> ( Step );
            for ( int lY = 0; lY <= Step; lY++ ) {
              float alphaY = 1.0f - static_cast<float> ( lY ) / static_cast<float> ( Step );
              float v = ( alphaX * alphaY * t11 + alphaX * ( 1.0f - alphaY ) * t12
                          + ( 1.0f - alphaX ) * alphaY * t21 + ( 1.0f - alphaX ) * ( 1.0f - alphaY ) * t22 );
              TimeSlice->set ( X + lX, Y + lY, v );
            }
          }
        }
      }
    }
  }
}


void Eikonal3d::cont_dist( ) {
  vector<TrialNodeFloat>::iterator it;

  float min = trialNodeHeap[ 0 ].val, max = trialNodeHeap[ 0 ].val;
  for ( it = trialNodeHeap.begin( ); it != trialNodeHeap.end( ); ++it ) {
    if ( ( *it ).val > max ) max = ( *it ).val;
  }

  for ( it = trialNodeHeap.begin( ); it != trialNodeHeap.end( ); ++it ) {
    ( *it ).val = ( ( *it ).val - min ) / ( max - min ) * 0.01f + min;
  }
}

void Eikonal3d::accept ( float T ) {
  //! @todo: implement this more efficiently

#ifdef SUCCESSIVE_SEGMENT
  if ( !tmpTimeField )
    tmpTimeField = new qc::ScalarArray<float, qc::QC_3D> ( numX + 1, numY + 1, numZ + 1 );
#else
  cerr << "This function has been turned off.\n";
  return;
#endif

  for ( unsigned int level = 0; level <= ocTree->getMaxLevel( ); level++ ) {
    int Step = 1 << ( GRID_DEPTH - level );

    for ( int X = 0; X < numX; X += Step ) {
      for ( int Y = 0; Y < numY; Y += Step ) {
        for ( int Z = 0; Z < numZ; Z += Step ) {
          if ( ocTree->getElement ( X, Y, Z, level ) &&
               checkElementDataSet ( X, Y, Z, level ) ) {

            float fll, flr, ful, fur, bll, blr, bul, bur;

            fll = timeField->get ( X, Y, Z );
            flr = timeField->get ( X + Step, Y, Z );
            ful = timeField->get ( X, Y + Step, Z );
            fur = timeField->get ( X + Step, Y + Step, Z );
            bll = timeField->get ( X, Y, Z + Step );
            blr = timeField->get ( X + Step, Y, Z + Step );
            bul = timeField->get ( X, Y + Step, Z + Step );
            bur = timeField->get ( X + Step, Y + Step, Z + Step );

            for ( int lX = 0; lX < Step; lX++ ) {
              for ( int lY = 0; lY < Step; lY++ ) {
                for ( int lZ = 0; lZ < Step; lZ++ ) {
                  const float alphaX = float ( lX ) / float ( Step ),
                                       alphaY = float ( lY ) / float ( Step ),
                                                alphaZ = float ( lZ ) / float ( Step );

                  const float ll = ( 1.0f - alphaZ ) * fll + alphaZ * bll;
                  const float lr = ( 1.0f - alphaZ ) * flr + alphaZ * blr;
                  const float ul = ( 1.0f - alphaZ ) * ful + alphaZ * bul;

                  const float ur = ( 1.0f - alphaZ ) * fur + alphaZ * bur;

                  const float l = ( 1.0f - alphaY ) * ll + alphaY * ul;
                  const float r = ( 1.0f - alphaY ) * lr + alphaY * ur;

                  const float v = ( ( 1.0f - alphaX ) * l + alphaX * r );

                  if ( v <= T ) {
                    tagField->set ( X + lX, Y + lY, Z + lZ, 255 );
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  exportToArray3d ( tmpTimeField );
  for ( int i = 0; i < numX*numY*numZ; i++ ) {
    finalTimeField->set ( i, aol::Min ( tmpTimeField->get ( i ) - T, finalTimeField->get ( i ) ) );
  }

#ifdef SUCCESSIVE_SEGMENT
  if ( tmpTimeField ) delete tmpTimeField;
  tmpTimeField = NULL;
#endif
}

void Eikonal3d::erase ( float T ) {

#ifdef SUCCESSIVE_SEGMENT
  if ( !tmpTimeField )
    tmpTimeField = new qc::ScalarArray<float, qc::QC_3D> ( numX + 1, numY + 1, numZ + 1 );
#else
  cerr << "This function has been turned off.\n";
  return;
#endif

  for ( unsigned int level = 0; level <= ocTree->getMaxLevel( ); level++ ) {
    int Step = 1 << ( GRID_DEPTH - level );

    for ( int X = 0; X < numX; X += Step ) {
      for ( int Y = 0; Y < numY; Y += Step ) {
        for ( int Z = 0; Z < numZ; Z += Step ) {
          if ( ocTree->getElement ( X, Y, Z, level ) &&
               checkElementDataSet ( X, Y, Z, level ) ) {

            float fll, flr, ful, fur, bll, blr, bul, bur;

            fll = timeField->get ( X, Y, Z );
            flr = timeField->get ( X + Step, Y, Z );
            ful = timeField->get ( X, Y + Step, Z );
            fur = timeField->get ( X + Step, Y + Step, Z );
            bll = timeField->get ( X, Y, Z + Step );
            blr = timeField->get ( X + Step, Y, Z + Step );
            bul = timeField->get ( X, Y + Step, Z + Step );
            bur = timeField->get ( X + Step, Y + Step, Z + Step );

            for ( int lX = 0; lX < Step; lX++ ) {
              for ( int lY = 0; lY < Step; lY++ ) {
                for ( int lZ = 0; lZ < Step; lZ++ ) {
                  const float alphaX = float ( lX ) / float ( Step ),
                                       alphaY = float ( lY ) / float ( Step ),
                                                alphaZ = float ( lZ ) / float ( Step );

                  const float ll = ( 1.0f - alphaZ ) * fll + alphaZ * bll;
                  const float lr = ( 1.0f - alphaZ ) * flr + alphaZ * blr;
                  const float ul = ( 1.0f - alphaZ ) * ful + alphaZ * bul;
                  const float ur = ( 1.0f - alphaZ ) * fur + alphaZ * bur;

                  const float l = ( 1.0f - alphaY ) * ll + alphaY * ul;
                  const float r = ( 1.0f - alphaY ) * lr + alphaY * ur;

                  const float v = ( ( 1.0f - alphaX ) * l + alphaX * r );

                  if ( v <= T ) {
                    tagField->set ( X + lX, Y + lY, Z + lZ, 0 );
                  }

                }
              }
            }
          }
        }
      }
    }
  }

  exportToArray3d ( tmpTimeField );
  for ( int i = 0; i < numX*numY*numZ; i++ ) {
    finalTimeField->set ( i, aol::Max ( T - tmpTimeField->get ( i ), finalTimeField->get ( i ) ) );
  }

#ifdef SUCCESSIVE_SEGMENT
  if ( tmpTimeField ) delete tmpTimeField;
  tmpTimeField = NULL;
#endif

}

// originates from file SegEstimator3d.cpp:

void SegEstimator3d::makeGradients( ) {
  qc::ScalarArray<unsigned char, qc::QC_3D> tempArray ( *this );

  float dx, dy, dz;

  cerr << "- computing gradients..\n";
  for ( int X = 0; X < numX; X++ ) {
    cerr << int ( float ( X ) / float ( numX - 1 ) *100.0f ) << "% done.   \r";
    for ( int Y = 0; Y < numY; Y++ )
      for ( int Z = 0; Z < numZ; Z++ ) {
        if ( X == 0 ) {
          dx = 2.0f * ( tempArray.get ( X + 1, Y, Z )
                        - tempArray.get ( X, Y, Z ) );
        } else if ( X == numX - 1 ) {
          dx = 2.0f * ( tempArray.get ( X, Y, Z )
                        - tempArray.get ( X - 1, Y, Z ) );
        } else {
          dx = static_cast<float> ( ( tempArray.get ( X + 1, Y, Z )
                                      - tempArray.get ( X - 1, Y, Z ) ) );
        }

        if ( Y == 0 ) {
          dy = 2.0f * ( tempArray.get ( X, 1, Z )
                        - tempArray.get ( X, 0, Z ) );
        } else if ( Y == numY - 1 ) {
          dy = 2.0f * ( tempArray.get ( X, Y, Z )
                        - tempArray.get ( X, Y - 1, Z ) );
        } else {
          dy = static_cast<float> ( ( tempArray.get ( X, Y + 1, Z )
                                      - tempArray.get ( X, Y - 1, Z ) ) );
        }

        if ( Z == 0 ) {
          dz = 2.0f * ( tempArray.get ( X, Y, 1 )
                        - tempArray.get ( X, Y, 0 ) );
        } else if ( Z == numZ - 1 ) {
          dz = 2.0f * ( tempArray.get ( X, Y, Z )
                        - tempArray.get ( X, Y, Z - 1 ) );
        } else {
          dz = static_cast<float> ( ( tempArray.get ( X, Y, Z + 1 )
                                      - tempArray.get ( X, Y, Z - 1 ) ) );
        }

        unsigned char grad = static_cast<unsigned char> ( sqrt ( dx * dx + dy * dy + dz * dz ) );
        set ( X, Y, Z, grad );
      }
  }
  cerr << "\n";
}

void SegEstimator3d::retrieveOriginalImage( ) {
  cerr << "retrieveOriginalImage called" << endl;
  if ( origImage != NULL ) {
    if ( numX == origImage->getNumX()
         && numY == origImage->getNumY()
         && numZ == origImage->getNumZ() ) {
      int m = numX * numY * numZ;
      for ( int i = 0; i < m; i++ ) {
        set ( i, static_cast <unsigned char> ( 255.0f * origImage->get ( i ) ) );
      }
    } else {
      for ( int i = 0; i < numX; i++ )
        for ( int j = 0; j < numY; j++ )
          for ( int k = 0; k < numZ; k++ ) {
            set ( i, j, k, static_cast <unsigned char> ( 255.0f * origImage->get ( i, j, k ) ) );
          }
    }
  }
  return;
}

}   // end namespace eik

