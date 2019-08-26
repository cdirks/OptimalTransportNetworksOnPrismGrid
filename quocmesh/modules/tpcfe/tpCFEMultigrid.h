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

#ifndef __TPCFEMULTIGRID_H
#define __TPCFEMULTIGRID_H

#include <tpCFEElement.h>
#include <tpCFEGrid.h>
#include <tpCFEMatrices.h>

#include <multigrid.h>

#include <preconditioner.h>
#include <smoother.h>

#include <arrayExtensions.h>
#include <matrixInverse.h>


/* Contents:
 * ---------
 * structure coarsening
 * infrastructure for slope bookkeeping
 * prolongation
 * restriction
 * operator (matrix) coarsening
 * multigrid solver
 */

namespace tpcfe {

// === some utility functions

inline bool CFEMG_is_inside ( const int x, const int y, const int z, const int n ) {
  return ( ( x >= 0 ) && ( x < n ) && ( y >= 0 ) && ( y < n ) && ( z >= 0 ) && ( z < n ) );
}

inline bool CFEMG_is_inside ( const qc::CoordType& pos, const int n ) {
  return ( ( pos[0] >= 0 ) && ( pos[0] < n ) && ( pos[1] >= 0 ) && ( pos[1] < n ) && ( pos[2] >= 0 ) && ( pos[2] < n ) );
}

 inline bool CFEMG_is_inside ( const qc::CoordType& pos, const aol::Vec3<int>& n ) {
  return ( ( pos[0] >= 0 ) && ( pos[0] < n[0] ) && ( pos[1] >= 0 ) && ( pos[1] < n[1] ) && ( pos[2] >= 0 ) && ( pos[2] < n[2] ) );
}

// \todo: if grid known, this can be checked on the grid; adapt to rect domain case

// === Structure Coarsening ===

//! Coarsen structure of domainNodes from level fine to level fine-1.
//! Node in coarse grid is domainNode if parent of fine domainNode
//! This probably only makes sense for CFE_CD
template< typename GridType >
void coarsenDomainNodeMask ( const GridType &fineGrid, GridType &coarseGrid ) {

  typedef typename GridType::RealType RealType;

  if ( fineGrid.getGridDepth() != coarseGrid.getGridDepth() + 1 )
    throw aol::Exception ( "Incompatible grid depths for coarsenDomainNodeMask", __FILE__, __LINE__ );

  qc::BitArray<qc::QC_3D> coarseDomainNodeMask ( qc::GridSize<qc::QC_3D>::createFrom ( coarseGrid ) );
  coarseDomainNodeMask.setAll ( false );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( coarseGrid ); bit.notAtEnd(); ++bit ) {
    bool res = false;
    for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
      const qc::CoordType NBPosf = static_cast<short> ( 2 ) * ( *bit ) + CFELookup<RealType>::_tetNbOffsets[i];
      if ( fineGrid.isAdmissibleNode ( NBPosf ) ) {
        res |= fineGrid.isDomainNode ( aol::Vec3<short> ( NBPosf ) );
      }
    }
    coarseDomainNodeMask.set ( *bit, res );
  }

  coarseGrid.setDomainNodeMask ( coarseDomainNodeMask );
}

/** Coarsen structure of Dirichlet nodes from level fine to level fine-1.
 *  Node in coarse grid is Dirichlet node if corresponding node in fine grid is Dirichlet node.
 *  This does not respect the  domainNodeMask, so it is possible to specify non-domain Dirichlet nodes on the finest grid.
 *  If only domainNodes are Dirichlet nodes, the Dirichlet boundary may become smaller on the coarsened grids, if there are non-domain Dirichlet nodes, the Dirichlet boundary may become bigger on the coarsened grids.
 *  This probably makes sense for all constraint types
 */

template< typename GridType >
void coarsenDirichletMask ( const GridType &fineGrid, GridType &coarseGrid ) {

  if ( fineGrid.getGridDepth() != coarseGrid.getGridDepth() + 1 )
    throw aol::Exception ( "Incompatible grid depths for coarsenDirichletMask", __FILE__, __LINE__ );

  qc::BitArray<qc::QC_3D> coarseDirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( coarseGrid ) );
  coarseDirichletMask.setAll ( false );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( coarseGrid ); bit.notAtEnd(); ++bit ) {
    // no check for admissibility of node necessary
    if ( fineGrid.isDirichletNode ( static_cast<short> ( 2 ) * ( *bit ) ) ) {
      coarseDirichletMask.set ( *bit, true );
    }
  }

  coarseGrid.setDirichletMask ( coarseDirichletMask );

}


//! Coarsen structure of degrees of freedom
//! Node in coarse grid is DOF if it is parent of DOF and not Dirichlet node. Being parent of DOF node implies being parent of domainNode.
template< typename GridType >
void coarsenDOFMask ( const GridType &fineGrid, GridType &coarseGrid ) {

  typedef typename GridType::RealType RealType;

  qc::BitArray<qc::QC_3D> coarseDOFMask ( qc::GridSize<qc::QC_3D>::createFrom ( coarseGrid ) );
  coarseDOFMask.setAll ( false );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( coarseGrid ); bit.notAtEnd(); ++bit ) {
    if ( !coarseGrid.isDirichletNode ( *bit ) ) {
      bool res = false;
      for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
        const qc::CoordType NBPosf = static_cast<short> ( 2 ) * ( *bit ) + CFELookup<RealType>::_tetNbOffsets[i];
        if ( fineGrid.isAdmissibleNode ( NBPosf ) ) {
          res |= fineGrid.isDOFNode ( aol::Vec3<short> ( NBPosf ) );
        }
      }

      coarseDOFMask.set ( *bit, res );

    }
  }

  coarseGrid.setDOFMask ( coarseDOFMask );
}

// === Slope (of Basis Functions) Management ===

//! Basis class for template specialization
template< typename GridType >
class SlopeInterface;


template< typename GridType >
class SlopeInterfaceBase {
friend class SlopeInterface< CFEGrid< typename GridType::RealType, GridType::CT > >; // want to call protected members from child class on foreign instance of parent class

protected:
  typedef uint64_t                      EdgeIndexType;
  typedef uint64_t                      CoarsenIndexType;
  typedef typename GridType::RealType   RealType;
  typedef aol::Vec2< RealType >         SlopesType;

protected:
  qc::AArray< std::map< EdgeIndexType, SlopesType >, qc::QC_3D >  _slopeStorage;          //!< for each node, store all nonzero slopes, even standard ones, identified by the edge (in sorted global node indices) on which they are computed
  qc::AArray< std::map< qc::CoordType, RealType >, qc::QC_3D >    _coarseningWeights;     //!< store all coarsening weights, even standard ones
  //  std::map< CoarsenIndexType, RealType >                          _coarseningWeights;
  qc::FastILexMapper<qc::QC_3D>                                   _iMap;                  //!< index mapper for the (fine) grid
  std::map< qc::CoordType, short >                                _neighborOffsetToIndex; //!< maps relative position of neighbor to index of neighbor

public:
  SlopeInterfaceBase ( ) {
    setUpLookup();
  }

  explicit SlopeInterfaceBase ( const GridType& Grid ) : _slopeStorage ( Grid ), _coarseningWeights ( Grid ), _iMap ( Grid ), _neighborOffsetToIndex( ) {
    setUpLookup();
  }

  virtual ~SlopeInterfaceBase ( ) {
  }

  virtual void determineSlopes ( const GridType& /*Grid*/ ) = 0;
  virtual void determineCoarseningWeights ( const GridType& /*Grid*/ ) = 0;
  virtual void coarsenSlopeInterface ( SlopeInterfaceBase<GridType> & /*coarseSlI*/, const GridType &/*fineGrid*/, const GridType &/*coarseGrid*/ ) = 0;

  void determineSlopesAndCoarseningWeights ( const GridType &Grid ) {
    determineSlopes ( Grid );
    determineCoarseningWeights ( Grid );
  }

  //! get slope of the basis function at BFPos on the edge [init point, term point] (in that direction) at initPoint
  inline RealType getOutSlope ( const qc::CoordType &BFPos, const qc::CoordType &initPoint, const qc::CoordType &termPoint ) const {
    // cerr << "getOutSlope: " << BFPos << " " << initPoint << " " << termPoint << " " << getSlope ( BFPos, initPoint, termPoint, 0 ) << endl;
    return ( getSlope ( BFPos, initPoint, termPoint, 0 ) );
  }

  //! get slope of the basis function at BFPos on the edge [init point, term point] (in that direction) at termPoint
  inline RealType getArrSlope ( const qc::CoordType &BFPos, const qc::CoordType &initPoint, const qc::CoordType &termPoint ) const {
    // cerr << "getArrSlope: " << BFPos << " " << initPoint << " " << termPoint << " " << getSlope ( BFPos, initPoint, termPoint, 1 ) << endl;
    return ( getSlope ( BFPos, initPoint, termPoint, 1 ) );
  }

  //! \todo think about what happens for nodes that are present on the coarse grid: storage of 1/2 not really necessary ...
  RealType getCoarseningWeight ( const qc::CoordType &fNPos, const qc::CoordType &cPPos ) const {
    const typename std::map< qc::CoordType, RealType>::const_iterator mit = _coarseningWeights.getRef ( cPPos ).find ( fNPos );
    if ( mit == _coarseningWeights.getRef ( cPPos ).end() )
      return ( aol::NumberTrait<RealType>::zero );
    else
      return ( mit->second );
  }

  const std::map< qc::CoordType, RealType >& getCoarseningWeightMapRef ( const qc::CoordType &cPPos ) const {
    return ( _coarseningWeights.getRef ( cPPos ) );
  }

  void dumpSlopes ( ostream &os = cout ) const;

  void dumpCoarseningWeights ( ostream &os = cout ) const;

  void dump ( ostream &os = cout ) const {
    dumpSlopes ( os );
    dumpCoarseningWeights ( os );
  }

  void checkSumOfCoarseningWeights ( ostream &os = cout ) const;

  //! call this method only if you know what you are doing and if you wish to use standard coarsening (1, 0.5) everywhere (which is incorrect for both jumping coefficients and complicated domains)
  //! compile with --expert :-)
  virtual void storeAllStandardSlopesOnly ( const GridType & ) = 0;

protected:
  void setUpLookup ( ) {
    for ( short i = 0; i < NUM_NEIGHBORS; ++i ) {
      _neighborOffsetToIndex[ CFETopoLookup::_tetNbOffsets[i] ] = i;
    }
  }

  inline EdgeIndexType getEdgeIndex ( const int node1, const int node2 ) const {
#ifdef DEBUG
    if ( node1 >= node2 )
      throw aol::Exception ( "SlopeStorage::getEdgeIndex: wrong orientation of edge", __FILE__, __LINE__ );
#endif
    const EdgeIndexType result = ( ( static_cast<EdgeIndexType> ( node1 ) << 32 ) + static_cast<EdgeIndexType> ( node2 ) );
    return ( result );
  }

  inline void splitEdgeIndex ( const EdgeIndexType edgeIndex, int &node1, int &node2 ) const {
    node1 = static_cast<int> ( edgeIndex >> 32 );
    node2 = static_cast<int> ( edgeIndex );
#ifdef DEBUG
    if ( node1 >= node2 )
      throw aol::Exception ( "SlopeStorage::splitEdgeIndex: wrong orientation of edge", __FILE__, __LINE__ );
#endif
  }

  inline CoarsenIndexType getCoarsenIndex ( const int fNIndex, const int cPIndex ) const {
    const CoarsenIndexType result = ( ( static_cast<CoarsenIndexType> ( fNIndex ) << 32 ) + static_cast<CoarsenIndexType> ( cPIndex ) );
    return ( result );
  }

  inline void splitCoarsenIndex ( const CoarsenIndexType restrictIndex, int &fineIndex, int &coarseIndex ) const {
    fineIndex = static_cast<int> ( restrictIndex >> 32 );
    coarseIndex = static_cast<int> ( restrictIndex );
  }

  //! Saves slopes of basis function at BFPos on edge [initPoint, termPoint]: outSlope at initPoint, arrSlope at termPoint. Internal storage is with respect to positive orientation.
  void setSlopes ( const qc::CoordType &BFPos, const qc::CoordType &initPoint, const qc::CoordType &termPoint, const RealType outSlope, const RealType arrSlope ) {
    const int initIndex = _iMap ( initPoint ), termIndex = _iMap ( termPoint );
    if ( initIndex < termIndex ) {
      // already have positive orientation
      _slopeStorage.getRef ( BFPos ) [ getEdgeIndex ( initIndex, termIndex ) ] = SlopesType ( outSlope, arrSlope );
    } else {
      // need to invert edge, flip and change sign of slopes
      _slopeStorage.getRef ( BFPos ) [ getEdgeIndex ( termIndex, initIndex ) ] = SlopesType ( - arrSlope, - outSlope );
    }
  }

  inline RealType getSlope ( const qc::CoordType &BFPos, const qc::CoordType &initPoint, const qc::CoordType &termPoint, const int which ) const {
    const int initIndex = _iMap ( initPoint ), termIndex = _iMap ( termPoint );
    if ( initIndex < termIndex ) {
      // asked for edge in positive orientation
      typename std::map< EdgeIndexType, SlopesType >::const_iterator mappos = _slopeStorage.getRef ( BFPos ).find ( getEdgeIndex ( initIndex, termIndex ) );
      return ( mappos != _slopeStorage.getRef ( BFPos ).end() ? + mappos->second[   which   ] : aol::NumberTrait<RealType>::zero );
    } else {
      // asked for edge in negative orientation
      typename std::map< EdgeIndexType, SlopesType >::const_iterator mappos = _slopeStorage.getRef ( BFPos ).find ( getEdgeIndex ( termIndex, initIndex ) );
      return ( mappos != _slopeStorage.getRef ( BFPos ).end() ? - mappos->second[ 1 - which ] : aol::NumberTrait<RealType>::zero );
    }
  }

  void setCoarseningWeight ( const qc::CoordType fNPos, const qc::CoordType cPPos, const RealType cWeight ) {
    _coarseningWeights.getRef ( cPPos ) [ fNPos ] = cWeight;
  }

  void removeStoredData ( ) {
    _slopeStorage.resetEntries();
    _coarseningWeights.resetEntries();
  }

  // end class SlopeInterfaceBase
};


//! Basis class for slope interface (compatibility only) classes for standard coarsening.
template< typename GridType >
class StandardSlopeInterface: public SlopeInterfaceBase< GridType > {
public:
  StandardSlopeInterface ( ) : SlopeInterfaceBase<GridType> ( ) {}

  StandardSlopeInterface ( const GridType& ) : SlopeInterfaceBase<GridType> ( ) {
    throw aol::Exception ( "StandardSlopeInterface for CFE_CD not useful", __FILE__, __LINE__ );
  }

  void determineSlopes ( const GridType& ) {
    throw aol::Exception ( "StandardSlopeInterface for CFE_CD not useful", __FILE__, __LINE__ );
  }

  virtual void determineCoarseningWeights ( const GridType& /*Grid*/ ) {
    throw aol::Exception ( "StandardSlopeInterface for CFE_CD not useful", __FILE__, __LINE__ );
  }

  virtual void coarsenSlopeInterface ( SlopeInterfaceBase<GridType> & /*coarseSlI*/, const GridType &/*fineGrid*/, const GridType &/*coarseGrid*/ ) {
    throw aol::Exception ( "StandardSlopeInterface for CFE_CD not useful", __FILE__, __LINE__ );
  }

  virtual void storeAllStandardSlopesOnly ( const GridType & ) {
    throw aol::Exception ( "StandardSlopeInterface for CFE_CD not useful", __FILE__, __LINE__ );
  }
};


//! Compatibility slope interface for CFE_CD
template< typename RealType, typename NodalCoeffType >
class SlopeInterface< CFEGrid< RealType, CFE_CD, NodalCoeffType > > : public StandardSlopeInterface< CFEGrid< RealType, CFE_CD, NodalCoeffType > > {
  typedef CFEGrid< RealType, CFE_CD, NodalCoeffType > GridType;

public:
  SlopeInterface ( ) : StandardSlopeInterface<GridType> () {}

  SlopeInterface ( const GridType &Grid ) : StandardSlopeInterface<GridType> ( Grid ) {}
};


//! Compatibility slope interface for CFE_TPOSELAST
template< typename RealType >
class SlopeInterface< CFEGrid< RealType, CFE_TPOSELAST, IsotropicElasticityCoefficient<RealType> > > : public StandardSlopeInterface< CFEGrid< RealType, CFE_TPOSELAST, IsotropicElasticityCoefficient<RealType> > > {
  typedef CFEGrid< RealType, CFE_TPOSELAST, IsotropicElasticityCoefficient<RealType> > GridType;

public:
  SlopeInterface ( ) : StandardSlopeInterface<GridType> () {}

  SlopeInterface ( const GridType &Grid ) : StandardSlopeInterface<GridType> ( Grid ) {}
};

template< typename RealType >
class SlopeInterface< CFEGrid< RealType, CFE_TPOSELAST, VoigtElasticityCoefficient<RealType> > > : public StandardSlopeInterface< CFEGrid< RealType, CFE_TPOSELAST, VoigtElasticityCoefficient<RealType> > > {
  typedef CFEGrid< RealType, CFE_TPOSELAST, VoigtElasticityCoefficient<RealType> > GridType;

public:
  SlopeInterface ( ) : StandardSlopeInterface<GridType> () {}

  SlopeInterface ( const GridType &Grid ) : StandardSlopeInterface<GridType> ( Grid ) {}
};


template< typename RealType >
class SlopeInterface< CFEGrid< RealType, CFE_LIEHR > > : public SlopeInterfaceBase< CFEGrid< RealType, CFE_LIEHR > > {
protected:
  typedef CFEGrid< RealType, CFE_LIEHR >                            GridType;
  typedef typename SlopeInterfaceBase<GridType>::CoarsenIndexType   CoarsenIndexType;

public:
  //! Standard constructor, probably not useful. Make sure you properly set _iMap before using it.
  SlopeInterface ( ) : SlopeInterfaceBase< GridType > ( ) {}

  explicit SlopeInterface ( const GridType &Grid ) : SlopeInterfaceBase<GridType> ( Grid ) {}

  void determineSlopes ( const GridType &Grid ) ;

  void determineCoarseningWeights ( const GridType &Grid );

  void coarsenSlopeInterface ( SlopeInterfaceBase<GridType> &coarseSlI, const GridType &fineGrid, const GridType &coarseGrid ) ;

  void removeData ( ) {
    this->removeStoredData();
  }

  virtual void storeAllStandardSlopesOnly ( const GridType & ) {
    throw aol::UnimplementedCodeException ( "storeAllStandardSlopesOnly for CFE_LIEHR not implemented yet.", __FILE__, __LINE__ );
  }

  // end class SlopeInterface<RealType, CFE_LIEHR>
};


template< typename RealType >
class SlopeInterface< CFEGrid< RealType, CFE_TPOS > > : public SlopeInterfaceBase< CFEGrid< RealType, CFE_TPOS > > {

protected:
  class RelevantEdgeIdentifier {
  protected:
    short _offsetI, _directionI; // can use unsigned char, then need to cast

  public:
    RelevantEdgeIdentifier ( ) : _offsetI ( 0 ), _directionI ( 0 ) {}

    // create from neighbor index: offset is index, direction is negative offset
    RelevantEdgeIdentifier ( const short neighborI ) : _offsetI ( NEIGHBOR_SELF ), _directionI ( neighborI ) {
      if ( _directionI == NEIGHBOR_SELF )
        throw aol::Exception ( "RelevantEdgeIdentifier ( const short neighborI ): neighborI must not be tpcfe::NEIGHBOR_SELF", __FILE__, __LINE__ );
    }

    // create from offset to neighbor and direction from neighbor
    RelevantEdgeIdentifier ( const short offsetI, const short directionI ) : _offsetI ( offsetI ), _directionI ( directionI ) {
      if ( _directionI == NEIGHBOR_SELF )
        throw aol::Exception ( "RelevantEdgeIdentifier: _directionI must not be tpcfe::NEIGHBOR_SELF", __FILE__, __LINE__ );
    }

    inline bool hasPositiveOrientation ( ) const {
      return ( _directionI > 7 );
    }

    inline short getOffsetI ( ) const {
      return ( _offsetI );
    }

    inline short getDirectionI ( ) const {
      return ( _directionI );
    }

    inline qc::CoordType getOffset ( ) const {
      return ( CFETopoLookup::_tetNbOffsets[ _offsetI ] );
    }

    inline qc::CoordType getDirection ( ) const {
      return ( CFETopoLookup::_tetNbOffsets[ _directionI ] );
    }

    //! radial edges must be given in the form _offsetI = self, _directionI = direction of edge; all other edges are circumferential edges.
    inline bool isRadialEdge ( ) const {
      return ( _offsetI == NEIGHBOR_SELF );
    }

    //! node of edge with smaller global index
    inline qc::CoordType getInitialPointInPositiveOrientation ( const qc::CoordType &pos ) const {
      qc::CoordType ret ( pos );
      ret += getOffset();
      if ( ! hasPositiveOrientation() )
        ret += getDirection();

      return ( ret );
    }

    //! node of edge with larger global index
    inline qc::CoordType getTerminalPointInPositiveOrientation ( const qc::CoordType &pos ) const {
      qc::CoordType ret ( pos );
      ret += getOffset();
      if ( hasPositiveOrientation() )
        ret += getDirection();

      return ( ret );
    }

    inline void getNAKedges ( const qc::CoordType &pos, qc::CoordType& E0, qc::CoordType& E1, qc::CoordType& E2 ) const {
      QUOC_ASSERT ( _directionI != NEIGHBOR_SELF );
      E0 = pos + getOffset() + getOffset();
      E1 =              E0                 + getDirection();
      E2 =                       E1                         + getDirection();
    }

    inline qc::CoordType getCOAneighbor ( const qc::CoordType &pos ) const {
      return ( pos + static_cast<short> ( 2 ) * getOffset() + getDirection() );
    }

    //! some less-than comparison for being able to use this class as key in an std::map
    inline bool operator< ( const RelevantEdgeIdentifier &other ) const {
      return ( ( this->_offsetI < other._offsetI ) || ( ( this->_offsetI == other._offsetI ) && ( this->_directionI < other._directionI ) ) );
    }

    // end nested class RelevantEdgeIdentifier
  };

  typedef CFEGrid< RealType, CFE_TPOS >                                  GridType;
  typedef typename SlopeInterfaceBase<GridType>::CoarsenIndexType        CoarsenIndexType;

  qc::AArray< std::vector< RelevantEdgeIdentifier >, qc::QC_3D >         _relevantEdges;         //!< for each node, store vector of relevant edges

public:
  static bool _useIncorrectStandardCoarsening;
  static bool _performNAKRelaxationAtCircumferentialEdges, _performNAKRelaxationAtAllRemainingEdges;
  static const char* _checkNAKQualityAfterCorrection; // file name mask, can be set to something like "out/nakQ_%02d.out"

public:
  //! Standard constructor, probably not useful. Make sure you properly set _iMap before using it.
  SlopeInterface ( ) : SlopeInterfaceBase<GridType> ( ) {}

  // copy constructor and assignment operator should do the right thing


  //! Constructor that only sets the index mapper; slopes must be set via coarsening
  explicit SlopeInterface ( const GridType& Grid ) : SlopeInterfaceBase<GridType> ( Grid ), _relevantEdges ( Grid ) {}

  void determineSlopes ( const GridType& Grid );

  void determineCoarseningWeights ( const GridType &Grid );

  void coarsenSlopeInterface ( SlopeInterfaceBase<GridType> &coarseSlIntf, const GridType &fineGrid, const GridType &coarseGrid );

  //! remove data stored
  void removeData ( ) {
    _relevantEdges.resetEntries();
    this->removeStoredData();
  }

  virtual void storeAllStandardSlopesOnly ( const GridType &Grid ) {
    this->removeStoredData();
    for ( qc::RectangularIterator<qc::QC_3D> bit ( Grid ); bit.notAtEnd(); ++bit ) {
      const qc::CoordType pos = *bit;
      for ( short nbi = 0; nbi < NUM_NEIGHBORS; ++nbi ) {
        if ( nbi != NEIGHBOR_SELF ) { // node to node itself is not an edge
          const qc::CoordType neigh = pos + CFETopoLookup::_tetNbOffsets[ nbi ];
          if ( Grid.isAdmissibleNode ( neigh )  ) { // edge present (node is present and neighbor needs to be)
            storeRelevantEdge ( pos, NEIGHBOR_SELF, nbi );
            this->setSlopes ( pos, pos, neigh, - aol::NumberTrait<RealType>::one, - aol::NumberTrait<RealType>::one );
          }
        }
      }
    }
  }

  void dumpRelevantEdges ( ) const {
    for ( qc::RectangularIterator<qc::QC_3D> ait ( _relevantEdges ); ait.notAtEnd(); ++ait ) {
      cerr << *ait << ": " << endl;
      const std::vector< RelevantEdgeIdentifier > &vc = _relevantEdges.getRef ( *ait );
      for ( typename std::vector< RelevantEdgeIdentifier >::const_iterator it = vc.begin(); it != vc.end(); ++it ) {
        cerr << CFETopoLookup::_tetNbOffsets[ it->getOffsetI() ] << " " << CFETopoLookup::_tetNbOffsets[ it->getDirectionI() ] << ";   ";
      }
      cerr << endl;
    }
  }

protected:
  inline void storeRelevantEdge ( const qc::CoordType &pos, const short offset, const short direction ) {
    _relevantEdges.getRef ( pos ).push_back ( RelevantEdgeIdentifier ( offset, direction ) );
  }

  //! pos1 must be on same side of the interface as pos (without grid, this cannot be verified here)
  inline void storeRelevantEdge ( const qc::CoordType &pos, const qc::CoordType &pos1, const qc::CoordType &pos2 ) {
    qc::CoordType offsetCo, directionCo;
    offsetCo = pos1 - pos; // offset does not cross interface, direction does
    directionCo = pos2 - pos1;

    const std::map< qc::CoordType, short >::const_iterator offsetIt = this->_neighborOffsetToIndex.find ( offsetCo ), directionIt = this->_neighborOffsetToIndex.find ( directionCo );
#ifdef DEBUG
    if ( offsetIt == this->_neighborOffsetToIndex.end() || directionIt == this->_neighborOffsetToIndex.end() ) {
      cerr << pos << "; " << pos1 << "; " << pos2 << endl;
      throw aol::Exception ( "SlopeInterface::storeRelevantEdge: illegal offset or direction", __FILE__, __LINE__ );
    }
#endif

    storeRelevantEdge ( pos, offsetIt->second, directionIt->second );
  }
  // end class SlopeInterface<RealType, CFE_TPOS>
};


// === Restriction and Prolongation Weight Management ===


/** Class used for the template specialization for cases where we want to use standard weights.
 *  \author Schwen
 */
template< typename GridType >
class MGStandardWeightProvider {
public:
  typedef SlopeInterface< GridType > SlopeInterfaceType;
  typedef typename GridType::RealType RealType;

protected:
  std::map< qc::CoordType, RealType > dummy;

public:
  // standard constructor, copy constructor and assignment operator do the correct thing

  //! get restriction weight for the nontrivial case (not fine node = coarse node)
  inline RealType getProperRestrictionWeight ( const qc::CoordType &finePos, const qc::CoordType &coarsePos ) const {
#ifdef DEBUG
    if ( finePos == static_cast<short>(2) * coarsePos )
      throw aol::Exception ( "MGStandardWeightProvider< RealType, CType >::getProperRestrictionWeight must not be called for finePos == 2 * coarsePos", __FILE__, __LINE__ );
#else
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( finePos );    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( coarsePos );
#endif
    return ( aol::NumberTrait<RealType>::one / 2 );
  }

  //! get prolongation weight for the nontrivial case (not fine node = coarse node)
  inline RealType getProperProlongationWeight ( const qc::CoordType &finePos, const qc::CoordType &coarsePos ) const {
#ifdef DEBUG
    if ( finePos == static_cast<short>(2) * coarsePos )
      throw aol::Exception ( "MGStandardWeightProvider< RealType, CType >::getProperRestrictionWeight must not be called for finePos == 2 * coarsePos", __FILE__, __LINE__ );
#else
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( finePos );    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( coarsePos );
#endif
    return ( aol::NumberTrait<RealType>::one / 2 );
  }

  const std::map< qc::CoordType, RealType >& getCoarseningWeightMapRef ( const qc::CoordType & ) const {
    throw aol::Exception ( "MGStandardWeightProvider<RealType, CType>: getCoarseningWeightMapRef must not be called.", __FILE__, __LINE__ );
    return ( dummy );
  }

  explicit MGStandardWeightProvider ( const SlopeInterface< GridType > & ) {
    throw aol::Exception ( "MGStandardWeightProvider<RealType, CType>: Constructor ( &SlopeInterface ) not useful.", __FILE__, __LINE__ );
    // but necessary for compatibility
  }

  MGStandardWeightProvider ( ) { }

  // methods not to be used:
  MGStandardWeightProvider< GridType >& operator= ( const MGStandardWeightProvider< GridType >& ) {
    throw aol::Exception ( "MGStandardWeightProvider< GridType >::operator=  must not be called.", __FILE__, __LINE__ );
    return ( *this );
  }

  void getCWeightMapIncludingImpropers ( const qc::CoordType &, std::map< qc::CoordType, RealType > & ) const {
    throw aol::Exception ( "MGStandardWeightProvider< GridType >::getCWeightMapIncludingImpropers should not be called.", __FILE__, __LINE__ ); // but may be implemented when needed
  }

};


/** Class used for the template specialization for jumping coefficient cases.
 *  \author Schwen
 */
template< typename GridType >
class MGWeightProviderJumpCoeff {

protected:
  typedef SlopeInterface< GridType > SlopeInterfaceType;
  typedef typename GridType::RealType RealType;
  const SlopeInterfaceType &_slopeInterface;

public:
  explicit MGWeightProviderJumpCoeff ( const SlopeInterfaceType &SlopeInterface ) : _slopeInterface ( SlopeInterface ) {}

  // copy constructor does the correct thing

  //! get restriction weight for the nontrivial case (not fine node = coarse node)
  inline RealType getProperRestrictionWeight ( const qc::CoordType &finePos, const qc::CoordType &coarsePos ) const {
#ifdef DEBUG
    if ( finePos == static_cast<short> ( 2 ) * coarsePos )
      throw aol::Exception ( "MGWeightProviderJumpCoeff< GridType >::getProperRestrictionWeight must not be called for finePos == 2 * coarsePos", __FILE__, __LINE__ );
#endif
    return ( _slopeInterface.getCoarseningWeight ( finePos, coarsePos ) );
  }

  //! get prolongation weight for the nontrivial case (not fine node = coarse node)
  inline RealType getProperProlongationWeight ( const qc::CoordType &finePos, const qc::CoordType &coarsePos ) const {
#ifdef DEBUG
    if ( finePos == static_cast<short> ( 2 ) * coarsePos )
      throw aol::Exception ( "MGWeightProviderJumpCoeff< GridType >::getProperProlongationWeight must not be called for finePos == 2 * coarsePos", __FILE__, __LINE__ );
#endif
    return ( _slopeInterface.getCoarseningWeight ( finePos, coarsePos ) );
  }

  const std::map< qc::CoordType, RealType >& getCoarseningWeightMapRef ( const qc::CoordType &cPPos ) const {
    return ( _slopeInterface.getCoarseningWeightMapRef ( cPPos ) );
  }

  //! \warning slow, copies map
  void getCWeightMapIncludingImpropers ( const qc::CoordType &cPos, std::map< qc::CoordType, RealType > &CWMap ) const {
    CWMap = _slopeInterface.getCoarseningWeightMapRef ( cPos );
    CWMap[ cPos + cPos ] = 1.0;
  }

// methods not to be used: standard constructor, assignment operator
  MGWeightProviderJumpCoeff ( ) : _slopeInterface ( * reinterpret_cast< const SlopeInterfaceType* > ( NULL ) ) {
    throw aol::Exception ( "MGWeightProviderJumpCoeff< GridType > standard constructor must not be called.", __FILE__, __LINE__ );
  }

  MGWeightProviderJumpCoeff< GridType >& operator= ( const MGWeightProviderJumpCoeff< GridType >& ) {
    throw aol::Exception ( "MGWeightProviderJumpCoeff< RealType, CT > assignment operator  must not be called.", __FILE__, __LINE__ );
  }
};


/** Class that provides proper restriction and prolongation weights
 *  (i. e. from neighboring fine grid nodes to coarse grid nodes, not
 *  from fine grid node coinciding with coarse grid node) for CFE
 *  multigrid.
 *  \author Schwen
 */
template< typename GridType >
class MGWeightProvider {};

/** Template specialization for complicated domain case (CFE_CD).
 */
template< typename RealType, typename NodalCoeffType >
class MGWeightProvider < CFEGrid< RealType, CFE_CD, NodalCoeffType > > : public MGStandardWeightProvider < CFEGrid< RealType, CFE_CD, NodalCoeffType > > {
  typedef CFEGrid< RealType, CFE_CD, NodalCoeffType > GridType;
public:
  MGWeightProvider ( ) {}
  explicit MGWeightProvider ( const typename MGStandardWeightProvider< GridType >::SlopeInterfaceType &SlopeInterface ) : MGStandardWeightProvider< GridType > ( SlopeInterface ) {}
};

/** Template specialization for jumping coefficient elasticty case (CFE_TPOSELAST) with incorrect coarsening weights!
 */
template< typename RealType >
class MGWeightProvider < CFEGrid < RealType, CFE_TPOSELAST, IsotropicElasticityCoefficient < RealType > > > : public MGStandardWeightProvider < CFEGrid < RealType, CFE_TPOSELAST, IsotropicElasticityCoefficient < RealType > > > {
  typedef CFEGrid< RealType, CFE_TPOSELAST, IsotropicElasticityCoefficient<RealType> > GridType;
public:
  MGWeightProvider ( ) {}
  explicit MGWeightProvider ( const typename MGStandardWeightProvider< GridType >::SlopeInterfaceType &SlopeInterface ) : MGStandardWeightProvider< GridType > ( SlopeInterface ) {}
};

template< typename RealType >
class MGWeightProvider < CFEGrid < RealType, CFE_TPOSELAST, VoigtElasticityCoefficient < RealType > > > : public MGStandardWeightProvider < CFEGrid < RealType, CFE_TPOSELAST, VoigtElasticityCoefficient < RealType > > > {
  typedef CFEGrid< RealType, CFE_TPOSELAST, VoigtElasticityCoefficient<RealType> > GridType;
public:
  MGWeightProvider ( ) {}
  explicit MGWeightProvider ( const typename MGStandardWeightProvider< GridType >::SlopeInterfaceType &SlopeInterface ) : MGStandardWeightProvider< GridType > ( SlopeInterface ) {}
};


/** Template specialization for jumping coefficient case (CFE_LIEHR).
 */
template< typename RealType >
class MGWeightProvider< CFEGrid< RealType, CFE_LIEHR > > : public MGWeightProviderJumpCoeff< CFEGrid< RealType, CFE_LIEHR > > {
public:
  MGWeightProvider ( ) {}
  explicit MGWeightProvider ( const typename MGWeightProviderJumpCoeff< CFEGrid< RealType, CFE_LIEHR > >::SlopeInterfaceType &SlopeInterface ) : MGWeightProviderJumpCoeff< CFEGrid< RealType, CFE_LIEHR > > ( SlopeInterface ) {}
};

/** Template specialization for jumping coefficient (CFE_TPOS) case
 */
template< typename RealType >
class MGWeightProvider< CFEGrid< RealType, CFE_TPOS > > : public MGWeightProviderJumpCoeff< CFEGrid< RealType, CFE_TPOS > > {
public:
  MGWeightProvider ( ) {}
  explicit MGWeightProvider ( const typename MGWeightProviderJumpCoeff< CFEGrid< RealType, CFE_TPOS > >::SlopeInterfaceType &SlopeInterface ) : MGWeightProviderJumpCoeff< CFEGrid< RealType, CFE_TPOS > > ( SlopeInterface ) {}
};



// === Prolongation Operator ===

/** Prolongation operator for tpCFE problems (for multigrid)
 *  Implemented for CFE_CD and CFE_LIEHR.
 *  For CFE_CD, (complicated domains) ignoring non-DOF nodes and respecting zero Dirichlet boundary conditions (the corresponding Dirichlet mask is stored in the tpcfe-grid).
 *  For CFE_LIEHR, obtaining the prolongation weights from a WeightProvider that is initialized with SlopeInterface.
 *  The operator can be assembled, even though this probably won't make sense.
 *  \author Schwen
 */
template< typename GridType >
class CFEProlongOp : public aol::BiOp< aol::Vector< typename GridType::RealType > > {

protected:
  typedef typename GridType::RealType   RealType;
  typedef MGWeightProvider< GridType >  MGWeightProviderType;

  const GridType                       &_coarseGrid, &_fineGrid;
  const aol::OperatorType              _opType;
  mutable aol::SparseMatrix<RealType>  *_pmat;
  const MGWeightProviderType           _mgWeightProvider;

public:
  //! Constructor in the complicated domain case
  CFEProlongOp ( const GridType &coarseGrid, const GridType &fineGrid, const aol::OperatorType OpType = aol::ONTHEFLY ) :
      _coarseGrid ( coarseGrid ),       _fineGrid ( fineGrid ),
      _opType ( OpType ),
      _pmat ( NULL ),
      _mgWeightProvider ( ) {

    if ( ( _coarseGrid.CT != CFE_CD ) && ( _coarseGrid.CT != CFE_TPOSELAST ) )
      throw aol::Exception ( "tpcfe::CFEProlongOp: This constructor may only be called for CFE_CD or CFE_TPOSELAST grids", __FILE__, __LINE__ );

    init();

  }

  //! Constructor in the jumping coefficient case
  CFEProlongOp ( const GridType &coarseGrid, const GridType &fineGrid, const SlopeInterface< GridType > &slopeInterface, const aol::OperatorType OpType = aol::ASSEMBLED ) :
      _coarseGrid ( coarseGrid ),       _fineGrid ( fineGrid ),
      _opType ( OpType ),
      _pmat ( NULL ),
      _mgWeightProvider ( slopeInterface ) {

    if ( ( _coarseGrid.CT != CFE_LIEHR ) && ( _coarseGrid.CT != CFE_TPOS ) ) {
      throw aol::Exception ( "tpcfe::CFEProlongOp: This constructor may only be called for CFE_LIEHR and CFE_TPOS grids", __FILE__, __LINE__ );
    }

    init();

  }


protected:
  void init ( ) {
    if ( _coarseGrid.getGridDepth() != _fineGrid.getGridDepth() - 1 ) {
      throw aol::Exception ( "Cannot construct CFE Prolongation operator: incorrect grid dephts", __FILE__, __LINE__ );
    }
    if ( _opType == aol::ASSEMBLED ) {
      assembleMatrix();
    }
  }

public:
  ~CFEProlongOp ( ) {
    if ( _pmat ) {
      delete ( _pmat );
    }
  }

public:
  virtual void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const;

  using aol::BiOp< aol::Vector<RealType> >::apply;

  virtual void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    cerr << "Are you sure you want to applyAdd the CFEProlongOp?" << endl; // in this case, delete this output ...
    aol::Vector<RealType> tmp ( Dest, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  using aol::BiOp< aol::Vector<RealType> >::applyAdd;

protected:
  virtual void assembleMatrix ( ) const {
    if ( _pmat ) delete _pmat ;
    _pmat = new aol::SparseMatrix<RealType> ( _fineGrid.getNumberOfNodes() , _coarseGrid.getNumberOfNodes() );
    assembleAddMatrix ( *_pmat );
  }

public:
  virtual void assembleAddMatrix ( aol::SparseMatrix<RealType> &mat ) const;

private:
  CFEProlongOp ( ); // do not implement
  CFEProlongOp ( const CFEProlongOp<GridType>& ); // do not implement
  CFEProlongOp<GridType>& operator= ( const CFEProlongOp<GridType>& ); // do not implement
  // end class CFEProlongOp
};


// === Restriction Operator ===

/** Restriction operator for tpCFE problems (for multigrid)
 *  Implemented for CFE_CD (tpcfe for complicated domains) ignoring non-DOF nodes
 *  and respecting zero Dirichlet boundary conditions (the corresponding Dirichlet mask is stored in the tpcfe-grid.
 *  The operator can be assembled, even though this probably won't make sense.
 *  \author Schwen
 */
template< typename GridType >
class CFERestrictOp : public aol::BiOp< aol::Vector< typename GridType::RealType > > {

protected:
  typedef typename GridType::RealType   RealType;
  typedef MGWeightProvider< GridType >  MGWeightProviderType;

  const GridType                       &_coarseGrid;
  const GridType                       &_fineGrid;
  const aol::OperatorType              _opType;
  mutable aol::SparseMatrix<RealType>  *_pmat;
  const MGWeightProviderType           _mgWeightProvider;

public:
  //! Constructor in the complicated domain case: assembling not necessary
  CFERestrictOp ( const GridType &coarseGrid, const GridType &fineGrid, const aol::OperatorType opType = aol::ONTHEFLY ) :
      _coarseGrid ( coarseGrid ), _fineGrid ( fineGrid ),
      _opType ( opType ),
      _pmat ( NULL ),
      _mgWeightProvider ( ) {

    if ( ( _coarseGrid.CT != CFE_CD ) && ( _coarseGrid.CT != CFE_TPOSELAST ) )
      throw aol::Exception ( "tpcfe::CFERestrictOp: This constructor may only be called for CFE_CD grids", __FILE__, __LINE__ );

    init ();
  }

  //! Constructor in the jumping coefficient case: assembling may be better because nonstandard weights are stored in map and need to be looked up.
  CFERestrictOp ( const GridType &coarseGrid, const GridType &fineGrid, const SlopeInterface< GridType > &slopeInterface, const aol::OperatorType OpType = aol::ASSEMBLED ) :
      _coarseGrid ( coarseGrid ), _fineGrid ( fineGrid ),
      _opType ( OpType ),
      _pmat ( NULL ),
      _mgWeightProvider ( slopeInterface ) {

    if ( ( _coarseGrid.CT != CFE_LIEHR ) && ( _coarseGrid.CT != CFE_TPOS ) ) {
      throw aol::Exception ( "tpcfe::CFERestrictOp: This constructor may only be called for CFE_LIEHR grids", __FILE__, __LINE__ );
    }

    init ();
  }


protected:
  void init ( ) {
    if ( _coarseGrid.getGridDepth() != _fineGrid.getGridDepth() - 1 ) {
      throw aol::Exception ( "Cannot construct CFE Restriction operator: incorrect grid dephts", __FILE__, __LINE__ );
    }
    if ( _opType == aol::ASSEMBLED ) {
      assembleMatrix();
    }
  }

public:
  ~CFERestrictOp ( ) {
    if ( _pmat ) {
      delete ( _pmat );
    }
  }

public:
  virtual void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const;
  using aol::BiOp< aol::Vector< RealType > >::apply;


  virtual void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    cerr << "Are you sure you want to applyAdd the CFERestrictOp?" << endl; // in this case, delete this output ...
    aol::Vector<RealType> tmp ( Dest, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  using aol::BiOp< aol::Vector< RealType > >::applyAdd;

protected:
  virtual void assembleMatrix ( ) const { // obviously, constness is nonsense. but necessary ...
    if ( _pmat ) delete _pmat ;
    _pmat = new aol::SparseMatrix<RealType> ( _coarseGrid.getNumberOfNodes(), _fineGrid.getNumberOfNodes() );
    assembleAddMatrix ( *_pmat );
  }

public:
  virtual void assembleAddMatrix ( aol::SparseMatrix<RealType> &mat ) const;

private:
  CFERestrictOp ( ); // do not implement
  CFERestrictOp ( const CFERestrictOp<GridType>& ); // do not implement
  CFERestrictOp<GridType>& operator= ( const CFERestrictOp<GridType>& ); // do not implement
  // end class CFERestrictOp
};


// === Operator (Matrix) Coarsening ===

/** Operator coarsening for CFE matrices with 15-band structure (CFE_CD, CFE_LIEHR)
 *  \param dropSmallEntries drop entries in the coarsened matrices that are smaller than this value (addition and subtraction may lead to such ghost entries)
 */
template< typename ConfiguratorType >
void coarsenCFEMatrix ( const typename ConfiguratorType::MatrixType &fineMat, typename ConfiguratorType::MatrixType &coarseMat,
                        const typename ConfiguratorType::GridType &fineGrid, const typename ConfiguratorType::GridType &coarseGrid,
                        const typename ConfiguratorType::RealType diagEntry,
                        const typename ConfiguratorType::RealType dropSmallEntries = aol::NumberTrait< typename ConfiguratorType::RealType >::zero,
                        const MGWeightProvider< typename ConfiguratorType::GridType > mgWeightProvider = MGWeightProvider< typename ConfiguratorType::GridType > ( ) ) {

  if ( ConfiguratorType::CT == tpcfe::CFE_TPOS ) {
    cerr << aol::color::brown << "on-the-fly coarsening of CFE_TPOS matrices is slow" << aol::color::reset;
  }

  typedef typename ConfiguratorType::RealType RealType;

#ifdef VERBOSE
  aol::ProgressBar<> pb ( "Coarsening matrix" );
  pb.start ( coarseGrid.getNumberOfNodes() );
#endif

  for ( qc::RectangularIterator<qc::QC_3D> bit ( coarseGrid ); bit.notAtEnd(); ++bit ) {
#ifdef VERBOSE
    pb++;
#endif
    const int coarse_i = coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit ); // global node index in coarse grid
    const qc::CoordType iPos = *bit;

    if ( coarseGrid.isDOFNode ( coarse_i ) ) { // automatically inside (loops)

      switch ( ConfiguratorType::CT ) {
        case CFE_CD:
        case CFE_LIEHR: {

          for ( int _l = 0; _l < NUM_NEIGHBORS; ++_l ) {
            const qc::CoordType lPos = iPos + CFELookup<RealType>::_tetNbOffsets[_l];
            int coarse_l;

            if ( coarseGrid.isAdmissibleNode ( lPos ) &&  coarseGrid.isDOFNode ( coarse_l = coarseGrid.getIndexMapperRef().getGlobalIndex ( lPos ) ) ) { // only set l to this value if inside.

              for ( int _j = 0; _j < NUM_NEIGHBORS; ++_j ) {
                const qc::CoordType jPos = static_cast<short> ( 2 ) * iPos + CFELookup<RealType>::_tetNbOffsets[_j];
                int fine_j;

                if ( fineGrid.isAdmissibleNode ( jPos ) && fineGrid.isDOFNode ( fine_j = fineGrid.getIndexMapperRef().getGlobalIndex ( jPos ) ) ) {

                  for ( int _k = 0; _k < NUM_NEIGHBORS; ++_k ) {
                    const qc::CoordType kPos = static_cast<short> ( 2 ) * lPos + CFELookup<RealType>::_tetNbOffsets[_k];
                    int fine_k;

                    if ( fineGrid.isAdmissibleNode ( kPos ) && fineGrid.isDOFNode ( fine_k = fineGrid.getIndexMapperRef().getGlobalIndex ( kPos ) ) ) {

#ifdef VERBOSE
                      if ( mgWeightProvider.getProperRestrictionWeight ( jPos, iPos ) != 0.5 || mgWeightProvider.getProperProlongationWeight ( kPos, lPos ) != 0.5 )
                        cerr << fine_j << " " << coarse_i << " " << mgWeightProvider.getProperRestrictionWeight ( jPos, iPos ) << "     "
                        << fine_k << " " << coarse_l << " " << mgWeightProvider.getProperProlongationWeight ( kPos, lPos ) << endl;
#endif

                      const RealType re_ij = ( _j == NEIGHBOR_SELF ? 1.0 : mgWeightProvider.getProperRestrictionWeight ( jPos, iPos ) ); // 1, if coarse and fine point coincide, 1/2 if interpolation point nearby
                      const RealType fi_jk = fineMat.get ( fine_j, fine_k );
                      const RealType pr_kl = ( _k == NEIGHBOR_SELF ? 1.0 : mgWeightProvider.getProperProlongationWeight ( kPos, lPos ) );

                      coarseMat.add ( coarse_i, coarse_l, re_ij * fi_jk * pr_kl );

                    } // k dof
                  } // _k
                } // j dof
              } // _j

              if ( fabs ( coarseMat.get ( coarse_i, coarse_l ) ) < dropSmallEntries ) {
#ifdef VERBOSE
                cerr << "Dropping entry: " << coarse_i << " " << coarse_l << " " << coarseMat.get ( coarse_i, coarse_l ) << endl;
#endif
                coarseMat.set ( coarse_i, coarse_l, aol::NumberTrait<RealType>::zero );
              }

            } // l dof
          } // _l

        }
        break;

        case CFE_TPOS: {
          for ( qc::LocalLInfBoxIterator< qc::QC_3D > lit ( iPos, 3, aol::Vec3<int> ( 0, 0, 0 ), coarseGrid.getSize() ); lit.notAtEnd(); ++lit ) {
            const qc::CoordType lPos = *lit;
            const int coarse_l = coarseGrid.getIndexMapperRef().getGlobalIndex ( lPos );
            if ( coarseGrid.isDOFNode ( coarse_l ) ) {

              typedef std::map< qc::CoordType, RealType > CWMapType;
              typedef typename CWMapType::const_iterator CWMapItType;

              CWMapType jmap;
              mgWeightProvider.getCWeightMapIncludingImpropers ( iPos, jmap );
              for ( CWMapItType jit = jmap.begin(); jit != jmap.end(); ++jit ) {
                const qc::CoordType jPos = jit->first;
                const int fine_j = fineGrid.getIndexMapperRef().getGlobalIndex ( jPos );
                if ( fineGrid.isDOFNode ( fine_j ) ) {

                  CWMapType kmap;
                  mgWeightProvider.getCWeightMapIncludingImpropers ( lPos, kmap );
                  for ( CWMapItType kit = kmap.begin(); kit != kmap.end(); ++kit ) {
                    const qc::CoordType kPos = kit->first;
                    const int fine_k = fineGrid.getIndexMapperRef().getGlobalIndex ( kPos );
                    if ( fineGrid.isDOFNode ( fine_k ) ) {

                      if ( !fineGrid.isDOFNode ( fine_k ) )
                        cout << "Trying nonDOF k " << kPos << endl;

                      const RealType re_ij = jit->second;
                      const RealType fi_jk = fineMat.get ( fine_j, fine_k );
                      const RealType pr_kl = kit->second;

                      coarseMat.add ( coarse_i, coarse_l, re_ij * fi_jk * pr_kl );

                    } // k DOF
                  } // kit loop
                } // j DOF
              } // jit loop

              if ( fabs ( coarseMat.get ( coarse_i, coarse_l ) ) < dropSmallEntries ) {
#ifdef VERBOSE
                cerr << "Dropping entry: " << coarse_i << " " << coarse_l << " " << coarseMat.get ( coarse_i, coarse_l ) << endl;
#endif
                coarseMat.set ( coarse_i, coarse_l, aol::NumberTrait<RealType>::zero );
              }
            }
          } // lit loop

        }
        break;

        default:
          throw aol::UnimplementedCodeException ( "tpcfe::coarsenCFEMatrix: only implemented for CFE_CD, CFE_LIEHR_, CFE_TPOS", __FILE__, __LINE__ );
      }

    } else  { // coarse non-DOF node
      coarseMat.set ( coarse_i, coarse_i, diagEntry );  // zero for off-diagonal blocks (vector valued case)

    } // i dof
  } // coarse grid nodes

#ifdef VERBOSE
  pb.finish();
#endif

}


// === The Multigrid Solver Classes (scalar and vector-valued case ) ===

/** Multigrid solver for complicated domain CFE_CD and jumping coefficients CFE_LIEHR, CFE_TPOS
 *  \author Schwen
 */
template< typename CFEOpType, aol::OperatorType RePrOpType, mg::CoarseningMode CMode >
class CFEMultigrid : public mg::ExplicitOperatorHierarchyMultigrid < typename CFEOpType::VectorType,
                                                                     typename CFEOpType::ConfiguratorType::MatrixType,
                                                                     mg::GaussSeidelSmoother< typename CFEOpType::VectorType, typename CFEOpType::ConfiguratorType::MatrixType >,
                                                                     CFERestrictOp< typename CFEOpType::ConfiguratorType::GridType >,
                                                                     CFEProlongOp< typename CFEOpType::ConfiguratorType::GridType >,
                                                                     typename CFEOpType::ConfiguratorType::GridType > {

protected:
  typedef typename CFEOpType::ConfiguratorType::RealType      RealType;
  typedef typename CFEOpType::VectorType                      VectorType;
  typedef typename CFEOpType::ConfiguratorType::MatrixType    MatrixType;
  typedef typename CFEOpType::ConfiguratorType::GridType      GridType;

  aol::GaussSeidelSweepingMode              _gsmode;
  RealType                                  _relax;
  RealType                                  _coarsenDropSmallEntries;

  std::vector< SlopeInterface<GridType> >   _slopeInterfaces; // note: if restriction/prolongation operators are assembled (which should be the standard case), this wastes some memory

public:

  //! Constructor for both complicated domain and jumping coefficient cases.
  CFEMultigrid ( const GridType &Grid,
                 const MatrixType &fineMatrix,
                 const int ExplicitLevel = 2,
                 const int PreSmooth = 3, const int PostSmooth = 3,
                 const aol::GaussSeidelSweepingMode GSmode = aol::GAUSS_SEIDEL_FORWARD,
                 const RealType Relax = 1.0,
                 const mg::MGCycleMode CycleMode = mg::MULTIGRID_V_CYCLE,
                 const RealType CycleEps = 1.0e-16,
                 const int MaxCycle = 100,
                 const aol::StoppingMode stop = aol::STOPPING_UNSET,
                 const RealType coarsenDropSmallEntries = 0.0f )
      : mg::ExplicitOperatorHierarchyMultigrid < VectorType,
                                                 MatrixType,
                                                 mg::GaussSeidelSmoother< VectorType, MatrixType >,
                                                 CFERestrictOp< GridType >,
                                                 CFEProlongOp< GridType >, GridType > ( fineMatrix,
                                                                                        Grid,
                                                                                        ExplicitLevel,
                                                                                        PreSmooth, PostSmooth,
                                                                                        Relax,
                                                                                        CycleMode, CycleEps, MaxCycle, stop ),
        _gsmode ( GSmode ),
        _relax ( Relax ),
        _coarsenDropSmallEntries ( coarsenDropSmallEntries ),
        _slopeInterfaces ( Grid.getGridDepth() + 1 ) {

    switch ( Grid.CT ) {
      case CFE_CD:
        // do nothing
        break;
      case CFE_LIEHR:
      case CFE_TPOS:
        _slopeInterfaces[ Grid.getGridDepth() ] = SlopeInterface<GridType> ( Grid ); // create and copy temporary object slightly inefficient, but currently not avoidable
        _slopeInterfaces[ Grid.getGridDepth() ].determineSlopes ( Grid );

        break;

      default:
        throw aol::Exception ( "tpcfe::CFEMultigrid: illegal constraint type", __FILE__, __LINE__ );
    }

    this->init();

    this->setCoarseSolverSteps ( 5000 );
    this->setCoarseSolverEps ( 1.0e-16 );

  }


public:
  virtual void coarsenGrid ( const int coarseLevel ) {
    GridType *pCoarseGrid = new GridType ( coarseLevel );

    switch ( GridType::CT ) {
      case CFE_CD:
        coarsenDomainNodeMask ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid ); // only for CFE_CD
        coarsenDirichletMask  ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid ); // also in other cases
        coarsenDOFMask        ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid );
        break;

      case CFE_LIEHR:
      case CFE_TPOS:
        coarsenDirichletMask  ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid ); // also in other cases
        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEMultigrid::coarsenGrid: illegal constraint type.", __FILE__, __LINE__ );

    } // end switch

    pCoarseGrid->setCoarsenedFlag( true );

    this->setpGrid ( coarseLevel, pCoarseGrid );
  }

protected:
  virtual void init ( ) {
    cerr << "Setting up list of grids and slopeInterfaces ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      coarsenGrid ( i );

      switch ( GridType::CT ) {
        case CFE_CD:
          break; // do nothing
        case CFE_LIEHR:
        case CFE_TPOS:
          _slopeInterfaces[ i ] = SlopeInterface<GridType> ( this->getGridRef ( i ) ); // create but do not fill
          break;
        default:
          throw aol::UnimplementedCodeException ( "CFEMultigrid::init: illegal constraint type.", __FILE__, __LINE__ );
      }

      cerr << i << " ";
    }

    switch ( GridType::CT ) {
      case CFE_CD:
        break; // do nothing
      case CFE_LIEHR:
      case CFE_TPOS:
        cerr << " coarsening slopeInterfaces ... ";
        for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
          coarsenSlopeInterface ( i );                                               // fill SlopeInterfaces here
          cerr << i << " ";
        }
        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEMultigrid::init: illegal constraint type.", __FILE__, __LINE__ );

    } // end switch

    cerr << " setting up restriction and prolongation operators ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      setProlongationOperator ( i );
      setRestrictionOperator ( i );
      cerr << i << " ";
    }

    cerr << " coarsening operators ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      coarsenOperator ( i );
      cerr << i << " ";
    }

    cerr << " setting up smoothOps ... ";
    for ( int i = ( this->_depth - 0 ); i > this->_explicitLevel; --i ) {
      // need smoothOp on finest level but not on _explicitLevel
      setSmoother ( i );
      cerr << i << " ";
    }

    setCoarsePreconditioner();

    // \todo Think about clearing the SlopeInterfaces at this point to potentially free memory.

    cerr << " done." << endl;

  }

  virtual void setRestrictionOperator ( const int level ) {
    switch ( GridType::CT ) {
      case CFE_CD:
        this->setpRestrictionOperator ( level, new CFERestrictOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), RePrOpType ) );
        break;

      case CFE_LIEHR:
      case CFE_TPOS:
        this->setpRestrictionOperator ( level, new CFERestrictOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), _slopeInterfaces[ level + 1 ], RePrOpType ) );
        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEMultigrid::setRestrictionOperator: illegal constraint type.", __FILE__, __LINE__ );

    } // end switch
  }

  virtual void setProlongationOperator ( const int level ) {
    switch ( GridType::CT ) {
      case CFE_CD:
        this->setpProlongationOperator ( level, new CFEProlongOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), RePrOpType ) );
        break;

      case CFE_LIEHR:
      case CFE_TPOS:
        this->setpProlongationOperator ( level, new CFEProlongOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), _slopeInterfaces[ level + 1 ], RePrOpType ) );
        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEMultigrid::setProlongationOperator: illegal constraint type.", __FILE__, __LINE__ );

    } // end switch
  }


  virtual void coarsenOperator ( const int coarseLevel ) {
    const int fineLevel = coarseLevel + 1;

    MatrixType *newpMat = new MatrixType ( this->getGridRef ( coarseLevel ) );

    switch ( CMode ) {
      case mg::MATRIXMULT_COARSENING: {
        // "assembled" version of coarsening
        // if we have assembled prolongation and restriction operators, assembling again may be unnecessary here

        const int coarseSize  = this->getGridRef ( coarseLevel ).getNumberOfNodes();
        const int fineSize    = this->getGridRef ( fineLevel ).getNumberOfNodes();

        typedef aol::SparseMatrix<RealType> SparseMatrixType;
        SparseMatrixType *reMat = new SparseMatrixType ( coarseSize, fineSize );
        this->getRestrictionOpRef ( coarseLevel ).assembleAddMatrix ( ( *reMat ) );

        cerr << "l";
        SparseMatrixType *reXfineop = new SparseMatrixType ( coarseSize, fineSize );
        ( *reXfineop ).addMatrixProduct ( ( *reMat ), this->getOperatorRef ( fineLevel ) ); // i. e. get matrix

        delete ( reMat );

        SparseMatrixType *prMat  = new SparseMatrixType ( fineSize, coarseSize );
        this->getProlongationOpRef ( coarseLevel ).assembleAddMatrix ( ( *prMat ) );

        cerr << "r";
        ( *newpMat ).addMatrixProduct ( ( *reXfineop ), ( *prMat ) );

        delete ( reXfineop );
        delete ( prMat );

        for ( int i = 0; i < coarseSize; ++i ) {
          if ( ! ( this->getGridRef ( coarseLevel ).isDOFNode ( i ) ) ) {
            newpMat->set ( i, i, aol::NumberTrait< RealType >::one );
          }
        }

      } break;

      case mg::ONTHEFLY_COARSENING:
        // "on-the-fly" coarsening (seems to be faster for CFE_CD, what about CFE_LIEHR?)

        switch ( GridType::CT ) {
            // \todo think about whether this distinction is necessary: first call creates "empty" weightProvider which is no better than the one we have.
            // same above!
          case CFE_CD:
            coarsenCFEMatrix<typename CFEOpType::ConfiguratorType> ( this->getOperatorRef ( coarseLevel + 1 ), ( *newpMat ),
                                                                     this->getGridRef ( coarseLevel + 1 ), this->getGridRef ( coarseLevel ),
                                                                     1, _coarsenDropSmallEntries );
            break;

          case CFE_LIEHR:
          case CFE_TPOS:
            coarsenCFEMatrix<typename CFEOpType::ConfiguratorType> ( this->getOperatorRef ( coarseLevel + 1 ), ( *newpMat ),
                                                                     this->getGridRef ( coarseLevel + 1 ), this->getGridRef ( coarseLevel ),
                                                                     1.0, _coarsenDropSmallEntries,
                                                                     MGWeightProvider< GridType > ( _slopeInterfaces[ fineLevel ] ) );
            break;

          default:
            throw aol::UnimplementedCodeException ( "CFEMultigrid::coarsenOperator: illegal constraint type.", __FILE__, __LINE__ );

        } // end switch


        break;

      default:
        throw aol::Exception ( "tpcfe::CFEMultigrid: illegal RePrCoOpType", __FILE__, __LINE__ );
    } // end switch

    this->setpOperator ( coarseLevel, newpMat );

  }

  void coarsenSlopeInterface ( const int coarseLevel ) {
    _slopeInterfaces[ coarseLevel + 1 ].coarsenSlopeInterface ( _slopeInterfaces[ coarseLevel ], this->getGridRef ( coarseLevel + 1 ), this->getGridRef ( coarseLevel ) );
  }

  virtual void setSmoother ( const int level ) {
    this->setpSmoother ( level, new mg::GaussSeidelSmoother< VectorType, MatrixType > ( this->getOperatorRef ( level ), -1, _relax, _gsmode ) ); // -1 iterations will be changed before smoothing
  }

  virtual void setCoarsePreconditioner ( ) {
    this->_coarsePreconditioner = new aol::DiagonalPreconditioner< aol::Vector<RealType> > ( this->getOperatorRef ( this->_explicitLevel ) );
  }

  // end class CFEMultigrid
};


/** Multigrid solver for CFE_CD vector valued problems that requires a block matrix on the finest grid
 *  \author Schwen
 */
template< typename CFEBlockOpType, aol::OperatorType RePrOpType, mg::CoarseningMode CMode >
class CFEBlockMultigrid : public mg::ExplicitOperatorHierarchyMultigrid < typename CFEBlockOpType::VectorType,
                                                                          aol::BlockMatrix< typename CFEBlockOpType::ConfiguratorType::MatrixType >,
                                                                          mg::BlockGaussSeidelSmoother< typename CFEBlockOpType::VectorType, aol::BlockMatrix< typename CFEBlockOpType::ConfiguratorType::MatrixType > > ,
                                                                          CFERestrictOp< typename CFEBlockOpType::ConfiguratorType::GridType >,
                                                                          CFEProlongOp< typename CFEBlockOpType::ConfiguratorType::GridType >, typename CFEBlockOpType::ConfiguratorType::GridType > {

protected:
  typedef typename CFEBlockOpType::ConfiguratorType::RealType     RealType;
  typedef typename CFEBlockOpType::VectorType                     VectorType; // note: VectorType is MultiVector in Block case!
  typedef typename CFEBlockOpType::ConfiguratorType::MatrixType   MatrixType;
  typedef typename CFEBlockOpType::ConfiguratorType::GridType     GridType;
  typedef aol::BlockMatrix< MatrixType >                          BlockMatrixType;

  aol::GaussSeidelSweepingMode   _gsmode;
  RealType                       _relax;
  RealType                       _coarsenDropSmallEntries;

  std::vector< SlopeInterface<GridType> >   _slopeInterfaces; // note: if restriction/prolongation operators are assembled (which should be the standard case), this wastes some memory

public:
  CFEBlockMultigrid ( const GridType &Grid,
                      const BlockMatrixType &fineBlockMatrix,
                      const int ExplicitLevel = 2,
                      const int PreSmooth = 3, const int PostSmooth = 3,
                      const aol::GaussSeidelSweepingMode GSmode = aol::GAUSS_SEIDEL_FORWARD,
                      const RealType Relax = 1.0,
                      const mg::MGCycleMode CycleMode = mg::MULTIGRID_V_CYCLE,
                      const RealType CycleEps = 1.0e-16,
                      const int MaxCycle = 100,
                      const aol::StoppingMode stop = aol::STOPPING_UNSET,
                      const RealType coarsenDropSmallEntries = 0 )
      : mg::ExplicitOperatorHierarchyMultigrid < VectorType,
                                                 aol::BlockMatrix< MatrixType >,
                                                 mg::BlockGaussSeidelSmoother< VectorType, aol::BlockMatrix< MatrixType > > ,
                                                 CFERestrictOp< GridType >,
                                                 CFEProlongOp< GridType >, GridType > ( fineBlockMatrix,
                                                                                        Grid,
                                                                                        ExplicitLevel,
                                                                                        PreSmooth, PostSmooth,
                                                                                        Relax,
                                                                                        CycleMode, CycleEps, MaxCycle, stop ),
        _gsmode ( GSmode ),
        _relax ( Relax ),
        _coarsenDropSmallEntries ( coarsenDropSmallEntries ),
        _slopeInterfaces ( Grid.getGridDepth() + 1 ) {

    switch ( Grid.CT ) {
      case CFE_CD:
      case CFE_TPOSELAST:
        // do nothing
        break;
      case CFE_LIEHR:
      case CFE_TPOS:
        _slopeInterfaces[ Grid.getGridDepth() ] = SlopeInterface<GridType> ( Grid ); // create and copy temporary object slightly inefficient, but currently not avoidable
        _slopeInterfaces[ Grid.getGridDepth() ].determineSlopes ( Grid );
        break;

      default:
        throw aol::Exception ( "tpcfe::CFEBlockMultigrid: illegal constraint type", __FILE__, __LINE__ );
    }

    this->init();

    this->setCoarseSolverSteps ( 5000 );
    this->setCoarseSolverEps ( 1.0e-16 );

  }

protected:

  virtual void coarsenGrid ( const int coarseLevel ) {
    GridType *pCoarseGrid = new GridType ( coarseLevel );

    switch ( GridType::CT ) {
      case CFE_CD:
        coarsenDomainNodeMask ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid ); // only for CFE_CD
        coarsenDirichletMask  ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid ); // also in other cases
        coarsenDOFMask        ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid );
        break;

      case CFE_LIEHR:
      case CFE_TPOS:
      case CFE_TPOSELAST:
        coarsenDirichletMask  ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid ); // coarsening DOF mask not necessary, here DOF <=> !Dirichlet
        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEMultigrid::coarsenGrid: illegal constraint type.", __FILE__, __LINE__ );

    } // end switch

    pCoarseGrid->setCoarsenedFlag( true );

    this->setpGrid ( coarseLevel, pCoarseGrid );
  }

  virtual void init ( ) {
    cerr << "Setting up list of grids and slopeInterfaces ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      coarsenGrid ( i );

      switch ( GridType::CT ) {
        case CFE_CD:
        case CFE_TPOSELAST:
          break; // do nothing
        case CFE_LIEHR:
        case CFE_TPOS:
          _slopeInterfaces[ i ] = SlopeInterface<GridType> ( this->getGridRef ( i ) ); // create but do not fill
          break;
        default:
          throw aol::UnimplementedCodeException ( "CFEMultigrid::init: illegal constraint type.", __FILE__, __LINE__ );
      }

      cerr << i << " ";
    }

    switch ( GridType::CT ) {
      case CFE_CD:
      case CFE_TPOSELAST:
        break; // do nothing
      case CFE_LIEHR:
      case CFE_TPOS:
        cerr << " coarsening slopeInterfaces ... ";
        for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
          coarsenSlopeInterface ( i );                                               // fill SlopeInterfaces here
          cerr << i << " ";
        }
        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEMultigrid::init: illegal constraint type.", __FILE__, __LINE__ );

    } // end switch

    cerr << " setting up restriction and prolongation operators ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      setProlongationOperator ( i );
      setRestrictionOperator ( i );
      cerr << i << " ";
    }

    cerr << " coarsening operators ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      coarsenOperator ( i );
      cerr << i << " ";
    }

    cerr << " setting up smoothOps ... ";
    for ( int i = ( this->_depth - 0 ); i > this->_explicitLevel; --i ) {
      // need smoothOp on finest level but not on _explicitLevel
      setSmoother ( i );
      cerr << i << " ";
    }

    setCoarsePreconditioner();

    // \todo Think about clearing the SlopeInterfaces at this point to potentially free memory.

    cerr << " done." << endl;

  }

  virtual void setRestrictionOperator ( const int level ) {
    switch ( GridType::CT ) {
      case CFE_CD:
      case CFE_TPOSELAST:
        this->setpRestrictionOperator ( level, new CFERestrictOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), RePrOpType ) );
        break;

      case CFE_LIEHR:
      case CFE_TPOS:
        this->setpRestrictionOperator ( level, new CFERestrictOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), _slopeInterfaces[ level + 1 ], RePrOpType ) );
        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEMultigrid::setRestrictionOperator: illegal constraint type.", __FILE__, __LINE__ );

    } // end switch
  }

  virtual void setProlongationOperator ( const int level ) {
    switch ( GridType::CT ) {
      case CFE_CD:
      case CFE_TPOSELAST:
        this->setpProlongationOperator ( level, new CFEProlongOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), RePrOpType ) );
        break;

      case CFE_LIEHR:
      case CFE_TPOS:
        this->setpProlongationOperator ( level, new CFEProlongOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), _slopeInterfaces[ level + 1 ], RePrOpType ) );
        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEMultigrid::setProlongationOperator: illegal constraint type.", __FILE__, __LINE__ );

    } // end switch
  }

  virtual void coarsenOperator ( const int coarseLevel ) {
    const int fineLevel = coarseLevel + 1;
    BlockMatrixType *newpBlockMat = new BlockMatrixType ( this->getGridRef ( coarseLevel ) );

    switch ( CMode ) {
      case mg::MATRIXMULT_COARSENING: {
        // "assembled" version of coarsening
        // if we have assembled prolongation and restriction operators, assembling again may be unnecessary here

        const int coarseSize  = this->getGridRef ( coarseLevel ).getNumberOfNodes();
        const int fineSize    = this->getGridRef ( fineLevel ).getNumberOfNodes();

        typedef aol::SparseMatrix<RealType> SparseMatrixType;
        SparseMatrixType reMat ( coarseSize, fineSize ), prMat ( fineSize, coarseSize );
        this->getRestrictionOpRef ( coarseLevel ).assembleAddMatrix ( reMat );
        this->getProlongationOpRef ( coarseLevel ).assembleAddMatrix ( prMat );

        // was ist hier mit den Diagonaleintraegen (0/1) auf der feinen Matrix, die zu nicht-DOFs korrespondieren??

        for ( int i = 0; i < 3; ++i ) {
          for ( int j = 0; j < 3; ++j ) {
            SparseMatrixType reXfineop ( coarseSize, fineSize );
            cerr << "l";
            reXfineop.addMatrixProduct ( reMat, this->getOperatorRef ( fineLevel ).getReference ( i, j ) );
            cerr << "r";
            newpBlockMat->getReference ( i, j ).addMatrixProduct ( reXfineop, prMat );

            for ( int k = 0; k < coarseSize; ++k ) {
              if ( ! ( this->getGridRef ( coarseLevel ).isDOFNode ( k ) ) ) {
                newpBlockMat->getReference ( i, j ).set ( k, k, ( i == j ? aol::NumberTrait<RealType>::one : aol::NumberTrait<RealType>::zero ) );
              }
            }
          }
        }

      } break;

      case mg::ONTHEFLY_COARSENING: {

          for ( int i = 0; i < 3; ++i ) {
            for ( int j = 0; j < 3; ++j ) {
              switch ( GridType::CT ) {
                  // \todo think about whether this distinction is necessary: first call creates "empty" weightProvider which is no better than the one we have.
                  // same above!
                case CFE_CD:
                  coarsenCFEMatrix<typename CFEBlockOpType::ConfiguratorType> ( this->getOperatorRef ( coarseLevel + 1 ).getReference ( i, j ),  newpBlockMat->getReference ( i, j ),
                                                                                this->getGridRef ( coarseLevel + 1 ), this->getGridRef ( coarseLevel ),
                                                                                ( i == j ? aol::NumberTrait<RealType>::one : aol::NumberTrait<RealType>::zero ), _coarsenDropSmallEntries );
                  break;

                case CFE_LIEHR:
                case CFE_TPOS:
                  coarsenCFEMatrix<typename CFEBlockOpType::ConfiguratorType> ( this->getOperatorRef ( coarseLevel + 1 ).getReference ( i, j ),  newpBlockMat->getReference ( i, j ),
                                                                                this->getGridRef ( coarseLevel + 1 ), this->getGridRef ( coarseLevel ),
                                                                                ( i == j ? aol::NumberTrait<RealType>::one : aol::NumberTrait<RealType>::zero ), _coarsenDropSmallEntries,
                                                                                MGWeightProvider< GridType > ( _slopeInterfaces[ fineLevel ] ) );
                  break;

                default:
                  throw aol::UnimplementedCodeException ( "CFEBlockMultigrid::coarsenOperator: illegal constraint type.", __FILE__, __LINE__ );

              } // end switch

            }
          }
        }
        break;

      default:
        throw aol::Exception ( "CFEBlockMultigrid::coarsenOperator: illegal Coarsening Mode", __FILE__, __LINE__ );
    }


    this->setpOperator ( coarseLevel, newpBlockMat );
  }

  void coarsenSlopeInterface ( const int coarseLevel ) {
    _slopeInterfaces[ coarseLevel + 1 ].coarsenSlopeInterface ( _slopeInterfaces[ coarseLevel ], this->getGridRef ( coarseLevel + 1 ), this->getGridRef ( coarseLevel ) );
  }

  virtual VectorType* createNewVector ( const int Level ) const {
    // creation of MultiVectors works differently ...
    VectorType* pNewVector = new VectorType ( this->getGridRef ( Level ) );
    return ( pNewVector );
  }

  virtual int getVectorSize ( int const Level ) const {
    // total vector size (for L2-norms) is dimension times vector size
    return ( this->getGridRef ( Level ).getDimOfWorld() *  this->getGridRef ( Level ).getNumberOfNodes() );
  }

  virtual void setSmoother ( const int level ) {
    // block-GaussSeidel requires BlockMatrix
    this->setpSmoother ( level, new mg::BlockGaussSeidelSmoother< VectorType, aol::BlockMatrix< MatrixType > > ( this->getOperatorRef ( level ), -1, _relax, _gsmode ) ); // -1 iterations will be changed before smoothing
  }

  virtual void setCoarsePreconditioner ( ) {
    this->_coarsePreconditioner = new aol::BlockDiagonalPreconditioner< RealType, qc::QC_3D > ( this->getOperatorRef ( this->_explicitLevel ) );
  }

  // end class CFEBlockMultigrid
};


template< typename CFEOpType, aol::OperatorType RePrOpType, mg::CoarseningMode CMode >
class CFEMultigridPreconditioner : public CFEMultigrid < CFEOpType, RePrOpType, CMode > {
  typedef typename CFEOpType::ConfiguratorType::GridType      GridType;
  typedef typename CFEOpType::ConfiguratorType::MatrixType    MatrixType;

public:
  CFEMultigridPreconditioner ( const GridType &Grid,
                               const MatrixType &fineMatrix,
                               const int ExplicitLevel = 2,
                               const int PreSmooth = 3, const int PostSmooth = 3,
                               const aol::GaussSeidelSweepingMode GSmode = aol::GAUSS_SEIDEL_FORWARD )
  : CFEMultigrid  < CFEOpType, RePrOpType, CMode > ( Grid, fineMatrix, ExplicitLevel, PreSmooth, PostSmooth, GSmode, 1.0, mg::MULTIGRID_V_CYCLE, 1.0 - 1.0e-16, 1, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM ) {
    this->setVerboseMode ( 0 );
  }
};

template< typename CFEBlockOpType, aol::OperatorType RePrOpType, mg::CoarseningMode CMode >
class CFEBlockMultigridPreconditioner : public CFEBlockMultigrid < CFEBlockOpType, RePrOpType, CMode > {
  typedef typename CFEBlockOpType::ConfiguratorType::GridType    GridType;
  typedef typename CFEBlockOpType::ConfiguratorType::MatrixType  MatrixType;
  typedef aol::BlockMatrix< MatrixType >                         BlockMatrixType;

public:
  CFEBlockMultigridPreconditioner ( const GridType &Grid,
                                    const BlockMatrixType &fineBlockMatrix,
                                    const int ExplicitLevel = 2,
                                    const int PreSmooth = 3, const int PostSmooth = 3,
                                    const aol::GaussSeidelSweepingMode GSmode = aol::GAUSS_SEIDEL_FORWARD )
  : CFEBlockMultigrid < CFEBlockOpType, RePrOpType, CMode > ( Grid, fineBlockMatrix, ExplicitLevel, PreSmooth, PostSmooth, GSmode, 1.0, mg::MULTIGRID_V_CYCLE, 1.0 - 1.0e-16, 1, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM ) {
    this->setVerboseMode ( 0 );
  }
};


}

#endif
