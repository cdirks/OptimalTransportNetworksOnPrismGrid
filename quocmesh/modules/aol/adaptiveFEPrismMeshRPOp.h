#ifndef __ADAPTIVEFEPRISMMESHRPOP_H
#define __ADAPTIVEFEPRISMMESHRPOP_H

#include <set>

#if defined ( USE_CPP11 ) || defined ( _MSC_VER )
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

#include <op.h>
#include <qmException.h>
#include <multiVector.h>
#include <vectorExtensions.h>
#include <sparseMatrices.h>
#include <sparseMatrixRowIterator.h>
#include <pointerClasses.h>


/**
 * Restriction/prolongation operator class for adaptive three-dimensional prism grid used for problem arising from functional lifting 
 * \author Dirks 
 */


template < typename RealType >
struct AdaptiveOpRPInfo {
  bool _hanging;                                  //!< Is the node a hanging node?
  unsigned int _node;                             //!< Indices of nodes for prolongation/restriction.
  RealType _prolongationWeight;                   //!< Weights for prolongation/restriction.
};


template < typename ConfiguratorType >
class AdaptiveFEPrismMeshRPOp {

public: 
    
  typedef typename ConfiguratorType::MapType MapType;
  typedef typename ConfiguratorType::InitType GridType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef aol::Op < typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType > OpType;
  typedef std::unordered_map < uint64_t, AdaptiveOpRPInfo < RealType > > IndexMapType;
  
  
  class DofIterator {
      
    typedef AdaptiveFEPrismMeshRPOp < ConfiguratorType > RPType;
    
  protected: 
      
    RPType _rp;
    unsigned int _iter;
    
  public: 
      
    DofIterator ( AdaptiveFEPrismMeshRPOp < ConfiguratorType > & rp ) : _rp ( rp ), _iter ( 0 ) {}
    
    bool atEnd () {
      return static_cast < int > (_iter) == _rp.getGridRef().getNumberOfNodes();
    }
    
    bool notAtEnd () {
      return !this->atEnd();   
    }
      
    DofIterator & operator++ () {
      do {
        ++_iter;
        // If iterator is at end, do not compute anything 
        if ( this->atEnd() ) {
          return *this; 
        }
      } while ( _rp.getGridRef().isHangingNode ( _iter ) );
      return *this;        
    }  
    
    inline int getNodeIndex () const {
      return static_cast < int > ( _iter );   
    }
      
    inline int getDofIndex () const {        
      return _rp.getRestrictedGlobalNodeIndex ( _iter );
    }
    
    const aol::Vec3 < RealType > getCoords() {
      return _rp.getGridRef().getNodeCoords ( _iter );
    }
    
    
      
  };  
  

  AdaptiveFEPrismMeshRPOp ( const GridType & grid, const ConfiguratorType & configurator )
  : _grid ( grid ), _configurator ( configurator ), _hangingNodesBefore ( _grid.getNumberOfNodes() ), _op ( NULL ) {
      
    // Get hanging node map   
    MapType & hangingNodeMap = _grid.getHangingNodeMap();

    // Set hangingNodesBefore 
    int numHangingNodes = 0;
    for ( int i = 0; i < _grid.getNumberOfNodes(); ++i ) {
      _hangingNodesBefore[i] = numHangingNodes;  
      if ( _grid.isHangingNode ( i ) ) ++numHangingNodes;
    }

    // Set index map    
    this->generateIndexMapZConst ( grid, hangingNodeMap, _indexMap );    
    
  }
  
  // Copy constructor 
  AdaptiveFEPrismMeshRPOp ( AdaptiveFEPrismMeshRPOp & other ) :
    _grid ( other._grid ),
    _configurator ( other._configurator ),
    _hangingNodesBefore ( other._hangingNodesBefore ),
    _indexMap ( other._indexMap ) {}
  
  const GridType& getGridRef ( ) {
    return _grid;
  }

  const ConfiguratorType& getConfiguratorRef ( ) {
    return _configurator;
  }

  const IndexMapType& getIndexMap () const {
    return _indexMap;
  }

  const aol::Vector < unsigned int > & getHangingNodesBefore () const {
    return _hangingNodesBefore;
  }
  
  //! Generate the mapping from indices to AdaptiveOpRPInfos
  static void generateIndexMapZConst ( const GridType & grid, const MapType & hangingNodeMap, IndexMapType & indexMap )  { 
    RealType zCoord1, zCoord2;  
    // Go through all grid nodes 
    for ( int nodeit = 0; nodeit < grid.getNumberOfNodes(); ++nodeit ) {
      AdaptiveOpRPInfo < RealType > rpInfo;  
      auto got = hangingNodeMap.find ( nodeit );
      // Node is non-hanging
      if ( got == hangingNodeMap.end() ) {
        rpInfo._hanging = false;
        rpInfo._node = nodeit;
        rpInfo._prolongationWeight = 1.0;
        indexMap[nodeit] = rpInfo;
      }
      // Node is hanging 
      else {
        rpInfo._hanging = true; 
        // Insert lower interpolation node for z constant prolongation 
        zCoord1 = grid.getNodeCoord ( got->second._interpolationNodes.first, 2 );
        zCoord2 = grid.getNodeCoord ( got->second._interpolationNodes.second, 2 );
        if ( zCoord1 < zCoord2 ) {
          rpInfo._node = got->second._interpolationNodes.first;
        }
        else {
          rpInfo._node = got->second._interpolationNodes.second;   
        }
        rpInfo._prolongationWeight = 1.0;
        indexMap[nodeit] = rpInfo;
      }      
    }
  }  
  
  //! \brief Returns the global index of a node, when hanging nodes are not included in the numbering.
  //! \warning Only works on non-hanging nodes!
  inline int getRestrictedGlobalNodeIndex ( int nodeIndex ) const {
    return nodeIndex - _hangingNodesBefore[nodeIndex];
  }
  
  void extendVectorZConst ( VectorType &arg ) const {
    VectorType origArg ( arg, aol::DEEP_COPY );
    arg.resize ( _grid.getNumberOfNodes () );
    this->extendVectorZConst ( origArg, arg );
  }
  
  void extendVectorZConst ( const VectorType &Arg, VectorType& Dest ) const {
    Dest.setZero ();
    for ( auto it = _indexMap.begin (); it != _indexMap.end (); ++it ) {
      // Combine hanging node values from neighbors...
      if ( it->second._hanging == true ) {
        Dest[it->first] += it->second._prolongationWeight * Arg[it->second._node - _hangingNodesBefore[it->second._node]];    
      }
      // ...and keep dof values.
      else {
        Dest[it->first] = Arg[it->first - _hangingNodesBefore[it->first]];   
      }
    }
  }
  
  void restrictVectorZConst ( VectorType &arg ) const {
    VectorType origArg ( arg, aol::DEEP_COPY );
    arg.resize ( _grid.getNumberOfDofs () );
    arg.setZero ();
    this->restrictVectorZConst ( origArg, arg );
  }
  
  void restrictVectorZConst ( const VectorType& Arg, VectorType& Dest ) const {
    for ( auto it = _indexMap.begin (); it != _indexMap.end (); ++it ) {
      // Throw away hanging nodes.
      if ( it->second._hanging == false ) {
        Dest[it->second._node - _hangingNodesBefore[it->second._node]] = Arg[it->first];  
      }
    }  
  }
   
  // Apply PTP in place 
  template < typename MatrixType > 
  void applyPTPInPlaceZConstSparse ( MatrixType &m ) {
    for ( typename IndexMapType::iterator it = _indexMap.begin (); it != _indexMap.end (); ++it ) {
      if ( it->second._hanging == true ) {
        for ( int j = 0; j < m.getNumCols (); ++j ) {
          m.add ( it->second._node, j, it->second._prolongationWeight * m.get ( it->first, j ) );  
          m.set ( it->first, j, 0.0 ); 
        }
      }
    }
    for ( typename IndexMapType::iterator it = _indexMap.begin (); it != _indexMap.end (); ++it ) {
      if ( it->second._hanging == true ) {
        for ( int i = 0; i < m.getNumRows (); ++i ) {
          m.add ( i, it->second._node, it->second._prolongationWeight * m.get ( i, it->first ) );
          m.set ( i, it->first, 0.0 );
        }
      }
    }  
  }
  
  // Apply PTP in place for three matrices at once (matrices must have the same size) 
  template < typename MatrixType > 
  void applyPTPInPlaceZConstSparse ( MatrixType &m1, MatrixType &m2, MatrixType &m3 ) {
    if ( m1.getNumRows() != m2.getNumRows() || m2.getNumRows() != m3.getNumRows() || m1.getNumRows() != m3.getNumRows() ) 
      throw aol::Exception ( "AdaptiveFEPrismMeshRPOp::applyPTPInPlaceZConstSparse ( m1, m2, m3 ): Matrices do not have the same size", __FILE__, __LINE__ );
    for ( typename IndexMapType::iterator it = _indexMap.begin (); it != _indexMap.end (); ++it ) {
      if ( it->second._hanging == true ) {
        for ( int j = 0; j < m1.getNumCols (); ++j ) {
          m1.add ( it->second._node, j, it->second._prolongationWeight * m1.get ( it->first, j ) );  
          m1.set ( it->first, j, 0.0 ); 
          m2.add ( it->second._node, j, it->second._prolongationWeight * m2.get ( it->first, j ) );  
          m2.set ( it->first, j, 0.0 ); 
          m3.add ( it->second._node, j, it->second._prolongationWeight * m3.get ( it->first, j ) );  
          m3.set ( it->first, j, 0.0 ); 
        }
      }
    }
    for ( typename IndexMapType::iterator it = _indexMap.begin (); it != _indexMap.end (); ++it ) {
      if ( it->second._hanging == true ) {
        for ( int i = 0; i < m1.getNumRows (); ++i ) {
          m1.add ( i, it->second._node, it->second._prolongationWeight * m1.get ( i, it->first ) );
          m1.set ( i, it->first, 0.0 );
          m2.add ( i, it->second._node, it->second._prolongationWeight * m2.get ( i, it->first ) );
          m2.set ( i, it->first, 0.0 );
          m3.add ( i, it->second._node, it->second._prolongationWeight * m3.get ( i, it->first ) );
          m3.set ( i, it->first, 0.0 );
        }
      }
    }  
  }
  
  
  void fastRestrictSparseMatrix ( aol::SparseMatrix < RealType > &m ) const {
      
    // Remove rows
    std::vector < typename IndexMapType::key_type > hangingNodesOrdered;
    for ( typename IndexMapType::iterator it = _indexMap.begin (); it != _indexMap.end (); ++it ) {
      if ( it->second._hanging == true ) {
        hangingNodesOrdered.push_back ( it->first );
      }
    }
    std::sort ( hangingNodesOrdered.begin (), hangingNodesOrdered.end () );
    for ( unsigned int i = 0; i < hangingNodesOrdered.size (); ++i )
      m.destroyRow ( hangingNodesOrdered[i] - _hangingNodesBefore[hangingNodesOrdered[i]] );

    // Remove columns
    for ( int i = 0; i < m.getNumRows (); ++i ) {
      for ( aol::SparseMatrixRowIterator < RealType > rowIt ( m, i ); rowIt.notAtEnd (); ++rowIt ) {
        // If not hanging, set i ~> i - hanging nodes before i.
        if ( _indexMap[rowIt->col]._hanging == false ) {
          rowIt->col -= _hangingNodesBefore[rowIt->col];
        }
      }
    }
    m.destructiveResize ( _grid.getNumberOfDofs (), _grid.getNumberOfDofs () );
    
  }
  
  
    
protected:
    
  const GridType &_grid;
  const ConfiguratorType &_configurator;
  aol::Vector < unsigned int > _hangingNodesBefore;  //!< _hangingNodesBefore[i]: The number of hanging nodes before node i.
  mutable OpType const* _op;  //!< The operator on the full grid (including hanging nodes).
  mutable IndexMapType _indexMap;  //!< Maps index in between prolongated and non-prolongated vectors.    

};

#endif