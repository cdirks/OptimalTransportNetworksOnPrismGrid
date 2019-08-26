#ifndef __ADAPTIVEPRISMMESHVALUEMAP_H
#define __ADAPTIVEPRISMMESHVALUEMAP_H


template < typename RealType >
struct IndexInfo {
  bool _hanging;                     //!< Is the node a hanging node?
  unsigned int _node;                //!< Indices of nodes for prolongation/restriction.
  RealType _prolongationWeight;      //!< Weights for prolongation/restriction.
};

   
//! \brief Helper class for transferring a function from one grid to another.
//! \todo 64 bit hash key length is hard coded
//! \author Rossmanith
template < typename RealType >
class AdaptivePrismMeshValueMap {
    
  typedef unsigned long int KeyFullType; 
  typedef std::unordered_map < KeyFullType, IndexInfo < RealType > > MapTypeNew;
    
public:
    
  //! \brief Transfer values from a vector to the value map.
  //! \param[in] The restriction/prolongation operator of the grid.
  //! \param[in] vector Vector of nodal values of the function.
  template < typename RPOpType >
  unsigned int getIndexFromNodes ( RPOpType &rp ) { 
    _indexMap.clear();
    // Go through all nodes in current grid 
    for ( typename RPOpType::GridType::NodeIterator it ( rp.getGridRef() ); it.notAtEnd(); ++it ) {
      // Create new indexInfo 
      IndexInfo < RealType > indexInfo;
      auto mapit = rp.getIndexMap().find ( it.getIndex() );
      indexInfo._hanging = mapit->second._hanging;
      //indexInfo._node = rp.getRestrictedGlobalNodeIndex ( mapit->second._node ); 
      indexInfo._node = mapit->second._node; // ATTENTION: Sets node indices, NOT dof indices 
      indexInfo._prolongationWeight = mapit->second._prolongationWeight; 
      _indexMap [ rp.getGridRef().convertCoord3DToHash ( it.getCoords() ) ] = indexInfo;     
      //_indexMap[it.getIndex()] = indexInfo;
    }
    return rp.getGridRef().getNumberOfNodes();
  }
  
  bool indexMapSet () const {
    return _indexMap.size() != 0;   
  }

  //! \brief Interpolation matrix 
  template < typename RPOpType, typename MatrixType > 
  void makeInterpolationMatrixZConst ( RPOpType &rp, MatrixType &matrix ) { 
    aol::BitVector visited ( rp.getGridRef().getNumberOfDofs() );
    visited.setAll ( false );  
    // Go through all elements in new grid 
    for ( typename RPOpType::GridType::ElementIterator it ( rp.getGridRef() ); it.notAtEnd(); ++it ) {
      // Go through all nodes in current element 
      aol::Vec < 6, GlobalIndex > nodeInds =  it.getNodeIndices();  
      for ( int i = 0; i < 6; ++i ) {  
            
        // Check if node is non-hanging and was not visited 
        if ( !rp.getGridRef().isHangingNode ( nodeInds[i] ) && !visited.get ( rp.getRestrictedGlobalNodeIndex ( nodeInds[i] ) ) ) {
          // Check if node existed before in old grid 
          auto got = _indexMap.find ( rp.getGridRef().convertCoord3DToHash ( rp.getGridRef().getNodeCoords ( nodeInds[i] ) ) );  
          // Node existed before (might have been hanging before): Interpolate from old interpolation nodes 
          if ( got != _indexMap.end() ) {  
            //matrix.add ( rp.getRestrictedGlobalNodeIndex ( nodeInds[i] ), got->second._node, got->second._prolongationWeight ); 
            matrix.add ( nodeInds[i], got->second._node, got->second._prolongationWeight );  // ATTENTION: Sets node indices, NOT dof indices 
          }
          // Node did not exist before: Interpolate from grid interpolation map   
          else {             
            auto interpol = rp.getGridRef().getInterpolationMap().find ( nodeInds[i] ); 
            if ( interpol == rp.getGridRef().getInterpolationMap().end() ) { 
              throw aol::Exception ( "AdaptivePrismMeshValueMap<>::makeInterpolationMatrixZConst(): New node is not in interpolation map, you did something stupid!", __FILE__, __LINE__ );   
            }
            for ( unsigned int j = 0; j < interpol->second._nodes.size(); ++j ) {
              //matrix.add ( rp.getRestrictedGlobalNodeIndex ( nodeInds[i] ), interpol->second._nodes[j], interpol->second._weights[j] );   
              matrix.add ( nodeInds[i], interpol->second._nodes[j], interpol->second._weights[j] );   // ATTENTION: Sets node indices, NOT dof indices 
            }
          }
          // Mark node as visited 
          visited.set ( rp.getRestrictedGlobalNodeIndex ( nodeInds[i] ), true );  
        }        
        
      }
    }
  }
    
        
protected: 
    
  MapTypeNew _indexMap;

};

#endif