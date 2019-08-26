#ifndef __ADAPTIVEPRISMMESH_H
#define __ADAPTIVEPRISMMESH_H

#include <aol.h>
#include <vec.h>
#include <math.h>
#include <multiVector.h>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <array>

// TODO: Change every int to GlobalIndex or LocalIndex 


/**
 * Adaptive three-dimensional prism grid class used for problem arising from functional lifting 
 * \author Dirks 
 */


template < typename GlobalIndex > 
struct HangingNodeInfo {
  std::unordered_set < GlobalIndex > _neighbors; 
  std::pair < GlobalIndex, GlobalIndex > _interpolationNodes; 
};

template < typename GlobalIndex, typename RealType > 
struct InterpolationInfo {
  std::vector < GlobalIndex > _nodes;
  std::vector < RealType > _weights;
};


template  < typename RealType, typename GlobalIndex > 
class PrismBaseElement {
  
protected: 

  GlobalIndex _globIdx;                     // Global indices of element
  aol::Vec < 6, GlobalIndex > _globNodeIdx; // Global indices of nodes 
  aol::Matrix22 < RealType > _ginv;         // Inverse transformation matrix of theta:T_ref->T (=Inverse Metric)
  RealType _vol;                            // Volume
  aol::Vec3 < RealType > _coords[6];        // Node coords 
  RealType _baseArea;                       // Triangle area of element base 
  unsigned short _levelz;                   // z level of element 
    
public:
    
  PrismBaseElement () : _globIdx ( -1 ) {}

  template < typename MeshType >
  PrismBaseElement ( const MeshType & Mesh, int globalIdx ) : _globIdx ( globalIdx ) {
    fillElement ( Mesh, globalIdx );
  }

  ~PrismBaseElement(){}  
    
  // Fill Element function
  template < typename MeshType >
  void fillElement ( const MeshType & Mesh, const int globalIdx ) {
    _globIdx = globalIdx;  
    // Fill geometry corresponding to StorageType
    for ( int i = 0; i < 6; ++i ) {
      _globNodeIdx[i] =  Mesh.getNodeIndicesOfElement ( _globIdx )[i];
      _coords[i] = Mesh.getNodeCoords ( _globNodeIdx[i] );
    }    
    _levelz = std::log2 ( 1.0/this->getHeight() + 0.5 );
    // Compute Inverse Metric and save
    setInverseMetric ();
  }  

  const aol::Matrix22 < RealType > & ginv() const {
    return _ginv;
  }  
  
  inline RealType vol ( ) const { 
    return _vol;
  }
  
  inline GlobalIndex getGlobalElementIndex () const {
    return _globIdx;   
  }
  
  inline RealType getHeight ( ) const {
    return _coords[3][2] - _coords[0][2];   
  }
  
  inline RealType getBaseArea ( ) const {
    return _baseArea;    
  }
  
  inline unsigned short getLevelz ( ) const {
    return _levelz;    
  }
  
  const GlobalIndex getGlobalNodeIndex ( const int localIndex ) const {
    return _globNodeIdx[localIndex];
  }
  
  aol::Vec3 < RealType > edge ( const int locNode1, const int locNode2 ) const {
    return aol::Vec3 < RealType > ( _coords[locNode2] - _coords[locNode1] );
  }
  
protected: 
    
  void setInverseMetric() {
    aol::Vec3 < RealType > dx01 = this->edge ( 0, 1 );
    aol::Vec3 < RealType > dx02 = this->edge ( 0, 2 );
    RealType g11 = ( dx01[0] * dx01[0] + dx01[1] * dx01[1] );
    RealType g12 = ( dx01[0] * dx02[0] + dx01[1] * dx02[1] );
    RealType g22 = ( dx02[0] * dx02[0] + dx02[1] * dx02[1] );
    RealType detg = g11 * g22 - aol::Sqr ( g12 );
    if ( detg < 1e-15 ) throw aol::Exception ( "aol::PrismBaseElement::setInverseMetric: determinant not positive.\n", __FILE__, __LINE__ );
    _ginv.set ( 0, 0, g22/detg );
    _ginv.set ( 1, 1, g11/detg);
    _ginv.set ( 1, 0, -g12/detg);
    _ginv.set ( 0, 1, -g12/detg);
    _baseArea = 0.5 * sqrt ( detg );
    _vol = _baseArea * this->getHeight();
  }
    
    
};


template < typename RealType, typename PrismType = PrismBaseElement < RealType, int > >
class AdaptiveFEPrismMesh {
    
  typedef std::vector < aol::Vec3 < int > > CellArray;  
  typedef std::vector < aol::Vec < 6, int > > CellArray3D;  

public: 
    
  class ElementIterator {
      
  public:
      
    typedef AdaptiveFEPrismMesh < RealType > MeshType;
    typedef ElementIterator EndType;
    typedef PrismType IteratedType;
      
  protected: 
      
    const MeshType &  _mesh;  
    unsigned int _iter; 
    mutable IteratedType _currentEl;
    mutable bool _filled;
    
  public:
      
    ElementIterator ( const MeshType & mesh ) : _mesh ( mesh ), _iter ( 0 ), _filled ( false ) {}
    
    ElementIterator ( const MeshType & mesh, const int iter ) : _mesh ( mesh ), _iter ( iter ), _filled ( false ) {}
    
    void set ( const int iter ) {
      _iter = iter; 
      _filled = false; 
    }
    
    bool operator!= ( const ElementIterator &otherIt ) const {
      return ( this->_iter != otherIt._iter );
    }
    
    ElementIterator & operator++ () {
      ++_iter;
      _filled = false;
      return *this;
    }
    
    IteratedType & operator*() {
      if (!_filled) fillElement();
      return _currentEl;
    }
    
    RealType getBaseArea () const {
      if (!_filled) fillElement();  
      return _currentEl.getBaseArea();   
    }
    
    unsigned short getLevelz () const {
      if (!_filled) fillElement();  
      return _currentEl.getLevelz();    
    }
    
    bool atEnd () const {
      return static_cast < int > (_iter) == _mesh.getNumberOfElements();
    }
    
    bool notAtEnd () const {
      return ! atEnd ();
    }
    
    int getIndex () const {
      return static_cast < int > ( _iter );
    }
    
    aol::Vec < 6, int > getNodeIndices() const {
      return _mesh.getNodeIndicesOfElement ( _iter );
    }
    
  protected:
      
    void fillElement () const {
      _currentEl.fillElement ( _mesh, _iter );
      _filled = true;
    }
    
  };
  
  class NodeIterator { 
  
    typedef AdaptiveFEPrismMesh < RealType > MeshType;
    
  protected: 
      
    const MeshType & _mesh;
    unsigned int _iter;
    
  public:
    
    NodeIterator ( const MeshType & mesh ) : _mesh ( mesh ), _iter ( 0 ) {}
    
    NodeIterator & operator++ () {
      ++_iter;
      return *this;
    }
    
    bool atEnd () const {
      return static_cast < int > (_iter) == _mesh.getNumberOfNodes();
    }
    
    bool notAtEnd () const {
      return ! atEnd ();
    }
    
    int getIndex () const {
      return static_cast < int > ( _iter );
    }
    
    aol::Vec3 < RealType > getCoords() const {
      return _mesh.getNodeCoords ( _iter );
    }
      
  };
  
  typedef ElementIterator ElementIteratorType;
  typedef NodeIterator NodeIteratorType;
  typedef short LocalIndex;
  typedef int GlobalIndex;
  typedef std::unordered_map < GlobalIndex, HangingNodeInfo < GlobalIndex > > MapType; 

  // Empty constructor 
  AdaptiveFEPrismMesh() : _maxlevelz ( 0 ), _nodes3D ( 0 ), _elements3D ( 0 ), _markedForXYRefinement ( 0 ), _markedForZRefinement ( 0 ), _neighborList ( 0 ), _visited ( 0 ), _neighborListTmp ( 0 ), begin_it ( *this ), end_it ( *this ) {}
  
  // Level constructor
  AdaptiveFEPrismMesh ( int levelxy, int levelz ) :
    begin_it ( *this ), 
    end_it ( *this, ( 1 << ( 2*levelxy + 1 ) ) * ( 1 << levelz ) ),
    _maxlevelz ( levelz ),
    _nodes3D ( 3, ( ( 1 << levelxy ) + 1 ) * ( ( 1 << levelxy ) + 1 ) * ( ( 1 << levelz ) + 1) ), 
    _markedForXYRefinement ( ( 1 << ( 2*levelxy + 1 ) ) * ( 1 << levelz ) ), 
    _markedForZRefinement ( ( 1 << ( 2*levelxy + 1 ) ) * ( 1 << levelz ) ),
    _neighborList ( ( 1 << ( 2*levelxy + 1 ) ) * ( 1 << levelz ) ) {
        
    int numGroundNodes = ( ( 1 << levelxy ) + 1 ) * ( ( 1 << levelxy ) + 1 );
    int numGroundElements = 1 << ( 2*levelxy + 1 );
    
    // Node list of 3D nodes 
    int index;
    int n = ( 1 << levelxy ) + 1;
    RealType stepZ = 1.0 / static_cast < RealType > ( 1 << levelz );
    aol::Vec3 < RealType > coords;
    RealType stepXY = 1.0 / static_cast < RealType > ( 1 << levelxy );
    for ( int i = 0; i < ( 1 << levelz ) + 1; ++i ) {
      for ( int j = 0; j < numGroundNodes; ++j ) {
        coords[0] = static_cast < RealType > ( ( j % n ) ) * stepXY; 
        coords[1] = static_cast < RealType > ( ( j - ( j % n ) ) / n ) * stepXY;
        coords[2] = static_cast < RealType >(i) * stepZ;
        index = i * numGroundNodes + j;
        _nodes3D.set ( index, coords ); 
        std::pair < unsigned long int, int > newHash ( convertCoord3DToHash ( coords ), index );
        _nodeMap.insert ( newHash ); 
      }
    }
    
    // Element list of 3D elements 
    CellArray groundTriang ( numGroundElements );
    CellArray3D elements3D ( numGroundElements * ( 1 << levelz ) );
    aol::Vec < 6, int > el;
    aol::Vector < int > aux ( 1 << ( levelxy+1 ) );
    aux[0] = 1; aux[1] = n;
    for ( int i = 2; i < ( 1 << ( levelxy+1 ) ); ++i ) {
      aux[i] = aux[i-2] + 1;
    }
    unsigned int numEls = 1 << ( 2*levelxy + 1 );
    int m = 1 << ( levelxy + 1 );
    for ( int i = 0; i < ( 1 << levelz ); ++i ) {
      for ( unsigned int j = 0; j < numEls; ++j ) {
        int p = j % m;
        int r = (j-p)/m;
        index = i * numGroundElements + j;
        if ( index % 2 == 0 ) {
          el[0] = r*n + floor((p+1)/2) + i*numGroundNodes; el[1] = r*n + aux[p] + i*numGroundNodes; el[2] = (r+1)*n + floor((p+1)/2) + i*numGroundNodes; 
          el[3] = r*n + floor((p+1)/2) + (i+1)*numGroundNodes; el[4] = r*n + aux[p] + (i+1)*numGroundNodes; el[5] = (r+1)*n + floor((p+1)/2) + (i+1)*numGroundNodes;
        }
        else {
          el[0] = (r+1)*n + floor((p+1)/2) + i*numGroundNodes; el[1] = r*n + aux[p] + i*numGroundNodes; el[2] = r*n + floor((p+1)/2) + i*numGroundNodes;
          el[3] = (r+1)*n + floor((p+1)/2) + (i+1)*numGroundNodes; el[4] = r*n + aux[p] + (i+1)*numGroundNodes; el[5] = r*n + floor((p+1)/2) + (i+1)*numGroundNodes;
        } 
        elements3D[index] = el;
      }
    }
    
    // Set neighbor list 
    int pTilde, xTilde, yTilde; int mTilde = 1 << ( levelxy + 1 ); 
    aol::Vec < 8, int > neighbor1; 
    neighbor1[0] = 1; neighbor1[1] = -1; neighbor1[2] = 1-mTilde; neighbor1[3] = neighbor1[0]; neighbor1[4] = neighbor1[1]; neighbor1[5] = neighbor1[2]; neighbor1[6] = -2; neighbor1[7] = numEls; 
    aol::Vec < 8, int > neighbor2; 
    neighbor2[0] = 0; neighbor2[1] = 2; neighbor2[2] = mTilde; neighbor2[3] = neighbor2[0]; neighbor2[4] = neighbor2[1]; neighbor2[5] = neighbor2[2]; neighbor2[6] = -2; neighbor2[7] = numEls+1;
    _neighborList[0] = neighbor1;
    _neighborList[1] = neighbor2;    
    for ( unsigned int i = 2; i < elements3D.size(); ++i ) {
      for ( int j = 0; j < 6; ++j ) {
        _neighborList[i][j] = _neighborList[i-2][j] + 2;
      }
      // Element below 
      if ( i < numEls ) {
        _neighborList[i][6] = -2;    
      }
      else { 
        _neighborList[i][6] = i - numEls;  
      }
      // Element above 
      if ( i >= elements3D.size() - numEls ) {
        _neighborList[i][7] = -2;    
      }
      else {
        _neighborList[i][7] = i + numEls;   
      }          
    }
    if ( levelz == 0 ) {
      _neighborList[0][7] = -2; _neighborList[1][7] = -2;
    }
    // Set indices outside of domain to -2
    // ATTENTION: _neighborList = -2 means that neighbor is not in domain, _neighborList = -1 means that edge does not appear in neighbor since neighbor.zlevel is greater than current element zlevel
    for ( unsigned int i = 0; i < elements3D.size(); ++i ) {
      if ( i % 2 == 0 ) {  
        pTilde = i % mTilde;
        xTilde = (i-pTilde)/mTilde;
        yTilde = xTilde % (mTilde/2);
        // Second and fifth column 
        if ( pTilde == 0 || pTilde == mTilde - 1 ) {
          _neighborList[i][1] = -2; _neighborList[i][4] = -2; 
        }
        // Third and sixth column   
        if ( yTilde == 0 ) { 
          _neighborList[i][2] = -2; _neighborList[i][5] = -2;
        }      
      }
      else {
        pTilde = i % mTilde;
        xTilde = (i-pTilde)/mTilde;
        yTilde = xTilde % (mTilde/2);
        // Third and sixth column 
        if ( yTilde == mTilde/2 - 1 ) { 
          _neighborList[i][2] = -2; _neighborList[i][5] = -2; 
        }   
        // Second and fifth column 
        if ( pTilde == 0 || pTilde == mTilde - 1 ) {
          _neighborList[i][1] = -2; _neighborList[i][4] = -2; 
        }  
      }
    }

    //_elements3D.reset ( new CellArray3D ( elements3D ), true );
    _elements3D = elements3D; 
    
    this->unmarkAll();
    
  }
  
  // Copy constructor 
  AdaptiveFEPrismMesh ( const AdaptiveFEPrismMesh & other ) : 
    begin_it ( *this ),
    end_it ( *this, other.getNumberOfElements() ),
    _maxlevelz ( other._maxlevelz ),
    _nodes3D ( other._nodes3D ), 
    _nodeMap ( other._nodeMap ), 
    _elements3D ( other._elements3D ), 
    _markedForXYRefinement ( other._markedForXYRefinement ), 
    _markedForZRefinement ( other._markedForZRefinement ), 
    _interpolationMap ( other._interpolationMap ), 
    _neighborList ( other._neighborList ), 
    _hangingNodeMap ( other._hangingNodeMap ), 
    _subelementList ( other._subelementList ),  
    _sizeBeforeRefinement ( other._sizeBeforeRefinement ) {      
    }
  
  // Operator = 
  AdaptiveFEPrismMesh & operator= ( const AdaptiveFEPrismMesh & other ) {
    _maxlevelz = other._maxlevelz;
    _nodes3D = other._nodes3D; 
    _nodeMap = other._nodeMap; 
    _elements3D = other._elements3D; 
    _markedForXYRefinement = other._markedForXYRefinement; 
    _markedForZRefinement = other._markedForZRefinement; 
    _interpolationMap = other._interpolationMap;
    _neighborList = other._neighborList; 
    _hangingNodeMap = other._hangingNodeMap; 
    _subelementList = other._subelementList; 
    _sizeBeforeRefinement = other._sizeBeforeRefinement;
    begin_it = other.begin_it; 
    end_it = other.end_it;
    return *this;
  }

  const ElementIterator & begin () const {
    return begin_it;
  }

  const ElementIterator & end () const {
    return end_it;
  }
  
  // 3D coordinate hashing 
  unsigned long int convertCoord3DToHash ( const aol::Vec3 < RealType > & coord, const unsigned long int& maxLevel = 10UL ) const {
    return ( ( 1UL << (maxLevel) ) + 1UL ) * ( ( 1UL << (maxLevel) ) + 1UL ) * static_cast<unsigned long int>( coord[2] * static_cast<RealType>( 1UL << (maxLevel) ) + static_cast<RealType>( 0.5 ) ) + ( ( 1UL << (maxLevel) ) + 1UL ) * static_cast<unsigned long int>( coord[1] * static_cast<RealType>( 1UL << (maxLevel) ) + static_cast<RealType>( 0.5 ) ) + static_cast<unsigned long int>( coord[0] * static_cast<RealType>( 1UL << (maxLevel) ) + static_cast<RealType>( 0.5 ) );
  }
  
  // 2D coordinate hashing 
  unsigned long int convertCoord2DToHash ( const aol::Vec3 < RealType > & coord, const unsigned long int& maxLevel = 10UL ) const {
    return ( ( 1UL << (maxLevel) ) + 1UL ) * static_cast<unsigned long int>( coord[1] * static_cast<RealType>( 1UL << (maxLevel) ) + static_cast<RealType>( 0.5 ) ) + static_cast<unsigned long int>( coord[0] * static_cast<RealType>( 1UL << (maxLevel) ) + static_cast<RealType>( 0.5 ) );
  }
  
  inline int getMaxLevelz () const {
    return _maxlevelz;
  }
  
  // Unmark all elements
  void unmarkAll () {
    for ( unsigned int i = 0; i < _markedForXYRefinement.size(); ++i ) {
      _markedForXYRefinement[i] = false;
    }
    for ( unsigned int i = 0; i < _markedForZRefinement.size(); ++i ) {
      _markedForZRefinement[i] = false;
    }
  }
  
  inline aol::MultiVector < RealType > & getNodes3D () {
    return _nodes3D;
  }
  
  const inline std::unordered_map < unsigned long int, GlobalIndex > & getNodeMap () const {
    return _nodeMap;   
  }
  
  CellArray3D getElements3DConst () const {
    return _elements3D;   
  }
  
  const CellArray3D & getElements3D () const {
    return _elements3D;   
  }
  
  const std::map < int, InterpolationInfo < GlobalIndex, RealType > > & getInterpolationMap () const {
    return _interpolationMap;
  }
  
  const std::map < GlobalIndex, std::vector < GlobalIndex > > & getSubelementList () const {
    return _subelementList;
  }
    
  int getNumberOfElements () const {
    return _elements3D.size();   
  }
  
  int getNumberOfNodes () const {
    return _nodes3D[0].size();
  }
  
  int getNumberOfDofs () const {
    return this->getNumberOfNodes() - _hangingNodeMap.size();  
  }
  
  int getNumberOfHangingNodes () const {
    return _hangingNodeMap.size();   
  }
  
  std::vector < aol::Vec < 8, int > > & getNeighborList () {
    return _neighborList;
  }
  
  inline MapType & getHangingNodeMap () const {
    return _hangingNodeMap;    
  }
  
  aol::Vec < 8, GlobalIndex > & getNeighborsOfElement ( const GlobalIndex elIndex ) {
    return _neighborList[elIndex];   
  }
  
  const aol::Vec < 6, int > & getNodeIndicesOfElement ( const int num ) const {
    return ( _elements3D[num] );
  }
  
  const aol::Vec3 < int > getGroundNodeIndicesOfElement ( const int num ) const {
    return aol::Vec3 < int > ( _elements3D[num][0], _elements3D[num][1], _elements3D[num][2] );
  }
  
  RealType getNodeCoord ( const int nodeIndex, const int numCoord ) const {
    return _nodes3D[numCoord][nodeIndex];
  }
  
  aol::Vec3 < RealType > getNodeCoords ( const int nodeIndex ) const {
    aol::Vec3 < RealType > coords ( this->getNodeCoord(nodeIndex,0), this->getNodeCoord(nodeIndex,1), this->getNodeCoord(nodeIndex,2) );
    return coords;   
  }
  
  GlobalIndex getElementNodeIndex ( GlobalIndex element, LocalIndex localInd ) const {
    return ( _elements3D[element][localInd] );
  }
  
  // Mark for xy refinement
  void markXY ( int element ) {
    if( element >= this->getNumberOfElements() )
      throw aol::Exception ( "AdaptiveFEPrismMesh<>::markXY(): out of range!", __FILE__, __LINE__ );    
    _markedForXYRefinement[ element ] = true;
  }
  
  // Mark for z refinement
  void markZ ( int element ) {
    if( element >= this->getNumberOfElements() )
      throw aol::Exception ( "AdaptiveFEPrismMesh<>::markZ(): out of range!", __FILE__, __LINE__ );    
    _markedForZRefinement[ element ] = true;
  }
  
  // Unmark for xy refinement
  void unmarkXY ( int element ) {
    _markedForXYRefinement[ element ] = false;
  }
  
  // Unmark for z refinement 
  void unmarkZ ( int element ) {
    _markedForZRefinement[ element ] = false;  
  }
  
  // Check if element is marked for xy refinement 
  bool isMarkedForXYRefinement( int element ) const {
    return _markedForXYRefinement[element];
  }
  
  // Check if element is marked for z refinement 
  bool isMarkedForZRefinement( int element ) const {
    return _markedForZRefinement[element];
  }
  
  // Check if node is a hanging node.
  bool isHangingNode ( const GlobalIndex elIndex ) const {
    auto got = _hangingNodeMap.find ( elIndex );
    if ( got == _hangingNodeMap.end () )
      return false;
    else
      return true;
  }
  
  // Save in legacy VTK format
  void saveAsLegacyVTK ( string filename ) const {
      
    std::ofstream out ( filename ); 
    
    // Write the file header:
    out << "<VTKFile type=\"UnstructuredGrid\">\n";
    out << "<UnstructuredGrid>\n";
    out << "<Piece NumberOfPoints=\"" << this->getNumberOfNodes () << "\" NumberOfCells=\"" << this->getNumberOfElements () << "\">\n";
    
    // Write point data 
    out << "<Points>\n";
    out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
    std::vector < std::pair < unsigned int, aol::Vec3 < RealType > > > indices;  
    for ( NodeIteratorType nodeit ( *this ); nodeit.notAtEnd(); ++nodeit ) {
      indices.push_back  ( std::make_pair ( nodeit.getIndex(), this->getNodeCoords ( nodeit.getIndex() ) ) );
      for ( int i = 0; i < 3; ++i ) out << this->getNodeCoord ( nodeit.getIndex(), i ) << ' ';
      out << '\n';
    }
    out << "</DataArray>\n";
    out << "</Points>\n";
    
    // Write cell data 
    out << "<Cells>\n";
    out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
       out << elit.getNodeIndices()[0] << " ";
       out << elit.getNodeIndices()[1] << " ";
       out << elit.getNodeIndices()[2] << " ";
       out << elit.getNodeIndices()[3] << " ";
       out << elit.getNodeIndices()[4] << " ";
       out << elit.getNodeIndices()[5] << " ";
       out << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    
    // Offset 
    unsigned int counter = 1;
    const unsigned int off = 6;
    // Set offsets for hexahedron cells written above:
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
      out << off * counter++ << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    
    // Cell type 
    unsigned int type = 13;
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
      out << type << "\n";
    }
    out << "</DataArray>\n";
    out << "</Cells>\n";
    out << "</Piece>\n";
    out << "</UnstructuredGrid>\n";
    out << "</VTKFile>\n";

    out.close ();
    
  }

  // Save in legacy VTK format with point data
  void saveAsLegacyVTK ( string filename, const aol::Vector < RealType > & pointData ) const {
      
    std::ofstream out ( filename ); 
    
    // Write the file header:
    out << "<VTKFile type=\"UnstructuredGrid\">\n";
    out << "<UnstructuredGrid>\n";
    out << "<Piece NumberOfPoints=\"" << this->getNumberOfNodes () << "\" NumberOfCells=\"" << this->getNumberOfElements () << "\">\n";    
    
    // Extract values from pointData
    std::unordered_map < std::string, const aol::Vector < RealType >* > _pointDataMap; 
    _pointDataMap["sol"] = &pointData;
    out << "<PointData Scalars=\"";
    for ( auto it = _pointDataMap.begin (); it != _pointDataMap.end (); ++it )
      out << ( it == _pointDataMap.begin () ? "" : ", " ) << it->first;
    out << "\">\n";

    for ( auto it = _pointDataMap.begin (); it != _pointDataMap.end (); ++it ) {
      out << "<DataArray Name=\"" << it->first << "\" type=\"Float32\">\n";
      for ( int i = 0; i < it->second->size(); i++ )
        out << it->second->operator[] (i) << ' ';
      out << '\n';
      out << "</DataArray>\n";
    }
    out << "</PointData>\n";
    
    // Write point data 
    out << "<Points>\n";
    out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
    std::vector < std::pair < unsigned int, aol::Vec3 < RealType > > > indices;  
    for ( NodeIteratorType nodeit ( *this ); nodeit.notAtEnd(); ++nodeit ) {
      indices.push_back  ( std::make_pair ( nodeit.getIndex(), this->getNodeCoords ( nodeit.getIndex() ) ) );
      for ( int i = 0; i < 3; ++i ) out << this->getNodeCoord ( nodeit.getIndex(), i ) << ' ';
      out << '\n';
    }
    out << "</DataArray>\n";
    out << "</Points>\n";
    
    // Write cell data 
    out << "<Cells>\n";
    out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
       out << elit.getNodeIndices()[0] << " ";
       out << elit.getNodeIndices()[1] << " ";
       out << elit.getNodeIndices()[2] << " ";
       out << elit.getNodeIndices()[3] << " ";
       out << elit.getNodeIndices()[4] << " ";
       out << elit.getNodeIndices()[5] << " ";
       out << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    
    // Offset 
    unsigned int counter = 1;
    const unsigned int off = 6;
    // Set offsets for hexahedron cells written above:
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
      out << off * counter++ << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    
    // Cell type 
    unsigned int type = 13;
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
      out << type << "\n";
    }
    out << "</DataArray>\n";
    out << "</Cells>\n";
    out << "</Piece>\n";
    out << "</UnstructuredGrid>\n";
    out << "</VTKFile>\n";

    out.close ();  
      
  }
  
  // Save in legacy VTK format with vector data 
  void saveAsLegacyVTKVectorField ( string filename, const aol::Vector < RealType > & pointDataX, const aol::Vector < RealType > & pointDataY, const aol::Vector < RealType > & pointDataZ ) const {
      
    std::ofstream out ( filename ); 
    
    // Write the file header:
    out << "<VTKFile type=\"UnstructuredGrid\">\n";
    out << "<UnstructuredGrid>\n";
    out << "<Piece NumberOfPoints=\"" << this->getNumberOfNodes () << "\" NumberOfCells=\"" << this->getNumberOfElements () << "\">\n";    
    
    // Extract values from pointData
    std::unordered_map < std::string, const aol::Vector < RealType >* > _pointDataMapX, _pointDataMapY, _pointDataMapZ; 
    _pointDataMapX["sol"] = &pointDataX;
    _pointDataMapY["sol"] = &pointDataY;
    _pointDataMapZ["sol"] = &pointDataZ;
    out << "<PointData Vectors=\"";
    for ( auto it = _pointDataMapX.begin (); it != _pointDataMapX.end (); ++it )  
      out << ( it == _pointDataMapX.begin () ? "" : ", " ) << it->first;
    out << "\">\n";
    
    auto itY = _pointDataMapY.begin();
    auto itZ = _pointDataMapZ.begin();   
    for ( auto itX = _pointDataMapX.begin (); itX != _pointDataMapX.end (); ++itX ) {
      out << "<DataArray Name=\"" << itX->first << "\" type=\"Float32\" NumberOfComponents=\"3\" >\n";
      for ( int i = 0; i < itX->second->size(); i++ ) {
        out << itX->second->operator[] (i) << ' ' << itY->second->operator[] (i) << ' ' << itZ->second->operator[] (i);
        out << '\n'; 
      }
      out << "</DataArray>\n";
    }
    out << "</PointData>\n";
    
    // Write point data 
    out << "<Points>\n";
    out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
    std::vector < std::pair < unsigned int, aol::Vec3 < RealType > > > indices;  
    for ( NodeIteratorType nodeit ( *this ); nodeit.notAtEnd(); ++nodeit ) {
      indices.push_back  ( std::make_pair ( nodeit.getIndex(), this->getNodeCoords ( nodeit.getIndex() ) ) );
      for ( int i = 0; i < 3; ++i ) out << this->getNodeCoord ( nodeit.getIndex(), i ) << ' ';
      out << '\n';
    }
    out << "</DataArray>\n";
    out << "</Points>\n";
    
    // Write cell data 
    out << "<Cells>\n";
    out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
       out << elit.getNodeIndices()[0] << " ";
       out << elit.getNodeIndices()[1] << " ";
       out << elit.getNodeIndices()[2] << " ";
       out << elit.getNodeIndices()[3] << " ";
       out << elit.getNodeIndices()[4] << " ";
       out << elit.getNodeIndices()[5] << " ";
       out << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    
    // Offset 
    unsigned int counter = 1;
    const unsigned int off = 6;
    // Set offsets for hexahedron cells written above:
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
      out << off * counter++ << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    
    // Cell type 
    unsigned int type = 13;
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
      out << type << "\n";
    }
    out << "</DataArray>\n";
    out << "</Cells>\n";
    out << "</Piece>\n";
    out << "</UnstructuredGrid>\n";
    out << "</VTKFile>\n";

    out.close ();  
      
  }
  
  // Save in legacy VTK format with cell data 
  void saveAsLegacyVTKCell ( string filename, const aol::Vector < RealType > & cellData ) const { 
      
    std::ofstream out ( filename ); 
    
    // Write the file header:
    out << "<VTKFile type=\"UnstructuredGrid\">\n";
    out << "<UnstructuredGrid>\n";
    out << "<Piece NumberOfPoints=\"" << this->getNumberOfNodes () << "\" NumberOfCells=\"" << this->getNumberOfElements () << "\">\n";    
    
    // Extract values from cell data
    std::unordered_map < std::string, const aol::Vector < RealType >* > _cellDataMap; 
    _cellDataMap["sol"] = &cellData;
    out << "<CellData Scalars=\"";
    for ( auto it = _cellDataMap.begin (); it != _cellDataMap.end (); ++it )
      out << ( it == _cellDataMap.begin () ? "" : ", " ) << it->first;
    out << "\">\n";

    for ( auto it = _cellDataMap.begin (); it != _cellDataMap.end (); ++it ) {
      out << "<DataArray Name=\"" << it->first << "\" type=\"Float32\">\n";
      for ( int i = 0; i < it->second->size(); i++ )
        out << it->second->operator[] (i) << ' ';
      out << '\n';
      out << "</DataArray>\n";
    }
    out << "</CellData>\n";
    
    // Write cell data 
    out << "<Points>\n";
    out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
    std::vector < std::pair < unsigned int, aol::Vec3 < RealType > > > indices;  
    for ( NodeIteratorType nodeit ( *this ); nodeit.notAtEnd(); ++nodeit ) {
      indices.push_back  ( std::make_pair ( nodeit.getIndex(), this->getNodeCoords ( nodeit.getIndex() ) ) );
      for ( int i = 0; i < 3; ++i ) out << this->getNodeCoord ( nodeit.getIndex(), i ) << ' ';
      out << '\n';
    }
    out << "</DataArray>\n";
    out << "</Points>\n";
    
    // Write cells  
    out << "<Cells>\n";
    out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
       out << elit.getNodeIndices()[0] << " ";
       out << elit.getNodeIndices()[1] << " ";
       out << elit.getNodeIndices()[2] << " ";
       out << elit.getNodeIndices()[3] << " ";
       out << elit.getNodeIndices()[4] << " ";
       out << elit.getNodeIndices()[5] << " ";
       out << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    
    // Offset 
    unsigned int counter = 1;
    const unsigned int off = 6;
    // Set offsets for hexahedron cells written above:
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
      out << off * counter++ << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    
    // Cell type 
    unsigned int type = 13;
    for ( ElementIteratorType elit ( *this ); elit.notAtEnd(); ++elit ) {
      out << type << "\n";
    }
    out << "</DataArray>\n";
    out << "</Cells>\n";
    out << "</Piece>\n";
    out << "</UnstructuredGrid>\n";
    out << "</VTKFile>\n";

    out.close ();  
      
  }
  
  // Refine all marked elements 
  void refineMarkedElements () {
      
    // Clear additional node list, interpolationMap and subelementList, since new refinement round starts 
    _addNodes.clear();
    _interpolationMap.clear();
    _subelementList.clear();
    _visited.resize ( this->getNumberOfElements() );
    _visited.setAll ( false );
    _sizeBeforeRefinement = this->getNumberOfElements(); 
    
    // Set copy of _neighborList to _neighborListTmp 
    _neighborListTmp.clear();
    _neighborListTmp.resize ( _neighborList.size() );
    for ( unsigned int i = 0; i < _neighborList.size(); ++i ) {
      for ( int j = 0; j < 8; ++j ) {
        _neighborListTmp[i][j].push_back ( _neighborList[i][j] );   
      }
    }
    
    // Go through all elements and mark all elements below and above for refinement, since the refinement does not take care of such neighbors 
    int numEls = this->getNumberOfElements();
    for ( int elit = 0; elit < numEls; ++elit ) {
      if ( _markedForXYRefinement[elit] ) {
        // Find elements below 
        GlobalIndex neighborBelow = _neighborList[elit][6];
        while ( neighborBelow >= 0 ) {
          this->markXY ( neighborBelow );
          neighborBelow = _neighborList[neighborBelow][6]; 
        }
        // Find elements above 
        GlobalIndex neighborAbove = _neighborList[elit][7]; 
        while ( neighborAbove >= 0 ) {
          this->markXY ( neighborAbove );
          neighborAbove = _neighborList[neighborAbove][7];
        }
      }
    }
    
    // Go through all elements and refine XY direction 
    for ( int elit = 0; elit < numEls; ++elit ) {
      // Check if element is marked and was not visited 
      if ( _markedForXYRefinement[elit] && !_visited[elit] ) {
        refineXY ( elit );        
      }
    }
    
    // Go through all elements and refine Z direction 
    numEls = this->getNumberOfElements();
    _visited.resize ( numEls );
    _visited.setAll ( false );
    for ( int elit = 0; elit < numEls; ++elit ) {
      // Check if element is marked and was not visited 
      if ( _markedForZRefinement[elit] && !_visited[elit] ) {
        refineZ ( elit );   
      }
    }
    
    // Set _neighborList from _neighborListTmp 
    _neighborList.resize ( _neighborListTmp.size() );
    for ( unsigned int i = 0; i < _neighborListTmp.size(); ++i ) {
      //cout << "Element " << i << endl;
      for ( int j = 0; j < 8; ++j ) { 
        _neighborList[i][j] = _neighborListTmp[i][j][0];   
      }
    }
    
    this->updateIterators(); 
    
    _markedForXYRefinement.resize ( this->getNumberOfElements() );
    _markedForZRefinement.resize ( this->getNumberOfElements() );
    
    this->unmarkAll();
      
  }
  
  // Update iterators 
  void updateIterators() {
    begin_it.set ( 0 );
    end_it.set ( this->getNumberOfElements() );
  }
  
  // Mark and refine all elements in current grid
  void refineAll () {
    // Mark all elements for XY and Z refinement 
    int numEls = this->getNumberOfElements(); 
    for ( int elit = 0; elit < numEls; ++elit ) {
      this->markXY ( elit );
      this->markZ ( elit );
    }
    this->refineMarkedElements();
  }
  
  // Mark and refine all elements in current grid in xy direction 
  void refineAllXY () {
    int numEls = this->getNumberOfElements(); 
    for ( int elit = 0; elit < numEls; ++elit ) {
      this->markXY ( elit );
    }
    this->refineMarkedElements();   
  }
  
  // Mark and refine all elements in current grid in z direction
  void refineAllZ () {
    int numEls = this->getNumberOfElements();
    for ( int elit = 0; elit < numEls; ++elit ) {
      this->markZ ( elit );
    }
    this->refineMarkedElements();  
  }

  ElementIterator begin_it;
  ElementIterator end_it;  
    
protected: 

  int _maxlevelz;
  aol::MultiVector < RealType >  _nodes3D;                                         //!< list of all node coordinates 
  std::unordered_map < unsigned long int, GlobalIndex > _nodeMap;                  //!< hashmap with all 3D nodes 
  CellArray3D _elements3D;                                                         //!< list of all elements containing their node indices 
  std::vector < bool > _markedForXYRefinement;
  std::vector < bool > _markedForZRefinement;
  std::map < int, InterpolationInfo < GlobalIndex, RealType > > _interpolationMap;
  std::vector < aol::Vec < 8, GlobalIndex > > _neighborList;                       //!< Neighbor[ind][i] returns global element index of neighbor element along edge i (the last two entries are the element below and above)
  mutable MapType _hangingNodeMap;                                                 //!< Hanging node map 
  std::map < GlobalIndex, std::vector < GlobalIndex > > _subelementList;           //!< List of subelements for each refined element 
  
  // For internal refinement 
  aol::BitVector _visited; 
  std::unordered_map < unsigned long int, GlobalIndex > _addNodes;                 //!< Additional nodes added in refinement routine 
  std::vector < std::array < std::vector < GlobalIndex >, 8 > > _neighborListTmp; 
  GlobalIndex _sizeBeforeRefinement;                                               //!< Number of element before refinement, needed in refinement routine 
  
  
  // Returns local index of node opposite of longest edge (ground local index)
  LocalIndex getLongestEdgeIndex ( GlobalIndex element ) const {
      
    // Get global node indices of triangle
    aol::Vec3 < int > groundNodeInds = this->getGroundNodeIndicesOfElement ( element );
    
    // Assume that edge12 is longest 
    LocalIndex localInd = 0;
    aol::Vec3 < RealType > edge = this->getNodeCoords( groundNodeInds[1] );
    edge -= this->getNodeCoords( groundNodeInds[2] );  
    RealType maxLength = edge.normSqr();
      
    // Check if edge02 is longer 
    edge = this->getNodeCoords ( groundNodeInds[2] );
    edge -= this->getNodeCoords ( groundNodeInds[0] );
    RealType candidate = edge.normSqr();
    if( candidate > maxLength ){
        localInd = 1;
        maxLength = candidate;
    }
    
    // Check if edge01 is longer 
    edge = this->getNodeCoords ( groundNodeInds[0] );
    edge -= this->getNodeCoords ( groundNodeInds[1] );
    candidate = edge.normSqr();
    
    if( candidate > maxLength )
        localInd = 2;
    
    return localInd;
      
  }
  
  // Create new node at edge midpoint 
  aol::Vec3 < RealType > setEdgeMidpoint ( const GlobalIndex nodeIndex1, const GlobalIndex nodeIndex2 ) {
    // Create new coords
    aol::Vec3 < RealType > newNodeCoords = this->getNodeCoords ( nodeIndex1 );
    newNodeCoords += this->getNodeCoords ( nodeIndex2 );
    newNodeCoords /= 2.;
    return newNodeCoords;
  }

  // Insert new element to element list 
  GlobalIndex insertNewElement ( const GlobalIndex elIndex, const aol::Vec < 6, int > & newElement1, const aol::Vec < 6, int > & newElement2 ) {
     _elements3D.push_back ( newElement2 ); 
     _elements3D[elIndex] = newElement1;
     return this->getNumberOfElements() - 1;
  }
  
  // Get lower neighbor of an element along edge 
  std::vector < GlobalIndex > getLowerNeighbor ( const GlobalIndex elIndex, const LocalIndex edge ) const {
    return _neighborListTmp[elIndex][edge];
  }
  
  // Get upper neighbor of an element along edge 
  std::vector < GlobalIndex > getUpperNeighbor ( const GlobalIndex elIndex, const LocalIndex edge ) const {
    return _neighborListTmp[elIndex][edge+3];
  }  
  
  // Check if element contains two nodes 
  bool elementContainsNodes ( const GlobalIndex elementIndex, const aol::Vec2 < GlobalIndex > & sharedEdgeNodeInds ) {
    bool containsBoth = false;
    aol::Vec < 6, GlobalIndex > elementNodeInds = this->getNodeIndicesOfElement ( elementIndex );  
    for ( int i = 0; i < 6; ++i ) {
      if ( elementNodeInds[i] == sharedEdgeNodeInds[0] ) {
        containsBoth = true; 
        break; 
      }
    }
    if ( containsBoth ) {
      containsBoth = false;
      for ( int i = 0; i < 6; ++i ) {
        if ( elementNodeInds[i] == sharedEdgeNodeInds[1] ) {
          containsBoth = true; 
          break; 
        }
      }
    }
    return containsBoth;    
  }
  
  // Check if element contains node 
  bool elementContainsNode ( const GlobalIndex elementIndex, const GlobalIndex nodeIndex ) {
    bool containsNode = false;
    aol::Vec < 6, GlobalIndex > elementNodeInds = this->getNodeIndicesOfElement ( elementIndex );  
    for ( int i = 0; i < 6; ++i ) {
      if ( elementNodeInds[i] == nodeIndex ) {
        containsNode = true; 
        break; 
      }
    }
    return containsNode;    
  }
    
  // Update neighbor list after xy refinement   
  void updateNeighborListXY ( const GlobalIndex elIndex, const GlobalIndex newElementIndex, const LocalIndex refinementEdge ) { 
    
    // Get old neighbor of element elIndex and all edge indices   
    std::array < std::vector < GlobalIndex >, 8 > neighbors = _neighborListTmp [ elIndex ];
    LocalIndex node1Lower = 1; 
    LocalIndex node2Lower = 2; 
    LocalIndex refinementEdgeUpper = 3; 
    LocalIndex node1Upper = 4; 
    LocalIndex node2Upper = 5;    

    //--------------------------------------------------------
    // Update neighbors of elIndex = E1 and newElementIndex = E2 
     
    // Clear old list of elIndex 
    for ( int edge = 0; edge < 8; ++edge ) {
      _neighborListTmp[elIndex][edge].clear();
    }
    
    // Add new entry for newElementIndex 
    _neighborListTmp.resize ( _neighborListTmp.size() + 1 );
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 0 
        
    // E1: Old neighbors at node2Lower 
    for ( unsigned int j = 0; j < neighbors[node2Lower].size(); ++j ) {
      _neighborListTmp[elIndex][0].push_back ( neighbors[node2Lower][j] );
    }    
        
    // E2: Old neighbors at node1Lower 
    for ( unsigned int j = 0; j < neighbors[node1Lower].size(); ++j ) {
      _neighborListTmp[newElementIndex][0].push_back ( neighbors[node1Lower][j] );
      // Update neighbors of neighbor 
      if ( neighbors[node1Lower][j] >= 0 ) {
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[node1Lower][j]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[node1Lower][j]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[node1Lower][j]][k][l] = newElementIndex; 
            }
          }
        }
      }
    }  
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 1
    
    // E1: Old neighbors at refinement edge 
    if ( neighbors[refinementEdge].size() > 1 ) {
       
      GlobalIndex sharedNode1_E1 = this->getElementNodeIndex ( elIndex, 0 ); // new node 
      GlobalIndex sharedNode2_E1 = this->getElementNodeIndex ( elIndex, 2 ); // node1 
      
      bool foundNode1_E1, foundNode2_E1;
    
      // Check E1
      for ( unsigned int k = 0; k < neighbors[refinementEdge].size(); ++k ) {
          
        foundNode1_E1 = false; 
        foundNode2_E1 = false;
        // Check if neighbor k index is >= 0
        if ( neighbors[refinementEdge][k] >= 0 ) {
          // Check if neighbor k contains sharedNode1_E1 and sharedNode2_E1 
          for ( int l = 0; l < 6; ++l ) {
            if ( this->getElementNodeIndex ( neighbors[refinementEdge][k], l ) == sharedNode1_E1 ) {
              foundNode1_E1 = true; 
            }
            if ( this->getElementNodeIndex ( neighbors[refinementEdge][k], l ) == sharedNode2_E1 ) {
              foundNode2_E1 = true;    
            }
            if ( foundNode1_E1 && foundNode2_E1 ) break;
          }
          if ( foundNode1_E1 && foundNode2_E1 ) {
            _neighborListTmp[elIndex][1].push_back ( neighbors[refinementEdge][k] );
            break; 
          }  
        }
      }
      // If no neighbor contains both nodes, neighbor must be -1 
      if ( !foundNode1_E1 && !foundNode2_E1 ) {
        _neighborListTmp[elIndex][1].push_back ( -1 );   
      }
      
    }
    else {
      _neighborListTmp[elIndex][1].push_back ( neighbors[refinementEdge][0] );  
    }
    
    // E2: E1 
    _neighborListTmp[newElementIndex][1].push_back ( elIndex ); 
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 2
    
    // E1: E2 
    _neighborListTmp[elIndex][2].push_back ( newElementIndex ); 
    
    // E2: Old neighbors at refinement edge 
    if ( neighbors[refinementEdge].size() > 1 ) {
        
      GlobalIndex sharedNode1_E2 = this->getElementNodeIndex ( newElementIndex, 0 ); // new node
      GlobalIndex sharedNode2_E2 = this->getElementNodeIndex ( newElementIndex, 1 ); // node2 
      
      bool foundNode1_E2, foundNode2_E2;

      for ( unsigned int k = 0; k < neighbors[refinementEdge].size(); ++k ) {
          
        foundNode1_E2 = false; 
        foundNode2_E2 = false;
        // Check if neighbor k index is >= 0
        if ( neighbors[refinementEdge][k] >= 0 ) {
          // Check if neighbor k contains sharedNode1_E2 and sharedNode2_E2 
          for ( int l = 0; l < 6; ++l ) {
            if ( this->getElementNodeIndex ( neighbors[refinementEdge][k], l ) == sharedNode1_E2 ) { 
              foundNode1_E2 = true; 
            }
            if ( this->getElementNodeIndex ( neighbors[refinementEdge][k], l ) == sharedNode2_E2 ) {
              foundNode2_E2 = true;    
            }
            if ( foundNode1_E2 && foundNode2_E2 ) break;
          }
          if ( foundNode1_E2 && foundNode2_E2 ) {
            _neighborListTmp[newElementIndex][2].push_back ( neighbors[refinementEdge][k] );
            // Update neighbor list of neighbor
            for ( int i = 0; i < 8; ++i ) {
              for ( unsigned int j = 0; j < _neighborListTmp[neighbors[refinementEdge][k]][i].size(); ++j ) {  
                if ( _neighborListTmp[neighbors[refinementEdge][k]][i][j] == elIndex ) { 
                  _neighborListTmp[neighbors[refinementEdge][k]][i][j] = newElementIndex;
                }
              }
            }
            break; 
          }  
        }
      }
      // If no neighbor contains both nodes, neighbor must be -1 
      if ( !foundNode1_E2 && !foundNode2_E2 ) {
        _neighborListTmp[newElementIndex][2].push_back ( -1 );   
      } 
      
    }
    else {
      _neighborListTmp[newElementIndex][2].push_back ( neighbors[refinementEdge][0] );
      
      // Update neighbor list of neighbor 
      if ( neighbors[refinementEdge][0] >= 0 ) {
        for ( int i = 0; i < 3; ++i ) {
          for ( unsigned int j = 0; j < _neighborListTmp[neighbors[refinementEdge][0]][i].size(); ++j ) {  
            if ( _neighborListTmp[neighbors[refinementEdge][0]][i][j] == elIndex ) { 
              _neighborListTmp[neighbors[refinementEdge][0]][i].push_back ( newElementIndex );              
            }
          }
        }
      }
    
    }  
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 3
    
    // E1: Old neighbors at node2Upper 
    for ( unsigned int j = 0; j < neighbors[node2Upper].size(); ++j ) {
      _neighborListTmp[elIndex][3].push_back ( neighbors[node2Upper][j] );
    }    
        
    // E2: Old neighbors at node1Upper 
    for ( unsigned int j = 0; j < neighbors[node1Upper].size(); ++j ) {
      _neighborListTmp[newElementIndex][3].push_back ( neighbors[node1Upper][j] );
      // Update neighbors of neighbor 
      if ( neighbors[node1Upper][j] >= 0 ) {
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[node1Upper][j]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[node1Upper][j]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[node1Upper][j]][k][l] = newElementIndex; 
            }
          }
        }
      }
    } 
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 4
    
    // E1: Old neighbors at refinement edge 
    if ( neighbors[refinementEdgeUpper].size() > 1 ) {
        
      GlobalIndex sharedNode1_E1 = this->getElementNodeIndex ( elIndex, 3 );
      GlobalIndex sharedNode2_E1 = this->getElementNodeIndex ( elIndex, 5 );
      
      bool foundNode1_E1, foundNode2_E1; 
    
      // Check E1
      for ( unsigned int k = 0; k < neighbors[refinementEdgeUpper].size(); ++k ) {
          
        foundNode1_E1 = false; 
        foundNode2_E1 = false;
        // Check if neighbor k index is >= 0
        if ( neighbors[refinementEdgeUpper][k] >= 0 ) {
          // Check if neighbor k contains sharedNode1_E1 and sharedNode2_E1 
          for ( int l = 0; l < 6; ++l ) {
            if ( this->getElementNodeIndex ( neighbors[refinementEdgeUpper][k], l ) == sharedNode1_E1 ) {
              foundNode1_E1 = true; 
            }
            if ( this->getElementNodeIndex ( neighbors[refinementEdgeUpper][k], l ) == sharedNode2_E1 ) {
              foundNode2_E1 = true;    
            }
            if ( foundNode1_E1 && foundNode2_E1 ) break;
          }
          if ( foundNode1_E1 && foundNode2_E1 ) {
            _neighborListTmp[elIndex][4].push_back ( neighbors[refinementEdgeUpper][k] );
            break; 
          }  
        }
      }
      // If no neighbor contains both nodes, neighbor must be -1 
      if ( !foundNode1_E1 && !foundNode2_E1 ) {
        _neighborListTmp[elIndex][4].push_back ( -1 );   
      }
      
    }
    else {        
      _neighborListTmp[elIndex][4].push_back ( neighbors[refinementEdgeUpper][0] );  
    }
    
    // E2: E1 
    _neighborListTmp[newElementIndex][4].push_back ( elIndex ); 
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 5
    
    // E1: E2 
    _neighborListTmp[elIndex][5].push_back ( newElementIndex ); 
    
    // E2: Old neighbors at refinement edge  
    if ( neighbors[refinementEdgeUpper].size() > 1 ) {
        
      GlobalIndex sharedNode1_E2 = this->getElementNodeIndex ( newElementIndex, 3 );
      GlobalIndex sharedNode2_E2 = this->getElementNodeIndex ( newElementIndex, 4 );  

      bool foundNode1_E2, foundNode2_E2;
      
      for ( unsigned int k = 0; k < neighbors[refinementEdgeUpper].size(); ++k ) {
          
        foundNode1_E2 = false; 
        foundNode2_E2 = false;
        // Check if neighbor k index is >= 0
        if ( neighbors[refinementEdgeUpper][k] >= 0 ) {
          // Check if neighbor k contains sharedNode1_E2 and sharedNode2_E2 
          for ( int l = 0; l < 6; ++l ) {
            if ( this->getElementNodeIndex ( neighbors[refinementEdgeUpper][k], l ) == sharedNode1_E2 ) { 
              foundNode1_E2 = true; 
            }
            if ( this->getElementNodeIndex ( neighbors[refinementEdgeUpper][k], l ) == sharedNode2_E2 ) {
              foundNode2_E2 = true;    
            }
            if ( foundNode1_E2 && foundNode2_E2 ) break;
          }
          if ( foundNode1_E2 && foundNode2_E2 ) {
            _neighborListTmp[newElementIndex][5].push_back ( neighbors[refinementEdgeUpper][k] );
            // Update neighbor list of neighbor
            for ( int i = 0; i < 8; ++i ) {
              for ( unsigned int j = 0; j < _neighborListTmp[neighbors[refinementEdgeUpper][k]][i].size(); ++j ) {  
                if ( _neighborListTmp[neighbors[refinementEdgeUpper][k]][i][j] == elIndex ) { 
                  _neighborListTmp[neighbors[refinementEdgeUpper][k]][i][j] = newElementIndex;
                }
              }
            }
            break; 
          }  
        }
      }
      // If no neighbor contains both nodes, neighbor must be -1 
      if ( !foundNode1_E2 && !foundNode2_E2 ) {
        _neighborListTmp[newElementIndex][5].push_back ( -1 );   
      } 
      
    }
    else {
      _neighborListTmp[newElementIndex][5].push_back ( neighbors[refinementEdgeUpper][0] );
      
      // Update neighbor list of neighbor 
      if ( neighbors[refinementEdgeUpper][0] >= 0 ) {
        for ( int i = 3; i < 6; ++i ) {
          for ( unsigned int j = 0; j < _neighborListTmp[neighbors[refinementEdgeUpper][0]][i].size(); ++j ) {  
            if ( _neighborListTmp[neighbors[refinementEdgeUpper][0]][i][j] == elIndex ) { 
              _neighborListTmp[neighbors[refinementEdgeUpper][0]][i].push_back ( newElementIndex );              
            }
          }
        }
      }
    
    }  
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 6

    if ( neighbors[6].size() == 1 ) {
      _neighborListTmp[elIndex][6].push_back ( neighbors[6][0] );  
      _neighborListTmp[newElementIndex][6].push_back ( neighbors[6][0] );  
      
      // Update neighbors of neighbor 
      if ( neighbors[6][0] >= 0 ) {
        for ( unsigned int l = 0; l < _neighborListTmp[neighbors[6][0]][7].size(); ++l ) {
          if ( _neighborListTmp[neighbors[6][0]][7][l] == elIndex ) {
            _neighborListTmp[neighbors[6][0]][7].push_back ( newElementIndex ); 
          }
        }
      }
      
    }
    else if ( neighbors[6].size() == 2 ) {
        
      aol::Vec3 < GlobalIndex > lowerNodesE1 ( this->getElementNodeIndex ( elIndex, 0 ), this->getElementNodeIndex ( elIndex, 1 ), this->getElementNodeIndex ( elIndex, 2 ) );
      aol::Vec3 < GlobalIndex > upperNodesNeighbor1 ( this->getElementNodeIndex ( neighbors[6][0], 3 ), this->getElementNodeIndex ( neighbors[6][0], 4 ), this->getElementNodeIndex ( neighbors[6][0], 5 ) );
      bool shareFace = this->shareFace ( upperNodesNeighbor1, lowerNodesE1 );
      
      if ( shareFace ) {
          
        _neighborListTmp[elIndex][6].push_back ( neighbors[6][0] );
        _neighborListTmp[newElementIndex][6].push_back ( neighbors[6][1] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[6][1]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[6][1]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[6][1]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      else {
        _neighborListTmp[elIndex][6].push_back ( neighbors[6][1] );   
        _neighborListTmp[newElementIndex][6].push_back ( neighbors[6][0] );  
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[6][0]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[6][0]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[6][0]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
        
    }    
    else if ( neighbors[6].size() == 3 ) {

      aol::Vec3 < GlobalIndex > lowerNodesE1 ( this->getElementNodeIndex ( elIndex, 0 ), this->getElementNodeIndex ( elIndex, 1 ), this->getElementNodeIndex ( elIndex, 2 ) );
      aol::Vec3 < GlobalIndex > lowerNodesE2 ( this->getElementNodeIndex ( newElementIndex, 0 ), this->getElementNodeIndex ( newElementIndex, 1 ), this->getElementNodeIndex ( newElementIndex, 2 ) );
      aol::Vec3 < GlobalIndex > upperNodesNeighbor1 ( this->getElementNodeIndex ( neighbors[6][0], 3 ), this->getElementNodeIndex ( neighbors[6][0], 4 ), this->getElementNodeIndex ( neighbors[6][0], 5 ) );   
      aol::Vec3 < GlobalIndex > upperNodesNeighbor2 ( this->getElementNodeIndex ( neighbors[6][1], 3 ), this->getElementNodeIndex ( neighbors[6][1], 4 ), this->getElementNodeIndex ( neighbors[6][1], 5 ) ); 
      aol::Vec3 < GlobalIndex > upperNodesNeighbor3 ( this->getElementNodeIndex ( neighbors[6][2], 3 ), this->getElementNodeIndex ( neighbors[6][2], 4 ), this->getElementNodeIndex ( neighbors[6][2], 5 ) ); 
       
      bool shareFaceE1n1 = this->shareFace ( upperNodesNeighbor1, lowerNodesE1 );
      bool shareFaceE1n2 = this->shareFace ( upperNodesNeighbor2, lowerNodesE1 );
      bool shareFaceE1n3 = this->shareFace ( upperNodesNeighbor3, lowerNodesE1 );
      bool shareFaceE2n1 = this->shareFace ( upperNodesNeighbor1, lowerNodesE2 );
      bool shareFaceE2n2 = this->shareFace ( upperNodesNeighbor2, lowerNodesE2 );
      bool shareFaceE2n3 = this->shareFace ( upperNodesNeighbor3, lowerNodesE2 );
      
      if ( shareFaceE1n1 ) _neighborListTmp[elIndex][6].push_back ( neighbors[6][0] );
      if ( shareFaceE1n2 ) _neighborListTmp[elIndex][6].push_back ( neighbors[6][1] );
      if ( shareFaceE1n3 ) _neighborListTmp[elIndex][6].push_back ( neighbors[6][2] );
      if ( shareFaceE2n1 ) {
        _neighborListTmp[newElementIndex][6].push_back ( neighbors[6][0] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[6][0]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[6][0]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[6][0]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      if ( shareFaceE2n2 ) {
        _neighborListTmp[newElementIndex][6].push_back ( neighbors[6][1] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[6][1]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[6][1]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[6][1]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      if ( shareFaceE2n3 ) {
        _neighborListTmp[newElementIndex][6].push_back ( neighbors[6][2] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[6][2]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[6][2]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[6][2]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      
    }
    else { 
        
      aol::Vec3 < GlobalIndex > lowerNodesE1 ( this->getElementNodeIndex ( elIndex, 0 ), this->getElementNodeIndex ( elIndex, 1 ), this->getElementNodeIndex ( elIndex, 2 ) );
      aol::Vec3 < GlobalIndex > lowerNodesE2 ( this->getElementNodeIndex ( newElementIndex, 0 ), this->getElementNodeIndex ( newElementIndex, 1 ), this->getElementNodeIndex ( newElementIndex, 2 ) );
      aol::Vec3 < GlobalIndex > upperNodesNeighbor1 ( this->getElementNodeIndex ( neighbors[6][0], 3 ), this->getElementNodeIndex ( neighbors[6][0], 4 ), this->getElementNodeIndex ( neighbors[6][0], 5 ) );   
      aol::Vec3 < GlobalIndex > upperNodesNeighbor2 ( this->getElementNodeIndex ( neighbors[6][1], 3 ), this->getElementNodeIndex ( neighbors[6][1], 4 ), this->getElementNodeIndex ( neighbors[6][1], 5 ) ); 
      aol::Vec3 < GlobalIndex > upperNodesNeighbor3 ( this->getElementNodeIndex ( neighbors[6][2], 3 ), this->getElementNodeIndex ( neighbors[6][2], 4 ), this->getElementNodeIndex ( neighbors[6][2], 5 ) ); 
      aol::Vec3 < GlobalIndex > upperNodesNeighbor4 ( this->getElementNodeIndex ( neighbors[6][3], 3 ), this->getElementNodeIndex ( neighbors[6][3], 4 ), this->getElementNodeIndex ( neighbors[6][3], 5 ) );
       
      bool shareFaceE1n1 = this->shareFace ( upperNodesNeighbor1, lowerNodesE1 );
      bool shareFaceE1n2 = this->shareFace ( upperNodesNeighbor2, lowerNodesE1 );
      bool shareFaceE1n3 = this->shareFace ( upperNodesNeighbor3, lowerNodesE1 );
      bool shareFaceE1n4 = this->shareFace ( upperNodesNeighbor4, lowerNodesE1 );
      bool shareFaceE2n1 = this->shareFace ( upperNodesNeighbor1, lowerNodesE2 );
      bool shareFaceE2n2 = this->shareFace ( upperNodesNeighbor2, lowerNodesE2 );
      bool shareFaceE2n3 = this->shareFace ( upperNodesNeighbor3, lowerNodesE2 );
      bool shareFaceE2n4 = this->shareFace ( upperNodesNeighbor4, lowerNodesE2 );
      
      if ( shareFaceE1n1 ) _neighborListTmp[elIndex][6].push_back ( neighbors[6][0] );
      if ( shareFaceE1n2 ) _neighborListTmp[elIndex][6].push_back ( neighbors[6][1] );
      if ( shareFaceE1n3 ) _neighborListTmp[elIndex][6].push_back ( neighbors[6][2] );
      if ( shareFaceE1n4 ) _neighborListTmp[elIndex][6].push_back ( neighbors[6][3] );
      if ( shareFaceE2n1 ) {
        _neighborListTmp[newElementIndex][6].push_back ( neighbors[6][0] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[6][0]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[6][0]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[6][0]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      if ( shareFaceE2n2 ) {
        _neighborListTmp[newElementIndex][6].push_back ( neighbors[6][1] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[6][1]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[6][1]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[6][1]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      if ( shareFaceE2n3 ) {
        _neighborListTmp[newElementIndex][6].push_back ( neighbors[6][2] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[6][2]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[6][2]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[6][2]][k][l] = newElementIndex; 
            }
          }
        }
        
      }  
      if ( shareFaceE2n4 ) {
        _neighborListTmp[newElementIndex][6].push_back ( neighbors[6][3] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[6][3]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[6][3]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[6][3]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
        
    }
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 7
    if ( neighbors[7].size() == 1 ) {
        
      _neighborListTmp[elIndex][7].push_back ( neighbors[7][0] );  
      _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][0] ); 
      
      // Update neighbors of neighbor 
      if ( neighbors[7][0] >= 0 ) {
        for ( unsigned int l = 0; l < _neighborListTmp[neighbors[7][0]][6].size(); ++l ) {
          if ( _neighborListTmp[neighbors[7][0]][6][l] == elIndex ) {
            _neighborListTmp[neighbors[7][0]][6].push_back ( newElementIndex ); 
          }
        }
      }
      
    }
    else if ( neighbors[7].size() == 2 ) {
            
      aol::Vec3 < GlobalIndex > upperNodesE1 ( this->getElementNodeIndex ( elIndex, 3 ), this->getElementNodeIndex ( elIndex, 4 ), this->getElementNodeIndex ( elIndex, 5 ) );
      aol::Vec3 < GlobalIndex > lowerNodesNeighbor1 ( this->getElementNodeIndex ( neighbors[7][0], 0 ), this->getElementNodeIndex ( neighbors[7][0], 1 ), this->getElementNodeIndex ( neighbors[7][0], 2 ) );
      bool shareFace = this->shareFace ( lowerNodesNeighbor1, upperNodesE1 );
      if ( shareFace ) {
        _neighborListTmp[elIndex][7].push_back ( neighbors[7][0] );
        _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][1] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[7][1]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[7][1]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[7][1]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      else {
        _neighborListTmp[elIndex][7].push_back ( neighbors[7][1] );   
        _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][0] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[7][0]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[7][0]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[7][0]][k][l] = newElementIndex; 
            }
          }
        }        
        
      }
        
    }    
    else if ( neighbors[7].size() == 3 ) {
        
      aol::Vec3 < GlobalIndex > upperNodesE1 ( this->getElementNodeIndex ( elIndex, 3 ), this->getElementNodeIndex ( elIndex, 4 ), this->getElementNodeIndex ( elIndex, 5 ) );
      aol::Vec3 < GlobalIndex > upperNodesE2 ( this->getElementNodeIndex ( newElementIndex, 3 ), this->getElementNodeIndex ( newElementIndex, 4 ), this->getElementNodeIndex ( newElementIndex, 5 ) );
      aol::Vec3 < GlobalIndex > lowerNodesNeighbor1 ( this->getElementNodeIndex ( neighbors[7][0], 0 ), this->getElementNodeIndex ( neighbors[7][0], 1 ), this->getElementNodeIndex ( neighbors[7][0], 2 ) );   
      aol::Vec3 < GlobalIndex > lowerNodesNeighbor2 ( this->getElementNodeIndex ( neighbors[7][1], 0 ), this->getElementNodeIndex ( neighbors[7][1], 1 ), this->getElementNodeIndex ( neighbors[7][1], 2 ) ); 
      aol::Vec3 < GlobalIndex > lowerNodesNeighbor3 ( this->getElementNodeIndex ( neighbors[7][2], 0 ), this->getElementNodeIndex ( neighbors[7][2], 1 ), this->getElementNodeIndex ( neighbors[7][2], 2 ) ); 
       
      bool shareFaceE1n1 = this->shareFace ( lowerNodesNeighbor1, upperNodesE1 );
      bool shareFaceE1n2 = this->shareFace ( lowerNodesNeighbor2, upperNodesE1 );
      bool shareFaceE1n3 = this->shareFace ( lowerNodesNeighbor3, upperNodesE1 );
      bool shareFaceE2n1 = this->shareFace ( lowerNodesNeighbor1, upperNodesE2 );
      bool shareFaceE2n2 = this->shareFace ( lowerNodesNeighbor2, upperNodesE2 );
      bool shareFaceE2n3 = this->shareFace ( lowerNodesNeighbor3, upperNodesE2 );
      
      if ( shareFaceE1n1 ) _neighborListTmp[elIndex][7].push_back ( neighbors[7][0] );
      if ( shareFaceE1n2 ) _neighborListTmp[elIndex][7].push_back ( neighbors[7][1] );
      if ( shareFaceE1n3 ) _neighborListTmp[elIndex][7].push_back ( neighbors[7][2] );
      if ( shareFaceE2n1 ) {
        _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][0] ); 
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[7][0]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[7][0]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[7][0]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      if ( shareFaceE2n2 ) { 
        _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][1] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[7][1]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[7][1]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[7][1]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      if ( shareFaceE2n3 ) {
        _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][2] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[7][2]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[7][2]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[7][2]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      
    }
    else {  
        
      aol::Vec3 < GlobalIndex > upperNodesE1 ( this->getElementNodeIndex ( elIndex, 3 ), this->getElementNodeIndex ( elIndex, 4 ), this->getElementNodeIndex ( elIndex, 5 ) );
      aol::Vec3 < GlobalIndex > upperNodesE2 ( this->getElementNodeIndex ( newElementIndex, 3 ), this->getElementNodeIndex ( newElementIndex, 4 ), this->getElementNodeIndex ( newElementIndex, 5 ) );
      aol::Vec3 < GlobalIndex > lowerNodesNeighbor1 ( this->getElementNodeIndex ( neighbors[7][0], 0 ), this->getElementNodeIndex ( neighbors[7][0], 1 ), this->getElementNodeIndex ( neighbors[7][0], 2 ) );   
      aol::Vec3 < GlobalIndex > lowerNodesNeighbor2 ( this->getElementNodeIndex ( neighbors[7][1], 0 ), this->getElementNodeIndex ( neighbors[7][1], 1 ), this->getElementNodeIndex ( neighbors[7][1], 2 ) ); 
      aol::Vec3 < GlobalIndex > lowerNodesNeighbor3 ( this->getElementNodeIndex ( neighbors[7][2], 0 ), this->getElementNodeIndex ( neighbors[7][2], 1 ), this->getElementNodeIndex ( neighbors[7][2], 2 ) ); 
      aol::Vec3 < GlobalIndex > lowerNodesNeighbor4 ( this->getElementNodeIndex ( neighbors[7][3], 0 ), this->getElementNodeIndex ( neighbors[7][3], 1 ), this->getElementNodeIndex ( neighbors[7][3], 2 ) ); 
      
      bool shareFaceE1n1 = this->shareFace ( lowerNodesNeighbor1, upperNodesE1 );
      bool shareFaceE1n2 = this->shareFace ( lowerNodesNeighbor2, upperNodesE1 );
      bool shareFaceE1n3 = this->shareFace ( lowerNodesNeighbor3, upperNodesE1 );
      bool shareFaceE1n4 = this->shareFace ( lowerNodesNeighbor4, upperNodesE1 );
      bool shareFaceE2n1 = this->shareFace ( lowerNodesNeighbor1, upperNodesE2 );
      bool shareFaceE2n2 = this->shareFace ( lowerNodesNeighbor2, upperNodesE2 );
      bool shareFaceE2n3 = this->shareFace ( lowerNodesNeighbor3, upperNodesE2 );
      bool shareFaceE2n4 = this->shareFace ( lowerNodesNeighbor4, upperNodesE2 );
      
      if ( shareFaceE1n1 ) _neighborListTmp[elIndex][7].push_back ( neighbors[7][0] );
      if ( shareFaceE1n2 ) _neighborListTmp[elIndex][7].push_back ( neighbors[7][1] );
      if ( shareFaceE1n3 ) _neighborListTmp[elIndex][7].push_back ( neighbors[7][2] );
      if ( shareFaceE1n4 ) _neighborListTmp[elIndex][7].push_back ( neighbors[7][3] ); 
      if ( shareFaceE2n1 ) {
        _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][0] ); 
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[7][0]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[7][0]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[7][0]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      if ( shareFaceE2n2 ) { 
        _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][1] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[7][1]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[7][1]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[7][1]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      if ( shareFaceE2n3 ) {
        _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][2] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[7][2]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[7][2]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[7][2]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      if ( shareFaceE2n4 ) {
        _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][3] );
        
        // Update neighbors of neighbor 
        for ( unsigned int k = 0; k < 8; ++k ) {
          for ( unsigned int l = 0; l < _neighborListTmp[neighbors[7][3]][k].size(); ++l ) {
            if ( _neighborListTmp[neighbors[7][3]][k][l] == elIndex ) {
              _neighborListTmp[neighbors[7][3]][k][l] = newElementIndex; 
            }
          }
        }
        
      }
      
    }
    
    
    //--------------------------------------------------------
    
  }
  
  // Update neighbor list after z refinement 
  void updateNeighborListZ ( const GlobalIndex elIndex, const GlobalIndex newElementIndex ) { 
        
    // Get old neighbors of element elIndex 
    std::array < std::vector < GlobalIndex >, 8 > neighbors = _neighborListTmp [ elIndex ];
  
    // Add new entry for newElementIndex
    _neighborListTmp.resize ( _neighborListTmp.size() + 1 );
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Edge neighbors 
    for ( int k = 0; k < 3; ++k ) {
      // Lower and upper neighbor are the same = neighbor is not refined  
      if ( neighbors[k][0] == neighbors[k+3][0] ) {
        if ( neighbors[k][0] == -2 ) {
          _neighborListTmp[elIndex][k+3][0] = -2;   
          _neighborListTmp[newElementIndex][k].push_back ( -2 );
        }
        else {
          _neighborListTmp[elIndex][k+3][0] = -1; 
          _neighborListTmp[newElementIndex][k].push_back ( -1 );
        }
        _neighborListTmp[newElementIndex][k+3].push_back ( neighbors[k][0] ); 
        // Update neighbors of neighbor 
        if ( neighbors[k][0] >= 0 ) {
          for ( int j = 3; j < 6; ++j ) {  
            if ( _neighborListTmp[neighbors[k][0]][j][0] == elIndex ) {
              _neighborListTmp[neighbors[k][0]][j][0] = newElementIndex;    
            }
          }
        }
      }
      // Lower and upper neighbor are different = neighbor is refined 
      else {
          
        _neighborListTmp[elIndex][k+3][0] = neighbors[k][0];
        _neighborListTmp[newElementIndex][k].push_back ( neighbors[k+3][0] ); 
        _neighborListTmp[newElementIndex][k+3].push_back ( neighbors[k+3][0] );
        // Update neighbors of neighbor 
        if ( neighbors[k][0] >= 0 ) {
          for ( int j = 0; j < 3; ++j ) {
            if ( _neighborListTmp[neighbors[k][0]][j][0] == elIndex ) {
              _neighborListTmp[neighbors[k][0]][j+3][0] = elIndex;   
            }
          }
        }
        if ( neighbors[k+3][0] >= 0 ) {
            
          for ( int j = 3; j < 6; ++j ) {
              
            if ( _neighborListTmp[neighbors[k+3][0]][j][0] == elIndex ) {
              _neighborListTmp[neighbors[k+3][0]][j][0] = newElementIndex;   
              _neighborListTmp[neighbors[k+3][0]][j-3][0] = newElementIndex;    
            }
          }
        }
      }    
    }
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 6 
    _neighborListTmp[newElementIndex][6].push_back ( elIndex ); 
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 7
    _neighborListTmp[elIndex][7][0] = newElementIndex; 
    _neighborListTmp[newElementIndex][7].push_back ( neighbors[7][0] );
    if ( neighbors[7][0] >= 0 ) {
      _neighborListTmp[neighbors[7][0]][6][0] = newElementIndex;   
    }
      
  }
  
  // Check if two elements share a face 
  bool shareFace ( const aol::Vec3 < GlobalIndex > & nodes1, const aol::Vec3 < GlobalIndex > & nodes2 ) { 
    
    aol::Vec3 < RealType > P1 = this->getNodeCoords ( nodes2[0] );
    aol::Vec3 < RealType > P2 = this->getNodeCoords ( nodes2[1] );
    aol::Vec3 < RealType > P3 = this->getNodeCoords ( nodes2[2] );
      
    // Get convex combination of nodes2 
    aol::Matrix22 < RealType > A ( P1[0]-P3[0], P2[0]-P3[0], P1[1]-P3[1], P2[1]-P3[1] );
    A.invert();
    aol::Vec2 < RealType > sol; 
    
    // Check all points of nodes1 
    bool check = true;
    for ( int i = 0; i < 3; ++i ) {
      aol::Vec2 < RealType > b ( this->getNodeCoords(nodes1[i])[0]-P3[0], this->getNodeCoords(nodes1[i])[1]-P3[1] );
      A.apply ( b, sol );
      
      if ( sol[0] < 0.0 || sol[0] > 1.0 || sol[1] < 0.0 || sol[1] > 1.0 || sol[0]+sol[1] < 0.0 || sol[0]+sol[1] > 1.0 ) { 
        check = false; 
        break;
      }
    }
      
    return check;  
  }
  
  // Get lower and upper neighbor of an element 
  std::array < std::vector < GlobalIndex >, 2 > getElementBelowAndAbove ( const GlobalIndex elIndex ) const { 
    return std::array < std::vector < GlobalIndex >, 2 > ( {_neighborListTmp[elIndex][6], _neighborListTmp[elIndex][7]} ); 
  }  
  
  // Get element below 
  std::vector < GlobalIndex > getElementBelow ( const GlobalIndex elIndex ) const {
    return _neighborListTmp[elIndex][6];   
  }
  
  // Get element above 
  std::vector < GlobalIndex > getElementAbove ( const GlobalIndex elIndex ) const {
    return _neighborListTmp[elIndex][7];   
  }

  // Update interpolation map after xy refinement 
  void updateInterpolationMap ( const GlobalIndex elIndex, const GlobalIndex newNodeIndex, const LocalIndex interpNode1, const LocalIndex interpNode2 ) {
      
    InterpolationInfo < GlobalIndex, RealType > interpolInfo; 
    // Check if interpolation nodes are new nodes as well 
    std::stack < std::pair < GlobalIndex, RealType > > interpolStack; 
    interpolStack.push ( std::pair < GlobalIndex, RealType > ( this->getElementNodeIndex ( elIndex, interpNode1 ), 0.5 ) );
    interpolStack.push ( std::pair < GlobalIndex, RealType > ( this->getElementNodeIndex ( elIndex, interpNode2 ), 0.5 ) );
    int counter = 0; 
    while ( !interpolStack.empty() ) {
      std::pair < GlobalIndex, RealType > current = interpolStack.top();
      interpolStack.pop();
      auto gotInterpolNode = _addNodes.find ( this->convertCoord3DToHash ( this->getNodeCoords ( current.first ) ) );  
      // Node not found: Node existed before, push node to interpolInfo 
      if ( gotInterpolNode == _addNodes.end() ) {
        //if ( newNodeIndex == 3149 ) cout << "Node existed before" << endl; 
        interpolInfo._nodes.push_back ( current.first );
        interpolInfo._weights.push_back ( current.second );
      }          
      // Node found: Node did not exist before, insert interpolation nodes of this node 
      else {
        //if ( newNodeIndex == 3149 ) cout << "Node did not exist before" << endl; 
        // Find current.first in interpolation map and push interpolation nodes to interpolStack   
        auto interpol = _interpolationMap.find ( current.first ); 
        for ( unsigned int j = 0; j < interpol->second._nodes.size(); ++j ) {
          //if ( newNodeIndex == 3149 ) cout << interpol->second._nodes[j] << ", coords " << this->getNodeCoords ( interpol->second._nodes[j] ) << endl;  
          interpolStack.push ( std::pair < GlobalIndex, RealType > ( interpol->second._nodes[j], current.second*interpol->second._weights[j] ) );
        }
      }
      ++counter;
      if ( counter > 100 ) {
        cout << "EMERGENCY BREAK!!!" << endl;  
        throw aol::Exception ( "AdaptiveFEPrismMesh::updateInterpolationMap: Interpolation stack is filling up!", __FILE__, __LINE__ );
      }
    }   
    _interpolationMap.insert ( std::pair < int, InterpolationInfo < GlobalIndex, RealType > > ( newNodeIndex, interpolInfo ) );

  }
  
  // Update interpolation map after z refinement 
  void updateInterpolationMapZ ( const GlobalIndex nodeInd, const GlobalIndex newNodeInd ) {                          
        
    InterpolationInfo < GlobalIndex, RealType > interpolInfo; 
    // Check if interpolation nodes are new nodes as well 
    std::stack < std::pair < GlobalIndex, RealType > > interpolStack; 
    interpolStack.push ( std::pair < GlobalIndex, RealType > ( nodeInd, 1.0 ) ); // ATTENTION: For linear interpolation, use lower and upper node with weight 1/2 instead 
    int counter = 0; 
    while ( !interpolStack.empty() ) {
      std::pair < GlobalIndex, RealType > current = interpolStack.top();
      interpolStack.pop();
      auto gotInterpolNode = _addNodes.find ( this->convertCoord3DToHash ( this->getNodeCoords ( current.first ) ) );  
      // Node not found: Node existed before, push node to interpolInfo 
      if ( gotInterpolNode == _addNodes.end() ) {
        //if ( nodeInd == 3149 ) cout << "Node existed before" << endl;   
        interpolInfo._nodes.push_back ( current.first );
        interpolInfo._weights.push_back ( current.second );
      }          
      // Node found: Node did not exist before, insert interpolation nodes of this node 
      else {
        //if ( nodeInd == 3149 ) cout << "Node " << current.first << " is additional node as well, check interpolation nodes" << endl; 
        // Find current.first in interpolation map and push interpolation nodes to interpolStack   
        auto interpol = _interpolationMap.find ( current.first ); 
        for ( unsigned int j = 0; j < interpol->second._nodes.size(); ++j ) {
          //if ( nodeInd == 3149 ) cout << this->getNodeCoords ( interpol->second._nodes[j] ) << endl;   
          interpolStack.push ( std::pair < GlobalIndex, RealType > ( interpol->second._nodes[j], current.second*interpol->second._weights[j] ) );
        }
      }
      ++counter; 
      if ( counter > 100 ) {
        cout << "EMERGENCY BREAK!!!" << endl;  
        throw aol::Exception ( "AdaptiveFEPrismMesh::updateInterpolationMapZ: Interpolation stack is filling up!", __FILE__, __LINE__ );
      }
    }   
    _interpolationMap.insert ( std::pair < int, InterpolationInfo < GlobalIndex, RealType > > ( newNodeInd, interpolInfo ) ); 

  }
  
  // Update subelement list 
  void updateSubelementList ( const GlobalIndex elIndex, const GlobalIndex newElementIndex ) {
  
    auto gotCurrentEl = _subelementList.find ( elIndex );
    // If element is already in list, check if element appeared in original grid 
    if ( gotCurrentEl != _subelementList.end() ) {
      if ( elIndex >= _sizeBeforeRefinement ) { 
        auto gotOriginalEl = _subelementList.find ( gotCurrentEl->second[0] );
        _subelementList[gotOriginalEl->first].push_back ( newElementIndex );
      }
      else {
        _subelementList[gotCurrentEl->first].push_back ( newElementIndex );  
      }
    }
    // If element is not in list, create new entry for element and new element 
    else {
      std::vector < GlobalIndex > newSubelement; newSubelement.push_back ( elIndex ); newSubelement.push_back ( newElementIndex );
      std::vector < GlobalIndex > newSubelement2; newSubelement2.push_back ( elIndex );
      _subelementList.insert ( std::pair < GlobalIndex, std::vector < GlobalIndex > > ( elIndex, newSubelement ) );
      _subelementList.insert ( std::pair < GlobalIndex, std::vector < GlobalIndex > > ( newElementIndex, newSubelement2 ) );
    }
    
  }
  
  // Refine element in xy direction 
  GlobalIndex refineXY ( const GlobalIndex elIndex ) {
      
      aol::Vec < 6, GlobalIndex > nodeInds =  this->getNodeIndicesOfElement ( elIndex );      

      // Find longest edge of element 
      LocalIndex longestEdgeIndex = this->getLongestEdgeIndex ( elIndex );
      
      // Set coordinates of new points 
      LocalIndex localIndTmp = (longestEdgeIndex+1)%3;
      aol::Vec3 < RealType > newNodeLower = this->setEdgeMidpoint ( nodeInds[localIndTmp], nodeInds[3 - (localIndTmp+longestEdgeIndex)] );
      aol::Vec3 < RealType > newNodeUpper = this->setEdgeMidpoint ( nodeInds[localIndTmp+3], nodeInds[6 - (localIndTmp+longestEdgeIndex)] );
      
      // Get new node index for lower node (first check if node is already in addNodes) 
      GlobalIndex newNodeIndexLower;
      std::unordered_map < long unsigned int, int >::const_iterator got = _addNodes.find ( this->convertCoord3DToHash ( newNodeLower ) );
      if ( got == _addNodes.end() ) {
        newNodeIndexLower = this->getNumberOfNodes();
        // Add node to _nodes3D, _nodeMap and _addNodes 
        for ( short i = 0; i < 3; ++i ) {
          _nodes3D[i].pushBack ( newNodeLower[i] );
        }
        std::pair < long unsigned int, int > newHash ( this->convertCoord3DToHash ( newNodeLower ), newNodeIndexLower );
        _nodeMap.insert ( newHash );
        _addNodes.insert ( newHash ); 
        
        // Update interpolation map 
        this->updateInterpolationMap ( elIndex, newNodeIndexLower, localIndTmp, 3 - (localIndTmp+longestEdgeIndex) );        
        
      }
      else {
        newNodeIndexLower = got->second;   
      }
      
      // Get new node index for upper node (first check if node is already in addNodes) 
      GlobalIndex newNodeIndexUpper;
      got = _addNodes.find ( this->convertCoord3DToHash ( newNodeUpper ) );
      if ( got == _addNodes.end() ) {
        newNodeIndexUpper = this->getNumberOfNodes();
        // Add node to _nodes3D, _nodeMap and _addNodes 
        for ( short i = 0; i < 3; ++i ) {
          _nodes3D[i].pushBack ( newNodeUpper[i] );
        }
        std::pair < long unsigned int, int > newHash ( this->convertCoord3DToHash ( newNodeUpper ), newNodeIndexUpper );
        _nodeMap.insert ( newHash );
        _addNodes.insert ( newHash ); 
        
        // Update interpolation map 
        this->updateInterpolationMap ( elIndex, newNodeIndexUpper, localIndTmp+3, 6 - (localIndTmp+longestEdgeIndex) );        
        
      }
      else {
        newNodeIndexUpper = got->second;   
      }
      
      // Create new elements 
      aol::Vec < 6, GlobalIndex > newElement1; 
      newElement1[0] = newNodeIndexLower; newElement1[1] = this->getElementNodeIndex ( elIndex, longestEdgeIndex ); newElement1[2] = this->getElementNodeIndex ( elIndex, localIndTmp );
      newElement1[3] = newNodeIndexUpper; newElement1[4] = this->getElementNodeIndex ( elIndex, longestEdgeIndex + 3 ); newElement1[5] = this->getElementNodeIndex ( elIndex, localIndTmp + 3 );
      aol::Vec < 6, GlobalIndex > newElement2;
      newElement2[0] = newNodeIndexLower; newElement2[1] = this->getElementNodeIndex ( elIndex, 3 - (localIndTmp+longestEdgeIndex) ); newElement2[2] = this->getElementNodeIndex ( elIndex, longestEdgeIndex );
      newElement2[3] = newNodeIndexUpper; newElement2[4] = this->getElementNodeIndex ( elIndex, 6 - (localIndTmp+longestEdgeIndex) ); newElement2[5] = this->getElementNodeIndex ( elIndex, longestEdgeIndex + 3 );
      
      // Change current element to newElement1 and push newElement2 to elements3D 
      GlobalIndex newElementIndex = this->insertNewElement ( elIndex, newElement1, newElement2 ); // ATTENTION: Changes _elements3D[elIndex] to first new element and pushes the second new element to the list 
      _visited.set ( elIndex, true ); 
      _visited.resize ( _visited.size() + 1 );
      _visited.set ( newElementIndex, true ); 
      
      // Update subelementList 
      this->updateSubelementList ( elIndex, newElementIndex );
      
      //_markedForXYRefinement[elIndex] = true; 
      _markedForXYRefinement.push_back ( false );
      
      // Mark both elements for z refinement if old element was marked 
      if ( _markedForZRefinement[elIndex] ) _markedForZRefinement.push_back ( true );   
      else _markedForZRefinement.push_back ( false );   
      
      // Find edge-sharing neighbors
      std::vector < GlobalIndex > neighbor1 = this->getLowerNeighbor ( elIndex, longestEdgeIndex ); 
      std::vector < GlobalIndex > neighbor2 = this->getUpperNeighbor ( elIndex, longestEdgeIndex );
      std::vector < GlobalIndex > neighborBelow = this->getElementBelow ( elIndex ); 
      std::vector < GlobalIndex > neighborAbove = this->getElementAbove ( elIndex ); 
      
      // Update hanging node map       
      this->updateHangingNodeMap ( nodeInds, newNodeLower, newNodeUpper, neighbor1, neighbor2, neighborBelow, neighborAbove, newNodeIndexLower, newNodeIndexUpper, elIndex, newElementIndex, longestEdgeIndex );
      
      // Update _neighborList 
      this->updateNeighborListXY ( elIndex, newElementIndex, longestEdgeIndex );
    
      // Check if neighbor is in domain 
      if ( neighbor1[0] != -2 ) {
         
        GlobalIndex neighborToRefine = -1; 
        
        //--------------------------------------------------------------------------------------------------- 
        // Check lower neighbor 
        if ( neighbor1.size() == 1 ) {
          if ( neighbor1[0] >= 0 ) {
            neighborToRefine = neighbor1[0]; 
            if ( !_visited[neighborToRefine] ) {
                this->refineNeighbor ( neighborToRefine, longestEdgeIndex, nodeInds ); 
            }
            else {
                
              // If neighbor was already visited, check if this neighbor has to be refined another time
              std::vector < GlobalIndex > currentNeighborsE1 = _neighborListTmp[elIndex][1];  
              std::vector < GlobalIndex > currentNeighborsE2 = _neighborListTmp[newElementIndex][2]; 
              if ( currentNeighborsE1[0] == currentNeighborsE2[0] ) {
                  refineXY ( currentNeighborsE1[0] ); 
              }
            }
          }
        }
        //--------------------------------------------------------------------------------------------------- 
        // Check upper neighbor 
        if ( neighbor2.size() == 1 ) {  
          if ( neighbor2[0] >= 0 ) {
            // If lower and upper neighbor coincide, neighbor2 was already refined 
            if ( neighbor1[0] != neighbor2[0] ) {
              neighborToRefine = neighbor2[0];
              
              if ( !_visited[neighborToRefine] ) {
                this->refineNeighbor ( neighborToRefine, longestEdgeIndex, nodeInds );     
              }
              else {  
                // If neighbor was already visited, check if this neighbor has to be refined another time
                std::vector < GlobalIndex > currentNeighborsE1 = _neighborListTmp[elIndex][4];  
                std::vector < GlobalIndex > currentNeighborsE2 = _neighborListTmp[newElementIndex][5]; 
                if ( currentNeighborsE1[0] == currentNeighborsE2[0] ) {
                  refineXY ( currentNeighborsE1[0] ); 
                }
              }
            }
          }
        }
        //---------------------------------------------------------------------------------------------------   
          
      }
      
      return newElementIndex;
  }
    
  // Refine element in z direction   
  GlobalIndex refineZ ( const GlobalIndex elIndex ) {
      
    aol::Vec < 6, GlobalIndex > nodeInds =  this->getNodeIndicesOfElement ( elIndex );  
    
    // Set coordinates and indices of new points (first check if nodes are already in nodeMap)
    std::vector < aol::Vec3 < RealType > > newNodeCoords; 
    std::vector < GlobalIndex > newNodeInds;
    for ( int i = 0; i < 3; ++i ) {        
      newNodeCoords.push_back ( this->setEdgeMidpoint ( nodeInds[i], nodeInds[i+3] ) );
      std::unordered_map < long unsigned int, int >::const_iterator got = _nodeMap.find ( this->convertCoord3DToHash ( newNodeCoords[i] ) );
      // If node is new, add it to interpolation map and hangingNodeMap
      if ( got == _nodeMap.end() ) {
        newNodeInds.push_back ( this->getNumberOfNodes() );   
        // Add node to _nodes3D and _nodeMap 
        for ( int j = 0; j < 3; ++j ) {
          _nodes3D[j].pushBack ( newNodeCoords[i][j] );    
        }
        std::pair < long unsigned int, int > newHash ( this->convertCoord3DToHash ( newNodeCoords[i] ), newNodeInds[i] );
        _nodeMap.insert ( newHash ); 
        _addNodes.insert ( newHash ); 
        
        // Update interpolation map 
        this->updateInterpolationMapZ ( nodeInds[i], newNodeInds[i] ); 
        
      }
      // If node already exists, get its index 
      else {
        newNodeInds.push_back ( got->second );  
      }
      
    }
    
    // Go through all new nodes 
    for ( int i = 0; i < 3; ++i ) {
        
      // Check for 2-1-rule   
      // Check if one interpolation node is hanging: If yes, recursively refine all neighbors of this hanging node  
      auto gotInterpolationNode1 = _hangingNodeMap.find ( nodeInds[i] );
      auto gotInterpolationNode2 = _hangingNodeMap.find ( nodeInds[i+3] );
      if ( gotInterpolationNode1 != _hangingNodeMap.end() ) {
        //cout << "Node 1 is hanging " << endl;  
        std::unordered_set < GlobalIndex > node1Neighbors = gotInterpolationNode1->second._neighbors;
        for ( auto j = node1Neighbors.begin(); j != node1Neighbors.end(); ++j ) {
          if ( !_visited[*j] ) refineZ ( *j );
        }
      }
      else if ( gotInterpolationNode2 != _hangingNodeMap.end() ) {
          //cout << "Node 2 is hanging " << endl; 
          std::unordered_set < GlobalIndex > node2Neighbors = gotInterpolationNode2->second._neighbors;
          for ( auto j = node2Neighbors.begin(); j != node2Neighbors.end(); ++j ) {
            if ( !_visited[*j] ) refineZ ( *j );  
          }
      }
    }
    
    // Create new elements 
    aol::Vec < 6, GlobalIndex > newElement1, newElement2;
    newElement1[0] = nodeInds[0]; newElement1[1] = nodeInds[1]; newElement1[2] = nodeInds[2]; newElement1[3] = newNodeInds[0]; newElement1[4] = newNodeInds[1]; newElement1[5] = newNodeInds[2]; 
    newElement2[0] = newNodeInds[0]; newElement2[1] = newNodeInds[1]; newElement2[2] = newNodeInds[2]; newElement2[3] = nodeInds[3]; newElement2[4] = nodeInds[4]; newElement2[5] = nodeInds[5]; 
    
    // Change current element to newElement1 and push newElement2 to elements3D 
    GlobalIndex newElementIndex = this->insertNewElement ( elIndex, newElement1, newElement2 ); // ATTENTION: Changes _elements3D[elIndex] to first new element and pushes the second new element to the list 
    _visited.set ( elIndex, true ); 
    _markedForZRefinement.push_back ( false );
    
    // Update subelementList 
    this->updateSubelementList ( elIndex, newElementIndex );
    
    // Go through all new nodes 
    for ( int i = 0; i < 3; ++i ) {
        
      // Update hanging node map 
      auto gotNewNode = _hangingNodeMap.find ( newNodeInds[i] );
      // If new node already is in hanging node map, remove this element from neighbor list 
      if ( gotNewNode != _hangingNodeMap.end() ) {
          
        auto gotNeighbor = gotNewNode->second._neighbors.find ( elIndex );
        gotNewNode->second._neighbors.erase ( gotNeighbor );
        
        // If list is empty now, this node is a dof now 
        if ( gotNewNode->second._neighbors.empty() ) {
          _hangingNodeMap.erase ( gotNewNode );
        }
      }
      // If new node is not in _hangingNodeMap, create new entry 
      else {
       
        std::unordered_set < GlobalIndex > newNeighbors = this->findElementsContainingNodes ( nodeInds[i], nodeInds[i+3] );  
        
        if ( !newNeighbors.empty() ) {
          HangingNodeInfo < GlobalIndex > newHangingNodeInfo; 
          std::pair < GlobalIndex, GlobalIndex > newInterpolationNodes ( nodeInds[i], nodeInds[i+3] );
          newHangingNodeInfo._interpolationNodes = newInterpolationNodes; 
          newHangingNodeInfo._neighbors = newNeighbors;
          _hangingNodeMap [ newNodeInds[i] ] = newHangingNodeInfo;
        }
        
      }
      
    }
    
    // Update neighborList 
    this->updateNeighborListZ ( elIndex, newElementIndex ); 
    
    // Update _maxlevelz
    RealType elHeight = this->getNodeCoord ( nodeInds[3], 2 ) - this->getNodeCoord ( nodeInds[0], 2 ); 
    int elLevel = std::log2 ( 1.0/elHeight + 0.5 );
    _maxlevelz = aol::Max ( _maxlevelz, elLevel + 1 );
    
    return newElementIndex;  
  }
  
  // Update hanging node map 
  void updateHangingNodeMap ( const aol::Vec < 6, GlobalIndex > & nodeInds, const aol::Vec3 < RealType > & newNodeLower, const aol::Vec3 < RealType > & newNodeUpper, const std::vector < GlobalIndex > & neighbor1, 
                              const std::vector < GlobalIndex > & neighbor2, const std::vector < GlobalIndex > & neighborBelow, const std::vector < GlobalIndex > & neighborAbove, const GlobalIndex newNodeIndexLower, 
                              const GlobalIndex newNodeIndexUpper, const GlobalIndex elIndex, const GlobalIndex newElementIndex, const LocalIndex longestEdgeIndex ) { 
      
    LocalIndex localIndTmp = (longestEdgeIndex+1)%3; 
      
    // If lower or upper edge neighbor are -1, neighbor z level is higher than of this element -> new hanging node will appear! 
    if ( neighbor1[0] == -1 ) {
        
      // Check if there are one or two elements below 
      if ( neighborBelow.size() > 1 ) {
        // New hanging node appears -> lower new node 
        HangingNodeInfo < GlobalIndex > newHangingNodeInfo;  
        aol::Vec3 < RealType > auxNodeCoords = this->getNodeCoords ( newNodeIndexLower );
        auxNodeCoords[2] = this->getNodeCoord ( this->getElementNodeIndex ( neighborBelow[0], 0 ), 2 );
        // Find index of auxNodeCoords (node must exist) 
        auto auxNode = _nodeMap.find ( this->convertCoord3DToHash ( auxNodeCoords ) );
        newHangingNodeInfo._interpolationNodes = std::pair < GlobalIndex, GlobalIndex > ( auxNode->second, newNodeIndexUpper );
        // Set neighbors of new hanging node = neighbor2 
        for ( unsigned int j = 0; j < neighbor2.size(); ++j ) {
          newHangingNodeInfo._neighbors.insert ( neighbor2[j] ); 
        }
        _hangingNodeMap [ newNodeIndexLower ] = newHangingNodeInfo;
      }
      
      // Update hanging node map for node over 0 
      RealType newZCoord = ( newNodeLower[2] + newNodeUpper[2] ) / 2.0;
      aol::Vec3 < RealType > nodeOverNode0 = this->getNodeCoords ( nodeInds[longestEdgeIndex] );
      nodeOverNode0[2] = newZCoord; 
      // Find points in nodeMap 
      auto gotNode0 = _nodeMap.find ( this->convertCoord3DToHash ( nodeOverNode0 ) );   
      if ( gotNode0 != _nodeMap.end() ) {          
        auto hangingNode0 = _hangingNodeMap.find ( gotNode0->second );  
        // If node exists and is hanging, update neighbor information 
        if ( hangingNode0 != _hangingNodeMap.end() ) {
          hangingNode0->second._neighbors.insert ( newElementIndex );
        }
      }
      
    }
    else if ( neighbor2[0] == -1 ) {
        
      // Check if there are one or two elements below 
      if ( neighborAbove.size() > 1 ) {
        // New hanging node appears -> lower new node 
        HangingNodeInfo < GlobalIndex > newHangingNodeInfo;  
        aol::Vec3 < RealType > auxNodeCoords = this->getNodeCoords ( newNodeIndexLower );
        auxNodeCoords[2] = this->getNodeCoord ( this->getElementNodeIndex ( neighborAbove[0], 3 ), 2 );
        // Find index of auxNodeCoords (node must exist) 
        auto auxNode = _nodeMap.find ( this->convertCoord3DToHash ( auxNodeCoords ) );
        newHangingNodeInfo._interpolationNodes = std::pair < GlobalIndex, GlobalIndex > ( newNodeIndexLower, auxNode->second );
        // Set neighbors of new hanging node = neighbor1
        for ( unsigned int j = 0; j < neighbor1.size(); ++j ) {
          newHangingNodeInfo._neighbors.insert ( neighbor1[j] ); 
        }
        _hangingNodeMap [ newNodeIndexUpper ] = newHangingNodeInfo;
      }
      
      // Update hanging node map for node over 0 
      RealType newZCoord = ( newNodeLower[2] + newNodeUpper[2] ) / 2.0;
      aol::Vec3 < RealType > nodeOverNode0 = this->getNodeCoords ( nodeInds[longestEdgeIndex] );
      nodeOverNode0[2] = newZCoord; 
      // Find points in nodeMap 
      auto gotNode0 = _nodeMap.find ( this->convertCoord3DToHash ( nodeOverNode0 ) );   
      if ( gotNode0 != _nodeMap.end() ) {          
        auto hangingNode0 = _hangingNodeMap.find ( gotNode0->second );  
        // If node exists and is hanging, update neighbor information 
        if ( hangingNode0 != _hangingNodeMap.end() ) {
          hangingNode0->second._neighbors.insert ( newElementIndex );
        }
      }
      
    }
    // If no edge neighbor is -1, no new hanging node appears, update hangingNodeMap 
    else {
            
      // Create node coordinates for possible hanging nodes 
      RealType newZCoord = ( newNodeLower[2] + newNodeUpper[2] ) / 2.0;
      aol::Vec3 < RealType > nodeOverNode0 = this->getNodeCoords ( nodeInds[longestEdgeIndex] );
      nodeOverNode0[2] = newZCoord;
      aol::Vec3 < RealType > nodeOverNode2 = this->getNodeCoords ( nodeInds[3-(localIndTmp+longestEdgeIndex)] );
      nodeOverNode2[2] = newZCoord;   
      // Find points in nodeMap 
      auto gotNode0 = _nodeMap.find ( this->convertCoord3DToHash ( nodeOverNode0 ) );   
      if ( gotNode0 != _nodeMap.end() ) {          
        auto hangingNode0 = _hangingNodeMap.find ( gotNode0->second );  
        // If node exists and is hanging, update neighbor information 
        if ( hangingNode0 != _hangingNodeMap.end() ) {
          hangingNode0->second._neighbors.insert ( newElementIndex );
        }
      }
      
      auto gotNode2 = _nodeMap.find ( this->convertCoord3DToHash ( nodeOverNode2 ) );   
      if ( gotNode2 != _nodeMap.end() ) {
        auto hangingNode2 = _hangingNodeMap.find ( gotNode2->second );  
        // If node exists and is hanging, update neighbor information 
        if ( hangingNode2 != _hangingNodeMap.end() ) {
          hangingNode2->second._neighbors.erase ( elIndex );
          hangingNode2->second._neighbors.insert ( newElementIndex );
        }
      }
      
        
    }  
      
  }
  
  // Find element containing two nodes 
  std::unordered_set < GlobalIndex > findElementsContainingNodes ( const GlobalIndex node1, const GlobalIndex node2 ) { 
    std::unordered_set < GlobalIndex > newSet; 
    
    // Go through all elements 
    bool found = false; 
    for ( int i = 0; i < this->getNumberOfElements(); ++i ) {
      found = false; 
      // Check if element contains node1 
      for ( int j = 0; j < 6; ++j ) {
        if ( this->getElementNodeIndex ( i, j ) == node1 ) {
          found = true; 
          break; 
        }
      }
      // If node1 was found, check if element contains node2 
      if ( found ) {
        found = false;
        for ( int j = 0; j < 6; ++j ) {
          if ( this->getElementNodeIndex ( i, j ) == node2 ) {
            found = true; 
            break;
          }
        }
      }
      if ( found ) {
        newSet.insert ( i ); 
      }
    }
    
    return newSet; 
  }
  
  // Refine neighbor of an element 
  void refineNeighbor ( const GlobalIndex neighborToRefine, const LocalIndex longestEdgeIndex, aol::Vec < 6, GlobalIndex > & nodeInds ) {
      
    LocalIndex localIndTmp = (longestEdgeIndex+1)%3;
              
    // Find shared edge and longest edge of neighbor 
    aol::Vec2 < GlobalIndex > sharedEdgeNodeInds ( nodeInds[localIndTmp], nodeInds[3-(localIndTmp+longestEdgeIndex)] );
    LocalIndex localIndTmp2 = this->getLongestEdgeIndex ( neighborToRefine );
    aol::Vec2 < GlobalIndex > longestEdgeNodeInds ( this->getElementNodeIndex ( neighborToRefine, (localIndTmp2+1)%3 ), this->getElementNodeIndex ( neighborToRefine, 3-((localIndTmp2+1)%3+localIndTmp2) ) );
    
    aol::Vec3 < RealType > P1 = this->getNodeCoords ( sharedEdgeNodeInds[0] );
    aol::Vec3 < RealType > P2 = this->getNodeCoords ( sharedEdgeNodeInds[1] );
    aol::Vec3 < RealType > P3 = this->getNodeCoords ( longestEdgeNodeInds[0] );
    aol::Vec3 < RealType > P4 = this->getNodeCoords ( longestEdgeNodeInds[1] );
    
    bool edgesCoincide1 = ( aol::Abs(P1[0]-P3[0])<1e-8 && aol::Abs(P1[1]-P3[1])<1e-8 && aol::Abs(P2[0]-P4[0])<1e-8 && aol::Abs(P2[1]-P4[1])<1e-8 );
    bool edgesCoincide2 = ( aol::Abs(P1[0]-P4[0])<1e-8 && aol::Abs(P1[1]-P4[1])<1e-8 && aol::Abs(P2[0]-P3[0])<1e-8 && aol::Abs(P2[1]-P3[1])<1e-8 );
    bool edgesCoincide = edgesCoincide1 | edgesCoincide2;
    
    // If shared edge = longest edge, only refine neighbor 
    if ( edgesCoincide ) {
      if ( !_visited[neighborToRefine] ) refineXY ( neighborToRefine );
    }
    // If shared edge != longest edge, refine neighbor and refine correct child 
    else {
      
      GlobalIndex newChildIndex = refineXY ( neighborToRefine );

      // Find child which contains sharedEdgeNodeInds and refine (can either be neighborToRefine or newChildIndex)
      aol::Vec3 < RealType > coordsSharedNode0 = this->getNodeCoords ( sharedEdgeNodeInds[0] );
      aol::Vec3 < RealType > coordsSharedNode1 = this->getNodeCoords ( sharedEdgeNodeInds[1] );
      
      bool child1ContainsNode1 = false, child1ContainsNode2 = false, child2ContainsNode1 = false, child2ContainsNode2 = false;
      // Check if child 1 contains sharedEdgeNodes
      for ( int i = 0; i < 3; ++i ) {
        aol::Vec3 < RealType > coords = this->getNodeCoords ( this->getNodeIndicesOfElement ( neighborToRefine )[i] );
        if ( coords[0] == coordsSharedNode0[0] ) {
          if ( coords[1] == coordsSharedNode0[1] ) {   
            child1ContainsNode1 = true; 
          }
        }
        if ( coords[0] == coordsSharedNode1[0] ) {
          if ( coords[1] == coordsSharedNode1[1] ) {   
            child1ContainsNode2 = true; 
          }
        } 
      }
      // Check if child 2 contains sharedEdgeNodes
      for ( int i = 0; i < 3; ++i ) {
        aol::Vec3 < RealType > coords = this->getNodeCoords ( this->getNodeIndicesOfElement ( newChildIndex )[i] );
        if ( coords[0] == coordsSharedNode0[0] ) {
          if ( coords[1] == coordsSharedNode0[1] ) {   
            child2ContainsNode1 = true; 
          }
        }
        if ( coords[0] == coordsSharedNode1[0] ) {
          if ( coords[1] == coordsSharedNode1[1] ) {   
            child2ContainsNode2 = true; 
          }
        } 
      }
      
      bool checkChild1 = child1ContainsNode1 & child1ContainsNode2;
      bool checkChild2 = child2ContainsNode1 & child2ContainsNode2;
      
      if ( !(checkChild1^checkChild2) ) throw aol::Exception ( "AdaptiveFEPrismMesh<>::refineNeighbor(): Both children contain refinement edge nodes, you did something stupid!", __FILE__, __LINE__ );    
      if ( checkChild1 ) {
        _markedForXYRefinement[neighborToRefine] = true; 
        refineXY ( neighborToRefine );
      }
      if ( checkChild2 ) {
        _markedForXYRefinement[newChildIndex] = true; 
        refineXY ( newChildIndex );  
      }
        
    }
  
  }
  
    
};

#endif
