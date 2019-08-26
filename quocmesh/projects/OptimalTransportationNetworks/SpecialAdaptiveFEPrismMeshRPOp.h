#ifndef __SPECIALADAPTIVEFEPRISMMESHRPOP_H
#define __SPECIALADAPTIVEFEPRISMMESHRPOP_H

#include<aol.h>
#include<adaptiveFEPrismMeshRPOp.h>

template < typename ConfiguratorType > 
class SpecialAdaptiveFEPrismMeshRPOp : public AdaptiveFEPrismMeshRPOp < ConfiguratorType > { 

public:
    
    typedef typename ConfiguratorType::InitType GridType;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::InitType::MapType MapType;
    typedef aol::SparseMatrix < int > MatrixType;
    typedef aol::SparseMatrix < RealType > RealMatrixType;
    typedef typename AdaptiveFEPrismMeshRPOp < ConfiguratorType >::IndexMapType IndexMapType;
    typedef aol::Vector < RealType > VectorType;
    typedef aol::Vector < int > IntVectorType;
    
    SpecialAdaptiveFEPrismMeshRPOp ( const GridType &grid, const ConfiguratorType &configurator, const int stopLevelxy, const int stopLevelz ) : AdaptiveFEPrismMeshRPOp < ConfiguratorType > ( grid, configurator ) {
      _stopLevelxy = stopLevelxy;
      _stopLevelz = stopLevelz;
    }
    
    SpecialAdaptiveFEPrismMeshRPOp ( const GridType &grid, const ConfiguratorType &configurator ) : AdaptiveFEPrismMeshRPOp < ConfiguratorType > ( grid, configurator ) {}
    
    inline int getStopLevelxy () const {
      return _stopLevelxy;
    }
    
    inline int getStopLevelz () const {
      return _stopLevelz;
    }
    
    void setLineSegmentList () {    
        // Go through all elements of grid 
        aol::Vec3 < RealType > coords3D;
        int zCoordInt;
        unsigned long int place, nodeIndex, dofIndex;
        std::unordered_map < unsigned long int, int > hashmap; 
        for ( typename GridType::NodeIterator nodeit ( this->getGridRef() ); nodeit.notAtEnd(); ++nodeit ) { 
            nodeIndex = nodeit.getIndex(); 
            if ( !this->getGridRef().isHangingNode ( nodeIndex ) ) {
                // Get coordinates of new point and find corresponding place in current lineSegmentList
                dofIndex = this->getRestrictedGlobalNodeIndex ( nodeIndex );
                coords3D = nodeit.getCoords();
                place = this->getGridRef().convertCoord2DToHash ( coords3D );
                auto mapit = hashmap.find ( place );    
                zCoordInt = static_cast < int > ( coords3D[2] * ( 1 << _stopLevelz ) );
                
                // If place exists: Add new point to lineSegmentList[place]
                if ( mapit != hashmap.end() ) {  
                    _lineSegmentList [ mapit->second ].insert ( std::pair < int, int > ( zCoordInt, dofIndex ) );
                }
                // If not: Push back to lineSegmentList and create new entry
                else {
                    hashmap.insert ( { place, _lineSegmentList.size() } );
                    std::map < int, int > newEntry;
                    newEntry [ zCoordInt ] = dofIndex;
                    _lineSegmentList.push_back ( newEntry );
                }
            }
        }
    }
    
    // ATTENTION: Clears existing line segment list 
    void setExtendedLineSegmentList () {
        _lineSegmentList.clear();
         // Go through all elements of grid 
        aol::Vec3 < RealType > coords3D;
        int zCoordInt;
        unsigned long int place, nodeIndex;
        std::unordered_map < unsigned long int, int > hashmap; 
        for ( typename GridType::NodeIterator nodeit ( this->getGridRef() ); nodeit.notAtEnd(); ++nodeit ) { 
            // Get coordinates of new point and find corresponding place in current lineSegmentList
            nodeIndex = nodeit.getIndex(); 
            coords3D = nodeit.getCoords();
            place = this->getGridRef().convertCoord2DToHash ( coords3D );
            auto mapit = hashmap.find ( place );    
            zCoordInt = static_cast < int > ( coords3D[2] * ( 1 << _stopLevelz ) );
            
            // If place exists: Add new point to lineSegmentList[place]
            if ( mapit != hashmap.end() ) {  
                _lineSegmentList [ mapit->second ].insert ( std::pair < int, int > ( zCoordInt, nodeIndex ) );
            }
            // If not: Push back to lineSegmentList and create new entry
            else {
                hashmap.insert ( { place, _lineSegmentList.size() } );
                std::map < int, int > newEntry;
                newEntry [ zCoordInt ] = nodeIndex;
                _lineSegmentList.push_back ( newEntry );
            }
        }
    }
    
    // Set constraint list for a uniform grid in RUN=0
    void setConstraintList () {
        int lineSegmentLength;
        std::vector < std::pair < int, int > > constraintVector;
        for ( unsigned int i = 0; i < _lineSegmentList.size(); ++i ) {
            constraintVector.clear();
            lineSegmentLength = _lineSegmentList[i].size();
            for ( int t1 = 0; t1 < lineSegmentLength-1; ++t1 ) {
                for ( int t2 = t1+1; t2 < lineSegmentLength; ++t2 ) {
                    std::pair < int, int > constraint ( t1, t2 );
                    constraintVector.push_back ( constraint );  
                }
            }
            _constraintList.push_back ( constraintVector );
        }  
    }
    
    // Set full constraint list for a uniform grid in RUN=0
    void setConstraintListFull () {
        _constraintListFull.clear();
        int lineSegmentLength;
        std::vector < std::pair < int, int > > constraintVector;
        for ( unsigned int i = 0; i < _lineSegmentList.size(); ++i ) {
            constraintVector.clear();
            lineSegmentLength = _lineSegmentList[i].size();
            for ( int t1 = 0; t1 < lineSegmentLength-1; ++t1 ) {
                for ( int t2 = t1+1; t2 < lineSegmentLength; ++t2 ) {
                    std::pair < int, int > constraint ( t1, t2 );
                    constraintVector.push_back ( constraint );  
                }
            }
            _constraintListFull.push_back ( constraintVector );
        }  
    }
    
    std::vector < std::map < int, int > > & getLineSegmentList () {
      return _lineSegmentList;
    }
    
    std::vector < std::vector < std::pair < int, int > > > & getConstraintList () {
      return _constraintList;
    }
    
    std::vector < std::vector < std::pair < int, int > > > & getConstraintListFull () {
      return _constraintListFull;    
    }
    
    // Update constraint list after grid was refined or within iteration 
    void updateConstraintListBT ( const aol::Vector < RealType > & powerLookuptable, aol::Vector < RealType > &phi1, aol::Vector < RealType > &phi2, const RealType tol, const bool start ) {
    
        if ( start ) {
            _constraintList.clear();
            _constraintList.resize ( _lineSegmentList.size() );
        }
  
        RealType A, B, intervalLength, valA, valB;
        RealType hz = 1.0 / ( static_cast < RealType > ( 1 << _stopLevelz ) );
        int totalIntervalLength;
        std::map < int, int >::const_iterator mapit; 
        std::map < int, int >::const_iterator mapit2;
        int constraintsAdded = 0;
  
        // Go through all line segments
        for ( unsigned int i = 0; i < _lineSegmentList.size(); ++i ) {
      
            int lineSegmentLength = _lineSegmentList[i].size();
            std::map < int, int >::const_iterator mapit_t1 = _lineSegmentList[i].begin();
            std::map < int, int >::const_iterator mapit_t2; 
        
            for ( int t1 = 0; t1 < lineSegmentLength-1; ++t1, ++mapit_t1 ) {
                
                mapit_t2 = mapit_t1;
                ++mapit_t2;
        
                for ( int t2 = t1+1; t2 < lineSegmentLength; ++t2, ++mapit_t2 ) {
                    
                    std::pair < int, int > constraint ( t1, t2 );
            
                    if ( start ) {
                    
                        // Check if constraint is active 
                        A = 0.0;
                        B = 0.0;
                        mapit = mapit_t1;
                        mapit2 = mapit_t1;  
                        ++mapit2;
                        for ( int k = t1; k < t2; ++k, ++mapit, ++mapit2 ) {
                            intervalLength = static_cast < RealType > ( mapit2->first - mapit->first );
                            A += phi1[mapit->second] * intervalLength;
                            B += phi2[mapit->second] * intervalLength;
                        }
                        A *= hz;
                        B *= hz;
                        totalIntervalLength = mapit_t2->first - mapit_t1->first;
                        
                        if ( sqrt ( A*A + B*B ) >= powerLookuptable[totalIntervalLength]*static_cast<RealType>(totalIntervalLength)*hz - tol ) {
                            _constraintList[i].push_back ( constraint );   
                            ++constraintsAdded;
                        }
                    
                    }
                    else {    
                    
                        auto it = std::find ( _constraintList[i].begin(), _constraintList[i].end(), constraint );
                        
                        if ( it == _constraintList[i].end() ) {
                            
                            // Check if constraint is active 
                            A = 0.0;
                            B = 0.0;
                            mapit = mapit_t1;
                            mapit2 = mapit_t1;
                            ++mapit2;
                            for ( int k = t1; k < t2; ++k, ++mapit, ++mapit2 ) {
                                intervalLength = static_cast < RealType > ( mapit2->first - mapit->first );
                                A += phi1[mapit->second] * intervalLength;
                                B += phi2[mapit->second] * intervalLength;
                            }
                            A *= hz;
                            B *= hz;
                            totalIntervalLength = mapit_t2->first - mapit_t1->first;
                            
                            if ( sqrt ( A*A + B*B ) >= powerLookuptable[totalIntervalLength]*static_cast<RealType>(totalIntervalLength)*hz - tol ) {
                                _constraintList[i].push_back ( constraint );   
                                ++constraintsAdded;
                            }
                        
                        }
                    
                    }
          
                }
      
            }
    
        }
  
    }  
    
    // Update constraint list after grid was refined or within iteration 
    void updateConstraintListUP ( const RealType epsilon, const RealType a, aol::Vector < RealType > &phi1, aol::Vector < RealType > &phi2, const RealType tol, const bool start ) {
    
        if ( start ) {
            _constraintList.clear();
            _constraintList.resize ( _lineSegmentList.size() );
        }
  
        this->extendVectorZConst ( phi1 );
        this->extendVectorZConst ( phi2 );
  
        RealType A, B, intervalLength, valA, valB;
        RealType hz = 1.0 / ( static_cast < RealType > ( 1 << _stopLevelz ) );
        int totalIntervalLength;
        std::map < int, int >::const_iterator mapit; 
        std::map < int, int >::const_iterator mapit2;
        int constraintsAdded = 0;
  
        // Go through all line segments
        for ( unsigned int i = 0; i < _lineSegmentList.size(); ++i ) {
      
            int lineSegmentLength = _lineSegmentList[i].size();
            std::map < int, int >::const_iterator mapit_t1 = _lineSegmentList[i].begin();
            std::map < int, int >::const_iterator mapit_t2; 
        
            for ( int t1 = 0; t1 < lineSegmentLength-1; ++t1, ++mapit_t1 ) {
                
                mapit_t2 = mapit_t1;
                ++mapit_t2;
        
                for ( int t2 = t1+1; t2 < lineSegmentLength; ++t2, ++mapit_t2 ) {
                    
                    std::pair < int, int > constraint ( t1, t2 );
            
                    if ( start ) {
                    
                        // Check if constraint is active 
                        A = 0.0;
                        B = 0.0;
                        mapit = mapit_t1;
                        mapit2 = mapit_t1;  
                        ++mapit2;
                        for ( int k = t1; k < t2; ++k, ++mapit, ++mapit2 ) {
                            intervalLength = static_cast < RealType > ( mapit2->first - mapit->first );
                            A += phi1[mapit->second] * intervalLength;
                            B += phi2[mapit->second] * intervalLength;
                        }
                        A *= hz;
                        B *= hz;
                        totalIntervalLength = mapit_t2->first - mapit_t1->first;
                        
                        if ( sqrt ( A*A + B*B ) >= aol::Min ( a*static_cast<RealType>(totalIntervalLength)*hz, static_cast<RealType>(totalIntervalLength)*hz + epsilon ) - tol ) { 
                            _constraintList[i].push_back ( constraint );   
                            ++constraintsAdded;
                        }
                    
                    }
                    else {    
                    
                        std::vector < std::pair < int, int > >::const_iterator it = std::find ( _constraintList[i].begin(), _constraintList[i].end(), constraint );
                        
                        if ( it == _constraintList[i].end() ) {
                            
                            // Check if constraint is active 
                            A = 0.0;
                            B = 0.0;
                            mapit = mapit_t1;
                            mapit2 = mapit_t1;
                            ++mapit2;
                            for ( int k = t1; k < t2; ++k, ++mapit, ++mapit2 ) {
                                intervalLength = static_cast < RealType > ( mapit2->first - mapit->first );
                                A += phi1[mapit->second] * intervalLength;
                                B += phi2[mapit->second] * intervalLength;
                            }
                            A *= hz;
                            B *= hz;
                            totalIntervalLength = mapit_t2->first - mapit_t1->first;
                            
                            if ( sqrt ( A*A + B*B ) >= aol::Min ( a*static_cast<RealType>(totalIntervalLength)*hz, static_cast<RealType>(totalIntervalLength)*hz + epsilon ) - tol ) { 
                                _constraintList[i].push_back ( constraint );   
                                ++constraintsAdded;
                            }
                        
                        }
                    
                    }
          
                }
      
            }
    
        }

        this->restrictVectorZConst ( phi1 );
        this->restrictVectorZConst ( phi2 );  
  
    }   

    
    
public: 
    
    int _stopLevelxy;
    int _stopLevelz;
    std::vector < std::map < int, int > > _lineSegmentList;
    std::vector < std::vector < std::pair < int, int > > > _constraintList;
    std::vector < std::vector < std::pair < int, int > > > _constraintListFull;
    
};

#endif