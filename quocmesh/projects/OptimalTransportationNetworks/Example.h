#ifndef __EXAMPLE_H
#define __EXAMPLE_H

#include<aol.h>
#include<adaptiveFEPrismMesh.h>
#include<SpecialAdaptiveFEPrismMeshRPOp.h>
#include<adaptiveTriangMesh.h>


enum boundaryType { quadratic, circular, circularWithMidpoint, middle1To1, middle1To2, middle1To3, middle1To4, middle1To5, middle2To2 };

template < typename ConfiguratorType2D, typename ConfiguratorType3D > 
class Example {
    
    typedef typename ConfiguratorType3D::RealType RealType;
    typedef aol::Vector < RealType > VectorType;
    typedef typename ConfiguratorType2D::InitType GridType2D;
    typedef typename ConfiguratorType3D::InitType GridType3D;
    typedef SpecialAdaptiveFEPrismMeshRPOp < ConfiguratorType3D > RPType;
    typedef AdaptiveFEPrismMesh < RealType > GridType;
    typedef typename GridType2D::VertexIterator VertexIterator;
    typedef std::vector < aol::Vec3 < int > > CellArray;  
    
public:
    
    Example ( const string example, int initialGridLevelxy, const std::string problemType ) {
        _INITIALGRIDLEVELXY = initialGridLevelxy;
        _exampleType = example;
        _problemType = problemType; 
        if ( example == "1To1_quadratic_boundary" ) _boundaryType = quadratic;
        else if ( example == "1To1_quadratic_boundaryAlternative" ) _boundaryType = quadratic;
        else if ( example == "1To1_quadratic_boundary_oblique" ) _boundaryType = quadratic; 
        else if ( example == "1To1_quadratic_boundary_oblique2" ) _boundaryType = quadratic; 
        else if ( example == "1To1_quadratic_boundary_arc" ) _boundaryType = quadratic; 
        else if ( example == "1To2_quadratic_boundary" ) _boundaryType = quadratic;
        else if ( example == "2To2_quadratic_boundary" ) _boundaryType = quadratic;
        else if ( example == "2To2_quadratic_boundary_shiftLeft" ) _boundaryType = quadratic; 
        else if ( example == "2To2_quadratic_boundary_shiftRight" ) _boundaryType = quadratic; 
        else if ( example == "4To4_quadratic_boundary" ) _boundaryType = quadratic;
        else if ( example == "4To4_quadratic_boundary_TEST" ) _boundaryType = quadratic;
        else if ( example == "3To3_quadratic_boundary" ) _boundaryType = quadratic;
        else if ( example == "8To8_quadratic_boundary" ) _boundaryType = quadratic;
        else if ( example == "16To16_quadratic_boundary" ) _boundaryType = quadratic;
        else if ( example == "16To16_quadratic_boundary_asym" ) _boundaryType = quadratic;
        else if ( example == "16To16_quadratic_boundary_asym2" ) _boundaryType = quadratic;
        else if ( example == "16To16_quadratic_boundary_asym3" ) _boundaryType = quadratic;
        else if ( example == "1To1_circular_boundary" ) _boundaryType = circular;
        else if ( example == "1To2_circular_boundary" ) _boundaryType = circular;
        else if ( example == "3To3_circular_boundary" ) _boundaryType = circular;
        else if ( example == "1To4_circular_midpoint" ) _boundaryType = circularWithMidpoint;
        else if ( example == "1To8_circular_midpoint" ) _boundaryType = circularWithMidpoint;
        else if ( example == "1To16_circular_midpoint" ) _boundaryType = circularWithMidpoint; 
        else if ( example == "1To32_circular_midpoint" ) _boundaryType = circularWithMidpoint; 
        else if ( example == "1To32_circular_midpoint_asym" ) _boundaryType = circularWithMidpoint;
        else if ( example == "1To32_circular_midpoint_asym2" ) _boundaryType = circularWithMidpoint;
        else if ( example == "1To1_middle" ) _boundaryType = middle1To1;
        else if ( example == "1To2_middle" ) _boundaryType = middle1To2;
        else if ( example == "1To3_middle" ) _boundaryType = middle1To3;
        else if ( example == "1To4_middle" ) _boundaryType = middle1To4;
        else if ( example == "1To5_middle" ) _boundaryType = middle1To5;
        else if ( example == "2To2_middle" ) _boundaryType = middle2To2;
    }
    
    // Copy constructor
    Example ( Example & other ) {
        _exampleType = other._exampleType;
        _boundaryType = other._boundaryType;
        _problemType = other._problemType; 
    }
    
    
public: 
    
    int _INITIALGRIDLEVELXY; 
    std::string _exampleType;
    std::string _problemType; 
    boundaryType _boundaryType;
    aol::BitVector _boundaryMask;
    VectorType _boundaryVals;
     
    void setInitialImage ( RPType & rp, VectorType & image ) {
      // Go through all grid dofs 
      for ( typename RPType::DofIterator it ( rp ); it.notAtEnd(); ++it ) {
        int i = it.getDofIndex(); 
        aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( it.getNodeIndex() );  
        image.set ( i, this->setValue3D ( coords ) );    
      }
    }
    
    void setBoundaryMask ( RPType & rp, const int run = -1 ) {  
        _boundaryMask.resize ( rp.getGridRef().getNumberOfDofs() );
        _boundaryMask.setAll ( false );
        _boundaryVals.resize ( rp.getGridRef().getNumberOfDofs() );
        _boundaryVals.setAll ( 0.0 );
        VectorType boundaryImage ( rp.getGridRef().getNumberOfDofs() );
        // Quadratic boundaries 
        if ( _boundaryType == quadratic ) {
          // Go through all grid dofs 
          for ( typename RPType::DofIterator it ( rp ); it.notAtEnd(); ++it ) {
            int i = it.getDofIndex(); 
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( it.getNodeIndex() );
            // Set boundary mask and values 
            //if ( coords[0] == 0.0 || coords[0] == 1.0 || coords[1] == 0.0 || coords[1] == 1.0 || coords[2] == 0.0 || coords[2] == 1.0 ) { // OLD VERSION
            if ( coords[0] == 0.0 || coords[0] == 1.0 || coords[1] == 0.0 || coords[1] == 1.0 || coords[2] == 1.0 ) { // NEW VERSION
              _boundaryMask.set ( i, true );
              _boundaryVals.set ( i, this->setValue3D ( coords ) );
              boundaryImage.set ( i, 1.0 );
            }
          }
        }
        else if ( _boundaryType == circular ) { 
          // Go through all grid dofs 
          for ( typename RPType::DofIterator it ( rp ); it.notAtEnd(); ++it ) {
            int i = it.getDofIndex(); 
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( it.getNodeIndex() );
            if ( sqrt ( ( coords[0] - 0.5 ) * ( coords[0] - 0.5 ) + ( coords[1] - 0.5 ) * ( coords[1] - 0.5 ) ) >= 0.45 ) {
              _boundaryMask.set ( i, true ); 
              _boundaryVals.set ( i, this->setValue3D ( coords ) );
              boundaryImage.set ( i, 1.0 );
            }
          }
        }
        else if ( _boundaryType == circularWithMidpoint ) { 
          // Go through all grid dofs 
          for ( typename RPType::DofIterator it ( rp ); it.notAtEnd(); ++it ) {  
            int i = it.getDofIndex(); 
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( it.getNodeIndex() );
            if ( sqrt ( ( coords[0] - 0.5 ) * ( coords[0] - 0.5 ) + ( coords[1] - 0.5 ) * ( coords[1] - 0.5 ) ) >= 0.45 ) {
              _boundaryMask.set ( i, true ); 
              _boundaryVals.set ( i, this->setValue3D ( coords ) ); 
              boundaryImage.set ( i, 1.0 ); 
            }
            RealType runTilde = std::floor(0.5*static_cast<RealType>(run)+0.1);
            RealType diff = 1.0/32.0*pow(0.5,runTilde);
            //if ( coords[0] >= 0.5 && ( coords[1] >= 31.0/64.0 && coords[1] <= 33.0/64.0 ) ) {
            //if ( coords[0] >= 0.5 && ( coords[1] >= 30.0/64.0 && coords[1] <= 34.0/64.0 ) ) {
            if ( coords[0] >= 0.5 && ( coords[1] >= 0.5-diff && coords[1] <= 0.5+diff ) ) {
              _boundaryMask.set ( i, true ); 
              _boundaryVals.set ( i, this->setValue3D ( coords ) ); 
              boundaryImage.set ( i, 1.0 ); 
            }                
          } 
        }
        else if ( _boundaryType == middle1To1 ) {
          // Go through all grid dofs 
          bool check1, check2, check3; 
          for ( typename RPType::DofIterator it ( rp ); it.notAtEnd(); ++it ) {  
            int i = it.getDofIndex(); 
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( it.getNodeIndex() );
            check1 = ( coords[0] == 0.0 || coords[0] == 1.0 || coords[1] == 0.0 || coords[1] == 1.0 || coords[2] == 0.0 || coords[2] == 1.0 );
            check2 = ( ( coords[0] > 0.8 && coords[0] < 0.85 ) && ( coords[1] > 0.2 && coords[1] < 0.8 ) );
            check3 = ( ( coords[0] > 0.5 && coords[0] < 0.85 ) && ( ( coords[1] > 0.2 && coords[1] < 0.3 ) || ( coords[1] > 0.7 && coords[1] < 0.8 ) ) );
            if ( check1 || check2 || check3 ) { 
              _boundaryMask.set ( i, true ); 
              _boundaryVals.set ( i, this->setValue3D ( coords ) );
              boundaryImage.set ( i, 1.0 );
            } 
          }
        }
        else if ( _boundaryType == middle1To2 ) {
          // Go through all grid dofs 
          bool check1, check2, check3, check4; 
          for ( typename RPType::DofIterator it ( rp ); it.notAtEnd(); ++it ) {  
            int i = it.getDofIndex(); 
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( it.getNodeIndex() );  
            check1 = ( coords[0] == 0.0 || coords[0] == 1.0 || coords[1] == 0.0 || coords[1] == 1.0 || coords[2] == 1.0 );
            check2 = ( coords[0] > 0.125 && coords[0] < 0.175 && coords[1] > 0.075 && coords[1] < 0.925 );
            check3 = ( coords[0] > 0.125 && coords[0] < 0.825 && coords[1] > 0.875 && coords[1] < 0.925 );
            check4 = ( coords[0] > 0.125 && coords[0] < 0.825 && coords[1] > 0.075 && coords[1] < 0.125 );
            if ( check1 || check2 || check3 || check4 ) { 
              _boundaryMask.set ( i, true ); 
              _boundaryVals.set ( i, this->setValue3D ( coords ) );
              boundaryImage.set ( i, 1.0 );
            } 
          }
        }
        else if ( _boundaryType == middle1To3 ) { 
          // Go through all grid dofs 
          bool check1, check2, check3, check4; 
          for ( typename RPType::DofIterator it ( rp ); it.notAtEnd(); ++it ) {  
            int i = it.getDofIndex(); 
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( it.getNodeIndex() );  
            check1 = ( coords[0] == 0.0 || coords[0] == 1.0 || coords[1] == 0.0 || coords[1] == 1.0 || coords[2] == 1.0 );
            check2 = ( coords[0] > 0.125 && coords[0] <= 0.225 && coords[1] > 0.125 && coords[1] <= 0.875 );
            check3 = ( coords[0] > 0.125 && coords[0] <= 0.875 && coords[1] > 0.125 && coords[1] <= 0.25 );
            check4 = ( coords[0] > 0.775 && coords[0] <= 0.875 && coords[1] > 0.125 && coords[1] <= 0.875 );
            if ( check1 || check2 || check3 || check4 ) { 
              _boundaryMask.set ( i, true ); 
              _boundaryVals.set ( i, this->setValue3D ( coords ) );
              boundaryImage.set ( i, 1.0 );
            } 
          }  
        }
        else if ( _boundaryType == middle1To4 ) { 
          // Go through all grid dofs 
          bool check1, check2, check3, check4, check5; 
          RealType pi = aol::NumberTrait < RealType >::getPi();
          for ( typename RPType::DofIterator it ( rp ); it.notAtEnd(); ++it ) {  
            int i = it.getDofIndex(); 
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( it.getNodeIndex() );  
            RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 ); 
            check1 = ( coords[0] == 0.0 || coords[0] == 1.0 || coords[1] == 0.0 || coords[1] == 1.0 || coords[2] == 1.0 );
            check2 = ( coords[0] > 0.8 && coords[0] <= 0.9 && coords[1] > 0.1 && coords[1] <= 0.9 );
            check3 = ( coords[0] > 0.1 && coords[0] <= 0.9 && coords[1] > 0.1 && coords[1] <= 0.2 );
            check4 = ( coords[0] > 0.1 && coords[0] <= 0.9 && coords[1] > 0.8 && coords[1] <= 0.9 );
            check5 = ( coords[0] > 0.1 && coords[0] <= 0.2 && coords[1] > 0.1 && coords[1] <= 0.9 );
            if ( check1 || check2 || check3 || check4 || check5 ) { 
              _boundaryMask.set ( i, true ); 
              _boundaryVals.set ( i, this->setValue3D ( coords ) );
              boundaryImage.set ( i, 1.0 );
            }             
          }
        }
        else if ( _boundaryType == middle1To5 ) { 
          // Go through all grid dofs 
          bool check1, check2, check3, check4, check5; 
          RealType pi = aol::NumberTrait < RealType >::getPi();
          for ( typename RPType::DofIterator it ( rp ); it.notAtEnd(); ++it ) {  
            int i = it.getDofIndex(); 
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( it.getNodeIndex() );  
            RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 ); 
            check1 = ( coords[0] == 0.0 || coords[0] == 1.0 || coords[1] == 0.0 || coords[1] == 1.0 || coords[2] == 1.0 );
            check2 = ( coords[0] > 0.8 && coords[0] <= 0.9 && coords[1] > 0.1 && coords[1] <= 0.9 );
            check3 = ( coords[0] > 0.1 && coords[0] <= 0.9 && coords[1] > 0.1 && coords[1] <= 0.2 );
            check4 = ( coords[0] > 0.1 && coords[0] <= 0.9 && coords[1] > 0.8 && coords[1] <= 0.9 );
            check5 = ( coords[0] > 0.1 && coords[0] <= 0.2 && coords[1] > 0.1 && coords[1] <= 0.9 );
            if ( check1 || check2 || check3 || check4 || check5 ) { 
              _boundaryMask.set ( i, true ); 
              _boundaryVals.set ( i, this->setValue3D ( coords ) );
              boundaryImage.set ( i, 1.0 );
            } 
          }
        }
        else if ( _boundaryType == middle2To2 ) {
          bool check1, check2, check3, check4;
          for ( typename RPType::DofIterator it ( rp ); it.notAtEnd(); ++it ) {  
            int i = it.getDofIndex(); 
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( it.getNodeIndex() );  
            check1 = ( coords[0] <= 0.25 );
            check2 = ( coords[0] >= 0.75 );
            check3 = ( coords[1] <= 0.25 );
            check4 = ( coords[1] >= 0.75 );
            if ( check1 || check2 || check3 || check4 ) {
              _boundaryMask.set ( i, true ); 
              _boundaryVals.set ( i, this->setValue3D ( coords ) );
              boundaryImage.set ( i, 1.0 );   
            }
          }
        }
    }
    
    
    RealType setValue3D ( const aol::Vec3 < RealType > coords ) { 
      RealType value2D = 0.0, value3D = 0.0;  
      if ( _exampleType == "1To1_quadratic_boundary" ) {
        value2D = 0.75;
        if ( coords[0] > 0.5 ) value2D = 0.25;  
      }
      else if ( _exampleType == "1To1_quadratic_boundaryAlternative" ) {
        value2D = 1.0;
        if ( coords[0] > 0.5 ) value2D = 0.0;  
      }
      else if ( _exampleType == "1To1_quadratic_boundary_oblique" ) {
        value2D = 1.0;
        if ( coords[0] > 0.0 && coords[1] < 1.0 ) value2D = 0.0;
        if ( coords[0] == 1.0 ) value2D = 0.0;
      } 
      else if ( _exampleType == "1To1_quadratic_boundary_oblique2" ) {
        value2D = 1.0;
        if ( coords[0] > 5.0/8.0 ) value2D = 0.0; 
        if ( coords[0] > 3.0/8.0 && coords[1] == 0.0 ) value2D = 0.0;
      }
      else if ( _exampleType == "1To1_quadratic_boundary_arc" ) { 
        value2D = 1.0;
        if ( coords[0] > 0.25+sqrt(0.25-(coords[1]-0.5)*(coords[1]-0.5)) ) value2D = 0.0;
      }
      else if ( _exampleType == "1To2_quadratic_boundary" ) {
        value2D = 0.75; 
        if ( coords[0] > 0.5 ) value2D = 0.25;
        if ( coords[0] > 0.0 && coords[0] < 1.0 && coords[1] == 0.0 ) value2D = 0.5; 
      }
      else if ( _exampleType == "2To2_quadratic_boundary" ) { 
        value2D = 0.75; 
        if ( coords[0] > 0.25 && coords[0] <= 0.75 ) value2D = 0.5;
        if ( coords[0] > 0.75 ) value2D = 0.25;
      }
      else if ( _exampleType == "2To2_quadratic_boundary_shiftLeft" ) {
        value2D = 0.75; 
        if ( coords[0] > 0.03125 && coords[0] <= 0.5+0.03125 ) value2D = 0.5;
        if ( coords[0] > 0.5+0.03125 ) value2D = 0.25;
      }
      else if ( _exampleType == "2To2_quadratic_boundary_shiftRight" ) {
        value2D = 0.75; 
        if ( coords[0] > 0.5-0.03125 && coords[0] <= 1.0-0.03125 ) value2D = 0.5;
        if ( coords[0] > 1.0-0.03125 ) value2D = 0.25;
      }   
      else if ( _exampleType == "2To2_quadratic_boundaryAlternative" ) { 
        value2D = 1.0; 
        if ( coords[0] > 0.25 && coords[0] <= 0.75 ) value2D = 0.5;
        if ( coords[0] > 0.75 ) value2D = 0.0;
      }
      else if ( _exampleType == "4To4_quadratic_boundary" ) {  
        value2D = 0.75;
        if ( coords[0] > 0.2 && coords[0] <= 0.4 ) value2D = 0.625;
        if ( coords[0] > 0.4 && coords[0] <= 0.6 ) value2D = 0.5; 
        if ( coords[0] > 0.6 && coords[0] <= 0.8 ) value2D = 0.375; 
        if ( coords[0] > 0.8 ) value2D = 0.25;
      }
      else if ( _exampleType == "4To4_quadratic_boundaryAlternative" ) {  
        value2D = 1.0;
        if ( coords[0] > 0.2 && coords[0] <= 0.4 ) value2D = 0.75;
        if ( coords[0] > 0.4 && coords[0] <= 0.6 ) value2D = 0.5; 
        if ( coords[0] > 0.6 && coords[0] <= 0.8 ) value2D = 0.25; 
        if ( coords[0] > 0.8 ) value2D = 0.0;
      }
      else if ( _exampleType == "4To4_quadratic_boundary_TEST" ) {
        value2D = 1.0;
        if ( coords[0] > 1.0/5.0 && coords[0] <= 2.0/5.0 ) value2D = 0.875;
        if ( coords[0] > 2.0/5.0 && coords[0] <= 3.0/5.0 ) value2D = 0.5;
        if ( coords[0] > 3.0/5.0 && coords[0] <= 4.0/5.0 ) value2D = 0.125;
        if ( coords[0] > 4.0/5.0 ) value2D = 0.0;
      }
      else if ( _exampleType == "3To3_quadratic_boundary" ) { 
        value2D = 0.75;
        if ( coords[0] > 0.25 && coords[0] <= 0.5 ) value2D = 0.625;
        if ( coords[0] > 0.5 && coords[0] <= 0.75 ) value2D = 0.375;
        if ( coords[0] > 0.75 ) value2D = 0.25;
      }
      else if ( _exampleType == "8To8_quadratic_boundary" ) {
        value2D = 1.0;
        if ( coords[0] > 1.0/9.0 && coords[0] <= 2.0/9.0 ) value2D = 0.875; 
        if ( coords[0] > 2.0/9.0 && coords[0] <= 3.0/9.0 ) value2D = 0.75; 
        if ( coords[0] > 3.0/9.0 && coords[0] <= 4.0/9.0 ) value2D = 0.625; 
        if ( coords[0] > 4.0/9.0 && coords[0] <= 5.0/9.0 ) value2D = 0.5; 
        if ( coords[0] > 5.0/9.0 && coords[0] <= 6.0/9.0 ) value2D = 0.375; 
        if ( coords[0] > 6.0/9.0 && coords[0] <= 7.0/9.0 ) value2D = 0.25; 
        if ( coords[0] > 7.0/9.0 && coords[0] <= 8.0/9.0 ) value2D = 0.125; 
        if ( coords[0] > 8.0/9.0 ) value2D = 0.0;    
      }
      else if ( _exampleType == "16To16_quadratic_boundary" ) {
        value2D = 1.0;
        if ( coords[0] > 1.0/17.0 && coords[0] <= 2.0/17.0 ) value2D = 15.0/16.0; 
        if ( coords[0] > 2.0/17.0 && coords[0] <= 3.0/17.0 ) value2D = 14.0/16.0; 
        if ( coords[0] > 3.0/17.0 && coords[0] <= 4.0/17.0 ) value2D = 13.0/16.0; 
        if ( coords[0] > 4.0/17.0 && coords[0] <= 5.0/17.0 ) value2D = 12.0/16.0; 
        if ( coords[0] > 5.0/17.0 && coords[0] <= 6.0/17.0 ) value2D = 11.0/16.0; 
        if ( coords[0] > 6.0/17.0 && coords[0] <= 7.0/17.0 ) value2D = 10.0/16.0; 
        if ( coords[0] > 7.0/17.0 && coords[0] <= 8.0/17.0 ) value2D = 9.0/16.0;   
        if ( coords[0] > 8.0/17.0 && coords[0] <= 9.0/17.0 ) value2D = 8.0/16.0; 
        if ( coords[0] > 9.0/17.0 && coords[0] <= 10.0/17.0 ) value2D = 7.0/16.0; 
        if ( coords[0] > 10.0/17.0 && coords[0] <= 11.0/17.0 ) value2D = 6.0/16.0; 
        if ( coords[0] > 11.0/17.0 && coords[0] <= 12.0/17.0 ) value2D = 5.0/16.0; 
        if ( coords[0] > 12.0/17.0 && coords[0] <= 13.0/17.0 ) value2D = 4.0/16.0; 
        if ( coords[0] > 13.0/17.0 && coords[0] <= 14.0/17.0 ) value2D = 3.0/16.0; 
        if ( coords[0] > 14.0/17.0 && coords[0] <= 15.0/17.90 ) value2D = 2.0/16.0;  
        if ( coords[0] > 15.0/17.0 && coords[0] <= 16.0/17.0 ) value2D = 1.0/16.0;  
        if ( coords[0] > 16.0/17.0 ) value2D = 0.0;  
      }
      else if ( _exampleType == "16To16_quadratic_boundary_asym" ) {
        value2D = 1.0;
        if ( coords[0] > 0.057 && coords[0] <= 0.124 ) value2D = 15.0/16.0; 
        if ( coords[0] > 0.124 && coords[0] <= 0.181 ) value2D = 14.0/16.0; 
        if ( coords[0] > 0.181 && coords[0] <= 0.238 ) value2D = 13.0/16.0; 
        if ( coords[0] > 0.238 && coords[0] <= 0.305 ) value2D = 12.0/16.0; 
        if ( coords[0] > 0.305 && coords[0] <= 0.362 ) value2D = 11.0/16.0; 
        if ( coords[0] > 0.362 && coords[0] <= 0.419 ) value2D = 10.0/16.0; 
        if ( coords[0] > 0.419 && coords[0] <= 0.476 ) value2D = 9.0/16.0;   
        if ( coords[0] > 0.476 && coords[0] <= 0.543 ) value2D = 8.0/16.0; 
        if ( coords[0] > 0.543 && coords[0] <= 0.6 ) value2D = 7.0/16.0; 
        if ( coords[0] > 0.6 && coords[0] <= 0.657 ) value2D = 6.0/16.0; 
        if ( coords[0] > 0.657 && coords[0] <= 0.714 ) value2D = 5.0/16.0; 
        if ( coords[0] > 0.714 && coords[0] <= 0.781 ) value2D = 4.0/16.0; 
        if ( coords[0] > 0.781 && coords[0] <= 0.838 ) value2D = 3.0/16.0; 
        if ( coords[0] > 0.838 && coords[0] <= 0.895 ) value2D = 2.0/16.0;  
        if ( coords[0] > 0.895 && coords[0] <= 0.952 ) value2D = 1.0/16.0;    
        if ( coords[0] > 0.952 ) value2D = 0.0;  
      }
      else if ( _exampleType == "16To16_quadratic_boundary_asym2" ) {
        value2D = 1.0;
        if ( coords[0] > 0.0568 && coords[0] <= 0.1136 ) value2D = 15.0/16.0; 
        if ( coords[0] > 0.1136 && coords[0] <= 0.1704 ) value2D = 14.0/16.0; 
        if ( coords[0] > 0.1704 && coords[0] <= 0.2272 ) value2D = 13.0/16.0; 
        if ( coords[0] > 0.2272 && coords[0] <= 0.2954 ) value2D = 12.0/16.0; 
        if ( coords[0] > 0.2954 && coords[0] <= 0.3522 ) value2D = 11.0/16.0; 
        if ( coords[0] > 0.3522 && coords[0] <= 0.409 ) value2D = 10.0/16.0; 
        if ( coords[0] > 0.409 && coords[0] <= 0.4658 ) value2D = 9.0/16.0;   
        if ( coords[0] > 0.4658 && coords[0] <= 0.534 ) value2D = 8.0/16.0; 
        if ( coords[0] > 0.534 && coords[0] <= 0.5908 ) value2D = 7.0/16.0; 
        if ( coords[0] > 0.5908 && coords[0] <= 0.6476 ) value2D = 6.0/16.0; 
        if ( coords[0] > 0.6476 && coords[0] <= 0.7044 ) value2D = 5.0/16.0; 
        if ( coords[0] > 0.7044 && coords[0] <= 0.7726 ) value2D = 4.0/16.0; 
        if ( coords[0] > 0.7726 && coords[0] <= 0.8294 ) value2D = 3.0/16.0; 
        if ( coords[0] > 0.8294 && coords[0] <= 0.8862 ) value2D = 2.0/16.0;  
        if ( coords[0] > 0.8862 && coords[0] <= 0.943 ) value2D = 1.0/16.0;    
        if ( coords[0] > 0.943 ) value2D = 0.0;  
      }
      else if ( _exampleType == "16To16_quadratic_boundary_asym3" ) {
        value2D = 1.0;
        if ( coords[0] > 1.0/32.0 && coords[0] <= 3.0/32.0 ) value2D = 15.0/16.0; 
        if ( coords[0] > 3.0/32.0 && coords[0] <= 5.0/32.0 ) value2D = 14.0/16.0; 
        if ( coords[0] > 5.0/32.0 && coords[0] <= 7.0/32.0 ) value2D = 13.0/16.0; 
        if ( coords[0] > 7.0/32.0 && coords[0] <= 9.0/32.0 ) value2D = 12.0/16.0; 
        if ( coords[0] > 9.0/32.0 && coords[0] <= 11.0/32.0 ) value2D = 11.0/16.0; 
        if ( coords[0] > 11.0/32.0 && coords[0] <= 13.0/32.0 ) value2D = 10.0/16.0; 
        if ( coords[0] > 13.0/32.0 && coords[0] <= 15.0/32.0 ) value2D = 9.0/16.0;   
        if ( coords[0] > 15.0/32.0 && coords[0] <= 17.0/32.0 ) value2D = 8.0/16.0; 
        if ( coords[0] > 17.0/32.0 && coords[0] <= 19.0/32.0 ) value2D = 7.0/16.0; 
        if ( coords[0] > 19.0/32.0 && coords[0] <= 21.0/32.0 ) value2D = 6.0/16.0; 
        if ( coords[0] > 21.0/32.0 && coords[0] <= 23.0/32.0 ) value2D = 5.0/16.0; 
        if ( coords[0] > 23.0/32.0 && coords[0] <= 25.0/32.0 ) value2D = 4.0/16.0; 
        if ( coords[0] > 25.0/32.0 && coords[0] <= 27.0/32.0 ) value2D = 3.0/16.0; 
        if ( coords[0] > 27.0/32.0 && coords[0] <= 29.0/32.0 ) value2D = 2.0/16.0;  
        if ( coords[0] > 29.0/32.0 && coords[0] <= 31.0/32.0 ) value2D = 1.0/16.0;    
        if ( coords[0] > 31.0/32.0 ) value2D = 0.0;    
      }
      else if ( _exampleType == "1To1_circular_boundary" ) {
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 );
        if ( ( psi >= 0.0 && psi <= 3.0/4.0*pi ) || ( psi < 0 && psi > -1.0/4.0*pi ) ) value2D = 0.25;
        else if ( psi > 3.0/4.0*pi || psi <= 1.0/4.0*pi ) value2D = 0.75; 
      }
      else if ( _exampleType == "1To2_circular_boundary" ) {
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 ); 
        if ( ( psi > 0.0 && psi <= 2.0/3.0*pi ) || psi <= -2.0/3.0*pi ) value2D = 0.75;
        else if ( psi > 2.0/3.0*pi || psi <= -2.0/3.0*pi ) value2D = 0.5;
        else if ( psi > -2.0/3.0*pi && psi <= 0 ) value2D = 0.25;
      }
      else if ( _exampleType == "3To3_circular_boundary" ) {
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 );   
        if ( psi > -5.0/12.0*pi+pi/4.0+pi/30.0 && psi <= 7.0/12.0*pi-pi/4.0+pi/30.0 ) value2D = 0.375;  
        if ( psi > 7.0/12.0*pi-pi/4.0+pi/30.0 && psi <= 7.0/12.0*pi ) value2D = 0.25; 
        if ( psi > 7.0/12.0*pi && psi <= 7.0/12.0*pi+pi/4.0+pi/30.0 ) value2D = 0.75;
        if ( psi > 7.0/12.0*pi+pi/4.0+pi/30.0 || psi <= -5.0/12.0*pi-pi/4.0+pi/30.0 ) value2D = 0.625; 
        if ( psi > -5.0/12.0*pi-pi/4.0+pi/30.0 && psi <= -5.0/12.0*pi ) value2D = 0.75; 
        if ( psi > -5.0/12.0*pi && psi <= -5.0/12.0*pi+pi/4.0+pi/30.0 ) value2D = 0.25; 
      }
      else if ( _exampleType == "1To4_circular_midpoint" ) {
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 );   
        if ( psi > 0.0 && psi <= pi/4.0-pi/20.0 ) value2D = 0.75;
        if ( psi > pi/4.0-pi/20.0 && psi <= 3.0/4.0*pi+pi/20.0 ) value2D = 0.625; 
        if ( psi > 3.0/4.0*pi+pi/20.0 || psi <= -3.0/4.0*pi-pi/20.0 ) value2D = 0.5;
        if ( psi > -3.0/4.0*pi-pi/20.0 && psi <= -pi/4.0+pi/20.0 ) value2D = 0.375; 
        if ( psi > -pi/4.0+pi/20.0 && psi <= 0.0 ) value2D = 0.25; 
      }
      else if ( _exampleType == "1To8_circular_midpoint" ) { 
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 );  
        if ( psi > 0.0 && psi <= pi/8.0 ) value2D = 1.0;  
        if ( psi > pi/8.0 && psi <= 3.0/8.0*pi ) value2D = 0.875; 
        if ( psi > 3.0/8.0*pi && psi <= 5.0/8.0*pi ) value2D = 0.75; 
        if ( psi > 5.0/8.0*pi && psi <= 7.0/8.0*pi ) value2D = 0.625;  
        if ( psi > 7.0/8.0*pi || psi <= -7.0/8.0*pi ) value2D = 0.5;  
        if ( psi > -7.0/8.0*pi && psi <= -5.0/8.0*pi ) value2D = 0.375;            
        if ( psi > -5.0/8.0*pi && psi <= -3.0/8.0*pi ) value2D = 0.25;  
        if ( psi > -3.0/8.0*pi && psi <= -pi/8.0 ) value2D = 0.125;  
        if ( psi > -pi/8.0 && psi <= 0.0 ) value2D = 0.0;
      }
      else if ( _exampleType == "1To16_circular_midpoint" ) { 
        cout << "WE ARE HERE" << endl; 
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 );    
        if ( psi > 0.0 && psi <= pi/16.0 ) value2D = 1.0;  
        if ( psi > pi/16.0 && psi <= 3.0/16.0*pi ) value2D = 15.0/16.0;; 
        if ( psi > 3.0/16.0*pi && psi <= 5.0/16.0*pi ) value2D = 14.0/16.0;; 
        if ( psi > 5.0/16.0*pi && psi <= 7.0/16.0*pi ) value2D = 13.0/16.0;;  
        if ( psi > 7.0/16.0*pi && psi <= 9.0/16.0*pi ) value2D = 12.0/16.0;; 
        if ( psi > 9.0/16.0*pi && psi <= 11.0/16.0*pi ) value2D = 11.0/16.0;;  
        if ( psi > 11.0/16.0*pi && psi <= 13.0/16.0*pi ) value2D = 10.0/16.0;; 
        if ( psi > 13.0/16.0*pi && psi <= 15.0/16.0*pi ) value2D = 9.0/16.0;; 
        if ( psi > 15.0/16.0*pi || psi <= -15.0/16.0*pi ) value2D = 8.0/16.0;; 
        if ( psi > -15.0/16.0*pi && psi <= -13.0/16.0*pi ) value2D = 7.0/16.0;; 
        if ( psi > -13.0/16.0*pi && psi <= -11.0/16.0*pi ) value2D = 6.0/16.0;; 
        if ( psi > -11.0/16.0*pi && psi <= -9.0/16.0*pi ) value2D = 5.0/16.0;; 
        if ( psi > -9.0/16.0*pi && psi <= -7.0/16.0*pi ) value2D = 4.0/16.0;; 
        if ( psi > -7.0/16.0*pi && psi <= -5.0/16.0*pi ) value2D = 3.0/16.0;; 
        if ( psi > -5.0/16.0*pi && psi <= -3.0/16.0*pi ) value2D = 2.0/16.0;; 
        if ( psi > -3.0/16.0*pi && psi <= -1.0/16.0*pi ) value2D = 1.0/16.0;; 
        if ( psi > -1.0/16.0*pi && psi <= 0.0 ) value2D = 0.0; 
        cout << "Psi = " << psi << endl; 
        cout << "Value2D = " << value2D << endl; 
      }
      else if ( _exampleType == "1To32_circular_midpoint" ) { 
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 );    
        if ( psi > 0.0 && psi <= pi/32.0 ) value2D = 1.0;  
        if ( psi > pi/32.0 && psi <= 3.0/32.0*pi ) value2D = 31.0/32.0; 
        if ( psi > 3.0/32.0*pi && psi <= 5.0/32.0*pi ) value2D = 30.0/32.0; 
        if ( psi > 5.0/32.0*pi && psi <= 7.0/32.0*pi ) value2D = 29.0/32.0;  
        if ( psi > 7.0/32.0*pi && psi <= 9.0/32.0*pi ) value2D = 28.0/32.0;  
        if ( psi > 9.0/32.0*pi && psi <= 11.0/32.0*pi ) value2D = 27.0/32.0;  
        if ( psi > 11.0/32.0*pi && psi <= 13.0/32.0*pi ) value2D = 26.0/32.0;  
        if ( psi > 13.0/32.0*pi && psi <= 15.0/32.0*pi ) value2D = 25.0/32.0;  
        if ( psi > 15.0/32.0*pi && psi <= 17.0/32.0*pi ) value2D = 24.0/32.0;   
        if ( psi > 17.0/32.0*pi && psi <= 19.0/32.0*pi ) value2D = 23.0/32.0;   
        if ( psi > 19.0/32.0*pi && psi <= 21.0/32.0*pi ) value2D = 22.0/32.0;
        if ( psi > 21.0/32.0*pi && psi <= 23.0/32.0*pi ) value2D = 21.0/32.0;
        if ( psi > 23.0/32.0*pi && psi <= 25.0/32.0*pi ) value2D = 20.0/32.0;
        if ( psi > 25.0/32.0*pi && psi <= 27.0/32.0*pi ) value2D = 19.0/32.0;
        if ( psi > 27.0/32.0*pi && psi <= 29.0/32.0*pi ) value2D = 18.0/32.0;
        if ( psi > 29.0/32.0*pi && psi <= 31.0/32.0*pi ) value2D = 17.0/32.0;
        if ( psi > 31.0/32.0*pi || psi <= -31.0/32.0*pi ) value2D = 16.0/32.0;
        if ( psi > -31.0/32.0*pi && psi <= -29.0/32.0*pi ) value2D = 15.0/32.0;
        if ( psi > -29.0/32.0*pi && psi <= -27.0/32.0*pi ) value2D = 14.0/32.0; 
        if ( psi > -27.0/32.0*pi && psi <= -25.0/32.0*pi ) value2D = 13.0/32.0;
        if ( psi > -25.0/32.0*pi && psi <= -23.0/32.0*pi ) value2D = 12.0/32.0;
        if ( psi > -23.0/32.0*pi && psi <= -21.0/32.0*pi ) value2D = 11.0/32.0;
        if ( psi > -21.0/32.0*pi && psi <= -19.0/32.0*pi ) value2D = 10.0/32.0;
        if ( psi > -19.0/32.0*pi && psi <= -17.0/32.0*pi ) value2D = 9.0/32.0;
        if ( psi > -17.0/32.0*pi && psi <= -15.0/32.0*pi ) value2D = 8.0/32.0;
        if ( psi > -15.0/32.0*pi && psi <= -13.0/32.0*pi ) value2D = 7.0/32.0;
        if ( psi > -13.0/32.0*pi && psi <= -11.0/32.0*pi ) value2D = 6.0/32.0;
        if ( psi > -11.0/32.0*pi && psi <= -9.0/32.0*pi ) value2D = 5.0/32.0;
        if ( psi > -9.0/32.0*pi && psi <= -7.0/32.0*pi ) value2D = 4.0/32.0;
        if ( psi > -7.0/32.0*pi && psi <= -5.0/32.0*pi ) value2D = 3.0/32.0;
        if ( psi > -5.0/32.0*pi && psi <= -3.0/32.0*pi ) value2D = 2.0/32.0;
        if ( psi > -3.0/32.0*pi && psi <= -1.0/32.0*pi ) value2D = 1.0/32.0;
        if ( psi > -1.0/32.0*pi && psi <= 0.0 ) value2D = 0.0; 
      }
      else if ( _exampleType == "1To32_circular_midpoint_asym" ) { 
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 );   
        RealType epsilon = pi/64.0;
        RealType h = 0.25*(pi/4.0-2.0*epsilon);
        RealType hTilde = h/2.0+epsilon;
        
        if ( psi > 0.0 && psi <= hTilde ) value2D = 1.0;
        if ( psi > hTilde && psi <= hTilde+h ) value2D = 31.0/32.0;
        if ( psi > hTilde+h && psi <= hTilde+2.0*h ) value2D = 30.0/32.0;
        if ( psi > hTilde+2.0*h && psi <= hTilde+3.0*h ) value2D = 29.0/32.0;
        
        if ( psi > 1.0*hTilde+3.0*h && psi <= 3.0*hTilde+3.0*h ) value2D = 28.0/32.0;
        if ( psi > 3.0*hTilde+3.0*h && psi <= 3.0*hTilde+4.0*h ) value2D = 27.0/32.0;
        if ( psi > 3.0*hTilde+4.0*h && psi <= 3.0*hTilde+5.0*h ) value2D = 26.0/32.0;
        if ( psi > 3.0*hTilde+5.0*h && psi <= 3.0*hTilde+6.0*h ) value2D = 25.0/32.0;
        
        if ( psi > 3.0*hTilde+6.0*h && psi <= 5.0*hTilde+6.0*h ) value2D = 24.0/32.0;
        if ( psi > 5.0*hTilde+6.0*h && psi <= 5.0*hTilde+7.0*h ) value2D = 23.0/32.0;
        if ( psi > 5.0*hTilde+7.0*h && psi <= 5.0*hTilde+8.0*h ) value2D = 22.0/32.0;
        if ( psi > 5.0*hTilde+8.0*h && psi <= 5.0*hTilde+9.0*h ) value2D = 21.0/32.0;
        
        if ( psi > 5.0*hTilde+9.0*h && psi <= 7.0*hTilde+9.0*h ) value2D = 20.0/32.0;
        if ( psi > 7.0*hTilde+9.0*h && psi <= 7.0*hTilde+10.0*h ) value2D = 19.0/32.0;
        if ( psi > 7.0*hTilde+10.0*h && psi <= 7.0*hTilde+11.0*h ) value2D = 18.0/32.0;
        if ( psi > 7.0*hTilde+11.0*h && psi <= 7.0*hTilde+12.0*h ) value2D = 17.0/32.0;
        
        if ( psi > 7.0*hTilde+12.0*h || psi <= -7.0*hTilde-12.0*h ) value2D = 16.0/32.0;
        if ( psi > -7.0*hTilde-12.0*h && psi <= -7.0*hTilde-11.0*h ) value2D = 15.0/32.0;
        if ( psi > -7.0*hTilde-11.0*h && psi <= -7.0*hTilde-10.0*h ) value2D = 14.0/32.0;
        if ( psi > -7.0*hTilde-10.0*h && psi <= -7.0*hTilde-9.0*h ) value2D = 13.0/32.0;
        
        if ( psi > -7.0*hTilde-9.0*h && psi <= -5.0*hTilde-9.0*h ) value2D = 12.0/32.0;
        if ( psi > -5.0*hTilde-9.0*h && psi <= -5.0*hTilde-8.0*h ) value2D = 11.0/32.0;
        if ( psi > -5.0*hTilde-8.0*h && psi <= -5.0*hTilde-7.0*h ) value2D = 10.0/32.0;
        if ( psi > -5.0*hTilde-7.0*h && psi <= -5.0*hTilde-6.0*h ) value2D = 9.0/32.0;
        
        if ( psi > -5.0*hTilde-6.0*h && psi <= -3.0*hTilde-6.0*h ) value2D = 8.0/32.0;
        if ( psi > -3.0*hTilde-6.0*h && psi <= -3.0*hTilde-5.0*h ) value2D = 7.0/32.0;
        if ( psi > -3.0*hTilde-5.0*h && psi <= -3.0*hTilde-4.0*h ) value2D = 6.0/32.0;
        if ( psi > -3.0*hTilde-4.0*h && psi <= -3.0*hTilde-3.0*h ) value2D = 5.0/32.0;
        
        if ( psi > -3.0*hTilde-3.0*h && psi <= -1.0*hTilde-3.0*h ) value2D = 4.0/32.0;
        if ( psi > -1.0*hTilde-3.0*h && psi <= -1.0*hTilde-2.0*h ) value2D = 3.0/32.0;
        if ( psi > -1.0*hTilde-2.0*h && psi <= -1.0*hTilde-1.0*h ) value2D = 2.0/32.0;
        if ( psi > -1.0*hTilde-1.0*h && psi <= -1.0*hTilde ) value2D = 1.0/32.0;
        if ( psi > -1.0*hTilde && psi <= 0.0 ) value2D = 0.0;
      }
      else if ( _exampleType == "1To32_circular_midpoint_asym2" ) { 
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 );   
        RealType delta = pi/64.0;
        RealType epsilon = pi/128.0;
        RealType h = pi/16.0 - delta/4.0 - epsilon/8.0;
        RealType hTilde = h+delta;
        RealType hHat = hTilde+epsilon;
        
        if ( psi > 0.0 && psi <= hHat/2.0 ) value2D = 1.0;
        if ( psi > hHat/2.0 && psi <= hHat/2.0+1.0*h ) value2D = 31.0/32.0;
        if ( psi > hHat/2.0+1.0*h && psi <= hHat/2.0+2.0*h ) value2D = 30.0/32.0;
        if ( psi > hHat/2.0+2.0*h && psi <= hHat/2.0+3.0*h ) value2D = 29.0/32.0;
        
        if ( psi > hHat/2.0+3.0*h && psi <= hHat/2.0+3.0*h+hTilde ) value2D = 28.0/32.0;
        if ( psi > hHat/2.0+3.0*h+hTilde && psi <= hHat/2.0+4.0*h+hTilde ) value2D = 27.0/32.0;
        if ( psi > hHat/2.0+4.0*h+hTilde && psi <= hHat/2.0+5.0*h+hTilde ) value2D = 26.0/32.0;
        if ( psi > hHat/2.0+5.0*h+hTilde && psi <= hHat/2.0+6.0*h+hTilde ) value2D = 25.0/32.0;
        
        
        if ( psi > hHat/2.0+6.0*h+hTilde && psi <= 3.0*hHat/2.0+6.0*h+hTilde ) value2D = 24.0/32.0;
        if ( psi > 3.0*hHat/2.0+6.0*h+hTilde && psi <= 3.0*hHat/2.0+7.0*h+hTilde ) value2D = 23.0/32.0;
        if ( psi > 3.0*hHat/2.0+7.0*h+hTilde && psi <= 3.0*hHat/2.0+8.0*h+hTilde ) value2D = 22.0/32.0;
        if ( psi > 3.0*hHat/2.0+8.0*h+hTilde && psi <= 3.0*hHat/2.0+9.0*h+hTilde ) value2D = 21.0/32.0;
        
        if ( psi > 3.0*hHat/2.0+9.0*h+hTilde && psi <= 3.0*hHat/2.0+9.0*h+2.0*hTilde ) value2D = 20.0/32.0;
        if ( psi > 3.0*hHat/2.0+9.0*h+2.0*hTilde && psi <= 3.0*hHat/2.0+10.0*h+2.0*hTilde ) value2D = 19.0/32.0;
        if ( psi > 3.0*hHat/2.0+10.0*h+2.0*hTilde && psi <= 3.0*hHat/2.0+11.0*h+2.0*hTilde ) value2D = 18.0/32.0;
        if ( psi > 3.0*hHat/2.0+11.0*h+2.0*hTilde && psi <= 3.0*hHat/2.0+12.0*h+2.0*hTilde ) value2D = 17.0/32.0;
        
        
        if ( psi > 3.0*hHat/2.0+12.0*h+2.0*hTilde || psi <= -3.0*hHat/2.0-12.0*h-2.0*hTilde ) value2D = 16.0/32.0;
        if ( psi > -3.0*hHat/2.0-12.0*h-2.0*hTilde && psi <= -3.0*hHat/2.0-11.0*h-2.0*hTilde ) value2D = 15.0/32.0;
        if ( psi > -3.0*hHat/2.0-11.0*h-2.0*hTilde && psi <= -3.0*hHat/2.0-10.0*h-2.0*hTilde ) value2D = 14.0/32.0; 
        if ( psi > -3.0*hHat/2.0-10.0*h-2.0*hTilde && psi <= -3.0*hHat/2.0-9.0*h-2.0*hTilde ) value2D = 13.0/32.0;
        
        if ( psi > -3.0*hHat/2.0-9.0*h-2.0*hTilde && psi <= -3.0*hHat/2.0-9.0*h-1.0*hTilde ) value2D = 12.0/32.0;
        if ( psi > -3.0*hHat/2.0-9.0*h-1.0*hTilde && psi <= -3.0*hHat/2.0-8.0*h-1.0*hTilde ) value2D = 11.0/32.0;
        if ( psi > -3.0*hHat/2.0-8.0*h-1.0*hTilde && psi <= -3.0*hHat/2.0-7.0*h-1.0*hTilde ) value2D = 10.0/32.0;
        if ( psi > -3.0*hHat/2.0-7.0*h-1.0*hTilde && psi <= -3.0*hHat/2.0-6.0*h-1.0*hTilde ) value2D = 9.0/32.0;
        
        
        if ( psi > -3.0*hHat/2.0-6.0*h-1.0*hTilde && psi <= -1.0*hHat/2.0-6.0*h-1.0*hTilde ) value2D = 8.0/32.0;
        if ( psi > -1.0*hHat/2.0-6.0*h-1.0*hTilde && psi <= -1.0*hHat/2.0-5.0*h-1.0*hTilde ) value2D = 7.0/32.0;
        if ( psi > -1.0*hHat/2.0-5.0*h-1.0*hTilde && psi <= -1.0*hHat/2.0-4.0*h-1.0*hTilde ) value2D = 6.0/32.0;
        if ( psi > -1.0*hHat/2.0-4.0*h-1.0*hTilde && psi <= -1.0*hHat/2.0-3.0*h-1.0*hTilde ) value2D = 5.0/32.0;
        
        if ( psi > -1.0*hHat/2.0-3.0*h-1.0*hTilde && psi <= -1.0*hHat/2.0-3.0*h ) value2D = 4.0/32.0;
        if ( psi > -1.0*hHat/2.0-3.0*h && psi <= -1.0*hHat/2.0-2.0*h ) value2D = 3.0/32.0;
        if ( psi > -1.0*hHat/2.0-2.0*h && psi <= -1.0*hHat/2.0-1.0*h ) value2D = 2.0/32.0;
        if ( psi > -1.0*hHat/2.0-1.0*h && psi <= -1.0*hHat/2.0 ) value2D = 1.0/32.0;
        if ( psi > -1.0*hHat/2.0 && psi <= 0.0 ) value2D = 0.0;
      }
      else if ( _exampleType == "1To1_middle" ) {
        value2D = 0.75;
        if ( coords[0] > 0.125 && coords[0] < 0.875 && coords[1] > 0.125 && coords[1] < 0.875 ) value2D = 0.5;  
      } 
      else if ( _exampleType == "1To2_middle" ) {
        value2D = 0.5;
        if ( coords[0] > 0.125 && coords[0] < 0.825 && coords[1] >= 0.075 && coords[1] < 0.5 ) value2D = 1.0;
        if ( coords[0] > 0.125 && coords[0] < 0.825 && coords[1] >= 0.5 && coords[1] <= 0.925 ) value2D = 0.0;
      }
      else if ( _exampleType == "1To3_middle" ) { 
        value2D = 0.75;
        if ( coords[1] <= 0.25 ) value2D = 0.25;
        if ( coords[0] <= 0.225 && coords[1] > 0.25 && coords[1] <= 0.75 ) value2D = 0.5;
        if ( coords[0] > 0.775 && coords[1] > 0.25 && coords[1] <= 0.75 ) value2D = 1.0;
      }
      else if ( _exampleType == "1To4_middle" ) { 
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 ); 
        RealType tol = 0.025;
        value2D = 0.5;
        if ( psi >= pi/5.0+tol && psi < 3.0*pi/5.0 ) value2D = 0.75;
        if ( psi >= 3.0*pi/5.0 && psi < pi ) value2D = 1.0;
        if ( psi >= -pi && psi < -3.0*pi/5.0 ) value2D = 0.0;
        if ( psi >= -3.0*pi/5.0 && psi < -pi/5.0-tol ) value2D = 0.25;
      }
      else if ( _exampleType == "1To5_middle" ) { 
        RealType pi = aol::NumberTrait < RealType >::getPi();
        RealType psi = std::atan2 ( coords[1]-0.5, coords[0]-0.5 );   
        RealType tol = 0.1;
        value2D = 0.625;
        if ( psi > pi/3.0+tol && psi <= 2.0*pi/3.0+tol ) value2D = 0.75;
        if ( psi > 2.0*pi/3.0+tol && psi <= pi ) value2D = 0.875;
        if ( psi > -pi && psi <= -2.0*pi/3.0-tol ) value2D = 0.25;
        if ( psi > -2.0*pi/3.0-tol && psi <= -pi/3.0-tol ) value2D = 0.375;
        if ( psi > -pi/3.0-tol && psi <= 0.0+tol ) value2D = 0.5;
      }
      else if ( _exampleType == "2To2_middle" ) {
        value2D = 0.75;
        if ( coords[0] >= 0.375 && coords[0] < 0.625 ) value2D = 0.5;
        if ( coords[0] >= 0.625 ) value2D = 0.25;
      }
      // Lift value 
      if ( value2D > coords[2] ) value3D = 1.0; 
      return value3D;        
    }
    
    void deliftImage ( RPType & rp, VectorType & image3D, VectorType & image2D, const int run ) { 
 
      // Build adaptive 2D grid 
      int nodeCounter = 1;
      aol::Vec3 < int > nodeInds2D;
      aol::MultiVector < RealType > nodes2D ( 3, 1 );
      aol::Vec3 < RealType > startCoords ( 0.0, 0.0, 0.0 );
      nodes2D.set ( 0, startCoords );
      std::unordered_map < long unsigned int, int > nodeHashMap;
      std::pair < long unsigned int, int > startHash ( rp.getGridRef().convertCoord3DToHash ( startCoords ), 0 );
      nodeHashMap.insert ( startHash );
      CellArray elements2D; 
      for ( typename RPType::GridType::ElementIterator elit ( rp.getGridRef() ); elit.notAtEnd(); ++elit ) {
        aol::Vec < 6, int > nodeInds = rp.getGridRef().getNodeIndicesOfElement ( elit.getIndex() );
        // If first node has z coord 0, element is ground element  
        if ( rp.getGridRef().getNodeCoords ( nodeInds[0] )[2] == 0.0 ) {
          // Go through nodes, check if they already exist and push them to nodes2D 
          for ( int i = 0; i < 3; ++i ) {
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( nodeInds[i] );  
            auto gotNode = nodeHashMap.find ( rp.getGridRef().convertCoord3DToHash ( coords ) );    
            // Node found: Get index 
            if ( gotNode != nodeHashMap.end() ) {
              nodeInds2D[i] = gotNode->second;
            }
            // Node not found: Push to node list and create new index 
            else {
              nodeInds2D[i] = nodeCounter; 
              ++nodeCounter; 
              std::pair < long unsigned int, int > newHash ( rp.getGridRef().convertCoord3DToHash ( coords ), nodeInds2D[i] );
              nodeHashMap.insert ( newHash );
              for ( short j = 0; j < 3; ++j ) {
                nodes2D[j].pushBack ( coords[j] );   
              }
            }
          }  
          // Push element to elements2D 
          elements2D.push_back ( nodeInds2D );
        }
      }
      
      // Save 2D grid 
      AdaptiveFETriangMesh < RealType > grid2D ( nodes2D, elements2D );
      
      // Get extended line segment list (containing all nodes) 
      rp.setExtendedLineSegmentList();      
      
      // Go through 2D grid and get corresponding line segment 
      rp.extendVectorZConst ( image3D );
      image2D.resize ( grid2D.getNumVertices() );
      RealType value2D; 
      RealType hz = 1.0 / static_cast < RealType > ( 1 << rp.getStopLevelz() );
      for ( typename AdaptiveFETriangMesh < RealType >::NodeIteratorType nodeit ( grid2D ); nodeit.notAtEnd(); ++nodeit ) {
          aol::Vec3 < RealType > refCoords = grid2D.getVertex ( nodeit.getIndex() );
          std::map < int, int >::const_iterator lineSegment; 
          // Find corresponding line segment 
          for ( unsigned int i = 0; i < rp.getLineSegmentList().size(); ++i ) {
            std::map < int, int >::const_iterator lineSegStart = rp.getLineSegmentList()[i].begin();  
            aol::Vec3 < RealType > coords = rp.getGridRef().getNodeCoords ( lineSegStart->second ); 
            if ( aol::Abs(coords[0]-refCoords[0]) < 1e-8 && aol::Abs(coords[1]-refCoords[1]) < 1e-8 ) {
              lineSegment = rp.getLineSegmentList()[i].begin();
              break;   
            }
          }
          
          // Delift and set value of image2D
          value2D = 0.0;
          
          while ( image3D[lineSegment->second] > 0.9 ) {
            ++lineSegment;
            value2D = hz*static_cast<RealType>(lineSegment->first);
          }        
          image2D.set ( nodeit.getIndex(), value2D );
      }
      
      stringstream filename;
      filename << "Results/solution" << "_" << run << ".vtu"; 
      std::string filenameStr = filename.str();
      grid2D.saveAsVTK ( filenameStr, image2D ); 
        
    }
    

};

#endif