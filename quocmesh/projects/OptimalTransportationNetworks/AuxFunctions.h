#ifndef __AUXFUNCTIONS_H
#define __AUXFUNCTIONS_H


#include<aol.h>
#include<iostream>
#include<algorithm>
#include<vector>
#include<scalarArray.h>
#include<FEOpInterface.h>
#include<configurators.h>
#include<Newton.h>
#include<omp.h>
#include<time.h>
#include<smallMat.h>
#include<math.h>
#include<tiffio.h>
//#include<suiteSparseSolver.h>

#include<SpecialAdaptiveFEPrismMeshRPOp.h>
#include<PDMinimizationSolver.h>
#include<Example.h>

using namespace std;
using namespace qc;
using namespace aol;
typedef double RealType;

//-----------------------------------------------------------------------------------------------------------------------------------------

/*
template < typename ConfiguratorType >
class StiffMatrixAssembler2D : public aol::FELinScalarWeightedStiffInterface < ConfiguratorType, StiffMatrixAssembler2D < ConfiguratorType > > {
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  GridType & _grid;    
    
public:
      
  StiffMatrixAssembler2D ( ConfiguratorType & configurator, GridType & grid ) 
    : aol::FELinScalarWeightedStiffInterface < ConfiguratorType, StiffMatrixAssembler2D < ConfiguratorType > > ( configurator, grid ), _grid ( grid ) { }  
   
  RealType getCoeff(const typename ConfiguratorType::ElementType& El, int QuadPoint, const typename ConfiguratorType::DomVecType & RefCoord) const {
    return 1.0;
  } 
   
};

template < typename ConfiguratorType >
class MixedMassStiffMatrixAssembler2DX : public aol::FELinScalarWeightedSemiDiffInterface < ConfiguratorType, MixedMassStiffMatrixAssembler2DX < ConfiguratorType > > {
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  GridType & _grid;    
    
public:
      
  MixedMassStiffMatrixAssembler2DX ( ConfiguratorType & configurator, GridType & grid ) 
    : aol::FELinScalarWeightedSemiDiffInterface < ConfiguratorType, MixedMassStiffMatrixAssembler2DX < ConfiguratorType > > ( configurator, grid, 0, true ), _grid ( grid ) { }  
   
  RealType getCoeff(const typename ConfiguratorType::ElementType& El, int QuadPoint, const typename ConfiguratorType::DomVecType & RefCoord) const {
    return 1.0;
  } 
   
};

template < typename ConfiguratorType >
class MixedMassStiffMatrixAssembler2DY : public aol::FELinScalarWeightedSemiDiffInterface < ConfiguratorType, MixedMassStiffMatrixAssembler2DY < ConfiguratorType > > {
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  GridType & _grid;    
    
public:
      
  MixedMassStiffMatrixAssembler2DY ( ConfiguratorType & configurator, GridType & grid ) 
    : aol::FELinScalarWeightedSemiDiffInterface < ConfiguratorType, MixedMassStiffMatrixAssembler2DY < ConfiguratorType > > ( configurator, grid, 1, true ), _grid ( grid ) { }  
   
  RealType getCoeff(const typename ConfiguratorType::ElementType& El, int QuadPoint, const typename ConfiguratorType::DomVecType & RefCoord) const {
    return 1.0;
  } 
   
};

void writeTIFF ( const std::string & filename, const qc::ScalarArray < RealType, qc::QC_2D > & image ) {

   TIFF *tif= TIFFOpen(filename.c_str(), "w");

   TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, image.getNumX() );
   TIFFSetField (tif, TIFFTAG_IMAGELENGTH, image.getNumY() );
   TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
   TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);
   TIFFSetField (tif, TIFFTAG_ROWSPERSTRIP, 1);
   TIFFSetField (tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
   TIFFSetField (tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
   TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 32);
   TIFFSetField (tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
   TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

   tsize_t strip_size = TIFFStripSize (tif);
   tstrip_t strips_num = TIFFNumberOfStrips (tif);

   float* strip_buf=static_cast<float*>( _TIFFmalloc(strip_size) );
   for (unsigned int s=0; s<strips_num; s++) {
     for (int col=0; col<image.getNumX(); col++) {
       strip_buf[col]=static_cast<float>( image.get(col,s) );
     }
     TIFFWriteEncodedStrip (tif, s, strip_buf, strip_size);
   }
   _TIFFfree(strip_buf);
   TIFFClose(tif);

}

aol::Vec2 < int > convertCoordsToMatrixIndex ( const aol::Vec3 < RealType > & coords, const RealType dx ) {
  aol::Vec2 < int > inds ( 0, 0 );
  inds[1] = static_cast < int > ( coords[1]/dx + 0.5 );
  inds[0] = static_cast < int > ( coords[0]/dx + 0.5 );
  return inds;
}

template < typename GridType, typename VectorType >
void createImageFromVectorFieldLaplace ( const GridType & grid, const VectorType & phi1, const VectorType & phi2, const int levelxy ) {
 
  // Create 2D grid   
  typedef FETriangMesh < RealType > GridType2D; 
  typedef TriangleIntegration < RealType, qc::QC_2D, 2 > QuadRule; 
  typedef TriangMeshConfigurator < RealType, GridType2D, QuadRule > ConfiguratorType2D;
  GridType2D grid2D ( levelxy );
  ConfiguratorType2D configurator2D ( grid2D );
  int gridDofs2D = grid2D.getNumVertices();
    
  // Create 2D psi = phi^\perp   
  VectorType psi1_2D ( gridDofs2D );
  VectorType psi2_2D ( gridDofs2D );
  int index = 0;
  int index2D = 0;
  for ( typename GridType::NodeIterator it ( grid ); it.notAtEnd(); ++it, ++index ) { 
    aol::Vec3 < RealType > coords = it.getCoords();
    // Get slice s = 1/2
    if ( abs(coords[2]-0.5) <= 1e-5 ) {
      // Find corresponding ground node index 
      index2D = 0;
      for ( GridType2D::VertexIterator it ( grid2D ); it.notAtEnd(); ++it, ++index2D ) {
        aol::Vec3 < RealType > coords2D = it.getCoords();
        if ( aol::Abs ( coords2D[0] - coords[0] ) <= 1e-5 && aol::Abs ( coords2D[1] - coords[1] ) <= 1e-5 ) {
          psi1_2D[index2D] = phi2[index];
          psi2_2D[index2D] = -phi1[index];
        }
      }
    } 
  }
  grid2D.saveAsVTK ( "Results/psi1Upper.vtu", psi1_2D );
  grid2D.saveAsVTK ( "Results/psi2Upper.vtu", psi2_2D );
  
  // Create 2D FE matrices 
  StiffMatrixAssembler2D < ConfiguratorType2D > stiffMatrixAssembler ( configurator2D, grid2D );
  MixedMassStiffMatrixAssembler2DX < ConfiguratorType2D > mixedMassStiffMatrixAssemblerX ( configurator2D, grid2D );
  MixedMassStiffMatrixAssembler2DY < ConfiguratorType2D > mixedMassStiffMatrixAssemblerY ( configurator2D, grid2D );
  SparseMatrix < RealType > stiffMat ( gridDofs2D, gridDofs2D );
  SparseMatrix < RealType > mixedMatXT ( gridDofs2D, gridDofs2D );
  SparseMatrix < RealType > mixedMatYT ( gridDofs2D, gridDofs2D );
  stiffMatrixAssembler.assembleAddMatrix ( stiffMat );
  mixedMassStiffMatrixAssemblerX.assembleAddMatrix ( mixedMatXT );
  mixedMassStiffMatrixAssemblerY.assembleAddMatrix ( mixedMatYT );
  
  // Solve 2D problem to obtain image u 
  Vector < RealType > u ( gridDofs2D );
  Vector < RealType > RHS ( gridDofs2D ); 
  mixedMatXT.apply ( psi1_2D, RHS ); 
  mixedMatYT.applyAdd ( psi2_2D, RHS );
  UMFPACKInverseOp < RealType, SparseMatrix < RealType > > op;
  op.setMatrix ( stiffMat );
  op.apply ( RHS, u );
  grid2D.saveAsVTK ( "Results/u.vtu", u );
    
  // Save image as matrix 
  int numNodesXY = ( 1 << levelxy ) + 1 ;
  RealType dx = 1.0/static_cast<RealType>(1<<levelxy);
  ScalarArray < RealType, qc::QC_2D > uMat ( numNodesXY, numNodesXY );
  index2D = 0;
  for ( GridType2D::VertexIterator it ( grid2D ); it.notAtEnd(); ++it, ++index2D ) { 
    // Get node coords and convert to matrix indices (i,j) 
    Vec3 < RealType > coords = it.getCoords();  
    Vec2 < int > matrixInds = convertCoordsToMatrixIndex ( coords, dx );
    uMat.set ( matrixInds[0], matrixInds[1], u[index2D] );
  }
  
  // Set min(uMat) to 0
  RealType minVal = uMat.getMinValue();
  for ( int i = 0; i < numNodesXY; ++i ) {
    for ( int j = 0; j < numNodesXY; ++j ) {
      uMat.set ( i, j, uMat.get(i,j) - minVal );
    }
  }
  writeTIFF ( "Results/imageFromVectorField.tiff", uMat ); 
  
}
*/

//-----------------------------------------------------------------------------------------------------------------------------------------
template < typename VectorType >
void setPowerLookuptable ( VectorType & powerLookuptable, const int MAXLEVELZ, const RealType epsilon ) {
    int num = 1 << MAXLEVELZ;
    RealType hz = 1.0 / ( static_cast < RealType > ( num ) );
    powerLookuptable.resize ( num + 1 );
    powerLookuptable[0] = 0.0;
    #pragma omp parallel for 
    for ( int i = 1; i < num; ++i ) {
        powerLookuptable[i] = pow ( hz*static_cast<RealType>(i), -epsilon );
    }
    powerLookuptable[num] = 1.0; 
}

//-----------------------------------------------------------------------------------------------------------------------------------------
template < typename ConfiguratorType > 
void saveImageAsVTK ( AdaptiveFEPrismMesh < RealType > & grid, SpecialAdaptiveFEPrismMeshRPOp < ConfiguratorType > & rp, aol::Vector < RealType > & vector, const std::string savename, const int run = -1, const bool start = false ) {
    if ( !start ) {
        rp.extendVectorZConst ( vector );
    }
    stringstream filename;
    if ( run == -1 ) {
        filename << "Results/" << savename << ".vtu"; 
    }
    else {
        filename << "Results/" << savename << "_" << run << ".vtu"; 
    }
    std::string filenameStr = filename.str();
    grid.saveAsLegacyVTK ( filenameStr, vector );   
    if ( !start ) {
        rp.restrictVectorZConst ( vector );
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------
template < typename ConfiguratorType > 
void saveImageAsVTKVectorField ( AdaptiveFEPrismMesh < RealType > & grid, SpecialAdaptiveFEPrismMeshRPOp < ConfiguratorType > & rp, aol::Vector < RealType > & vectorX, aol::Vector < RealType > & vectorY, aol::Vector < RealType > & vectorZ, const std::string savename, const int run = -1, const bool start = false ) {
    if ( !start ) {
        rp.extendVectorZConst ( vectorX );
        rp.extendVectorZConst ( vectorY );
        rp.extendVectorZConst ( vectorZ );
    }
    stringstream filename;
    if ( run == -1 ) {
        filename << "Results/" << savename << ".vtu"; 
    }
    else {
        filename << "Results/" << savename << "_" << run << ".vtu"; 
    }
    std::string filenameStr = filename.str();
    grid.saveAsLegacyVTKVectorField ( filenameStr, vectorX, vectorY, vectorZ );   
    if ( !start ) {
        rp.restrictVectorZConst ( vectorX );
        rp.restrictVectorZConst ( vectorY );
        rp.restrictVectorZConst ( vectorZ );
    }
}


//-----------------------------------------------------------------------------------------------------------------------------------------
void saveImageAsVTKCell ( AdaptiveFEPrismMesh < RealType > & grid, aol::Vector < RealType > & vector, const std::string savename, const int run = -1 ) {
    stringstream filename;
    if ( run == -1 ) {
        filename << "Results/" << savename << ".vtu"; 
    }
    else {
        filename << "Results/" << savename << "_" << run << ".vtu"; 
    }
    std::string filenameStr = filename.str();
    grid.saveAsLegacyVTKCell ( filenameStr, vector );   
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// Function that creates a matrix P that interpolates a vector corresponding to an adaptive grid to a vector of another adaptive grid.
// It goes through all elements of all levels and checks if a node of the new grid exists in the old grid -> takes the old value.
// If a node is new it is interpolated by its neighboured nodes which are already set with values. This is hard coded for all possible positions of nodes.
template < typename ConfiguratorType >
void makeInterpolationMatrixOntoAdaptiveGridZConstant ( SpecialAdaptiveFEPrismMeshRPOp < ConfiguratorType > &rp, SpecialAdaptiveFEPrismMeshRPOp < ConfiguratorType > &newRp, typename ConfiguratorType::MatrixType &P ) {
  typedef SpecialAdaptiveFEPrismMeshRPOp < ConfiguratorType > RPOpType;
  // Go through all elements in new grid
  aol::BitVector visited ( newRp.getGridRef().getNumberOfDofs() );
  visited.setAll ( false );
  for ( typename RPOpType::GridType::ElementIterator it ( newRp.getGridRef() ); it.notAtEnd(); ++it ) { 
    // Go through all nodes in current element 
    aol::Vec < 6, GlobalIndex > nodeInds =  it.getNodeIndices();  
    for ( int i = 0; i < 6; ++i ) {
      // Check if node is non-hanging and was not visited 
      if ( !newRp.getGridRef().isHangingNode ( nodeInds[i] ) && !visited.get ( newRp.getRestrictedGlobalNodeIndex ( nodeInds[i] ) ) ) {  
        // Check if node existed before in old rp's indexMap 
        auto got = rp.getIndexMap().find ( nodeInds[i] ); 
        // Node existed before (might have been hanging before)
        if ( got != rp.getIndexMap().end() ) {
          // Node was hanging before, interpolate from old interpolation nodes 
          P.set ( newRp.getRestrictedGlobalNodeIndex ( nodeInds[i] ), got->second._node, got->second._prolongationWeight ); 
        }  
        // Node did not exist before: Interpolate from grid interpolation map  
        else {
          auto interpol = newRp.getGridRef().getInterpolationMap().find ( nodeInds[i] );
          P.set ( newRp.getRestrictedGlobalNodeIndex ( nodeInds[i] ), interpol->second[0], 0.5 );  
          P.set ( newRp.getRestrictedGlobalNodeIndex ( nodeInds[i] ), interpol->second[1], 0.5 );  
        }  
        // Mark node as visited 
        visited.set ( newRp.getRestrictedGlobalNodeIndex ( nodeInds[i] ), true );  
      }
    }
  }        
}  

#endif