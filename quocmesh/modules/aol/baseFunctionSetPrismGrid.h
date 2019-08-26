#ifndef __BASEFUNCTIONSETPRISMGRID_H
#define __BASEFUNCTIONSETPRISMGRID_H

#include <aol.h>
#include <quoc.h>
#include <vec.h>
#include <smallMat.h>

template < typename RealType, typename QuadRuleType, typename PrismType >
class BaseFunctionSetPrismGrid  {
    
  typedef aol::Vec3 < RealType > DomVecType;   
  
  static RealType _b1    ( const aol::Vec3<RealType> &RefCoord ) { return 1.0 - RefCoord[0] - RefCoord[1] - RefCoord[2] + RefCoord[0]*RefCoord[2] + RefCoord[1]*RefCoord[2]; }
  static RealType _b2    ( const aol::Vec3<RealType> &RefCoord ) { return RefCoord[0] - RefCoord[0]*RefCoord[2]; }
  static RealType _b3    ( const aol::Vec3<RealType> &RefCoord ) { return RefCoord[1] - RefCoord[1]*RefCoord[2]; }
  static RealType _b4    ( const aol::Vec3<RealType> &RefCoord ) { return RefCoord[2] - RefCoord[0]*RefCoord[2] - RefCoord[1]*RefCoord[2]; }
  static RealType _b5    ( const aol::Vec3<RealType> &RefCoord ) { return RefCoord[0]*RefCoord[2]; }
  static RealType _b6    ( const aol::Vec3<RealType> &RefCoord ) { return RefCoord[1]*RefCoord[2]; }
  
  static RealType _dx_b1 ( const aol::Vec3<RealType> &RefCoord ) { return -1.0 + RefCoord[2]; }
  static RealType _dx_b2 ( const aol::Vec3<RealType> &RefCoord ) { return 1.0 - RefCoord[2]; }
  static RealType _dx_b3 ( const aol::Vec3<RealType> &/*RefCoord*/ ) { return 0.0; }
  static RealType _dx_b4 ( const aol::Vec3<RealType> &RefCoord ) { return -RefCoord[2]; }
  static RealType _dx_b5 ( const aol::Vec3<RealType> &RefCoord ) { return RefCoord[2]; }
  static RealType _dx_b6 ( const aol::Vec3<RealType> &/*RefCoord*/ ) { return 0.0; }

  static RealType _dy_b1 ( const aol::Vec3<RealType> &RefCoord ) { return -1.0 + RefCoord[2]; }
  static RealType _dy_b2 ( const aol::Vec3<RealType> &/*RefCoord*/ ) { return 0.0; }
  static RealType _dy_b3 ( const aol::Vec3<RealType> &RefCoord ) { return 1.0 - RefCoord[2]; }
  static RealType _dy_b4 ( const aol::Vec3<RealType> &RefCoord ) { return -RefCoord[2]; }
  static RealType _dy_b5 ( const aol::Vec3<RealType> &/*RefCoord*/ ) { return 0.0; }
  static RealType _dy_b6 ( const aol::Vec3<RealType> &RefCoord ) { return RefCoord[2]; }

  static RealType _dz_b1 ( const aol::Vec3<RealType> &RefCoord ) { return -1.0 + RefCoord[0] + RefCoord[1]; }
  static RealType _dz_b2 ( const aol::Vec3<RealType> &RefCoord ) { return -RefCoord[0]; }
  static RealType _dz_b3 ( const aol::Vec3<RealType> &RefCoord ) { return -RefCoord[1]; }
  static RealType _dz_b4 ( const aol::Vec3<RealType> &RefCoord ) { return 1.0 - RefCoord[0] - RefCoord[1]; }
  static RealType _dz_b5 ( const aol::Vec3<RealType> &RefCoord ) { return RefCoord[0]; }
  static RealType _dz_b6 ( const aol::Vec3<RealType> &RefCoord ) { return RefCoord[1]; }
  
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const aol::Vec3 < RealType > & RefCoord );
  BASIS_FUNC_TYPE _deriv_basis[3][6];
  BASIS_FUNC_TYPE _basis[6];
  
  const PrismType *_prism;
  RealType _h;
  
  
public:
    
  BaseFunctionSetPrismGrid () : _prism ( NULL ), _h ( 1.0 ) {
      
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;
    _deriv_basis[0][4] = _dx_b5;
    _deriv_basis[0][5] = _dx_b6;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    _deriv_basis[1][4] = _dy_b5;
    _deriv_basis[1][5] = _dy_b6;

    _deriv_basis[2][0] = _dz_b1;
    _deriv_basis[2][1] = _dz_b2;
    _deriv_basis[2][2] = _dz_b3;
    _deriv_basis[2][3] = _dz_b4;
    _deriv_basis[2][4] = _dz_b5;
    _deriv_basis[2][5] = _dz_b6; 
      
  }
  
  BaseFunctionSetPrismGrid ( const RealType H ) : _prism ( NULL ), _h ( H ) {
      
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;
    _deriv_basis[0][4] = _dx_b5;
    _deriv_basis[0][5] = _dx_b6;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    _deriv_basis[1][4] = _dy_b5;
    _deriv_basis[1][5] = _dy_b6;

    _deriv_basis[2][0] = _dz_b1;
    _deriv_basis[2][1] = _dz_b2;
    _deriv_basis[2][2] = _dz_b3;
    _deriv_basis[2][3] = _dz_b4;
    _deriv_basis[2][4] = _dz_b5;
    _deriv_basis[2][5] = _dz_b6;

  }
  
  enum { numBaseFuncs = 6 };
  
  void setPrism ( const PrismType &P ) {
    _prism = &P;
    _h = P.getHeight();
  }
  
  int numQuadPoints () const {
    return QuadRuleType::numQuadPoints;
  }
    
  inline RealType getWeight ( int QuadPoint ) const {
    return _quadRule.getWeight ( QuadPoint );
  }
  
  inline const DomVecType & getRefCoord ( int QuadPoint ) const {
    return _quadRule.getRefCoord ( QuadPoint );
  }
 
  RealType evaluate ( int BaseFuncNum, const aol::Vec3 < RealType > & RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }
 
  void evaluateGradient ( int BaseFuncNum, const aol::Vec3 < RealType > & RefCoord, aol::Vec3 < RealType > & Gradient ) const { 
    // Initialize vectors
    aol::Vec2 < RealType > tmp, tmp2;  
    // Gradient at quad point in barycentric coords
    tmp[0] = _deriv_basis[0][BaseFuncNum] ( RefCoord );
    tmp[1] = _deriv_basis[1][BaseFuncNum] ( RefCoord );  
    // Change coordinate system 
    _prism->ginv().mult ( tmp, tmp2 );
    // Compute tangent vectors
    const aol::Vec3 < RealType > dir0 = _prism->edge(0,1);
    const aol::Vec3 < RealType > dir1 = _prism->edge(0,2);
    for ( int i = 0; i < 2; i++ ) {
      Gradient[i] = tmp2[0] * dir0[i] + tmp2[1] * dir1[i];
    }
    // Set z gradient 
    Gradient[2] = _deriv_basis[2][BaseFuncNum] ( RefCoord ) / _h;
  }

  inline aol::Vec3 < RealType > evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    aol::Vec3 < RealType > g;
    this->evaluateGradient ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), g );
    return g;
  }
  

  
protected: 
  
  QuadRuleType _quadRule;

};

#endif