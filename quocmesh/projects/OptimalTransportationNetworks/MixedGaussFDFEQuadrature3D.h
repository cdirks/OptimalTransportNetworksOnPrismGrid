#ifndef __MIXEDGAUSSFDQUADRATURE3DNONUNIFORM_H
#define __MIXEDGAUSSFDQUADRATURE3DNONUNIFORM_H

#include<aol.h>


// Gauss quadrature in x/y direction and FD quadrature in z direction on prism element
template < typename _RealType >
class MixedGaussFDFEQuadrature3D {
  
public:
  enum { numQuadPoints = 3 };
  
  MixedGaussFDFEQuadrature3D() {
    
    _points[0][0] = 0.1666666666666667;
    _points[0][1] = 0.1666666666666667;
    _points[0][2] = 0.0;
    
    _points[1][0] = 0.6666666666666667;
    _points[1][1] = 0.1666666666666667;
    _points[1][2] = 0.0;
    
    _points[2][0] = 0.1666666666666667;
    _points[2][1] = 0.6666666666666667;
    _points[2][2] = 0.0;
    
    _weights[0] = 0.3333333333333333;
    _weights[1] = 0.3333333333333333;
    _weights[2] = 0.3333333333333333;
    
  }
  
  inline const aol::Vec3 < _RealType > &getRefCoord ( int QuadPoint ) const {
    return _points[QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    return _weights[QuadPoint];
  }

protected:
  _RealType _weights[numQuadPoints];
  aol::Vec3 < _RealType > _points[numQuadPoints];
  
};


#endif