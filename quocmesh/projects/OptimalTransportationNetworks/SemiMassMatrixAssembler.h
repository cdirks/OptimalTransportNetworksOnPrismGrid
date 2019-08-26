#ifndef __SEMIMASSMATRIXASSEMBLER_H
#define __SEMIMASSMATRIXASSEMBLER_H


#include<aol.h>
#include<quoc.h>
#include<iostream>
#include<scalarArray.h>
#include<FEOpInterface.h>
#include<configurators.h>
#include<Newton.h>
#include<omp.h>

//------------------------------------------------------------------------------------------------------------------------------------------
// Class for matrix assembling: L^1_ij = \int Basis_i*gradx(Basis_j) dxds
// Derives from FELinScalarWeightedSemiDiffInterface -> new class for interface!
template < typename ConfiguratorType >
class SemiMassMatrixAssemblerX : public aol::FELinScalarWeightedSemiDiffInterface < ConfiguratorType, SemiMassMatrixAssemblerX < ConfiguratorType > > 
{
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  GridType & _grid;
  
  public:
  SemiMassMatrixAssemblerX ( ConfiguratorType & Configurator, GridType & Grid ) 
    : aol::FELinScalarWeightedSemiDiffInterface < ConfiguratorType, SemiMassMatrixAssemblerX < ConfiguratorType > > ( Configurator, Grid, 0, false, aol::ASSEMBLED ), _grid(Grid) { }
  
  // Define coefficients w(x) = 1
  RealType getCoeff(const typename ConfiguratorType::ElementType& /*El*/, int /*QuadPoint*/, const typename ConfiguratorType::DomVecType & /*RefCoord*/) const {
    return 1.0;
  }
  
};
//------------------------------------------------------------------------------------------------------------------------------------------

// Class for matrix assembling: L^2_ij = \int Basis_i*grady(Basis_j) dxds
// Derives from FELinScalarWeightedSemiDiffInterface -> new class for interface!
template < typename ConfiguratorType >
class SemiMassMatrixAssemblerY : public aol::FELinScalarWeightedSemiDiffInterface < ConfiguratorType, SemiMassMatrixAssemblerY < ConfiguratorType > > 
{
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  GridType & _grid;
  
  public:
  SemiMassMatrixAssemblerY ( ConfiguratorType & Configurator, GridType & Grid )
    : aol::FELinScalarWeightedSemiDiffInterface< ConfiguratorType, SemiMassMatrixAssemblerY < ConfiguratorType > > ( Configurator, Grid, 1, false, aol::ASSEMBLED ), _grid(Grid) { }  

  // Define coefficients w(x) = 1
  RealType getCoeff(const typename ConfiguratorType::ElementType& /*El*/, int /*QuadPoint*/, const typename ConfiguratorType::DomVecType & /*RefCoord*/) const {
    return 1.0;
  }
};
//------------------------------------------------------------------------------------------------------------------------------------------

// Class for matrix assembling: L^2_ij = \int Basis_i*gradz(Basis_j) dxds
// Derives from FELinScalarWeightedSemiDiffInterface -> new class for interface!
template < typename ConfiguratorType >
class SemiMassMatrixAssemblerZ : public aol::FELinScalarWeightedSemiDiffInterface < ConfiguratorType, SemiMassMatrixAssemblerZ < ConfiguratorType > > 
{
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  GridType & _grid;
  
  public:
  SemiMassMatrixAssemblerZ ( ConfiguratorType & Configurator, GridType & Grid )
    : aol::FELinScalarWeightedSemiDiffInterface< ConfiguratorType, SemiMassMatrixAssemblerZ < ConfiguratorType > > ( Configurator, Grid, 2, false, aol::ASSEMBLED), _grid(Grid) { }  

  // Define coefficients w(x) = 1
  RealType getCoeff(const typename ConfiguratorType::ElementType& /*El*/, int /*QuadPoint*/, const typename ConfiguratorType::DomVecType & /*RefCoord*/) const {
    return 1.0;
  }
};
//------------------------------------------------------------------------------------------------------------------------------------------


#endif

