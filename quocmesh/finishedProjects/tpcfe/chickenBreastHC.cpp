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

#include <indexMapper.h>
#include <parameterParser.h>
#include <scalarArray.h>

#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>


template< typename DataType >
void generateEllipsoidLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, aol::Vec3<DataType> radii ) {
  const int width = levelset.getNumX();

  const DataType ctr_x  = 0.5, ctr_y  = 0.5, ctr_z  = 0.5;

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {

    const DataType x = ( 1.0 * (*bit)[0] ) / ( width - 1 ), y = ( 1.0 * (*bit)[1] ) / ( width - 1 ), z = ( 1.0 * (*bit)[2] ) / ( width - 1 );
    const DataType value = sqrt ( aol::Sqr ( x - ctr_x ) / aol::Sqr ( radii[0] ) + aol::Sqr ( y - ctr_y ) / aol::Sqr ( radii[1] )  + aol::Sqr ( z - ctr_z ) / aol::Sqr ( radii[2] ) )  - 1 ;

    levelset.set ( *bit, value );
  }

}


template< typename GridType, typename ArrayType >
void writeTemperatureSlicePLYs ( const GridType &grid, const ArrayType &temperature, const char* fnmask ) {
  // cerr << "Saving PLYs to " << fnmask;
  typedef typename GridType::RealType RealType;
  // const int numXYZ = grid.getNumXYZ();
  for ( short diri = 0; diri < 3; ++diri ) {
    const qc::Comp dir = static_cast<qc::Comp>( diri );
    // for ( short pos = 0; pos < grid.getNumXYZ(); ++ pos ) {
    // short pos = grid.getNumXYZ() / 2; {
    for ( short pos = grid.getNumXYZ()/2 - 2; pos < grid.getNumXYZ()/2 + 3; ++ pos ) {

      tpcfe::CFEInterfaceTriangulationWithScalarDataGenerator< GridType > itg ( grid, temperature );
      itg.determineSliceTriangulation ( dir, pos );
      // not bad so far, but want more than the slice intersected with the object

      char filename[1024];
      sprintf ( filename, fnmask, dir, pos );
      itg.saveToPLYFile ( filename );
    }
  }
}


typedef double RealType;

static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;
typedef tpcfe::CFEGrid < RealType, AT >                           GridType;
typedef tpcfe::CFEHybridMatrix<GridType>                        MatrixType;

typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
typedef tpcfe::CFEMassOp  < ConfiguratorType >                  MassOpType;
typedef tpcfe::CFEStiffOpWI < ConfiguratorType >              WStiffOpType;
typedef tpcfe::CFEMassOpWI < ConfiguratorType >                WMassOpType;

int main ( int argc, char** argv ) {
  try {

    tpcfe::CFEVirtualNode<RealType,AT,RealType>::_maximalInversionError = 1e-15;

    aol::ParameterParser params ( argc < 2 ? "par/chickenBreastHC.par" : argv[1] );

    params.dump();

    // "numerics parameters"
    const int
      level         = params.getInt ( "level" ),
      numTimesteps  = params.getInt ( "numTimesteps" ),
      savePlyEvery  = params.getInt ( "savePlyEvery" ),
      explicitLevel = 3,
      nSmooth       = 3;

    bool useLumpedMassMat = ( params.getInt ( "useLumpedMassMat" ) != 0 );

    aol::Vector<RealType> MeatTempTime ( numTimesteps  + 1 ), WaterTempTime ( numTimesteps + 1 );

    // meat: negative = interior
    GridType grid ( level );

    qc::ScalarArray<RealType, qc::QC_3D> dataset ( grid );

    generateEllipsoidLevelset ( dataset, aol::Vec3<RealType> ( 1.25/12.5, 3.25/12.5, 5.00/12.5 ) ); // with sigma = 0.25 m, this corresponds to an ellipse with diameters 10, 6.5, 2.5 cm

    grid.addStructureFrom ( dataset );

    // "material parameters"
    const RealType
      sigma         = params.getDouble ("sigma"),                          // size of bounding box; in m
      volume_Meat   = aol::Cub ( sigma ) * grid.getTotalInnerVolume(),
      volume_RWater = params.getDouble ( "volume_RWater" ),                // volume of real water; in m^3
      volume_VWater = aol::Cub ( sigma ) - volume_Meat,
      lambda_Meat   = params.getDouble ( "lambda_Meat" ),                  // heat conductivity of meat; this is the actual unknown to be determined from time scaling; in W K^{-1} m^{-1}
      lambda_VWater = params.getDouble ( "lambda_VWater" ),                // virtual heat conductivity of virtual water modelling continuous stirring; in W K^{-1} m^{-1}
      rhoC_Meat     = params.getDouble ( "rhoC_Meat" ),                    // heat capacity times density of meat, determined by experiment; in J K^{-1} m^{-3}
      rhoC_RWater   = params.getDouble ( "rhoC_RWater" ),                  // heat capacity times density of real water; in J K^{-1} m^{-3}
      rhoC_VWater   = ( volume_RWater / volume_VWater ) * rhoC_RWater,     // rhoC for virtual water so that both volumes have same total capacity; in J K^{-1} m^{-2}
      temperatureTolerance = params.getDouble ( "temperatureTolerance" );  // stop time stepping when relative temperature difference is smaller than this threshold

    cerr << "volume_Meat   = " << volume_Meat << endl
         << "volume_RWater = " << volume_RWater << endl
         << "volume_VWater = " << volume_VWater << endl
         << "lambda_Meat   = " << lambda_Meat << endl
         << "lambda_VWater = " << lambda_VWater << endl
         << "rhoC_Meat     = " << rhoC_Meat << endl
         << "rhoC_VWater   = " << rhoC_VWater << endl;

    const qc::CoordType
      waterSensorPos ( 1, 1, 1 ),
      meatSensorPos ( grid.getNumX() / 2, grid.getNumY() / 2, grid.getNumZ() / 2 );


    //     const short numSensors = 6;
    //     aol::MultiVector<RealType> MeatMultiTempTime ( numSensors, numTimesteps + 1 );
    //     std::vector< qc::CoordType > meatMultiSensorPos ( numSensors );
    //     const int xctr = grid.getNumX() / 2, yctr = xctr, nz = grid.getNumZ();
    //     meatMultiSensorPos[0] = qc::CoordType ( xctr, yctr, ( 60 * nz ) / 100 );
    //     meatMultiSensorPos[1] = qc::CoordType ( xctr, yctr, ( 65 * nz ) / 100 );
    //     meatMultiSensorPos[2] = qc::CoordType ( xctr, yctr, ( 70 * nz ) / 100 );
    //     meatMultiSensorPos[3] = qc::CoordType ( xctr, yctr, ( 75 * nz ) / 100 );
    //     meatMultiSensorPos[4] = qc::CoordType ( xctr, yctr, ( 80 * nz ) / 100 );
    //     meatMultiSensorPos[5] = qc::CoordType ( xctr, yctr, ( 85 * nz ) / 100 );


    // "initial conditions"
    const RealType
      initialTemperatureWater = params.getDouble ( "initialTemperatureWater" ),  // in K
      initialTemperatureMeat  = params.getDouble ( "initialTemperatureMeat" ),  // in K
      tau                     = params.getDouble ( "tau" ); // in s


    // set up CFE
    qc::AArray<RealType, qc::QC_3D> coeffLambda ( grid );
    for ( int i = 0; i < dataset.size(); ++i )
      coeffLambda[i] = ( dataset[i] < 0 ? lambda_Meat : lambda_VWater );
    grid.detectAndInitVirtualNodes ( coeffLambda );

    // set up yet another correction term: initial temperature profile is discontinuous and under-represented by CFE interpolation
    RealType initTempCorrFactor = 0;
    {
      MassOpType unweightedMassOp ( grid, aol::ASSEMBLED );
      qc::ScalarArray<RealType, qc::QC_3D> chiMeat ( grid ), dwarf ( grid );
      for ( int i = 0; i < chiMeat.size(); ++i )
        if ( dataset[i] < 0 )
          chiMeat[i] = 1;
      unweightedMassOp.apply ( chiMeat, dwarf );
      initTempCorrFactor = volume_Meat / ( aol::Cub ( sigma ) * dwarf.sum() );
      cerr << initTempCorrFactor << endl;
    }

    qc::AArray<RealType, qc::QC_3D> coeffRhoC ( grid );
    for ( int i = 0; i < dataset.size(); ++i )
      coeffRhoC[i] = ( dataset[i] < 0 ? initTempCorrFactor * rhoC_Meat : rhoC_VWater );


    aol::DiagonalMatrix<RealType> massMatDiag ( grid ); // this will become the lumped mass matrix
    MatrixType stiffMat ( grid ), sysMat ( grid ), massMatNondiag ( grid ); // full mass matrix, not merely diagonal (lumped) one
    {
      WMassOpType massOp ( coeffRhoC, grid, aol::ASSEMBLED );
      WStiffOpType stiffOp ( coeffLambda, grid, aol::ASSEMBLED );
      massOp.assembleAddMatrix ( massMatNondiag );

      massMatNondiag *= aol::Cub ( sigma );

      // poor man's mass lumping:
      massMatDiag.setToLumpedCopyOf ( massMatNondiag );

      stiffOp.assembleAddMatrix ( stiffMat );
      stiffMat *= sigma;

      sysMat += stiffMat;
      sysMat *= tau;
      if ( useLumpedMassMat ) {
        sysMat += massMatDiag;
      } else {
        sysMat += massMatNondiag;
      }
    }


    qc::ScalarArray<RealType, qc::QC_3D> temperature ( grid ), rhs ( grid );


#if 1

    tpcfe::CFEMultigrid< /*misuse: */ tpcfe::CFEStiffOp<ConfiguratorType>, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > heatSolver ( grid, sysMat, explicitLevel, nSmooth, nSmooth, aol::GAUSS_SEIDEL_FORWARD, 1.0,
                                                                                                                                    mg::MULTIGRID_V_CYCLE, 1.0e-16, 500, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
#else
    { RealType dummy = explicitLevel; dummy *= nSmooth; } // prevent unused variable warnings
    aol::DiagonalPreconditioner<RealType> prec ( sysMat );
    aol::PCGInverse< aol::Vector<RealType> > heatSolver ( sysMat, prec, 1.0e-16, 2000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
#endif

    // initial values for temperature
    temperature.setAll ( initialTemperatureWater );

    for ( qc::RectangularIterator<qc::QC_3D> pit ( grid ); pit.notAtEnd(); ++pit ) {
      if ( dataset.get ( *pit ) < 0 ) {
        temperature.set ( *pit, initialTemperatureMeat );
      }
    }

    MeatTempTime[0]  = temperature.get ( meatSensorPos );
    WaterTempTime[0] = temperature.get ( waterSensorPos );
    cerr << aol::shortFormat ( tau * 0 ) << " " << aol::longScientificFormat ( MeatTempTime[0] ) << " " << aol::longScientificFormat ( WaterTempTime[0] ) << endl;

    //     for ( short i = 0; i < numSensors; ++i ) {
    //       MeatMultiTempTime[i][0] = temperature.get( meatMultiSensorPos[i] );
    //     }

    char filenamemask[1024];
    sprintf ( filenamemask, "out/chickenBreastTemperature_%03d_%%d_%%03d.ply.bz2", 0 );
    writeTemperatureSlicePLYs ( grid, temperature, filenamemask );

    int numTimestepsActuallyPerformed = 0;
    // THE loop
    for ( short t = 0; t < numTimesteps; ++t, ++numTimestepsActuallyPerformed ) {

      if ( useLumpedMassMat ) {
        massMatDiag.apply ( temperature, rhs );
      } else {
        massMatNondiag.apply ( temperature, rhs );
      }

      temperature.setZero();
      heatSolver.apply ( rhs, temperature );

      MeatTempTime[t+1]  = temperature.get ( meatSensorPos );
      WaterTempTime[t+1] = temperature.get ( waterSensorPos );
      cerr << aol::shortFormat ( tau * (t+1) ) << " " << aol::longScientificFormat ( MeatTempTime[t+1] ) << " " << aol::longScientificFormat ( WaterTempTime[t+1] ) << endl;

      if ( (t+1) % savePlyEvery == 0 ) {
        sprintf ( filenamemask, "out/chickenBreastTemperature_%03d_%%d_%%03d.ply.bz2", t+1 );
        writeTemperatureSlicePLYs ( grid, temperature, filenamemask );
      }

      //       for ( short i = 0; i < numSensors; ++i ) {
      //         MeatMultiTempTime[i][t+1] = temperature.get( meatMultiSensorPos[i] );
      //       }

      if ( fabs ( MeatTempTime[t+1] - WaterTempTime[t+1] ) < temperatureTolerance * fabs ( initialTemperatureMeat - initialTemperatureWater ) )
        break;

      {
        aol::Vector<RealType> dummy ( temperature, aol::STRUCT_COPY );
        if ( useLumpedMassMat ) {
          massMatDiag.apply ( temperature, dummy );
        } else {
          massMatNondiag.apply ( temperature, dummy );
        }
        cerr << "Total Energy: " << dummy.sum() << endl;
      }
    }

    // output data for displaying:
    dataset.save ( "out/chickenBreastLevelset.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

    // ofstream avgTempOut ( "out/avgTemp.dat" );
    cout << endl;
    for ( int t = 0; t <= numTimestepsActuallyPerformed; ++t ) {
      cout << aol::detailedFormat ( tau * t ) << " " << aol::longScientificFormat ( MeatTempTime[t] ) << " " << aol::longScientificFormat ( WaterTempTime[t] ) << endl;
    }

    //     for ( short i = 0; i < numSensors; ++i ) {
    //       cout << endl << "Outputting additional sensor data no " << i << endl;
    //       for ( int t = 0; t <= numTimestepsActuallyPerformed; ++t ) {
    //         cout << aol::detailedFormat ( tau * t ) << " " << aol::longScientificFormat ( MeatMultiTempTime[i][t] ) << " " << aol::longScientificFormat ( aol::NumberTrait<RealType>::NaN ) << endl;
    //       }
    //       cout << endl;
    //     }

    // parameter fitting and tau rescaling in separate program.


  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}

#if 0
template< typename DataType >
class ColorLookupTable {
protected:
  std::vector< aol::Vec3<DataType> > _discreteColors;
  DataType _min, _max;

public:
  ColorLookupTable ( ) : _discreteColors(), _min ( 0 ), _max ( 0 ) {
  }

  void pushBackColor ( const aol::Vec3<DataType> &newColor ) {
    _discreteColors.push_back ( newColor );
  }

  void setMin ( const DataType min ) {
    _min = min;
  }

  void setMax ( const DataType max ) {
    _max = max;
  }

  aol::Vec3<DataType> getColor ( const DataType value ) const {
    if ( ( value < _min ) || ( value > _max ) ) {
      cerr << value << " " << _min << " " << _max << endl;
      throw aol::Exception("ColorLookupTable: value out of range [_min,_max]", __FILE__, __LINE__ );
    }

    const DataType bary = ( value - _min ) / ( _max - _min );
    const DataType baryS = bary * ( _discreteColors.size() - 1 );

    if ( static_cast<unsigned int>( baryS ) == _discreteColors.size() - 1 ) {
      return ( _discreteColors[ _discreteColors.size() - 1 ] );
    } else {
      const unsigned int colA = static_cast<unsigned int> ( floor ( baryS ) ), colB = colA + 1;
      const DataType weiB = baryS - colA, weiA = 1 - weiB;
      return ( weiA * _discreteColors[ colA ] + weiB * _discreteColors[ colB ] );
    }
  }

};
#endif

#if 0
template< typename RealType >
void writeColorImage ( const qc::ScalarArray<RealType, qc::QC_3D> &values, const char* filenamemask, ColorLookupTable<RealType> &CLT ) {
  for ( int z = 0; z < values.getNumZ(); ++z ) {
    qc::MultiArray<RealType,2,3> imageData ( values.getNumX(), values.getNumY() );
    for ( short y = 0; y < values.getNumY(); ++y ) {
      for ( short x = 0; x < values.getNumX(); ++x ) {
        imageData.set ( aol::Vec2<short>(x,y), 255.0 * CLT.getColor ( values.get ( x, y, z ) ) );
      }
    }
    char filename[1024];
    sprintf ( filename, filenamemask, z );
    imageData.savePNG ( filename );
  }
}
#endif

#if 0
template< typename GridType, typename RealType >
RealType integrateFunctionOverPart ( const GridType& grid, const qc::ScalarArray<RealType, qc::QC_3D> &values, const short which ) {
  RealType integral = 0;
  const RealType h = grid.H();

  tpcfe::CFEElementIterator<RealType> it ( grid );
  for ( it.start(); ! ( it.atEnd() ); ++it ) {
    const tpcfe::CFEElement<RealType> &element = *it;
    const int cfetype = element.pureCFEType();
    element.computeAssembleData ( grid );

    tpcfe::CFETetrahedronIterator<RealType> tit;
    for ( tit.start ( cfetype, which ); ! ( tit.atEnd() ); ++tit ) {
      const tpcfe::CFETetra<RealType> &tetra = * ( *tit );
      aol::Vec<4,RealType> vertexValues;
      for ( short i = 0; i < 4; ++i ) {
        if ( tetra.isVirtualNode(i) ) {
          vertexValues[i] = grid.getVirtualNodeRef ( element.getGlobalIndex( tetra(i,0) ), element.getGlobalIndex( tetra(i,1) ) ).extrapolate ( values );
        } else {
          vertexValues[i] = values.get ( element.getGlobalIndex( tetra.getNode ( i ) ) );
        }
      }
      integral += tetra.getVolume() * aol::Cub ( h ) * vertexValues.sum() / 4;
    }
  }

  return ( integral );
}
#endif
