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

/** \file
 *  \brief generates and modifies 3D level arrays
 *
 *  The program is called with one or more parameter files as console argument(s).
 *  Each parameter file has to specify an "action" from the following list:
 *
 *     setConst
 *     addConst
 *     linComb
 *     maxOfTwo
 *     minOfTwo
 *     product
 *     smooth
 *     move
 *     reflect
 *     reorderComps
 *     redistance
 *     cleanLevelFct
 *     generateLevelFunction3D
 *
 *  and an "output" file. Depending on the selected action, further parameters are
 *  read from the parameter file:
 *
 *  (a) setConst
 *
 *      requires "value" and "gridLevel" and produces a constant image with this value
 *
 *  (b) addConst
 *
 *      requires "input" and "value" and adds the value to the given image
 *
 *  (c) linComb, maxOfTwo, minOfTwo, product
 *
 *      They all require "input1" and "input2", linComb additionally searches the
 *      parameter file for "coeff1" and "coeff2". They write into the output file
 *      what you would expect from their names.
 *
 *  (d) smooth
 *
 *      needs "input" data array and smoothing parameter "sigma". With higher sigma,
 *      smoothing is stronger.
 *
 *  (e) move
 *
 *      needs "input", "compToMove" and "moveAmount". The compToMove'th coordinate is
 *      moved by moveAmound \in [0;1], at the other end of the 3D array, the last
 *      slice is repeatedly added to fill the emerging empty part of the array.
 *
 *  (f) reflect
 *
 *      reflects the "input" array's entries in the "compToReflect" 'th component.
 *
 *  (g) reorderComps
 *
 *      reorders the components. Reads a int vector "newOrder" of length 3, say: p,
 *      and re-arranges the array such that
 *
 *             output( a ) = input ( b )       where b_i = a_p(i).
 *
 *  (e) redistance
 *
 *      writes a level array to "output" that has the same zero-level-set as
 *      "input", but is a discretized signed distance function
 *
 *  (i) cleanLevelFct
 *
 *      calls the class SmallComponentsPoempelRemover which removes all
 *      connectivity components of the sub-zero level set but the largest
 *      one. This is done by finding the largest sub-zero level set of "input"
 *      setting negative level values to "newValue" everywhere else.
 *
 *  (j) generateLevelFunction3D
 *
 *      uses a qc::DataGenerator object to produce level arrays. The general type of
 *      level set is given as "type". Nearly all types need "center" and "radius", so
 *      these parameters are not mentioned explicitely further on. All types will
 *      generate a ScalarArray of grid depth "gridLevel" and save is as "output".
 *
 *      (1) sphere
 *          creates a sphere with given center and radius
 *
 *      (2) l1sphere
 *          creates an l_1 metric sphere (pyramid) with given center and radius
 *
 *      (3) linftysphere
 *          creates an l_infty metric sphere (cube) with given center and radius
 *
 *      (4) ellipsoid
 *          takes additional parameter "scaling". With center p, radius r and scaling s,
 *          the zero level set is a discretization of
 *
 *                 \f$ \{ x | \sqrt (\sum [s_i x_i - c_i]^2 ) < r \} \f$
 *
 *      (5) plane
 *          creates a plane with given "normal" and "point" lying on it
 *          (needs no "center" nor "radius").
 *
 *      (6) cylinder
 *          needs "center", "scaling" and "radius".
 *
 *  \author Stefan W. von Deylen
 */

#include <generator.h>
#include <configurators.h>
#include <parameterParser.h>
#include <linearSmoothOp.h>
#include <signedDistanceOp.h>
#include <tpCFELevelsets.h>

using namespace aol;
using namespace qc;

#include <cstdlib>
#include <iostream>

using namespace std;

//------------------------------------------------------------------------------------------------

typedef long double                                          RealType;

typedef QuocConfiguratorTraitMultiLin<RealType,
                                      QC_3D,
                                      GaussQuadrature
                                      <RealType, QC_3D, 3> > ConfType;

typedef FastUniformGridMatrix<RealType, ConfType::Dim>       MatrixType;
typedef ConfType::ArrayType                                  ArrayType;

const SaveType saveType = PGM_DOUBLE_BINARY;

//------------------------------------------------------------------------------------------------

#include "dataGenerator2D3D.h"

//------------------------------------------------------------------------------------------------

void generateLevelFunction3D ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  int level = parser.getInt ( "gridLevel" );
  ConfType::InitType grid ( level, QC_3D );
  DataGenerator<ConfType> generator ( grid );
  ArrayType levelValues ( grid );

  string type = parser.getString ( "type" );
  if ( type == "sphere" ) {
    Vector<RealType> centerVector;
    parser.getRealVec ( "center", centerVector );
    Vec3<RealType> center ( centerVector[0], centerVector[1], centerVector[2] );
    RealType radius = parser.getDouble ( "radius" );
    clog << "generating sphere with center " << centerVector << endl
         << "and radius " << radius << endl;
    generator.generateSphereLevelset ( levelValues, center, radius );
  }
  else if ( type == "l1sphere" ) {
    Vector<RealType> centerVector;
    parser.getRealVec ( "center", centerVector );
    Vec3<RealType> center ( centerVector[0], centerVector[1], centerVector[2] );
    RealType radius = parser.getDouble ( "radius" );
    clog << "generating l1 sphere with center " << centerVector << endl
         << "and radius " << radius << endl;
    generator.generateL1SphereLevelset ( levelValues, center, radius );
  }
  else if ( type == "linfinitysphere" ) {
    Vector<RealType> centerVector;
    parser.getRealVec ( "center", centerVector );
    Vec3<RealType> center ( centerVector[0], centerVector[1], centerVector[2] );
    RealType radius = parser.getDouble ( "radius" );
    clog << "generating l_infty sphere with center " << centerVector << " and radius " << radius << endl;
    generator.generateLInfinitySphereLevelset ( levelValues, center, radius );
  }
  else if ( type == "ellipsoid" ) {
    Vector<RealType> centerVector, scalingVector;
    parser.getRealVec ( "center", centerVector );
    parser.getRealVec ( "scaling", scalingVector );
    Vec3<RealType> center ( centerVector.getData() ), scaling ( scalingVector.getData() );
    RealType radius = parser.getDouble ( "radius" );
    clog << "generating ellipsoid with center = " << centerVector << "," << endl
          << "scaling = " << scalingVector << " and radius " << radius << endl;
    generator.generateEllipsoidLevelset ( levelValues, center, scaling, radius );
  }
  else if ( type == "plane" ) {
    Vector<RealType> normalVector, pointVector;
    parser.getRealVec ( "normal", normalVector );
    parser.getRealVec ( "point", pointVector );
    Vec3<RealType> normal ( normalVector.getData() ), point  ( pointVector.getData() );
    clog << "generating plane with normal " << normalVector << endl
         << "and point " << pointVector << endl;
    generator.generatePlaneLevelset ( levelValues, normal, point );
  }
  else if (type == "cylinder" ) {
    Vector<RealType> centerVector, scalingVector;
    parser.getRealVec ( "center", centerVector );
    parser.getRealVec ( "scaling", scalingVector );
    Vec3<RealType> center ( centerVector.getData() ), scaling ( scalingVector.getData() );
    RealType radius = parser.getDouble ( "radius" );
    clog << "generating cylinder with center = " << centerVector << "," << endl
         << ", scaling = " << scalingVector << " and radius " << radius << endl;
    generator.generateCylindricalLevelset ( levelValues, center, scaling, radius );
  }
  else {
    cerr << aol::color::red << "level set type \"" << type << "\" not supported." << aol::color::reset << endl;
    return;
  }

  clog << "saving as " << output << " ... ";
  levelValues.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

void move ( const ArrayType & input, ArrayType & output, int compToMove, RealType moveAmount ) {

  int d = static_cast<int> ( input.getSize()[compToMove] * moveAmount );

  short n_toMove = input.getSize()[compToMove],
        n0 = input.getSize()[0],
        n1 = input.getSize()[1],
        n2 = input.getSize()[2];

  CoordType index;
  for ( index[0] = 0; index[0] < n0; ++index[0] )
    for ( index[1] = 0; index[1] < n1; ++index[1] )
      for ( index[2] = 0; index[2] < n2; ++index[2] ) {
        CoordType movedIndex = index;
        short notMoved = index[compToMove];
        movedIndex[compToMove] = notMoved + d >= 0  // is moved index still inside?
                                 ? notMoved + d < n_toMove
                                   ? notMoved + d   // take index i + d if inside
                                   : n_toMove - 1   // else highest
                                 : 0;         // or lowest index
        output.set ( index, input.get ( movedIndex ) );
      }
}

void move ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  string input  = parser.getString ( "input" );

  int compToMove = parser.getInt ( "compToMove" );
  RealType moveAmount = parser.getDouble ( "moveAmount" );

  clog << "moving " << input << " about " << moveAmount << " in coordinate direction " << compToMove << endl;

  ArrayType in  ( input.c_str() ),
            out ( in, STRUCT_COPY );

  move ( in, out, compToMove, moveAmount );

  clog << "saving as " << output << " ... ";
  out.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

void reflect ( ArrayType & output, const ArrayType & input, int compToReflect ) {

  int r = compToReflect;
  GridSize3d size ( input );

  CoordType index;
  for ( index[0] = 0; index[0] < size[0]; ++index[0] )
    for ( index[1] = 0; index[1] < size[1]; ++index[1] )
      for ( index[2] = 0; index[2] < size[2]; ++index[2] ) {
        CoordType rIndex = index;
        rIndex[r] = size[r] - 1 - index[r];
        output.set ( rIndex, input.get ( index ) );
      }
}

void reflect ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  string input  = parser.getString ( "input" );

  int compToReflect = parser.getInt ( "compToReflect" );

  if (compToReflect < 0 || compToReflect > 2) {
    cerr << aol::color::red << "parameter \"compToReflect\" has to be between 0 and 2, "
            "you passed " << compToReflect << aol::color::reset << endl;
    return;
  }

  string numeral = (   compToReflect == 0 ? "0th"
                   : ( compToReflect == 1 ? "1st" : "2nd" ) );

  clog << "reflecting " << numeral << " component of " << input << endl;

  ArrayType in ( input.c_str() ),
            out ( in, STRUCT_COPY );

  reflect ( out, in, compToReflect );
  clog << "saving as " << output << " ... ";
  out.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

void reorderComps ( ArrayType & output, const ArrayType & input, Vec3<int> newOrder ) {

  GridSize3d size ( input );

  CoordType index;
  for ( index[0] = 0; index[0] < size[0]; ++index[0] )
    for ( index[1] = 0; index[1] < size[1]; ++index[1] )
      for ( index[2] = 0; index[2] < size[2]; ++index[2] ) {
        CoordType newIndex;
        for (int i = 0; i < 3; ++i)
          newIndex[i] = index[newOrder[i]];
        output.set ( newIndex, input.get ( index ) );
      }
}

void reorderComps ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  string input  = parser.getString ( "input" );

  Vector<int> newOrderVector;
  parser.getIntVec ( "newOrder", newOrderVector );
  Vec3<int> newOrder ( newOrderVector.getData() );

  clog << "reordering components of " << input << ":" << endl
       << "0 -> " << newOrder[0] << endl
       << "1 -> " << newOrder[1] << endl
       << "2 -> " << newOrder[2] << endl;

  ArrayType in ( input.c_str() ),
            out ( in, STRUCT_COPY );

  reorderComps ( out, in, newOrder );
  clog << "saving as " << output << " ... ";
  out.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

void cleanLevelFct ( ArrayType & inOut, RealType newValue ) {
  tpcfe::SmallComponentsPoempelRemover<RealType> remover ( newValue );
  remover.applySingle ( inOut );
}

void cleanLevelFct ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  string input  = parser.getString ( "input" );
  RealType newValue = parser.getDouble ( "newValue" );
  bool reflectInsideOut = (parser.getString ( "reflectInsideOut" ) == "yes");

  clog << "loading input array " << input << " ... ";
  ArrayType inOut ( input.c_str() );
  clog << "done." << endl;

  clog << "cleaning level function (newValue = " << newValue << ") ... " << endl;
  if ( reflectInsideOut )
    inOut *= -1.;
  cleanLevelFct ( inOut, newValue );
  if ( reflectInsideOut )
    inOut *= -1.;
  clog << "done." << endl;

  clog << "saving as " << output << " ... ";
  inOut.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

int main( int argc, char ** argv )
{
  if (argc < 2) {
    cerr << aol::color::error << "Usage: " << argv[0] << " <parameter files>" << endl
         << aol::color::reset;
    return -1;
  }

  try {
    for ( int i = 1; i < argc; ++i ) {
      aol::ParameterParser parser ( argv[i] );
      string action = parser.getString ( "action" );

      clog << aol::color::blue;

      if ( action == "generateLevelFunction3D" )      generateLevelFunction3D ( parser );
      else if ( action == "setConst" )                setConst ( parser );
      else if ( action == "addConst" )                addConst ( parser );
      else if ( action == "linComb" )                 linComb ( parser );
      else if ( action == "maxOfTwo" )                maxOfTwo ( parser );
      else if ( action == "minOfTwo" )                minOfTwo ( parser );
      else if ( action == "product" )                 product ( parser );
      else if ( action == "smooth" )                  smooth ( parser );
      else if ( action == "move" )                    move ( parser );
      else if ( action == "reflect" )                 reflect ( parser );
      else if ( action == "reorderComps" )            reorderComps ( parser );
      else if ( action == "redistance" )              redistance<ConfType> ( parser );
      else if ( action == "cleanLevelFct" )           cleanLevelFct ( parser );
      else
        clog << aol::color::red << "Unknown action \"" << action << "\"." << endl;

    }
  }
  catch (Exception & exc) {
    exc.dump();
  }
  catch (...) {
  }

  clog << aol::color::reset;
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
//------------------------------------------------------------------------------------------------
