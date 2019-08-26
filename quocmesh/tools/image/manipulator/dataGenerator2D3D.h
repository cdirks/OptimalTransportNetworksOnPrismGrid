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

//------------------------------------------------------------------------------------------------

void setConst ( ArrayType & output, RealType value ) {
  output.setAll( value );
}

void setConst ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  RealType value = parser.getDouble ( "value" );
  int level = parser.getInt ( "gridLevel" );

  clog << "setting " << output << " to value " << value << endl;

  ConfType::InitType grid ( level, ConfType::Dim );
  ArrayType z ( grid );

  setConst ( z, value );
  clog << "saving as " << output << " ... ";
  z.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

void addConst ( ArrayType & output, const ArrayType & input, RealType value ) {
  output = input;
  output.addToAll ( value );
}

void addConst ( ParameterParser & parser ) {

  string output  = parser.getString ( "output" );
  string input   = parser.getString ( "input" );
  RealType value = parser.getDouble ( "value" );

  clog << "adding " << value << " to " << input << endl;

  ArrayType x ( input.c_str() ),
            z ( x, STRUCT_COPY );

  addConst ( z, x, value );
  clog << "saving as " << output << " ... ";
  z.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

//! writes coeff1 input1 + coeff2 input2 into output
void linComb ( ArrayType & output, const ArrayType & input1, const ArrayType & input2,
               RealType coeff1, RealType coeff2 ) {
  output.setZero();
  output.addMultiple ( input1, coeff1 );
  output.addMultiple ( input2, coeff2 );
}

void linComb ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  string input1 = parser.getString ( "input1" );
  string input2 = parser.getString ( "input2" );

  RealType coeff1 = parser.getDouble ( "coeff1" );
  RealType coeff2 = parser.getDouble ( "coeff2" );

  clog << "computing " << coeff1 << " * " << input1 << " + " << coeff2 << " * " << input2 << endl;

  ArrayType x ( input1.c_str() ),
            y ( input2.c_str() ),
            z ( x, STRUCT_COPY );

  linComb ( z, x, y, coeff1, coeff2 );
  clog << "saving as " << output << " ... ";
  z.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

void product ( ArrayType & output, const ArrayType & input1, const ArrayType & input2 ) {

  for ( int i = 0; i < output.size(); ++i )
    output[i] = input1[i] * input2[i];
}

void product ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  string input1 = parser.getString ( "input1" );
  string input2 = parser.getString ( "input2" );

  clog << "computing product of " << input1 << " and " << input2 << endl;

  ArrayType x ( input1.c_str() ),
            y ( input2.c_str() ),
            z ( x, STRUCT_COPY );

  product ( z, x, y );
  clog << "saving as " << output << " ... ";
  z.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

void maxOfTwo ( ArrayType & output, const ArrayType & input1, const ArrayType & input2 ) {

  for ( int i = 0; i < output.size(); ++i )
    output[i] = Max ( input1[i], input2[i] );
}

void maxOfTwo ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  string input1 = parser.getString ( "input1" );
  string input2 = parser.getString ( "input2" );

  clog << "computing max of " << input1 << " and " << input2 << endl;

  ArrayType x ( input1.c_str() ),
            y ( input2.c_str() ),
            z ( x, STRUCT_COPY );

  maxOfTwo ( z, x, y );
  clog << "saving as " << output << " ... ";
  z.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

void minOfTwo ( ArrayType & output, const ArrayType & input1, const ArrayType & input2 ) {

  for ( int i = 0; i < output.size(); ++i )
    output[i] = Min ( input1[i], input2[i] );
}

void minOfTwo ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  string input1 = parser.getString ( "input1" );
  string input2 = parser.getString ( "input2" );

  clog << "computing min of " << input1 << " and " << input2 << endl;

  ArrayType x ( input1.c_str() ),
            y ( input2.c_str() ),
            z ( x, STRUCT_COPY );

  minOfTwo ( z, x, y );
  clog << "saving as " << output << " ... ";
  z.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

void smooth ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  string input  = parser.getString ( "input" );
  RealType sigma = parser.getDouble ( "sigma" );

  clog << "smoothing " << input << " with sigma = " << sigma << endl;

  ArrayType in ( input.c_str() );
  ConfType::InitType grid ( in.getSize() );
  ArrayType out ( in, STRUCT_COPY );

  LinearSmoothOp<RealType> smoothOp ( grid );
  smoothOp.setSigma ( sigma );
  smoothOp.apply ( in, out );

  clog << "saving as " << output << " ... ";
  out.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
//------------------------------------------------------------------------------------------------

template <typename ConfiguratorType>
void redistance ( typename ConfiguratorType::ArrayType & output, const typename ConfiguratorType::ArrayType & input ) {
  typename ConfiguratorType::InitType grid ( qc::GridSize<ConfiguratorType::Dim>::createFrom ( input ) );
  typename qc::SignedDistanceOpTrait<ConfiguratorType, ConfiguratorType::Dim>::OpType signedDistanceOp ( grid );
  signedDistanceOp.apply ( input, output );
}

template <typename ConfiguratorType>
void redistance ( ParameterParser & parser ) {

  string output = parser.getString ( "output" );
  string input  = parser.getString ( "input" );

  clog << "loading input array " << input << " ... ";
  typename ConfiguratorType::ArrayType in ( input.c_str() ),
            out ( in, STRUCT_COPY );
  clog << "done." << endl;

  clog << "performing redistancing ... " << endl;
  redistance<ConfiguratorType> ( out, in );
  clog << "done." << endl;

  clog << "saving as " << output << " ... ";
  out.save ( output.c_str(), saveType );
  clog << "done." << endl;
}
