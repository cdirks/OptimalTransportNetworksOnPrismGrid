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

#ifndef __TYPEDPARSER_H
#define __TYPEDPARSER_H

#include <aol.h>
#include <qmException.h>

namespace aol {

/******************************************************************************************************/
//! A Simple Parameter Parser that can read parameter files of the format "type identifier value".
//! It can easily be extended for new types by adding them to the table in TypedParameterParser::read
//! Each type must understand the operator >> that is used to interpret the value
//! Run-Time Type safety is guaranteed via typeid ()
class TypedParameterParser {

private:

  //! Type information that can be compared
  typedef const type_info* TypePtr;

  /**********************************************************************************/
  //! Storage for a variable type information, identifier, and (still-to-decode) value
  class Variable {

  private:

    //! Type information.
    //! Will not be deleted, since memory is owned by run time system
    TypePtr _type;
    //! Name identifier
    string _identifier;
    //! Value as a string, will be decoded to appropriate type by operator >>
    string _value;
    //! Internal type name
    string _typeName;

  public:

    //! Default constructor for STL compatibility
    Variable ();

    //! Copy Constructor
    Variable ( const Variable& variable );

    //! Assignment Operator
    Variable& operator = ( const Variable& variable );

    //! Destructor
    virtual ~Variable ();

    //! Construct from its members
    Variable ( TypePtr type, string identifier, string value, string typeName );

    //! Decode value to appropriate type.
    //! Type safety guaranteed at run time by typeid comparison
    template <class DataType> void getValue ( DataType& value ) const;

    //! Get value as a string
    const string & getValueString () const;

    //! Get type name as a string
    const string & getTypeString () const;
  };

  //! Get Variable for specified identifier
  Variable getVariable ( string identifier ) const;

public:

  //! Default Constructor
  TypedParameterParser ();

  //! Copy Constructor
  TypedParameterParser ( const TypedParameterParser& parser );

  //! Assignment Operator
  TypedParameterParser& operator = ( const TypedParameterParser& parser );

  //! Destructor
  virtual ~TypedParameterParser () ;

  //! Read Parameter Values from File
  explicit TypedParameterParser ( istream& stream );

  //! Read values from specified file or default
  TypedParameterParser ( int argc, char* argv[] );

  //! Read values from a file given by filename, possibly give another parameter file with
  //! default parameters, whichs values will be overwritten by the values in filename,
  //! if they are specified in both files.
  TypedParameterParser ( string filename, string defaultParameterFilename = "" );

  //! Decode value of given identifier to appropriate type.
  //! Type safety is guaranteed by Variable::get<DataType>
  template <class DataType> void get ( string identifier, DataType& value ) const;

  //! Read data from file, do not delete old entries
  //! \param[in] allowDoubleParameters Allows parameters to be defined more than once (needed for default parameters).
  void read ( istream& stream, const bool allowDoubleParameters = false );

  //! Check if Variable for specified identifier exists
  bool hasVariable ( string identifier ) const;

  //! Save current parameters
  void dump ( string filename ) const;

private:

  //! Stores Identifier--Variable Pairs
  map <string, Variable> _values;
};

// Template functions (remainder in parser.cpp)

template <class DataType> void TypedParameterParser::get ( string identifier, DataType& value ) const {
    Variable variable = getVariable ( identifier );
    variable.getValue ( value );
  }

template <class DataType> void TypedParameterParser::Variable::getValue ( DataType& value ) const {
  if ( *_type != typeid ( DataType ) )
    throw aol::TypeException ( "TypedParameterParser::Variable::getValue, for Variable " + _identifier, __FILE__, __LINE__ );
  stringstream temp ( _value );
  temp >> value;
}

}

#endif
