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

#include <elasticTensor.h>
#include <typedParser.h>
#include <qmException.h>

namespace aol {

// Class TypedParameterParser

TypedParameterParser::TypedParameterParser ()
    : _values () {}

TypedParameterParser::TypedParameterParser ( const TypedParameterParser& parser )
    : _values ( parser._values ) {}

TypedParameterParser& TypedParameterParser::operator = ( const TypedParameterParser& parser ) {
  // Beware of self-assignment
  if ( this == &parser ) return *this;

  _values = parser._values;

  return *this;
}

TypedParameterParser::~TypedParameterParser () {}

TypedParameterParser::TypedParameterParser ( int argc, char* argv[] )
    : _values () {
  string parameterFileName;

  switch ( argc ) {
  case 1:
    parameterFileName = "parameter";
    break;
  case 2:
    parameterFileName = argv [1];
    break;
  default:
    cerr << "Usage: " << argv [0] << " [parameterfilename]" << endl;
    throw ParameterException ( "TypedParameterParser::TypedParameterParser, too many parameters", __FILE__, __LINE__ );
  }

  ifstream parameterFile ( parameterFileName.c_str () );
  read ( parameterFile );
}

TypedParameterParser::TypedParameterParser ( istream& stream )
    : _values () {
  read ( stream );
}

TypedParameterParser::TypedParameterParser ( const string filename, const string defaultParameterFilename )
    : _values () {
  ifstream parameterFile ( filename.c_str () );

  // First read default parameters...
  if ( defaultParameterFilename.length () > 0 ) {
    ifstream defaultParameterFile ( defaultParameterFilename.c_str () );
    read ( defaultParameterFile );
  }

  // ... s.t. they can be overwritten by parameters in the parameter file.
  read ( parameterFile, true );
}

void TypedParameterParser::read ( istream& stream, const bool allowDoubleParameters ) {
  if ( !stream ) throw FileDoesNotExistException ( "TypedParameterParser::read", __FILE__, __LINE__ );

  // Extend this for vadditional types
  map <string, TypePtr> types;
  types["int"]           = &typeid ( int );
  types["uint"]          = &typeid ( unsigned int );
  types["double"]        = &typeid ( double );
  types["string"]        = &typeid ( string );
  types["bool"]          = &typeid ( bool );
  types["Elast"]         = &typeid ( aol::ElasticTensor<double> );
  types["Vec2"]          = &typeid ( aol::Vec2<double> );
  types["Vec3"]          = &typeid ( aol::Vec3<double> );
  types["Mat22"]         = &typeid ( aol::Matrix22<double> );
  types["Mat33"]         = &typeid ( aol::Matrix33<double> );
  types["VecInt"]        = &typeid ( aol::Vector<int> );
  types["VecDouble"]     = &typeid ( aol::Vector<double> );

  string line;

  while ( getline ( stream, line ) ) {

    std::string::size_type hash = line.find ( "#" ); // Cut comments
    if ( hash != string::npos ) line = line.substr ( 0, hash );

    stringstream temp ( line );
    string type;
    string identifier;
    string value;

    temp >> type;
    if ( type == "" ) continue;

    // Check if type is a known type
    if ( types.find ( type ) == types.end () )
      throw FileFormatException ( "TypedParameterParser::read, type " + type + " not recognized.", __FILE__, __LINE__ );

    temp >> identifier;
    getline ( temp, value );

    if ( !allowDoubleParameters ) {
      if ( _values.find ( identifier ) != _values.end () )
        throw FileFormatException ( "TypedParameterParser::read, variable " + identifier + " was found more than once in parameter file", __FILE__, __LINE__ );
    }

    _values [identifier] = Variable ( types [type], identifier, string ( value ), type );
  }
}

TypedParameterParser::Variable TypedParameterParser::getVariable ( string identifier ) const {
  map<string, Variable>::const_iterator found = _values.find ( identifier );
  if ( found == _values.end () ) throw FileFormatException ( "TypedParameterParser::getVariable, variable " + identifier + " not found in parameter file", __FILE__, __LINE__ );
  return ( * ( found ) ).second;
}

bool TypedParameterParser::hasVariable ( string identifier ) const {
  map<string, Variable>::const_iterator found = _values.find ( identifier );
  return ( found != _values.end () );
}

void TypedParameterParser::dump ( string filename ) const {
  ofstream outputFile ( filename.c_str () );
  if ( !outputFile.is_open () ) {
    throw aol::Exception ( "Failed to open output file", __FILE__, __LINE__, __FUNCTION__ );
  }
  for ( map<string, Variable>::const_iterator iter = _values.begin(); iter != _values.end(); ++iter )  {
    string tmp;
    tmp = iter->second.getValueString();
    // Trim leading whitespaces
    tmp.erase ( tmp.begin (), std::find_if ( tmp.begin (), tmp.end (), std::not1 ( std::ptr_fun<int, int> ( std::isspace ) ) ) );
    outputFile << setw(15) << left << iter->second.getTypeString() << setw(30) << left << iter->first << tmp << endl;
  }
  outputFile.close();
}

// Class TypedParameterParser::Variable

TypedParameterParser::Variable::Variable ()
  : _type ( 0 ), _identifier (), _value (), _typeName () {}

TypedParameterParser::Variable::Variable ( const TypedParameterParser::Variable& variable )
  : _type ( variable._type ), _identifier ( variable._identifier ), _value ( variable._value ), _typeName ( variable._typeName ) {}

TypedParameterParser::Variable& TypedParameterParser::Variable::operator = ( const TypedParameterParser::Variable& variable ) {
  // Beware of self-assignment
  if ( this == &variable ) return *this;

  _type = variable._type;
  _identifier = variable._identifier;
  _value = variable._value;
  _typeName = variable._typeName;

  return *this;
}

TypedParameterParser::Variable::~Variable () {}

TypedParameterParser::Variable::Variable ( TypePtr type, string identifier, string value, string typeName )
    : _type ( type ), _identifier ( identifier ), _value ( value ), _typeName ( typeName ) {}

const string & TypedParameterParser::Variable::getValueString () const {
  return _value;
}

const string & TypedParameterParser::Variable::getTypeString () const {
  return _typeName;
}

}
