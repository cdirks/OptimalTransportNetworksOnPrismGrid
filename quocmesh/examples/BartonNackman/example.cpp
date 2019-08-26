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

/**
 * \file
 * \brief example for the so-called Barton--Nackman-Trick.
 *
 * The so-called Barton--Nackman-Trick is a way to implement polymerphism in a static way, i.e. write abstract basis classes
 * without virtual methods, where the decision which method is called is made at compile time.
 *
 * The example contains an abstract vector class SomeVector that just knows that entries can be read and written. One can write algorithms
 * ( in this example: SomeVector::print() and operator==() ) on this class that do not know about the actual implementation, the
 * implementation-dependent access functions ( here: SomeVector::get() and SomeVector::set() ) just delegate work down to the actual implementation
 * (of which two versions, a fully populated and a sparse vector, are supplied in this example). Thus, the basis class has to know
 * the type of the implementation via a template parameter (here: \code Imp \endcode) and must be able to cast itself into an instance of
 * the implementation type ( here: via SomeVector::asImp() ).
 *
 * \author Lenz
 */


#include <iostream>
#include <map>

#include "aol.h"

// Abstract vector, does not contain implementation
template <class Imp> class SomeVector {
public:
  // Constructor
  SomeVector (int size) : _size (size) {}
  int size () const { return _size; }
  // Convert to implementation type
  Imp& asImp () { return * static_cast<Imp*> (this); }
  const Imp& asImp () const { return * static_cast<const Imp*> (this); }
  // Access functions, just delegates access to implemtation
  double get (int i) const { return asImp ().get (i); }
  double set (int i, double v) { return asImp ().set (i, v); }
  // Algorithm written for all vectors
  void print () const { for (int i = 0; i < _size; ++i) std::cout << aol::mixedFormat ( this->get (i) ) << " "; std::cout << std::endl; }
protected:
  int _size;
};

class FullVector : public SomeVector<FullVector> {
public:
  // Constructor
  FullVector (int size) : SomeVector<FullVector> (size), _data (new double [size]) {}
  // Destructor
  ~FullVector () { if (_data) delete [] _data; };
  // Implementation of access functions
  double get (int i) const { return _data [i]; }
  void set (int i, double v) { _data [i] = v; }
private:
  double* _data;
};

class SparseVector : public SomeVector<SparseVector> {
public:
  // Constructor
  SparseVector (int size) : SomeVector<SparseVector> (size) {}
  // Implementation of access function
  double get (int i) const { std::map<int,double>::const_iterator it = _data.find (i); if (it != _data.end ()) return it->second; else return 0; }
  void set (int i, double v) { if (v) _data [i] = v; else _data.erase (i); }
private:
  std::map<int,double> _data;
};

template <class lImp, class rImp>
bool operator == (const SomeVector<lImp>& lv, const SomeVector<rImp>& rv)
{
  if (lv.size () != rv.size ()) return false;
  for (int i = 0; i < lv.size (); ++i) if (lv.get(i) != rv.get(i)) return false;
  return true;
}

int main ()
{
  FullVector fv (10);
  SparseVector sv (10);

  fv.set (4, 4.14); // Uses access functions in implementation
  sv.set (4, 4.14); // Uses access functions in implementation

  fv.print (); // Uses algorithm from abstract base class, which then uses the access functions in the implementation
  sv.print (); // Uses algorithm from abstract base class, which then uses the access functions in the implementation

  if (fv == sv) std::cout << "Vectors are equal." << endl;

  return 0;
}
