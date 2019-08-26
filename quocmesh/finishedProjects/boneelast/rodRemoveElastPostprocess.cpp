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

#include <aol.h>


int main ( int argc, char** argv ) {
  if ( argc != 2 ) {
    cerr << "usage: " << argv[0] << " infile " << endl;
    return ( EXIT_FAILURE );
  }
  try {

    std::vector< int > percentages;
    std::vector< double > Emoduli, stddevis, linErrors;

    double Emax = 0.0, Emin = 0.0, avgLinError = 0.0;
    const int min_p = 0, max_p = 30;

    ifstream datafile ( argv[1] );
    if ( !datafile ) {
      cerr << "could not open file " << endl;
      return ( EXIT_FAILURE );
    }

    while (!datafile.eof() ) {

      int percentage = 0;
      double Emodulus = 0.0, stddevi = 0.0;
      datafile >> percentage;
      datafile >> Emodulus;
      datafile >> stddevi;

      if ( Emodulus != 0.0 ) {
        percentages.push_back( percentage );
        Emoduli.push_back( Emodulus );
        stddevis.push_back( stddevi );
      }

    }

    { // compute average minimal and maximal E moduli
      double dwarf_min = 0.0, dwarf_max = 0.0; int n_avg_min = 0, n_avg_max = 0;
      for ( unsigned int i = 0; i < percentages.size(); ++i ) {
        if ( percentages[i] == max_p ) {
          dwarf_min += Emoduli[i];
          ++n_avg_min;
        }
        if ( percentages[i] == min_p ) {
          dwarf_max += Emoduli[i];
          ++n_avg_max;
        }
      }
      Emin = dwarf_min / n_avg_min;
      Emax = dwarf_max / n_avg_max;
    }

    linErrors.resize( percentages.size() );
    double lin_error_sum = 0.0;
    for ( unsigned int i = 0; i < percentages.size(); ++i ) {
      linErrors[i] = aol::Abs ( Emoduli[i] - ( Emax + ( ( 1.0 / ( max_p - min_p ) ) * ( percentages[i] - min_p ) ) * ( Emin - Emax ) ) );
      lin_error_sum += linErrors[i];
    }

    avgLinError = lin_error_sum / percentages.size();
    avgLinError /= ( ( Emin + Emax ) / 2 );

    for ( unsigned int i = 0; i < percentages.size(); ++i ) {
      //      cerr << percentages[i] << "  " << Emoduli[i] << "  " << stddevis[i] << "  " << linErrors[i] << endl;
    }

    cerr << "Emax for full structure: " << Emax << ",   Emin for reduced structure: " << Emin << ", Emin relative to Emax: " << Emin / Emax << ",   average absolute difference to secant: " << avgLinError << endl;
    cout << ", " << Emax << " + " << Emax - Emin << " / 30 *  x w l " << endl;
    cout << Emax << " " << Emin <<  " " << Emin / Emax << " " << avgLinError << endl;

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}


#if 0
// multiply AN to FN standard deviations by sqrt (129)
int main ( int argc, char** argv ) {
  if ( argc != 3 ) {
    cerr << "usage: " << argv[0] << " infile outfile" << endl;
    return ( EXIT_FAILURE );
  }
  try {

    ifstream ins  ( argv[1] );
    ofstream outs ( argv[2] );
    if ( !ins || !outs) {
      cerr << "could not open file " << endl;
      return ( EXIT_FAILURE );
    }

    while (!ins.eof() ) {

      double percentage, value, error;
      ins >> percentage;
      if ( ins.eof() )
        break;
      ins >> value;
      ins >> error;

      error *= sqrt ( 129. );
      aol::MixedFormat pcF ( 2, 3 ), vlF ( 5, 9 );

      outs << pcF( percentage ) << "  " << vlF ( value ) << "  " << vlF ( error ) << endl;
    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
#endif

#if 0
// multiply OSMR01 to 07 (except 04a,b) by sqrt(129)
int main ( int argc, char** argv ) {
  if ( argc != 3 ) {
    cerr << "usage: " << argv[0] << " infile outfile" << endl;
    return ( EXIT_FAILURE );
  }
  try {

    ifstream ins  ( argv[1] );
    ofstream outs ( argv[2] );
    if ( !ins || !outs) {
      cerr << "could not open file " << endl;
      return ( EXIT_FAILURE );
    }

    while (!ins.eof() ) {

      double percentage1, percentage2, value, error;
      ins >> percentage1;
      if ( ins.eof() )
        break;
      ins >> percentage2;
      ins >> value;
      ins >> error;

      error *= sqrt ( 129. );
      aol::MixedFormat pcF ( 2, 3 ), vlF ( 5, 9 );

      outs << pcF( percentage1 ) << "  " << pcF( percentage2 ) << "  " << vlF ( value ) << "  " << vlF ( error ) << endl;
    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
#endif
