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
 *  \brief Small utility classes in the aol that are not related to numerics
 *
 *  This is a program that shows how to use "small utility classes" in
 *  the aol. These are tools that are not for image processing or
 *  numerics and should be in alphabetical order.
 *
 *  \author Schwen
 */

#include<epswriter.h>
#include<rgbColorMap.h>
#include<parameterParser.h>
#include<progressBar.h>
#include<randomGenerator.h>


using namespace std;

int main ( int, char** ) {
  try {

    {
      cerr << "Example for colored shell output (defined in platformDependent.h)" << endl;
      cerr << "**********************************************************" << endl;

      cerr << "Printing text in " << aol::color::red << " color " << aol::color::reset << " but not the rest of the line" << endl;

      cerr << "More colors: "
      << aol::color::black << " black "
      << aol::color::red << " red "
      << aol::color::green << " green "
      << aol::color::brown << " brown "
      << aol::color::blue << " blue "
      << aol::color::purple << " purple "
      << aol::color::cyan << " cyan "
      << aol::color::light_grey << " light_grey "
      << aol::color::dark_grey << " dark_grey "
      << aol::color::light_red << " light_red "
      << aol::color::light_green << " light_green "
      << aol::color::yellow << " yellow "
      << aol::color::light_blue << " light_blue "
      << aol::color::pink << " pink "
      << aol::color::light_cyan << " light_cyan "
      << aol::color::white << " white "  << aol::color::reset << endl
      << "And some special orders: "
      << aol::color::beep << " beep "
      << aol::color::error << " error "
      << aol::color::ok << " ok "
      << aol::color::residuum << " residuum "
      << aol::color::invert << " invert "
      << aol::color::reset << " reset gives black again! " ;
    }

    {
      cerr << endl << endl;
      cerr << "Example for EpsWriter, a class to write eps files in color" << endl;
      cerr << "**********************************************************" << endl;

      aol::EpsWriter epsWriter ( "fourty-two.eps", aol::HSV_BLUE_TO_RED );

      /*                    x0      y0      x1      y1      thickness color */
      epsWriter.writeLine ( 100000, 300000, 100000, 450000, 50000,  0 );
      epsWriter.writeLine ( 100000, 450000, 300000, 450000, 50000, 30 );
      epsWriter.writeLine ( 300000, 300000, 300000, 600000, 50000, 60 );

      epsWriter.writeLine ( 400000, 600000, 600000, 600000, 50000,  90 );
      epsWriter.writeLine ( 400000, 600000, 400000, 450000, 50000, 120 );
      epsWriter.writeLine ( 400000, 450000, 600000, 450000, 50000, 150 );
      epsWriter.writeLine ( 600000, 450000, 600000, 300000, 50000, 180 );
      epsWriter.writeLine ( 400000, 300000, 600000, 300000, 50000, 210 );
      // output file closed by destructor of EpsWriter
    }

    {
      cerr << endl << endl;
      cerr << "Example for nicely formatted output (defined in aol.cpp)" << endl;
      cerr << "**********************************************************" << endl;

      const double exactNumber = 1.234567891011121314151617181920;

      cout  << "scientificFormat    : " << aol::scientificFormat ( exactNumber ) << endl;
      cout  << "longScientificFormat: " << aol::longScientificFormat ( exactNumber ) << endl;
      cout  << "shortFormat         : " << aol::shortFormat ( exactNumber ) << endl;
      cout  << "intFormat           : " << aol::intFormat ( exactNumber ) << endl;
      cout  << "fillFormat          : " << aol::fillFormat ( exactNumber ) << endl;
      cout  << "mixedFormat         : " << aol::mixedFormat ( exactNumber ) << endl;
      cout  << "detailedFormat      : " << aol::detailedFormat ( exactNumber ) << endl;

    }

    {
      cerr << "Example for ParameterParser, a parameter parser" << endl;
      cerr << "***********************************************" << endl;

      aol::ParameterParser parser ( "parameters.par" );                     // read parameters from given filename

      parser.dump();                                                        // print all variable names and values read from parameter file

      int num_steps = parser.getInt ( "num_steps" );                        // get integer
      typedef float RealType;
      float alpha = static_cast<RealType> ( parser.getDouble ( "alpha" ) ); // sometimes, casts can be necessary
      char filename[1024];
      parser.getString ( "filename", filename );                            // get a string
      bool has_beta = parser.hasVariable ( "beta" );                        // check whether some variable is in the parameter file
      cerr << num_steps << " " << alpha << " " << filename << " " << has_beta << endl;
    }

    {
      cerr << endl << endl;
      cerr << "Example for ProgressBar, a text console progress bar" << endl;
      cerr << "****************************************************" << endl;

      const int i = 400, j = 400;
      double values[i][j];

      aol::ProgressBar<> pb ( "Doing something" );   // Text that is shown
      pb.start ( i * j );                            // exact total number of steps (for being able to compute progess in percent)
      for ( int m = 0; m < i; ++m ) {
        for ( int n = 0; n < j; ++n, pb++ ) {        // pb only has postfix increment, for fine steps, it should be incremented in an inner loop
          values[m][n] = 2.818 * m + 1.618 * n;      // there should be no console output while the progress bar is printing
        }
      }
      pb.finish();                                   // essentially to continue console output in next line
      cerr << "Some result: " << values[42][23] << endl;
    }

    {
      cerr << endl << endl;
      cerr << "Example for RandomGenerator, a class that produces random data of various types " << endl;
      cerr << "******************************************************************************* " << endl;

      aol::RandomGenerator rg ( 31415 ); // specify a random seed: this way, different but reproducible sequences of random data can be generated

      cerr << "Random Integer in the range [23, 42) " << rg.rInt ( 23, 42 ) << endl;
      cerr << "Random double in the range [1.414, 3.141) " << rg.rReal ( 1.414, 3.141 ) << endl;
      cerr << "Random bool " << rg.rBool() << endl;
    }

    {
      cerr << endl << endl;
      cerr << "Example for StopWatch, a class to measure computation time" << endl;
      cerr << "**********************************************************" << endl;

      const int i = 800, j = 800;
      aol::Vector<double> vec ( i * j );    // this is just some computations that takes some time and is not optimized away

      aol::StopWatch timer;
      timer.start();
      for ( int m = 0; m < i; ++m ) {
        for ( int n = 0; n < j; ++n ) {
          vec[i*m+n] = 2.818 * m + 1.618 * n;
        }
      }
      timer.stop();
      cerr << "Doing someting else took " << timer.elapsedCpuTime() << " seconds, that is " << timer.elapsedCpuTimeString() << endl;
      cerr << "It started at " << timer.startedAt() << " and finished at " << timer.stoppedAt() << ", thus took " << timer.elapsedCpuTimeString() << "." << endl;
    }

    return ( EXIT_SUCCESS );

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
  return ( EXIT_FAILURE );
}
