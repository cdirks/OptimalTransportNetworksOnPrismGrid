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

#include <configurators.h>
#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define RealType double

// HACK: For printing
#define X 32
#define Y 32
#define Z 32


void bubble_sort (RealType *list, int n)
{
  RealType tmp;
  int i, j;

  for (i=0; i<n-1; i++)
    {
    for (j=0; j<n-1-i; j++)
      if (list[j+1] < list[j])
        {
        tmp = list[j];
        list[j] = list[j+1];
        list[j+1] = tmp;
        }
    }
}


RealType getVariance( RealType *list, RealType E, int n )
{
  RealType v = 0.;

  for (int i=0; i<n; i++)
    v += aol::Sqr( E-list[i] );

  return v;
}


int MaxFn (int a, int b)
{
  if (a>b) return a; else return b;
}

int MinFn (int a, int b)
{
  if (a>=b) return b; else return a;
}

int main( int argc, char **argv)
{
  try
    {
      if ( argc < 5 )
  {
   cerr <<aol::color::red<< "usage: "<<argv[0]<<" input output max_radius step_radius (variance_output) \n";
          cerr<<aol::color::reset;
    return EXIT_FAILURE;
  }
      qc::ScalarArray<RealType, qc::QC_3D> inputimg( argv[1]);
      int N = inputimg.getNumX ();
      int d = qc::logBaseTwo (N);
      qc::GridDefinition grid( d, qc::QC_3D );
      qc::ScalarArray<RealType, qc::QC_3D> outputimg(grid);
      qc::ScalarArray<RealType, qc::QC_3D> varianceImg(grid);      // to evaluate the radii estimation

      int cx, cy, cz, x, y, z, i, M0, ri;
      RealType maxr, stepr, r, r2;
      maxr = atof(argv[3]);
      stepr = atof(argv[4]);
      int nr = 1 + static_cast<int>(floor((maxr-1.0)*2.0));
      int nrmed = static_cast<int>(rint(0.5*static_cast<RealType>(nr)));
      RealType *vrl = static_cast<RealType*>(calloc(nr,sizeof(RealType)));


      // just for displaying the progress
      RealType percent = 0;

      for (cx = 0; cx<N; cx++)
      {
        // display progress
        percent = static_cast<RealType> ( cx * 100. / N );
        cerr.precision(4);
        cerr<<"\rCalculating..."<<percent<<" %            ";

        for (cy = 0; cy<N; cy++)
          for (cz = 0; cz<N; cz++)
            {

//               if (inputimg.get(cx,cy,cz)>0.)
//                     cerr<<"("<<cx<<","<<cy<<","<<cz<<"); ";

              // HACK: For printing the values
//               if (cx == X && cy == Y && cz == Z) cerr<<"\n\n";


              r = 1.;
              i = 0;
              while (r<=maxr)
                {
                  M0 = 0;
                  ri = static_cast<int>(ceil(r));
                  r2 = r*r;
                  for (x = MaxFn(cx-ri,0); x<= MinFn(cx+ri,N-1); x++)
                    for (y = MaxFn(cy-ri,0); y<= MinFn(cy+ri,N-1); y++)
                      for (z = MaxFn(cz-ri,0); z<= MinFn(cz+ri,N-1); z++)
                        {
                          if ((static_cast<RealType>((x-cx)*(x-cx)+(y-cy)*(y-cy)+(z-cz)*(z-cz)))<=r2)
                            {
                              if (inputimg.get(x,y,z)>0.5) M0++;
                            }
                        }
                  vrl[i] = sqrt((static_cast<RealType>(M0))/(2.0*r*M_PI));

                  // HACKED: Print the values
//                   if (cx == X && cy == Y && cz == Z) cerr<<r<<" "<<vrl[i]<<endl;


                  r += stepr;
                  i++;
                }
              bubble_sort(vrl,nr);
              outputimg.set(cx,cy,cz,vrl[nrmed]);

              // ------------ compute the variance of the radius computation -------------
              if (argc == 6) varianceImg.set( cx,cy,cz, getVariance( vrl, vrl[nrmed], nr ) );

              // HACK: For printing the values
//               if (cx == X && cy == Y && cz == Z) cerr<<"\n\n";

            }
      }
      free(vrl);
      outputimg.save(argv[2], qc::PGM_DOUBLE_BINARY);
      if (argc == 6) varianceImg.save(argv[5], qc::PGM_DOUBLE_BINARY);
    }
  catch ( aol::Exception &ex )
    {
      ex.dump();
    }
}
