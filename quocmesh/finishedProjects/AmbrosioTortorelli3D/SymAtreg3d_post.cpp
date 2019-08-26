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

//#define FIX_PHASE_REFERENCE
//#define INCLUDE_WEIGHTS
//#define INCLUDE_WEIGHTS2
//#define  __SIX_FRAME__
#define __MIDDLE_FRAME__


#include <aol.h>
#include <parameterParser.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <AmbrosioTortorelli.h>
#include "ambrotelli.h"
#include "SymAtreg3d_atom.h"
#include "utilities.h"


typedef float RType;
//typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_3D, aol::GaussQuadrature<RType,qc::QC_3D,3> > ConfType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType;



void loadAndPrepareImage( const char* Filename,
                          ConfType::ArrayType &Dest,
                          const bool NoScaling)
{
  Dest.load( Filename );
  if( !NoScaling ){
    ConfType::RealType minValue = Dest.getMinValue();
    if( minValue < 0.){ Dest.addToAll( - minValue );}

    ConfType::RealType maxValue = Dest.getMaxValue();
    if( maxValue != 0.){ Dest /= maxValue;}
  }
}

int main( int argc, char **argv ) {

  try {
    char parameterfilename[1024];

    if ( argc > 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
      return EXIT_FAILURE;
    }
    if ( argc == 2 ) {
      sprintf( parameterfilename, "%s",  argv[1] );
    }
    if ( argc == 1 ) {
      if( ConfType::Dim == 3)
        sprintf( parameterfilename, "ParameterFiles/SymRegPost_3D.par" );
      else
        sprintf( parameterfilename, "ParameterFiles/SymRegPost_2D.par" );
    }
    cerr << "Reading parameters from " << parameterfilename << endl;
    aol::ParameterParser parser( parameterfilename );

    //load transformation
    aol::MultiVector<ConfType::RealType>* pPhi=NULL;
    int iDataLevel  = parser.getInt("dataLevel");
    if(qc::LoadMultiVector<ConfType>(iDataLevel,pPhi,"Results/Transform_R2T.dat")==false)
    {
       cerr<<"Loading data fail"<<endl;
       return 0;
    }
    aol::MultiVector<ConfType::RealType>* pPsi=NULL;
    if(qc::LoadMultiVector<ConfType>(iDataLevel,pPsi,"Results/Transform_T2R.dat")==false)
    {
       cerr<<"Loading data fail"<<endl;
       return 0;
    }

    aol::MultiVector<ConfType::RealType> TransNull(ConfType::Dim, (*pPhi)[0].size() );
    TransNull.setZero();

#ifdef __FOUR_FRAME__
    //load images.
    const qc::GridDefinition grid(iDataLevel, ConfType::Dim);

    ConfType::ArrayType Im0( grid);
    ConfType::ArrayType Im1( grid);
    ConfType::ArrayType Im2( grid);
    ConfType::ArrayType Im3( grid);

    char buffer[1024];

    parser.getString("Image0",buffer);
    loadAndPrepareImage(buffer,Im0,false);

    parser.getString("Image1",buffer);
    loadAndPrepareImage(buffer,Im1,false);

    parser.getString("Image2",buffer);
    loadAndPrepareImage(buffer,Im2,false);

    parser.getString("Image3",buffer);
    loadAndPrepareImage(buffer,Im3,false);


    //transform the data.
    ConfType::ArrayType Im03( grid);
    ConfType::ArrayType Im30( grid);

    qc::deformImageWithCoarseDeformation<ConfType>(Im0,grid,grid,Im03,*pPhi);
    qc::deformImageWithCoarseDeformation<ConfType>(Im3,grid,grid,Im30,*pPsi);

    //interpolation.
    ConfType::ArrayType ImEst1(grid);
    ConfType::ArrayType ImEst2(grid);

    (*pPhi)*= 1./3.;
    (*pPsi)*= 1./3.;

    qc::deformImageWithCoarseDeformation<ConfType>(Im0,grid,grid,ImEst1,*pPhi);
    qc::deformImageWithCoarseDeformation<ConfType>(Im3,grid,grid,ImEst2,*pPsi);



    //write the data
    qc::writeImage<ConfType::RealType> (grid, Im0, "results/Image0.pgm");
    qc::writeImage<ConfType::RealType> (grid, Im1, "results/Image1.pgm");
    qc::writeImage<ConfType::RealType> (grid, Im2, "results/Image2.pgm");
    qc::writeImage<ConfType::RealType> (grid, Im3, "results/Image3.pgm");



    //transformed data
    qc::writeImage<ConfType::RealType> (grid, ImEst1, "results/Image1E.pgm");
    qc::writeImage<ConfType::RealType> (grid, ImEst2, "results/Image2E.pgm");


    //diff data
    qc::writeDiffImage<ConfType>(grid,grid,"results/ImageRDiff10.pgm", Im1,Im0,TransNull);
    qc::writeDiffImage<ConfType>(grid,grid,"results/ImageRDiff1.pgm", Im1,Im0,*pPhi);
    qc::writeDiffImage<ConfType>(grid,grid,"results/ImageRDiff2.pgm", Im2,Im3,*pPsi);
    qc::writeDiffImage<ConfType>(grid,grid,"results/ImageRDiff20.pgm", Im2,Im3,TransNull);

#endif

#ifdef __SIX_FRAME__
    //load images.
    const qc::GridDefinition grid(iDataLevel, ConfType::Dim);

    const int interFrameNum = 6;
    const int allFrameNum = interFrameNum+2;

    //load Array for images
    ConfType::ArrayType* pIm[allFrameNum];
    char buffer[1024];
    char strname[1024];
    for(int iFrame = 0; iFrame<allFrameNum ; iFrame++)
    {
       pIm[iFrame] = new ConfType::ArrayType( grid);
       sprintf(strname,"Image%d",iFrame);
       parser.getString(strname,buffer);
       loadAndPrepareImage(buffer,*pIm[iFrame],false);

    }

    //output absolute diff between first and last frame
    qc::writeDiffImage<ConfType>(grid,grid,"results/ImageDiff.pgm", *pIm[0],*pIm[allFrameNum-1],TransNull);

    //forward
    ConfType::ArrayType ImTrans(grid);
    double factor = 0;
    double DiffForward[allFrameNum-1];
    for(int iFrame = 1; iFrame<=(interFrameNum+1); iFrame++)
    {
        //current transform
       factor = (iFrame)/(interFrameNum+1.);
       (*pPhi)*= factor;

       //estimated frame
       //qc::deformImageWithCoarseDeformation<ConfType>(*pIm[0],grid,grid,ImTrans,*pPhi);

       //difference frame
       sprintf(strname,"results/ImageFDiff%d.pgm",iFrame);
       DiffForward[iFrame-1] = qc::writeDiffImage<ConfType>(grid,grid,strname, *pIm[iFrame],*pIm[0],*pPhi);

       //recover phi
       (*pPhi)*= 1./factor;
    }

    //inward
    double DiffBackward[allFrameNum-1];
    for(int iFrame = 0; iFrame<=interFrameNum ; iFrame++)
    {
        //current transform
       factor = (interFrameNum-iFrame+1)/(interFrameNum+1.);
       (*pPsi)*= factor;

       //estimated frame
       //qc::deformImageWithCoarseDeformation<ConfType>(*pIm[allFrameNum-1],grid,grid,ImTrans,*pPsi);

       //difference frame
       sprintf(strname,"results/ImageBDiff%d.pgm",iFrame);
       DiffBackward[iFrame] = qc::writeDiffImage<ConfType>(grid,grid,strname, *pIm[iFrame],*pIm[allFrameNum-1],*pPsi);

       //recover phi
       (*pPsi)*= 1./factor;
    }

    //print out diff

    std::cerr<<std::endl<<"Forward:\t";
    for(int iFrame = 0; iFrame<allFrameNum-1 ; iFrame++)
       std::cerr<<"\t"<<DiffForward[iFrame];

    std::cerr<<std::endl<<"Backward:\t";
    for(int iFrame = 0; iFrame<allFrameNum-1 ; iFrame++)
    {
       std::cerr<<"\t"<<DiffBackward[iFrame];
    }
    std::cerr<<std::endl;


    std::ofstream file_diff;
    file_diff.open("results/diff_file.txt");
    file_diff<<0<<"\t"<<0<<"\t"<<DiffBackward[0]<<std::endl;
    for(int iFrame = 0; iFrame<interFrameNum ; iFrame++)
    {
       file_diff<<iFrame+1<<"\t"<<DiffForward[iFrame]<<"\t"<<DiffBackward[iFrame+1]<<std::endl;
    }
    file_diff<<interFrameNum+1<<"\t"<<DiffForward[interFrameNum]<<"\t"<<0<<std::endl;

    //delete images
     for(int iFrame = 0; iFrame<allFrameNum-1 ; iFrame++)
    {
       delete pIm[iFrame];
    }


#endif //__SIX_FRAME__

#ifdef __MIDDLE_FRAME__
    //load images.
    const qc::GridDefinition grid(iDataLevel, ConfType::Dim);



    //load Array for images
    ConfType::ArrayType ImI(grid), ImP(grid), ImB(grid);
    char buffer[1024];

    parser.getString("ImageI",buffer);
    loadAndPrepareImage(buffer,ImI,false);
    parser.getString("ImageP",buffer);
    loadAndPrepareImage(buffer,ImP,false);
    parser.getString("ImageB",buffer);
    loadAndPrepareImage(buffer,ImB,false);


    *pPhi *= 0.5;
    *pPsi *= 0.5;
    //forward estimation and reverse sestimation
    ConfType::ArrayType ImEstForward(grid),ImEstReverse(grid), ImEstBidirectional(grid);
    qc::deformImageWithCoarseDeformation<ConfType>(ImI,grid,grid,ImEstForward,*pPhi);
    qc::deformImageWithCoarseDeformation<ConfType>(ImP,grid,grid,ImEstReverse,*pPsi);



    //output
    qc::writeImage<ConfType::RealType> (grid, ImEstForward, "results/ImageEstForward.pgm");
    qc::writeImage<ConfType::RealType> (grid, ImEstReverse, "results/ImageEstReverse.pgm");
    ImEstForward += ImEstReverse;
    ImEstForward *=0.5;
     qc::writeImage<ConfType::RealType> (grid, ImEstForward , "results/ImageEstBidirectional.pgm");
#endif
    //clear up
    if(pPhi){delete pPhi;}
    if(pPsi){delete pPsi;}
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
