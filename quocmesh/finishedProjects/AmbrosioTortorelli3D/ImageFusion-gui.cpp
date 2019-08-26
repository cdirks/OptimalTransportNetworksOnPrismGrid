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

/*********************************************
*  File:   ImageFusion-gui.cpp
*  Date:   18.04.2006
*  Brief:  Combine the images for visualization
*
***********************************************/

#ifdef USE_EXTERNAL_FOX

#include <foxIncludes.h>

#define PATTEN_2D "*.pgm"
#define PATTEN_3D "*.raw"
#define PATTEN_TRANSFORM "*.dat"

#define GUI_PRESENT

#include <scalarArray.h>
#include <configurators.h>
#include <AmbrosioTortorelli.h>
#include "ambrotelli.h"
#include "utilities.h"

typedef float RType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType;
// Event Handler Object

class ImageFusionPanel: public FXMainWindow {

  // Macro for class hierarchy declarations
  FXDECLARE(ImageFusionPanel)

  typedef  ConfType::ArrayType ArrayType;
  typedef  aol::MultiVector< ConfType::RealType> MultiArrayType;
private:

  const int _numberOfUpdateableImages;
  int _levelOfDeformation;
  std::vector<FXColor*> _imageDataVector;
  //FXint _width, _height, _npixels;
  FXImageUpdater<RType, ConfType::Dim> *_imageUpdater;

  FXTextField* m_TextInput1 ;
  FXTextField* m_TextInput2;
  FXTextField* m_TextDeform;
  FXTextField* m_TextResult;

  FXSpinner* m_Gridspin;

  enum {  NUMFUSIONSTYLE=4 };
  enum STYLE {  STYLE_CHECKBOARD, STYLE_STRIP, STYLE_ABS, STYLE_TRANFORM };


  FXRadioButton* m_pRBStyle[NUMFUSIONSTYLE];
  FXTabBook* tabbookPara;
  FXTabItem* m_pTabPara[NUMFUSIONSTYLE];

  FXCheckButton* m_CKState;
  FXTextField* m_TextCKSize;
  FXTextField* m_TextStripSize;

  STYLE m_Style;

  qc::GridDefinition* m_pGrid;
  ArrayType *m_Image1;
  ArrayType *m_Image2;
  ArrayType *m_Output;


  MultiArrayType* m_Trans;


protected:
  ImageFusionPanel():
  _numberOfUpdateableImages(0),
  _levelOfDeformation(0),
  m_Style(STYLE_CHECKBOARD),
  m_pGrid(NULL),
  m_Image1(NULL),
  m_Image2(NULL),
  m_Output(NULL),
  m_Trans(NULL) {}

public:

  // Message handlers
  long onRunClick(FXObject*,FXSelector,void*);
  long onAboutClick(FXObject*,FXSelector,void*);
  long onStyleRBClick(FXObject*,FXSelector,void*);

  long onImage1OpenClick(FXObject*,FXSelector,void*);
  long onImage2OpenClick(FXObject*,FXSelector,void*);
  long onDeformOpenClick(FXObject*,FXSelector,void*);
  long onResultOpenClick(FXObject*,FXSelector,void*);

  long onLoad1Click(FXObject*,FXSelector,void*);
  long onLoad2Click(FXObject*,FXSelector,void*);
  long onLoadDeformClick(FXObject*,FXSelector,void*);
  long onSaveResultClick(FXObject*,FXSelector,void*);

  long onGridChange(FXObject*,FXSelector,void*);
  long onCKDeform(FXObject*,FXSelector,void*);
  long onUpdateImages(FXObject*,FXSelector,void*);
public:

  // Messages for our class
  enum{
    ID_RUN=FXMainWindow::ID_LAST,
    ID_OPENIN1,
    ID_OPENIN2,
    ID_OPENDEFORM,
    ID_OPENRESULT,
    ID_LOAD1,
    ID_LOAD2,
    ID_LOADDEFORM,
    ID_SAVERESULTS,
    ID_RKCHANGE,
    ID_TABBOOKINTERACTIVE,
    ID_GRIDCHANGE,
    ID_ABOUT,
    ID_SWITCH,
    ID_UPDATE_IMAGES,
    ID_LAST
  };

public:

  // constructor
  ImageFusionPanel(FXApp* a);

  // Initialize
  virtual void create();

  virtual ~ImageFusionPanel();

public:

protected:
  bool openFileFromDiag(const FXString& initDir, const FXString& pattern, FXString& outFile);
  FXString getInitDir(FXString curDir, FXString defaultDir1, FXString defaultDir2 );
  bool FusionUpdate();
};


// Message Map for the Scribble Window class
FXDEFMAP(ImageFusionPanel) ImageFusionPanelMap[]={
                                                   //________Message_Type_____________________ID____________Message_Handler_______
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_RUN,  ImageFusionPanel::onRunClick),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_RKCHANGE,  ImageFusionPanel::onStyleRBClick),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_ABOUT,  ImageFusionPanel::onAboutClick),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_OPENIN1,  ImageFusionPanel::onImage1OpenClick),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_OPENIN2,  ImageFusionPanel::onImage2OpenClick),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_OPENDEFORM,  ImageFusionPanel::onDeformOpenClick),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_OPENRESULT,  ImageFusionPanel::onResultOpenClick),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_LOAD1,  ImageFusionPanel::onLoad1Click),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_LOAD2,  ImageFusionPanel::onLoad2Click),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_LOADDEFORM,  ImageFusionPanel::onLoadDeformClick),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_SAVERESULTS,  ImageFusionPanel::onSaveResultClick),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_GRIDCHANGE,  ImageFusionPanel::onGridChange),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_RKCHANGE,  ImageFusionPanel::onCKDeform),
                                                   FXMAPFUNC(SEL_COMMAND,            ImageFusionPanel::ID_UPDATE_IMAGES,  ImageFusionPanel::onUpdateImages)
                                                 };



// Macro for the ScribbleApp class hierarchy implementation
FXIMPLEMENT(ImageFusionPanel,FXMainWindow,ImageFusionPanelMap,ARRAYNUMBER(ImageFusionPanelMap))


// Create and initialize
void ImageFusionPanel::create() {

  // Create the windows
  FXMainWindow::create();

  // Make the main window appear
  show(PLACEMENT_SCREEN);
}

// Create and initialize
ImageFusionPanel::ImageFusionPanel(FXApp *a)
    : FXMainWindow(a,"Image Fusion GUI",NULL,NULL,DECOR_ALL),
    _numberOfUpdateableImages(3),
    _levelOfDeformation(0),
    _imageDataVector(_numberOfUpdateableImages),
    m_Style(STYLE_CHECKBOARD),
    m_pGrid(NULL),
    m_Image1(NULL),
    m_Image2(NULL),
    m_Output(NULL),
    m_Trans(NULL) {

  FXHorizontalFrame *TopFrame = new FXHorizontalFrame(this);
  // top frame
  FXVerticalFrame *contentsFrame=new FXVerticalFrame(TopFrame);


  // loading frame
  FXMatrix* fileMatrix = new FXMatrix(contentsFrame, 4);

  new FXButton(fileMatrix,"...",NULL,this,ID_OPENIN1,LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X);
  new FXButton(fileMatrix,"...",NULL,this,ID_OPENIN2,LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X);
  new FXButton(fileMatrix,"...",NULL,this,ID_OPENDEFORM,LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X);
  new FXButton(fileMatrix,"...",NULL,this,ID_OPENRESULT,LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X);

  m_TextInput1 = new FXTextField(fileMatrix,35,NULL,0,FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP|LAYOUT_FIX_HEIGHT,0,0, 150,20);
  m_TextInput2 = new FXTextField(fileMatrix,35,NULL,0,FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP|LAYOUT_FIX_HEIGHT,0,0,150,20);
  m_TextDeform = new FXTextField(fileMatrix,35,NULL,0,FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP|LAYOUT_FIX_HEIGHT,0,0,150,20);
  m_TextResult = new FXTextField(fileMatrix,35,NULL,0,FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP|LAYOUT_FIX_HEIGHT,0,0,150,20);

  new FXButton(fileMatrix,"&Load ImR",NULL,this,ID_LOAD1,LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X);
  new FXButton(fileMatrix,"&Load ImT",NULL,this,ID_LOAD2,LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X);
  new FXButton(fileMatrix,"&Load Phi",NULL,this,ID_LOADDEFORM,LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X);
  new FXButton(fileMatrix,"&Save R*T",NULL,this,ID_SAVERESULTS,LAYOUT_CENTER_X | BUTTON_NORMAL | LAYOUT_FILL_X);

  new FXLabel(fileMatrix,"Grid:",NULL,JUSTIFY_CENTER_X|LAYOUT_FILL_X);
  m_Gridspin=new FXSpinner(fileMatrix,3,this,ID_GRIDCHANGE,FRAME_SUNKEN|FRAME_THICK|SPIN_NOMAX|LAYOUT_CENTER_Y|LAYOUT_RIGHT);

  //new FXLabel(fileMatrix,"On",NULL,JUSTIFY_CENTER_X|LAYOUT_FILL_X);
  //m_Gridspin=new FXSpinner(fileMatrix,3,InputType_tgt,ID_SWITCH,FRAME_SUNKEN|FRAME_THICK|SPIN_NOMAX|LAYOUT_CENTER_Y|LAYOUT_RIGHT);
  m_CKState = new FXCheckButton(fileMatrix,FXString("Deform"),this,ID_SWITCH);
  FXHorizontalFrame *optionFrame = new FXHorizontalFrame(contentsFrame);

  FXGroupBox* pGroupBox = new FXGroupBox(optionFrame,"Fusion Style",GROUPBOX_TITLE_RIGHT|FRAME_RIDGE|LAYOUT_FILL_X|LAYOUT_FILL_Y);
  m_pRBStyle[0] = new FXRadioButton(pGroupBox,("CheckBox"), this, ID_RKCHANGE);
  m_pRBStyle[1] = new FXRadioButton(pGroupBox,("Strip"), this, ID_RKCHANGE);
  m_pRBStyle[2] = new FXRadioButton(pGroupBox,("Difference"), this, ID_RKCHANGE);
  m_pRBStyle[3] = new FXRadioButton(pGroupBox,("No Fusion"), this, ID_RKCHANGE);
  FXHorizontalFrame *tabFrame = new FXHorizontalFrame(optionFrame,FRAME_SUNKEN);
  tabbookPara = new FXTabBook(tabFrame,this);
  //parameter of check box
  m_pTabPara[0] = new FXTabItem(tabbookPara,("CK"),NULL);
  FXVerticalFrame* pFrameCheckboard = new FXVerticalFrame(tabbookPara, FRAME_SUNKEN);
  new FXLabel(pFrameCheckboard,"Setting:",NULL,JUSTIFY_CENTER_X|LAYOUT_FILL_X);
  // loading frame
  FXMatrix* CheckboardMatrix = new FXMatrix(pFrameCheckboard, 1);
  new FXLabel(CheckboardMatrix,"Width of box:",NULL,JUSTIFY_CENTER_X|LAYOUT_FILL_X);
  m_TextCKSize = new FXTextField(CheckboardMatrix,9,this,ID_UPDATE_IMAGES,FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP|LAYOUT_FIX_HEIGHT,0,0, 150,20);


  m_pTabPara[1] = new FXTabItem(tabbookPara,("St"),NULL);
  FXVerticalFrame* pFrameStrip = new FXVerticalFrame(tabbookPara, FRAME_SUNKEN);
  new FXLabel(pFrameStrip,"Setting:",NULL,JUSTIFY_CENTER_X|LAYOUT_FILL_X);
  FXMatrix* StripMatrix = new FXMatrix(pFrameStrip, 1);
  new FXLabel(StripMatrix,"Width of Strip:",NULL,JUSTIFY_CENTER_X|LAYOUT_FILL_X);
  m_TextStripSize = new FXTextField(StripMatrix,9,this,ID_UPDATE_IMAGES,FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP|LAYOUT_FIX_HEIGHT,0,0, 150,20);

  m_pTabPara[2] = new FXTabItem(tabbookPara,("Di"),NULL);
  FXVerticalFrame* pFrameABS = new FXVerticalFrame(tabbookPara, FRAME_SUNKEN);
  new FXLabel(pFrameABS,"Setting:",NULL,JUSTIFY_CENTER_X|LAYOUT_FILL_X);

  m_pTabPara[3] = new FXTabItem(tabbookPara,("No"),NULL);
  FXVerticalFrame* pFrameTrans = new FXVerticalFrame(tabbookPara, FRAME_SUNKEN);
  new FXLabel(pFrameTrans,"Setting:",NULL,JUSTIFY_CENTER_X|LAYOUT_FILL_X);


  //initialization
  m_pRBStyle[0]->setCheck(true);
  m_pRBStyle[1]->setCheck(false);
  m_pRBStyle[2]->setCheck(false);
  m_pRBStyle[3]->setCheck(false);


  m_TextCKSize->setText("32");
  m_TextStripSize->setText("32");

  tabbookPara->setCurrent(0);

  m_Gridspin->setIncrement(1);
  m_Gridspin->setRange(1,30);
  m_Gridspin->setValue(8);

  m_CKState->getCheck();

  m_pGrid = new qc::GridDefinition(m_Gridspin->getValue(),ConfType::Dim);
  m_Image1 = new ArrayType(*m_pGrid);
  m_Image2 = new ArrayType(*m_pGrid);
  m_Output = new ArrayType(*m_pGrid);


  //m_TextInput1->setText("P:/local/project/input/2D/center3.pgm");
  //m_TextInput2->setText("P:/local/project/input/2D/center4.pgm");
  //m_TextDeform->setText("P:/local/project/quocmesh/projects/ambrosio_tortorelli3D/results/Transform_T2R.dat");
  //m_TextResult->setText("P:/local/tmp/out.pgm");

  //m_Trans->setAll(.1);
  //qc::SaveMultiVector<ConfType>(*this->m_pGrid,*this->m_Trans,m_TextDeform->getText().text());



  for( int i = 0; i < _numberOfUpdateableImages; i++) {
    FXCALLOC(&_imageDataVector[i],FXColor,m_pGrid->getNumberOfNodes());
  }


  FXMatrix* imageMatrix = new FXMatrix(TopFrame, 1);
  std::vector<FXVerticalFrame*> columns(_numberOfUpdateableImages);
  for( int i = 0; i < _numberOfUpdateableImages; i++) {
    columns[i] = new FXVerticalFrame(imageMatrix);
  }


  std::vector<FXImage*> imageVector(_numberOfUpdateableImages);
  std::vector<FXImageFrame*> imageFrameVector(_numberOfUpdateableImages);
  for( int i = 0; i < _numberOfUpdateableImages; i++) {
    imageVector[i] = new FXImage(getApp(), _imageDataVector[i], 0, m_pGrid->getWidth(), m_pGrid->getWidth());
    imageFrameVector[i] = new FXImageFrame(columns[i], imageVector[i]);
  }

  _imageUpdater = new FXImageUpdater<ConfType::RealType, ConfType::Dim>( getApp(), imageFrameVector, imageVector, _imageDataVector, m_pGrid->getWidth() );

}
ImageFusionPanel::~ImageFusionPanel() {
  for( int i = 0; i < _numberOfUpdateableImages; i++) {
    FXFREE(&_imageDataVector[i]);
  }
  if(_imageUpdater)
    delete _imageUpdater;

  if(m_Image1) {
    delete m_Image1;
    m_Image1=NULL;
  }

  if(m_Image2) {
    delete m_Image2;
    m_Image2=NULL;
  }
  if(m_Trans) {
    delete m_Trans;
    m_Trans=NULL;
  }

  if(m_Output) {
    delete m_Output;
    m_Output=NULL;
  }

  if(m_pGrid) {
    delete m_pGrid;
    m_pGrid=NULL;
  }
}
long ImageFusionPanel::onGridChange(FXObject*,FXSelector,void*) {
  int newGridDepth= m_Gridspin->getValue();
  if(m_pGrid->getGridDepth()!=newGridDepth) {
    cerr<<"Grid level changeg to " << newGridDepth << endl;
    if(m_pGrid) {
      delete m_pGrid;
      m_pGrid=NULL;
    }
    m_pGrid=new qc::GridDefinition(newGridDepth,ConfType::Dim);
    if(m_Image1) {
      delete m_Image1;
      m_Image1=NULL;
    }
    m_Image1 = new ArrayType(*m_pGrid);


    if(m_Image2) {
      delete m_Image2;
      m_Image2=NULL;
    }
    m_Image2 = new ArrayType(*m_pGrid);

    if(m_Output) {
      delete m_Output;
      m_Output=NULL;
    }
    m_Output = new ArrayType(*m_pGrid);

    _imageUpdater->updateImage( *this->m_pGrid, *m_Image1, 0 );
    _imageUpdater->updateImage( *this->m_pGrid, *m_Image2, 1 );
    _imageUpdater->updateImage( *this->m_pGrid, *m_Output, 2 );
  }
  return 1;
}

long ImageFusionPanel::onImage1OpenClick(FXObject*,FXSelector,void*) {
  FXString iniDir = "";
  FXString outFile = "";
  bool fileFound=false;
  if(ConfType::Dim==qc::QC_3D)
    fileFound = openFileFromDiag(iniDir, PATTEN_3D, outFile);
  else
    fileFound = openFileFromDiag(iniDir, PATTEN_2D, outFile);

  if(fileFound)
    m_TextInput1->setText(outFile);

  if(m_TextInput2->getText().trim()=="")
    m_TextInput2->setText(outFile);

  if(m_TextResult->getText().trim()=="") {
    FXString resultFile = m_TextInput2->getText().trim().rbefore(FXchar('\\'));
    if(ConfType::Dim==qc::QC_3D)
      resultFile = resultFile+"\\out.raw";
    else
      resultFile =resultFile+"\\out.pgm";

    m_TextResult->setText(resultFile);
  }
  return 1;
}
long ImageFusionPanel::onImage2OpenClick(FXObject*,FXSelector,void*) {

  FXString iniDir = this->getInitDir(m_TextInput2->getText().trim(), m_TextInput1->getText().trim(), "" );
  FXString outFile = "";
  bool fileFound=false;
  if(ConfType::Dim==qc::QC_3D)
    fileFound = openFileFromDiag(iniDir, PATTEN_3D, outFile);
  else
    fileFound = openFileFromDiag(iniDir, PATTEN_2D, outFile);

  if(fileFound)
    m_TextInput2->setText(outFile);


  return 1;
}
long ImageFusionPanel::onDeformOpenClick(FXObject*,FXSelector,void*) {
  FXString iniDir = this->getInitDir(m_TextDeform->getText().trim(), m_TextInput1->getText().trim(), "" );
  FXString outFile = "";
  bool fileFound=false;
  fileFound = openFileFromDiag(iniDir, PATTEN_TRANSFORM, outFile);

  if(fileFound)
    m_TextDeform->setText(outFile);
  return 1;
}
long ImageFusionPanel::onResultOpenClick(FXObject*,FXSelector,void*) {
  FXString iniDir =  this->getInitDir(m_TextResult->getText().trim(), m_TextInput1->getText().trim(), "" );
  FXString outFile = "";
  bool fileFound=false;
  if(ConfType::Dim==qc::QC_3D)
    fileFound = openFileFromDiag(iniDir, PATTEN_3D, outFile);
  else
    fileFound = openFileFromDiag(iniDir, PATTEN_2D, outFile);

  if(fileFound)
    m_TextResult->setText(outFile);
  return 1;
}


long ImageFusionPanel::onLoad1Click(FXObject*,FXSelector,void*) {
  FXString imageDir=m_TextInput1->getText().trim();
  if(imageDir=="")
    return 1;

  m_Image1->load(imageDir.text());

  ConfType::RealType minValue = m_Image1->getMinValue();
  if( minValue < 0.)
    (*m_Image1).addToAll( -minValue );

  ConfType::RealType maxValue = m_Image1->getMaxValue();
  if( maxValue != 0.)
    (*m_Image1) /= maxValue;

  _imageUpdater->updateImage( *this->m_pGrid, *m_Image1, 0 );
  return 1;
}
long ImageFusionPanel::onLoad2Click(FXObject*,FXSelector,void*) {
  FXString imageDir=m_TextInput2->getText().trim();
  if(imageDir=="")
    return 1;

  m_Image2->load(imageDir.text());

  ConfType::RealType minValue = m_Image2->getMinValue();
  if( minValue < 0.)
    (*m_Image2).addToAll( -minValue );

  ConfType::RealType maxValue = m_Image2->getMaxValue();
  if( maxValue != 0.)
    (*m_Image2) /= maxValue;
  _imageUpdater->updateImage( *this->m_pGrid, *m_Image2, 1 );
  return 1;
}
long ImageFusionPanel::onLoadDeformClick(FXObject*,FXSelector,void*) {
  FXString imageDir=m_TextDeform->getText().trim();
  if(imageDir=="")
    return 1;
  if( m_Trans ) {
    delete m_Trans;
    m_Trans=NULL;
  }
  if(!qc::LoadMultiVector<ConfType>(_levelOfDeformation, m_Trans, imageDir.text())) {
    cerr<<"Loading of deformation failed."<<endl;
  }
  return 1;
}
long ImageFusionPanel::onCKDeform(FXObject*,FXSelector,void*) {
  if(this->m_CKState->getCheck()) {
    this->m_CKState->setCheck(false);
  } else {
    this->m_CKState->setCheck(true);
  }
  FusionUpdate();
  return 1;
}
long ImageFusionPanel::onSaveResultClick(FXObject*,FXSelector,void*) {
  FXString outFile=     this->m_TextResult->getText().trim().rbefore(FXchar('.'));

  qc::writeImage<ConfType::RealType> ( *this->m_pGrid, *this->m_Output, outFile.text() );
  return 1;
}

long ImageFusionPanel::onRunClick(FXObject*,FXSelector,void*) {

  return 1;
}
long ImageFusionPanel::onStyleRBClick(FXObject* sender,FXSelector,void*) {

  for(int i=0;i<NUMFUSIONSTYLE;i++) {
    m_pRBStyle[i]->setCheck(false);
  }
  sender->handle(this,FXSEL(SEL_COMMAND,ID_CHECK),NULL);
  for(int i=0;i<NUMFUSIONSTYLE;i++) {
    if( m_pRBStyle[i]->getCheck()) {
      m_Style = static_cast<STYLE>(i);
      break;
    }
  }
  m_pTabPara[int(m_Style)]->setFocus();
  tabbookPara->setCurrent(int(m_Style));
  FusionUpdate();
  return 1;
}

long ImageFusionPanel::onAboutClick(FXObject*,FXSelector,void*) {
  FXDialogBox about(this,tr("About FUSION GUI"),DECOR_TITLE|DECOR_BORDER,0,0,0,0, 0,0,0,0, 0,0);

  FXColor* data;
  FXuchar *pp;
  FXint width, height, npixels;

  data = NULL;
  width = 257;
  height = 257;
  npixels=width*height;

  FXCALLOC(&data,FXColor,npixels);
  pp = (FXuchar*)data;

  for( int i = 0; i < height; i++) {
    for( int j = 0; j < width; j++,pp+=4) {
      unsigned char temp = static_cast<unsigned char>( static_cast<double>(i*j)/static_cast<double>((width-1)*(height-1))*255.0);
      pp[0] = temp;
      pp[1] = temp;
      pp[2] = temp;
      pp[3] = 255;
    }
  }
  FXImage test(getApp(), data, 0, width, height);

  FXHorizontalFrame aboutContents(&about,LAYOUT_SIDE_TOP|LAYOUT_FILL_X|LAYOUT_FILL_Y,0,0,0,0,0,0,0,0);
  new FXImageFrame(&aboutContents,&test);
  FXVerticalFrame side (&aboutContents,LAYOUT_SIDE_RIGHT|LAYOUT_FILL_X|LAYOUT_FILL_Y,0,0,0,0, 10,10,10,10, 0,0);
  new FXLabel(&side,"Ambrosio Tortorelli Registration GUI",NULL,JUSTIFY_LEFT|ICON_BEFORE_TEXT|LAYOUT_FILL_X);
  new FXHorizontalSeparator(&side,SEPARATOR_LINE|LAYOUT_FILL_X);
  new FXLabel(&side,"This tool registers two images by an extension of the Ambrosio--Tortorelli approximation of the Mumford--Shah Model",NULL,JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_Y);
  FXButton button(&side,tr("&OK"),NULL,&about,FXDialogBox::ID_ACCEPT,BUTTON_DEFAULT|FRAME_RAISED|FRAME_THICK|LAYOUT_RIGHT,0,0,0,0,32,32,2,2);
  button.setFocus();
  about.execute(PLACEMENT_OWNER);
  FXFREE(&data);
  return 1;
}
bool ImageFusionPanel::openFileFromDiag(const FXString& initDir, const FXString& pattern, FXString& outFile) {
  FXFileDialog FileOpen(this,"Open");
  FileOpen.setDirectory (initDir);
  FileOpen.setPattern(pattern);
  if(FileOpen.execute()) {
    outFile=FileOpen.getFilename ();
    return true;
  } else {
    outFile="";
    return false;
  }

  return true;
}

FXString ImageFusionPanel::getInitDir(FXString curDir, FXString defaultDir1, FXString defaultDir2 ) {
  FXString iniDir;
  if(curDir=="") {
    if(defaultDir1=="") {
      iniDir=defaultDir2;
    } else {
      iniDir=defaultDir1.rbefore(FXchar('\\'));
    }
  } else {
    iniDir = curDir.rbefore(FXchar('\\'));
  }
  return iniDir;
}
bool ImageFusionPanel::FusionUpdate() {
  FXString strDeform= m_TextDeform->getText().trim();

  if(m_Image1==NULL||m_Image2==NULL||m_Output==NULL||m_pGrid==NULL) {
    cerr<<"Invalid member."<<endl;
    return false;
  }


  qc::Array<ConfType::RealType> im1(*this->m_Image1, *this->m_pGrid, aol::FLAT_COPY);
  qc::Array<ConfType::RealType> im2(*this->m_Image2, *this->m_pGrid, aol::FLAT_COPY);
  qc::Array<ConfType::RealType> imOut(*this->m_Output, *this->m_pGrid, aol::FLAT_COPY);
  qc::Array<ConfType::RealType> imDeform(*this->m_pGrid);

  if(m_CKState->getCheck()) {
    if(m_Trans==NULL) {
      cerr<<"No Transform."<<endl;
      return false;
    }
    qc::GridDefinition gridOfDeformation(_levelOfDeformation, ConfType::Dim);
    qc::deformImageWithCoarseDeformation<ConfType>(im2, *this->m_pGrid, gridOfDeformation, imDeform, *this->m_Trans);
  } else {
    imDeform = im2;
  }

  qc::DataGenerator<ConfType> generator ( *this->m_pGrid );

  switch (m_Style) {
  case STYLE_CHECKBOARD:
    generator.generateCheckView ( imOut, im1, imDeform, FXIntVal(m_TextCKSize->getText()), false );
    break;
  case STYLE_STRIP:
    generator.generateCheckView ( imOut, im1, imDeform, FXIntVal(m_TextStripSize->getText()), true );
    break;
  case STYLE_ABS:
    qc::ComputeDifference<ConfType>( *this->m_pGrid,im1,imDeform,imOut);
    break;
  case STYLE_TRANFORM:
    imOut = imDeform;
    //qc::deformImageWithCoarseDeformation<ConfType>(im2, *this->m_pGrid,*this->m_pGrid, imOut, *this->m_Trans);
    break;
  default:
    break;
  }


  _imageUpdater->updateImage( *this->m_pGrid, *m_Output, 2 );
  return 1;
}

long ImageFusionPanel::onUpdateImages(FXObject*,FXSelector,void*){
  FusionUpdate();
  return 1;
}

int main( int argc, char **argv ) {
  try {
    // Make application
    FXApp application("ImageFusion","Fox");

    // Start app
    application.init(argc,argv);

    // Scribble window
    new ImageFusionPanel(&application);

    // Create the application's windows
    application.create();

    // Run the application
    application.run();

  } catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  return 0;
}

#else

#include <aol.h>

int main ( int, char** ) {
  cerr << "Needs to be compiled with -DUSE_EXTERNAL_FOX\n";
  return 0;
}

#endif