#ifndef __QCIMVIEWAPP_H
#define __QCIMVIEWAPP_H

#include <qtIncludes.h>
#include <scalarArray.h>
#include <multiArray.h>
#include <gnuplotter.h>

#include "ui_qcImViewApp.h"

/**
 * \brief Very simple QT example class that opens and displays a 2D Quoc array file.
 *
 * \author Berkels
 */
class QuocImageViewApp : public QMainWindow, private Ui::QuocImageViewer {
  Q_OBJECT

  std::string _dirName;
  std::vector<std::string> _imageList;
  unsigned int _currentImageIndex;
  qc::ScalarArray<double, qc::QC_2D> _image;
  int _scaleExponent;

  QLabel *_sbarFileName;
  QLabel *_sbarFileNum;

  void loadImage ( const char *Filename ) {
    imageLabel->setText ( "" );
    QImage qimage;
    if ( aol::fileNameEndsWith ( Filename, ".png" ) ) {
      qc::MultiArray<unsigned char, qc::QC_2D, 3> image;
      // Be carefuly when loading the file. Possibly it can't be read or is not even an image.
      try {
        image.loadPNG( Filename );
      }
      catch ( aol::Exception &el ) {
        el.dump();
        imageLabel->setText( el.getMessage().c_str() );
        image.reallocate ( 1, 1 );
      }
      const int numX = image[0].getNumX();
      const int numY = image[0].getNumY();

      qimage = QImage ( numX, numY, QImage::Format_RGB32 );
      for ( int j = 0; j < numY; ++j )
        for ( int i = 0; i < numX; ++i )
          qimage.setPixel ( i, j, qRgb(image[0].get(i,j),image[1].get(i,j),image[2].get(i,j)) );

      _image.reallocate ( numX, numY );
      for ( int i = 0; i < _image.size();  ++i ) {
        aol::Vec3<double> color;
        for ( int c = 0; c < 3; ++c )
          color[c] = image[c][i];
        _image[i] = color.getMeanValue();
      }
    }
    else {
      // Be carefuly when loading the file. Possibly it can't be read or is not even an image.
      try {
        _image.load ( Filename );
      }
      catch ( aol::Exception &el ) {
        el.dump();
        imageLabel->setText( el.getMessage().c_str() );
        _image.reallocate ( 1, 1 );
      }

      const bool rangeIs8Bit = ( aol::fileNameEndsWith ( Filename, ".pgm" ) );

      const double minVal = rangeIs8Bit ? 0 : _image.getMinValue();
      const double maxVal = rangeIs8Bit ? 255 : _image.getMaxValue();
      quocArrayToQImage ( qimage, minVal, maxVal );
    }

    if ( imageLabel->text().isEmpty() ) {
      imageLabel->setPixmap ( QPixmap::fromImage(qimage) );
      if ( actionAutoContrast->isChecked() )
        on_actionEnhanceContrast_triggered();
    }
  }

  void quocArrayToQImage ( QImage &Qimage, const double MinVal, const double MaxVal ) const {
    const double valScaleFactor = 255 / ( MaxVal - MinVal );
    const int numX = _image.getNumX();
    const int numY = _image.getNumY();

    Qimage = QImage( numX, numY, QImage::Format_Indexed8);

    // QImage doesn't know grayscale as format, so we have to create a grayscale palette.
    for ( int i = 0; i < 256; ++i )
      Qimage.setColor ( i, qRgb(i,i,i) );

    // Convert the values in image to [0,255] and store them in the QImage.
    for ( int j = 0; j < numY; ++j )
      for ( int i = 0; i < numX; ++i )
        Qimage.setPixel ( i, j, static_cast<unsigned char> ( aol::Clamp<double> ( valScaleFactor * ( _image.get ( i, j ) - MinVal ), 0, 255 ) ) );

  }
public:
  QuocImageViewApp ( const char * InputFilename = NULL )
    : _scaleExponent ( 0 ),
      _sbarFileName ( new QLabel ),
      _sbarFileNum ( new QLabel( "0/0" ) ) {
    setupUi ( this );
    statusbar->addPermanentWidget( _sbarFileName, 1 );
    statusbar->addPermanentWidget( _sbarFileNum );
    imageLabel->setAlignment ( Qt::AlignCenter );
    scrollArea->setBackgroundRole ( QPalette::Dark );
    scrollArea->setAlignment ( Qt::AlignCenter );
    scrollArea->setWidget ( imageLabel );
    scrollArea->setWidgetResizable ( true );
    connect ( actionAboutQt, SIGNAL ( triggered() ), qApp, SLOT ( aboutQt() ) );
    loadAndShowImage ( InputFilename );

// Workaround for a Qt5 bug under OS X: https://bugreports.qt-project.org/browse/QTBUG-38256
#if defined ( __APPLE__ ) && ( QT_VERSION >= 0x050000 )
    foreach ( QAction* a, menuFile->actions())
      QObject::connect ( new QShortcut(a->shortcut(), a->parentWidget()), SIGNAL(activated()), a, SLOT(trigger()) );
    foreach ( QAction* a, menuView->actions())
      QObject::connect ( new QShortcut(a->shortcut(), a->parentWidget()), SIGNAL(activated()), a, SLOT(trigger()) );
#endif
  }

  void loadAndShowImage ( const char * InputFilename ) {
    if ( InputFilename != NULL ) {
      loadImage ( InputFilename );
      const QString path = InputFilename;
      initDirectory ( path );
    }
  }

  ~QuocImageViewApp ( ) {
    delete _sbarFileName;
    delete _sbarFileNum;
  }

  bool hasImageLoaded ( ) const {
    return ( _image.size() != 0 );
  }

protected:
  void initDirectory ( const QString Filename ) {
    _dirName = QFileInfo ( Filename ).absolutePath().toLatin1().constData();
    _dirName += "/";
    _imageList.clear();

    std::vector<std::string> dirList;
    aol::createDirectoryListing ( _dirName.c_str(), dirList );

    std::vector<std::string> supportedSuffixes;
    supportedSuffixes.push_back( ".png" );
    supportedSuffixes.push_back( ".dat.bz2" );
    supportedSuffixes.push_back( ".q2bz" );
    supportedSuffixes.push_back( ".pgm" );
    supportedSuffixes.push_back( ".tif" );
    supportedSuffixes.push_back( ".dm3" );
    supportedSuffixes.push_back( ".dm4" );

    for ( unsigned int i = 0; i < dirList.size(); ++i ) {
      // Ignore OS X meta data files.
      if ( ( dirList[i].size() > 1 ) && ( dirList[i].compare ( 0, 2, "._" ) == 0 ) )
        continue;

      for ( unsigned int j = 0; j < supportedSuffixes.size(); ++j ) {
        if ( aol::fileNameEndsWith ( dirList[i].c_str(), supportedSuffixes[j].c_str() ) ) {
          _imageList.push_back ( dirList[i] );
          break;
        }
      }
    }

    _currentImageIndex = 0;

    // Make sure that the filename uses '/' as dir seperator.
    std::string filename = QDir::fromNativeSeparators ( Filename ).toLatin1().constData ();

    for ( unsigned int i = 0; i < _imageList.size(); ++i ) {
      if ( ( _dirName + _imageList[i] ).compare ( filename ) == 0 )
        _currentImageIndex = i;
    }
    updateStatusBar();
  }

  void updateStatusBar ( ) {
    _sbarFileName->setText( ( _dirName + ( ( _imageList.size() > 0 ) ? _imageList [ _currentImageIndex ] : "unknown" ) ).c_str() );
    _sbarFileNum->setText( aol::strprintf ( "%d/%d", _currentImageIndex + 1, _imageList.size() ).c_str() );
  }

  void keyPressEvent ( QKeyEvent * event )  {
    if( event->key() == Qt::Key_F ) {
      if ( this->windowState() & Qt::WindowMaximized )
        this->showNormal();
      else
        this->showMaximized();
    }
    else if( event->key() == Qt::Key_Enter ) {
      if ( this->windowState() & Qt::WindowFullScreen )
        this->showNormal();
      else
        this->showFullScreen();
    }
  }

  void scaleImage ( const int ExponentChange ) {
    if ( actionFitToWindow->isChecked() )
      actionFitToWindow->trigger ( );
    _scaleExponent += ExponentChange;
    const double scaleFactor = std::pow ( 2., _scaleExponent );
    imageLabel->resize ( scaleFactor * _image.getNumX(), scaleFactor * _image.getNumY() );
  }

protected slots:
  void on_actionShowInfo_triggered () {
    std::stringstream message;
    message << "Size " << _image.getNumX() << "x" << _image.getNumY() << endl;
    message << "Value range " << _image.getMinValue() << " to " << _image.getMaxValue() << endl;
    QMessageBox::information(NULL, "Info", message.str ().c_str());
  }

  void on_actionEnhanceContrast_triggered () {
    if ( hasImageLoaded() ) {
      QImage qimage;
      aol::Vec2<double> minMax = _image.getSaturatedMinMaxValue( 0.15 );
      quocArrayToQImage ( qimage, minMax[0], minMax[1] );
      imageLabel->setPixmap ( QPixmap::fromImage(qimage) );
    }
  }

  void on_actionHSVColormap_triggered () {
    if ( hasImageLoaded() ) {
      const aol::RGBColorMap<double> hsvMap ( 0., 1., aol::HSV_BLUE_TO_RED );
      aol::Vec2<double> minMax = _image.getMinMaxValue( );

      const int numX = _image.getNumX();
      const int numY = _image.getNumY();

      QImage qimage = QImage ( numX, numY, QImage::Format_RGB32 );

      for ( int j = 0; j < numY; ++j ) {
        for ( int i = 0; i < numX; ++i ) {
          aol::Vec3<double> color;
          hsvMap.scalarToColor ( ( _image.get ( i, j ) - minMax[0] ) / ( minMax[1] - minMax[0] ), color );
          color *= 255;
          qimage.setPixel ( i, j, qRgb( color[0], color[1], color[2]) );
        }
      }
      imageLabel->setPixmap ( QPixmap::fromImage(qimage) );
    }
  }

  void on_actionNextFile_triggered () {
    _currentImageIndex = ( _currentImageIndex + 1 ) % _imageList.size();
    loadImage ( ( _dirName + _imageList[_currentImageIndex] ).c_str() );
    updateStatusBar();
  }

  void on_actionPrevFile_triggered () {
    _currentImageIndex = ( _currentImageIndex + _imageList.size() - 1 ) % _imageList.size();
    loadImage ( ( _dirName + _imageList[_currentImageIndex] ).c_str() );
    updateStatusBar();
  }

  void on_actionLoad_triggered () {
    const QString path = QFileDialog::getOpenFileName ( this, "Select 2D Quoc array file" );
    if ( path.isEmpty() == false )
      loadAndShowImage ( path.toLatin1().constData() );
  };

  void on_actionPlotHistogram_triggered () {
    if ( aol::runGnuplot ( "--version" ) == false )
      QMessageBox::information(NULL, "Error", "Calling gnuplot failed! Make sure gnuplot is in your search path.");
    else {
      aol::Vector<int> histo;
      _image.createHistogramOfValues ( histo, 256 );
      aol::plotHistogram<double> ( histo, "", true, true );
    }
  };

  void on_actionPlotCenterLine_triggered () {
    if ( aol::runGnuplot ( "--version" ) == false )
      QMessageBox::information(NULL, "Error", "Calling gnuplot failed! Make sure gnuplot is in your search path.");
    else {
      aol::Plotter<double> plotter;
      aol::PlotDataFileHandler<double> plotHandler;
      plotHandler.generateCenterLinePlot ( _image );
      plotter.addPlotCommandsFromHandler ( plotHandler );
      plotter.plotToScreen();
    }
  };

  void on_actionFitToWindow_triggered () {
    if ( actionFitToWindow->isChecked() )
      scrollArea->setWidgetResizable ( true );
    else {
      scrollArea->setWidgetResizable ( false );
      imageLabel->resize ( _image.getNumX(), _image.getNumY() );
    }
  };

  void on_actionAutoContrast_triggered () {
    if ( actionAutoContrast->isChecked() )
      on_actionEnhanceContrast_triggered();
  };

  void on_actionZoomIn_triggered () {
    scaleImage ( 1 );
  };

  void on_actionZoomOut_triggered () {
    scaleImage ( -1 );
  };
};

class QuocViewerApp : public QApplication {
  Q_OBJECT

  std::vector<QuocImageViewApp *> _dialogs;
public:
  QuocViewerApp ( int & argc, char **argv ) : QApplication ( argc, argv ) {
    QuocImageViewApp *dialog = new QuocImageViewApp ( ( argc > 1 ) ? argv[1] : NULL );
    dialog->show();
    _dialogs.push_back ( dialog );
  }
  virtual ~QuocViewerApp() {
    for ( unsigned int i = 0; i < _dialogs.size(); ++i )
      delete _dialogs[i];
  }

  bool notify ( QObject *Receiver, QEvent *Event ) {
    try {
      return QApplication::notify( Receiver, Event );
    }
    catch ( aol::Exception &el ) {
      QMessageBox::information(NULL, "Info", el.getMessage().c_str());
      return false;
    }
  }
protected:
  // This is necessary because OS X doesn't tell an app which file to open with command line arguments, but with an event.
  bool event ( QEvent *Event ) {
    switch ( Event->type() ) {
      case QEvent::FileOpen: {
        QString filename = static_cast<QFileOpenEvent *>(Event)->file();
        if ( ( _dialogs.size() == 1 ) && ( _dialogs[0]->hasImageLoaded() == false ) ) {
          _dialogs[0]->loadAndShowImage ( filename.toLatin1().constData () );
        }
        else {
          QuocImageViewApp *dialog = new QuocImageViewApp( filename.toLatin1().constData () );
          dialog->show();
          _dialogs.push_back ( dialog );
        }
        return true;
      }
      default:
        return QApplication::event ( Event );
    }
  }
};

#endif // __QCIMVIEWAPP_H
