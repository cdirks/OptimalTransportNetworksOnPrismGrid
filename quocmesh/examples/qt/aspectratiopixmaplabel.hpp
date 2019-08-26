#ifndef __ASPECTRATIOPIXMAPLABEL_HPP
#define __ASPECTRATIOPIXMAPLABEL_HPP

#include <qtIncludes.h>

// Based on AspectRatioPixmapLabel from http://stackoverflow.com/questions/8211982/qt-resizing-a-qlabel-containing-a-qpixmap-while-keeping-its-aspect-ratio
// author http://stackoverflow.com/users/999943/phyatt
class AspectRatioPixmapLabel : public QLabel {
  Q_OBJECT
public:
  explicit AspectRatioPixmapLabel ( QWidget *parent = 0 );
  virtual int heightForWidth ( int width ) const;
  virtual QSize sizeHint() const;
  const QPixmap &getPixmap () const;
signals:

public slots:
  void setPixmap ( const QPixmap & );
  void resizeEvent ( QResizeEvent * );
private:
  QPixmap _pix;
};

#endif // __ASPECTRATIOPIXMAPLABEL_HPP
