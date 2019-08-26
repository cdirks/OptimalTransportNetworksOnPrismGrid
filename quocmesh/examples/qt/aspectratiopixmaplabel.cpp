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

#include "aspectratiopixmaplabel.hpp"

// Based on AspectRatioPixmapLabel from http://stackoverflow.com/questions/8211982/qt-resizing-a-qlabel-containing-a-qpixmap-while-keeping-its-aspect-ratio
// author http://stackoverflow.com/users/999943/phyatt
AspectRatioPixmapLabel::AspectRatioPixmapLabel ( QWidget *parent ) : QLabel ( parent ) {
  this->setMinimumSize ( 1, 1 );
}

void AspectRatioPixmapLabel::setPixmap ( const QPixmap & p ) {
  _pix = p;
  // This takes care of rescaling the new pixmap to display it.
  resizeEvent ( NULL );
}

int AspectRatioPixmapLabel::heightForWidth ( int width ) const {
  return static_cast<qreal> ( _pix.height() * width ) / _pix.width();
}

QSize AspectRatioPixmapLabel::sizeHint() const {
  const int w = this->width();
  return QSize ( w, heightForWidth ( w ) );
}

void AspectRatioPixmapLabel::resizeEvent ( QResizeEvent * /*e*/ ) {
  if ( text().isEmpty() ) {
    QSize size = this->size();
#if QT_VERSION >= 0x050000
    // Make use of the full resolution the display offers.
    _pix.setDevicePixelRatio ( this->devicePixelRatio() );
    size *= this->devicePixelRatio();
#endif
    QLabel::setPixmap ( _pix.scaled ( size, Qt::KeepAspectRatio, Qt::FastTransformation ) );
  }
}

const QPixmap &AspectRatioPixmapLabel::getPixmap () const {
  return _pix;
}

