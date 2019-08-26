#!/bin/sh
cd ..
wget http://downloads.sourceforge.net/project/vtkfox/vtkfox/06064/vtkfox-06064.tar.bz2
tar xvf vtkfox-06064.tar.bz2
cd vtkfox
patch -p0 -i vtkfox.patch
