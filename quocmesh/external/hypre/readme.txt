standard external
interface to the hypre preconditioner library

Note that MPI may cause problems, so hypre needs to be compiled by hand (do not use the one provided elsewhere).

The following is known to work:

unpack hypre-2.0.0
patch the file hypre-2.0.0/src/utilities/hypre_memory.c using 
####
126c126
<    int   size = count*elt_size;
---
>    //   int   size = count*elt_size;
128c128
<    if (size > 0)
---
>    if ((count > 0) && (elt_size > 0))
141c141
<         hypre_OutOfMemory(size);
---
>         hypre_OutOfMemory(count*elt_size);
####
and then compile using
####
cd hypre-2.0.0/src/
CFLAGS=-fPIC CXXFLAGS=-fPIC LDFLAGS=-fPIC ./configure --prefix=/home/prog/hypre-2.0.0-64 --without-MPI --without-blas --without-lapack

sed -i "s/\${CC} -c dlamch.c/\${CC} -fPIC -c dlamch.c/" lapack/Makefile

make
make install
cd /home/prog/hypre-2.0.0-64/lib
for F in lib*.a ; do ar x $F ; done
g++ -shared *.o -o libHYPRE.so -Xlinker -soname=`pwd`/libHYPRE.so
rm -f *.o *.a
chown root:root /home/prog/hypre-2.0.0-64 -R
chmod a+rX /home/prog/hypre-2.0.0-64 -R
cd /home/prog/hypre-2.0.0-64
mv hypre-2.0.0/docs ./
rm -rf hypre-2.0.0 hypre-2.0.0.tar.bz2
####