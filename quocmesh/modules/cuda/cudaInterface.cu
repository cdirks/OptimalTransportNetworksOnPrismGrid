#include <cudaInterface.h>

namespace cuda {

void printCurrentDevice ( ) {
  cudaDeviceProp deviceProp;
  int devID = 0;
  cudaGetDevice ( &devID );
  CUDA_SAFE_CALL( cudaGetDeviceProperties(&deviceProp, devID) );
  printf("Current CUDA device is %s (%d MB total global memory).\n", deviceProp.name, ( deviceProp.totalGlobalMem ) / ( 1024 * 1024 ) );
}

template <>
void cublasScale ( const int Size, const float Value, float *PData ) {
  cublasSscal ( Size, Value, PData, 1 );
}

template <>
void cublasScale ( const int Size, const double Value, double *PData ) {
  cublasDscal ( Size, Value, PData, 1 );
}

template <>
float cublasScalarProduct ( const int Size, const float *PDataA, const float *PDataB ) {
 return cublasSdot ( Size, PDataA, 1, PDataB, 1 );
}

template <>
double cublasScalarProduct ( const int Size, const double *PDataA, const double *PDataB ) {
 return cublasDdot ( Size, PDataA, 1, PDataB, 1 );
}

template <>
float cublasNormSqr ( const int Size, const float *PData ) {
  return cublasSnrm2 ( Size, PData, 1 ); 
}

template <>
double cublasNormSqr ( const int Size, const double *PData ) {
  return cublasDnrm2 ( Size, PData, 1 ); 
}

template <>
void cublasAddMultiple ( const int Size, const float Factor, const float *PDataArg, float *PDataDest ) {
  cublasSaxpy ( Size, Factor, PDataArg, 1, PDataDest, 1); 
}

template <>
void cublasAddMultiple ( const int Size, const double Factor, const double *PDataArg, double *PDataDest ) {
  cublasDaxpy ( Size, Factor, PDataArg, 1, PDataDest, 1); 
}

}
