#include <cudaInterface.h>
#include <cudaVector.h>

namespace cuda {

template <typename DataType>
GPUMemoryBlock<DataType>::GPUMemoryBlock ( int Size )
  : _pData ( cuda::allocate<DataType> ( Size ) ),
    _size ( Size ) {
}

template <typename DataType>
GPUMemoryBlock<DataType>::~GPUMemoryBlock () {
  cuda::deallocate ( _pData );
}

template <typename DataType>
void GPUMemoryBlock<DataType>::setZero () {
  cuda::memset ( _pData, 0, _size );
}

template <typename DataType>
void GPUMemoryBlock<DataType>::copyFromHost ( const DataType * const PHostMem ) {
  cuda::memcpyHostToGPU ( _pData, PHostMem, _size );
}

template <typename DataType>
void GPUMemoryBlock<DataType>::copyToHost ( DataType *PHostMem ) const {
  cuda::memcpyGPUToHost ( PHostMem, _pData, _size );
}

template class GPUMemoryBlock<double>;
template class GPUMemoryBlock<float>;
template class GPUMemoryBlock<int>;

template <typename RealType>
GPUVector<RealType>::GPUVector ( const GPUVector<RealType> &Vec, const int CopyFlag )
  : GPUMemoryBlock<RealType> ( Vec._size ) {
  assert ( ( CopyFlag == 1 ) || ( CopyFlag == 2 ) );

  if ( CopyFlag == 1 )
    *this = Vec;
  else if ( CopyFlag == 2 )
    this->setZero();
  else
    printf ( "Unexpected CopyFlag value specified.\n" );
}

template <typename RealType>
GPUVector<RealType>& GPUVector<RealType>::operator= ( const GPUVector<RealType> &Vec ) {
  assert ( Vec.size() == this->_size );

  // Beware of self-assignment
  if ( this->_pData != Vec._pData )
    cuda::memcpyGPUToGPU ( this->_pData, Vec._pData, this->_size );

  return *this;
}
template <typename RealType>

GPUVector<RealType>& GPUVector<RealType>::operator+= ( const GPUVector<RealType> &Vec ) {
  return addMultiple ( Vec, 1 );
}
template <typename RealType>

GPUVector<RealType>& GPUVector<RealType>::operator-= ( const GPUVector<RealType> &Vec ) {
  return addMultiple ( Vec, -1 );
}

template <typename RealType>
GPUVector<RealType>& GPUVector<RealType>::operator*= ( const RealType Value ) {
  cublasScale<RealType> ( this->_size, Value, this->_pData );
  return *this;
}

template <typename RealType>
RealType GPUVector<RealType>::operator* ( const GPUVector<RealType> &Vec ) const {
  return cublasScalarProduct<RealType> ( this->_size, this->_pData, Vec._pData );
}

template <typename RealType>
RealType GPUVector<RealType>::norm ( ) const {
  return cublasNormSqr<RealType> ( this->_size, this->_pData ); 
}

template <typename RealType>
GPUVector<RealType>& GPUVector<RealType>::addMultiple ( const GPUVector<RealType> &Vec, RealType Factor ) {
  cublasAddMultiple<RealType> ( this->_size, Factor, Vec._pData, this->_pData ); 
  return *this;
}

template class GPUVector<double>;
template class GPUVector<float>;

void initLibCublas ( ) {
  cublasStatus stat = cublasInit(); 
  if ( stat != CUBLAS_STATUS_SUCCESS ) {
    printf ( "Failed to initialize the CUBLAS library\n" );
  }
}

}
