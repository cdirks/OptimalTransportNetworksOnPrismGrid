#include <cudaInterface.h>
#include <cudaMatrix.h>

namespace cuda {

template <typename RealType>
GPUFullMatrix<RealType>::GPUFullMatrix ( const int NumRows, const int NumCols )
  : _data ( NumRows*NumCols ),
    _numRows ( NumRows ),
    _numCols ( NumCols ) {
}

template <typename RealType>
void GPUFullMatrix<RealType>::setZero () {
  _data.setZero();
}

template <typename RealType>
void GPUFullMatrix<RealType>::apply ( const GPUVector<RealType> &Arg, GPUVector<RealType> &Dest ) const {
  cublasSgemv ( 'n', _numRows, _numCols, 1, _data.getDataPointer(), _numRows, Arg.getDataPointer(), 1, 0, Dest.getDataPointer(), 1);
}

template <typename RealType>
void GPUFullMatrix<RealType>::applyAdd ( const GPUVector<RealType> &Arg, GPUVector<RealType> &Dest ) const {
  cublasSgemv ( 'n', _numRows, _numCols, 1, _data.getDataPointer(), _numRows, Arg.getDataPointer(), 1, 1, Dest.getDataPointer(), 1);
}

template <typename RealType>
void GPUFullMatrix<RealType>::copyFromHost ( const RealType * const PHostMem, const bool Flip ) {
  if ( Flip == false )
    _data.copyFromHost ( PHostMem );
  else {
    RealType *pBuffer = new RealType [ _numRows * _numCols ];
    for ( int i = 0; i < _numRows; ++i )
      for ( int j = 0; j < _numCols; ++j )
        pBuffer[i + j*_numRows] = PHostMem[i*_numCols + j];
    _data.copyFromHost ( pBuffer );
    delete[] pBuffer;
  }
}

template class GPUFullMatrix<float>;

template <typename RealType>
void GPUDiagonalMatrix<RealType>::apply ( const GPUVector<RealType> &Arg, GPUVector<RealType> &Dest ) const {
  // In order to use cubals to apply our diagonal vector as matrix, we interpret this vector as band matrix
  // with zero sub- and superdiagonals.
  cublasSgbmv ('n', getNumRows(), getNumCols(), 0, 0, 1, this->getDataPointer(), 1, Arg.getDataPointer(), 1, 0, Dest.getDataPointer(), 1);
}

template class GPUDiagonalMatrix<float>;

/**
 * \author Berkels
 */
template <typename RealType>
__device__  void cache2DArray ( const RealType* Arg, RealType* ArgCached, const int ix, const int iy, const int Width, const int CacheSizeX, const int row ) {

  const int tx = threadIdx.x;
  const int ty = threadIdx.y;

  ArgCached[getGlobalIndex ( tx + 1, ty + 1, CacheSizeX )] = Arg[row];

  if ( tx == blockDim.x-1 )
    ArgCached[ getGlobalIndex ( tx + 2, ty + 1, CacheSizeX ) ] = ( ix < Width - 1 ) ? Arg [ getGlobalIndex ( ix + 1, iy, Width ) ] : 0;
  else if ( tx == 0 )
    ArgCached[ getGlobalIndex ( tx, ty + 1, CacheSizeX ) ] = ( ix > 0 ) ? Arg [ getGlobalIndex ( ix - 1, iy, Width ) ] : 0;
  if ( ty == blockDim.y-1 )
    ArgCached[ getGlobalIndex ( tx + 1, ty + 2, CacheSizeX ) ] = ( iy < Width - 1 ) ? Arg [ getGlobalIndex ( ix, iy + 1, Width ) ] : 0;
  else if ( ty == 0 )
    ArgCached[ getGlobalIndex ( tx + 1, ty, CacheSizeX ) ] = ( iy > 0 ) ? Arg [ getGlobalIndex ( ix, iy - 1, Width ) ] : 0;

  if ( ( tx == 0 ) && ( ty == 0 ) )
    ArgCached[ getGlobalIndex ( tx, ty, CacheSizeX ) ] = ( ( ix > 0 ) && ( iy > 0 ) ) ? Arg [ getGlobalIndex ( ix - 1, iy - 1, Width ) ] : 0;
  else if ( ( tx == 0 ) && ( ty == blockDim.y-1 ) )
    ArgCached[ getGlobalIndex ( tx, ty + 2, CacheSizeX ) ] = ( ( ix > 0 ) && ( iy < Width - 1 ) ) ? Arg [ getGlobalIndex ( ix - 1, iy + 1, Width ) ] : 0;
  else if ( ( tx == blockDim.x-1 ) && ( ty == 0 ) )
    ArgCached[ getGlobalIndex ( tx + 2, ty, CacheSizeX ) ] = ( ( ix < Width - 1 ) && ( iy > 0 ) ) ? Arg [ getGlobalIndex ( ix + 1, iy - 1, Width ) ] : 0;
  else if ( ( tx == blockDim.x-1 ) && ( ty == blockDim.y-1 ) )
    ArgCached[ getGlobalIndex ( tx + 2, ty + 2, CacheSizeX ) ] = ( ( ix < Width - 1 ) && ( iy < Width - 1 ) ) ? Arg [ getGlobalIndex ( ix + 1, iy + 1, Width ) ] : 0;
}

/**
 * \author Berkels
 */
template <typename RealType>
__device__  void applyRowOf2DQuadraticGridMatrix ( const RealType* Rows, const RealType* ArgCached, RealType* Dest, const int CacheSizeX, const int row ) {
  int k = 0;

  const int tx = threadIdx.x;
  const int ty = threadIdx.y;

  Dest[row] = 0;
  for ( int j = -1; j <= 1; ++j ) {
    int Y = ty + 1 + j;
    for ( int i = -1; i <= 1; ++i ) {
      int X = tx + 1 + i;
      Dest[row] += Rows[row*9+k] * ArgCached[ getGlobalIndex ( X, Y, CacheSizeX ) ];
      ++k;
    }
  }
}

/**
 * A kernel to apply cuda::GPU2DQuadraticGridMatrix to a vector while respecting that
 * the vector actually is a 2D array. Assumes to be used with 2D blocks.
 *
 * \author Berkels
 */
__global__  void apply2DQuadraticGridMatrix ( const float* Rows, const float* Arg, float* Dest, const int NumDofs, const int Width ) {
  // calculate point coordinates in the 2D grid
  unsigned int ix = blockDim.x*blockIdx.x + threadIdx.x;
  unsigned int iy = blockDim.y*blockIdx.y + threadIdx.y;

  // for bigger thread numbers there isn't a corresponding grid point
  if ( ( ix >= Width ) || ( iy >= Width ) )
    return;

  const int row = getGlobalIndex ( ix, iy, Width );

  extern __shared__ float argCached[];
  const int cacheSizeX = blockDim.x + 2;

  cache2DArray<float> ( Arg, argCached, ix, iy, Width, cacheSizeX, row );
  __syncthreads();
  applyRowOf2DQuadraticGridMatrix<float> ( Rows, argCached, Dest, cacheSizeX, row );
}

/**
 * Double version of cuda::GPU2DQuadraticGridMatrix. Unfortunately CUDA doesn't seem to allow to
 * templetize the type of shared memory (and doesn't like two versions of a kernel to use the same
 * name for their shared memory, so this code duplication can't be completely prevented.
 *
 * \author Berkels
 */
__global__  void apply2DQuadraticGridMatrix ( const double* Rows, const double* Arg, double* Dest, const int NumDofs, const int Width ) {
// Older architectures don't support double precision.
#if __CUDA_ARCH__ >= 130
  // calculate point coordinates in the 2D grid
  unsigned int ix = blockDim.x*blockIdx.x + threadIdx.x;
  unsigned int iy = blockDim.y*blockIdx.y + threadIdx.y;

  // for bigger thread numbers there isn't a corresponding grid point
  if ( ( ix >= Width ) || ( iy >= Width ) )
    return;

  const int row = getGlobalIndex ( ix, iy, Width );

  extern __shared__ double argCachedDouble[];
  const int cacheSizeX = blockDim.x + 2;

  cache2DArray<double> ( Arg, argCachedDouble, ix, iy, Width, cacheSizeX, row );
  __syncthreads();
  applyRowOf2DQuadraticGridMatrix<double> ( Rows, argCachedDouble, Dest, cacheSizeX, row );
#endif
}

/**
 * Straightforward implementation of a kernel to apply cuda::GPU2DQuadraticGridMatrix to a vector.
 * Assumes to be used with 1D blocks.
 *
 * \author Berkels
 */
__global__  void apply2DQuadraticGridMatrixSimple ( const float* Rows, const float* Arg, float* Dest, const int NumDofs, const int Width ) {
  int row = blockIdx.x * blockDim.x + threadIdx.x;

  if ( row < NumDofs ) {
    Dest[row] = 0;
    int startIndex = row - Width - 1;
    int k = 0;

    int y = row / Width;
    int x = row % Width;

    for ( int j = -1; j <= 1; ++j ) {
      int Y = y + j;
      for ( int i = -1; i <= 1; ++i ) {
        int X = x + i;
        if ( X >= 0 && X < Width && Y >= 0 && Y < Width ) {
          Dest[row] += Rows[row*9+k] * Arg[ startIndex ];
        }
        startIndex++;
        ++k;
      }
      startIndex += Width - 3;
    }
  }
}

template <typename RealType>
void GPU2DQuadraticGridMatrix<RealType>::apply ( const GPUVector<RealType> &Arg, GPUVector<RealType> &Dest ) const {
  //! \todo Make the block size variable.
#if 1
  cuda::GridBlockConfig2D gridConf( 16, 16, _gridWidth, _gridWidth );
  apply2DQuadraticGridMatrix <<<gridConf.dimGrid, gridConf.dimBlock, (gridConf.dimBlock.x+2) * (gridConf.dimBlock.y+2) * sizeof (RealType) >>>( this->getDataPointer(), Arg.getDataPointer(), Dest.getDataPointer(), _numRowAndCols, _gridWidth );
#else
  cuda::GridBlockConfig2D gridConf( 512, 1, _numRowAndCols, 1 );
  apply2DQuadraticGridMatrixSimple <<<gridConf.dimGrid, gridConf.dimBlock>>>( this->getDataPointer(), Arg.getDataPointer(), Dest.getDataPointer(), _numRowAndCols, _gridWidth );
#endif
}

template class GPU2DQuadraticGridMatrix<double>;
template class GPU2DQuadraticGridMatrix<float>;

// The CUSPARSE library needs CUDA 3.2 or later.
#if CUDA_VERSION >= 3020

#include <cusparse.h>

CUDASparseHandle::CUDASparseHandle ( )
  : _pHandle ( new cusparseHandle_t ) {
  if ( cusparseCreate ( static_cast<cusparseHandle_t*>( _pHandle ) ) != CUSPARSE_STATUS_SUCCESS ) {
    fprintf( stderr, "Error while initializing the CUSPARSE library.\n" );
    abort();
  }
  else
    fprintf( stderr, "Successfully initialized the CUSPARSE library.\n" );
}

CUDASparseHandle::~CUDASparseHandle ( ) {
  cusparseHandle_t *pHandle = static_cast<cusparseHandle_t*>( _pHandle );
  if ( cusparseDestroy ( *pHandle ) != CUSPARSE_STATUS_SUCCESS )
    fprintf( stderr, "Error while deinitializing the CUSPARSE library.\n" );

  delete pHandle;
}

template <typename RealType>
GPUCSRMatrix<RealType>::GPUCSRMatrix ( const int NumRows,
                                       const int NumCols,
                                       const int NumEntries,
                                       const RealType *Values,
                                       const int *RowStartingIndices,
                                       const int *ColumnIndices,
                                       const CUDASparseHandle &SparseHandle )
  : _values ( NumEntries ),
    _rowStartingIndices ( NumRows + 1 ),
    _columnIndices ( NumEntries ),
    _numRows ( NumRows ),
    _numCols ( NumCols ),
    _pMatDescr ( new cusparseMatDescr_t ),
    _sparseHandle ( SparseHandle ) {

  _values.copyFromHost ( Values );
  _rowStartingIndices.copyFromHost ( RowStartingIndices );
  _columnIndices.copyFromHost ( ColumnIndices );

  cusparseMatDescr_t &descr = *(static_cast<cusparseMatDescr_t *> ( _pMatDescr ));
  if ( cusparseCreateMatDescr(&descr) != CUSPARSE_STATUS_SUCCESS ) {
    fprintf( stderr, "Error calling cusparseCreateMatDescr!\n" );
    abort();
  }

  cusparseSetMatType ( descr, CUSPARSE_MATRIX_TYPE_GENERAL );
  cusparseSetMatIndexBase( descr, CUSPARSE_INDEX_BASE_ZERO );
}

template <typename RealType>
GPUCSRMatrix<RealType>::~GPUCSRMatrix ( ) {
  cusparseMatDescr_t *pMatDescr = static_cast<cusparseMatDescr_t*>( _pMatDescr );
  if ( cusparseDestroyMatDescr ( *pMatDescr ) != CUSPARSE_STATUS_SUCCESS )
    fprintf( stderr, "Error calling cusparseDestroyMatDescr!\n" );

  delete ( pMatDescr );
}

template <typename RealType>
void GPUCSRMatrix<RealType>::apply ( const GPUVector<RealType> &Arg, GPUVector<RealType> &Dest ) const {
  cusparseScsrmv( *static_cast<const cusparseHandle_t*>( _sparseHandle.getHandlePointer() ), CUSPARSE_OPERATION_NON_TRANSPOSE,
                  _numRows, _numCols, 1, *static_cast<const cusparseMatDescr_t*>( _pMatDescr ),
                  _values.getDataPointer(), _rowStartingIndices.getDataPointer(), _columnIndices.getDataPointer(),
                  Arg.getDataPointer(), 0, Dest.getDataPointer() );
}

template class GPUCSRMatrix<float>;

#endif 

}
