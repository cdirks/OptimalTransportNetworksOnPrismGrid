#include <cudaInterface.h>
#include <cudaVector.h>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//divergenz
__device__ float div ( float* p1, float* p2, int i, int j, int nx, int ny ) {
  float div = 0;

  if ( i == 0 )
    div += p1[ cuda::getGlobalIndex ( i, j, nx ) ];
  if( i == nx-1 )
    div += -p1[ cuda::getGlobalIndex ( i-1, j, nx ) ];
  if ( 0 < i && i < nx-1 )
    div += p1[ cuda::getGlobalIndex ( i, j, nx ) ] - p1[ cuda::getGlobalIndex ( i-1, j, nx ) ];

  if ( j == 0 )
    div += p2[ cuda::getGlobalIndex ( i, j, nx ) ];
  if ( j == ny-1 )
    div += -p2[ cuda::getGlobalIndex ( i, j-1, nx ) ];
  if ( 0 < j && j < ny-1 )
    div += p2[ cuda::getGlobalIndex ( i, j, nx ) ]-p2[ cuda::getGlobalIndex ( i, j-1, nx ) ];

  return div;
}

__device__ float dXFD2D ( float* image, int x, int y, int nx, int ny ) {
  return image[ cuda::getGlobalIndex ( x + 1, y, nx ) ] - image[ cuda::getGlobalIndex ( x, y, nx ) ];
}

__device__ float dYFD2D ( float* image, int x, int y, int nx, int ny ) {
  return image[ cuda::getGlobalIndex ( x, y + 1, nx ) ] - image[ cuda::getGlobalIndex ( x, y, nx ) ];
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class TVChambolleHelper {
public:
  static __device__ float func ( float* p1, float* p2, float *image, float *, int x, int y, int nx, int ny ) {
    return image[ cuda::getGlobalIndex ( x, y, nx ) ] - div ( p1, p2, x, y, nx, ny );
  };
};

class MSSegHelper {
public:
  static __device__ float func ( float* p1, float* p2, float *twiceIndicator2OverLambda, float *twiceIndicator1Plus2, int x, int y, int nx, int ny ) {
    const int globalIndex = cuda::getGlobalIndex ( x, y, nx );
    return ( twiceIndicator2OverLambda[globalIndex] - div ( p1, p2, x, y, nx, ny ) ) / twiceIndicator1Plus2[globalIndex];
  };
};

template <typename T>
__global__ void update_duals ( float* p1, float* p2, float* p1New, float* p2New,
                               float* in1, float* in2, float timestep,
                               int nx, int ny) {
  // calculate point in colormap discretization
  unsigned int ix = blockDim.x*blockIdx.x + threadIdx.x;
  unsigned int iy = blockDim.y*blockIdx.y + threadIdx.y;

  // for bigger thread numbers there isn't a corresponding point in colormap discretization
  if ( ( ix >= nx ) || ( iy >= ny ) )
    return;

  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  const int globalIdx = cuda::getGlobalIndex ( ix, iy, nx );

  // calculate ( image - div p ) and store it in shared memory, so that all threads can access it.
  // The amound of memory necessary for this is specified in the execution configuration of update_duals.
  extern __shared__ float imageMinusDivP[];
  imageMinusDivP[ cuda::getGlobalIndex ( tx, ty, blockDim.x + 1 ) ] =  T::func ( p1, p2, in1, in2, ix, iy, nx, ny );
  if ( tx == blockDim.x-1 )
    imageMinusDivP[ cuda::getGlobalIndex ( tx + 1, ty, blockDim.x + 1 ) ] = ( ix < nx-1 ) ? T::func ( p1, p2, in1, in2, ix + 1, iy, nx, ny ) : 0;
  if ( ty == blockDim.y-1 )
    imageMinusDivP[ cuda::getGlobalIndex ( tx, ty + 1, blockDim.x + 1 ) ] = ( iy < ny-1 ) ? T::func ( p1, p2, in1, in2, ix, iy + 1, nx, ny ) : 0;
  __syncthreads();

  //-\nabla colormap + \nabla ind *div p + ind*\nabla div p
  const float buffer1 = ( ix < nx-1 ) ? - dXFD2D ( imageMinusDivP, tx, ty, blockDim.x + 1, blockDim.y + 1 ) : 0;
  const float buffer2 = ( iy < ny-1 ) ? - dYFD2D ( imageMinusDivP, tx, ty, blockDim.x + 1, blockDim.y + 1 ) : 0;

  const float norm = sqrt(buffer1*buffer1+buffer2*buffer2);

  //update p1, p2
  p1New[globalIdx] = ( p1[globalIdx] + timestep*buffer1 ) / ( 1 + timestep*norm);
  p2New[globalIdx] = ( p2[globalIdx] + timestep*buffer2 ) / ( 1 + timestep*norm);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename T>
__global__ void constructImageFromDuals ( float* in1, float* in2, float* p1, float* p2,
                                          float lambda, int nx, int ny) {
  // calculate point in colormap discretization
  unsigned int ix = blockDim.x*blockIdx.x + threadIdx.x;
  unsigned int iy = blockDim.y*blockIdx.y + threadIdx.y;

  // for bigger thread numbers there isn't a corresponding point in colormap discretization
  if ( ( ix >= nx ) || ( iy >= ny ) )
    return;

  in1[ cuda::getGlobalIndex ( ix, iy, nx ) ] = lambda * T::func( p1, p2, in1, in2, ix, iy, nx, ny );
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename T>
void ChambolleTVAlgo ( cuda::GPUVector<float> &DevIn1, cuda::GPUVector<float> &DevIn2, int nx, int ny, float lambda, float tol, float timeStep, int maxIter ) {
  // Size of a thread block, should be optimized for the GPU the code is run on.
  // The size of 16x16 seems to work well for a GeForce 8800 GT.
  const int blockSizeX = 16;
  const int blockSizeY = 16;

  ////////////////////////////////////////////////////////////////////////////////
  // compute number of Blocks & Grid
  cuda::GridBlockConfig2D gridConf(  blockSizeX, blockSizeY, nx, ny );

  cuda::GPUMultiVector<float> dev_p ( 2, nx*ny ); 
  dev_p.setZero();
  cuda::GPUMultiVector<float> dev_pNew ( 2, nx*ny );
  dev_pNew.setZero();

  ////////////////////////////////////////////////////////////////////////////////

  cublasStatus stat = cublasInit();
  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf ("Failed to initialize the CUBLAS library\n" );
  }

  cuda::GPUMultiVector<float> errorMVec ( 2, nx*ny );
  errorMVec.setZero();

  float error = tol + 1;
  //calculate dual variables p1,p2 through fix point iteration scheme
  for ( int i = 1; i < maxIter; ++i ) {
    //calculate p^n+1= (p^n + timestep*buffer) / (1+timestep*norm)
    // with     buffer =  \nabla (div p^n - image / lambda)
    // and      norm = ||buffer||_2
    update_duals<T><<<gridConf.dimGrid, gridConf.dimBlock, (gridConf.dimBlock.x+1) * (gridConf.dimBlock.y+1) * sizeof (float) >>>(dev_p[0].getDataPointer(), dev_p[1].getDataPointer(), dev_pNew[0].getDataPointer(), dev_pNew[1].getDataPointer(), DevIn1.getDataPointer(), DevIn2.getDataPointer(), timeStep, nx, ny);
    CUT_CHECK_ERROR("ERROR: dual Kernel Failed !");
    dev_p = dev_pNew;
    if ( i % 100 == 0 ) {
      //save old variables
      errorMVec = dev_p;
    }
    if ( i % 100 == 1 ) {
      //difference of old and new
      errorMVec -= dev_p;
      error = errorMVec[0].norm() + errorMVec[1].norm();
    }

    if ( abs(error) < tol ) {
      break;
    }
  }

  //reconstruct colormap from duals
  constructImageFromDuals<T><<<gridConf.dimGrid, gridConf.dimBlock>>>(DevIn1.getDataPointer(), DevIn2.getDataPointer(), dev_p[0].getDataPointer(), dev_p[1].getDataPointer(), lambda, nx, ny);
  CUT_CHECK_ERROR("ERROR: construct Kernel Failed !");
}

void tv_smooth_gpu(float* image, int nx, int ny, float lambda, float tol, float timeStep, int maxIter) {

  ////////////////////////////////////////////////////////////////////////////////
  //allocation & init

  cuda::GPUVector<float> dev_image ( nx*ny );
  dev_image.copyFromHost ( image );
  dev_image *= 1.f / lambda;

  ChambolleTVAlgo<TVChambolleHelper>( dev_image, dev_image, nx, ny, lambda, tol, timeStep, maxIter);

  ////////////////////////////////////////////////////////////////////////////////
  // copy back to host
  dev_image.copyToHost ( image );
}

/**
 * Result is stored in indicator1.
 *
 * \author Berkels
 */
void MSSegGPU(float* indicator1, float* indicator2, int nx, int ny, float lambda, float tol, float timeStep, int maxIter) {
  ////////////////////////////////////////////////////////////////////////////////
  //allocation & init
  const float lambdaHDependent = lambda * ( ( nx > ny ) ? nx - 1 : ny - 1 );

  cuda::GPUVector<float> twiceIndicator1Plus2 ( nx*ny );
  twiceIndicator1Plus2.copyFromHost ( indicator1 );
  cuda::GPUVector<float> twiceIndicator2OverLambda ( nx*ny );
  twiceIndicator2OverLambda.copyFromHost ( indicator2 );
  twiceIndicator1Plus2 += twiceIndicator2OverLambda;
  twiceIndicator1Plus2 *= 2.f;
  twiceIndicator2OverLambda *= 2.f / lambdaHDependent;

  ChambolleTVAlgo<MSSegHelper>( twiceIndicator2OverLambda, twiceIndicator1Plus2, nx, ny, lambdaHDependent, tol, timeStep, maxIter);

  ////////////////////////////////////////////////////////////////////////////////
  // copy back to host
  twiceIndicator2OverLambda.copyToHost ( indicator1 );
}
