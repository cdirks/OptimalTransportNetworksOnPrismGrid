#include <cudaInterface.h>
#include <cudaVector.h>

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

__global__  void jacobi_kernel(float* device_A, float* device_f, float* device_u_new, float* device_u_old, int nnu) {
  /* get global thread id */
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  /* check if id corresponds to a variable */
  if (i<nnu) {
    float temp=device_f[i];
    for (int j=0;j < i  ;++j)
      temp -= device_A[i*nnu+j]*device_u_old[j];
    for (int j=i+1;j<nnu;++j)
      temp -= device_A[i*nnu+j]*device_u_old[j];
    device_u_new[i]=temp/device_A[i*nnu+i];
  }
}

void cudaJacobiSolver ( float *pMatrix, float *pRhs, float *pSolution, const int Length ) {
  /* Allocation for Memory on device */
  cuda::GPUVector<float> device_A ( Length*Length );
  cuda::GPUVector<float> device_f ( Length );
  cuda::GPUVector<float> device_u_new ( Length );
  cuda::GPUVector<float> device_u_old ( Length );

  /* Copy from host to device */
  device_A.copyFromHost ( pMatrix );
  device_f.copyFromHost ( pRhs );
  device_u_new.setZero();
  device_u_old.setZero();

  /* use one thread for every variable
  organize line of blocks with 512 threads each in a grid
  until number of variables reached */
  int block_size = 512;
  dim3 dimBlock(block_size,1,1);
  int dimenx = Length/block_size;
  if (dimenx*block_size < Length)
    dimenx++;
  dim3 dimGrid(dimenx,1,1);

  for ( int i = 0; i < 1000; ++i ) {
    /* start kernel code on device */
    jacobi_kernel<<<dimGrid, dimBlock>>>(device_A.getDataPointer(), device_f.getDataPointer(), device_u_new.getDataPointer(), device_u_old.getDataPointer(), Length);
    device_u_old = device_u_new;
  }

  /* copy back to host */
  device_u_new.copyToHost( pSolution );
}
