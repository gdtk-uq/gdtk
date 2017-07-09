#include <stdio.h>

__global__ void hello_kernel(char *a, int *b)
{
 a[threadIdx.x] += b[threadIdx.x];	
}

extern "C" {
       void hello(char *ad, int *bd, int blocksize)
       {
           dim3 dimBlock( blocksize, 1 );
     	   dim3 dimGrid( 1, 1 );
           hello_kernel<<<dimGrid, dimBlock>>>(ad, bd);
}
}	 