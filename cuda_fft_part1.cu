#include<iostream>
#include<complex>
#include<stdlib.h>
#include <cuda.h>
#include <cuComplex.h>
#include <math.h>
#include<curand.h>
#include<curand_kernel.h>
#include "stopwatch.hpp"
stopwatch<std::milli, float> sw;

#define RADIUS 3

#define FRAME_SIZE 4096*8*7

#define NBPSC 2
#define NSC 64 
#define NCBPS  128
typedef std::complex<double> Complex;
using namespace std;
__global__ void init(unsigned int seed, curandState_t* states) {

     /* we have to initialize the state */
     curand_init(seed, /* the seed can be the same for each core, here we pass the time in from the CPU */
                blockIdx.x*blockDim.x+threadIdx.x, /* the sequence number should be different for each core (unless you want all
                cores to get the same sequence of numbers for some reason - use thread id! */
                0, /* the offset is how much extra we advance in the sequence for each call, can be 0 */
                &states[blockIdx.x]);
}


/* this GPU kernel takes an array of states, and an array of ints, and puts a random int into each */
__global__ void randoms(curandState_t* states, int* numbers) {
    /* curand works like rand - except that it takes a state as a parameter */
    numbers[blockIdx.x*blockDim.x+threadIdx.x] = curand(&states[blockIdx.x*blockDim.x+threadIdx.x]) % 2;
}

__global__ void scramble(int* numbers,int *scrambler_bits) {
    /* curand works like rand - except that it takes a state as a parameter */
     numbers[blockIdx.x*blockDim.x+threadIdx.x] ^=scrambler_bits[ (blockIdx.x*blockDim.x+threadIdx.x) % 128];
}


int checkResults(int startElem, int endElem, float* cudaRes, float* res)
{
    int nDiffs=0;
    const float smallVal = 0.0001f;
    for(int i=startElem; i<endElem; i++) {
        if(fabs(cudaRes[i]-res[i])>smallVal) {
            nDiffs++;
            std::cout << i << std::endl;
        }
    }
    return nDiffs;
}

void initializeWeights(float* weights, int rad)
{
    // for now hardcoded for RADIUS=3
    weights[0] = 0.50f;
    weights[1] = 0.75f;
    weights[2] = 1.25f;
    weights[3] = 2.00f;
    weights[4] = 1.25f;
    weights[5] = 0.75f;
    weights[6] = 0.50f;
}

void initializeArray(FILE* fp,float* arr, int nElements)
{
    for( int i=0; i<nElements; i++){
        	int r=fscanf(fp,"%f",&arr[i]);
		if(r == EOF){
		  rewind(fp);
		}
    }
}

void applyStencil1D_SEQ(int sIdx, int eIdx, const float *weights, float *in, float *out) {
  
  for (int i = sIdx; i < eIdx; i++) {   
    out[i] = 0;
    //loop over all elements in the stencil
    for (int j = -RADIUS; j <= RADIUS; j++) {
      out[i] += weights[j + RADIUS] * in[i + j]; 
    }
    out[i] = out[i] / (2 * RADIUS + 1);
  }
}


__global__ void interleave(int *in, int *out) {

    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    const int symbol_no = tid/NCBPS;
    const int k = tid%NCBPS;
    const int i = (NCBPS/16) * (k % 16) + ((k/16)) ;
    const int s =1; //for QPSK
    const int j = s * (i/s) + ((int)(i + NBPSC- ((16 * i)/NCBPS)) % s);
    out[symbol_no*NCBPS + j] = in[symbol_no*NCBPS + k];
    
}

__global__ void modulate(int *in, cuDoubleComplex *out) {
 
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    out[tid]=make_cuDoubleComplex(2*in[2*tid]-1,2*in[2*tid+1]-1);
}

__global__ void encode_s(int sIdx, int eIdx, int *in, int *out) {
    int i = sIdx + blockIdx.x*blockDim.x + threadIdx.x;
    volatile __shared__ int s_in[512+6];
    int offset = blockIdx.x*blockDim.x;
    s_in[threadIdx.x+3] = in[i] ;
    if(threadIdx.x<=3) {
      s_in[sIdx-threadIdx.x] = in [sIdx+offset-threadIdx.x];
      s_in[sIdx+blockDim.x + threadIdx.x] = in [ sIdx+ offset+ blockDim.x +threadIdx.x];
    }
    __syncthreads(); 
    if( i < eIdx ) {
        int result0 = 0;
        int result1 = 0;
        result0 ^= s_in[-offset+i-2];
        result0 ^= s_in[-offset+i-3];
        result0 ^= s_in[-offset+i-5];
        result0 ^= s_in[-offset+i-6];
        result0 ^= s_in[-offset+i];

        result1 ^= s_in[-offset+i-1];
        result1 ^= s_in[-offset+i-2];
        result1 ^= s_in[-offset+i-3];
        result1 ^= s_in[-offset+i-6];
        result1 ^= s_in[-offset+i];
        out[2*i] = result0;
        out[2*i+1] = result1;
    }
}


__global__ void applyStencil1D(int sIdx, int eIdx, const float *weights, float *in, float *out) {
    int i = sIdx + blockIdx.x*blockDim.x + threadIdx.x;
    if( i < eIdx ) {
        float result = 0.f;
        result += weights[0]*in[i-3];
        result += weights[1]*in[i-2];
        result += weights[2]*in[i-1];
        result += weights[3]*in[i];
        result += weights[4]*in[i+1];
        result += weights[5]*in[i+2];
        result += weights[6]*in[i+3];
        result /=7.f;
        out[i] = result;
    }
}

int main(int argc, char* argv[]) {

    int num_frames=10;

    if(argc>1) num_frames = atoi ( argv[1]);
    float time = 0.f;

    float generate_time= 0, scramble_time=0, encode_time=0, encode_parallel_time=0, interleave_time=0, modulate_time = 0, fft_time=0 ; 
    int frame_size = FRAME_SIZE;
    int *frame;
    int *scrambler_bits;
    int *encoded_frame;
    int *interleaved_frame;
    cuDoubleComplex *modulated_frame;

    cudaEvent_t startevent, stopevent;
    cudaEventCreate(&startevent);
    cudaEventCreate(&stopevent);
    cudaEventRecord(startevent,0);  

    curandState_t* states;
    /* allocate space on the GPU for the random states */
    cudaMalloc((void**) &states, FRAME_SIZE * sizeof(curandState_t));
    /* invoke the GPU to initialize all of the random states */
    init<<<FRAME_SIZE/128, 128>>>(1000, states);

    cudaMalloc((void**) &scrambler_bits, 128 * sizeof(int));
    randoms<<<1, 128>>>(states, scrambler_bits);
    /* allocate an array of unsigned ints on the CPU and GPU */
    cudaMalloc((void**) &frame, FRAME_SIZE * sizeof(int));
    cudaMalloc((void**) &encoded_frame, 2*FRAME_SIZE * sizeof(int));
    cudaMalloc((void**) &interleaved_frame, 2*FRAME_SIZE * sizeof(int));
    cudaMalloc((void**) &modulated_frame, FRAME_SIZE * sizeof(cuDoubleComplex));
    /* invoke the kernel to get some random numbers */
    for(int ii=0; ii<num_frames;ii++) {
      sw.start();
      randoms<<<FRAME_SIZE/128, 128>>>(states, frame);
      cudaDeviceSynchronize();  sw.stop(); generate_time += sw.count(); sw.start();
      scramble<<<FRAME_SIZE/128, 128>>>( frame,scrambler_bits);
      cudaDeviceSynchronize();  sw.stop(); scramble_time += sw.count(); sw.start();
      encode_s<<<FRAME_SIZE/512, 512>>>(6,FRAME_SIZE, frame,encoded_frame);
      cudaDeviceSynchronize();  sw.stop(); encode_time += sw.count(); sw.start();
      interleave<<<FRAME_SIZE/512, 512>>> (encoded_frame,interleaved_frame);
      cudaDeviceSynchronize();  sw.stop(); interleave_time += sw.count(); sw.start();
      modulate<<<FRAME_SIZE/512, 512>>> (interleaved_frame,modulated_frame);
      cudaDeviceSynchronize();  sw.stop(); modulate_time += sw.count();
    }
      cout<< " generate_time " << generate_time << endl;
      cout<< " scramble_time " << scramble_time << endl;
      cout<< " encode_time " << encode_time << endl;
      cout<< " interleave_time " << interleave_time << endl;
      cout<< " modulate_time " << modulate_time << endl;
    cudaEventRecord(stopevent,0);  //ending timing for inclusive
     cudaEventSynchronize(stopevent);   
     cudaEventElapsedTime(&time, startevent, stopevent);
    cout << "Total time minus fft" << time << endl;

    
    return 0;
}
