#include <complex>
#include <iostream>
#include <valarray>
#include <vector>
#include <bitset>
#include "stopwatch.hpp"
#include <cuda.h>
//#include <omp.h>
#include <cufft.h>
 
#define NBPSC 2
#define NSC 64 
#define NCBPS  128


const double PI = 3.141592653589793238460;
using namespace std;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> Complex_array;
int NUM_THREADS = 2;
#define FRAME_SIZE 4096*8*7
stopwatch<std::milli, float> sw;

void print_array( int x[],int size)  {
  for(int ii=0 ; ii<size; ii++) cout<<x[ii];
}
//substitute cuda here ( optionally vectorized fft)
void fft(Complex x[])
{
	// DFT
	unsigned int N = NSC, k = N, n;
	double thetaT = 3.141592653589793238462/ N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				Complex t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
}
 
void generate_frame(int data[],int size)  {
  //parallel for
  //#pragma omp parallel for num_threads(NUM_THREADS) schedule(static) 
  for(int ii=0; ii<size; ii++)  {
    //one block per 128 numbers
    //num_blocks=size/128
    #ifdef OMP
    unsigned int myseed = omp_get_thread_num();
    data[ii]=rand_r(&myseed)%2;
    #else
    data[ii]=rand()%2;
    #endif
  }
}

void scramble(int data[], int size, bool init_scrambler=false )  {
  static int xor_sequence[NCBPS];
  //one block per 128 numbers  
  if(init_scrambler)  {
    for(int ii=0; ii<NCBPS; ii++)  {
      xor_sequence[ii]=rand()%2;//hardcode for cuda, input ot kernel
    } 
  } else  {
  
  //parallel for
  //#pragma omp parallel for num_threads(NUM_THREADS) schedule(static) 
    for(int ii=0; ii<size; ii++)  {
      data[ii] ^= xor_sequence[ii%NCBPS];//merge with generate_frame kernel
    } 
  }
}
//encoder
void encode(int data[], int size, int encoded_data[], int shift_reg[], int start=0)  {
  //convolution
  //can be sped up by performing encoding on a symbol by symbol basis
  //int b0; int b1;
  for(int ii=start; ii<size;ii++)  {
    //convolution
    encoded_data[2*ii]= data[ii] ^ shift_reg[1] ^ shift_reg[2] ^ shift_reg[4] ^ shift_reg[5];
    encoded_data[2*ii+1]= data[ii] ^ shift_reg[0] ^ shift_reg[1] ^ shift_reg[2] ^ shift_reg[5];
    for(int jj=6; jj>0; jj--)  {
      shift_reg[jj] = shift_reg [jj-1];
    }
    
    shift_reg[0] = data[ii];
  }

  
}

void encode_parallel(int data[], int size, int encoded_data_t0[])  {
    
//   encode(data+ENCODE_CHUNK, ENCODE_CHUNK , encoded_data_t0);
//  int encoded_data_t0[ENCODE_CHUNK*2];
  //#pragma omp parallel for num_threads(NUM_THREADS) schedule(static) default(shared) 
  for(int ii=0;ii<NUM_THREADS;ii++)  {
    int shift_reg_c[]={0,0,0,0,0,0};
    int *shift_reg=shift_reg_c;
    if(ii>0)
      shift_reg = data + ii*size/NUM_THREADS;
    encode(data, size/NUM_THREADS , encoded_data_t0,shift_reg, ii*size/NUM_THREADS);
  }

}
  
//interleaver

//interleaving process is a two permutation block interleaving. The first permutation is defined by the following:
//
//i = (NCBPS/16) · (k mod 16) + (floor(k/16) i=0,1,....,NCBPS-1
//
//where i is the index after the permutation of the k'th input bit in the block of NCBPS bits, and NCBPS is defined by the data rate.
//
//The second permutation is defined by the following:
//
//j = s * floor(i/s) + (i + NBPSC- floor(16 * i/NCBPS)) mod s
//
//s=max(NBPSC/2,1)
//
//where j is the index after the permutation of the i'th bit from the first permutation. The values for NCBPS and NBPSC for the different data rates are:
//
//Data Rate (Mbits/s) NCBPS NBPSC
//6 48  1
//9 48  1
//12  96  2
//18  96  2
//24  192 4
//36  192 4
//48  288 6
//54  288 6
void interleave( int data[], int size, int interleaved_data[])  {
  //int s=((NBPSC/2)>1) ?  (NBPSC/2):1 
  int s = 1;
  #ifdef debug_print
  cout<<"Interleaving";
  for(int ii=0; ii<size; ii++) cout<<data[ii];
  #endif
  //parallel for
  //#pragma omp parallel for num_threads(NUM_THREADS) schedule(static) 
  for (int symbol_no =0; symbol_no< size/NCBPS; symbol_no++)  {
    //one block per symbol. 128 bits get arranged in rows of 16 each, and are read out column wise, perhaps gather
    for(int k=0 ; k < NCBPS; k++)  {  
      int i = (NCBPS/16) * (k % 16) + (floor(k/16)) ;
      int j = s * floor(i/s) + ((int)(i + NBPSC- floor(16 * i/NCBPS)) % s);
      interleaved_data[symbol_no*NCBPS + j] = data[symbol_no*NCBPS + k];
    }
  }
}
//modulator
void modulate ( int data[], int size, Complex modulated_data[][NSC] )  {
  #ifdef debug_print
  cout<<"Modulating";
  #endif
  //parallel for
  //#pragma omp parallel for num_threads(NUM_THREADS) schedule(static) 
  for(int symbol_num=0; symbol_num< size/NCBPS; symbol_num++)  {
    for(int ii=0; ii<  NSC*NBPSC ; ii+=1)  {
      modulated_data[symbol_num][ii] =  Complex (2*data[2*ii+symbol_num*NSC*NBPSC] -1 ,2 * data[2*ii + 1 + symbol_num*NSC*NBPSC] -1 );
    }  
  }
}
//fft
 
void perform_fft_per_symbol ( Complex modulated_data[][NSC], int size )  {
  //parallel for
  
  //#pragma omp parallel for num_threads(NUM_THREADS) schedule(static) 
  for(int ii =0 ;ii< size/NSC/NBPSC; ii++) { //for each symbol
    #ifdef debug_print
    for(int jj=0; jj<NSC; jj++) cout << modulated_data[ii][jj]<<endl;
    cout << "pre Symbol complete";
    #endif
    // cuda implementation

     //cuda implementation ended
    //cufftExecC2C(plan, (cufftComplex *)d_modulatedFrame+ii*64,(cufftComplex *)d_modulatedFrame+ii*64,CUFFT_FORWARD);

    
    
    
    
    //fft(modulated_data[ii]);
    #ifdef debug_print
    for(int jj=0; jj<NSC; jj++) cout << modulated_data[ii][jj]<<endl;
    cout << "Symbol complete";
    #endif
  }

}
int main(int argc, char* argv[])
{ 

    int frame_size = FRAME_SIZE;
    int frame[FRAME_SIZE];
    int encoded_frame[FRAME_SIZE];
    int interleaved_frame[FRAME_SIZE*2];
    Complex modulated_frame[FRAME_SIZE/NSC][NSC];
    
    int shift_reg[] = {0,0,0,0,0,0};
    float generate_time= 0, scramble_time=0, encode_time=0, encode_parallel_time=0, interleave_time=0, modulate_time = 0, fft_time=0 ; 
    int num_frames = 10;
    if(argc>1) num_frames = atoi ( argv[1]);
    if(argc>2) NUM_THREADS = atoi ( argv[2]);
  //parallel for
    scramble(frame,frame_size,true); //initialize scrambler
    for(int ii =0 ; ii <num_frames; ii++)  {
      sw.start();
      generate_frame(frame,frame_size);
      sw.stop(); generate_time += sw.count(); sw.start();
      scramble(frame,frame_size,false); 
      sw.stop(); scramble_time += sw.count(); sw.start();
      #ifdef OMP_ENCODE_PARALLEL
      encode_parallel(frame,frame_size,encoded_frame);
      sw.stop(); encode_parallel_time += sw.count(); sw.start();
      #endif
      encode(frame,frame_size,encoded_frame,shift_reg);
      sw.stop(); encode_time += sw.count(); sw.start();
      interleave(encoded_frame,frame_size*2,interleaved_frame);
      sw.stop(); interleave_time += sw.count(); sw.start();
      modulate(interleaved_frame,frame_size*2,modulated_frame);
      sw.stop(); modulate_time += sw.count();  
     
        Complex *d_modulatedFrame;
        cudaMalloc((void**)&d_modulatedFrame,sizeof(Complex)*frame_size);
        cudaMemcpy(d_modulatedFrame,modulated_frame,sizeof(Complex)*frame_size,cudaMemcpyHostToDevice);
        cufftHandle plan;
        cufftPlan1d(&plan,NSC,CUFFT_C2C,frame_size/64);


        sw.start();


        cufftExecC2C(plan, (cufftComplex *)d_modulatedFrame,(cufftComplex *)d_modulatedFrame,CUFFT_FORWARD);
      //perform_fft_per_symbol(modulated_frame, frame_size);
        sw.stop(); fft_time += sw.count();

        cudaMemcpy(modulated_frame,d_modulatedFrame,sizeof(Complex)*frame_size,cudaMemcpyDeviceToHost);
        cudaFree(d_modulatedFrame);
        cufftDestroy(plan);
    }
    //cout<< " generate_time " << generate_time << endl;
    //cout<< " scramble_time " << scramble_time << endl;
    //cout<< " encode_time " << encode_time << endl;
    //cout<< " encode_parallel_time " << encode_parallel_time << endl;
    //cout<< " interleave_time " << interleave_time << endl;
    //cout<< " modulate_time " << modulate_time << endl;
    cout<< " fft_time " << fft_time << endl;
  return 0;
}

//dynamic scheduling is disastrous
//rand_r
