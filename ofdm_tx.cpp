#include <complex>
#include <iostream>
#include <valarray>
#include <vector>
 
#define NBPSC 2
#define NSC 64 
#define NCBPS  128


const double PI = 3.141592653589793238460;
using namespace std;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> Complex_array;

void print_array( int x[],int size)  {
  for(int ii=0 ; ii<size; ii++) cout<<x[ii];
}
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
  for(int ii=0; ii<size; ii++)  {
    data[ii]=rand()%2;
  }
}

void scramble(int data[], int size, bool init_scrambler=false )  {
  static int xor_sequence[NCBPS];
  if(init_scrambler)  {
    for(int ii=0; ii<NCBPS; ii++)  {
      xor_sequence[ii]=rand()%2;
    } 
  }
  
  for(int ii=0; ii<size; ii++)  {
    data[ii] ^= xor_sequence[ii];
  }
}
//encoder
void encode(int data[], int size, int encoded_data[])  {
  int shift_reg[] = {0,0,0,0,0,0};
  int b0; int b1;
  for(int ii=0; ii<size;ii++)  {
    encoded_data[2*ii]= data[ii] ^ shift_reg[1] ^ shift_reg[2] ^ shift_reg[4] ^ shift_reg[5];
    encoded_data[2*ii+1]= data[ii] ^ shift_reg[0] ^ shift_reg[1] ^ shift_reg[2] ^ shift_reg[5];
    for(int jj=6; jj>0; jj--)  {
      shift_reg[jj] = shift_reg [jj-1];
    }
    
    shift_reg[0] = data[ii];
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
  cout<<"Interleaving";
  for(int ii=0; ii<size; ii++) cout<<data[ii];
  for (int symbol_no =0; symbol_no< size/NCBPS; symbol_no++)  {
    for(int k=0 ; k < NCBPS; k++)  {  
      int i = (NCBPS/16) * (k % 16) + (floor(k/16)) ;
      int j = s * floor(i/s) + ((int)(i + NBPSC- floor(16 * i/NCBPS)) % s);
      interleaved_data[symbol_no*NCBPS + j] = data[symbol_no*NCBPS + k];
    }
  }
}
//modulator
void modulate ( int data[], int size, Complex modulated_data[][NSC] )  {
  cout<<"Modulating";
  for(int symbol_num=0; symbol_num< size/NCBPS; symbol_num++)  {
    for(int ii=0; ii<  NSC*NBPSC ; ii+=1)  {
      modulated_data[symbol_num][ii] =  Complex (2*data[2*ii+symbol_num*NSC*NBPSC] -1 ,2 * data[2*ii + 1 + symbol_num*NSC*NBPSC] -1 );
    }  
  }
}
//fft
 
void perform_fft_per_symbol ( Complex modulated_data[][NSC], int size )  {
  for(int ii =0 ;ii< size/NSC/NBPSC; ii++) {
    for(int jj=0; jj<NSC; jj++) cout << modulated_data[ii][jj]<<endl;
    cout << "pre Symbol complete";
    fft(modulated_data[ii]);
    for(int jj=0; jj<NSC; jj++) cout << modulated_data[ii][jj]<<endl;
    cout << "Symbol complete";
  }
}
int main()
{ 
  int frame_size = 1024;
  int frame[1024];
  int encoded_frame[1024*2];
  int interleaved_frame[1024*2];
  Complex modulated_frame[1024/NSC][NSC];
  generate_frame(frame,frame_size);
  scramble(frame,frame_size,true);
  encode(frame,frame_size,encoded_frame);
  interleave(encoded_frame,frame_size*2,interleaved_frame);
  modulate(interleaved_frame,frame_size*2,modulated_frame);
  perform_fft_per_symbol(modulated_frame, frame_size);
  return 0;
}
