#include <complex>
#include <iostream>
#include <valarray>
 
const double PI = 3.141592653589793238460;
using namespace std;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
 
//// Cooley–Tukey FFT (in-place, divide-and-conquer)
//// Higher memory requirements and redundancy although more intuitive
//void fft(CArray& x)
//{
//    const size_t N = x.size();
//    if (N <= 1) return;
// 
//    // divide
//    CArray even = x[std::slice(0, N/2, 2)];
//    CArray  odd = x[std::slice(1, N/2, 2)];
// 
//    // conquer
//    fft(even);
//    fft(odd);
// 
//    // combine
//    for (size_t k = 0; k < N/2; ++k)
//    {
//        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
//        x[k    ] = even[k] + t;
//        x[k+N/2] = even[k] - t;
//    }
//}
 
// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// Better optimized but less intuitive
// !!! Warning : in some cases this code make result different from not optimased version above (need to fix bug)
// The bug is now fixed @2017/05/30 
void fft(CArray &x)
{
	// DFT
	unsigned int N = x.size(), k = N, n;
	double thetaT = 3.14159265358979323846264338328L / N;
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
	//// Normalize (This section make it not working correctly)
	//Complex f = 1.0 / sqrt(N);
	//for (unsigned int i = 0; i < N; i++)
	//	x[i] *= f;
}
 
// inverse fft (in-place)
void ifft(CArray& x)
{
    // conjugate the complex numbers
    x = x.apply(std::conj);
 
    // forward fft
    fft( x );
 
    // conjugate the complex numbers again
    x = x.apply(std::conj);
 
    // scale the numbers
    x /= x.size();
}
void generate_frame(int data[],int size)  {
  for(int ii=0; ii<size; ii++)  {
    data[ii]=rand()%2;
  }
}

void scramble(int data[], int size, bool init_scrambler=false )  {
  static int xor_sequence[128];
  if(init_scrambler)  {
    for(int ii=0; ii<128; ii++)  {
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
  int NCBPS = 96;
  int NBPSC = 2; 
  //int s=((NBPSC/2)>1) ?  (NBPSC/2):1 ;
  int s = 1;
  for(int k=0 ; k < size; k++)  {  
    int i = (NCBPS/16) * (k % 16) + (floor(k/16)) ;
    int j = s * floor(i/s) + ((int)(i + NBPSC- floor(16 * i/NCBPS)) % s);
    interleaved_data[j] = data[k];
  }
}
//modulator
void modulate ( int data[], int size, Complex modulated_data[] )  {
  int NBPSC = 2; 
  for(int ii=0; ii< size/NBPSC ; ii++)  {
    
  modulated_data[ii] =  Complex (data[2*ii],data[2*ii + 1] );
  }
}
//fft
 
void perform_fft_per_symbol ( Complex modulated_data[], int size )  {
  for(int ii =0 ;ii< size/64; ii++) {
    CArray data(modulated_data[ii*64], 64);
    std::cout << modulated_data[ii] << std::endl;
    fft(data);
    for (int i = 0; i < 64; ++i)  {
        std::cout << data[i] << std::endl;
     }
  }
}
int main()
{ 
  int frame_size = 1024;
  int frame[1024];
  int encoded_frame[1024*2];
  int interleaved_frame[1024*2];
  Complex modulated_frame[1024];
  generate_frame(frame,frame_size);
  encode(frame,frame_size,encoded_frame);
  interleave(encoded_frame,frame_size*2,interleaved_frame);
  modulate(interleaved_frame,frame_size*2,modulated_frame);
  perform_fft_per_symbol(modulated_frame, frame_size);
//   const Complex test[] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
//
//    
//    CArray data(test, 8);
// 
//    // forward fft
//    fft(data);
// 
//    std::cout << "fft" << std::endl;
//    for (int i = 0; i < 8; ++i)
//    {
//        std::cout << data[i] << std::endl;
//    }
// 
//    // inverse fft
//    ifft(data);
// 
//    std::cout << std::endl << "ifft" << std::endl;
//    for (int i = 0; i < 8; ++i)
//    {
//        std::cout << data[i] << std::endl;
//    }
    return 0;
}
