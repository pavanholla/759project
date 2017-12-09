# 759project

The code is set up to build four targets:

1. single_thread_out := single thread implementation of OFDM

2. openmp_out := openmp version of OFDM code

3. cuda_fft_part1 := cuda version of OFDM code (scrambler,encoder,interleaver,mapper and modulator)

4. cuda_fft := cuda version of OFDM code(fft only)


For openMP part:

The number of threads is decided through the argument in the user-input and the number of frames is also passed through user input

For the CUDA part:

It takes in one argument to fix the number of frames for the ofdm calculation

The number of threads is fixed at 128 per block and the number of blocks is fixed at FRAME_SIZE/128


