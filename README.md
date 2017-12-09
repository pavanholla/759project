# 759project

The code is set up to build four targets:

1. single_thread_out := single thread implementation of OFDM

2. openmp_out := openmp version of OFDM code

3. cuda_fft_part1 := cuda version of OFDM code (scrambler,encoder,interleaver,mapper and modulator)

4. cuda_fft := cuda version of OFDM code(fft only)


For openMP part:

The number of threads is decided through the argument in the user-input and the number of frames is 512

For the CUDA part:

The number of frames for the ofdm calculation is 512

The number of threads is fixed at 128 per block and the number of blocks is fixed at FRAME_SIZE/128

RUN FILES:
To run the code on Euler, we need ".sh" files

Following four ".sh" files are generated to run each of the four codes:
run_single_thread.sh
run_openmp.sh
run_cuda_fft_part1.sh
run_cuda_fft.sh

The ".out" files contain the results for each run.
