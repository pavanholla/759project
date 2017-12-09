# Warnings
WFLAGS	:= -Wall -Wextra -Wsign-conversion -Wsign-compare

# Optimization and architecture
OPT		:= -O3
ARCH   	:= -march=native

# Language standard
CCSTD	:= -std=c99
CXXSTD	:= -std=c++11

EXEC := single_thread_out openmp_out #cuda_fft 

# Linker options
LDOPT 	:= $(OPT)
LDFLAGS := 
BIN = "/usr/local/gcc/6.4.0/bin/gcc"
.DEFAULT_GOAL := all

.PHONY: debug
debug : OPT  := -O0 -g -fno-omit-frame-pointer -fsanitize=address
debug : LDFLAGS := -fsanitize=address
debug : ARCH :=
debug : $(EXEC)

all : $(EXEC)


#module load cuda;nvcc -o problem2 $(OPT) vector_reduction.cu vector_reduction_gold.cpp -ccbin $(BIN)
#module load cuda;nvcc -o problem3 $(OPT) problem3.cu -ccbin $(BIN)

single_thread_out: ofdm_tx.cpp
	gcc -g $(CXXSTD) -lm  -fopt-info -lstdc++ ofdm_tx.cpp -o single_thread_out

openmp_out: ofdm_openmp.cpp
	gcc -g $(CXXSTD) -lm  -fopt-info -lstdc++ -fopenmp ofdm_openmp.cpp -o openmp_out -DOMP -DOMP_ENCODE_PARALLEL

cuda_fft: ofdm_tx.cu
	module load cuda;nvcc -o cuda_fft -O1 ofdm_tx.cu --ptxas-options=-v --use_fast_math -lcufft -ccbin $(BIN)

problem1a: problem1a.cu
	module load cuda;nvcc -o problem1a -O1 problem1a.cu --ptxas-options=-v --use_fast_math -lcufft -ccbin $(BIN)
# TODO: add targets for building executables

.PHONY: clean
clean:
	rm -f single_thread_out openmp_out cuda_fft
