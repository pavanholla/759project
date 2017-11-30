# Warnings
WFLAGS	:= -Wall -Wextra -Wsign-conversion -Wsign-compare

# Optimization and architecture
OPT		:= -O3
ARCH   	:= -march=native

# Language standard
CCSTD	:= -std=c99
CXXSTD	:= -std=c++11

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

all : problem1



problem1:  
	gcc -g $(CXXSTD) -lm  -fopt-info -lstdc++ ofdm_tx.cpp -o single_thread_out
	gcc -g $(CXXSTD) -lm  -fopt-info -lstdc++ -fopenmp ofdm_tx.cpp -o openmp_out -DOMP -DOMP_ENCODE_PARALLEL
	#module load cuda;nvcc -o problem1 $(OPT) problem1.cu -ccbin $(BIN)
	#module load cuda;nvcc -o problem2 $(OPT) vector_reduction.cu vector_reduction_gold.cpp -ccbin $(BIN)
	#module load cuda;nvcc -o problem3 $(OPT) problem3.cu -ccbin $(BIN)



# TODO: add targets for building executables

.PHONY: clean
clean:
	rm -f problem1
