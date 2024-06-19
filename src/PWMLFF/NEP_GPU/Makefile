# Compiler
CC = nvcc

# Flags
CFLAGS = -std=c++17 -O3 -arch=sm_60 -Xcompiler "-fPIC"
INC = -I./ -I./force -I./utilities
LDFLAGS = -shared
LIBS = -lcublas -lcusolver

# Source files
CUDA_SOURCES = force/nep3.cu utilities/error.cu utilities/gpu_vector.cu utilities/main_common.cu

# Object files
CUDA_OBJECTS = $(CUDA_SOURCES:.cu=.o)

# Shared library
SHARED_LIB = libnep_gpu.so

all: $(SHARED_LIB)

$(SHARED_LIB): $(CUDA_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

%.o: %.cu
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

clean:
	rm -f force/*.o utilities/*.o $(SHARED_LIB)

.PHONY: all clean

