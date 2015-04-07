CC = gcc
CPP = g++
FC = gfortran
LD = g++
CFLAGS += -std=c++0x #-DMKLRNG #-DDEBUG
#CFLAGS += -DUNROLL
CFLAGS += -g -m64 -Wall -O3
FFLAGS += -g -m64 -Wall -O3
LDFLAGS += -lgfortran

CPPFLAGS = $(CFLAGS)

BLAS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
#BLAS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#BLAS = -llapack -lgoto2 -lpthread
#BLAS = -framework Accelerate

LDFLAGS += $(BLAS)

