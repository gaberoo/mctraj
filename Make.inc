CC = gcc
CPP = g++
FC = gfortran
LD = g++
CFLAGS = -Imodels -I. -I.. -isystem ./include
#CFLAGS += -DDEBUG
#CFLAGS += -fopenmp 
#CFLAGS += -std=c++0x -DMKLRNG 
#CFLAGS += -DUNROLL
CFLAGS += -g -m64 -Wall -O3
FFLAGS += -g -m64 -Wall -O3
#LDFLAGS += -lgfortran

CPPFLAGS = $(CFLAGS)

BLAS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
#BLAS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#BLAS = -llapack -lgoto2 -lpthread
#BLAS = -framework Accelerate

LDFLAGS += $(BLAS)

MPICC = $(CC) `openmpicc -showme:compile`
MPICPP = $(CPP) `openmpic++ -showme:compile`
MPICFLAGS = -DUSE_MPI $(CFLAGS)
MPIFFLAGS += -DUSE_MPI $(CPPFLAGS)
MPILD = $(CPP)
MPILDFLAGS = $(LDFLAGS) `openmpic++ -showme:link`

