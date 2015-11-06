CC = gcc
CPP = g++
FC = gfortran
LD = g++

#CC = clang
#CPP = clang++
#LD = clang++

CFLAGS += -std=c++11
CFLAGS += -Imodels -I. -I.. -isystem ./include
#CFLAGS += -DDEBUG
#CFLAGS += -fopenmp 
#CFLAGS += -DMKLRNG 
#CFLAGS += -DUNROLL
CFLAGS += -g -m64 -Wall -O3
FFLAGS += -g -m64 -Wall -O3
#LDFLAGS += -lgfortran

CPPFLAGS = $(CFLAGS)

#BLAS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread
BLAS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#BLAS = -llapack -lgoto2 -lpthread
#BLAS = -framework Accelerate

LDFLAGS += $(BLAS)

MPICC = $(CC) `openmpicc -showme:compile`
MPICPP = $(CPP) `openmpic++ -showme:compile`
MPICFLAGS = -DUSE_MPI $(CFLAGS)
MPIFFLAGS += -DUSE_MPI $(CPPFLAGS)
MPILD = $(CPP)
MPILDFLAGS = $(LDFLAGS) `openmpic++ -showme:link`

