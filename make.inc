# Platform
ifeq ($(BUILD_TYPE),)
#BUILD_TYPE = ethz-gnu64
#BUILD_TYPE = icsmaster-gnu64
	BUILD_TYPE = kaust
endif

# Default values which might be overridden
MAKE     = make
AR       = ar
ARFLAGS  = ruv
AR_FLAGS = ruv
RANLIB   = ranlib

ifeq ($(BUILD_TYPE),kaust)
        FLAGS         = -DAdd_ -DSOLVE_CRAY -D__MINMAX -Wno-write-strings
        MKLHOME       = $(MKLROOT)/lib/intel64
        BLAHOME       = $(MKLHOME)
        LAPHOME       = $(MKLHOME)
        SCAHOME       = $(MKLHOME)
        CUDAHOME      = /usr/local/cuda
        MKL           = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
        LOADER        = mpic++ -Wl,-rpath -Wl,$(INTEL64)/compiler/lib/intel64
        MPICXX        = mpic++ -Wl,-rpath -Wl,$(INTEL64)/compiler/lib/intel64
        NVCC          = nvcc
        LOADFLAGS     = -O3 -ffast-math -funroll-loops -L$(GNU64)/lib64 -lgfortran -L$(INTEL64)/compiler/lib/intel64 -L$(MPIHOME)/lib/ #-lmpifort
        NVCCFLAGS     = -arch=compute_70 -code=sm_70
        CUDA_SDK      = 3
        COMMONFLAGS   = -DCUDA_SDK=$(CUDA_SDK)
        F90_LIBS      = -L$(GNU64)/lib64 -lgfortran -Wl,-rpath,L$(GNU64)/lib64
        LAPACK        = -L$(LAPHOME) -lmkl_gf_lp64 -lmkl_core -Wl,-rpath -Wl,$(MKLHOME)
        #SCALAPACK     = -L/lib/x86_64-linux-gnu -lscalapack-openmpi #-L$(SCAHOME) -lmkl_scalapack_lp64 -Wl,-rpath -Wl,$(MKLHOME)
        LINKS         = -Wl,-rpath -Wl,$(GNU64)/lib64
        LIBSCI        =
        BLAS          =
        OPENMP        = -fopenmp
        export LD_LIBRARY_PATH=$(GNU64)/lib64:$(INTEL64)/lib
endif

#INCCUD        = -I$(CUDAHOME)/include
#CUDA          = -L$(CUDAHOME)/lib64 -lcuda -L$(CUDAHOME)/lib64 -lcudart -lcublas -lcufft -lcusparse -Wl,-rpath,$(CUDAHOME)/lib64
INCMAG        = -I/home/x_pollakgr/RGF/applications/magma-2.5.4/include
MAGMA         = -L/home/x_pollakgr/applications/magma-2.5.4/lib -lmagma

INCCUD	      = -I$(CUDAHOME)/include
# -lcufft -lcuda (Cuda Driver API) not need we use (cuda runtime library)
CUDA          = -L$(CUDAHOME)/lib64 -lcudart -lcublas -lcusparse -lnvToolsExt -Wl,-rpath,$(CUDAHOME)/lib64
#INCMAG        ?= -I./magma-1.6.2p4/include/
#MAGMA         ?= -L./magma-1.6.2p4/lib/ -lmagma

# Common definitions
CD            = cd
ECHO          = echo
LN            = ln
LNFLAGS       = -s
RM            = rm
RMFLAGS       = -f
MV            = mv
MVFLAGS       = -f
SHELL         = sh
LIB           = -lm
LEX           = flex
YACC          = bison
