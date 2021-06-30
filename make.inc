# Platform
ifeq ($(BUILD_TYPE),)
# Set build type to kaust if non specified
# BUILD_TYPE = ethz-gnu64
# BUILD_TYPE = icsmaster-gnu64
	BUILD_TYPE = kaust
endif

# Default values which might be overridden
MAKE     = make
AR       = ar
ARFLAGS  = ruv
AR_FLAGS = ruv
RANLIB   = ranlib
# Local Directories
# Platform/compiler specific definitions

ifeq ($(BUILD_TYPE),ethz-gnu64)
	FLAGS         = -DAdd_ -DSOLVE_CRAY -D__MINMAX -Wno-write-strings
	GNU64	      = /usr/pack/gcc-9.2.0-af/linux-x64/
	MPIHOME	      = /usr/local/mpich-3.3.2/gcc/
	INTEL64	      = /usr/pack/intel_compiler-2015.0.090-af/Linux-x86_64/
	MKLHOME       = $(INTEL64)/mkl/lib/intel64
	BLAHOME	      = $(MKLHOME)
	LAPHOME       = $(MKLHOME)
	SCAHOME	      = $(MKLHOME)
	CUDAHOME      = /usr/local/cuda-10.0
	LOADER        = $(MPIHOME)/bin/mpic++ -Wl,-rpath -Wl,$(INTEL64)/compiler/lib/intel64
	MPICXX        = $(MPIHOME)/bin/mpic++ -Wl,-rpath -Wl,$(INTEL64)/compiler/lib/intel64
	NVCC	      = $(CUDAHOME)/bin/nvcc
	NVCC_LIBS     =
	CXXFLAGS      = -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK
	LOADFLAGS     = -O3 -ffast-math -funroll-loops -L$(GNU64)/lib64 -lgfortran -L$(INTEL64)/compiler/lib/intel64 -L$(MPIHOME)/lib/ -lmpifort
	NVCC_FLAGS    = -arch=compute_60 -code=sm_60
	CUDA_SDK      = 3
	COMMONFLAGS   = -DCUDA_SDK=$(CUDA_SDK)
	F90_LIBS      = -L$(GNU64)/lib64 -lgfortran -Wl,-rpath,L$(GNU64)/lib64
	LAPACK        = -L$(LAPHOME) -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath -Wl,$(MKLHOME)
	BLACS         = -L$(BLAHOME) -lmkl_blacs_sgimpt_lp64 -Wl,-rpath -Wl,$(MKLHOME)
	SCALAPACK     = -L$(SCAHOME) -lmkl_scalapack_lp64 -Wl,-rpath -Wl,$(MKLHOME)
	LINKS         = -Wl,-rpath -Wl,$(GNU64)/lib64
	LIBSCI        =
	BLAS          =
	OPENMP        = -fopenmp
	export LD_LIBRARY_PATH=$(GNU64)/lib64:$(INTEL64)/lib
endif

ifeq ($(BUILD_TYPE),icsmaster-gnu64)
	FLAGS         = -DAdd_ -DSOLVE_CRAY -D__MINMAX -Wno-write-strings
	GNU64	      = /apps/gcc/10.1.0
	INTEL64	      = /apps/spack_install/linux-centos8-haswell/gcc-10.1.0/intel-mkl-2020.1.217-qsctnr6bcemy7dlw6zav67jmdhge2fqi/compilers_and_libraries_2020.1.217/linux
	MKLHOME       = $(INTEL64)/mkl/lib/intel64
	BLAHOME	      = $(MKLHOME)
	LAPHOME       = $(MKLHOME)
	SCAHOME	      = $(MKLHOME)
	CUDAHOME      = /apps/cudatoolkit/11.1/gcc-10.1.0
	LOADER        = g++ -Wl,-rpath -Wl,$(INTEL64)/compiler/lib/intel64
	CXX           = g++ -Wl,-rpath -Wl,$(INTEL64)/compiler/lib/intel64
	NVCC	      = $(CUDAHOME)/bin/nvcc
	NVCC_LIBS     =
	CXXFLAGS      = -O3 -ffast-math -funroll-loops
	LOADFLAGS     = -O3 -ffast-math -funroll-loops -L$(GNU64)/lib64 -lgfortran -L$(INTEL64)/compiler/lib/intel64 #-lmpifort
	NVCC_FLAGS    = -arch=compute_60 -code=sm_60
	CUDA_SDK      = 3
	COMMONFLAGS   = -DCUDA_SDK=$(CUDA_SDK)
	F90_LIBS      = -L$(GNU64)/lib64 -lgfortran -Wl,-rpath,L$(GNU64)/lib64
	LAPACK        = -L$(LAPHOME) -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath -Wl,$(MKLHOME)
	BLACS         = #-L$(BLAHOME) -lmkl_blacs_sgimpt_lp64 -Wl,-rpath -Wl,$(MKLHOME)
	SCALAPACK     = -L/scratch/janalik/spack/opt/spack/linux-centos8-haswell/gcc-10.1.0/netlib-scalapack-2.1.0-qw6hespe7r5z3o4gxffrpvjq63pr2grh/lib -lscalapack #-L$(SCAHOME) -lmkl_scalapack_lp64 -Wl,-rpath -Wl,$(MKLHOME)
	LINKS         = -Wl,-rpath -Wl,$(GNU64)/lib64
	LIBSCI        =
	BLAS          =
	OPENMP        = -fopenmp
	INC_MAG        =
	MAGMA         = -lmagma
	export LD_LIBRARY_PATH=$(GNU64)/lib64:$(INTEL64)/lib
endif

ifeq ($(BUILD_TYPE),kaust)
	FLAGS         = -DAdd_ -DSOLVE_CRAY -D__MINMAX -Wno-write-strings
	MKLHOME       = $(MKLROOT)/lib/intel64
	BLAHOME       = $(MKLHOME)
	LAPHOME       = $(MKLHOME)
	SCAHOME       = $(MKLHOME)
	CUDAHOME      = /usr/local/cuda
#   $(INTEL64)/compiler/lib/intel64
	LOADER        = mpic++ -Wl,-rpath -Wl,$(INTEL64)/compiler/lib/intel64
	MPICXX        = mpic++ -Wl,-rpath -Wl,$(INTEL64)/compiler/lib/intel64
	NVCC          = nvcc
	NVCC_LIBS     =
	CXXFLAGS      = -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK
	LOADFLAGS     = -O3 -ffast-math -funroll-loops -L$(GNU64)/lib64 -lgfortran -L$(INTEL64)/compiler/lib/intel64 -L$(MPIHOME)/lib/ #-lmpifort
	NVCC_FLAGS     = -arch=compute_70 -code=sm_70
	CUDA_SDK      = 3
	COMMONFLAGS   = -DCUDA_SDK=$(CUDA_SDK)
	F90_LIBS      = -L$(GNU64)/lib64 -lgfortran -Wl,-rpath,L$(GNU64)/lib64
	LAPACK        = -L$(LAPHOME) -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath -Wl,$(MKLHOME)
	SCALAPACK     = -L/lib/x86_64-linux-gnu -lscalapack-openmpi #-L$(SCAHOME) -lmkl_scalapack_lp64 -Wl,-rpath -Wl,$(MKLHOME)
	LINKS         = -Wl,-rpath -Wl,$(GNU64)/lib64
	LIBSCI        =
	BLAS          =
	OPENMP        = -fopenmp
	export LD_LIBRARY_PATH=$(GNU64)/lib64:$(INTEL64)/lib
endif

INC_CC        = -Iinclude
#INC_CUD      = -I$(CUDAHOME)/include
#CUDA         = -L$(CUDAHOME)/lib64 -lcuda -L$(CUDAHOME)/lib64 -lcudart -lcublas -lcufft -lcusparse -Wl,-rpath,$(CUDAHOME)/lib64
INC_MAG       = -I$(MAGMA_DIR)/magma-2.5.4/include
MAGMA         = -L$(MAGMA_DIR)/magma-2.5.4/lib -lmagma

INC_CUD	      = -I$(CUDAHOME)/include
CUDA          = -L/usr/lib64/ -lcuda -L$(CUDAHOME)/lib64 -lcudart -lcublas -lcufft -lcusparse -Wl,-rpath,$(CUDAHOME)/lib64
#INC_MAG       ?= -I./magma-1.6.2p4/include/
#MAGMA        ?= -L./magma-1.6.2p4/lib/ -lmagma

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
