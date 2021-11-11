# compile into shared library and then include this in main makefile

MKLHOME       = $(MKLROOT)/lib/intel64
BLAHOME       = $(MKLHOME)
LAPHOME       = $(MKLHOME)
SCAHOME       = $(MKLHOME)

#include ./make.inc_kaust
#         $(CXX) -c $< $(CXXFLAGS) $(DEBUG) $(OPENMP) $(FLAGS) $(INCLUDEDIR)
include ./make.inc

DEBUG       =   # -Wall -g -fsanitize=address,signed-integer-overflow

# for mainEigen
INCEIGEN=-I/usr/include/eigen3

INCLUDEDIR  = $(INCCUD) $(INCMAG) $(INCEIGEN)
LIBS	    = $(SCALAPACK) $(BLACS) $(LAPACK) $(BLAS) $(LINKS) $(OPENMP) $(CUDA) $(MAGMA) $(F90_LIBS)


CC_FILES   = Utilities.o mainEigen.o #main.o
CU_FILES   = CWC_utility.o

mainEigen: $(CC_FILES) $(CU_FILES)
	$(RGF_LOADER) $(CC_FILES) $(CU_FILES) $(RGF_LOADFLAGS) $(LIBS) -lm -o $@ $(DEBUG)

#main: $(CC_FILES) $(CU_FILES)
#	$(LOADER) $(CC_FILES) $(CU_FILES) $(RGF_LOADFLAGS) $(LIBS) -lm -o $@ $(DEBUG)

#librgf.a: $(CC_FILES) $(CU_FILES)
#	ar rvs $@ $^

#librgf.so: $(CC_FILES) $(CU_FILES)
#	$(LOADER) --shared $(CC_FILES) $(CU_FILES) $(LOADFLAGS) $(LIBS) -lm -o $@ $(DEBUG) -fPIC

#main: $(CC_FILES) $(CU_FILES)
#	$(RGF_LOADER) $(CC_FILES) $(CU_FILES) $(RGF_LOADFLAGS) $(OPENMP) $(LIBS) -lm -o $@ $(DEBUG)

.C.o:
	$(CXX) -c $< $(RGF_CXXFLAGS) $(DEBUG) $(RGF_FLAGS) $(INCLUDEDIR) $(OPENMP) -fPIC

CWC_utility.o: CWC_utility.cu
	$(NVCC) -c CWC_utility.cu $(NVCCFLAGS) $(INCCUD) --compiler-options '-fPIC' 

clean:	
	rm -f *.o RGFSolver
	rm -f mainEigen
	rm -f main
	#rm -f librgf.so
	#rm -f librgf.a
