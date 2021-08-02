#include ./make.inc_kaust
#         $(CXX) -c $< $(CXXFLAGS) $(DEBUG) $(OPENMP) $(FLAGS) $(INCLUDEDIR)
include ./make.inc

DEBUG       =     #-g -fsanitize=address,signed-integer-overflow

INCLUDEDIR  = $(INCCUD) $(INCMAG)
LIBS	    = $(SCALAPACK) $(BLACS) $(LAPACK) $(BLAS) $(LINKS) $(OPENMP) $(CUDA) $(MAGMA) $(F90_LIBS)

CC_FILES   = Utilities.o main.o
CU_FILES   = CWC_utility.o

mainConstInd: $(CC_FILES) $(CU_FILES)
	$(LOADER) $(CC_FILES) $(CU_FILES) $(LOADFLAGS) $(LIBS) -lm -o $@ $(DEBUG)

.C.o:
	$(CXX) -c $< $(CXXFLAGS) $(DEBUG) $(FLAGS) $(INCLUDEDIR)

CWC_utility.o: CWC_utility.cu
	$(NVCC) -c CWC_utility.cu $(NVCCFLAGS) $(INCLUDEDIR)

clean:	
	rm -f *.o *.c *.h RGFSolver
