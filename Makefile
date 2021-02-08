include ./make.inc

DEBUG       =     

INCLUDEDIR  = $(INCCUD) $(INCMAG)
LIBS	    = $(SCALAPACK) $(BLACS) $(LAPACK) $(BLAS) $(LINKS) $(OPENMP) $(MPI) $(CUDA) $(MAGMA) $(F90_LIBS)

CC_FILES   = Utilities.o main.o
CU_FILES   = CWC_utility.o

RGFSolver: $(CC_FILES) $(CU_FILES)
	$(LOADER) $(CC_FILES) $(CU_FILES) $(LOADFLAGS) $(LIBS) -lm -o $@ -g

.C.o:
	$(MPICXX) -c $< $(CXXFLAGS) $(DEBUG) $(FLAGS) $(INCLUDEDIR) -g 

CWC_utility.o: CWC_utility.cu
	$(NVCC) -c CWC_utility.cu $(NVCCFLAGS) $(INCCUD) -g

clean:	
	rm -f *.o *.c *.h RGFSolver
