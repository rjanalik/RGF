include ./make.inc

DEBUG       =     

INCLUDEDIR  = $(INCCUD) $(INCMAG)
LIBS	    = $(SCALAPACK) $(BLACS) $(LAPACK) $(BLAS) $(LINKS) $(OPENMP) $(MPI) $(CUDA) $(MAGMA) $(F90_LIBS)

CC_FILES   = Utilities.o RGF.o main.o
CU_FILES   = CWC_utility.o

RGFSolver: $(CC_FILES) $(CU_FILES)
	$(LOADER) $(CC_FILES) $(CU_FILES) $(LOADFLAGS) $(LIBS) -lm -o $@$

.C.o:
	$(MPICXX) -c $< $(CXXFLAGS) $(DEBUG) $(FLAGS) $(INCLUDEDIR)

CWC_utility.o: CWC_utility.cu
	$(NVCC) -c CWC_utility.cu $(NVCCFLAGS) $(INCCUD)

clean:	
	rm -f *.o *.c *.h RGFSolver
