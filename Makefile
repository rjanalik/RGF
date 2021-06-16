#include ./make.inc_kaust
#         $(CXX) -c $< $(CXXFLAGS) $(DEBUG) $(OPENMP) $(FLAGS) $(INCLUDEDIR)
include ./make.inc

DEBUG       =     #-g -fsanitize=address,signed-integer-overflow

INCLUDEDIR  = $(INCCUD) $(INCMAG)
LIBS	    = $(SCALAPACK) $(BLACS) $(LAPACK) $(BLAS) $(LINKS) $(OPENMP) $(CUDA) $(MAGMA) $(F90_LIBS)

SRC_DIR := src
INC_DIR := include
BUILD_DIR := build
BIN_DIR := $(BUILD_DIR)/bin
OBJ_DIR := $(BUILD_DIR)/obj
EXEC := $(BIN_DIR)/main

# CC_FILES   = Utilities.o main.o
# CU_FILES   = CWC_utility.o
# TODO in order to have them as dependencies for all/mainConstInd
CC_SRC := $(wildcard $(SRC_DIR)/*.c)
CU_SRC := $(wildcard $(SRC_DIR)/*.cu)
CC_OBJ := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(CC_SRC))
CU_OBJ := $(patsubst $(SRC_DIR)/%.cu, $(OBJ_DIR)/%.o, $(CU_SRC))

.PHONY: all clean

all: $(EXEC)
$(EXEC): $(CC_OBJ) $(CU_OBJ) | $(BIN_DIR) $(OBJ_DIR)
	$(LOADER) $^ $(LOADFLAGS) $(LIBS) -lm -o $@ $(DEBUG)

# %.o: %.c:
# 	$(CXX) -c $< $(CXXFLAGS) $(DEBUG) $(FLAGS) $(INCLUDEDIR)

# Compile C++ source files to object files:
$(CC_OBJ): $(CC_SRC)
	$(CXX) $(CXXFLAGS) $(DEBUG) $(FLAGS) $(INCLUDEDIR) -c $< -o $@

# Compile CUDA source files to object files:
$(CU_OBJ): $(CU_SRC)
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS)

$(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR):
	mkdir -p $@

clean:	
	$(RM) -rv $(BUILD_DIR)
