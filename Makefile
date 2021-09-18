#include ./make.inc_kaust
#         $(CXX) -c $< $(CXXFLAGS) $(DEBUG) $(OPENMP) $(FLAGS) $(INCLUDEDIR)
include ./make.inc

SRC_DIR := src
INC_DIR := include
INC_CC  := -I$(INC_DIR)
C_EXTEN := cpp
CU_EXTEN := cu
BUILD_DIR := build
BIN_DIR := $(BUILD_DIR)/bin
OBJ_DIR := $(BUILD_DIR)/obj
EXEC := main

OBJ_FILES := ModelGenerator.o utilities.o main.o
CC_FILES := $(pathsubst %.o, %.$(C_EXTEN), $(OBJ_FILES))
CU_FILES := $(CWC_utilitiy.cu)

CC_SRC := $(addprefix $(SRC_DIR)/,$(CC_FILES))
$(info $(CC_SRC))
CU_SRC := $(addprefix $(SRC_DIR)/,$(CU_FILES))

CC_OBJECTS:= $(addprefix $(OBJ_DIR)/,$(OBJ_FILES))
CU_OBJECTS:= $(addprefix $(OBJ_DIR)/,CWC_utility.o)
$(info $(CC_OBJECTS))
$(info $(CU_OBJECTS))

INCLUDEDIR  = $(INCCUD) $(INCMAG)
LIBS	    = $(SCALAPACK) $(BLACS) $(LAPACK) $(BLAS) $(LINKS) $(OPENMP) $(CUDA) $(MAGMA) $(F90_LIBS)

CC_FILES   = utilities.o main.o ModelGenerator.o
CU_FILES   = CWC_utility.o

DEBUG ?= 1
ifeq ($(DEBUG), 1)
	DEBUG_FLAGS  +=-DDEBUG -g -Wall -fsanitize=address,signed-integer-overflow
	DEBUG_FLAGS_NVCC+=-DDEBUG -g
	CXXFLAGS     = -O0
else ifeq ($(DEBUG), 2)
	DEBUG_FLAGS  +=-DDEBUG -g -Wall #-fno-stack-protector
	CXXFLAGS     = -O0
	DEBUG_FLAGS_NVCC=-DDEBUG -g
else
	CXXFLAGS     = -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK
	NVCC_FLAGS   +=
	DEBUG_FLAGS=-DNDEBUG
	DEBUG_FLAGS_NVCC=-DNDEBUG
endif

.PHONY: clean

$(EXEC): $(CC_OBJECTS) $(CU_OBJECTS)
	$(LOADER) $(CC_OBJECTS) $(CU_OBJECTS) $(LOADFLAGS) $(LIBS) -lm -o $(BIN_DIR)/$@

# canned old version of:  %.o : %.c
# $(CC_OBJECTS): $(CC_SRC)
# 	@echo "Compiling .$(C_EXTEN) files"
# 	$(CXX) -c $< $(CXXFLAGS) $(INC_CC)  $(FLAGS) $(INCLUDEDIR)
build/obj/main.o: src/main.cpp | $(OBJ_DIR) $(BIN_DIR)
	@echo "Compiling .$(C_EXTEN) files"
	$(CXX) -c $< $(CXXFLAGS) $(INC_CC)  $(FLAGS) $(INCLUDEDIR) -o $@

build/obj/utilities.o: src/utilities.cpp include/utilities.hpp | $(OBJ_DIR) $(BIN_DIR)
	@echo "Compiling .$(C_EXTEN) files"
	$(CXX) -c $< $(CXXFLAGS) $(INC_CC)  $(FLAGS) $(INCLUDEDIR) -o $@

build/obj/ModelGenerator.o: src/ModelGenerator.cpp include/ModelGenerator.hpp | $(OBJ_DIR) $(BIN_DIR)
	@echo "Compiling .$(C_EXTEN) files"
	$(CXX) -c $< $(CXXFLAGS) $(INC_CC)  $(FLAGS) $(INCLUDEDIR) -o $@

build/obj/CWC_utility.o: src/CWC_utility.cu | $(OBJ_DIR) $(BIN_DIR)
	@echo "Compiling .$(CU_EXTEN) files"
	$(NVCC) -c $< $(INC_CC) $(NVCCFLAGS) $(INCLUDEDIR) -o $@


$(BIN_DIR):
	@mkdir -p $@

$(OBJ_DIR):
	@mkdir -p $@

clean:
	@$(RM) -rv $(BUILD_DIR)
