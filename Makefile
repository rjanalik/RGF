#include ./make.inc_kaust
include ./make.inc


LIBS	    = $(SCALAPACK) $(BLACS) $(LAPACK) $(BLAS) $(LINKS) $(OPENMP) $(CUDA) $(MAGMA) $(F90_LIBS)

SRC_DIR := src
INC_DIR := include
BUILD_DIR := build
BIN_DIR := $(BUILD_DIR)/bin
OBJ_DIR := $(BUILD_DIR)/obj
EXEC := $(BIN_DIR)/main
DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d

DEBUG ?= 1
ifeq ($(DEBUG), 1)
	DEBUG_FLAGS  =-DDEBUG -g -fsanitize=address,signed-integer-overflow -Wall
	DEBUG_FLAGS_NVCC=-DDEBUG -g
	CXXFLAGS     += -O0
else
	CXXFLAGS     += -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK
	DEBUG_FLAGS=-DNDEBUG
	DEBUG_FLAGS_NVCC=-DNDEBUG
endif



CC_SRC := $(wildcard $(SRC_DIR)/*.cpp)
CU_SRC := $(wildcard $(SRC_DIR)/*.cu)
CC_OBJ := $(CC_SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
# same as: CC_OBJ := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(CC_SRC))
CU_OBJ := $(CU_SRC:$(SRC_DIR)/%.cu=$(OBJ_DIR)/%.o)
# same as: CU_OBJ := $(patsubst $(SRC_DIR)/%.cu, $(OBJ_DIR)/%.o, $(CU_SRC))

.PHONY: all clean

all: $(EXEC)
$(EXEC): $(CC_OBJ) $(CU_OBJ)
	$(LOADER) $^ $(LINKFLAGS) $(LOADFLAGS) $(LIBS) -lm -o $@ $(DEBUG_FLAGS)

# Compile C++ source files to object files:
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR) $(BIN_DIR) $(DEPDIR)
	@echo "Compiling .cpp files"
	$(CXX) $(CXXFLAGS) $(FLAGS) $(OPEN_MP) $(INC_CC) $(INC_MAG) $(DEBUG_FLAGS) -c $< -o $@
# $(DEPFLAGS)

# Compile CUDA source files to object files:
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu | $(OBJ_DIR) $(BIN_DIR)
	@echo "Compiling .cu files"
	$(NVCC) $(NVCC_FLAGS) $(INC_CC) $(INC_CUD) $(INC_MAG) $(DEBUG_FLAGS_NVCC) -c $< -o $@ $(NVCC_LIBS)

$(BIN_DIR):
	@mkdir -p $@

$(OBJ_DIR):
	@mkdir -p $@

$(DEPDIR):
	@mkdir -p $@

DEPFILES := $(SRCS:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):
	include $(wildcard $(DEPFILES))

clean:	
	@$(RM) -rv $(BUILD_DIR)
