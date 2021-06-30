#!/bin/bash

source /opt/intel/oneapi/mkl/latest/env/vars.sh  intel64
export LD_LIBRARY_PATH=~/RGF/external/magma-2.5.4/lib:$LD_LIBRARY_PATH
export MAGMA_DIR=~/RGF/external
