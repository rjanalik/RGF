#!/bin/bash

source /opt/intel/oneapi/mkl/latest/env/vars.sh  intel64
export LD_LIBRARY_PATH=~/RGF/applications/magma-2.5.4/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=~/RGF/applications/magma-2.5.4/lib:$LIBRARY_PATH
export PATH=~/RGF/applications/magma-2.5.4/bin:$PATH
export CPATH=~/RGF/applications/magma-2.5.4/include:$CPATH
export MAGMA_DIR=~/RGF/applications
