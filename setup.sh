#!/bin/bash
#ml gcc/10.1.0 intel-mkl/2020.1.217-gcc-10.1.0-qsctnr6 cudatoolkit/11.1
ml gcc/10.1.0 intel-mkl/2020.1.217-gcc-10.1.0-qsctnr6
. /scratch/janalik/spack/share/spack/setup-env.sh
spack load cuda@11.2.1
spack load netlib-scalapack@2.1.0
spack load magma@2.5.4
#spack load armadillo@9.800.3
spack unload cuda@11.1.0
