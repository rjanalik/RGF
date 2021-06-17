# RGF
## Installation
## Running the code 
### KAUST
1. Run an the MKL script to load the required MKL envrioment variables into our shell:
    ``` sh
    source /opt/intel/mkl/bin/mklvars.sh intel64
    ```
2. Add the magma library to the `ld` linked environment path variable:

    ``` sh
    export LD_LIBRARY_PATH=<path to library>/magma-2.5.4/lib:$LD_LIBRARY_PATH
    ```
3. Run (still a problem with my makefile):
    ``` sh
    make
    ```

