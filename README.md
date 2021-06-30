# RGF
## TODO:
- Check version conflicts:  
`warning: libcusparse.so.11, needed by /home/x_pollakgr/applications/magma-2.5.4/lib/libmagma.so, may conflict with libcusparse.so.10`
- CUDA header dependencies not in makefile included
## Installation
## Fetching the data & Creating the Model
### Fetching the data 
1. Change to:
    ``` sh
    <project root>/scripts/R
    ```
2. Run: 
   ```sh
   Rscript fetch_ghcn_daily.R [--help]
   ```
### Creating The Models
1. Change to:
    ``` sh
    <project root>/scripts/R
    ```
2. Run: 
   ```sh
   Rscript generate_spatial_temporal_model.R [--help]
   ```
## Running the solver
### KAUST
1. Run an the MKL script to load the required MKL envrioment variables into our shell:
    ``` sh
    source /opt/intel/oneapi/mkl/latest/env/vars.sh  intel64
    ```
2. Add the magma library to the `ld` linked environment path variable:

    ``` sh
    export LD_LIBRARY_PATH=<path to library>/magma-2.5.4/lib:$LD_LIBRARY_PATH
    ```
3. Set the local MAGMA root directory path (needed for make.inc file):
   ```sh
     export MAGMA_DIR=/home/x_pollakgr/applications
   ```
3. Run the makefile:
    ``` sh
    make
    ```
