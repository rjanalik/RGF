# RGF
    .
    ├── docs                    # Documentation with `index.html`
    ├── build                   # Build files
    │   ├── obj                 # Object Files
    │   └── bin                 # Binary 
    ├── data                    # Data files 
    │   ├── input               # Input files/matrices
    │   │   └── ghcn            # Global Historical Climatology Network (GHCN)
    │   │        ├── 2018 
    │   │        ├── 2019
    │   │        └── ...
    │   └── output              # Output Files From ./main
    ├── external                # External Libraries s.a. magma
    ├── include                 # Header Files
    ├── lib                     # Libraries
    ├── scripts                 # Scripts i.e. for setting up the environment
    │   ├── R                   # R script, mainly for data generation
    │   └── ...
    └── src                    # Source files
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
### Compiling The Code
#### KAUST
1. Setup the enviroment:
- Run the setup script:
   ```shell
   source ~/RGF/scripts/setup_kaust.sh
   ```
- or alternatively:
    1. Run an the MKL script to load the required MKL envrioment variables into our shell:
        ``` sh
        source /opt/intel/oneapi/mkl/latest/env/vars.sh  intel64
        ```
    2. Add the magma library to the `ld` linked environment path variable:
        ``` sh
        export LD_LIBRARY_PATH=<path to library>/magma-2.5.4/lib:$LD_LIBRARY_PATH
        ```
    3. Set the local MAGMA root directory path i.e. where MAGMA directory resides (needed for make.inc file):
    ```sh
        export MAGMA_DIR=<path to library>/external
    ```
2. Run the makefile:
    ``` sh
    make
    ```
### Running the code 
#### KAUST
    ```sh
    CUDA_VISIBLE_DEVICES="nb" ./main <data_folder_path> ns nt nb no > <folder_path/output.txt>
    ```
    i.e.
    ```sh
    CUDA_VISIBLE_DEVICES=0 ./build/bin/main -p data/input/ghcn/2019/spatio_temporal/ns42_nt3 -s 42 -t 3 -f 2 -n 35542
    ```
